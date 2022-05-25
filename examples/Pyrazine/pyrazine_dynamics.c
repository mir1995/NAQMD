/* dynamics.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include <time.h> // srand(time(0))
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../src/Odeint/odeint.h"
#include "../../src/Potential/potential.h"
#include "../../setup.h" 


/*
# define EPS 0.074 / 27.211378682468403 
# define a1 1.703  / 27.211378682468403 
# define a2 1.0/ 27.211378682468403 
# define a3 1.595/ 27.211378682468403 
*/
# define EPS 0.074  
# define w1 1.703  
# define w2 1.0  
# define w3 1.595 
//# define w1 1.703 
//# define w2 1.0 
//# define w3 1.595 

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  long int npart;
  int dim, t, s;
  double dt, q[DIM], p[DIM], std[2*DIM];
  double param[1];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim = DIM; 
  t = 37, dt = 0.001;
  q[0] = 0, q[1] = 0, q[2] = 0;
  p[0]=0, p[1]=0, p[2] = 0;
  std[0] = sqrt(EPS/2.0/w1), std[1] =sqrt(EPS/2.0/w2), std[2] =sqrt(EPS/2.0/w3); 
  std[3] = sqrt(EPS/2.0*w1), std[4] =sqrt(EPS/2.0*w2), std[5] =sqrt(EPS/2.0*w3); 
  npart = 500; //pow(10,3);
  s = atoi(argv[2]);
  param[0] = EPS;
  char *rate = argv[1];
  /*
   *  PRINT SIMULATION PARAMETERS
  */
  printf("--->  w1    :  %.4f\n", w1);
  printf("--->  eps    :  %.4f\n", EPS);

  printf(" ***************************************************************\n");
  printf(" Pyrazine: Single Switch Surface Hopping  \n");
  printf("---> NPART     :  %ld\n", npart);
  printf("---> DIM       :  %d\n", dim);
  printf("---> T         :  %d\n", t);
  printf("---> dt        :  %.4f\n", dt);
  printf("---> POTENTIAL :  Pyrazine \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(npart, dim); // it's a weird data struct - to improve
  struct Odeint *solver = odeint_new(t, dt, dim, "lietrotter_symplectic");
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "Pyrazine", param);

  struct Observables *observables = sh_observables_new(npart, dim);
  struct Hopper *hopper = sh_hopper_new(rate);

  char filename[200];
  sprintf(filename, "./data/observables_npart%ld_rate%s_seed%d.txt", npart, rate, s); 
  file = fopen(filename, "w"); 


  /*
   *  SAMPLE PARTICLES FROM ITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
   */
  printf(" ***************************************************************\n");
  printf("    INITIALISE PARTICLES FROM WIGNER DISTRIBUTION   \n");
  printf("    INITIALISE PARTICLES' DATA: POS, MOM, POT, ...   \n");
  printf(" ***************************************************************\n");

  srand(s); // INITIALISE SEED // what is the actual variance...?
  sh_wigner_fill(particles, q, p, std, npart, dim);
  
  sh_particle_potential_init(particles, pot, npart, dim);// initialise particle values - potential, gradient, level ...

  /*
   *  SIMULATION: STEP - CHECK FOR AVOIDED CROSSING - UPDATE POTENTIAL - STEP
   *  AT LEAST ONE SIMULATION CHECKING SURFACE HOPPING TO WAVEPACKET DYNAMICS
   */
  
  printf(" ***************************************************************\n");
  printf("   Surface Particle Hopping Simulation  \n");
  printf(" ***************************************************************\n");
  
  fprintf(file, "t \t pos_up \t pos_down \t \
      mom_up \t mom_down \t ke_up \t ke_down \t \
      e_up \t e_down \t mass_up \t mass_down \n");

  for (int itr=0; itr < (int)((t*1.)/dt); itr++){

    struct Particle *part = particles; // come up with a better structure than a linked list

    for(long int i=0; i<npart; i++, part++){
      // 2) STEP SOLUTION IN TIME 
      solver->func_dostep(solver, part, pot);
      // 3) UPDATE PARTICLE INFORMATION 
      sh_particle_potential_update(part, pot, dim);
      // 4) CALL HOPPER 
      hopper->func_hop(part, hopper, pot, solver);
    }
    if(itr % (int)(((t*1.)/dt)/500) == 0){
      sh_observables_update(observables, particles);
      fprintf(file, "%.17g \t %.17g \t %.17g \t  %.17g \
                          \t %.17g \t %.17g \t %.17g \
                          \t %.17g \t %.17g \t %.17g \
                          \t %.17g \n",
          itr * dt,
          observables->x_up[0], observables->x_down[0],
          observables->p_up[0], observables->p_down[0],
          observables->ke_up, observables->ke_down,
          observables->e_up, observables->e_down,
          observables->mass_up, observables->mass_down
          );
    }
  }
  
  fclose(file);
  return 0;

}

