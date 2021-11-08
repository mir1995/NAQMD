/* dynamics.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../src/Odeint/odeint.h"
#include "../../src/Potential/potential.h"
#include <time.h> // srand(time(0))
#include "../../setup.h" 

#define GAMMA 0
#define EPS 0.01
#define DIM 2

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  unsigned int npart, dim;
  int t, s;
  double dt;
  double q[DIM], p[DIM];
  double param[4];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim = DIM; 
  t = 6, dt = 0.001;
  q[0]=5*sqrt(EPS), q[1]=0.5*sqrt(EPS), p[0]=0, p[1]=0;
  npart = atoi(argv[1]); // given that you know the convergence rate of the 
  // two methods, what is the corresponding number of points that 
  // accuracy
  s = 1;
  param[0] = EPS;
  param[3] = GAMMA;
  char *rate = argv[2];
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" Jahn Teller: Probabilistic Single Switch Surface Hopping  \n");
  printf("---> NPART     :  %d\n",npart);
  printf("---> DIM       :  %d\n",dim);
  printf("---> T         :  %d\n", t);
  printf("---> dt        :  %.4f\n", dt);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  Jahn Teller \n"); 
  printf("---> RATE      :  %s\n", rate);
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(npart, dim); // i think there is actually no need for the pointer 
  
  struct Odeint *solver = odeint_new(t, dt, dim, "lietrotter_symplectic");
  // thinking I do not need to pass this as pointer
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "Jahn Teller", param);
  struct Observables *observables = sh_observables_new(npart, dim);
  struct Hopper *hopper = sh_hopper_new(rate);

  char filename[200];
  sprintf(filename, "./data/observables_gamma%dseed%d_npart%d_rate%s.txt", GAMMA, s, npart, rate); 
  file = fopen(filename, "w"); 


  /*
   *  SAMPLE PARTICLES FROM ITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
   */
  printf(" ***************************************************************\n");
  printf("    INITIALISE PARTICLES FROM WIGNER DISTRIBUTION   \n");
  printf("    INITIALISE PARTICLES' DATA: POS, MOM, POT, ...   \n");
  printf(" ***************************************************************\n");

  srand(s); // INITIALISE SEED // what is the actual variance...?
  sh_wigner_fill(particles, q, p, sqrt(EPS/2), npart, dim);
  sh_particle_potential_init(particles, pot, npart, dim);// initialise particle values - potential, gradient, level ...

  /*
   *  SIMULATION: STEP - CHECK FOR AVOIDED CROSSING - UPDATE POTENTIAL - STEP
   *  AT LEAST ONE SIMULATION CHECKING SURFACE HOPPING TO WAVEPACKET DYNAMICS
   */
  
  printf(" ***************************************************************\n");
  printf("   Surface Particle Hopping Simulation  \n");
  printf(" ***************************************************************\n");
  
  fprintf(file, "itr \t x_up[0] \t x_up[1] \t x_down[0] \t x_down[1] \t \
      p_up[0] \t p_up[1] \t p_down[0] \t p_down[1] \t ke_up \t ke_down \t \
      e_up \t e_down \t mass_up \t mass_down \n");

  for (int itr=0; itr < (int)((t*1.)/dt); itr++){

    struct Particle *part = particles;
    for(int i=0; i<npart; i++, part++){
      // 2) STEP SOLUTION IN TIME 
      solver->func_dostep(solver, part, pot);
      // 3) UPDATE PARTICLE INFORMATION 
      sh_particle_potential_update(part, pot, dim);
      // 4) CALL HOPPER 
      hopper->func_hop(part, hopper, pot, solver);
    }
    if(itr % (int)(((t*1.)/dt)/80) == 0){
      sh_observables_update(observables, particles);
      fprintf(file, "%d \t %.17g \t %.17g \t %.17g \t %.17g \t  %.17g \
                         \t %.17g \t %.17g \t %.17g \t %.17g \t %.17g \
                          \t %.17g \t %.17g \t %.17g \
                          \t %.17g \n",
          itr,
          observables->x_up[0], observables->x_up[1],
          observables->x_down[0], observables->x_down[1],
          observables->p_up[0], observables->p_up[1],
          observables->p_down[0], observables->p_down[1],
          observables->ke_up, observables->ke_down,
          observables->e_up, observables->e_down,
          observables->mass_up, observables->mass_down
          );
    }
  }
  
  fclose(file);
  return 0;
}
