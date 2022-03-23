/* dynamics.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../src/Odeint/odeint.h"
#include "../../src/Potential/potential.h"
#include "../../setup.h"
#include <time.h> // srand(time(0))


#define ALPHA 0.5 // might need vary your parameter space
#define DELTA 0.6 
#define EPS 0.1

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  int npart;
  int dim, t, s;
  double dt;
  double q[DIM], p[DIM];
  double param[3];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim = DIM; 
  q[0] = -10, p[0] = 4;
  t = 20/p[0], dt = 1.0 / 100 / p[0];
  npart = pow(10,7);
  s = 1;
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  char *rate = argv[1];
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" Simple I (model study): Single Switch Surface Hopping  \n");
  printf("---> NPART     :  %d\n",npart);
  printf("---> DIM       :  %d\n",dim);
  printf("---> T         :  %d\n", t);
  printf("---> dt        :  %.4f\n", dt);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  Simple1 \n"); 
  printf("---> RATE      :  %s\n", rate);
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(npart, dim); // it's a weird data struct - to improve
  struct Odeint *solver = odeint_new(t, dt, dim, "lietrotter_symplectic");
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "Simple I", param);
  struct Observables *observables = sh_observables_new(npart, dim);
  struct Hopper *hopper = sh_hopper_new(rate);

  char filename[200];
  sprintf(filename, "./data/observables_alpha%.4f_delta%.4f_npart%d_rate%s.txt", ALPHA, DELTA, npart, rate); 
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
  
  fprintf(file, "itr \t pos_up \t pos_down \t \
      mom_up \t mom_down \t ke_up \t ke_down \t \
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
    if(itr % (int)(((t*1.)/dt)/40) == 0){
      sh_observables_update(observables, particles);
      fprintf(file, "%d \t %.17g \t %.17g \t  %.17g \
                          \t %.17g \t %.17g \t %.17g \
                          \t %.17g \t %.17g \t %.17g \
                          \t %.17g \n",
          itr,
          observables->x_up[0], observables->x_down[0],
          observables->p_up[0], observables->p_down[0],
          observables->ke_up, observables->ke_down,
          observables->e_up, observables->e_down,
          observables->mass_up, observables->mass_down
          );
    }
  }
#ifdef ALTO
#endif
  
  fclose(file);
  return 0;
}

