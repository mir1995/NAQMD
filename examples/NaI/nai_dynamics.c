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


// alpha is used in the approximation to the transition rate
#define ALPHA 0.002669794983146837
#define DELTA 0.002010562925688143 
#define EPS 0.014654629670711006

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  long int npart;
  int dim, t, s;
  double dt, q[DIM], p[DIM];
  double param[3];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim =DIM; 
  t = 90, dt = 0.0025;
  q[0] = 5, p[0] = 0;
  npart = atoi(argv[1]);
  s = 1;
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" NaI case study: Single Switch Surface Hopping  \n");
  printf("---> NPART     :  %ld\n",npart);
  printf("---> DIM       :  %d\n",dim);
  printf("---> T         :  %d\n", t);
  printf("---> dt        :  %.4f\n", dt);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  NaI \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Odeint *solver = odeint_new(t, dt, dim, "lietrotter_symplectic");
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "NaI", param);

  char filename[200];
  sprintf(filename, "./data/particles_crossing_npart%ldt%ddt%.17g.txt", npart, t, dt); 
  file = fopen(filename, "a"); 


  /*
   *  SAMPLE PARTICLES FROM ITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
   */
  printf(" ***************************************************************\n");
  printf("    INITIALISE PARTICLES FROM WIGNER DISTRIBUTION   \n");
  printf("    INITIALISE PARTICLES' DATA: POS, MOM, POT, ...   \n");
  printf(" ***************************************************************\n");

  srand(s); // INITIALISE SEED // what is the actual variance...?

  /*
   *  SIMULATION: STEP - CHECK FOR AVOIDED CROSSING - UPDATE POTENTIAL - STEP
   *  AT LEAST ONE SIMULATION CHECKING SURFACE HOPPING TO WAVEPACKET DYNAMICS
   */
  
  printf(" ***************************************************************\n");
  printf("   Surface Particle Hopping Simulation  \n");
  printf(" ***************************************************************\n");
  
  fprintf(file, "t_c \t x(t_c) \t p(t_c) \n");
    
  // the crossing has been pre-computed
  double x_c[1] = {13.27801894097567};
  struct Particle *particle = sh_particles_create(1, dim); // it's a weird data struct - to improve
  
  for (long int i=0; i<npart; i++){
      
      int itr = 0;

      sh_wigner_fill(particle, q, p, sqrt(EPS/2), 1, dim);
      
      sh_particle_potential_init(particle, pot, 1, dim);// initialise particle value
      
      while (particle->x_new[0] <= x_c[0]){
          // 2) STEP SOLUTION IN TIME 
          solver->func_dostep(solver, particle, pot);
          // 3) UPDATE PARTICLE INFORMATION 
          sh_particle_potential_update(particle, pot, dim);
          itr+=1;
      }
      
      fprintf(file, "%d \t %.17g \t %.17g \n",
              itr, particle->x_curr[0], particle->p_curr[0]);
  }

  fclose(file);
  return 0;
}

