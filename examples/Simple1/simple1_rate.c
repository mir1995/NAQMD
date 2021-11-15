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


#define ALPHA 0.05 //0.5 //0.05 // might need vary your parameter space
#define DELTA 0.5 
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
  npart = pow(10,6);
  s = 1;
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" Simple I (model study): Energy conservation  \n");
  printf("---> NPART     :  %d\n",npart);
  printf("---> DIM       :  %d\n",dim);
  printf("---> T         :  %d\n", t);
  printf("---> dt        :  %.4f\n", dt);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  Simple1 \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(npart, dim); // it's a weird data struct - to improve
  struct Odeint *solver = odeint_new(t, dt, dim, "lietrotter_symplectic");
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "Simple I", param);

  char filename[200];
  sprintf(filename, "./data/transition_rate_alphadelta%.5f_npart%d.txt", ALPHA / DELTA, npart); 
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
  
  fprintf(file, "SA - SA1 \t SA3 - SA13 \n");

  double x_c[1] = {0.0};
  struct Particle *part = particles;

  for(int i=0; i<npart; i++, part++){
    part->p_curr[0] = sqrt(pow(part->p[0], 2) + \
        2 * (pot->func_potup(pot, part->x_new, 1) - pot->func_potup(pot, x_c, 1)));
    part->x_curr[0] = x_c[0]; // the crossing is at zero
    part->rho_curr = DELTA;
    double sa = sh_transition_sa(part, pot, solver);
    double sa1 = sh_transition_sa1(part, pot, solver);
    double sa3 = sh_transition_sa3(part, pot, solver);
    double sa13 = sh_transition_sa13(part, pot, solver);
    
    printf("sa %.17g sa1 %.17g sa3 %.17g sa13 %.17g \n", sa, sa1 , sa3, sa13);

    fprintf(file, "%.17g \t %.17g \n",
        (sh_transition_sa(part, pot, solver) - sh_transition_sa1(part, pot, solver))/sh_transition_sa(part, pot, solver),
        (sh_transition_sa3(part, pot, solver) - sh_transition_sa13(part, pot, solver))/sh_transition_sa3(part, pot, solver));
  }
  
  fclose(file);
  #ifdef ALRO
  #endif  
  return 0;
}

