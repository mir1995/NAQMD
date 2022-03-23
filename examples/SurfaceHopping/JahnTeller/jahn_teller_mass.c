/* many.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../setup.h"
#include "../../src/Odeint/odeint.h"
#include "../../src/Potential/potential.h"
#include <time.h> // srand(time(0))


#define GAMMA 3
#define ALPHA 0.5 
#define DELTA 0.5
#define EPS 0.01
#define DIM 2

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  unsigned long int npart; 
  unsigned int dim;
  int s;
  double q[DIM], p[DIM];
  double param[4];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim =DIM; 
  q[0]=5*sqrt(EPS), q[1]=0.5*sqrt(EPS), p[0]=0, p[1]=0;
  npart = pow(10,4); // given that you know the convergence rate of the 
  // 
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  param[3] = GAMMA;
  s = 1;
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" Jahn Teller: No dynamics - estimating number of particles  \n");
  printf("---> NPART     :  %ld\n",npart);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  NaI \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(pow(10,4), dim); 
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "JahnTeller", param);
  struct Odeint *solver = odeint_new(20, 0.001, 1, "lietrotter_symplectic"); // not needed
  
  char filename[200];
  sprintf(filename, "data/mass_transitioned_seed%d_npart%ld.txt", s, npart);
  sprintf(filename, "data/testing_seed%d_npart%ld.txt", s, npart);
  file = fopen(filename, "a"); 
 
  fprintf(file, "rate \t mass \n");
  
  /*
   *  SAMPLE PARTICLES FROM INITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
   */
  printf(" ***************************************************************\n");
  printf("    INITIALISE PARTICLES FROM WIGNER DISTRIBUTION   \n");
  printf("    INITIALISE PARTICLES' DATA: POS, MOM, POT, ...   \n");
  printf(" ***************************************************************\n");

  srand(s); // INITIALISE SEED // what is the actual variance...?
  int count_lzdia = 0; // count number of particles which have transitioned
  int count_lzadia = 0;
  // the crossing has been pre-computed
  double x_c[2];
  x_c[0] = 0;
  x_c[1] = 0;
  

  for (int i=0; i< (int) (npart/pow(10,4)); i++){
    sh_wigner_fill(particles, q, p, sqrt(EPS/2), pow(10,4), dim);  
    sh_particle_potential_init(particles, pot, pow(10,4), dim);// initialise particle values - potential, gradient, level ...
    struct Particle *part = particles; // come up with a better structure than a linked list
    for(unsigned int i=0; i<pow(10,4); i++, part++){
      // energy conservation
      // the following is missing something
      part->p_curr[0] = sqrt(pow(part->p[0], 2) + \
          2 * (pot->func_potup(pot, part->x, dim) - pot->func_potup(pot, x_c, dim)));
      part->p_curr[1] = 0; 
      part->x_curr[0] = x_c[0];
      part->x_curr[1] = x_c[1];
      part->rho_curr = pot->func_rho(pot, x_c, dim);
      printf("%.17g\n", part->rho_curr); 
      double pr = ((double)rand() / RAND_MAX); 
      if (sh_transition_lzdia(part, pot, solver)>= pr){ 
        //printf("%.17g\n", sh_transition_lzdia(part, pot, solver)); 
        count_lzdia += 1;
      }
      /*
      if (sh_transition_lzadia(part, pot, solver)>= pr){ 
        count_lzadia += 1;
      }
      */
    }
  }
  fprintf(file, "%s \t %.17g \n", "lz_dia", count_lzdia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "lz_adia", count_lzadia * 1.0 / npart);
  return 0;
}

