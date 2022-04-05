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


// alpha is used in the approximation to the transition rate
#define ALPHA 0.002669794983146837
#define DELTA 0.002010562925688143 
#define EPS 0.014654629670711006

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  long int npart;
  int dim, s;
  double q[DIM], p[DIM];
  double param[3];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim =DIM; 
  q[0] = 5, p[0] = 0;
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  npart = pow(10,5);
  s = 1;
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" NaI case study: No dynamics - estimating number of particles  \n");
  printf("---> NPART     :  %ld\n",npart);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  NaI \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(pow(10,5), dim); 
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "NaI", param);
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
  int count_sa = 0;
  int count_sa1 = 0;
  int count_sa2 = 0;
  int count_sa3 = 0;
  int count_sa12 = 0;
  int count_sa13 = 0;
  int count_sa23 = 0;
  int count_sa123 = 0;
  // the crossing has been pre-computed
  double x_c[1] = {13.27801894097567};
  

  for (int i=0; i< (int) (npart/pow(10,5)); i++){
    sh_wigner_fill(particles, q, p, sqrt(EPS/2), pow(10,5), dim);  
    sh_particle_potential_init(particles, pot, pow(10,5), dim);// initialise particle values - potential, gradient, level ...

    struct Particle *part = particles; // come up with a better structure than a linked list
    for(unsigned int i=0; i<pow(10,5); i++, part++){
      // energy conservation
      part->p_curr[0] = sqrt(pow(part->p[0], 2) + \
          2 * (pot->func_potup(pot, part->x_new, 1) - pot->func_potup(pot, x_c, 1)));
      part->x_curr[0] = x_c[0];
      part->rho_curr = DELTA;
      double pr = ((double)rand() / RAND_MAX); 
      if (sh_transition_lzdia(part, pot, solver)>= pr){ 
        count_lzdia += 1;
      }
      if (sh_transition_lzadia(part, pot, solver)>= pr){ 
        count_lzadia += 1;
      }
      if (sh_transition_sa(part, pot, solver)>= pr){ 
        count_sa += 1;
      }
      if (sh_transition_sa1(part, pot, solver)>= pr){ 
        count_sa1 += 1;
      }
      if (sh_transition_sa2(part, pot, solver)>= pr){ 
        count_sa2 += 1;
      }
      if (sh_transition_sa3(part, pot, solver)>= pr){ 
        count_sa3 += 1;
      }
      if (sh_transition_sa12(part, pot, solver)>= pr){ 
        count_sa12 += 1;
      }
      if (sh_transition_sa13(part, pot, solver)>= pr){ 
        count_sa13 += 1;
      }
      if (sh_transition_sa23(part, pot, solver)>= pr){ 
        count_sa23 += 1;
      }
      if (sh_transition_sa123(part, pot, solver)>= pr){ 
        count_sa123 += 1;
      }
    }
  }
  fprintf(file, "%s \t %.17g \n", "lz_dia", count_lzdia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "lz_adia", count_lzadia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa", count_sa * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa1", count_sa1 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa2", count_sa2 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa3", count_sa3 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa12", count_sa12 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa13", count_sa13 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa23", count_sa23 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \n", "sa123", count_sa123 * 1.0 / npart);
  return 0;
}

