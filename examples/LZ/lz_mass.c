/* main.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../src/Odeint/odeint.h"
#include "../../src/Potential/potential.h"
#include <time.h> // srand(time(0))


#define ALPHA 0.5
#define EPS 0.05

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  long int npart;
  int dim, s;
  double q, p;
  double param[3];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim =1; 
  q = -10, p = 4;
  s = 1;
  param[0] = EPS;
  param[2] = ALPHA;
  npart = pow(10,11);
  param[1] = atof(argv[1]); // delta
  //char *rate = argv[2];
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" NaI case study: No dynamics - estimating number of particles  \n");
  printf("---> NPART     :  %ld\n",npart);
  printf("---> q         :  %.4f\n", q);
  printf("---> p         :  %.4f\n", p);
  printf("---> POTENTIAL :  LZ \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(dim, pow(10,7)); // it's a weird data struct - to improve
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &dd_v_up, &dd_v_down, &get_tau, "LZ", param);
  struct Odeint *solver = odeint_new(20, 0.001, 1, "lietrotter_symplectic"); 
  
  char filename[200];
  sprintf(filename, "data/mass_transitioned_npart%ld.txt", npart); 
  file = fopen(filename, "a"); 
  
  if (param[1] == 1){
    fprintf(file, "rate \t delta \t mass \n");
  }
  /*
   *  SAMPLE PARTICLES FROM INITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
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

  
  long int count_lzdia = 0; // count number of particles which have transitioned
  long int count_lzadia = 0;
  long int count_sa = 0;
  long int count_sa1 = 0;
  long int count_sa2 = 0;
  long int count_sa3 = 0;
  long int count_sa12 = 0;
  long int count_sa13 = 0;
  long int count_sa23 = 0;
  long int count_sa123 = 0;
  // the crossing is at x = 0
  double x_c[1] = {0};
  for (int i=0; i< (int) (npart/pow(10,7)); i++){
    sh_wigner_fill(particles, q, p, sqrt(EPS/2), pow(10,7), dim);  
    sh_particle_potential_init(particles, pot, dim);// initialise particle values - potential, gradient, level ...
    struct Particle *part = particles; // come up with a better structure than a linked list
    while(part != NULL){
      // energy conservation
      part->p_curr[0] = sqrt(pow(part->p[0], 2) + \
          2 * (pot->func_potup(pot, part->x, 1) - pot->func_potup(pot, x_c, 1)));
      part->x_curr[0] = x_c[0];
      double pr = (double)rand() / RAND_MAX; 
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
      part = part->next;
    }
  }
  fprintf(file, "%s \t %.17g \t %.17g \n", "lz_dia", param[1], count_lzdia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "lz_adia", param[1], count_lzadia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa", param[1], count_sa * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa1", param[1], count_sa1 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa2", param[1], count_sa2 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa3", param[1], count_sa3 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa12", param[1], count_sa12 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa13", param[1], count_sa13 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa23", param[1], count_sa23 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa123", param[1], count_sa123 * 1.0 / npart);
  return 0;
}
