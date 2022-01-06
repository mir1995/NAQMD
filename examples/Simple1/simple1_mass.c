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


#define ALPHA 0.5 //0.5 //0.05 // might need vary your parameter space

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
  
  dim = DIM; 
  q[0] = -10, p[0] = atof(argv[3]); // 2 4 6 8
  npart = pow(10,9);
  s = atoi(argv[4]);
  param[1] = atof(argv[2]); // DELTA
  param[0] = atof(argv[1]); // EPS
  param[2] = ALPHA;
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" Simple I (model study): Energy conservation  \n");
  printf("---> NPART     :  %ld\n",npart);
  printf("---> DIM       :  %d\n",dim);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  Simple1 \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(pow(10,7), dim); // it's a weird data struct - to improve
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "Simple I", param);
  struct Odeint *solver = odeint_new(20, 0.001, dim, "lietrotter_symplectic");

  char filename[200];
  sprintf(filename, "./data/transition_rate_seed%d_p%d_npart%ld.txt", s, (int) p[0], npart); 
  file = fopen(filename, "a"); 

  if (param[1] == 1){
    fprintf(file, "rate \t eps \t delta \t mass \n");
  }

  /*
   *  SAMPLE PARTICLES FROM ITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
   */
  printf(" ***************************************************************\n");
  printf("    INITIALISE PARTICLES FROM WIGNER DISTRIBUTION   \n");
  printf("    INITIALISE PARTICLES' DATA: POS, MOM, POT, ...   \n");
  printf(" ***************************************************************\n");

  
  srand(s); // INITIALISE SEED // what is the actual variance...?

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
  // the crossing has been pre-computed
  double x_c[1] = {0.0};
  
  for (int i=0; i< (int) (npart/pow(10,7)); i++){
    sh_wigner_fill(particles, q, p, sqrt(param[0]/2), pow(10,7), dim);  
    sh_particle_potential_init(particles, pot, pow(10,7), dim);// initialise particle values - potential, gradient, level ...

    struct Particle *part = particles; // come up with a better structure than a linked list
    for(unsigned int i=0; i<pow(10,7); i++, part++){
      // energy conservation
      part->p_curr[0] = sqrt(pow(part->p[0], 2) + \
              2 * (pot->func_potup(pot, part->x_new, 1) - pot->func_potup(pot, x_c, 1)));
      part->x_curr[0] = x_c[0];
      part->rho_curr = param[1];
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
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "lz_dia", param[0], param[1], count_lzdia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "lz_adia", param[0], param[1], count_lzadia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa", param[0], param[1], count_sa * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa1", param[0], param[1], count_sa1 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa2", param[0], param[1], count_sa2 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa3", param[0], param[1], count_sa3 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa12", param[0], param[1], count_sa12 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa13", param[0], param[1], count_sa13 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa23", param[0], param[1], count_sa23 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \t %.17g \n", "sa123", param[0], param[1], count_sa123 * 1.0 / npart);
  return 0;
}

