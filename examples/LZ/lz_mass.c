/* main.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../setup.h"
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
  double q[DIM], p[DIM];
  double param[3];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim = DIM; 
  q[0] = -10, p[0] = 4;
  s = atoi(argv[3]);
  param[0] = EPS;
  param[2] = ALPHA;
  npart = atoi(argv[1]);// pow(10,10); //pow(10,10);
  param[1] = atof(argv[2]); // delta
  //char *rate = argv[2];
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" NaI case study: No dynamics - estimating number of particles  \n");
  printf("---> NPART     :  %ld\n",npart);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  LZ \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(pow(10,5), dim); // it's a weird data struct - to improve
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "LZ", param);
  struct Odeint *solver = odeint_new(20, 0.001, 1, "lietrotter_symplectic"); 
  
  char filename[200];
  sprintf(filename, "data/new_code_mass_transitioned_seed%d_npart%ld.txt", s, npart); 
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

  
  long int count_lzadia = 0;
  long int count_sa = 0;
  long int count_sa2 = 0;
  long int count_sa3 = 0;
  // the crossing is at x = 0
  double x_c[1] = {0.0};
  for (int i=0; i< (int) (npart/pow(10,5)); i++){
    sh_wigner_fill(particles, q, p, sqrt(EPS/2), pow(10,5), dim);  
    sh_particle_potential_init(particles, pot, pow(10,5), dim);// initialise particle values - potential, gradient, level ...
    
    struct Particle *part = particles; // come up with a better structure than a linked list
    for(unsigned int j=0; j<pow(10,5); j++, part++){
      // energy conservation
      part->p_curr[0] = sqrt(pow(part->p[0], 2) + \
              2 * (pot->func_potup(pot, part->x_new, 1) - pot->func_potup(pot, x_c, 1)));
      part->x_curr[0] = x_c[0];
      part->rho_curr = param[1];
      double pr = (double)rand() / RAND_MAX; 
      if (sh_transition_lzadia(part, pot, solver)>= pr){ 
        count_lzadia += 1;
      }
      if (sh_transition_sa(part, pot, solver)>= pr){ 
        count_sa += 1;
      }
      if (sh_transition_sa2(part, pot, solver)>= pr){ 
        count_sa2 += 1;
      }
      if (sh_transition_sa3(part, pot, solver)>= pr){ 
        count_sa3 += 1;
      }
    }
  }
  fprintf(file, "%s \t %.17g \t %.17g \n", "lz_adia", param[1], count_lzadia * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa", param[1], count_sa * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa2", param[1], count_sa2 * 1.0 / npart);
  fprintf(file, "%s \t %.17g \t %.17g \n", "sa3", param[1], count_sa3 * 1.0 / npart);
  return 0;
}
