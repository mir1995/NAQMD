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
  int dim, s; 
  double param[3];
  FILE *infl, *oufl;

  dim =DIM; 
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  s=1;
  // read each line 
  // set particle curr pos and mom and gap curr state = 1 
  char filename[200];
  sprintf(filename, "data/dynamics_mass.txt");
  oufl = fopen(filename, "a"); 
 
  fprintf(oufl, "rate \t mass \n");
  
  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *part = sh_particles_create(1, dim); 
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "NaI", param);
  sh_particle_potential_init(part, pot, 1, dim);// initialise particle values - potential, gradient, level ...
  struct Odeint *solver = odeint_new(20, 0.001, 1, "lietrotter_symplectic"); // not needed
  

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

  /*
   *
   * READ INPUT FILE 
   */
  infl = fopen("data/particles_crossing_npart100000000000t80dt0.01.txt", "r");
  if (infl==NULL){
      perror("Unable to open file");
      exit(1);
  }
  char line[150];
  double itr, x, p;
  long int npart = -1;
  while(fgets(line, sizeof(line), infl) != NULL){
    npart += 1;
    if (npart!=0){
        // scan individual data
        sscanf(line, "%lf %lf %lf", &itr, &x, &p );
        // update particle position, momentum and gap
        part->p_curr[0] = p;
        part->x_curr[0] = x;
        part->rho_curr = DELTA; // evaluate the gap at x
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
      // energy conservation

  fprintf(oufl, "%s \t %.17g \n", "lz_dia", count_lzdia * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "lz_adia", count_lzadia * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa", count_sa * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa1", count_sa1 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa2", count_sa2 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa3", count_sa3 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa12", count_sa12 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa13", count_sa13 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa23", count_sa23 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa123", count_sa123 * 1.0 / npart);
#ifdef BLA
#endif 
  return 0;
}

