/* many.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../setup.h"
#include "../../src/Odeint/odeint.h"
#include "../../src/Potential/potential.h"
#include "../../src/Auxiliary/metrics.h"


#define GAMMA 0
#define EPS 0.01

double temp_transition_samultid(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  int dim = (int) odeint->dim;
  
  double p = norm_l2(part->p_curr, odeint->dim); 
  double sign = (part->state == 1) ? (1.0) : (-1.0);
  double k = sqrt(pow(p, 2) + sign * 4*part->rho_curr);
  
  return exp(- M_PI / 2 / pot->eps * part->rho_curr * fabs(k-p) );
}

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  int dim, s; 
  double param[4];
  FILE *infl, *oufl;
  

  dim =DIM; 
  s=1;
  param[0] = EPS;
  param[3] = GAMMA;
  // read each line 
  // set particle curr pos and mom and gap curr state = 1 
  char filename[200];
  sprintf(filename, "data/jahn_teller_dynamics_mass.txt");
  oufl = fopen(filename, "a"); 
 
  fprintf(oufl, "rate \t mass \n");
  
  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *part = sh_particles_create(1, dim); 
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "JahnTeller", param);
  sh_particle_potential_init(part, pot, 1, dim);// initialise particle values - potential, gradient, level ...
  struct Odeint *solver = odeint_new(20, 0.001, DIM, "lietrotter_symplectic"); // not needed
  

  srand(s); // INITIALISE SEED // what is the actual variance...?
  long int count_lzdia = 0; // count number of particles which have transitioned
  long int count_lzadia = 0;
  long int count_sa1 = 0;
  long int count_sa12 = 0;
  long int count_sa13 = 0;

  /*
   *
   * READ INPUT FILE 
   */
  infl = fopen("data/particles_crossing.txt", "r");
  if (infl==NULL){
      perror("Unable to open file");
      exit(1);
  }
  char line[150];
  double itr, x[DIM], p[DIM];
  long int npart = -1;
  printf("npart %ld\n",npart);
  while(fgets(line, sizeof(line), infl) != NULL){
    npart += 1;
    if (npart!=0){
        // scan individual data
        sscanf(line, "%lf %lf %lf %lf %lf", &itr, &x[0], &(x[1]), &(p[0]), &(p[1]));
        // update particle position, momentum and gap
        part->p_curr[0] = p[0];
        part->p_curr[1] = p[1];
        part->x_curr[0] = x[0];
        part->x_curr[1] = x[1];
        //printf("%.4f\n", part->x_curr[0]);
        part->rho_curr = sqrt(pow(x[0],2) + pow(x[1],2)); // evaluate the gap at x
        double pr = ((double)rand() / RAND_MAX); 
        if (sh_transition_lzdia(part, pot, solver)>= pr){ 
          count_lzdia += 1;
        }
        if (sh_transition_lzadia(part, pot, solver)>= pr){ 
          count_lzadia += 1;
        }
        //if (sh_transition_sa1(part, pot, solver)>= pr){ 
        //count_sa1 += 1;
        //}
        //if (sh_transition_sa12(part, pot, solver)>= pr){ 
        //count_sa12 += 1;
        //}
        if (temp_transition_samultid(part, pot, solver)>= pr){ 
          count_sa13 += 1;
        }
      }
    }
      // energy conservation

  printf(" %ld \n",  npart);
  fprintf(oufl, "%s \t %.17g \n", "lz_dia", count_lzdia * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "lz_adia", count_lzadia * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa1", count_sa1 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa12", count_sa12 * 1.0 / npart);
  fprintf(oufl, "%s \t %.17g \n", "sa13", count_sa13 * 1.0 / npart);
  return 0;
}

