/* main.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../SurfaceHopping/surface_hopping.h"
#include "../Odeint/odeint.h"
#include "../Potential/potential.h"
#include <time.h> // srand(time(0))


#define ALPHA 0.5
#define EPS 0.05

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  
  int npart, dim, s;
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
  npart = pow(10,9);
  param[1] = atof(argv[1]);
  char *rate = argv[2];
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" NaI case study: No dynamics - estimating number of particles  \n");
  printf("---> NPART     :  %d\n",npart);
  printf("---> q         :  %.4f\n", q);
  printf("---> p         :  %.4f\n", p);
  printf("---> POTENTIAL :  LZ \n"); 
  printf("---> RATE      :  %s\n", rate);
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(dim, npart); // it's a weird data struct - to improve
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &dd_v_up, &dd_v_down, &get_tau, "LZ", param);
  struct Odeint *solver = odeint_new(20, 0.001, 1, "lietrotter_symplectic"); // not needed
  struct Hopper *hopper = sh_hopper_new(rate);
  
  char filename[200];
  sprintf(filename, "data/mass_transitioned_npart%d_rate%s.txt", npart, rate); 
  file = fopen(filename, "a"); 
  
  if (param[1] == 1){
    fprintf(file, "delta \t mass \n");
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

  
  int count = 0; // count number of particles which have transitioned
  double x_c[1] = {0};
  for (int i=0; i< (int) npart/pow(10,7); i++){
    sh_wigner_fill(particles, q, p, sqrt(EPS/2), npart, dim);  
    sh_particle_potential_init(particles, pot, dim);// initialise particle values - potential, gradient, level ...
    struct Particle *part = particles; // come up with a better structure than a linked list
    while(part != NULL){
      // energy conservation
      part->p_curr[0] = sqrt(pow(part->p[0], 2) + \
          2 * (pot->func_potup(pot, part->x, 1) - pot->func_potup(pot, x_c, 1)));
      part->x_curr[0] = x_c[0];
      if (hopper->func_transition_probability(part, pot, solver)>= ((double)rand() / RAND_MAX)){ // why stochastic and not deterministic
        count += 1;
      }
      part = part->next;
    }
  }
  fprintf(file, "%.17g \t %.17g \n", param[1], count * 1.0 / npart);
  return 0;
}

