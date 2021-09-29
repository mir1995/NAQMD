/* main.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../SurfaceHopping/surface_hopping.h"
#include "../Odeint/odeint.h"
#include "../Potential/potential.h"
#include <time.h> // srand(time(0))


// alpha is used in the approximation to the transition rate
#define ALPHA 0.002669794983146837
#define DELTA 0.002010562925688143 
#define EPS 0.0019227199721312948

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
  q = 5, p = 0;
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  npart = atoi(argv[1]);
  s = atoi(argv[2]);
  char *rate = argv[3];
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" NaI case study: No dynamics - estimating number of particles  \n");
  printf("---> NPART     :  %d\n",npart);
  printf("---> q         :  %.4f\n", q);
  printf("---> p         :  %.4f\n", p);
  printf("---> POTENTIAL :  NaI \n"); 
  printf("---> RATE      :  %s\n", rate);
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  struct Particle *particles = sh_particles_create(dim, pow(10,7)); // it's a weird data struct - to improve
  struct Potential *pot = potential_construct(&rho, &v_zd, &v_v12d, &grad_v_up, &grad_v_down,
                                              &dd_v_up, &dd_v_down, &get_tau, "nai", param);
  struct Odeint *solver = odeint_new(20, 0.001, 1, "lietrotter_symplectic"); // not needed
  struct Hopper *hopper = sh_hopper_new(rate);
  
  char filename[200];
  sprintf(filename, "data/mass_transitioned_seed%d_rate%s.txt", s, rate);
  file = fopen(filename, "a"); 
 
  if (npart == pow(10, 7)){
    fprintf(file, "npart \t mass \n");
  }
  
  /*
   *  SAMPLE PARTICLES FROM INITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
   */
  printf(" ***************************************************************\n");
  printf("    INITIALISE PARTICLES FROM WIGNER DISTRIBUTION   \n");
  printf("    INITIALISE PARTICLES' DATA: POS, MOM, POT, ...   \n");
  printf(" ***************************************************************\n");

  srand(s); // INITIALISE SEED // what is the actual variance...?
  int count = 0; // count number of particles which have transitioned
  double x_c[1] = {13.27801894097567};
    
  for (int i=0; i< (int) npart/pow(10,7); i++){
    sh_wigner_fill(particles, q, p, sqrt(EPS/2), pow(10,7), dim);  
    sh_particle_potential_init(particles, pot, dim);// initialise particle values - potential, gradient, level ...

    /*
     *  SIMULATION: STEP - CHECK FOR AVOIDED CROSSING - UPDATE POTENTIAL - STEP
     *  AT LEAST ONE SIMULATION CHECKING SURFACE HOPPING TO WAVEPACKET DYNAMICS
     */

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
  fprintf(file, "%d \t %.17g \n", npart, count * 1.0 / npart);
  return 0;
}

