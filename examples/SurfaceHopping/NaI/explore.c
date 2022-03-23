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
#define EPS 0.0019227199721312948 // wrong - recompute

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
  npart = atoi(argv[1]);
  s = 1;
  param[0] = EPS;
  param[1] = DELTA;
  param[2] = ALPHA;
  char rate[] = "goddard";
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
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &dd_v_up, &dd_v_down, &get_tau,
                                              "NaI", param);
  struct Odeint *solver = odeint_new(20, 0.001, 1, "lietrotter_symplectic"); // not needed
  struct Hopper *hopper = sh_hopper_new(rate);
  
  /*
   * Print diabatic potential elements and see what differs - also compute alpha to see if it agrees 
   *
   */
  double x[1] = {13.27801894097567};
  double p0 = 0.24843052119644576;
  double rho = sqrt(pow(v_z(pot, x, dim), 2) + pow(v_v12(pot, x, dim), 2)); 
  double rhodd = (pow(v_zd(pot, x, 1), 2) + v_z(pot, x, 1)*v_zdd(pot, x, 1) + \
      pow(v_v12d(pot, x, dim), 2) + v_v12(pot, x, dim)*v_v12dd(pot, x, 1)) / rho - \
                 pow(v_z(pot, x, 1)*v_zd(pot, x, 1) + v_v12(pot, x, dim)*v_v12d(pot, x, dim), 2) / pow(rho, 3); 
  printf("v_vz = %.17g \n", v_z(pot, x, 1));
  printf("v_vzd = %.17g \n", v_zd(pot, x, 1));
  printf("v_vzdd = %.17g \n", v_zdd(pot, x, 1));
  printf("v_v12 = %.17g \n", v_v12(pot, x, 1));
  printf("v_v12d = %.17g \n", v_v12d(pot, x, 1));
  printf("v_v12dd = %.17g \n", v_v12dd(pot, x, 1));
  printf("v_vtrace = %.17g \n", v_trace(pot, x, 1));
  // alpha
  printf("alpha = %.17g \n", sqrt(DELTA * rhodd));
  // k in Belayev rate
  double k = sqrt(pow(p0,2) * (pow(pot->func_zd(pot, x, 1), 2) + \
        pow(pot->func_v12d(pot, x, 1), 2) )); 
  k =  pow(p0, 2) \
       * (pot->func_z(pot, x, 1) * pot->func_zdd(pot, x, 1) +\
         pot->func_v12(pot, x, 1) * pot->func_v12dd(pot, x, 1));
  printf("Belayev rate k = %.17g \n", k);

    
}

