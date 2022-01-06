/* dynamics.c */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "../../src/SurfaceHopping/surface_hopping.h"
#include "../../src/Odeint/odeint.h"
#include "../../src/Potential/potential.h"
#include <time.h> // srand(time(0))
#include "../../setup.h" 

#define GAMMA 0
#define EPS 0.01

int main(int argc, char *argv[]){
    
  /*
   * INITIALISE PARAMETERS
   */
  long int npart;
  int dim, t, s;
  double dt, q[DIM], p[DIM];
  double param[4];
  FILE *file;
  
  /*
   * SET PARAMETERS
   */
  
  dim = DIM; 
  t = 3, dt = 0.001;
  q[0]=5*sqrt(EPS), q[1]=0.5*sqrt(EPS), p[0]=0, p[1]=0;
  npart = pow(10,11); // given that you know the convergence rate of the 
  // two methods, what is the corresponding number of points that 
  // accuracy
  s = 1;
  param[0] = EPS;
  param[3] = GAMMA;
  /*
   *  PRINT SIMULATION PARAMETERS
  */

  printf(" ***************************************************************\n");
  printf(" Jahn Teller: Probabilistic Single Switch Surface Hopping  \n");
  printf("---> NPART     :  %ld\n",npart);
  printf("---> DIM       :  %d\n",dim);
  printf("---> T         :  %d\n", t);
  printf("---> dt        :  %.4f\n", dt);
  printf("---> q         :  %.4f\n", q[0]);
  printf("---> p         :  %.4f\n", p[0]);
  printf("---> POTENTIAL :  Jahn Teller \n"); 
  printf(" ***************************************************************\n");

  /*
   *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
   */
  
  struct Odeint *solver = odeint_new(t, dt, dim, "lietrotter_symplectic");
  // thinking I do not need to pass this as pointer
  struct Potential *pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                                              &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                                              &get_tau, "Jahn Teller", param);

  char filename[200];
  sprintf(filename, 
          "./data/particles_crossing_gamma%d_npart%ld_t%d_dt%.17g.txt", 
          GAMMA, npart, t, dt); 
  file = fopen(filename, "a"); 


  /*
   *  SAMPLE PARTICLES FROM ITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
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
  
  printf(" ***************************************************************\n");
  printf("   Surface Particle Hopping Simulation  \n");
  printf(" ***************************************************************\n");
  
  fprintf(file, "t_c \t x(t_c)[0] \t x(t_c)[1] \t p(t_c)[0] \t p(t_c)[1] \n");

  struct Particle *particle = sh_particles_create(1, dim); 
  
  for (long int i=0; i<npart; i++){
      int itr = 0;

      sh_wigner_fill(particle, q, p, sqrt(EPS/2), 1, dim);

      sh_particle_potential_init(particle, pot, 1, dim);// initialise particle values - potential, gradient, level ...
      
      while (!((particle->rho_new - particle->rho_curr) * (particle->rho_old - particle->rho_curr) > 0)){
          
          solver->func_dostep(solver, particle, pot);
          // 3) UPDATE PARTICLE INFORMATION 
          sh_particle_potential_update(particle, pot, dim);
          itr+=1;
      }
      fprintf(file, "%d \t %.17g \t %.17g \t %.17g \t %.17g \n",
              itr, particle->x_curr[0], particle->x_curr[1],
              particle->p_curr[0], particle->p_curr[1]);
  }
  
  fclose(file);
  return 0;

}
