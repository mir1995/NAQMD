/* main.c */
/* 0 copyright/licensing */
/* 1 includes */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include "surface_hopping.h"
#include "../Odeint/odeint.h"
#include "../Potential/potential.h"
#include <time.h> // srand(time(0))
/* 2 defines */
#define GetBit(A,k)         ( A[(k)/16] & (1 << ((k)%16) ) )
#define SetBit(A,k)         ( A[(k)/16] |= (1 << ((k)%16)) )
#define ClearBit(A,k)       ( A[(k)/16] &= ~(1 << ((k)%16)) )
#define TestBit(A,k)        ( A[(k)/16] & (1 << ((k)%16)) )
#define SUB(i,j,d)          ( (i)*(d) + (j) )
//#define Type double
/* 3 external declarations */
/* 4 typedefs */
/* 5 global variable declarations */
/* 6 function prototypes */
// surface hop
//void pos_mean(double** pos, double*  const int n, const int dim)


int main(int argc, char *argv[]){
    // you might want to consider having two .c files that you run using extern?

    /*
     * INITIALISE PARAMETERS
     */
    int npart, dim, t;
    double alpha, delta, eps, dt, q, p;
    double param[3];
    FILE *qfile;

    /*
     * SET PARAMETERS FROM COMMAND LINE
     */
    if (argc == 5){
      delta = atof(argv[1]);
      eps = atof(argv[2]);
      p = atof(argv[3]);
      npart = atoi(argv[4]);
      //
      alpha = 0.5;
      q = -10;
      dim = 1;
      param[0] = eps;
      param[1] = delta;
      param[2] = alpha;
      t = 20/(int)p;
      dt = 0.01/p;
    }
    else {
      npart = 10000;
      dim = 1;
      q = -10;
      p = 2;
      delta = 0.5;
      alpha = 0.5;
      eps = 0.1;
      param[0] = eps;
      param[1] = delta;
      param[2] = alpha;
      t = 20/(int)p;
      dt = 0.01/p;
      printf("Using default parameters. \n");
    }

    /*
     *  PRINT SIMULATION PARAMETERS
    */

    printf(" ***************************************************************\n");
    printf("  Surface Particle Hopping Simulation: parameters  \n");
    printf("---> NPART     :  %d\n",npart);
    printf("---> DIM       :  %d\n",dim);
    printf("---> T         :  %d\n", t);
    printf("---> dt        :  %.4f\n", dt);
    printf("---> q         :  %.4f\n", q);
    printf("---> p         :  %.4f\n", p);
    printf("---> POTENTIAL :  Simple1\n"); 
    printf("---> RATE      :  Goddard\n");
    printf(" ***************************************************************\n");

    /*
     *  GENERATE ARRAY FOR PARTICLES, INITIALISE SOLVER AND POTENTIAL
     */
    struct Particle *particles = sh_particles_create(dim, npart); // I would have just one functions returning the values - "hide what is not needed"
    struct Odeint *solver = odeint_new(t, dt, dim, "lietrotter_symplectic");
    struct Potential *pot = potential_construct(&v_up, &v_down, &grad_v_up, &grad_v_down,
                                                &dd_v_up, &dd_v_down, &get_tau, "simple1", param);
    struct Observables *observables = sh_observables_new(npart, dim);
    struct Hopper *hopper = sh_hopper_new("lasser");

    char filename[200];
    sprintf(filename, "./SurfaceHopping/data/observablesdelta%.4feps%.4falpha%.4fp%dfnpart%d%s.txt",
                      param[1], param[0], param[2], (int)p, npart, "lasser"); // fix transition rate name
    qfile = fopen(filename, "w"); // put some parameters on the file name


    /*
     *  SAMPLE PARTICLES FROM INITIAL WIGNER DISTRIBUTION (GAUSSIAN FOR THE MOMENT)
     */
    printf(" ***************************************************************\n");
    printf("    INITIALISE PARTICLES FROM WIGNER DISTRIBUTION   \n");
    printf("    INITIALISE PARTICLES' DATA: POS, MOM, POT, ...   \n");
    printf(" ***************************************************************\n");

    srand(1); // INITIALISE SEED // what is the actual variance...?
    sh_wigner_fill(particles, q, p, sqrt(eps/2), npart, dim);
    
    sh_particle_potential_init(particles, pot, dim);// initialise particle values - potential, gradient, level ...


    /*
     * TO DO:
     * convergence tests for numerical scheme, wigner sampling, .... ans so on
     */

    /*
     *  SIMULATION: STEP - CHECK FOR AVOIDED CROSSING - UPDATE POTENTIAL - STEP
     *  AT LEAST ONE SIMULATION CHECKING SURFACE HOPPING TO WAVEPACKET DYNAMICS
     */
    
    printf(" ***************************************************************\n");
    printf("   Surface Particle Hopping Simulation  \n");
    printf(" ***************************************************************\n");
    fprintf(qfile, "itr \t pos_up \t pos_down \t \
        mom_up \t mom_down \t ke_up \t ke_down \t \
        e_up \t e_down \t mass_up \t mass_down \n");

    
    for (int itr=0; itr < (int)((t*1.)/dt); itr++){

      struct Particle *part = particles; // come up with a better structure than a linked list

      while(part != NULL){
        // 2) STEP SOLUTION IN TIME //
        solver->func_dostep(solver, part, pot);
        // 3) UPDATE PARTICLE INFORMATION //
        sh_particle_potential_update(part, pot, dim);
        // 4) CALL HOPPER //
        hopper->func_hop(part, hopper, pot, solver);
        part = part->next;
      }
      
      if(itr % (int)(((t*1.)/dt)/40) == 0){
        sh_observables_update(observables, particles, dim);
        fprintf(qfile, "%d \t %.17g \t %.17g \t  %.17g \
                            \t %.17g \t %.17g \t %.17g \
                            \t %.17g \t %.17g \t %.17g \
                            \t %.17g \n",
            itr,
            observables->x_up[0], observables->x_down[0],
            observables->p_up[0], observables->p_down[0],
            observables->ke_up, observables->ke_down,
            observables->e_up, observables->e_down,
            observables->mass_up, observables->mass_down
            );
        // 1) COMPUTE OBSERVABLES AND PRINT TO THEM TO FILE //
        /*
        printf(" ***************************************************************\n");
        printf("       Iteration # %d       \n", itr);
        printf("       Observables          \n");
        printf("mass_up %.4f \n", observables->mass_up);
        printf("mass_down %.4f \n", observables->mass_down);
        printf("position up %.5f \n", observables->x_up[0]);
        printf("position down %.5f \n", observables->x_down[0]);
        printf(" ***************************************************************\n");
        */
      }
    }
    fclose(qfile);
    return 0;
}

// gcc -Wall -c foo.c , gcc -Wall -c main.c  - link them ? --> gcc -o name foo.o

/*
double mean = 0;
double mean2 = 0;
// test variance of wigner distribution
struct Particle *part = particles;
while(part != NULL){
  mean += part->x[0];
  mean2 += pow(part->x[0],2);
  part = part->next;
}
mean /= (double)npart;
mean2 /= (double)npart;
printf("mean is %f\n", mean);
printf("std is %f\n", sqrt(mean2 - pow(mean,2)));
printf("std given is %f\n", sqrt(eps/2));
test variance of wigner distribution 
*/
