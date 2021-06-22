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
    if (argc == 4){
      q = -10;
      dim = 1;
      delta = 0.5;
      alpha = 0.5;
      eps = atof(argv[1]);
      p = atof(argv[2]);
      npart = atoi(argv[3]);
      //
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
      p = 6; // 2 3 4 5 6
      delta = 0.5;
      alpha = 0.5;
      eps = 0.02; // 0.1 0.05 0.03 0.25 0.02 0.01
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

    char filename[200];
    sprintf(filename, "./SurfaceHopping/data/variancedelta%.4falpha%.4fnpart%d.txt",
                      param[1], param[2], npart); // fix transition rate name
    qfile = fopen(filename, "a"); // put some parameters on the file name


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


    printf(" ***************************************************************\n");
    printf("   Surface Particle Hopping Simulation  \n");
    printf(" ***************************************************************\n");

    fprintf(qfile, "eps \t p0 \t std_pos_up \n");
    
    for (int itr=0; itr < (int)((t*1.)/dt); itr++){

      struct Particle *part = particles; // come up with a better structure than a linked list

      while(part != NULL){
        // 2) STEP SOLUTION IN TIME //
        solver->func_dostep(solver, part, pot);
        // 3) UPDATE PARTICLE INFORMATION //
        sh_particle_potential_update(part, pot, dim);
        part = part->next;
      }
      
      sh_observables_update(observables, particles, dim);
      
    }
    // compute variance of random variable (a \circ Phi^t*)
    struct Particle *part = particles;
    double sum2 = 0;
    while (part != NULL){
      sum2 += pow(part->x[0], 2);  
      part = part->next;
    }
    sum2 /= npart;
    fprintf(qfile, "%.17g \t %.17g \t %.17g\n", 
                    eps, p, sqrt(sum2 - pow(observables->x_up[0],2)));
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
