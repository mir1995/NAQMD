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
//#define Type double
/* 3 external declarations */
/* 4 typedefs */
/* 5 global variable declarations */
/* 6 function prototypes */

int main(int argc, char *argv[]){
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
    if (argc == 7){
      npart = atoi(argv[1]);
      dim = atoi(argv[2]);
      delta = atof(argv[3]);
      eps = atof(argv[4]);
      q = atof(argv[5]);
      p = atof(argv[6]);
      //
      alpha = 0.5;
      param[0] = eps;
      param[1] = delta;
      param[2] = alpha;
      t = 20/(int)p;
      dt = 0.01/p;
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

    /*
     *  SIMULATION: STEP - CHECK FOR AVOIDED CROSSING - UPDATE POTENTIAL - STEP
     *  AT LEAST ONE SIMULATION CHECKING SURFACE HOPPING TO WAVEPACKET DYNAMICS
     */
    
    printf(" ***************************************************************\n");
    printf("   Surface Particle Hopping Simulation  \n");
    printf(" ***************************************************************\n");
    //fprintf(qfile, "pos_up \t pos_down \t \
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
    }
    sh_observables_update(observables, particles, dim);
    fprintf(qfile, "%.17g \t %.17g \t  %.17g \
                \t %.17g \t %.17g \t %.17g \
                \t %.17g \t %.17g \t %.17g \
                \t %.17g \n",
              observables->x_up[0], observables->x_down[0],
              observables->p_up[0], observables->p_down[0],
              observables->ke_up, observables->ke_down,
              observables->e_up, observables->e_down,
              observables->mass_up, observables->mass_down);

    fclose(qfile);
    return 0;
}

