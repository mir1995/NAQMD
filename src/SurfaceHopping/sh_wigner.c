#include <stdlib.h>
#include <math.h>
#include "surface_hopping.h"
#include <stdio.h>

void sh_wigner_fill(struct Particle *particles,
                    double *q, double *p, double std, 
                    const unsigned int npart, const unsigned int dim){
    // Box-Muller generates two samples
    // you might want to pass in two samples
    // the loop is a tod odd because of the way 
    // we are setting up the particle class
    //for (int i=0;; i < npart; i++){
    //void *end;
    //for (end = &ptr[npart]; ptr != end; ptr++){
    //linked lists are not ideal for what we want to do
    struct Particle *part = particles;
    for (unsigned int i=0; i<npart; i++, part++){
      for (unsigned int j=0; j < dim; j++){
        // i think jump conditions instructions are used in assembly
        // generate two from U[0,1]
        double x = (double)rand() / RAND_MAX;   
        double y = (double)rand() / RAND_MAX;   
        // https://www.geeksforgeeks.org/rand-and-srand-in-ccpp/
        // generate exp distribution
        // project onto cartesian co-ordinates
        // Box-Muller for independent distributions
        // the mean and variance are rather restrictive at the moment
        // are the two samples correlated though?
        part->x_new[j] = (sqrt(- 2 * log(x)) * cos(2*M_PI*y))*std + q[j]; 
        part->x_curr[j] = part->x_new[j]; 
        part->x_old[j] = part->x_new[j]; 
        part->p[j] = (sqrt(- 2 * log(x)) * sin(2*M_PI*y))*std + p[j]; 
      }
    }
}

// you should have a random library where you generate gaussian samples
