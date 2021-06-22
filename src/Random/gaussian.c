#include <stdlib.h>
#include <math.h>
#include "random.h"


double rand_gaussian_sample(double mean, double std){
  // Box-Muller generates two samples
  double x = (double)rand() / RAND_MAX;   
  double y = (double)rand() / RAND_MAX;   
  // https://www.geeksforgeeks.org/rand-and-srand-in-ccpp/
  // generate exp distribution
  // project onto cartesian co-ordinates
  // Box-Muller for independent distributions
  //ptr->p[j] = sqrt(- 2 * log(x)) * sin(2*M_PI*y); 
  return (sqrt(- 2 * log(x)) * cos(2*M_PI*y) )*std + mean; 
}

// you should have a random library where you generate gaussian samples
