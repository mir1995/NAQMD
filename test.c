#include "surface_hopping.h"
#include <stdlib.h>
#include <stdio.h>

int main(){
  int npart = 2; 
  struct Particle *particles = (struct Particle *)malloc(npart * sizeof(struct Particle));
  
  printf("%p\n", particles);
  for (int i=0; i<npart; i++){
    printf("%f ", particles[i].pot_old);
    printf("0 %p 1 %p \n", &particles[i].pot_old, &particles[i].pot_curr);
  }

  printf("%p\n", particles++);
  printf("%lu\n", sizeof(particles));
  // new idea -> create one particle, assigning memory, and 
  // then get the hopefully new size 
  // to create an array of those sized struct
  printf("%lu\n", sizeof(struct Particle));
  free(particles);
  printf("%lu\n", sizeof(double));
  printf("%lu\n", sizeof(unsigned int));
}

