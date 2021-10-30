#include <stdlib.h>
#include "surface_hopping.h"
#include <stdio.h>



struct Particle *sh_particles_create(const unsigned int npart, const unsigned int dim){
  
  // check pre-processing stage has worked -  correct dimension
  //if ()
  struct Particle *particles = malloc(npart * sizeof(struct Particle));    

  return particles;
}


void sh_particle_potential_update(struct Particle *part, struct Potential *pot,
                                  const unsigned int dim){
  part->pot_old = part->pot_curr; // how do I know it will not just change what is pointing to! TO FIGURE OUT
  part->pot_curr = part->pot_new;
  part->rho_old = part->rho_curr; // how do I know it will not just change what is pointing to! TO FIGURE OUT
  part->rho_curr = part->rho_new;
  if (part->state){
    part->pot_new = pot->func_potup(pot, part->x, dim);
    pot->func_gradup(pot, part->pot_grad, part->x, dim);
  }
  else{
    part->pot_new = pot->func_potdown(pot, part->x, dim);
    pot->func_graddown(pot, part->pot_grad, part->x, dim);
  }
  part->rho_new = pot->func_rho(pot, part->x, dim);
}

void sh_particle_potential_init(struct Particle *particles, struct Potential *pot, 
                                const unsigned int npart, const unsigned int dim){
  struct Particle *part = particles;
  for(unsigned int i=0; i<npart; i++, part++){
    // state
    part->state = true;
    // potential
    part->pot_old = pot->func_potup(pot, part->x, dim);
    part->pot_curr = part->pot_old;
    part->pot_new= part->pot_old;
    // gap
    part->rho_old = pot->func_rho(pot, part->x, dim);
    part->rho_curr = part->rho_old;
    part->rho_new= part->rho_old;
    // gradient
    pot->func_gradup(pot, part->pot_grad, part->x, dim);
  }
}

