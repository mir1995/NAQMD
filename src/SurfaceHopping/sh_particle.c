#include <stdlib.h>
#include "surface_hopping.h"
#include <stdio.h>



struct Particle *sh_particle_new(int dim){ // perhaps you want to do library_file_object_function_name
    // perhaps make this internal
  // i very much doubt this is a good structure...it is like a pointer to a pointer  
  // these are called self-referential structs
  struct Particle *ptr = malloc(sizeof(struct Particle));
  ptr->x            =   (double*)calloc(dim, sizeof(double));
  ptr->p            =   (double*)calloc(dim, sizeof(double));
  ptr->x_curr       =   (double*)calloc(dim, sizeof(double));
  ptr->p_curr       =   (double*)calloc(dim, sizeof(double));
  ptr->pot_grad     =   (double*)calloc(dim, sizeof(double));
  // the last three need not be pointers...
  ptr->next         =   NULL;
  //ptr->p            =   ptr->x + dim;
  return ptr;
}

struct Particle *sh_particle_add(int dim, struct Particle *particles){ // perhaps you want to do library_file_object_function_name
 
  struct Particle *particle = sh_particle_new(dim);
  particle->next = particles;
  
  return particle;
}

struct Particle *sh_particles_create(int dim, int npart){
  
  struct Particle *particles = sh_particle_new(dim);    
  for(int i=1; i<npart; i++){
    particles = sh_particle_add(dim, particles);
  }

  return particles;
}

void sh_particle_destroy(struct Particle *ptr){
        free(ptr->x);
        ptr->x = ptr->p = ptr->pot_grad = NULL;
}

void sh_particle_potential_update(struct Particle *ptr, struct Potential *pot,
                                  unsigned int dim){
  ptr->pot_old = ptr->pot_curr; // how do I know it will not just change what is pointing to! TO FIGURE OUT
  ptr->pot_curr = ptr->pot_new;
  ptr->rho_old = ptr->rho_curr; // how do I know it will not just change what is pointing to! TO FIGURE OUT
  ptr->rho_curr = ptr->rho_new;
  if (ptr->state){
    ptr->pot_new = pot->func_potup(pot, ptr->x, dim);
    pot->func_gradup(pot, ptr->pot_grad, ptr->x, dim);
  }
  else{
    ptr->pot_new = pot->func_potdown(pot, ptr->x, dim);
    pot->func_graddown(pot, ptr->pot_grad, ptr->x, dim);
  }
  ptr->rho_new = pot->func_rho(pot, ptr->x, dim);
}

void sh_particle_potential_init(struct Particle *ptr, struct Potential *pot, 
                                          unsigned int dim){
  struct Particle *part = ptr;
  while(part != NULL){
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
    part = part -> next;
  }
}

