#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "surface_hopping.h"
#include "../Potential/potential.h"

struct Hopper *sh_hopper_new(char* transition_name){

  struct Hopper *hopper = malloc(sizeof(struct Hopper));
  if (!strcmp(transition_name, "sa")){
    hopper->func_transition_probability = sh_transition_sa;
  }
  else if (!strcmp(transition_name, "sa1")){
    hopper->func_transition_probability = sh_transition_sa1;
  }
  else if (!strcmp(transition_name, "sa2")){
    hopper->func_transition_probability = sh_transition_sa2;
  }
  else if (!strcmp(transition_name, "sa3")){
    hopper->func_transition_probability = sh_transition_sa3;
  }
  else if (!strcmp(transition_name, "sa12")){
    hopper->func_transition_probability = sh_transition_sa12;
  }
  else if (!strcmp(transition_name, "sa13")){
    hopper->func_transition_probability = sh_transition_sa13;
  }
  else if (!strcmp(transition_name, "sa23")){
    hopper->func_transition_probability = sh_transition_sa23;
  }
  else if (!strcmp(transition_name, "sa123")){
    hopper->func_transition_probability = sh_transition_sa123;
  }
  else if (!strcmp(transition_name, "lzadia")){
    hopper->func_transition_probability = sh_transition_lzadia;
  }
  else if (!strcmp(transition_name, "lzdia")){
    hopper->func_transition_probability = sh_transition_lzdia;
  }
  hopper->func_hop = sh_hopper_hop; 

  return hopper;
}

/* Single Switch Surface Hopping */
void sh_hopper_hop(struct Particle *particle, struct Hopper *hopper, 
                    struct Potential *potential, struct Odeint *odeint){
  
  if ( (particle->rho_new - particle->rho_curr) * (particle->rho_old - particle->rho_curr) > 0){
    double p = hopper->func_transition_probability(particle, potential, odeint); // particle location and transition rate should suffice
    if (p >= ((double)rand() / RAND_MAX)){ // why stochastic and not deterministic
      /* change state of particle */
      particle->state = !(particle->state);
      /* update potential value new - do I have to for rho as well?? */
      if(particle->state){
        particle->pot_new = potential->func_potup(potential, particle->x, 1);
        potential->func_gradup(potential, particle->pot_grad, particle->x, 1);
      }
      else{
        particle->pot_new = potential->func_potdown(potential, particle->x, 1);
        potential->func_graddown(potential, particle->pot_grad, particle->x, 1);
        particle->p[0] = sqrt(pow(particle->p_curr[0],2) + 4*potential->delta); // to check
      }
    }
  }
}

