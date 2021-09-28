#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "surface_hopping.h"
#include "../Potential/potential.h"

struct Hopper *sh_hopper_new(char* transition_name){

  struct Hopper *hopper = malloc(sizeof(struct Hopper));
  if (!strcmp(transition_name, "goddard")){
    hopper->func_transition_probability = sh_transition_goddard;
  }
  else if (!strcmp(transition_name, "goddard1")){
    hopper->func_transition_probability = sh_transition_goddard1;
  }
  else if (!strcmp(transition_name, "goddard2")){
    hopper->func_transition_probability = sh_transition_goddard2;
  }
  else if (!strcmp(transition_name, "goddard3")){
    hopper->func_transition_probability = sh_transition_goddard3;
  }
  else if (!strcmp(transition_name, "lasser")){
    hopper->func_transition_probability = sh_transition_lasser;
  }
  else if (!strcmp(transition_name, "belayev")){
    hopper->func_transition_probability = sh_transition_belayev;
  }
  else if (!strcmp(transition_name, "landau")){
    hopper->func_transition_probability = sh_transition_landau;
  }
  hopper->func_hop = sh_hopper_hop; 

  return hopper;
}

/* Single Switch Surface Hopping */
// I would implement a different hopping function for Tully's algorithm called 
// FSSH
// check but it could also be the case that you do not have to 
// pass all these function parameters? as long as they are 
// included in the header files
void sh_hopper_hop(struct Particle *particle, struct Hopper *hopper, 
                    struct Potential *potential, struct Odeint *odeint){
  
  if ( (particle->rho_new - particle->rho_curr) * (particle->rho_old - particle->rho_curr) < 0){
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

