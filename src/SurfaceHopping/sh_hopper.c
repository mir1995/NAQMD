#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "surface_hopping.h"
#include "../Auxiliary/metrics.h"
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
  else if (!strcmp(transition_name, "sa3multid")){
    hopper->func_transition_probability = sh_transition_sa3multid;
  }
  hopper->func_hop = sh_hopper_hop; 

  return hopper;
}

/* Single Switch Surface Hopping */
void sh_hopper_hop(struct Particle *part, struct Hopper *hopper, 
                    struct Potential *pot, struct Odeint *odeint){
  
  if ( (part->rho_new - part->rho_curr) * (part->rho_old - part->rho_curr) > 0){
    // exit function if particle's momentum doesn't have enough momentum
    if (!(part->state) && pow(norm_l2(part->p_curr, odeint->dim),2) - 4*part->rho_curr < 0){return;}
    // call transition rate
    double p = hopper->func_transition_probability(part, pot, odeint); 
    // probabilistic SSSH
    if (p >= ((double)rand() / RAND_MAX)){ 
      /* change state of particle */
      part->state = !(part->state);
      // momentum adjustment - re-scaling 
      double k;
      // if hopped from the lower to the higher energy level
      if(part->state){
        part->pot_new = pot->func_potup(pot, part->x_new, odeint->dim);
        pot->func_gradup(pot, part->pot_grad, part->x_new, odeint->dim);
        k = sqrt(1 - 4 * part->rho_curr / pow(norm_l2(part->p, odeint->dim), 2) );
      }
      else{
        part->pot_new = pot->func_potdown(pot, part->x_new, odeint->dim);
        pot->func_graddown(pot, part->pot_grad, part->x_new, odeint->dim);
        k = sqrt(1 + 4 * part->rho_curr / pow(norm_l2(part->p, odeint->dim), 2) );
      }
      for (int i=0; i<odeint->dim; i++){
        part->p[i] *= k;
      }
    }
  }
}

