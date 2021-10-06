#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "surface_hopping.h"

double sh_transition_belayev(struct Particle *particle, struct Potential *potential,
                              struct Odeint *odeint){
  double k = sqrt(pow(particle->p_curr[0],2) * (pow(potential->func_zd(potential, particle->x_curr, 1), 2) + \
             pow(potential->func_v12d(potential, particle->x_curr, 1), 2) ) + pow(particle->p_curr[0], 2) *\
      (potential->func_z(potential, particle->x_curr, 1) * potential->func_zdd(potential, particle->x_curr, 1) + \
      potential->func_v12(potential, particle->x_curr, 1) * potential->func_v12dd(potential, particle->x_curr, 1)) + \
      - (potential->func_z(potential, particle->x_curr, 1) * potential->func_zd(potential, particle->x_curr, 1) + \
      potential->func_v12(potential, particle->x_curr, 1) * potential->func_v12d(potential, particle->x_curr, 1)) * \
	potential->func_traced(potential, particle->x_curr, 1) );
  return exp(- M_PI / potential->eps * pow(potential->delta,2) / k ); 
}

double sh_transition_landau(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  return 0.5;
}

double sh_transition_lasser(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double k = sqrt(pow(potential->func_zd(potential, particle->x_curr, 1), 2) + \
             pow(potential->func_v12d(potential, particle->x_curr, 1), 2)) \
             * particle->p_curr[0];
  return exp(- M_PI / potential->eps * pow(potential->delta,2) / k);
}

double sh_transition_goddard(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = potential->func_get_tau(potential);
  // is there a minus if we are transitioning from bottom up?
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  return exp(- tau/(potential->delta*potential->eps) * fabs(k - particle->p_curr[0]));
}

double sh_transition_goddard1(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = M_PI * pow(potential->delta,2) / 2 / potential->alpha; // this is an approximate computation of tau
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  return exp(- tau/(potential->delta*potential->eps) * fabs(k - particle->p_curr[0]));
}

double sh_transition_goddard2(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = potential->func_get_tau(potential);
  return exp(- tau/(potential->delta*potential->eps) * \
      fabs(2*potential->delta / particle->p_curr[0]));
}

double sh_transition_goddard3(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = M_PI * pow(potential->delta,2) / 2 / potential->alpha;
  return exp(- tau/(potential->delta*potential->eps) * fabs(2*potential->delta / particle->p_curr[0]));
}

