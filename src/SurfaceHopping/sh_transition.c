#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "surface_hopping.h"

// add a reference to the paper equation number
double sh_transition_lzadia(struct Particle *particle, struct Potential *potential,
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

double sh_transition_lzdia(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double k = sqrt(pow(potential->func_zd(potential, particle->x_curr, 1), 2) + \
             pow(potential->func_v12d(potential, particle->x_curr, 1), 2)) \
             * particle->p_curr[0];
  return exp(- M_PI / potential->eps * pow(potential->delta,2) / k);
}

double sh_transition_sa(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  // there will be a minus sign if transitioning from top to bottom
  // you could use the status of the particle for that
  // there could also be a sign inside of b
  double tau = potential->func_get_tau(potential);
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  double b = pow((k + particle->p_curr[0]) / 2 / particle->p_curr[0], 2);
  return b * exp(- tau/(potential->delta*potential->eps) * fabs(k - particle->p_curr[0]));
}

double sh_transition_sa1(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = M_PI * pow(potential->delta,2) / 2 / potential->alpha; // this is an approximate computation of tau
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  double b = pow((k + particle->p_curr[0]) / 2 / particle->p_curr[0], 2);
  return b * exp(- tau/(potential->delta*potential->eps) * fabs(k - particle->p_curr[0]));
}

double sh_transition_sa2(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = potential->func_get_tau(potential);
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  double b = pow((k + particle->p_curr[0]) / 2 / particle->p_curr[0], 2);
  return b * exp(- tau/(potential->delta*potential->eps) * \
      fabs(2*potential->delta / particle->p_curr[0]));
}

double sh_transition_sa3(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = potential->func_get_tau(potential);
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  return exp(- tau/(potential->delta*potential->eps) * fabs(k - particle->p_curr[0]));
}

double sh_transition_sa12(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = M_PI * pow(potential->delta,2) / 2 / potential->alpha;
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  double b = pow((k + particle->p_curr[0]) / 2 / particle->p_curr[0], 2);
  return b * exp(- tau/(potential->delta*potential->eps) * fabs(2*potential->delta / particle->p_curr[0]));
}

double sh_transition_sa13(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = M_PI * pow(potential->delta,2) / 2 / potential->alpha;
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  return exp(- tau/(potential->delta*potential->eps) * fabs(k - particle->p_curr[0]));
}

double sh_transition_sa23(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = potential->func_get_tau(potential);
  return exp(- tau/(potential->delta*potential->eps) * fabs(2*potential->delta / particle->p_curr[0]));
}

double sh_transition_sa123(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double tau = M_PI * pow(potential->delta,2) / 2 / potential->alpha;
  return exp(- tau/(potential->delta*potential->eps) * fabs(2*potential->delta / particle->p_curr[0]));
}

