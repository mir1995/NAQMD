#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "surface_hopping.h"

double sh_transition_belayev(struct Particle *particle, struct Potential *potential,
                              struct Odeint *odeint){
  double *dd_v = malloc(sizeof(double) * odeint->dim);
  double *qc = malloc(sizeof(double) * odeint->dim); 
  qc[0] = 0; 
  potential->func_dd_up(potential, dd_v, qc, odeint->dim); // not a good idea since already evolved?
  double ddz = 2*dd_v[0]; 
  potential->func_gradup(potential, dd_v, qc, odeint->dim); // not a good idea since already evolved?
  double dz = 2 * dd_v[0]; 
  double dq = (double)particle->p_curr[0]; // should do it more generally
  double ddq = (particle->p[0] - particle->p_curr[0]) / odeint->dt;
  double f = sqrt(pow(2*potential->delta,3) / (pow(dq,2) * ddz + dz * ddq));
  //double f = sqrt(pow(2*fabs(particle->pot_curr[0]),3) / (pow(dq,2) * ddz + dz * ddq));
  free(dd_v);
  free(qc);
  return exp(- M_PI / (2 * potential->eps) * f ); 
}

double sh_transition_landau(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  return 0.5;
}

double sh_transition_lasser(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  double z = 2*potential->delta; // could also consider but it should not matter
  double k = abs(potential->alpha * (1 - pow(tanh(particle->x_curr[0]),2)) * particle->p_curr[0]); // DO IT FOR GENERAL POTENTIAL
  return exp(- M_PI / potential->eps * pow(z,2) / 4 / k);
}

double sh_transition_goddard(struct Particle *particle, struct Potential *potential,
                            struct Odeint *odeint){
  // how would you implemented in general dimension?
  // what wounld you about the momentum of the particle...?
  double tau = potential->func_get_tau(potential);
  // is there a minus if we are transitioning from bottom up?
  double k = sqrt(pow(particle->p_curr[0], 2) + 4*potential->delta);
  return exp(- tau/(potential->delta*potential->eps) * fabs(k - particle->p_curr[0]));
  // this is best explained as the difference between the
  // momentum of classical particle between outgoing and incoming
  // why does more oscillations imply less transition...
}


