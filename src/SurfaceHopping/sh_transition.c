#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "surface_hopping.h"
#include "../Auxiliary/metrics.h"

// add a reference to the paper equation number
double sh_transition_lzadia(struct Particle *part, struct Potential *pot,
                              struct Odeint *odeint){
  
  int dim = (int) odeint->dim;
  double *grad_zd = (double *)malloc(dim * sizeof(double));
  double *grad_v12d = (double *)malloc(dim * sizeof(double));
  double *grad_zdd = (double *)malloc((dim * dim) * sizeof(double));
  double *grad_v12dd = (double *)malloc((dim * dim) * sizeof(double));
  double *grad_v = (double *)malloc(dim * sizeof(double));
  double a=0,b=0,c=0,d=0; 
  
  pot->func_zd(pot, part->x_curr, grad_zd, dim); 
  pot->func_v12d(pot, part->x_curr, grad_v12d, dim);
  pot->func_zdd(pot, part->x_curr, grad_zdd, dim); 
  pot->func_v12dd(pot, part->x_curr, grad_v12dd, dim);
  if (part->state){
  pot->func_gradup(pot, grad_v, part->x_curr, dim);
  }
  else{
  pot->func_graddown(pot, grad_v, part->x_curr, dim);
  }
 
  double v1[dim], v2[dim];
  for (int i=0; i<dim; i++){
    v1[i]=0;
    v2[i]=0;
    for (int j=0; j<dim; j++){
      v1[i] += grad_zdd[j + (i*dim)] * part->p_curr[j]; 
      v2[i] += grad_v12dd[j + (i*dim)] * part->p_curr[j]; 
    }
  }
  for (int i=0; i<dim; i++){
    a+= grad_zd[i] *  part->p_curr[i]; // r1 c1
    b+= grad_v12d[i] *  part->p_curr[i]; // r2 c1
    // check the following but it should be correct 
    // implementation
    c+= v1[i] * part->p_curr[i] - (grad_zd[i] * grad_v[i]); 
    d+= v2[i] * part->p_curr[i]- (grad_v12d[i] * grad_v[i]); 
  }
  double k = sqrt(pow(a, 2) + pow(b, 2) +\
      pot->func_z(pot, part->x_curr, dim) * c + \
      pot->func_v12(pot, part->x_curr, dim) * d);
  
  free(grad_zd); 
  free(grad_v12d); 
  free(grad_zdd); 
  free(grad_v12dd); 
  free(grad_v); 
  // the following is wrong
  return exp(- M_PI / pot->eps * pow(part->rho_curr,2) / k ); 
}


double sh_transition_lzdia(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  unsigned int dim = odeint->dim;
  double *grad_zd = (double *)malloc(dim * sizeof(double));
  double *grad_v12d = (double *)malloc(dim * sizeof(double));
  double a=0,b=0; 
  
  pot->func_zd(pot, part->x_curr, grad_zd, odeint->dim);
  pot->func_v12d(pot, part->x_curr, grad_v12d, odeint->dim);
  
  for (int i=0; i<odeint->dim; i++){
    a+= grad_zd[i] *  part->p_curr[i];
    b+= grad_v12d[i] *  part->p_curr[i];
  }
  double k = sqrt(pow(a, 2) + pow(b, 2));
  
  free(grad_zd); 
  free(grad_v12d); 
  
  return exp(- M_PI / pot->eps * pow(part->rho_curr,2) / k);
}

double sh_transition_sa(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = pot->func_get_tau(pot);
  double p = norm_l2(part->p_curr, odeint->dim); // even in 1d my p will be positive 
  double sign = (part->state == 1) ? (1.0) : (-1.0);
  double k = sqrt(pow(p, 2) + sign * 4*part->rho_curr);
  double b = pow((k + p) / 2 / p, 2); 
  
  return b * exp(- tau/(part->rho_curr*pot->eps) * fabs(k - p));
}

double sh_transition_sa1(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = M_PI * pow(part->rho_curr,2) / 2 / pot->alpha; // this is an approximate computation of tau
  double p = norm_l2(part->p_curr, odeint->dim); 
  double sign = (part->state == 1) ? (1.0) : (-1.0);
  double k = sqrt(pow(p, 2) + sign * 4*part->rho_curr);
  double b = pow((k + p) / 2 / p, 2);

  return b * exp(- tau/(part->rho_curr*pot->eps) * fabs(k - p));
}

double sh_transition_sa2(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = pot->func_get_tau(pot);
  double p = norm_l2(part->p_curr, odeint->dim); 
  double sign = (part->state == 1) ? (1.0) : (-1.0);
  double k = sqrt(pow(p, 2) + sign * 4*part->rho_curr);
  double b = pow((k + p) / 2 / p, 2);

  return b * exp(- tau/(part->rho_curr*pot->eps) * fabs(2*part->rho_curr / p));
}

double sh_transition_sa3(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = pot->func_get_tau(pot);
  double p = norm_l2(part->p_curr, odeint->dim); 
  double sign = (part->state == 1) ? (1.0) : (-1.0);
  double k = sqrt(pow(p, 2) + sign * 4*part->rho_curr);

  return exp(- tau/(part->rho_curr*pot->eps) * fabs(k - p));
}

double sh_transition_sa12(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = M_PI * pow(part->rho_curr,2) / 2 / pot->alpha;
  double p = norm_l2(part->p_curr, odeint->dim); 
  double sign = (part->state == 1) ? (1.0) : (-1.0);
  double k = sqrt(pow(p, 2) + sign * 4*part->rho_curr);
  double b = pow((k + p) / 2 / p, 2);

  return b * exp(- tau/(part->rho_curr*pot->eps) * fabs(2*part->rho_curr / p));
}

double sh_transition_sa13(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = M_PI * pow(part->rho_curr,2) / 2 / pot->alpha;
  double p = norm_l2(part->p_curr, odeint->dim); 
  double sign = (part->state == 1) ? (1.0) : (-1.0);
  double k = sqrt(pow(p, 2) + sign * 4*part->rho_curr);

  return exp(- tau/(part->rho_curr*pot->eps) * fabs(k - p));
}

double sh_transition_sa23(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = pot->func_get_tau(pot);
  double p = norm_l2(part->p_curr, odeint->dim); 
  return exp(- tau/(part->rho_curr*pot->eps) * fabs(2*part->rho_curr / p));
}

double sh_transition_sa123(struct Particle *part, struct Potential *pot,
                            struct Odeint *odeint){
  double tau = M_PI * pow(part->rho_curr,2) / 2 / pot->alpha;
  double p = norm_l2(part->p_curr, odeint->dim); 
  return exp(- tau/(part->rho_curr * pot->eps) * fabs(2*part->rho_curr / p));
}

