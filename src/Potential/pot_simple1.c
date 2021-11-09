#include <string.h>
#include <math.h>
#include "potential.h"

/* Handles construction of potential in the adiabatic representation */

double v_z(struct Potential *pot, double *x, unsigned int dim){
  return pot->alpha * tanh(x[0]);
}

void v_zd(struct Potential *pot, double *x, double *grad, unsigned int dim){
  grad[0] = pot->alpha * 1/pow(cosh(x[0]), 2);
}

void v_zdd(struct Potential *pot, double *x, double *hess, unsigned int dim){
  hess[0] = - pot->alpha * tanh(x[0]) / pow(cosh(x[0]), 2); 
}

double v_v12(struct Potential *pot, double *x, unsigned int dim){
  return pot->delta;
}

void v_v12d(struct Potential *pot, double *x, double *grad, unsigned int dim){
  grad[0] = 0;
}

void v_v12dd(struct Potential *pot, double *x, double *hess, unsigned int dim){
  hess[0] = 0; 
}

double v_trace(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

void v_traced(struct Potential *pot, double *x, double *grad, unsigned int dim){
  grad[0] = 0;
}

double get_tau(struct Potential *pot){
  // probably want to break the code if this is called
  return 2*sqrt(pow(pot->alpha,2) + pow(pot->delta,2)) * \
    (atan(pot->alpha / pot->delta) + atan(pot->delta / pot->alpha) ) - pot->alpha*M_PI;
}

