// should I assume the surfaces are analytically given 
// o.w. need to evaluate gradient numerically ?
#include <math.h>
#include "potential.h"



/* Handles construction of potential in the adiabatic representation */

double v_z(struct Potential *pot, double *x, unsigned int dim){
  return x[0];
}

void  v_zd(struct Potential *pot, double *x, double *grad, unsigned int dim){
  grad[0] = 1;
  grad[1] = 0;
}

void v_zdd(struct Potential *pot, double *x, double *hess, unsigned int dim){
  hess[0] = 0; // use some MACRO?
  hess[1] = 0; // use some MACRO?
  hess[2] = 0; // use some MACRO?
  hess[3] = 0; // use some MACRO?
}

double v_v12(struct Potential *pot, double *x, unsigned int dim){
  return x[1];
}

void v_v12d(struct Potential *pot, double *x, double *grad, unsigned int dim){
  grad[0] = 0; // i wonder if incrementing the pointer is faster - try on some large example
  grad[1] = 1;
}

void v_v12dd(struct Potential *pot, double *x, double *hess, unsigned int dim){
  hess[0] = 0; // use some MACRO?
  hess[1] = 0; // use some MACRO?
  hess[2] = 0; // use some MACRO?
  hess[3] = 0; // use some MACRO?
}

double v_trace(struct Potential *pot, double *x, unsigned int dim){
  return pot->gamma * (pow(x[0], 2) + pow(x[1], 2));
}

void v_traced(struct Potential *pot, double *x, double *grad, unsigned int dim){
  grad[0] = 2*pot->gamma*x[0];
  grad[1] = 2*pot->gamma*x[1];
}

double get_tau(struct Potential *pot){
  // probably want to break the code if this is called
  return 0.0;
}


