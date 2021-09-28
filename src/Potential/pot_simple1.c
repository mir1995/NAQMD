#include <string.h>
#include <math.h>
#include "potential.h"
// in this case, should I have a header file for pot_simple1.h?? it should be a reason for it
// it makes sense not to have a struct Particle as a parameter as you might want to evaluate the 
// potential on a grid and do not want to create particles in order to do this
// Have instead a middle function to take care of this
//
// THIS SIMPLE1 IS A ONE DIMENSIONAL EXAMPLE - CAN YOU EXTEND IT TO MORE THAN ONE DIMENSION?
/*
double v_up(struct Potential *pot, double *x, unsigned int dim){
  return sqrt(pow(pot->alpha * tanh(x[0]), 2) + pow(pot->delta,2));
}

double v_down(struct Potential *pot, double *x, unsigned int dim){
  return  - v_up(pot, x, dim);
}
*/
double rho(struct Potential *pot, double *x, unsigned int dim){
  return sqrt(pow(pot->alpha * tanh(x[0]), 2) + pow(pot->delta,2));
}

double d(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

double get_tau(struct Potential *pot){
  return 2*sqrt(pow(pot->alpha,2) + pow(pot->delta,2)) * \
    (atan(pot->alpha / pot->delta) + atan(pot->delta / pot->alpha) ) - pot->alpha*M_PI;
}

void grad_v_up(struct Potential *pot, double *grad_v, double *x, unsigned int dim){
  //double y[dim]; // i think the problem here is that you have created 
  // a variable whose scope is not global
  // "y is created with automatic storage duration and references to it will 
  // become invalid once it leaves its declaring scope, i.e. when the function 
  // returns"
  for(unsigned int i=0; i<dim; i++){
    double s = tanh(x[i]);
    grad_v[i] = pow(pot->alpha,2) * s * (1 - pow(s,2)) /v_up(pot, x, dim);
  }
}

void grad_v_down(struct Potential *pot, double *grad_v, double *x, unsigned int dim){

  for(unsigned int i=0; i<dim; i++){
    double s = tanh(x[i]);
    grad_v[i] = - pow(pot->alpha,2) * s * (1 - pow(s,2)) /v_up(pot, x, dim);
  }
}

void dd_v_up(struct Potential *pot, double *dd_v, double *x, unsigned int dim){
  for(unsigned int i=0; i<dim; i++){
    double s = tanh(x[i]);
    dd_v[i] = - pow(pot->alpha,4) * pow(s,2) * pow(1 - pow(s,2),2) / pow(v_up(pot, x, dim), 3) \
              + pow(pot->alpha,2) * pow(1 - pow(s,2), 2) / v_up(pot, x, dim) \
              - 2 * pow(pot->alpha,2) * pow(s,2) * (1 - pow(s, 2)) / v_up(pot, x, dim);
  }
}

void dd_v_down(struct Potential *pot, double *dd_v, double *x, unsigned int dim){
  for(unsigned int i=0; i<dim; i++){
    double s = tanh(x[i]);
    dd_v[i] =  pow(pot->alpha,4) * pow(s,2) * pow(1 - pow(s,2),2) / pow(v_up(pot, x, dim), 3) \
              - pow(pot->alpha,2) * pow(1 - pow(s,2), 2) / v_up(pot, x, dim) \
              + 2 * pow(pot->alpha,2) * pow(s,2) * (1 - pow(s, 2)) / v_up(pot, x, dim);
  }
}

/* --------------- Second approach you might want to try ----------------- */
/*
double *v_up(struct Potential *pot, double *x, unsigned int dim, char* level){
  double y[dim]; 
  if(!strcmp(level, "up")){
    y[0] = sqrt(pow(pot->alpha * tanh(x[0]), 2) + pow(pot->delta,2));
  }
  return y 
}
*/
