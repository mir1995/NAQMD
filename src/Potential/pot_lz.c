#include <string.h>
#include <math.h>
#include "potential.h"

/*
 * The Landau-Zener type potential is given in ?
 */



static double v12(double x, double delta){
  return delta;
}

static double v12d(double x){
  return 0;
}

/* Handles construction of potential in the adiabatic representation */

static double z(double x, double alpha){
  return alpha * x;
}

static double zd(double x, double alpha){
  return alpha;
}

double d(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

static double dd(double *x){
  return 0.0;
}

double rho(struct Potential *pot, double *x, unsigned int dim){
  return sqrt(pow(pot->alpha * x[0],2) + pow(pot->delta,2));
}

double v_zd(struct Potential *pot, double *x, unsigned int dim){
  return pot->alpha;
}

double v_v12d(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

void grad_v_up(struct Potential *pot, double *grad_v, double *x, unsigned int dim){
  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = (z(x[0], pot->alpha) * zd(x[0], pot->alpha) + v12(x[0], pot->delta) * v12d(x[0])) / rho(pot, x, dim) + dd(x);
  }
}

void grad_v_down(struct Potential *pot, double *grad_v, double *x, unsigned int dim){

  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = - (z(x[0], pot->alpha) * zd(x[0], pot->alpha) +\
        v12(x[0], pot->delta) * v12d(x[0])) / rho(pot, x, dim) + dd(x);
  }
}

double get_tau(struct Potential *pot){
  // the value of tau has been computed separately 
  return 0.;
}
// the following are still to compute
void dd_v_up(struct Potential *pot, double *dd_v, double *x, unsigned int dim){
  return 0;
}

void dd_v_down(struct Potential *pot, double *dd_v, double *x, unsigned int dim){
  return 0;
}

