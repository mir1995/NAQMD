#include <string.h>
#include <math.h>
#include "potential.h"

/*
 * The Landau-Zener type potential is given in ?
 */

/* Handles construction of potential in the adiabatic representation */

double v_z(struct Potential *pot, double *x, unsigned int dim){
  return pot->alpha * x[0];
}

double v_zd(struct Potential *pot, double *x, unsigned int dim){
  return pot->alpha;
}

double v_zdd(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

double v_v12(struct Potential *pot, double *x, unsigned int dim){
  return pot->delta;
}

double v_v12d(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

double v_v12dd(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

double v_trace(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

double v_traced(struct Potential *pot, double *x, unsigned int dim){
  return 0.0;
}

double get_tau(struct Potential *pot){
  return M_PI * pow(pot->delta,2) / 2 / pot->alpha;
}

