#include <string.h>
#include <math.h>
#include "potential.h"

/*
 * The NAI diabatic potential is given by Engel et Metiu 
 * in https://aip.scitation.org/doi/pdf/10.1063/1.456377
 */

// probably create a file of macros for the units -> import it and use it for the conversion
// if conversion then convert and so on
#define EVTOH 0.036749405469679  
#define ATOB 1.8897261254535  
#define E2 1.0000019559716125 
// (temporary) - the following have been precomputed
#define A1 0.02987726664684903 
#define BETA1 2.159043019538545
#define R0 5.045568754960845
#define A2 101.42835909631405
#define B2 2.998501950492449
#define C2 18.911325265329094
#define LAMBDAP 2.753320477414946
#define LAMBDAM 43.39853919180274
#define RO 0.6593254451707261
#define DELTAE0 0.07625501634958394
#define A12 0.002021217300832345
#define BETA12 0.19408776741560804
#define RX 13.095802049392756

/* potential in the diabatic representation */
static double v11(double *x){
  return A1 * exp(-BETA1 * (x[0] - R0));
}

static double v11d(double *x){
  return - BETA1 * v11(x);
}

static double v22(double *x){
  return (A2 + pow(B2 / x[0], 8)) * exp(-x[0]/RO) - E2 / x[0] \
    - E2 * (LAMBDAP + LAMBDAM) / (2 * pow(x[0], 4)) - C2 / pow(x[0], 6) \
    - 2 * E2 * LAMBDAP * LAMBDAM / pow(x[0], 7) + DELTAE0;
}

static double v22d(double *x){
  return - 8 * pow(B2, 8) / pow(x[0], 9) * exp(- x[0] / RO) + \
    - 1/RO * (A2 + pow(B2 / x[0] , 8)) * exp(-x[0] / RO) + E2 / pow(x[0],2) + \
    2 * E2 * (LAMBDAP + LAMBDAM) / pow(x[0],5) + 6 * C2 / pow(x[0], 7) + \
    14 * E2 * LAMBDAP * LAMBDAM / pow(x[0],8);
}

static double v12(double *x){
  return A12 * exp( - BETA12 * pow(x[0] - RX, 2));
}

static double v12d(double *x){
  return - 2 * BETA12  * (x[0] - RX) * v12(x);
}

/* Handles construction of potential in the adiabatic representation */

// the following functions should go in some common file, so that you 
// do not repeat code
static double z(double *x){
  return 0.5 * (v11(x) - v22(x));
}

static double zd(double *x){
  return 0.5 * (v11d(x) - v22d(x));
}

double d(struct Potential *pot, double *x, unsigned int dim){
  return 0.5 * (v11(x) + v22(x));
}

static double dd(double *x){
  return 0.5 * (v11d(x) + v22d(x));
}

double rho(struct Potential *pot, double *x, unsigned int dim){
  return sqrt(pow(z(x), 2) + pow(v12(x),2));
}
// to solve: can define v11, v12, v22 and then call into other functions? it would be much more general...
/*
double v_up(struct Potential *pot, double *x, unsigned int dim){
  return rho(x) + d(x);
}

double v_down(struct Potential *pot, double *x, unsigned int dim){
  return - rho(x) + d(x);
}
*/
double v_zd(struct Potential *pot, double *x, unsigned int dim){
  return zd(x);
}

double v_v12d(struct Potential *pot, double *x, unsigned int dim){
  return v12d(x);
}

void grad_v_up(struct Potential *pot, double *grad_v, double *x, unsigned int dim){
  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = (z(x) * zd(x) + v12(x) * v12d(x)) / rho(pot, x, dim) + dd(x);
  }
}

void grad_v_down(struct Potential *pot, double *grad_v, double *x, unsigned int dim){

  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = - (z(x) * zd(x) + v12(x) * v12d(x)) / rho(pot, x, dim) + dd(x);
  }
}

double get_tau(struct Potential *pot){
  // the value of tau has been computed separately 
  return 0.0023923147;
}
// the following are still to compute
void dd_v_up(struct Potential *pot, double *dd_v, double *x, unsigned int dim){
  return 0;
}

void dd_v_down(struct Potential *pot, double *dd_v, double *x, unsigned int dim){
  return 0;
}

