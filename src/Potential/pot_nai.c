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

/* 
 * Potential in the diabatic representation 
 *
 */

static double v11(double *x){
  return A1 * exp(-BETA1 * (x[0] - R0));
}

static double v11d(double *x){
  return - BETA1 * v11(x);
}

static double v11dd(double *x){
  return pow(BETA1,2) * v11(x);
}
// at some point write the following in terms of itself
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

static double v22dd(double *x){
  return 72 * pow(B2, 8) / pow(x[0], 10) * exp(- x[0] / RO) \
    + 1/pow(RO,2) * (A2 + pow(B2 / x[0] , 8)) * exp(-x[0] / RO) + \
    16 * pow(B2,8) / RO / pow(x[0], 9) * exp(- x[0] / RO) \
    - 2 * E2 / pow(x[0],3) + \
    - 10 * E2 * (LAMBDAP + LAMBDAM) / pow(x[0],6) - 42 * C2 / pow(x[0], 8) + \
    - 112 * E2 * LAMBDAP * LAMBDAM / pow(x[0],9);
}

static double v12(double *x){
  return A12 * exp( - BETA12 * pow(x[0] - RX, 2));
}

static double v12d(double *x){
  return - 2 * BETA12  * (x[0] - RX) * v12(x);
}

static double v12dd(double *x){
  return (- 2 * BETA12 + 4 * pow(BETA12, 2) * pow(x[0] - RX, 2) )* v12(x);
}


/* 
 * Handles construction of potential in the diabatic representation 
 *
 * */

double v_z(struct Potential *pot, double *x, unsigned int dim){
  return 0.5 * (v11(x) - v22(x));
}
double v_zd(struct Potential *pot, double *x, unsigned int dim){
  return 0.5 * (v11d(x) - v22d(x));
}
double v_zdd(struct Potential *pot, double *x, unsigned int dim){
  return 0.5 * (v11dd(x) - v22dd(x));
}

double v_v12(struct Potential *pot, double *x, unsigned int dim){
  return v12(x);
}

double v_v12d(struct Potential *pot, double *x, unsigned int dim){
  return v12d(x);
}

double v_v12dd(struct Potential *pot, double *x, unsigned int dim){
  return v12dd(x);
}

double v_trace(struct Potential *pot, double *x, unsigned int dim){
  return 0.5 * (v11(x) + v22(x));
}

double v_traced(struct Potential *pot, double *x, unsigned int dim){
  return 0.5 * (v11d(x) + v22d(x));
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

