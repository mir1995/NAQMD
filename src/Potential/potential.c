#include <string.h>
#include <stdlib.h>
#include "potential.h"
#include <math.h>


/* ------------ Construct potential from model systems -------------*/
double rho(struct Potential *pot, double *x, unsigned int dim){
  return sqrt(pow(v_z(pot, x, dim), 2) + pow(v_v12(pot, x, dim),2));
}

double v_up(struct Potential *pot, double *x, unsigned int dim){
  return rho(pot, x, dim) + v_trace(pot, x, dim);
}

double v_down(struct Potential *pot, double *x, unsigned int dim){
  return - rho(pot, x, dim) + v_trace(pot, x, dim);
}

void grad_v_up(struct Potential *pot, double *grad_v, double *x, unsigned int dim){
  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = (v_z(pot, x, dim) * v_zd(pot, x, dim) + v_v12(pot, x, dim) * \
        v_v12d(pot, x, dim)) / rho(pot, x, dim) + v_traced(pot, x, dim);
  }
}

void grad_v_down(struct Potential *pot, double *grad_v, double *x, unsigned int dim){

  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = -(v_z(pot, x, dim) * v_zd(pot, x, dim) + v_v12(pot, x, dim) * \
        v_v12d(pot, x, dim)) / rho(pot, x, dim) + v_traced(pot, x, dim);
  }
}

struct Potential    *potential_construct(
    // i want to pass in the minimum necessary functions that i need to 
    // construct my potential structure
    //double (*func_rho)(struct Potential *pot, double *x, unsigned int dim),
    //double (*func_potup)(struct Potential *pot, double *x, unsigned int dim), 
    //double (*func_potdown)(struct Potential *pot, double *x, unsigned int dim), 
    //void (*func_gradup)(struct Potential *pot, double *grad_v, double *x, unsigned int dim), 
    //void (*func_graddown)(struct Potential *pot, double *grad_v, double *x, unsigned int dim),
    double (*func_trace)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_z)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_v12)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_traced)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_zd)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_v12d)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_zdd)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_v12dd)(struct Potential *pot, double *x, unsigned int dim),
    void (*func_dd_up)(struct Potential *pot, double *dd_v, double *x, unsigned int dim), 
    void (*func_dd_down)(struct Potential *pot, double *dd_v, double *x, unsigned int dim),
    double (*func_get_tau)(struct Potential *pot),
    char* potential_name, double param[]){
  
  struct Potential *pot = malloc(sizeof(struct Potential));
  // parameters
  pot->delta = param[1];
  pot->alpha = param[2];
  pot->eps = param[0];
  // functions
  pot->func_rho = rho;
  pot->func_potup = v_up;
  pot->func_potdown = v_down;
  pot->func_trace = v_trace;
  pot->func_z = func_z;
  pot->func_v12 = func_v12;
  pot->func_traced = v_traced;
  pot->func_zd = func_zd;
  pot->func_v12d = func_v12d;
  pot->func_zdd = func_zdd;
  pot->func_v12dd = func_v12dd;
  pot->func_gradup = grad_v_up;
  pot->func_graddown = grad_v_down;
  pot->func_dd_up = func_dd_up;
  pot->func_dd_down = func_dd_down;
  pot->func_get_tau = func_get_tau;

  return pot;
}


