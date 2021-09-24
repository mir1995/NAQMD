#include <string.h>
#include <stdlib.h>
#include "potential.h"


/* ------------ Construct potential from model systems -------------*/
double v_up(struct Potential *pot, double *x, unsigned int dim){
  return rho(pot, x, dim) + d(pot, x, dim);
}

double v_down(struct Potential *pot, double *x, unsigned int dim){
  return - rho(pot, x, dim) + d(pot, x, dim);
}

struct Potential    *potential_construct(
    double (*func_rho)(struct Potential *pot, double *x, unsigned int dim),
    //double (*func_potup)(struct Potential *pot, double *x, unsigned int dim), 
    //double (*func_potdown)(struct Potential *pot, double *x, unsigned int dim), 
    double (*func_zd)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_v12d)(struct Potential *pot, double *x, unsigned int dim),
    void (*func_gradup)(struct Potential *pot, double *grad_v, double *x, unsigned int dim), 
    void (*func_graddown)(struct Potential *pot, double *grad_v, double *x, unsigned int dim),
    void (*func_dd_up)(struct Potential *pot, double *dd_v, double *x, unsigned int dim), 
    void (*func_dd_down)(struct Potential *pot, double *dd_v, double *x, unsigned int dim),
    double (*func_get_tau)(struct Potential *pot),
    char* potential_name, double param[]){
  
  struct Potential *pot = malloc(sizeof(struct Potential));
  pot->func_rho = func_rho;
  pot->func_potup = v_up;
  pot->func_potdown = v_down;
  pot->func_zd = func_zd;
  pot->func_v12d = func_v12d;
  pot->func_gradup = func_gradup;
  pot->func_graddown = func_graddown;
  pot->func_dd_up = func_dd_up;
  pot->func_dd_down = func_dd_down;
  pot->func_get_tau = func_get_tau;
  pot->eps = param[0];
  pot->delta = param[1];
  pot->alpha = param[2];
  //pot->delta = param[0];
  //pot->param = param;

  return pot;
}


