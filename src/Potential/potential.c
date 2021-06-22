#include <string.h>
#include <stdlib.h>
#include "potential.h"


/* ------------ Construct potential from model systems -------------*/
struct Potential    *potential_construct(
    double (*func_potup)(struct Potential *pot, double *x, unsigned int dim), 
    double (*func_potdown)(struct Potential *pot, double *x, unsigned int dim), 
    void (*func_gradup)(struct Potential *pot, double *grad_v, double *x, unsigned int dim), 
    void (*func_graddown)(struct Potential *pot, double *grad_v, double *x, unsigned int dim),
    void (*func_dd_up)(struct Potential *pot, double *dd_v, double *x, unsigned int dim), 
    void (*func_dd_down)(struct Potential *pot, double *dd_v, double *x, unsigned int dim),
    double (*func_get_tau)(struct Potential *pot),
    char* potential_name, double param[]){
  
  struct Potential *pot = malloc(sizeof(struct Potential));
  pot->func_potup = func_potup;
  pot->func_potdown = func_potdown;
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

