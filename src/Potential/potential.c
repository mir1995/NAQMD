#include <string.h>
#include <stdlib.h>
#include <stdio.h>
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

// i am thinking this will not work
void grad_v_up(struct Potential *pot, double *grad_v, double *x, unsigned int dim){

  double *grad_zd = (double *)malloc(dim * sizeof(double));
  double *grad_v12d = (double *)malloc(dim * sizeof(double));
  double *grad_traced = (double *)malloc(dim * sizeof(double));
  
  v_zd(pot, x, grad_zd, dim); 
  v_v12d(pot, x, grad_v12d, dim);
  v_traced(pot, x, grad_traced, dim);

  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = (v_z(pot, x, dim) * grad_zd[i] + v_v12(pot, x, dim) * \
        grad_v12d[i]) / rho(pot, x, dim) + grad_traced[i];
  }
  free(grad_zd); 
  free(grad_v12d); 
  free(grad_traced); 
}

void grad_v_down(struct Potential *pot, double *grad_v, double *x, unsigned int dim){

  double *grad_zd = (double *)malloc(dim * sizeof(double));
  double *grad_v12d = (double *)malloc(dim * sizeof(double));
  double *grad_traced = (double *)malloc(dim * sizeof(double));
  
  v_zd(pot, x, grad_zd, dim); 
  v_v12d(pot, x, grad_v12d, dim);
  v_traced(pot, x, grad_traced, dim);
  
  for(unsigned int i=0; i<dim; i++){
    grad_v[i] = -(v_z(pot, x, dim) * grad_zd[i] + v_v12(pot, x, dim) * \
        grad_v12d[i]) / rho(pot, x, dim) + grad_traced[i];
  }
  free(grad_zd); 
  free(grad_v12d); 
  free(grad_traced); 
}

struct Potential    *potential_construct(
    double (*func_trace)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_z)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_v12)(struct Potential *pot, double *x, unsigned int dim),
    void (*func_traced)(struct Potential *pot, double *x, double *grad, unsigned int dim),
    void (*func_zd)(struct Potential *pot, double *x, double *grad, unsigned int dim),
    void (*func_v12d)(struct Potential *pot, double *x, double *grad, unsigned int dim),
    void (*func_zdd)(struct Potential *pot, double *x, double *hess, unsigned int dim),
    void (*func_v12dd)(struct Potential *pot, double *x, double *hess, unsigned int dim),
    double complex (*func_get_tau)(struct Potential *pot),
    char* potential_name, double param[]){
  
  struct Potential *pot = malloc(sizeof(struct Potential));
  // parameters
  pot->delta = param[1];
  pot->alpha = param[2];
  pot->eps = param[0];
  if (param[3]){pot->gamma = param[3];}
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
  pot->func_get_tau = func_get_tau;

  return pot;
}


