#include <math.h>
#include "potential.h"


/* I am wondering whether I should
 * include the parameters from an 
 * input file */

# define C1 4.39

# define A1 -0.521
# define A2 0.081
# define A3 0.0

# define a1 1.703
# define a2 1.0
# define a3 1.595

# define C2 -0.45

# define B1 0.698 
# define B2 -0.467
# define B3 0.0

# define b 1.216

/* you could make dim an attribute of potential */

/* Handles construction of potential in the adiabatic representation */

double v_z(struct Potential *pot, double *x, unsigned int dim){
  double B[] = {B1, B2, B3};
  double sum = 0;
  for (unsigned int i=0; i<dim; i++){
    sum += B[i] * x[i];
  }
  return C2 + sum;
}

void v_zd(struct Potential *pot, double *x, double *grad, unsigned int dim){
  double B[] = {B1, B2, B3};
  for (unsigned int i=0; i<dim; i++){
    grad[i] = B[i];
  }
}

void v_zdd(struct Potential *pot, double *x, double *hess, unsigned int dim){
  for (unsigned int i=0; i<dim*dim; i++){
    hess[i] = 0;
  }
}

double v_v12(struct Potential *pot, double *x, unsigned int dim){
  return b * x[2];
}

void v_v12d(struct Potential *pot, double *x, double *grad, unsigned int dim){
  grad[0] = 0; 
  grad[1] = 0;
  grad[2] = b;
}

void v_v12dd(struct Potential *pot, double *x, double *hess, unsigned int dim){
  for (unsigned int i=0; i<dim*dim; i++){
    hess[i] = 0;
  }
}

double v_trace(struct Potential *pot, double *x, unsigned int dim){
  double A[] = {A1, A2, A3};
  double a[] = {a1, a2, a3};
  double sum = 0;
  for (unsigned int i=0; i<dim; i++){
    sum += A[i] * x[i] + 0.5 * pow( a[i]*x[i] , 2);
  }
  return C1 + sum;
}

void v_traced(struct Potential *pot, double *x, double *grad, unsigned int dim){
  double A[] = {A1, A2, A3};
  double a[] = {a1, a2, a3};
  for (unsigned int i=0; i<dim; i++){
    grad[i] = A[i] + pow( a[i], 2 )*x[i];
  }
}

double get_tau(struct Potential *pot){
  // probably want to break the code if this is called
  return 0.0;
}


