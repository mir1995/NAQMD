#include "odeint.h"
#include <stdlib.h>


/* ------- The odeint allows to retrieve scheme parameters
 * ------- while particle for value of potential*/

void odeint_solver_eulsym(struct Odeint *odeint, struct Particle *ptr, struct Potential *pot){
  for (unsigned int j=0; j < odeint->dim; j++){
    ptr->p[j] -= ptr->pot_grad[j] * odeint->dt;
    ptr->x[j] += ptr->p[j]        * odeint->dt;
  }
}

void odeint_solver_lietrotter(struct Odeint *odeint, struct Particle *ptr, struct Potential *pot){
  for (unsigned int j=0; j<odeint->dim; j++){
    ptr->p_curr[j] = ptr->p[j];
    ptr->x_curr[j] = ptr->x[j]; // you'd probably want these two lines more explicitly somewhere else
    ptr->p[j]     -= odeint->dt / 2 * ptr->pot_grad[j];
    ptr->x[j]     += odeint->dt * ptr->p[j];
  }
  double *grad = (double*)calloc(odeint->dim, sizeof(double));
  
  if (ptr->state){
    pot->func_gradup(pot, grad, ptr->x, odeint->dim);
  }
  else{
    pot->func_graddown(pot, grad, ptr->x, odeint->dim);
  }
  for (unsigned int j=0; j<odeint->dim; j++){
    ptr->p[j] -= odeint->dt / 2 * grad[j];
  }
  free(grad);
}
