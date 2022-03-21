#include <string.h>
#include <stdlib.h>
#include "odeint.h"

struct Odeint *odeint_new(double t, double dt, 
                          unsigned int dim, char* solver){
  struct Odeint *odeint = malloc(sizeof(struct Odeint));
  odeint->t = t;
  odeint->dt = dt;
  odeint->dim = dim;
  //if (!strcmp(solver, "euler_symplectic")){
  //  odeint->func_dostep = odeint_solver_eulsym;
  //}
  //else if (!strcmp(solver, "lietrotter_symplectic")){
  //  odeint->func_dostep = odeint_solver_lietrotter;
  //}
  return odeint;
}
