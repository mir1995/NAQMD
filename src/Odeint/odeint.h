#ifndef ODEINT_DOT_H
#define ODEINT_DOT_H 

#include "../SurfaceHopping/surface_hopping.h"

//struct Odeint;
struct Odeint{
  double        t;
  double        dt;
  unsigned int  dim;
  double        abstol;
  double        reltol;
  void (*func_dostep)(struct Odeint *odeint, struct Particle *particle, struct Potential *pot);
};

/* -------------- Generate l-------- */

struct Odeint    *odeint_new(double t, double dt,
                    unsigned int dim, char* solver);

/* ------------ Solvers for hamiltonian systems ------------ */

void odeint_solver_eulsym(struct Odeint *odeint, struct Particle *ptr, struct Potential *pot);
void odeint_solver_lietrotter(struct Odeint *odeint, struct Particle *ptr, struct Potential *pot);

#endif
