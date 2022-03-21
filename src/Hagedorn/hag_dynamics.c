# include <math.h>
# include <stdlib.h>
# include "hagedorn.h"
# include "../Potential/potential.h"
# include "../Odeint/odeint.h"
# include <stdio.h>

static void hag_dynamics_kinetic_step(struct HagedornWaves *params, double dt){

  for (unsigned int i=0; i<params->dim; i++){
    params->q[i] += params->p[i] * dt;
    // for (unsigne int j = 0; j<params->dim; j++)
    params->Q[i] += params->P[i] * dt;
    params->s += pow(params->p[i],2) * dt / 2;
  }

}

static void hag_dynamics_potential_step(struct HagedornWaves *params, struct Potential *pot, double dt){

  double v;
  double *grad = (double*)calloc( params->dim, sizeof(double) );
  double *hess = (double*)calloc( pow(params->dim, 2 ), sizeof(double) );
  
  if (params->state){
    v = pot->func_potup(pot, params->q, params->dim);
    pot->func_gradup(pot, grad, params->q, params->dim);
    pot->func_hessup(pot, hess, params->q, params->dim);
  }
  else{
    v = pot->func_potdown(pot, params->q, params->dim);
    pot->func_graddown(pot, grad, params->q, params->dim);
    pot->func_hessdown(pot, hess, params->q, params->dim);
  }
  
  for (unsigned int i=0; i<params->dim; i++){
    params->p[i] -= grad[i] * dt;
    //printf("grad %.4f\n", grad[i]);
    // for (unsigne int j = 0; j<params->dim; j++)
    params->P[i] -= hess[i] * dt * params->Q[i];
    params->s -= v * dt;
  }

  free(grad);
  free(hess);

}

/******************************************************************************/
void hag_dynamics_do_step(struct HagedornWaves *params, struct Potential *pot, 
                          struct Odeint *odeint){
/******************************************************************************/
/*
  Purpose:

   Evolve Hagedorn parameters over a small time step.

  Discussion:

    This program computes the projection of the wavepacket
    at the crossing using a Gaussian-Hermite quadrature rule
    and writes the set of parameters and Coefficients of the
    projected wavepacket to a file.

    The method was developed by Faou, Gradinaru & Lubich in
    ...

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 March 2022

  Author:

    Michael Redenti
*/

  double dt = odeint->dt; 

  hag_dynamics_kinetic_step(params, dt/2);
  hag_dynamics_potential_step(params, pot, dt);
  hag_dynamics_kinetic_step(params, dt/2);

}

