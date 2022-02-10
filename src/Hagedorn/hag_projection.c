# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include "hag_projection.h"
# include "hagedorn.h"
# include "hag_superadiabatic.h"
# include "../Quadrature/hermite_rule.h"
# include "../Potential/potential.h"


static void evaluate_f(struct Potential *pot,
                                struct HagedornWaves *params_in, 
                                struct HagedornWaves *params_out,
                                unsigned int index_in, unsigned int index_out, 
                                double *x, double complex *f, unsigned int n)
{
  double complex h; // clarify the terms
  
  // can think about splitting it into easier tasks
  hag_superadiabatic_formula(x, f, params_in, n, pot, "constant");
  
  for (unsigned int i=0; i<n; i++){
      h = hag_wavepackets_polynomial_evaluate((double complex) x[i], params_out, index_out); // evaluate Hagedorn polynomial at index_out
      f[i] = f[i] * h;
  }
}

// write this function in a linear algebra programme
static double complex inner_product(double *w, double complex *f, unsigned int n)
{
  double complex sum=0;
  for (unsigned int i=0; i<n; i++){
    sum += ((double complex) w[i]) * f[i];
  }
  return sum;
}


/******************************************************************************/

struct HagedornWaves *hag_project(struct HagedornWaves *params_in, unsigned int n, 
                                  unsigned int K, struct Potential *pot)

/******************************************************************************/
/*
  Purpose:

   Projection of superadiabatic transmitted wavepacket onto Hagedorn basis. 

  Discussion:

    This program computes the projection of the wavepacket 
    at the crossing using a Gaussian-Hermite quadrature rule.

    The input parameters:
    * N, the number of points in the rule
    * A, the center point;

    If SCALE is nonzero, then the factor C is sqrt ( B ) / sqrt ( PI ).
    which means that the function f(x)=1 will integrate to 1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 February 2021

  Author:

    Michael Redenti
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int kind;
  int scale;
  double *w;
  double *x;
  double complex *f;
  

  /******************************************************************************/
  // compute the new parameter set ?
  /******************************************************************************/
  // can you create a copy instead? using pointers perhaps?
  struct HagedornWaves *params_out = malloc(sizeof(struct HagedornWaves)); 
  params_out->state = 0;
  params_out->size = K;
  params_out->eps = params_in->eps;
  params_out->s = params_in->s;
  params_out->q[0] = sqrt(pow(params_in->q[0],2) + 4*pot->delta); // try and compute it exactly 
  params_out->p[0] = params_in->p[0];// this is position
  params_out->Q[0] = params_in->Q[0];
  params_out->P[0] = params_in->P[0];
  params_out->c  = ( double complex * ) malloc ( n * sizeof ( double complex ) ); // these are pointers so easy to pass

  /******************************************************************************/
  // Compute weights and nodes of Gauss-Hermite rule 
  /******************************************************************************/
  w = ( double * ) malloc ( n * sizeof ( double ) ); // your main function could pass these values to a function call
  x = ( double * ) malloc ( n * sizeof ( double ) ); // these are pointers so easy to pass
  f = ( double complex * ) malloc ( n * sizeof ( double complex ) ); // these are pointers so easy to pass
  
  a = creal(params_out->q[0]);
  b = 1.0 / 4 / creal(params_in->eps);
  alpha = 0.0;
  beta = 0.0;
  scale = 0;
  kind = 6; // Gauss Hermite Rule - a name would have been better 
  // compute set of coefficients until what ?? stop changing by some TOL?  
  // for the moment project onto the first 20 wavepackets...?
  

  cgqf ( n, kind, alpha, beta, a, b, x, w );

/*
  Normalize the rule.
  This way, the weights add up to 1.
*/
  if ( scale == 1 )
  {
    for (unsigned int i = 0; i < n; i++ )
    {
      w[i] = w[i] * sqrt ( b ) / sqrt ( M_PI );
    }
  }
  

  /******************************************************************************/
  // Compute Hagedorn coefficients of projected wavepacket
  /******************************************************************************/
  for (unsigned int index=0; index < params_out->size; index++){
    evaluate_f(pot, params_in, params_out, 0, index, x, f, n);  
    params_out->c[index] = pow(M_PI * params_out->eps, - 1.0/4) / params_out->Q[0] * inner_product(w, f, n);
  }
  // return set of parameters for transmitted wavepacket 
  return params_out;// this might not work as the scope end?
}

