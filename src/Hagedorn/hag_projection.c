# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include "hag_projection.h"
# include "hagedorn.h"
# include "hag_superadiabatic.h"
# include "../Quadrature/hermite_rule.h"
# include "../Potential/potential.h"
# include "../Hagedorn/hag_observables.h"
# include "../../src/Grid/grid.h"


static void evaluate_f(struct Potential *pot,
                                struct HagedornWaves *params_in,
                                struct HagedornWaves *params_out,
                                unsigned int index_in, unsigned int index_out,
                                double *x, double complex *f, long unsigned int n)
{
  double complex h; // clarify the terms

  // can think about splitting it into easier tasks
  hag_superadiabatic_formula(x, f, params_in, n, pot, "constant");

  for (long unsigned int i=0; i<n; i++){
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
  double a, b, alpha, beta;
  int kind, scale;
  double *x, *w;
  double complex *f;


  w = ( double * ) malloc ( n * sizeof ( double ) ); // your main function could pass these values to a function call
  x = ( double * ) malloc ( n * sizeof ( double ) ); // these are pointers so easy to pass
  f = ( double complex * ) malloc ( n * sizeof ( double complex ) ); // these are pointers so easy to pass
  /******************************************************************************/
  // compute the new parameter set ?
  /******************************************************************************/
  // can you create a copy instead? using pointers perhaps?
  struct HagedornWaves *params_out = malloc(sizeof(struct HagedornWaves));
  params_out->state = 0;
  params_out->size = K;
  params_out->eps = params_in->eps;
  params_out->s = params_in->s;
  // energy conservation for new momentum parameter
  params_out->q[0] = sqrt(pow(params_in->q[0],2) + 4*pot->delta);
  // exact computation of new momentum parameter
  params_out->q[0] = hag_observables_get_mean(x, w, f, params_in, pot, n);
  // keep the same position parameter
  params_out->p[0] = params_in->p[0];// this is position
  params_out->Q[0] = params_in->Q[0];
  params_out->P[0] = params_in->P[0];
  params_out->c  = ( double complex * ) malloc ( K * sizeof ( double complex ) ); // these are pointers so easy to pass

  /******************************************************************************/
  // Compute weights and nodes of Gauss-Hermite rule
  /******************************************************************************/

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
    params_out->c[index] = pow(M_PI * params_out->eps, - 1.0/4) \
                           * pow(params_out->Q[0], -0.5 ) \
                           * inner_product(w, f, n);
  }
  // return set of parameters for transmitted wavepacket
  return params_out;// this might not work as the scope end?
}




struct HagedornWaves *hag_project_exact(struct HagedornWaves *params_in, long unsigned int n,
                                  unsigned int K, struct Potential *pot, double complex mean)

/******************************************************************************/
/*
  Purpose:

  Compute exact integrals for coefficients.

*/
{
  double *x;
  double complex *f;


  x = ( double * ) malloc ( n * sizeof ( double ) ); // these are pointers so easy to pass
  f = ( double complex * ) malloc ( n * sizeof ( double complex ) ); // these are pointers so easy to pass
  /******************************************************************************/
  // compute the new parameter set ?
  /******************************************************************************/
  // can you create a copy instead? using pointers perhaps?
  struct HagedornWaves *params_out = malloc(sizeof(struct HagedornWaves));
  params_out->state = 0;
  params_out->size = K;
  params_out->eps = params_in->eps;
  params_out->s = params_in->s;
  // energy conservation for new momentum parameter
  params_out->q[0] = mean;
  // keep the same position parameter
  params_out->p[0] = params_in->p[0];// this is position
  params_out->Q[0] = params_in->Q[0];
  params_out->P[0] = params_in->P[0];
  params_out->c  = ( double complex * ) malloc ( K * sizeof ( double complex ) ); // these are pointers so easy to pass

  /******************************************************************************/
  // Compute Hagedorn coefficients of projected wavepacket
  /******************************************************************************/
  grid_uniform_fill(x, -2, 5, n);
  for (unsigned int index=0; index < params_out->size; index++){
    evaluate_f(pot, params_in, params_out, 0, index, x, f, n);
    for (long unsigned int i=0; i<n; i++){
        params_out->c[index] += ( f[i] * hag_wavepackets_gaussian_evaluate(x[i], params_out) );
    }
    params_out->c[index] *= fabs(x[1] - x[0]);
  }

  return params_out;// this might not work as the scope end?
}


struct HagedornWaves **hag_projection_partition_unity(struct HagedornWaves *params, unsigned int *F){
/******************************************************************************/
/*
  Purpose:

  Computes the parameters for the families
  of Hagedorn wavepackets and then projects
  the partition of unity onto these families
  to get the coefficients.

  The programs a pointer to the family 
  of Hagedorn wavepackets.

*/

  // number of quadrature points
  unsigned int n;
  double x_left, x_right;
  double delta;
  double *x, *w;
  double complex *f;

  /* Define partition of unity domain */
  x_left = params->q[0] - 10*hag_observables_get_variance(params);
  x_right = 2*params->q[0] - x_left;

  /* support of atomic function in the partition */
  delta = sqrt(params->eps); // we need a recipe

  /* generate means */
  unsigned int M = (int) ((x_right - x_left) / delta / 2);
  *F = M;
  printf("Number of nodes n = %d", M);
  // it'd be a good idea to separate the function above, 
  // nodes, variance
  /* array of */
  struct HagedornWaves **families; // (struct HagedornWaves*)malloc(M*sizeof(struct HagedornWaves));
  families = malloc( sizeof(struct HagedornWaves*) * M );   // allocate space for 10 int *
  *families = malloc( sizeof **families * M );   // allocate space for 10 int

  double q[params->dim], p[params->dim];
  double complex Q[params->dim * params->dim], P[params->dim * params->dim];
  int K = 10;
  double complex c[K];

  /* set number of quadrature points */
  n = params->size * K; 

  // for each family of Hagedorn wavepackets
  for (unsigned int i=0; i<M; i++){
    // set parameters for family[i]
    //families[i]->q[0] = x_left + i*delta/2;
    q[0] = x_left + i*delta/2;
    p[0] = params->p[0];
    Q[0] = 1; // variance delta // the variance is given by what...?
    P[0] = 1;

    /*
     *
     * Potentially, I could end up sampling only from
     * only one Gaussian
     */
    /*
     * A product of two complex Gaussian is given by
     * A * exp{i/eps (x-q_1)^TB(x-q_1) + (x-q_1)^Tb + c}
     * where
     * A = ...
     * B = C_1 - conj(C_2), b = (p_1 - p_2) - conj(C_2)(q_1 - q_2) )
     * c = - 1/2 (q_1 - q_2)^Tconj(C_2)(q_1 - q_2) - p_2^T(q_1 - q_2)
     */

    // quadrature points - match highest degree polynomials involved?
    x = ( double * ) malloc ( n * sizeof ( double ) ); // these are pointers so easy to pass
    // weight
    w = ( double * ) malloc ( n * sizeof ( double ) ); // your main function could pass these values to a function call
    // function values
    f = ( double complex * ) malloc ( n * sizeof ( double complex ) ); // these are pointers so easy to pass

    // for each Hagedorn family i, project onto the first K wavepackets
    for (unsigned int j=0; j<K; j++){
      for (unsigned int m=0; m<n; m++){
        c[j] += w[m] * 2;
      }
      c[j] *= 1;
    }


    // you want to re-scale the original variance
    *(families + i) = hag_wavepacket_create(params->dim, K,
                                          params->state, params->eps,
                                          params->s, q, params->p,
                                          Q, P, c);

    printf("q[%d]=%.17g\n", i, families[i]->q[0]);
    // hagedorn family as input and also a pointer to a function
    // compute coefficients by calling a projection function, which takes
    // hagedorn family as input and also a pointer to a function
    //hag_projec(*func_partition, params, families[i])
  }
  free(x);
  free(w);
  free(f);
  
  return families;
}
