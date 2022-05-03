# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include "partition_unity.h"
# include "../Hagedorn/hagedorn.h"
# include "../Quadrature/hermite_rule.h"
# include "../Hagedorn/hag_observables.h"


static double partition_unity_function(double x, double mean, double delta){
  
  double diff = fabs(x - mean);  
  if (diff < delta / 2){
    double a, b;
    a = ( delta / 2 ) / (- diff + delta / 2);
    b = ( delta / 2 ) / diff ;
    return exp(-a) / (exp(-a) + exp(-b));
  }

  return 0;
}

static double complex hag_polynomial_evaluate(double x, double eps, double q, 
                                              double complex Q, int index){
  double complex y = 1, z = 1, temp;
  for (int i=0; i<index; i++){
    temp = (sqrt(2/eps)/Q*(x - q)*y -\
           conj(Q)/Q*sqrt(i)*z)/sqrt(i + 1);
    z = y; 
    y = temp; // there could be a problem with pointers here
  }
  return y;
}

  struct HagedornWaves **hag_projection_partition_unity(struct HagedornWaves *params, 
                                                        unsigned int *F){
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

  double x_left, x_right;
  double delta, var;
  // number of quadrature points
  unsigned int n;
  double *x, *w;
  double complex *f;

  /* Pick grid endpoints for partition of unity */
  // they should also depend on the largest order of the Hagedorn wavepacket
  var = hag_observables_get_variance(params); 
  x_left = params->q[0] - 4*sqrt(var);
  x_right = 2*params->q[0] - x_left;
  printf("xl = %.4f, xr = %.4f\n", x_left, x_right);

  /* Define support of each function in the partition  */
  delta = 0.5*sqrt(var); // we need a recipe
  printf("delta=%.4f\n", delta);

  /* Number of nodes on the uniform grid */
  unsigned int M = (int) ((x_right - x_left) * 2 / delta);
  *F = M;
  printf("Number of nodes n = %d\n", M);

  // it'd be a good idea to separate the function above, 
  // nodes, variance
  /* array of */
  struct HagedornWaves **families; // (struct HagedornWaves*)malloc(M*sizeof(struct HagedornWaves));
  families = malloc( sizeof(struct HagedornWaves*) * M );   // allocate space for 10 int *
  *families = malloc( sizeof **families * M );   // allocate space for 10 int

  double q[params->dim], p[params->dim];
  double complex Q[params->dim * params->dim], P[params->dim * params->dim];
  int K =  50;
  double complex c[K];

  /* set number of quadrature points */
  n = 100; // params->size * K; 
  // quadrature points - match highest degree polynomials involved?
  x = ( double * ) malloc ( n * sizeof ( double ) ); // these are pointers so easy to pass
  // weight
  w = ( double * ) malloc ( n * sizeof ( double ) ); // your main function could pass these values to a function call
  // function values
  f = ( double complex * ) malloc ( n * sizeof ( double complex ) ); // these are pointers so easy to pass


  // for each family of Hagedorn wavepackets
  for (unsigned int i=0; i<M; i++){
    // set parameters for family[i]
    //families[i]->q[0] = x_left + i*delta/2;
    q[0] = x_left + i*delta/2;
    p[0] = params->p[0];
    // it makes sense to scale the variance by the number of nodes?
    //Q[0] = sqrt(2), P[0] = I / sqrt(2)
    double sigma = 2; //1 ,2 ,4
    Q[0] = 1 / sqrt(params->eps) * (delta/2/sigma); // variance delta // the variance is given by what...?
    /* need to satisfy symplectic constraints */
    P[0] = I / Q[0];

    /*
     * Now, we are ready to project onto this family of 
     * Hagedorn wavepackets
     * Potentially, I could end up sampling only from
     * only one Gaussian
     */
    /*
     * A product of two complex Gaussian is given by
     * A * exp{i/eps ( (x-q_1)^TB(x-q_1) + (x-q_1)^TC + D ) }
     * where
     * A = (pi*eps)^{-d/2} (det Q)
     * B = C_1 - conj(C_2), 
     * C = (p_1 - p_2) - conj(C_2)(q_1 - q_2) )
     * D = - 1/2 (q_1 - q_2)^Tconj(C_2)(q_1 - q_2) - p_2^T(q_1 - q_2)
     */

    double a, b;
    double complex C1, C2;
    a = q[0];
    printf("Q0 %.4f + 1j %.4f\n", creal(Q[0]), cimag(P[0]));
    C2 = P[0]/Q[0];
    printf("C2 %.4f + 1j %.4f\n", creal(C2), cimag(C2));
    C1 = params->P[0]/params->Q[0]; 
    printf("C1 %.17G + 1j %.17g\n", creal(C1), cimag(C1));
    b = cimag(C2 - conj(C1)) / 2 / params->eps;
    printf("b  %.4f\n", b);
    cgqf ( n, 6, 0, 0, a, b, x, w );
    
    double complex q1, p1, q2, p2;
    double complex q_diff, p_diff;
    double complex A, B, C, D;
    q1 = params->q[0], p1 = params->p[0];
    q2 = q[0], p2 = p[0];
    q_diff = q2 - q1; 
    p_diff = p2 - p1;
    // for each Hagedorn family i, project onto the first K wavepackets
    for (unsigned int j=0; j<K; j++){
      c[j] = 0;
      for (unsigned int m=0; m<n; m++){
        /* re-write more clearlt */
        // this term is highly oscillatory
        f[m] = cexp(I/params->eps * (x[m] - q2) * ( p_diff - conj(C1) * q_diff ) ); 
        //printf("value f[m]=%.4f , x[m]=%.4f\n", creal(f[m]), x[m]);
        f[m] *= cexp(I / 2.0 / params->eps * (x[m] - q2) * creal( C2 - conj(C1)) * (x[m] - q2)); 
        /* Multiply by j^th Hagedorn polynomial */
        f[m] *= hag_polynomial_evaluate(x[m], params->eps, q[0], Q[0], j);
        //printf("value f[m]=%.4f , x[m]=%.4f\n", creal(f[m]), x[m]);
        f[m] *= partition_unity_function(x[m], q[0], delta);
        printf("value f[m]=%.4f , x[m]=%.4f, w[m]=%.4f\n", creal(f[m]), x[m], w[m]);
        c[j] += w[m]*f[m];
      }
      // could potentially write this outside
      printf("c[%d]=%.17g\n", j, creal(c[j]));
      A = pow(M_PI * params->eps, - (params->dim*1.0) / 2 )*\
          cpow(cabs(params->Q[0]), - 0.5)*\
          cpow(cabs(Q[0]), - 0.5);
      D =  cexp(I/params->eps * ( - 1.0 / 2 * q_diff * conj(C1) * q_diff - p1*q_diff ) ); 
      c[j] *= A * D;
      printf("c[%d]=%.17g\n", j, creal(c[j]));
    }


    // you want to re-scale the original variance
    *(families + i) = hag_wavepacket_create(params->dim, K,
                                          params->state, params->eps,
                                          params->s, q, p,
                                          Q, P, c);

    printf("q[%d]=%.17g\n", i, families[i]->q[0]);
    printf("c[%d]=%.17g\n", i, creal(families[i]->c[0]));
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
