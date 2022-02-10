# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <complex.h>
# include "../src/Quadrature/hermite_rule.h" 
# include "../src/Hagedorn/hagedorn.h" 



/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HAGEDORN_LINEAR.C

  Discussion:

    This program computes the inner product 
    <varphi_k, varphi_l> where varphi_k is an Hagedorn 
    wavepacket of unit norm. 
    
    To evaluate the integral we use a Gaussian-Hermite 
    quadrature rule with N=floor( ( k+l-1 ) / 2 )
    points. This value of N should yield the exact 
    value. 

    Can you say anything about the convergence for the 
    number of points when the integrand is a polynomial.


    The user specifies:
    * N, the number of points in the rule
    * A, the center point;
    * B, a scale factor;
    * SCALE, is 1 if the weights are to be normalized;
    * FILENAME, the root name of the output files.

    If SCALE = 0, then the factor C in front of the integrand is 1.
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
  struct HagedornWaves params; 
  double a, b, alpha, beta, a_kl; 
  int kind, scale; 
  double *w, *x;
  double complex *h_k, *h_l;
  unsigned int n, K, N;
  

  K = 5;
  N = ceilf( ( 2.0 * K + 1.0 ) / 2 );
  printf("N=%d\n", N);
  
  params.dim = DIM;
  params.size = 1; // wavepacket at crossing is a gaussian
  params.s = 0;
  params.eps = 0.1;
  params.q[0] = 2;
  params.p[0] = 0;
  params.Q[0] = sqrt(2); // verify these are correct for the Fourier transform
  params.P[0] = I / sqrt(2);
  params.c = ( double complex * ) malloc ( params.size * sizeof ( double complex ) ); // assume initial condition is Gaussian
  params.c[0] = 1; // wavepacket at crossing is a gaussian
  
  a = creal(params.q[0]);
  b = 1.0 / 2.0 / creal(params.eps);
  alpha = 0.0;
  beta = 0.0;
  scale = 0;
  kind = 6; // Gauss Hermite Rule - a name would have been better 

  w = ( double * ) malloc ( N * sizeof ( double ) );
  x = ( double * ) malloc ( N * sizeof ( double ) );
  h_k = ( double complex * ) malloc ( N * sizeof ( double complex ) );
  h_l = ( double complex * ) malloc ( N * sizeof ( double complex ) );
  
  timestamp ( );
  printf ( "\n" );
  printf ( "TEST_HERMITE_HAGEDORN\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  printf ( "\n" );
  
  for (unsigned int k=0; k<K; k++){
    printf ( "\n" );
    for (unsigned int l=k; l<K; l++){
      a_kl = 0;
      // sufficient number of points
      n = ceilf( ( k + l + 1.0 ) / 2.0 );
      // generate N hermite points
      cgqf ( n, kind, alpha, beta, a, b, x, w ); // use the points generated from Chebishev polynomial and compare
      // evaluate H_{k},H_{l} at nodes X 
      //hag_wavepackets_polynomial_fill(x, h_k, &params, k, n);
      //hag_wavepackets_polynomial_fill(x, h_l, &params, l, n);
      // inner product with weights W 
      for (unsigned int i=0; i<n; i++){
        h_k[i] = hag_wavepackets_polynomial_evaluate(x[i], &params, k);
        h_l[i] = hag_wavepackets_polynomial_evaluate(x[i], &params, l); // we have detected the problem - it is in the function HAG_WAVEPACKETS_POLYNOMIAL_FILL
        h_l[i] *= conj(h_k[i]);
        a_kl += (h_l[i] * w[i]);
      }
      
      a_kl /= params.Q[0] ;
      a_kl *= pow(M_PI * params.eps, -0.5) ;
      printf("%.17g \t \t", creal(a_kl));  // try later without a_{%d,%d} 
    }
  }
  
  free(w);
  free(x);
  free(h_k);
  free(h_l);


  timestamp ( );


}
