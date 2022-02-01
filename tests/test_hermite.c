# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include "../src/Quadrature/hermite_rule.h"

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HERMITE_RULE.

  Discussion:

    This program computes a Gauss-Hermite quadrature rule 
    and writes it to a file.

    The integral to be approximated has the form

      C * Integral ( -oo < x < +oo ) f(x) rho(x) dx

    where the weight rho(x) is:

      rho(x) = exp ( - b * ( x - a )^2 ) * sqrt ( b / pi ) dx

    and A and B are parameters.

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

    06 February 2014

  Author:

    John Burkardt
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int kind;
  int n;
  double pi = M_PI;
  int scale;
  double *w;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  printf ( "\n" );
  printf ( "  Compute a Gauss-Hermite quadrature rule for\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) e^{-x^2} dx\n" );
  printf ( "  using N points.\n" );
  printf ( "\n" );
  printf ( "  The user specifies N, A, B, SCALE\n" );
  printf ( "\n" );
  printf ( "  N is the number of points;\n" );
  printf ( "  A is the center point, usually 0.0:\n" );
  printf ( "  B is a scale factor, usually 0.5 or 1.0;\n" );
  printf ( "  SCALE is 1 if the weights are to be normalized to unit sum.\n" );
/*
  Initialize parameters;
*/
  alpha = 0.0;
  beta = 0.0;
  a = 0;
  b = pow(M_PI , 1/3);
  scale = 0;
  kind = 6; // Gauss Hermite Rule - a name would have been better 
/*
  Get N.
*/
  if ( 1 < argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter N (1 or greater)\n" );
    scanf ( "%d", &n );
  }
/*
  Input summary.
*/
  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "  SCALE = %d\n", scale );
/*
  Construct the rule.
*/
  w = ( double * ) malloc ( n * sizeof ( double ) ); // your main function could pass these values to a function call
  x = ( double * ) malloc ( n * sizeof ( double ) ); // these are pointers so easy to pass
  
  cgqf ( n, kind, alpha, beta, a, b, x, w );
/*
  Normalize the rule.
  This way, the weights add up to 1.
*/
  if ( scale == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      w[i] = w[i] * sqrt ( b ) / sqrt ( M_PI );
    }
  }
/*
  compute integral
*/
  double f(double x){return 1;}
  double sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum += w[i] * f(x[i]);
    printf("\n x[i] = %.17g \n", x[i]);
    printf("   w[i] = %.17g \n", w[i]);
  }
  printf("Value of integral is %.17g", sum);
/*
  Free memory.
*/
  free ( w );
  free ( x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HERMITE_RULE:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}


