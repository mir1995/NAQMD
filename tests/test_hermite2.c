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
  char filename[80];
  int i;
  int kind;
  int n;
  double pi = 3.14159265358979323846264338327950;
  double *r;
  int scale;
  double *w;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "HERMITE_RULE\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  printf ( "\n" );
  printf ( "  Compute a Gauss-Hermite quadrature rule for\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) rho(x) dx\n" );
  printf ( "  where the weight rho(x) is:\n" );
  printf ( "    exp ( - b * ( x - a )^2 ) * sqrt ( b / pi ) dx\n" );
  printf ( "  using N points.\n" );
  printf ( "\n" );
  printf ( "  The user specifies N, A, B, SCALE, and FILENAME.\n" );
  printf ( "\n" );
  printf ( "  N is the number of points;\n" );
  printf ( "  A is the center point, usually 0.0:\n" );
  printf ( "  B is a scale factor, usually 0.5 or 1.0;\n" );
  printf ( "  SCALE is 1 if the weights are to be normalized to unit sum.\n" );
  printf ( "  FILENAME is used to generate 3 files:\n" );
  printf ( "    filename_w.txt - the weight file\n" );
  printf ( "    filename_x.txt - the abscissa file.\n" );
  printf ( "    filename_r.txt - the region file.\n" );
/*
  Initialize parameters;
*/
  alpha = 0.0;
  beta = 0.0;
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
  Get A.
*/
  if ( 2 < argc )
  {
    a = atof ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter A, the center point (default 0)\n" );
    scanf ( "%lf", &a );
  }
/*
  Get B.
*/
  if ( 3 < argc )
  {
    b = atof ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter B, the scale (default 1)\n" );
    scanf ( "%lf", &b );
  }
/*
  Get SCALE.
*/
  if ( 4 < argc )
  {
    scale = atoi ( argv[4] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter SCALE, the normalization option (0/1)\n" );
    scanf ( "%d", &scale );
  }
/*
  Get FILENAME:
*/
  if ( 5 < argc )
  {
    strcpy ( filename, argv[5] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter FILENAME, the \"root name\" of the quadrature files).\n" );
    scanf ( "%s", filename );
  }
/*
  Input summary.
*/
  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "  SCALE = %d\n", scale );
  printf ( "  FILENAME = \"%s\"\n", filename );
/*
  Construct the rule.
*/
  w = ( double * ) malloc ( n * sizeof ( double ) ); // your main function could pass these values to a function call
  x = ( double * ) malloc ( n * sizeof ( double ) ); // these are pointers so easy to pass
  
  kind = 6;
  cgqf ( n, kind, alpha, beta, a, b, x, w );
/*
  Normalize the rule.
  This way, the weights add up to 1.
*/
  if ( scale == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      w[i] = w[i] * sqrt ( b ) / sqrt ( pi );
    }
  }
/*
  Set the R values to suggest an infinite interval.
*/
  r = ( double * ) malloc ( 2 * sizeof ( double ) );
  r[0] = - r8_huge ( );
  r[1] =   r8_huge ( );
/*
  Write the rule to a file.
*/
  rule_write ( n, filename, x, w, r );
/*
  Free memory.
*/
  free ( r );
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


