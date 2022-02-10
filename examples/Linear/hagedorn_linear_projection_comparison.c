# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <complex.h>
# include "../../src/Quadrature/hermite_rule.h"
# include "../../src/Potential/potential.h"
# include "../../src/Potential/potential_io.h" 
# include "../../src/Grid/grid.h" 
# include "../../src/Hagedorn/hagedorn.h" 
# include "../../src/Hagedorn/hag_io.h" 
# include "../../src/Hagedorn/hag_superadiabatic.h" 



/******************************************************************************/

int main()

/******************************************************************************/
/*
  Purpose:
    MAIN is the main program for HAGEDORN_LINEAR.C

  Discussion:

    This program compares the transmitted wavepacket 
    to a hagedorn basis approximation.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 February 2021

  Author:

    Michael Redenti
*/
{
  char fname[80];
  unsigned int n_points;
  double *x;
  double complex *f, *f_hag; 
  struct Potential pot; 
  struct HagedornWaves params_in, params_out; 
  
  timestamp ( );
  printf ( "\n" );
  printf ( "HAGEDORN_LINEAR_PROJECTION_COMPARISON\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  printf ( "\n" );
  printf ( "  Compare transmitted wavepacket to the Hagedorn representation\n" );
  printf ( "\n" );
  printf ( "  The data (Hagedorn parameters and coefficients) will be saved in:\n" );
  printf ( "    filename_x.txt .\n" );
/*
  Get Hagedorn parameters for transmitted wavepacket.
*/

  n_points = pow(2,8);
  x = ( double * ) malloc( n_points * sizeof( double ) );
  f = ( double complex * ) malloc( n_points * sizeof( double complex ) );
  f_hag = ( double complex * ) malloc( n_points * sizeof( double complex ) );
  
  // if you have a hag_read function then you would only 
  // need to specify the filename
  
  // could do hag_read().. maybe can do sizeof(params)
  // you might want to choose a different file
  // pot_read(fname, pot);
  sprintf(fname, "hag_linear_projection_potential"); 
  pot_read(fname, &pot);
  sprintf(fname, "hag_linear_projection_params_in"); 
  hag_read(fname, &params_in); 
  sprintf(fname, "hag_linear_projection_params_out"); 
  hag_read(fname, &params_out); 
  
  printf("params_in->size=%.d\n", params_in.size);
  printf("c_in[0]=%.17g\n", creal(params_in.c[0]));
  printf("c_out[0]=%.17g\n", creal(params_out.c[0]));

#ifdef LATEO
  // build momentum grid 
  grid_uniform_fill(x, -2, 4, n_points);
  // evaluate superadiabatic wavepacket at a set of points
  hag_superadiabatic_formula(x, f, &params_in, n_points, &pot, "constant");
  // evaluate hagedorn wavepacket at a set of points
  hag_wavepackets_fill(x, f_hag, &params_out, n_points); // this turns out to be the problem
  // print grid values x and difference to a file 
#endif
  free(x);
  free(f);
  free(f_hag);
  
  timestamp ( );
}
