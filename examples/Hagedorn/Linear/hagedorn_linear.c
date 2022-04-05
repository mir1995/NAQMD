# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <complex.h>
# include "../../src/Quadrature/hermite_rule.h" 
# include "../../src/Hagedorn/hagedorn.h" 
# include "../../src/Hagedorn/hag_projection.h" 
# include "../../src/Hagedorn/hag_superadiabatic.h" 
# include "../../src/Hagedorn/hag_io.h" 
# include "../../src/Potential/potential.h"
# include "../../src/Potential/potential_io.h"
# include "../../src/Auxiliary/metrics.h"
# include "../../src/Grid/grid.h"



/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HAGEDORN_LINEAR.C

  Discussion:

    This program computes the projection of the wavepacket 
    at the crossing using a Gaussian-Hermite quadrature rule
    and writes the set of parameters and Coefficients of the 
    projected wavepacket to a file.

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

    06 February 2021

  Author:

    Michael Redenti
*/
{
  struct HagedornWaves params_in; 
  struct HagedornWaves *params_out;
  struct Potential *pot; 
  double *x;
  double complex *f, *f_hag;
  double pot_param[4];
  unsigned int n, K, power; // quadrature points, momentum grid points 
  long unsigned int n_points;
  FILE *file_error, *file_grid;
  char fname_error[200], fname_grid[200];

  
  n = atoi(argv[1]);
  K = atoi(argv[2]);
  power = 14;
  params_in.dim = DIM;
  params_in.eps = 0.1;
  params_in.s = 0;
  params_in.q[0] = 2;
  params_in.p[0] = 0;
  params_in.Q[0] = sqrt(2); // verify these are correct for the Fourier transform
  params_in.P[0] = I / sqrt(2);
  params_in.c = ( double complex * ) malloc ( 1 * sizeof ( double complex ) ); // assume initial condition is Gaussian
  params_in.size = 1; // wavepacket at crossing is a gaussian
  params_in.c[0] = 1; // wavepacket at crossing is a gaussian

  pot_param[0] = params_in.eps;
  pot_param[2] = 0.5; 
  pot_param[1] = 0.5; //  delta
  pot_param[3] = 3; // gamma
  pot = potential_construct(&v_trace, &v_z, &v_v12, &v_traced, 
                            &v_zd, &v_v12d, &v_zdd, &v_v12dd,
                            &get_tau, "Linear", pot_param);
  timestamp ( );
  printf ( "\n" );
  printf ( "HAGEDORN_LINEAR\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  printf ( "\n" );
  printf ( "  Compute projection of transmitted wavepacket \n" );
  printf ( "  for a potential with rho(x) = delta\n" );
  printf ( "  where the weight rho(x) is:\n" );
  printf ( "    exp ( - b * ( x - a )^2 ) * sqrt ( b / pi ) dx\n" );
  printf ( "  using N points.\n" );
  printf ( "\n" );
  printf ( "  The data (Hagedorn parameters and coefficients) will be saved in:\n" );
  printf ( "    filename_x.txt .\n" );
/*
  Get Hagedorn parameters for transmitted wavepacket.
*/
  params_out = malloc(sizeof(struct HagedornWaves)); 
  params_out = hag_project(&params_in, n, K, pot); 
  for (unsigned int i=0; i<K; i++){
    printf("c[%d]=%.17g\n", i, creal(params_out->c[i]) );
  }

// until we are able to save structs with pointers
  n_points = pow(2,power);
  x = (double *)malloc(n_points * sizeof(double));
  f = (double complex *)malloc(n_points * sizeof(double complex));
  f_hag = (double complex *)malloc(n_points * sizeof(double complex));
  
  // build momentum grid 
  grid_uniform_fill(x, 1, 4.5, n_points);
  // evaluate superadiabatic wavepacket at a set of points
  hag_superadiabatic_formula(x, f, &params_in, n_points, pot, "constant");
  // evaluate hagedorn wavepacket at a set of points
  hag_wavepackets_fill(x, f_hag, params_out, n_points); // this turns out to be the problem
  
  // print grid values x and difference to a file 
  
  sprintf(fname_error, "data/linear_hagedorn_projection_error_%d.txt", n); 
  sprintf(fname_grid, "data/difference_K%d.txt", K); 
  //sprintf(fname_grid, "data/wavepacket_crossing_transmitted.txt"); 
  
  file_grid = fopen(fname_grid, "w"); 
  file_error = fopen(fname_error, "a+"); 
  
  if (K == 1){
    fprintf(file_error, "K \t N \t L2_abs \t L2_rel \n");
  }
  
  //fprintf(file_grid, "p \t (f - f_hag).re  \t (f - f_hag).im\n");
  fprintf(file_grid, "p \t +.re  \t +.im \t -.re \t -.im \n");
  for (unsigned int i=0; i<n_points; i++){
    //fprintf(file_grid, "%.17g \t %.17g \t %.17g\n", x[i], creal(f[i] - f_hag[i]), cimag(f[i] - f_hag[i]) );
    fprintf(file_grid, "%.17g \t %.17g \t %.17g \t %.17g \t %.17g \n", x[i], creal(f_hag[i]), cimag(f_hag[i]), creal(f[i]), cimag(f[i]));
    //printf("f[i].re=%.17g \t f_hag.re[i]=%.17g\n", creal(f[i]), creal(f_hag[i]));
  }
  double abs_error = distance_L2(x, f, f_hag, n_points );
  double rel_error = abs_error / norm_L2(x, f, n_points);
  fprintf(file_error, "%d \t %d \t %.17g \t %.17g \n", K, n, abs_error, rel_error);
  
  free(x);
  free(f);
  free(f_hag);
  fclose(file_error);
  timestamp ( );
}
  /*
  sprintf(fname, "hag_linear_projection_potential"); 
  pot_write(fname, pot);
  sprintf(fname, "hag_linear_projection_params_in"); 
  hag_write(fname, &params_in);
  sprintf(fname, "hag_linear_projection_params_out"); 
  hag_write(fname, params_out);
  */
