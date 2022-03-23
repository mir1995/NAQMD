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

    This program computes the c_0l coefficient and 
    outputs the error as a function of 
    the number of quadrature points. 

    The user specifies:
    * l, the Hagedorn wavepacket to project onto

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 February 2021

  Author:

    Michael Redenti
*/
{
  struct HagedornWaves params_in; 
  struct HagedornWaves *params_out, *params_exact;
  struct Potential *pot; 
  double pot_param[4];
  unsigned int n, K; // quadrature points, momentum grid points 
  FILE *f;
  char fname[200];

  
  n = atoi(argv[1]);
  K = atoi(argv[2]);
  // create a set up file for these parameter common
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
  

/*
  Get Hagedorn parameters for transmitted wavepacket.
*/
  params_out = malloc(sizeof(struct HagedornWaves)); 
  params_exact = malloc(sizeof(struct HagedornWaves)); 
  params_out = hag_project(&params_in, n, K, pot); 
  params_exact = hag_project_exact(&params_in, pow(2,14), K, pot, params_out->q[0]);
  for (unsigned int i=0; i<K; i++){
    sprintf(fname, "data/linear_hagedorn_c0%d.txt", i); 
    f = fopen(fname, "a+"); 
    fprintf(f, "%d \t %.17g\n", n, (cabs(params_out->c[i] - params_exact->c[i]) )\
        / cabs(params_exact->c[i]) );
  }

  fclose(f);
  timestamp ( );
}
