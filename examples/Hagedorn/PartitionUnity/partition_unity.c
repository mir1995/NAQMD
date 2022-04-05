# include <stdlib.h>
# include <stdbool.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <complex.h>
# include "../../../src/Hagedorn/hagedorn.h" 

# define TIME 10
# define TIME_STEP 0.01
# define EPS 0.1
# define X0 0
# define K0 0 
# define SIZE 1


/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HAG_HARMONIC_OSCILLATOR.C

  Discussion:

    This program simulates the dynamics of a wavepacket 
    on a Harmonic oscillator.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Author:

    Michael Redenti
*/
{
  struct HagedornWaves params; 
  int n;
  double xl, xr, gamma, delta;
  char fname[100];
  FILE *f;


  printf ( "\n" );
  printf ( "HAG_HARMONIC_OSCILLATOR\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  printf ( "\n" );
  printf ( "\n" );
  printf("Specify parameters for initial Hagedorn wavepacket\n");
  printf ( "\n" );
  /* parameters for Hagedorn wavepacket at crossing */
  params.dim = DIM;
  params.eps = EPS;
  params.state = true;
  params.s = 0;
  params.q[0] = X0;
  params.p[0] = K0;
  params.Q[0] = sqrt(2); 
  params.P[0] = I / sqrt(2);
  params.size = SIZE; // wavepacket at crossing is a gaussian
  //params.c = ( double complex * ) malloc ( params.size * sizeof ( double complex ) ); // assume initial condition is Gaussian
  params.c[0] = 1; // wavepacket at crossing is a gaussian

  printf ( "\n" );
  printf("Parameters have been fixed.\n");
  printf ( "\n" );
  // allocate N Hagedorn data structures
  xl = params.q[0] - 5; 
  xr = params.q[0] + 5;
  delta = 10*params.eps;
  n = (int) ((xr - xl) / delta);
  gamma = 0.5;
  struct HagedornWaves *params_projection = (struct HagedornWaves *)malloc(n * sizeof(struct HagedornWaves));
  for (int i = 0; i < n; ++i) {
    (params_projection + i)->dim = params.dim;
    (params_projection + i)->eps = params.eps;
    (params_projection + i)->state = params.state;
    (params_projection + i)->s = params.s;
    (params_projection + i)->q[0] = xl + i*delta;
    (params_projection + i)->p[0] = params.p[0];
    (params_projection + i)->Q[0] = gamma * params.Q[0];
    (params_projection + i)->q[0] = 1 / gamma * params.P[0];
    // project onto lth Hagedorn wavepacket
    for (int j = 0; j < 5; j++){
      (params_projection + i)->c[j] = project();
    }
  }

  params_projection->c[0] = 5;
  params_projection->c[1] = 5;
  printf("c[1] = %.4f\n", creal(params_projection->c[1]));
  
  free(params_projection); 

}

