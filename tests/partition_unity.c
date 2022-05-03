# include <stdlib.h>
# include <stdbool.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <complex.h>
# include "../src/Hagedorn/hagedorn.h" 
# include "../src/PartitionUnity/partition_unity.h"
# include "../src/Hagedorn/hag_io.h"

# define TIME 10
# define TIME_STEP 0.01
# define EPS 0.1
# define X0 0
# define K0 0 
# define SIZE 1
# define STATE true


/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PARTITION_UNITY.C

  Discussion:

    This program tests the decomposition of a broad Hagedorn 
    wavepacket into a families of narrower Hagedorn wavepackets
    passing through a partition of unity. 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Author:

    Michael Redenti
*/
{
  struct HagedornWaves *params; 
  double s, q[DIM], p[DIM];
  double complex Q[DIM * DIM], P[DIM * DIM], c[SIZE];   
  int n;
  double xl, xr, gamma, delta;
  char fname[100];
  FILE *f;


  printf ( "\nPARTITION OF UNITY\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  /* parameters for Hagedorn wavepacket at crossing */
  q[0] = X0, p[0] = K0, Q[0] = sqrt(2), P[0] = I / sqrt(2);
  s = 0, c[0] = 1; // wavepacket at crossing is a gaussian

  /* HOW DO I CHANGE THE VARIANCE */

  /* Generate initial condition */
  params = hag_wavepacket_create(DIM, SIZE, STATE, EPS, s, q, p, Q, P, c);

  hag_write("partition_unity_hagedorn_wavepacket_parameters", params);
 
  unsigned int F;
  struct HagedornWaves **families = hag_projection_partition_unity(params, &F);
  for (unsigned int i=0; i<F; i++){
    hag_write("partition_unity_hagedorn_wavepacket_parameters_families_sigma2", families[i]);
  }
  // return variance of wavepacket assume it's gaussian
  // look up Caroline's paper for tails norm
  /* families of Hagedorn Wavepackets */
  
  
  /* Quadrature rule start */

  /* pointer to an array of families of Hagedorn Wavepackets*/
  //struct HagedornWaves *families = hag_projection_partition_unity(params);
  

  // build a routine to get x_l, x_r, nodes
  // no matter what the number of families 
  // should always at least be two otherwise 
  // there is no point
  

  // projection onto family of HagedornWaves 
  // using a partition of unity

  // allocate N Hagedorn data structures
  /*
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

  */

}

