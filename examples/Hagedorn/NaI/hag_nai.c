# include <stdlib.h>
# include <stdbool.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <complex.h>
# include "../../../src/Hagedorn/hagedorn.h" 
# include "../../../src/Hagedorn/hag_dynamics.h"
# include "../../../src/Potential/potential.h"
# include "../../../src/Odeint/odeint.h"

# define TIME 80
# define TIME_STEP 0.01
#define ALPHA 0.002669794983146837
#define DELTA 0.002010562925688143 
#define EPS 0.014654629670711006
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
  struct Potential *pot; 
  double pot_param[3];
  char fname[100];
  FILE *f;


  printf ( "\n" );
  printf ( "HAG_NAI\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__  );
  printf ( "\n" );
  // set_up(filename) - sets up all paramters for method
  // write a program that set up all the objects by reading an input textfile of parameters
  // kind of like MOLPRO
  printf ( "\n" );
  printf("Specify parameters for initial Hagedorn wavepacket\n");
  printf ( "\n" );
  params.dim = DIM;
  params.eps = EPS;
  params.state = true;
  params.s = 0;
  params.q[0] = 5;
  params.p[0] = 0;
  params.Q[0] = sqrt(2); // verify these are correct for the Fourier transform
  params.P[0] = I / sqrt(2);
  params.size = SIZE; // wavepacket at crossing is a gaussian
  params.c = ( double complex * ) malloc ( params.size * sizeof ( double complex ) ); // assume initial condition is Gaussian
  params.c[0] = 1; // wavepacket at crossing is a gaussian

  printf ( "\n" );
  printf("Construct Harmonic oscillator potential\n");
  printf ( "\n" );
  pot_param[0] = EPS;
  pot_param[1] = DELTA;
  pot_param[2] = ALPHA;
  pot = potential_construct(&v_trace, &v_z, &v_v12, 
                            &v_traced, &v_zd, &v_v12d, 
                            &v_tracedd, &v_zdd, &v_v12dd,
                            &get_tau, "NaI", 
                            pot_param);
  
  printf ( "\n" );
  printf("Set up time discretisation\n");
  printf ( "\n" );
  struct Odeint *odeint = odeint_new(TIME, TIME_STEP, DIM, "none");

  printf ( "\n" );
  printf("Evolve Hagedorn wavepacket for T=%f (a.u.)\n", odeint->t);
  printf ( "\n" );

  sprintf(fname, "data.txt");
  f = fopen(fname, "w");
  


  fprintf(f, "Time[a.u.] \t S \t q[0] \t p[0] \t Q[0].re \t Q[0].im \t P[0].re \t P[0].im \n");
  // i do not want to write for(int itr bla ...)) every time. 
  for (unsigned int itr=0; itr< odeint->t / odeint-> dt; itr++){
    if (itr % ( (int) (odeint->t / odeint-> dt / 40) ) == 0 ){
      fprintf(f, "%.17g \t %.17g \t %.17g \t %.17g \t %.17g \t %.17g \t %.17g \t %.17g \n",
                  itr * odeint->dt,
                  params.s, params.q[0], params.p[0], 
                  creal(params.Q[0]), cimag(params.Q[0]), 
                  creal(params.P[0]), cimag(params.P[0]));
    }
    hag_dynamics_do_step(&params, pot, odeint);
  }

  fclose(f);
}

