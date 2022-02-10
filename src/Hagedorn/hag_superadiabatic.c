# include <math.h>
# include <stdio.h>
# include <string.h>
# include "hagedorn.h"
# include "hag_superadiabatic.h"
# include "../Potential/potential.h" 


/******************************************************************************/

void hag_superadiabatic_formula(double *x, double complex *f,
                                struct HagedornWaves *params_in, long unsigned int n,
                                struct Potential *pot, char *version)

/******************************************************************************/
/*

  Discussion:

    This program computes the transmitted wavepacket using 
    the superadiabatic formula found in [reference]

    The user specifies:
    * f, the incoming wavepacket at the crossing 
    * probably replace this parameters hagedorn 
    * wavepackets
    *
     
    while the program modifies:
    * f, the value of the transmitted wavepacket at some values.

    Note that the Gaussian becomes part of the f(x) and not 
    the weight w(x) in the Gauss-Hermite integrnd.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 February 2021

  Author:

    Michael Redenti
*/
{ 
  // there is a different pre-factor when the 
  // eigenvalues are constant vs when they are not
  if (!strcmp(version, "constant")){
    //printf("hello");
    // apply formula for constant eigenvalues 
    double complex tau, a, v, v2, u, z, g1;
    double delta, gamma, eps; 

    delta = pot->delta; 
    eps = pot->eps;
    gamma = pot->gamma; 
    tau = pot->func_get_tau(pot);
  
    a = - sin(M_PI / 2 / gamma); 
    for (long unsigned int i=0; i<n; i++){
      v2 = pow(x[i],2) - 4*delta; 
      if ( creal(v2) > 0 ){
        v = x[i] / fabs(x[i]) * sqrt(v2); // there is a sign(function) missing 
        u = (x[i] + v) / 2.0 / fabs(creal(v));
        z = exp(I * tau * fabs(x[i] - creal(v)) / 2.0 / delta / eps);
        g1 = hag_wavepackets_gaussian_evaluate(v, params_in); // evaluate_gaussian(v) evaluate gaussian at shifted momentum value, v

        f[i] = a * u * z * g1; 
      }
      else{
        f[i] = 0;
      }
    }
  }
}
