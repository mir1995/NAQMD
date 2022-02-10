#include <stdio.h>
#include <complex.h>
#include <tgmath.h>
#include "hagedorn.h"




void grid_uniform_fill(double complex x[], double complex x_left, 
    double complex x_right, int n){
  double step = (x_right - x_left)/n;
  x[0] = x_left;
  for (int i=1; i<n; i++){
    x[i] = x[i-1] + step;
  }
}

void myarray_multiply(double complex x[], double complex y[], double complex z[], int n){
  for (int i=0; i<n; i++){
    z[i] = x[i]*y[i];
  }
}

/* Pass hamiltonian/potential so as to get transition rates  tau_c,... */
void supadia_transition_formula(double complex x[], double complex wave_in[], double complex wave_out[], 
    double complex tau, double complex delta, double complex eps, int n){
  
  double complex alpha = sin(M_PI * gamma / 2); 
  double complex step  = (x[1] - x[0]);
  for (int i=0; i<n; i++){
    if (x[i] != 0){
      double complex sign = (creal(x[i]) > 0) - (creal(x[i]) < 0);
      double complex nu = sign*sqrt(pow(x[i], 2) + 4*delta);  
      // is there a better way to multiply them?
      wave_out[(int) ((nu - x[0])/step)] =  alpha * sign * (1 + nu  / x[i]) * \
                                   exp(-tau/eps * cabs(nu - x[i])) * wave_in[i];
    }
  }
}

int main(){

  /*
   * Define parameters]
   */
  int n, M;
  double complex tau, delta, gamma, qc;
  HagedornParameters param;
 
  /*
   * Initialise structs 
   */

  // q is momentum, p is position
  param.q = 4, param.p = 0, param.Q = 1, param.P = 1*I, param.eps = 0.1;
  int index[4] = {0,1,3,10};
  n = pow(2,8);
  M = pow(2,6);
  delta = 0.5, alpha = 0.5;
  tau = 2*sqrt(pow(alpha, 2) + pow(delta,2)) * \
        (atan(alpha/delta) + atan(delta/alpha)) - alpha*M_PI;
  double complex x[n];
  double complex wave0[n], polyk[n], wavek[n], wave_out[n];
  /*
   *
   * Print parameters for information 
   */

  printf("q = %.3f, p=%.3f, Q = %.3f + %.3fi \n",
      creal(param.q), creal(param.p),
      creal(param.Q), cimag(param.Q));
  printf("epsilon = %.6f, n=%d \n", creal(param.eps), n);
  
  /*
   *
   * Generate wavepackets at crossing ---> transmit ---> save
   *
   */
  grid_uniform_fill(x, param.q - 20*param.eps, param.q + 20*param.eps, n);
  hag_wavepackets_gaussian_fill(x, wave0, param, n);
  
  
  for (int i=0; i< (size_t) (sizeof(index)/sizeof(index[0])); i++){
    double complex mass = 0;
    /* Apply transmission formula */
    hag_wavepackets_polynomial_fill(x, polyk, param, index[i], n);
    /* Apply transmission formula */
    myarray_multiply(wave0, polyk, wavek, n);
    /* Apply transmission formula */
    supadia_transition_formula(x, wavek, wave_out, tau, delta, param.eps, gamma, n); 
    
    char filename[100];
    sprintf(filename, "./data/hag_wave%d%.3f_transmitted.txt", index[i], creal(delta));
    FILE *qfile = fopen(filename, "w");
    fprintf(qfile, "xgrid, wave%d_r, wave%d_i waveout%d_r, waveout%d_i\n",
        index[i], index[i], index[i], index[i]);
  /*
   *
   * Print and mass and check other observables
   */
    for (int i=0; i<n; i++){
      fprintf(qfile, "%.12f, %.12f, %.12f, %.12f, %.12f \n", 
          creal(x[i]), creal(wavek[i]), cimag(wavek[i]),
          creal(wave_out[i]), cimag(wave_out[i]));
      mass += pow(cabs(wavek[i]), 2);
    }
    mass *= (x[1] - x[0]);
    printf("Mass=%.12f", creal(mass));
    fclose(qfile);
  }



  /*
   *
   *
   * Project transmitted wavepacket onto new basis - integrate
   *
   */


  /*
   *  
   *  Print Gaussian to tab-file for Gnuplot
   *
   */

  return 0;
}
