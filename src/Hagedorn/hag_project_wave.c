#include <stdio.h>
#include <stdlib.h>
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
    double complex tau, double complex delta, double complex eps,  double complex gamma, int n){
  
  double complex alpha = sin(M_PI * gamma / 2); 
  double complex step  = (x[1] - x[0]);
  int start = (int) ((- 2 - x[0])/step); // problem with accessing values
  for (int i=start; i<n; i++){
    if (x[i] != 0){
      double complex sign = (creal(x[i]) > 0) - (creal(x[i]) < 0);
      double complex nu = sign*sqrt(cpow(x[i], 2) + 4*delta);  
      // is there a better way to multiply them?
      wave_out[(int) ((nu - x[0])/step)] =  alpha * sign * (1 + nu  / x[i]) * \
                                   exp(-tau/eps * cabs(nu - x[i])) * wave_in[i];
    }
  }
}

int main(int argc, char *argv[]){
  
  /*
   * Define parameters]
   */
  int n, M;
  n = pow(2,14);
  double complex x[n], wave_0[n], wave_out[n];
  double complex tau, delta, gamma, qc;
  HagedornParameters param_0, param_l;
 
  /*
   * Initialise structs 
   */

  // q is momentum, p is position
  delta = 0.1;
  qc = 0.5; 
  tau = 2*delta*qc; 
  gamma = - qc/2;
  param_0.q = 4, param_0.p = 0, param_0.Q = 1, param_0.P = 1*I, param_0.eps = 0.1;
  param_l.q = 4, param_l.p = 0, param_l.Q = 1, param_l.P = 1*I, param_l.eps = 0.1;
  M = pow(2,10); // number of Monte Carlo samples
  /*
   * Print parameters for information 
   */
  printf("q = %.3f, p=%.3f, Q = %.3f + %.3fi \n",
      creal(param_0.q), creal(param_0.p),
      creal(param_0.Q), cimag(param_0.Q));
  printf("epsilon = %.6f, n=%d \n", creal(param_0.eps), n);
  
  
  /*
   *
   *
   * Project transmitted wavepacket onto new basis - integrate
   *
   */
  int L = 20; // number of coefficients to compute
  double complex nu, c[L]; 
  double complex eps = param_0.eps;
  double complex mass = 0;
  #define RIEMANN TRUE
  #ifdef RIEMANN
    /*
     * Generate complex Gaussian and apply transmission formula
     */
    grid_uniform_fill(x, param_0.q - 8, param_0.q + 12, n);
    hag_wavepackets_gaussian_fill(x, wave_0, param_0, n);
    supadia_transition_formula(x, wave_0, wave_out, tau, delta, eps, gamma, n);
    int start = (int) ( (2*sqrt(delta) - x[0])/(x[1] - x[0]) );
    //int stop = (int) ( (6 - x[0])/(x[1] - x[0]));
    /*
     *
     * Compute mass transmitted wavepacket to normalise things
     */
    //param_l.q = 0;
    for (int i = start; i<n; i++){
      mass += pow(cabs(wave_out[i]),2); 
      //param_l.q += pow(cabs(wave_out[i]),2)*x[i];
    }
    mass *= (x[1] - x[0]);
    //param_l.q *= ((x[1] - x[0])/mass);
    param_l.q = sqrt(pow(param_0.q,2) + 4*delta);
    //param_l.q = param_0.q;
    printf("mass %.5f\n", creal(mass));
    printf("p_prime %.5f\n", creal(param_l.q));
    for (int l=0; l< L; l++){
      c[l] = 0;
      double complex poly_l[n], wave_l[n];
      hag_wavepackets_gaussian_fill(x, wave_l, param_l, n);
      hag_wavepackets_polynomial_fill(x, poly_l, param_l, l, n);
      myarray_multiply(wave_l, poly_l, wave_l, n); 
      for (int i=start; i<n; i++){ // not all of them
        c[l] += conj(wave_l[i]) * wave_out[i];
      }
      c[l] *= ((x[1] - x[0])/sqrt(mass));
      // printf("Mass=%.12f", creal(mass));
    }
    printf("Result of Riemann Integration \n");
    for (int l=0; l<L; l++){
      printf("c_{%d} = %.12f + i%.12f \n", l, creal(c[l]), cimag(c[l]));
    }
    char filename[100];
    sprintf(filename, "./data/hag_coefficients0%.3f_conservation.txt", creal(delta));
    FILE *qfile = fopen(filename, "w");
    fprintf(qfile, "index, c_real, c_imag \n");
    for (int l=0; l<L; l++){
      fprintf(qfile, "%d, %.12f, %.12f \n", l, creal(c[l]), cimag(c[l]));
    }
    fclose(qfile);
  #else
    double complex alpha, beta, s; // s is a sample, c stores the coeffs
    double complex Im, Re;
    double complex p0, p1, std;

    Im = 1, Re = 1;
    p0 = param.q, p1 = p0 + 2*sqrt(delta) , std = sqrt(eps/2/Im);
    // constant prefactor
    alpha = sin(M_PI * gamma / 2) * \
            exp(-1/(2*param.eps) * (Im*(cpow(p1,2)/2 + cpow(p0,2) - 4*delta)) - \
            I * (Re*(cpow(p0,2) - cpow(p1,2) - 4*delta)));
    printf("%.3f", creal(alpha));
    printf("hello \n");
    for (int m=0; m < M; m++){
      s = -2;
      while (creal(s) <= 2*sqrt(creal(delta))) {
        s = p1/2 + std * rand_gaussian_sample();
      }
      printf("hello \n");
      nu = sqrt(cpow(s,2) - 4 * delta);
      beta = (1 + s/(sqrt(cpow(s,2) - 4*delta))) * \
             hag_wavepackets_polynomial_evaluate(nu, param, index[0]) * \
             exp(-1/eps * (s * (qc - I * Re * p1) - nu * (p0*Im + qc - I*Re*p0) ) );
      printf("%.3f", creal(beta));
      for (int l=0; l<L; l++){
        c[l] += alpha * beta * hag_wavepackets_polynomial_evaluate(s, param, l);
      }
    }
  
    printf("Result of Monte Carlo Integration \n");
    for (int i=0; i<L; i++){
      c[i] /= M; // average over MC samples
      printf("c_{%d} = %.9f \n", i, creal(c[i]));
    }
  #endif 
  // you want to also compare against not 
  // computing a different momentum value

  return 0;
}

