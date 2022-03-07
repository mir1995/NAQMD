# include <math.h>
# include "../Quadrature/hermite_rule.h"
# include "hag_observables.h"


static double complex inner_product(double *w, double complex *f, unsigned int n)
{
  double complex sum=0;
  for (unsigned int i=0; i<n; i++){
    sum += ((double complex) w[i]) * f[i];
  }
  return sum;
}

double complex hag_observables_get_mean(double *x, double *w, double complex *f,   
                                        struct HagedornWaves *params_in, 
                                        struct Potential *pot,
                                        unsigned int n){
  double a, alpha, b, beta;
  int kind, scale;
  double complex tau, d, v, v2, u, z, mass, mean;
  double delta, gamma, eps, s; 

  a = creal(params_in->q[0]); // 0? // depends how you do the CoV
  b = 1.0 / 2 / creal(params_in->eps);
  alpha = 0.0;
  beta = 0.0;
  scale = 0;
  kind = 6; 
  
  cgqf ( n, kind, alpha, beta, a, b, x, w );


  delta = pot->delta; 
  eps = pot->eps;
  gamma = pot->gamma; 
  tau = pot->func_get_tau(pot);

  d = pow(sin(M_PI / 2 / gamma), 2) * pow(2*M_PI*eps, - 0.5); 
  for (long unsigned int i=0; i<n; i++){
    s = x[i] / fabs(x[i]) * sqrt( pow(x[i], 2) + 4 * pot->delta);
    v2 = pow(s,2) - 4*delta; 
    if ( creal(v2) > 0 ){
      v = s / fabs(s) * sqrt(v2); // there is a sign(function) missing 
      u = pow( (s + v) / 2.0 / fabs(creal(v)) , 2);
      z = exp(I * tau * fabs(s - creal(v)) / delta / eps);

      f[i] = d * u * z * x[i] / s; // there might be a sign missing here 
    }
    else{
      f[i] = 0;
    }
    
    mass = inner_product(w, f, n);
  }
  
  for (long unsigned int i=0; i<n; i++){
    s = x[i] / fabs(x[i]) * sqrt( pow(x[i], 2) + 4 * pot->delta);
    f[i] *= s;
  }
  
  mean = inner_product(w, f, n) / mass;
  
  return mean;

}

