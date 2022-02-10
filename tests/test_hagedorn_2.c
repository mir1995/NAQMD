#include <stdio.h>
#include <complex.h>
#include <tgmath.h>


struct HagParameters{
  double complex q;
  double complex p;
  double complex Q;
  double complex P;
  double complex eps;
};


typedef struct HagParameters HagParameters;
double complex hag_gaussian_evaluate(
    double complex x, 
    HagParameters param)
{  
  return pow(M_PI * param.eps, -0.25) * pow(param.Q,-0.5) * \
    exp(I/(2*param.eps)*pow(x - param.q, 2)*param.P/param.Q + 2*param.p*(x - param.q));
}

void hag_polynomial_fill(double complex x[], double complex y[],
    HagParameters param, int index, int n)
{
  double complex z[n], temp;
  for (int i=0; i<index; i++){
    for (int j=0; j<n; j++){
      if (i==0){y[j]=1; z[j]=0;}
      temp = (sqrt(2/param.eps)/param.Q*(x[j] - param.q)*y[j] -\
             conj(param.Q)/param.Q*sqrt(index)*z[j])/sqrt(i + 1);
      z[j] = y[j];
      y[j] = temp;
    }
  }
}

void hag_gaussian_fill(double complex x[], double complex y[], 
    HagParameters param, int n){
  for (int i=0; i < n; i++){
    y[i] = hag_gaussian_evaluate(x[i], param);
  }
}

void grid_uniform_fill(double complex x[], double complex x_left, 
    double complex x_right, int n){
  double step = (x_right - x_left)/n;
  x[0] = x_left;
  for (int i=1; i<n; i++){
    x[i] = x[i-1] + step;
  }
}



int main(){
  
  struct HagParameters param;
  param.q = 4, param.p = 0, param.Q = 1, param.P = 1*I; // pass as type... eventually?
  param.eps = 0.05; 
  int n = pow(2,14), index = 10;
  double complex x[n], wave0[n], polyk[n], wavek[n];
  
  /* Create a class of parameters at some point
   * */

  grid_uniform_fill(x, param.q - 20*param.eps, param.q + 20*param.eps, n);
  hag_gaussian_fill(x, wave0, param, n);
  hag_polynomial_fill(x, polyk, param, index, n);
  
  /*
   *
   * Print parameters for information 
   * */

  printf("q = %.3f, p=%.3f, Q = %.3f + %.3fi \n",
      creal(param.q), creal(param.p),
      creal(param.Q), cimag(param.Q));
  printf("epsilon = %.6f, n=%d \n", creal(param.eps), n);

  
  /*
   *
   * Print and mass and check other observables
   */

  double complex sum;
  for (int i=0; i<n; i++){
    sum += pow(cabs(wave0[i]*polyk[i]), 2);
  }
  sum *= (x[1] - x[0]);
  printf("Mass=%.12f", creal(sum));
  /*
   *  
   *  Print Gaussian to tab-file for Gnuplot
   *
   */

  FILE *qfile = fopen("./hagedorn_wave0k.txt", "w");
  fprintf(qfile, "xgrid, wave0_r, wave0_i, \
      wave%d_r, wave%d_i", index, index);
  for (int i=0; i<n; i++){
    fprintf(qfile, "%.12f, %.12f, %.12f, %.12f, %.12f \n", 
        creal(x[i]), creal(wave0[i]), cimag(wave0[i]), 
        creal(polyk[i]), cimag(polyk[i]));
  }
  
  fclose(qfile);

  return 0;
}
