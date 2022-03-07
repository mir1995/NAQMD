# include <stdio.h> 
# include <stdlib.h>
# include <math.h> 
# include <complex.h> 
# include "../src/Quadrature/hermite_rule.h"

#define SIZE 4

double complex inner_product(double complex *y, double *w, unsigned int n);
void generate_quadrature(double *x, double *w, double a, double b, unsigned int n);
void f(double *x, double complex *y, double b, unsigned int n);

int main(){
  

  unsigned int i, j, n; 
  double a, b, c[SIZE];
  double *x, *w;
  double complex sum;
  double complex *y; 
  
  a = 0, b = 0.5;

  for (i=1; i<5; i++){
    n = pow(10,i);
    generate_quadrature(x, w, a, b, n);
    for (j = 0; j<SIZE; j++){
      c[j] = pow(10, - 1 + j);
      f(x, y, c[j], n);
      sum = inner_product(y, w, n);
      printf("M = %d, b = %.4f -----------> I ~ %.17g + I %.17g\n", n, c[j], creal(sum), cimag(sum));
    }
  }


}



double complex inner_product(double complex *y, double *w, unsigned int n){

  double complex sum = 0; 
  
  for (unsigned int i=0; i<n; i++){
    sum += y[i] * w[i];
  }
  
  return sum; 
}

void generate_quadrature(double *x, double *w, double a, double b, unsigned int n){

  double alpha, beta;
  int kind; 

  alpha = 0, beta =0; 
  kind = 6;

  x = (double *)realloc(x, n * sizeof(double));
  w = (double *)realloc(w, n * sizeof(double));
  

  cgqf ( n, kind, alpha, beta, a, b, x, w );

}

void f(double *x, double complex *y, double b, unsigned int n){
  
  for (unsigned int i = 0; i<n; i++){
    y[i] = exp( - pow(x[i], 2) / 2 + sqrt(-1)*sqrt(2)*b*x[i] );
  }
}
