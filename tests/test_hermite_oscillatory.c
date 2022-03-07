# include <stdio.h> 
# include <stdlib.h>
# include <math.h> 
# include <complex.h> 
# include "../src/Quadrature/hermite_rule.h"

#define SIZE 4

double complex inner_product(double complex *y, double *w, unsigned int n);
void f(double *x, double complex *y, double b, unsigned int n);

int main(){
  

  unsigned int i, j, n; 
  double a, b, c[SIZE];
  double *x, *w;
  double complex sum;
  double complex *y; 
  FILE *file; 
  
  a = 0.0, b = 0.5;
  c[0] = 0.1, c[1] = 1, c[2] = 3, c[3] = 6;
 
  file = fopen("hermite_oscillatory.txt", "w");
  fprintf(file, "b \t N \t |I(N) - I| \n");
  
  for (i=4; i<12; i++){
    n = pow(2,i);
    y = (double complex *)malloc(n * sizeof(double complex));
    x = (double *)malloc(n * sizeof(double));
    w = (double *)malloc(n * sizeof(double));
    cgqf ( n, 6, 0.0, 0.0, a, b, x, w );
    for (j = 0; j<SIZE; j++){
      //c[j] = pow(10, - 1.0 + j);
      f(x, y, c[j], n);
      sum = inner_product(y, w, n) / sqrt(M_PI) / exp(- pow( c[j], 2 ) / 2 ) ;
      printf("M = %d, b = %.4f -----------> I ~ %.17g + I %.17g\n", n, c[j], creal(sum), cimag(sum));
      fprintf(file, "%.4f \t %d \t %.17g \n", c[j], n, cabs(sum - 1 ) );
    }
    free(y);
    free(x);
    free(w);
  }
  fclose(file);
}



double complex inner_product(double complex *y, double *w, unsigned int n){

  double complex sum = 0; 
  
  for (unsigned int i=0; i<n; i++){
    sum += y[i] * w[i];
  }
  
  return sum; 
}

void f(double *x, double complex *y, double b, unsigned int n){
  
  for (unsigned int i = 0; i<n; i++){
    y[i] = cexp( - pow(x[i], 2)  / 2.0 + I * sqrt(2)*b*x[i] );
  }
}
