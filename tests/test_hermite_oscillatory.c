# include <stdio.h> 
# include <stdlib.h>
# include <math.h> 
# include <complex.h> 
# include "../src/Quadrature/hermite_rule.h"

#define SIZE 4

int main(){
  

  unsigned int i, j, n; 
  double a, b, c[SIZE];
  double *x, *w;
  double complex *y; 

  for (i=1; i<5; i++){
    n = pow(10,i);
    generate_quadrature(x, w, a, b, n);
    for (j = 0; j<SIZE; j++){
      c[j] = pow(10, - 1 + j);
    }
  }


}
