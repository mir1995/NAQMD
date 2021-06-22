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

void myarray_multiply_conjugate(double complex x[], double complex y[], double complex z[], int n){
  for (int i=0; i<n; i++){
    z[i] = x[i]*conj(y[i]);
  }
}


void W(double complex x[], double complex w[], int n){
  for (int i=0; i<n; i++){
    // could potentially be just multiplication
    w[i] = x[i] * pow(x[i], 3);
  }
}

double complex riemann_sum(double complex x[], double complex y[], int n){
  double complex dx = x[1] - x[0];
  double complex sum = 0;
  for (int i=0; i<n; i++){
    sum += y[i];
  }
  return sum*dx;
}

int main(int argc, char *argv[]){
  // 
  int n;
  HagedornParameters param;

  // this will be a real valued gaussian...?
  param.q=0; param.p=4; param.Q=1; param.P=1*I, param.eps=0.1;
  n=pow(2,10);
  double complex x[n];
  double complex w[n];
  double complex wave[n], polyk[n], wavek[n];
  // fix l (second index) to be zero
  int index[2] = {0,0};
  
  //fill grid 
  grid_uniform_fill(x, param.q - 10, param.q + 10, n);
  // complex gaussian
  hag_wavepackets_gaussian_fill(x, wave, param, n);
  // kth hagedorn wavepacket
  hag_wavepackets_polynomial_fill(x, polyk, param, index[0], n);
  myarray_multiply(wave, polyk, wavek, n);
  // get potential
  W(x, w, n);
  myarray_multiply(wave, w, wave, n);
  // integrate varphi_k varphi_l
  myarray_multiply_conjugate(wavek, wave, wave, n);
  //myarray_multiply_conjugate(wave, wave, wave, n);
  printf("Integral Real part %.17g\n", creal(riemann_sum(x, wave, n)));
  printf("Integral Imaginary part %.17g\n", cimag(riemann_sum(x, wave, n)));
  printf("Integral Real part %.17g\n", creal(sqrt(2/param.eps)*riemann_sum(x, wave, n)));
  printf("Integral Imaginary part %.17g\n", cimag(riemann_sum(x, wave, n)));
}
