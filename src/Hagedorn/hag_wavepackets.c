#include "hagedorn.h"
#include <complex.h>
#include <tgmath.h>
#include <stdio.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


//typedef struct HagedornParameters HagedornParameters;

double complex hag_wavepackets_gaussian_evaluate(
    double complex x, HagedornParameters params){  
  return pow(M_PI * params.eps, -0.25) * pow(params.Q,-0.5) * \
    exp(I/(2*params.eps)*(pow(x - params.q, 2)*params.P/params.Q + 2*params.p*(x - params.q)));
}

double complex hag_wavepackets_polynomial_evaluate(double complex x,
    HagedornParameters params, int index){
  double complex y = 1, z = 1, temp;
  for (int i=0; i<index; i++){
    temp = (sqrt(2/params.eps)/params.Q*(x - params.q)*y -\
           conj(params.Q)/params.Q*sqrt(i)*z)/sqrt(i + 1);
    z = y; 
    y = temp; // there could be a problem with pointers here
  }
  return y;
}

void hag_wavepackets_gaussian_fill(
    double complex x[], double complex y[], 
    HagedornParameters params, int n){
  for (int i=0; i < n; i++){
    y[i] = hag_wavepackets_gaussian_evaluate(x[i], params);
  }
}


void hag_wavepackets_polynomial_fill(double complex x[], double complex y[],
    HagedornParameters params, int index, int n){
  if (index == 0){
    for (int i=0; i<n; i++) {y[i] = 1;}
    return;
  }
  double complex z[n], temp;
  for (int i=0; i<index; i++){
    for (int j=0; j<n; j++){
      if (i==0){y[j]=1; z[j]=1;}
      temp = (sqrt(2/params.eps)/params.Q*(x[j] - params.q)*y[j] -\
             conj(params.Q)/params.Q*sqrt(i)*z[j])/sqrt(i + 1);
      z[j] = y[j];
      y[j] = temp;
    }
  }
}


