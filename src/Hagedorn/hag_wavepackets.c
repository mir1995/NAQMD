#include "hagedorn.h"
#include <complex.h>
#include <tgmath.h>
#include <stdio.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


//typedef struct HagedornParameters HagedornParameters;
// for the moment you are assuming one dimension
double complex hag_wavepackets_gaussian_evaluate(
    double complex x, HagedornWaves *params){  
  return pow(M_PI * params->eps, -0.25) * pow(params->Q[0],-0.5) * \
    exp(I/(2*params->eps)*(pow(x - params->q[0], 2)*params->P[0]/params->Q[0] + 2*params->p[0]*(x - params->q[0])));
}

double complex hag_wavepackets_polynomial_evaluate(double complex x,
    HagedornWaves *params, int index){
  double complex y = 1, z = 1, temp;
  for (int i=0; i<index; i++){
    temp = (sqrt(2/params->eps)/params->Q[0]*(x - params->q[0])*y -\
           conj(params->Q[0])/params->Q[0]*sqrt(i)*z)/sqrt(i + 1);
    z = y; 
    y = temp; // there could be a problem with pointers here
  }
  return y;
}

void hag_wavepackets_gaussian_fill(
    double complex x[], double complex y[], 
    HagedornWaves *params, int n){
  for (int i=0; i < n; i++){
    y[i] = hag_wavepackets_gaussian_evaluate(x[i], params);
  }
}


void hag_wavepackets_polynomial_fill(double complex x[], double complex y[],
    HagedornWaves *params, int index, int n){
  if (index == 0){
    for (int i=0; i<n; i++) {y[i] = 1;}
    return;
  }
  double complex z[n], temp;
  for (int i=0; i<n; i++){
    for (int j=0; j<index; j++){
      if (j==0){
        y[i]=1; z[i]=1;
      }
      else{
        temp = (sqrt(2/params->eps)/params->Q[0]*(x[i] - params->q[0])*y[i] -\
               conj(params->Q[0])/params->Q[0]*sqrt(j)*z[i])/sqrt(j + 1);
        z[i] = y[i];
        y[i] = temp;
      }
    }
  }
}

void hag_wavepackets_fill(double x[], double complex y[],
    HagedornWaves *params, long unsigned int n){
  for (long unsigned int i=0; i<n; i++){
    y[i] = 0;
    for (unsigned int k=0; k<params->size; k++){ 
      y[i] += params->c[k]*hag_wavepackets_polynomial_evaluate(x[i], params, k);
    }
    y[i] *= exp(I * params->s / params->eps) *\
            hag_wavepackets_gaussian_evaluate(x[i], params);
  }
}

