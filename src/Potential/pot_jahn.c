// should I assume the surfaces are analytically given 
// o.w. need to evaluate gradient numerically ?
#include <math.h>
#include "potential.h"


// perhaps should extend a struct whereby you fix the different 
// values as function pointers such as...
void h_jahn_teller_v1(double *pt, double* pos, double *delta){
    *pt = sqrt(pow(pos[0],2) + pow(pos[1],2) + pow(*delta,2));
}

void set_ptdwn(double *pt, double* pos, double *delta){
    *pt = - sqrt(pow(pos[0],2) + pow(pos[1],2) + pow(*delta,2));
}

void set_grdup(double* grdbfr, double* pos, double *delta){
    double den = sqrt(pow(pos[0],2) + pow(pos[1],2) + pow(*delta,2)); 
    grdbfr[0] = pos[0] / den;
    grdbfr[1] = pos[1] / den;
}

void set_grddwn(double* grdbfr, double* pos, double *delta){
    double den = sqrt(pow(pos[0],2) + pow(pos[1],2) + pow(*delta,2)); 
    grdbfr[0] = - pos[0] / den;
    grdbfr[1] = - pos[1] / den;
}

