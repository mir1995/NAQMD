#ifndef HAG_OBSERVABLES_DOT_H
#define HAG_OBSERVABLES_DOT_H

#include <complex.h>
#include "hagedorn.h"
#include "../Potential/potential.h" 


double complex hag_observables_get_mean(double *x, double *w, double complex *f,  
                                        struct HagedornWaves *params_in,
                                        struct Potential *pot,
                                        unsigned int n);

#endif 


