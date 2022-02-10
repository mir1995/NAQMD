#ifndef HAGEDORN_PROJECTION_DOT_H
#define HAGEDORN_PROJECTION_DOT_H

#include <complex.h>
#include "hagedorn.h"
#include "../Potential/potential.h" 


// now you are passing f(x) but perhaps you can even Construct it 
// inside the projection function ?
struct HagedornWaves *hag_project(struct HagedornWaves *params_in,
                                  unsigned int n, unsigned int K, struct Potential *pot);


#endif 

