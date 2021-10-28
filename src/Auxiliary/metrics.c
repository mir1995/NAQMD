#include <math.h>
#include "metrics.h"


double norm_l2(const double *v, unsigned int dim){
  double l = 0;
  for (int i=0; i<dim; i++){
    l += pow(v[i], 2);
  }
  return sqrt(l);
}

