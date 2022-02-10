# include <math.h>
# include <complex.h>
# include "metrics.h"


double norm_l2(const double *v, unsigned int dim){
  double l = 0;
  for (int i=0; i<dim; i++){
    l += pow(v[i], 2);
  }
  return sqrt(l);
}


double norm_L2(const double *x, const double complex *f, long unsigned int dim){
  double l = 0;
  for (unsigned int i=0; i<dim; i++){
    l += pow(cabs(f[i]),2);
  }
  return l * (x[1] - x[0]);
}


double distance_L2(const double *x, const double complex *g, const double complex *h, long unsigned int dim){
  double l = 0;
  for (unsigned int i=0; i<dim; i++){
    l += pow(cabs(g[i] - h[i]),2);
  }
  return l * (x[1] - x[0]);
}

