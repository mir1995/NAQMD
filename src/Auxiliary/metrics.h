#ifndef METRICS_DOT_H
#define METRICS_DOT_H

# include <complex.h>


double norm_l2(const double *v, unsigned int dim);
double norm_L2(const double *x, const double complex *f, long unsigned int dim);
double distance_L2(const double *x, const double complex *g, 
                    const double complex *h, long unsigned int dim);

#endif /* SURFACE_HOPPING_DOT_H */

