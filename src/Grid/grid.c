# include "grid.h"


/* Handles creation of a uniform grid */
void grid_uniform_fill(double x[], double x_left, double x_right, long unsigned int n){
  double step = (x_right - x_left) / (1.0 * n);
  x[0] = x_left;
  for (long unsigned int i=1; i<n+1; i++){
    x[i] = x[i-1] + step;
  }
}
