# ifndef HAG_SUPERADIABATIC_DOT_H
# define HAG_SUPERADIABATIC_DOT_H

# include "../Potential/potential.h" 
# include <complex.h>
# include <string.h>

void hag_superadiabatic_formula(double *x, double complex *f,
                                struct HagedornWaves *params_in, long unsigned int n,
                                struct Potential *pot, char *version);

#endif 

