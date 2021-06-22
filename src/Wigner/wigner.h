#ifndef WIGNER_DOT_H    /* This is an "include guard" */
#define WIGNER_DOT_H    /* prevents the file from being included twice. */
                               /* Including a header file twice causes all kinds */
                               /* of interesting problems.*/
// https://riptutorial.com/c/example/3250/calling-a-function-from-another-c-file 
/**
 * This is a function declaration.
 * It tells the compiler that the function exists somewhere.
 */
double gaussian_get_sample();
void wigner_set_samples(ParticleT *ptr, int npart, int dim);

#endif /* WIGNER_DOT_H */
