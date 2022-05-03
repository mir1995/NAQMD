# ifndef HAGEDORN_DOT_H
# define HAGEDORN_DOT_H

# include <complex.h>
# include <stdbool.h>
# include "../../setup.h"

struct HagedornWaves{
  unsigned int   dim;
  unsigned int   size; //  number of coefficients
  bool           state; // upper or lower level
  double eps;
  double s;
  double q[DIM]; // i am not yet sure whether we need it to be complex
  double p[DIM];
  double complex Q[DIM * DIM];
  double complex P[DIM * DIM];
  double complex c[]; // coefficients of Hagedorn wavepackets
};
typedef struct HagedornWaves HagedornWaves;


/*--------- Handles creation of a HagedornWaves struct */
struct HagedornWaves *hag_wavepacket_create(unsigned int dim, unsigned int size, 
                                            bool state, double eps, double s, 
                                            double q[], double p[], double complex Q[],
                                            double complex P[], double complex c[]);


/*--------- Handles generation of Hagedorn wavepackets -------------- */

double complex hag_wavepackets_gaussian_evaluate(
    double complex x, struct HagedornWaves *params);

double complex hag_wavepackets_polynomial_evaluate(double complex x,
    struct HagedornWaves *params, int index);

void hag_wavepackets_gaussian_fill(
  double complex x[], double complex y[], 
  struct HagedornWaves *params, int n);

void hag_wavepackets_polynomial_fill(
    double complex x[], double complex y[], 
    struct HagedornWaves *params, int index, int n);

void hag_wavepackets_fill(
    double x[], double complex y[], 
    struct HagedornWaves *params, long unsigned int n);

/*--------- Handles projection of transmitted wavepacket onto Hagedorn basis ---------- */

// i have a separate .h file for this function now
void hag_projection_get_coefficients(double complex c[] ,int N);



#endif 
