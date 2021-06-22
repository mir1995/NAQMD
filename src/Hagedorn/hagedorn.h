#ifndef HAGEDORN_DOT_H
#define HAGEDORN_DOT_H

#include <complex.h>

struct HagedornParameters{
  double complex q;
  double complex p;
  double complex Q;
  double complex P;
  double complex eps;
};
typedef struct HagedornParameters HagedornParameters;

/*--------- Handles generation of Hagedorn wavepackets -------------- */

double complex hag_wavepackets_gaussian_evaluate(
    double complex x, HagedornParameters params);

double complex hag_wavepackets_polynomial_evaluate(double complex x,
    HagedornParameters params, int index);

void hag_wavepackets_gaussian_fill(
  double complex x[], double complex y[], 
  HagedornParameters params, int n);

void hag_wavepackets_polynomial_fill(
    double complex x[], double complex y[], 
    HagedornParameters params, int index, int n);

/*--------- Handles projection of transmitted wavepacket onto Hagedorn basis ---------- */

void hag_projection_get_coefficients(double complex c[] ,int N);



#endif 
