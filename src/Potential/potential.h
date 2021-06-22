#ifndef POTENTIAL_DOT_H
#define POTENTIAL_DOT_H

#include "../SurfaceHopping/surface_hopping.h"

struct Potential;
struct Potential{
  double delta; 
  double alpha;
  double eps;
  double   (*func_potup)(struct Potential *pot, double *x, unsigned int dim);
  double   (*func_potdown)(struct Potential *pot, double *x, unsigned int dim);
  void   (*func_gradup)(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
  void   (*func_graddown)(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
  void   (*func_dd_up)(struct Potential *pot, double *dd_v, double *x, unsigned int dim);
  void   (*func_dd_down)(struct Potential *pot, double *dd_v, double *x, unsigned int dim);
  double (*func_get_tau)(struct Potential *pot);
  //double  param[]; // you'd probably want a dictionary
};

/* -------------- Construct potential from model systems -------- */

struct Potential    *potential_construct(
    double (*func_potup)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_potdown)(struct Potential *pot, double *x, unsigned int dim),
    void (*func_gradup)(struct Potential *pot, double *grad_v, double *x, unsigned int dim),
    void (*func_graddown)(struct Potential *pot, double *grad_v, double *x, unsigned int dim),
    void (*func_dd_up)(struct Potential *pot, double *dd_v, double *x, unsigned int dim),
    void (*func_dd_down)(struct Potential *pot, double *dd_v, double *x, unsigned int dim),
    double (*func_get_tau)(struct Potential *pot),
    char* potential_name, double param[]);

double v_up(struct Potential *pot, double *x, unsigned int dim);
double v_down(struct Potential *pot, double *x, unsigned int dim);
double get_tau(struct Potential *pot);
void grad_v_up(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
void grad_v_down(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
void dd_v_up(struct Potential *pot, double *dd_v, double *x, unsigned int dim);
void dd_v_down(struct Potential *pot, double *dd_v, double *x, unsigned int dim);

#endif

/*--------------- Secod approach - you might want to think about the implications of the two ------------ */
/*
struct Potential{
  double delta; 
  double alpha;
  void   (*func_potup)(struct Potential *pot, double *x, unsigned int dim);
  void   (*func_potdown)(struct Potential *pot, double *x, unsigned int dim);
  void   (*func_gradup)(struct Potential *pot, double *x, unsigned int dim);
  void   (*func_graddown)(struct Potential *pot, double *x, unsigned int dim);
  //double  param[]; // you'd probably want a dictionary
};
*/
