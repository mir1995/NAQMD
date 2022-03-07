#ifndef POTENTIAL_DOT_H
#define POTENTIAL_DOT_H

# include "complex.h"

struct Potential;
struct Potential{
  double delta; // half gap 
  double alpha; // approx tau
  double eps;
  double gamma;
  double   (*func_rho)(struct Potential *pot, double *x, unsigned int dim);
  double   (*func_potup)(struct Potential *pot, double *x, unsigned int dim);
  double   (*func_potdown)(struct Potential *pot, double *x, unsigned int dim);
  double (*func_trace)(struct Potential *pot, double *x, unsigned int dim);
  double (*func_z)(struct Potential *pot, double *x, unsigned int dim);
  double (*func_v12)(struct Potential *pot, double *x, unsigned int dim);
  void (*func_traced)(struct Potential *pot, double *x, double *grad, unsigned int dim);
  void (*func_zd)(struct Potential *pot, double *x, double *grad, unsigned int dim);
  void (*func_v12d)(struct Potential *pot, double *x, double *grad, unsigned int dim);
  void (*func_tracedd)(struct Potential *pot, double *x, double *hess, unsigned int dim);
  void (*func_zdd)(struct Potential *pot, double *x, double *hess, unsigned int dim);
  void (*func_v12dd)(struct Potential *pot, double *x, double *hess, unsigned int dim);
  void   (*func_gradup)(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
  void   (*func_graddown)(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
  void   (*func_hessup)(struct Potential *pot, double *hess_v, double *x, unsigned int dim);
  void   (*func_hessdown)(struct Potential *pot, double *hess_v, double *x, unsigned int dim);
  double complex (*func_get_tau)(struct Potential *pot);
  //double  param[]; // you'd probably want a dictionary
};

/* -------------- Construct potential from model systems -------- */

struct Potential    *potential_construct(
    double (*func_trace)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_z)(struct Potential *pot, double *x, unsigned int dim),
    double (*func_v12)(struct Potential *pot, double *x, unsigned int dim),
    void (*func_traced)(struct Potential *pot, double *x, double *grad, unsigned int dim),
    void (*func_zd)(struct Potential *pot, double *x, double *grad, unsigned int dim),
    void (*func_v12d)(struct Potential *pot, double *x, double *grad, unsigned int dim),
    void (*func_tracedd)(struct Potential *pot, double *x, double *hess, unsigned int dim),
    void (*func_zdd)(struct Potential *pot, double *x, double *hess, unsigned int dim),
    void (*func_v12dd)(struct Potential *pot, double *x, double *hess, unsigned int dim),
    double complex (*func_get_tau)(struct Potential *pot),
    char* potential_name, double param[]);

double rho(struct Potential *pot, double *x, unsigned int dim);
double v_trace(struct Potential *pot, double *x, unsigned int dim);
double v_z(struct Potential *pot, double *x, unsigned int dim);
double v_v12(struct Potential *pot, double *x, unsigned int dim);
void v_traced(struct Potential *pot, double *x, double *grad, unsigned int dim);
void v_zd(struct Potential *pot, double *x, double *grad, unsigned int dim);
void v_v12d(struct Potential *pot, double *x, double *grad, unsigned int dim);
void v_tracedd(struct Potential *pot, double *x, double *hess, unsigned int dim);
void v_zdd(struct Potential *pot, double *x, double *hess, unsigned int dim);
void v_v12dd(struct Potential *pot, double *x, double *hess, unsigned int dim);
double v_up(struct Potential *pot, double *x, unsigned int dim);
double v_down(struct Potential *pot, double *x, unsigned int dim);
double complex get_tau(struct Potential *pot);
void grad_v_up(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
void grad_v_down(struct Potential *pot, double *grad_v, double *x, unsigned int dim);
void hess_v_up(struct Potential *pot, double *hess_v, double *x, unsigned int dim);
void hess_v_down(struct Potential *pot, double *hess_v, double *x, unsigned int dim);

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
