#ifndef SURFACE_HOPPING_DOT_H
#define SURFACE_HOPPING_DOT_H

#include <stdbool.h>
#include "../Odeint/odeint.h"
#include "../Potential/potential.h"
#include "../../setup.h"

/* ------------ Structures (shared) -------------*/
struct Particle{
    double      pot_old;
    double      pot_curr;
    double      pot_new;
    double      rho_old;
    double      rho_curr;
    double      rho_new;
    bool        state;
    double      x_old[DIM];
    double      x_curr[DIM];
    double      x_new[DIM];
    double      p[DIM];
    double      p_curr[DIM];
    double      pot_grad[DIM];
};

struct Observables{
    unsigned int npart;
    unsigned int dim;
    double mass_up;
    double mass_down;
    double ke_up;
    double ke_down;
    double e_up;
    double e_down;
    double x_up[DIM];
    double x_down[DIM];
    double p_up[DIM];
    double p_down[DIM];
};

struct Hopper{
  double (*func_transition_probability)(struct Particle *part, 
                                        struct Potential *pot,
                                        struct Odeint *odeint);
  void (*func_hop)(struct Particle *part, 
                    struct Hopper *hopper,
                    struct Potential *pot, 
                    struct Odeint *odeint);
};


/* ---------- Allocate memory -------- to be improved */

struct Particle    *sh_particles_create(const unsigned int npart, const unsigned int dim); // allocate mem for n particles
struct Particle    *sh_sample_wigner(int dim, int npart); // allocate mem for n particles


/* ---------- Initialise Wigner samples --------- */

void        sh_wigner_fill(struct Particle *particles, double *q, double *p, double std, 
                          const unsigned int npart, const unsigned int dim);


/* ---------- Update particle data --------- */

void        sh_particle_potential_init(struct Particle *particles, struct Potential *pot,
                                       const unsigned int npart, const unsigned int dim);
void        sh_particle_potential_update(struct Particle *part, struct Potential *pot,
                                         const unsigned int dim);

/* ---------- Construct observables structure--------- */

struct Observables *sh_observables_new(unsigned int npart, unsigned int dim);
void sh_observables_update(struct Observables *observables, struct Particle *particles);


/* ---------- Hopping algorithms ---------
 * Handles the hopping of particles for different implementation. */

struct Hopper *sh_hopper_new(char* transition_name);
void        sh_hopper_hop(struct Particle *part, struct Hopper *hopper,
                          struct Potential *pot, struct Odeint *odeint);    // maybe pass transition rate as pointer to function for this one above


// there are not many transition rates so could use switch principle as with numerical solvers
/* ---------- Transition rates --------- */

double      sh_transition_lzadia(struct Particle *part, struct Potential *pot, struct Odeint *odeint); // probably not in this file then...?
double      sh_transition_lzdia(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa1(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa2(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa3(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa12(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa13(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa23(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_sa123(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_samultid_space(struct Particle *part, struct Potential *pot, struct Odeint *odeint);
double      sh_transition_samultid_time(struct Particle *part, struct Potential *pot, struct Odeint *odeint);

#endif /* SURFACE_HOPPING_DOT_H */

