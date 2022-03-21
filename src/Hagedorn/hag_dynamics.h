# ifndef HAG_DYNAMICS_DOT_H
# define HAG_DYNAMICS_DOT_H

# include "../Potential/potential.h"
# include "../Odeint/odeint.h"

void hag_dynamics_do_step(struct HagedornWaves *params, struct Potential *pot, 
                          struct Odeint *odeint);

# endif
