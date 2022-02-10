# ifndef POTENTIAL_IO_DOT_H
# define POTENTIAL_IO_DOT_H


# include "potential.h"


/*--------- Handles saving of data -------------- */
void pot_write(char *fname, struct Potential *pot);
void pot_read(char *fname, struct Potential *pot);

# endif
