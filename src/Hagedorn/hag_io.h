# ifndef HAG_IO_DOT_H
# define HAG_IO_DOT_H

# include "hagedorn.h"

/*--------- Handles saving of data -------------- */

void hag_write(char *fname, struct HagedornWaves *params);
void hag_read(char *fname, struct HagedornWaves *params);

#endif 

