# ifndef PARTITION_UNITY_DOT_H
# define PARTITION_UNITY_DOT_H

# include <complex.h>
# include <stdbool.h>
# include "../../setup.h"
#include "../Hagedorn/hagedorn.h"

struct HagedornWaves **hag_projection_partition_unity(struct HagedornWaves *params, unsigned int *F);


#endif 
