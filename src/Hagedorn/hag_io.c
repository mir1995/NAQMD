# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include "hagedorn.h"


/*--------- Handles saving of data -------------- */
void hag_write(char *fname, struct HagedornWaves *params){
  // open fname to write 
  
  FILE *f;
  char path[100];
  size_t i;
  
  strcpy(path, fname);
  strcat(path, ".txt");
  f = fopen(path, "a+");
 

  fseek (f, 0, SEEK_END);
  long int size = ftell(f);

  if (0 == size) {
    printf("file is empty\n");
    /* column names */
    fprintf(f, "DIM\tSIZE\tSTATE\tEPS\tS");
    for (i=0; i<params->dim; i++){
      fprintf(f, "\tq[%ld]", i);
    }
    for (i=0; i<params->dim; i++){
      fprintf(f, "\tp[%ld]", i);
    }
    for (i=0; i<pow(params->dim,2); i++){
      fprintf(f, "\tQ[%ld].re\tQ[%ld].im", i, i);
    }
    for (i=0; i<pow(params->dim,2); i++){
      fprintf(f, "\tP[%ld].re\tP[%ld].im", i, i);
    }
    for (i=0; i<params->size; i++){
      fprintf(f, "\tc[%ld].re\tc[%ld].im", i, i);
    }
    fprintf(f, "\n");
  }

  // write some metadata
  fprintf(f, "%d \t %d \t %d \t %.17g \t %.17g ",
              params->dim,
              params->size,
              params->state, 
              params->eps,
              params->s);
  for (i=0; i<params->dim; i++){
    fprintf(f, "\t %.17g ", params->q[i]);
  }
  for (i=0; i<params->dim; i++){
    fprintf(f, "\t %.17g ", params->p[i]);
  }
  for (i=0; i<pow(params->dim,2); i++){
    fprintf(f, "\t %.17g \t %.17g", creal(params->Q[i]), cimag(params->Q[i]));
  }
  for (i=0; i<pow(params->dim,2); i++){
    fprintf(f, "\t %.17g \t %.17g", creal(params->P[i]), cimag(params->P[i]));
  }
  for (i=0; i<params->size; i++){
    fprintf(f, "\t %.17g \t %.17g", creal(params->c[i]), cimag(params->c[i]));
  }
  fprintf(f, "\n");
  //fwrite(params, sizeof(*(params)), 1, file);
  //fwrite(params->c, sizeof(*(params->c)), params->size, file);
  fclose(f);
}

void hag_read(char *fname, struct HagedornWaves *params){
  
  FILE *file;
  char fname_rb[100];
  
  /*
  strcpy(fname_rb, fname);
  strcat(fname_rb, ".dat");
  file = fopen(fname_rb, "rb");
  
  fread(params, sizeof(*(params)), 1, file); 
  params->c = (double complex *)malloc(params->size * sizeof(double complex));
  
  fread(params->c, sizeof(*(params->c)), params->size , file);
  fclose(file);
  */
}
