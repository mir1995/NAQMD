# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include "hagedorn.h"


/*--------- Handles saving of data -------------- */
void hag_write(char *fname, struct HagedornWaves *params){
  // open fname to write 
  
  FILE *file;
  char fname_ab[100];
  
  strcpy(fname_ab, fname);
  strcat(fname_ab, ".dat");
  file = fopen(fname_ab, "ab");
  
  fwrite(params, sizeof(*(params)), 1, file);
  fwrite(params->c, sizeof(*(params->c)), params->size, file);
  fclose(file);
}

void hag_read(char *fname, struct HagedornWaves *params){
  
  FILE *file;
  char fname_rb[100];
  
  strcpy(fname_rb, fname);
  strcat(fname_rb, ".dat");
  file = fopen(fname_rb, "rb");
  
  fread(params, sizeof(*(params)), 1, file); 
  params->c = (double complex *)malloc(params->size * sizeof(double complex));
  
  fread(params->c, sizeof(*(params->c)), params->size , file);
  fclose(file);
}
