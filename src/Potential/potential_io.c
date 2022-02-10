# include <stdio.h>
# include <string.h>
# include "potential.h"


/*--------- Handles saving of data -------------- */
void pot_write(char *fname, struct Potential *pot){
  // open fname to write 
  
  FILE *file;
  char fname_ab[100];
  
  strcpy(fname_ab, fname);
  strcat(fname_ab, ".dat");
  file = fopen(fname_ab, "ab");
  
  fwrite(pot, sizeof(*(pot)), 1, file);
  fclose(file);
// to read 
// f=fopen(fname, "rb");
// fread(&struct, sizeof(strtct), 1, f)
//fclose(fname);


}
void pot_read(char *fname, struct Potential *pot){
  
  FILE *file;
  char fname_rb[100];
  
  strcpy(fname_rb, fname);
  strcat(fname_rb, ".dat");
  file = fopen(fname_rb, "rb");
  
  fread(pot, sizeof(*(pot)), 1, file); 
  fclose(file);
}

