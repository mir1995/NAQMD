/* test_gaussian.c */
/* 0 copyright/licensing */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <stdio.h> /* printf */
#include <time.h> // srand(time(0)) 
#include "../sample_wigner.h"
/* 3 external declarations */
/* 4 typedefs */
/* 5 global variable declarations */
/* 6 function prototypes */
int main(int argc, char *argv[]){
    int n;
    if (argc == 2){
            n = atoi(argv[1]);
            printf("%d samples from a Gaussian distribution \n", n);
    }
    else {
            printf("Error: one argument (# samples) expected \n");
            return 1;
    }
    
    srand(time(0));
    double *pos = malloc(n * sizeof(double));
    double *mom = malloc(n * sizeof(double));
    
    // standard gaussian samples
    sample_wigner(pos, mom, n);
    FILE *qfile = fopen("./position.txt", "w");
    for (int i=0; i<n; i++){
            fprintf(qfile, "%.12f \n", pos[i]*20 + 100);
    }
    fclose(qfile);
    return 0;
}
// gcc -Wall -c foo.c , gcc -Wall -c main.c  - link them ? --> gcc -o name foo.o
// main.o
/* 8 function declarations  https://www.tutorialspoint.com/cprogramming/c_typedef.htm*/
