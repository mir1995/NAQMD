/* test_wigner_sampling.c */
/* 0 copyright/licensing */
/* 1 includes */
#include <stdlib.h> /*supplies malloc(), calloc(), realloc() */
#include <unistd.h> /* EXIT_SUCCESS */
#include <stdio.h> /* printf */
#include <math.h>
#include <time.h> // srand(time(0))
#include "../src/Random/random.h"
/* 2 defines */
/* 3 external declarations */
/* 4 typedefs */
/* 5 global variable declarations */
/* 6 function prototypes */


int main(int argc, char *argv[]){
    /*
     * INITIALISE PARAMETERS
     */
    int npart, K;
    double eps, mean, std, sum, rmse;
    FILE *qfile;

    /*
     * SET PARAMETERS FROM COMMAND LINE
     */
    if (argc == 5){
            npart = atoi(argv[1]);
            eps = atof(argv[2]);
            mean = atof(argv[3]);
            std = atof(argv[4]);
    }
    else {
            npart = 1000;
            eps = 0.1;
            mean = 0;
            std = sqrt(eps/2);
            printf("Using default parameters. \n");
    }

    /*
     *  PRINT SIMULATION PARAMETERS
    */

    printf(" ***************************************************************\n");
    printf("  Testing Box-Muller sampling: parameters  \n");
    printf("---> NPART     :  %d\n",npart);
    printf("---> eps       :  %.4f\n", eps);
    printf("---> mean      :  %.4f\n", mean);
    printf("---> std       :  %.4f\n", std);
    printf(" ***************************************************************\n");

    char filename[200];
    sprintf(filename, "./meaneps%.4fq%dp%d.txt", eps, (int)mean, (int)std); 
    qfile = fopen(filename, "w"); 
    fprintf(qfile, "npart \t RMSE \n");

    srand(1); // INITIALISE SEED 
    K = 200;
    for (int i=1; i<8; i++){
      npart = pow(10,i);
      rmse = 0;
      //double sample[npart]; // possibly do not need this one
      for (int k=0; k<K; k++){
        sum = 0;
        for (int j=0; j<npart; j++){
          sum += rand_gaussian_sample(mean, std);
          //rmse += pow(rand_gaussian_sample(mean, std),2);
        }
        sum /= (double)npart; // sample mean
        rmse += pow(sum, 2);
      }
      // variance of the sampling distribution of the sample for fixed npart -> RMSE
      //fprintf(qfile, "%d \t %f \n", npart, sum);
      fprintf(qfile, "%d \t %f \n", npart, sqrt(rmse/K));
    }
    


    return 0;
}


/* 8 function declarations  https://www.tutorialspoint.com/cprogramming/c_typedef.htm*/

