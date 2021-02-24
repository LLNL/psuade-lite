#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs;
   double X[10], Y, T;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = X[0]*X[0] + X[1]*X[1] + X[0]*X[1] - 14*X[0] - 16*X[1] + 
       pow(X[2]-10,2) + 4 * pow(X[3]-5,2) + pow(X[4]-3,2) +
       2*pow(X[5]-1,2) + 5*X[6]*X[6] + 7 * pow(X[7]-11,2) + 
       2*pow(X[8]-10,2) + pow(X[9]-7,2) + 45; 
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   T = 105 - 4*X[0] - 5*X[1] + 3*X[6] - 9*X[7];
   fprintf(fOut, " %24.16e\n", T);
   T = - 10*X[0] + 8*X[1] + 17*X[6] - 2*X[7];
   fprintf(fOut, " %24.16e\n", T);
   T = 8*X[0] - 2*X[1] - 5*X[8] + 2*X[9] + 12;
   fprintf(fOut, " %24.16e\n", T);
   T = -3*pow(X[0]-2,2) - 4*pow(X[1]-3,2) - 2*pow(X[2],2) +
       7*X[3] + 120;
   fprintf(fOut, " %24.16e\n", T);
   T = -5*X[0]*X[0] - 8*X[1] - pow(X[2]-6,2) + 2*X[3] + 40;
   fprintf(fOut, " %24.16e\n", T);
   T = -0.5*pow(X[0]-8,2) - 2*pow(X[1]-4,2) - 3*X[4]*X[4] + X[5] + 30;
   fprintf(fOut, " %24.16e\n", T);
   T = -X[0]*X[0] - 2*pow(X[1]-2,2) + 2*X[0]*X[1] - 14*X[4] + 6*X[5];
   fprintf(fOut, " %24.16e\n", T);
   T = 3*X[0] - 6*X[1] - 12*pow(X[8]-8,2) + 7*X[9];
   fprintf(fOut, " %24.16e\n", T);
   fclose(fOut);   
}

