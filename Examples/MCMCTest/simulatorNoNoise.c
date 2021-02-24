#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, j, nInputs;
   double X[25], Y[25], Z[25], A = 0.25, B = 0.75, F, AB[2], error;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
                                                                                
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &AB[i]);
   for (i = 0; i < 5; i++)
   {
      for (j = 0; j < 5; j++)
      {
         X[i*5+j] = i * 0.25;
         Y[i*5+j] = j * 0.25;
      }
   }
   for (i = 0; i < 25; i++)
      Z[i] = A * X[i] * X[i] + B * Y[i] * Y[i];
                                                                                
   F = 0.0;
   for (i = 0; i < 25; i++)
   {
      error = Z[i] - AB[0]*X[i]*X[i] - AB[1]*Y[i]*Y[i];
      F += error * error;
   }
   F = F / 25.0;
   F = sqrt(F);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", F);
   fclose(fOut);
}

