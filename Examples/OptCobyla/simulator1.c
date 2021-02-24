#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs;
   double X[2], Y;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = X[0] * X[0] + X[1] * X[1];
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fprintf(fOut, " %24.16e\n", 1-X[0]-X[1]);
   fprintf(fOut, " %24.16e\n", 1-X[0]*X[0]-X[1]*X[1]);
   fprintf(fOut, " %24.16e\n", 9-9*X[0]*X[0]-X[1]*X[1]);
   fprintf(fOut, " %24.16e\n", X[1]-X[0]*X[0]);
   fprintf(fOut, " %24.16e\n", X[0]-X[1]*X[1]);
   fclose(fOut);   
}

