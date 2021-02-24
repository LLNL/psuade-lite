#include <math.h>
#include <stdio.h>
#include <stdlib.h>
main(int argc, char **argv)
{
   int    n, i;
   double *X, Y, Y2, Y3;

   FILE  *fIn  = fopen(argv[1], "r");
   FILE  *fOut;
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   for (i = 0; i < n; i++) fscanf(fIn, "%lg", &X[i]);
   fOut = fopen(argv[2], "w");

   Y = 80.0 * X[0] + 95.0 * X[1];
   Y2 = 10.0 * X[0] + 15.0 * X[1] - 100;
   if (Y2 >= 0) Y2 = 0;
   Y3 = 20.0 * X[0] + 15.0 * X[1] - 160;
   if (Y3 >= 0) Y3 = 0;
   Y = Y - 10000 * (Y2 + Y3);
   fprintf(fOut, "%24.16e\n", Y);
   free(X);
   fclose(fIn);   
   fclose(fOut);   
}

