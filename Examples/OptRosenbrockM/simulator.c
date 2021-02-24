#include <math.h>
#include <stdio.h>
#include <stdlib.h>
main(int argc, char **argv)
{
   int    n, i;
   double *X, Y;
   char   lineIn[100], stringPtr[100], equal[2];

   FILE  *fIn  = fopen(argv[1], "r");
   FILE  *fOut;
   if (fIn == NULL)
   {
      printf("RosenbrockM ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   for (i = 0; i < n; i++) fscanf(fIn, "%lg", &X[i]);
   Y = 0.0;
   for (i = 0; i < n-1; i++)
      Y += (pow(1.0 - X[i],2.0) + 100.0 * pow(X[i+1] - X[i]*X[i],2.0));
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   free(X);
   fclose(fIn);   
   fclose(fOut);   
}

