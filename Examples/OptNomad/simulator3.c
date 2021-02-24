#include <math.h>
#include <stdio.h>
#include <stdlib.h>
main(int argc, char **argv)
{
   int    n, i;
   double *X, Y, pi=3.1415828;

   FILE  *fIn  = fopen(argv[1], "r");
   FILE  *fOut;
   if (fIn == NULL)
   {
      printf("CRIT3 ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   for (i = 0; i < n; i++) fscanf(fIn, "%lg", &X[i]);
   fOut = fopen(argv[2], "w");

   Y = 0.0;
   for (i = 0; i < n/2; i++) 
      Y += (pow(X[i],2.0)-cos(0.2*pi*X[i])*10+
            pow(X[i+15],2.0)-cos(0.2*pi*X[i+15])*10);
   Y *= 150.0;
   fprintf(fOut, "%24.16e\n", Y);
   free(X);
   fclose(fIn);   
   fclose(fOut);   
}

