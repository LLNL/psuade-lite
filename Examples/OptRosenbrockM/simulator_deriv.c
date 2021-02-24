#include <math.h>
#include <stdio.h>
#include <stdlib.h>
main(int argc, char **argv)
{
   int    n, i;
   double *X, *YY, Y;

   FILE *fIn  = fopen(argv[1], "r");
   FILE *fOut;
   if (fIn == NULL)
   {
      printf("Rosenbrock ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   YY = (double *) malloc((n+1) * sizeof(double));
   for (i = 0; i < n; i++) fscanf(fIn, "%lg", &X[i]);
   Y = 0.0;
   for (i = 0; i < n-1; i++)
   {
      Y += (pow(1.0 - X[i],2.0) + 100.0 * pow(X[i+1] - X[i]*X[i],2.0));
      YY[i] = -2.0*(1.0-X[i]) + 200.0*(2.0*pow(X[i],3.0)-2.0*X[i]*X[i+1]);
      if (i > 0) YY[i] += 200.0 * (X[i] - X[i-1] * X[i-1]);
   }
   YY[n-1] = 200.0 * (X[n-1] - X[n-2] * X[n-2]);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, "%24.16e\n", Y);
   for (i = 0; i < n; i++) fprintf(fOut, "%24.16e\n", YY[i]);
   free(X);
   free(YY);
   fclose(fIn);   
   fclose(fOut);   
}

