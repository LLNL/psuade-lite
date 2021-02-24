#include <math.h>
#include <stdio.h>
#include <stdlib.h>
main(int argc, char **argv)
{
   int    n, i;
   double *X, Y;

   FILE  *fIn  = fopen(argv[1], "r");
   FILE  *fOut;
   if (fIn == NULL)
   {
      printf("Asaadi 1 ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d",  &n);
   X = (double *) malloc(n * sizeof(double));
   for (i = 0; i < n; i++) fscanf(fIn, "%lg", &X[i]);
   fOut = fopen(argv[2], "w");

   Y = pow(X[0]-10.0, 2.0) + 5.0 * pow(X[1]-12.0, 2.0) +
       pow(X[2], 3.0) + 3.0 * pow(X[3]-11.0, 2.0) +
       10.0 * pow(X[4], 6.0) + 7.0 * pow(X[5], 2.0) +
       pow(X[6], 4.0) - 4.0 * X[5] * X[6] - 10.0 * X[5] - 8.0 * X[6];
   fprintf(fOut, "%24.16e\n", Y);

   Y = - (-2.0 * X[0] * X[0] - 3.0 * pow(X[1], 4.0) - X[2] - 
          4.0 * X[3] * X[3] - 5.0 * X[4] + 127.0);
   fprintf(fOut, "%24.16e\n", Y);

   Y = - (7.0 * X[0] - 3.0 * X[1] - 10.0 * X[2] * X[2] - 
          X[3] + X[4] + 282.0);
   fprintf(fOut, "%24.16e\n", Y);

   Y = - (23.0 * X[0] - X[1] * X[1] - 6.0 * X[5] * X[5] + 
          8.0 * X[6] + 196.0);
   fprintf(fOut, "%24.16e\n", Y);

   Y = - (-4.0 * X[0] * X[0] - X[1] * X[1] + 3.0 * X[0] * X[1] - 
          2.0 * X[2] * X[2] - 5.0 * X[5] + 11.0 * X[6]);
   fprintf(fOut, "%24.16e\n", Y);

   free(X);
   fclose(fIn);   
   fclose(fOut);   
}

