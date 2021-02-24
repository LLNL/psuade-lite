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

   Y = X[0]*X[0] + X[1]*X[1] + 2.0*X[2]*X[2] + X[3]*X[3] - 
       5.0*X[0] - 5.0*X[1] - 21.0*X[2] + 7.0*X[3];
   fprintf(fOut, "%24.16e\n", Y);
   Y = -(-X[0]*X[0] - X[1]*X[1] - X[2]*X[2] - X[3]*X[3] - 
         X[0] + X[1] - X[2] + X[3] + 8.0);
   fprintf(fOut, "%24.16e\n", Y);

   Y = -(-X[0]*X[0] - 2.0*X[1]*X[1] - X[2]*X[2] - 2.0*X[3]*X[3] + 
          X[0] + X[3] + 10.0);
   fprintf(fOut, "%24.16e\n", Y);

   Y = -(-2.0*X[0]*X[0] - X[1]*X[1] - X[2]*X[2] - 
          2.0*X[0] + X[1] + X[3] + 5.0);
   fprintf(fOut, "%24.16e\n", Y);

   fprintf(fOut, "%24.16e\n", -X[0]);
   fprintf(fOut, "%24.16e\n", X[0]-4.0);
   fprintf(fOut, "%24.16e\n", -X[1]);
   fprintf(fOut, "%24.16e\n", X[1]-4.0);
   fprintf(fOut, "%24.16e\n", -X[2]);
   fprintf(fOut, "%24.16e\n", X[2]-4.0);
   fprintf(fOut, "%24.16e\n", -X[3]);
   fprintf(fOut, "%24.16e\n", X[3]-4.0);

   free(X);
   fclose(fIn);   
   fclose(fOut);   
}

