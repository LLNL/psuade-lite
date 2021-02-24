#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int    i, nInputs, count;
   double X[2], Y, pi=3.1415928;
   FILE   *fIn;
   FILE   *fOut;

   fIn  = fopen(argv[1], "r");
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = 0;
   if (X[0] > 0.8) 
      Y = Y + sin(5*(X[0]-0.8)*2*pi);
   if (X[1] > 0.8)
      Y = Y + sin(5*(X[1]-0.8)*2*pi);
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

