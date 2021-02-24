#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

main(int argc, char **argv)
{
   int    nInputs;
   double Y, X;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
                                                                                
   if (fIn == NULL)
   {
      printf("Simulator2 ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   fscanf(fIn, "%lg", &X);
   fclose(fIn);
   Y = X + X * X;
   fOut  = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);
}

