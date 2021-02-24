#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    nInputs;
   double X1, X2, X3, X4, X5, Y;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
                                                                                
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   fscanf(fIn, "%lg", &X1);
   fscanf(fIn, "%lg", &X2);
   fscanf(fIn, "%lg", &X3);
   fscanf(fIn, "%lg", &X4);
   fscanf(fIn, "%lg", &X5);
   Y = 1.0 + 2 * X1 + 3.5 * X2 + 4.5 * X3 + 5 * X4;
   Y += 6 * X1 * X1 + 7.5 * X1 * X2 + 8.5 * X1 * X3 + 9 * X1 * X4;
   Y += 10.5 * X2 * X2 + 11 * X2 * X3 + 13 * X2 * X4;
   Y += 24.5 * X3 * X3 + 15 * X3 * X4 + 15 * X4 * X4;
   Y += X5;
   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);
}

