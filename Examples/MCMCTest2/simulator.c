#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    nInputs;
   double X1, X2, X3, X4, Y;
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
   Y = 1.0 + 2 * X1 + 3.5 * X2 + 3.5 * X3 + 5 * X4;
   Y += 6 * X1 * X1 + 7.5 * X1 * X2 + 7.5 * X1 * X3 + 9 * X1 * X4;
   Y += 10.5 * X2 * X2 + 11 * X2 * X3 + 13 * X2 * X4;
   Y += 10.5 * X3 * X3 + 13 * X3 * X4 + 15 * X4 * X4;

   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y);
   Y = 10 + 11 * X1 + 12.5 * X2 + 12.5 * X3 + 15 * X4;
   Y += 14 * X1 * X1 + 7 * X1 * X2 + 7 * X1 * X3 + 3 * X1 * X4;
   Y += 4 * X2 * X2 + 8 * X2 * X3 + 5 * X2 * X4;
   Y += 4 * X3 * X3 + 5 * X3 * X4 + 2 * X4 * X4;
   fprintf(fOut, " %24.16e\n", Y);

   fclose(fOut);
}

