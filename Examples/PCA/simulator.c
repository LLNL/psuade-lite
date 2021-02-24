#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs;
   double X[12], Y1, Y2, Y3, Y4, Y5, pi2=2*3.14159;
   FILE   *fIn  = fopen(argv[1], "r");
   FILE   *fOut;
   char   lineIn[100];

   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in/out files.\n");
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);

   Y1 = Y2 = Y3 = Y4 = 0.0;
   Y1 = Y1 + sin(X[0] * pi2);
   Y2 = Y2 + sin(X[0] * pi2 * 2);
   Y3 = Y3 + sin(X[0] * pi2 * 4);
   Y4 = Y4 + sin(X[0] * pi2 * 8);
   Y5 = 0.0;

   fOut = fopen(argv[2], "w");
   fprintf(fOut, " %24.16e\n", Y1);
   fprintf(fOut, " %24.16e\n", Y2);
   fprintf(fOut, " %24.16e\n", Y3);
   fprintf(fOut, " %24.16e\n", Y4);
   fprintf(fOut, " %24.16e\n", Y5);
   fclose(fOut);   
}

