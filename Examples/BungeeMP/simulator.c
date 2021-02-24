#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main(int argc, char **argv)
{
   int    i, nInputs, count, id;
   double X[3], Y, kel=1.5, g=9.8;
   char   pString[500];
   FILE   *fIn;
   FILE   *fOut;

   sscanf(argv[1], "%d", &id);
   sprintf(pString, "psuadeApps_ct.in.%d", id);
   fIn = fopen(pString, "r");
   if (fIn == NULL)
   {
      printf("Simulator ERROR - cannot open in file %s.\n",pString);
      exit(1);
   }
   fscanf(fIn, "%d", &nInputs);
   for (i = 0; i < nInputs; i++) fscanf(fIn, "%lg", &X[i]);
   Y = X[0] - 2.0 * X[1] * g / (kel * X[2]);
   system("sleep 60");
   sprintf(pString, "psuadeApps_ct.out.%d", id);
   fOut = fopen(pString, "w");
   fprintf(fOut, " %24.16e\n", Y);
   fclose(fOut);   
}

