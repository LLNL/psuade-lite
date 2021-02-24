#include <stdio.h>

main()
{
   int n = 279936, i, j, ncnt, *part;
   FILE *fp;

   part = (int *) malloc(n*sizeof(int));
   fp = fopen("new.mgraph.part.270", "r"); 
   for ( i = 0; i < n; i++ ) fscanf(fp,"%d", &part[i]);
   fclose(fp);
   for ( i = 0; i < 270; i++ )
   {
      ncnt = 0;
      for ( j = 0; j < n; j++ ) if (part[j] == i) ncnt++;
      printf("part %3d has %5d nodes\n", i, ncnt);
   }
}
   
      
