#include <stdio.h>

main()
{
   int dim = 5, n1d = 4, ii, i1, i2, jj;
   int npts = 270, incrs[11], *ia, *conn, nnz, N;
   FILE *fp;

   N = 1;
   incrs[0] = N;
   for ( i1 = 1; i1 <= dim; i1++ ) {N *= n1d; incrs[i1] = N;}
   printf("N = %d\n", N);

   ia   = (int *) malloc((N+1)*sizeof(int));
   conn = (int *) malloc((N*dim*2)*sizeof(int));
   nnz = 0;
   ia[0] = nnz;
   for ( i1 = 0; i1 < N; i1++ )
   {
      ii = i1;
      for ( i2 = 0; i2 < dim; i2++ )
      {
         jj = ii % n1d;
         ii = ii / n1d;
         if ( jj > 0     ) conn[nnz++] = i1 - incrs[i2] + 1;
         if ( jj < n1d-1 ) conn[nnz++] = i1 + incrs[i2] + 1;
      }
      ia[i1+1] = nnz;
   }
   fp = fopen("new.mgraph", "w"); 
   fprintf(fp,"%d %d 0 0\n", N, nnz/2); 
   for ( i1 = 0; i1 < N; i1++ )
   {
      for ( i2 = ia[i1]; i2 < ia[i1+1]; i2++ )
         fprintf(fp, "%d ", conn[i2]);
      fprintf(fp, "\n");
   }
   fclose(fp);
   printf("nnz = %d\n", nnz);
}
   
      
