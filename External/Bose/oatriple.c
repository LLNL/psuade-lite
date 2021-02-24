/*

  These programs construct and manipulate orthogonal 
arrays.  They were prepared by

    Art Owen
    Department of Statistics
    Sequoia Hall
    Stanford CA 94305

  They may be freely used and shared.  This code comes
with no warranty of any kind.  Use it at your own
risk.

  I thank the Semiconductor Research Corporation and
the National Science Foundation for supporting this
work.

*/


#include <stdio.h>
#include "oa.h"
int *ivector(), **imatrix();

/* Count triple agreements among rows of an array */

main(argc,argv)
int  argc;
char *argv[];
{
int a3, q, nrow, ncol, **A;
int i1, i2, j1, j2, j3;
int num3;

OA_parsein( argc,argv, &q, &nrow, &ncol, &A );

num3 = 0;
for(  j1=0;    j1<ncol; j1++  )
for(  j2=j1+1; j2<ncol; j2++  )
for(  j3=j2+1; j3<ncol; j3++  ){
  a3 = 0;
  for( i1=0;    i1<nrow; i1++  )
  for( i2=i1+1; i2<nrow; i2++  )
    a3 += ( A[i1][j1]==A[i2][j1] )&&( A[i1][j2]==A[i2][j2] )&&( A[i1][j3]==A[i2][j3] );
  if( a3 ){
    printf("Cols %d %d %d match in %d distinct pairs of rows.\n",j1,j2,j3,a3 );
    num3++;
  }
} 
printf("There are %d distinct triples of columns that agree\n",num3);
printf("in at least two distinct rows.\n");
}
