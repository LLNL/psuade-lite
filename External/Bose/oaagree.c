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

/* Count Agreements among rows of an array */

#define ROWCHECK 50

int **A;

main(argc,argv)
int  argc;
char *argv[];
{
int q, nrow, ncol, **A;
int i, j, k;
int agree, maxagr;
int mrow1, mrow2;

OA_parsein( argc,argv, &q, &nrow, &ncol, &A );

maxagr = mrow1 = mrow2 = 0;

for( i=0; i<nrow; i++  ){
  for( j=i+1; j<nrow; j++  ){
    agree = 0;
    for( k=0; k<ncol; k++  )
      agree += (A[i][k]==A[j][k]);
    if(  agree>maxagr  ){
      maxagr = agree;
      mrow1 = i;
      mrow2 = j;
      printf("New max %d %d %d\n",i,j,agree);
    }
  }
  if(  i && i % ROWCHECK == 0  )
    printf("Checked rows <= %d vs all other rows.\n",i);
}
if(  maxagr == 0  )
  printf("No two distinct rows agree in any columns.\n");
else{
  printf("Maximum number of columns matching for two distinct rows is %d.\n",
	 maxagr);
  printf("This is attained by rows %d and %d.\n",mrow1,mrow2);
}

}
