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

/*  
 
   Use brute force to check the actual strength
of the input array.

*/

main(argc,argv)
int  argc;
char *argv[];
{
int q, nrow, ncol, **A;
int str;

double work;

OA_parsein( argc,argv, &q, &nrow, &ncol, &A );
printf("\nThe array has %d rows, %d columns and appears\n",nrow,ncol);
printf("to have %d symbols, since the largest symbol is %d.\n",q,q-1);
}
