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
OA_strength( q,nrow,ncol,A,&str,2 );

if(  str <0  ){
  printf("\nThe array does not even have strength 0, meaning that\n");
  printf("it is not composed of symbols 0 through %d.\n");
}
else
  printf("\nThe array has strength %d and no higher strength.\n",str);
}
