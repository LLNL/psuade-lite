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


#include <math.h>
#include <stdio.h>
#include "galois.h"

main(argc,argv)
int  argc;
char *argv[];
{
int       q, ncol, **A;
struct GF gf;

if(  argc==1  )
  scanf("%d %d",&q,&ncol);
else if( argc==2  ){
  sscanf(argv[1],"%d",&q);
  ncol = 2*q*q + 2*q +1;  /*  2(q^3-1)/(q-1) - 1  */
}else{
  sscanf(argv[1],"%d",&q);
  sscanf(argv[2],"%d",&ncol);
}

if(  !GF_getfield(q, &gf)  ){
  fprintf(stderr,"Could not construct the Galois field needed\n");
  fprintf(stderr,"for the Addelman-Kempthorne (n=3) design.\n");
  exit(1);
}

A = imatrix( 0, 2*q*q*q-1, 0, ncol-1  );
if(  !A  ){
  fprintf(stderr,"Could not allocate array for Addelman-Kempthorne (n=3) design.\n");
  exit(1);
}  

if(  addelkemp3( &gf, A, ncol )  ){
  OA_put( A, 2*q*q*q, ncol, q );
  exit(0);
}
else{
  fprintf(stderr,"Unable to construct Addelman-Kempthorne (n=3) design q=%d, ncol=%d.\n",
	  q,ncol);
  exit(1);
}
}

