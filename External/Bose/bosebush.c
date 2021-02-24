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
  ncol = 2*q;    /* Can get 2q without defect, 2q+1 with defect */
}else{
  sscanf(argv[1],"%d",&q);
  sscanf(argv[2],"%d",&ncol);
}

if(  q%2  ){
  fprintf(stderr,"This implementation of Bose-Bush only works for a number\n");
  fprintf(stderr,"q of levels equal to a power of 2.  q=%d was requested.\n",q);
  fprintf(stderr,"The Addelman-Kempthorne designs might be available.\n");
  exit(1);
}

if(  !GF_getfield(2*q, &gf)  ){
  fprintf(stderr,"Could not construct the Galois field needed\n");
  fprintf(stderr,"for the Bose Bush design.\n");
  exit(1);
}

A = imatrix( 0,2*q*q-1, 0,ncol-1  );
if(  !A  ){
  fprintf(stderr,"Could not allocate array for Bose design.\n");
  exit(1);
}  

if(  bosebush( &gf, A, ncol )  ){
  OA_put( A,2*q*q,ncol,q );
  exit(0);
}
else{
  fprintf(stderr,"Unable to construct Bose design q=%d, ncol=%d.\n",
	  q,ncol);
  exit(1);
}
}

