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
int       str, q, ncol, **A;
struct GF gf;

if(  argc==1  )
  scanf("%d %d %d",&str,&q,&ncol);
else if( argc==2  ){
  sscanf(argv[1],"%d",&str);
  scanf("%d %d",&q,&ncol);
}else if( argc==3  ){
  sscanf(argv[1],"%d",&str);
  sscanf(argv[2],"%d",&q);
  ncol = q+1;  /* can sometimes get 1 more with even q */
}else{
  sscanf(argv[1],"%d",&str);
  sscanf(argv[2],"%d",&q);
  sscanf(argv[3],"%d",&ncol);
}

if( str < 2 ){
  fprintf(stderr,"Bush designs not provided for strength %d.\n",str);
  exit(1);
}

if( ncol > q+1 ){
  fprintf(stderr,"Only q+1 = %d columns provided in Bush design.\n",q+1);
  fprintf(stderr,"Columns requested was %d.\n",ncol);
  exit(1);
}

if(  !GF_getfield(q, &gf)  ){
  fprintf(stderr,"Could not construct the Galois field needed\n");
  fprintf(stderr,"for the strength %d Bush design on %d levels.\n",str,q);
  exit(1);
}

A = imatrix( 0,ipow(q,str)-1, 0,ncol-1  );
if(  !A  ){
  fprintf(stderr,"Could not allocate array for Bush design.\n");
  exit(1);
}  

if(  bush( &gf, A, str, ncol )  ){
  OA_put( A,ipow(q,str),ncol,q );
  exit(0);
}
else{
  fprintf(stderr,"Unable to construct the strength %d Bush design q=%d, ncol=%d.\n",str, q,ncol);
  exit(1);
}
}

