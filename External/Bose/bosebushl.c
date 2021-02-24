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
int       lam, q, ncol, **A;
struct GF gf;
int       pq, nq, isppq,  pl, nl, isppl;

if(  argc==1  )
  scanf("%d %d %d",&lam,&q,&ncol);
else if( argc==2  ){
  sscanf(argv[1],"%d",&lam);
  scanf("%d %d",&q,&ncol);
}else if( argc==3 ){
  sscanf(argv[1],"%d",&lam);
  sscanf(argv[2],"%d",&q);
  ncol = q*lam;
}else{
  sscanf(argv[1],"%d",&lam);
  sscanf(argv[2],"%d",&q);
  sscanf(argv[3],"%d",&ncol);
}

primepow( lam, &pl, &nl, &isppl  );
primepow(  q , &pq, &nq, &isppq  );

if(  !isppq  ){
  fprintf(stderr,"The Bose-Bush construction requires that q be a prime\n");
  fprintf(stderr,"raised to a positive integral power. q=%d was requested\n",q);
  fprintf(stderr," and is not such a prime power.\n");
  exit(1);
}

if(  !isppl  ){
  fprintf(stderr,"The Bose-Bush construction requires that lambda be a prime\n");
  fprintf(stderr,"raised to a positive integral power. lambda=%d was requested\n",lam);
  fprintf(stderr," and is not such a prime power.\n");
  exit(1);
}

if(  pl != pq  ){
  fprintf(stderr,"The Bose-Bush construction requires lambda and q\n");
  fprintf(stderr,"to be powers of the same prime. So lambda = %d = %d^%d\n",
	  lam,pl,nl);
  fprintf(stderr,"and q = %d = %d^%d are not suitable.\n",q,pq,nq);
  exit(1);
}

if(  !GF_getfield(lam*q, &gf)  ){
  fprintf(stderr,"Could not construct the Galois field needed\n");
  fprintf(stderr,"for the Bose-Bush design.\n");
  exit(1);
}

A = imatrix( 0,lam*q*q-1, 0,ncol-1  );
if(  !A  ){
  fprintf(stderr,"Could not allocate array for Bose design.\n");
  exit(1);
}  

if(  bosebushl( &gf, lam, A, ncol )  ){
  OA_put( A,lam*q*q,ncol,q );
  exit(0);
}
else{
  fprintf(stderr,"Unable to construct Bose design lambda=%d, q=%d, ncol=%d.\n",
	  lam,q,ncol);
  exit(1);
}
}



