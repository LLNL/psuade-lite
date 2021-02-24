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


/*     Addelman-Kempthorne Constructions with general n.

       The article is quite vague on this.  Page 1173 states
   "When n>2 the same procedure will yield the desired plans
   if Lemma 5a is used in place of Lemma 5."  Page 1175
   provides the example n=3,q=3 which is OA( 54,25,3,2 ).
   Based on this example it is possible to make an educated
   guess as to how the construction generalizes.

*/

#include <math.h>
#include <stdio.h>

#include "galois.h"

addelkempncheck( q,p,akn,ncol  )
int q,p,akn,ncol;
{

if(  akn<2  ){
  fprintf(stderr,"This Addelman-Kempthorne OA(2q^n,ncol,q,2) is only\n");
  fprintf(stderr,"available for n >= 2.  n = %d was requested.\n",akn);
  return 0;
}

if(  p==2 && q>4 ){
  fprintf(stderr,"This Addelman-Kempthorne OA(2q^n,ncol,q,2) is only\n");
  fprintf(stderr,"available for odd prime powers q and for even prime\n");
  fprintf(stderr,"powers q<=4.\n");
  return 0;
}

if(  ncol > 2*(ipow(q,akn)-1)/(q-1) -1){
  fprintf(stderr,"The Addelman-Kempthorne construction needs\n");
   fprintf(stderr,"ncol <= 2(q^n-1)(q-1) -1. Can't have ncol = %d\n",ncol);
   fprintf(stderr,"with n=%d and q = %d,\n",akn,q);
  return 0;
}

return 1;
}


addelkempn( gf, akn, A, ncol )
/* Implement Addelman and Kempthorne's 1961 A.M.S. method with n=3 */
struct GF *gf;
int **A, ncol, akn;
{
int i,p,q;
int kay,*b,*c,*k;  /* A&K notation */
int row, col, square, ksquare;
int monic, numin, *x, *coef, *indx, *s, poly, elt, sub;

p=gf->p, q=gf->q;

if(  !addelkempncheck( q,p,akn,ncol  )  )return 0;

b = ivector( 0,q-1 );
c = ivector( 0,q-1 );
k = ivector( 0,q-1 );

x = ivector( 0,akn-1 );
s    = ivector( 0,akn-1 );
coef = ivector( 0,akn-1 );
indx = ivector( 0,akn-1 );

for(  i=0; i<akn; i++  )
  x[i] = 0;
for(  row=0; row<ipow(q,akn); row++  ){  /* First q^akn rows */
  /*printf("x: ");for(i=0;i<akn;i++)printf("%d ",x[i]);printf("\n");*/
  col = 0;
  s[0] = 1;          
  for(  i=1; i<akn; i++  )               /* first subset */
    s[i] = 0;           /* nonempty subsets of x indices */
  for(  sub=1; sub<ipow(2,akn) && col<ncol; sub++  ){
/*    if(!row)
      {printf("s: ");for(i=0;i<akn;i++)printf("%d ",s[i]);printf("\n");}*/
    monic = -1;
    numin = 0;
    for(  i=0; i<akn; i++  )
      if( s[i]  )
	if(  monic == -1  )
	  monic = i;
	else
	  indx[numin++] = i;
    for(  i=0; i<numin; i++  )
      coef[i] = 1;
    for(  poly=0; poly<ipow(q-1,numin) && col<ncol; poly++  ){
/*      if(  row==0  ){
	printf("  p: ");for(i=0;i<numin;i++)printf("%d ",coef[i]);printf("\n");
      }*/
      elt = x[monic];
      for(  i=0; i<numin; i++  )
	elt = gf->plus[elt][gf->times[coef[i]][x[indx[i]]]];
      A[row][col++] = elt;
      for( i=numin-1; i>=0; i--  ){
	coef[i] = (coef[i]+1) % q;
	if(  coef[i]  )
	  break;
	else
	  coef[i] = 1;
      }
    }
    for(  i=0; i<akn; i++  ){
      s[i] = (s[i]+1) % 2;
      if(  s[i]  )break;
    }      
  }

  square = gf->times[x[0]][x[0]];

  s[1] = 1;          
  for(  i=2; i<akn; i++  )               /* first subset */
    s[i] = 0;           /* nonempty subsets of x indices */
  for(  sub=1; sub<ipow(2,akn-1) && col<ncol; sub++  ){
/*    if(!row)
      {printf("s: ");for(i=0;i<akn;i++)printf("%d ",s[i]);printf("\n");}*/
    monic = -1;
    numin = 0;
    for(  i=1; i<akn; i++  )
      if( s[i]  )
	if(  monic == -1  )
	  monic = i;
	else
	  indx[numin++] = i;
    coef[0] = 0;
    for(  i=1; i<numin+1; i++  )
      coef[i] = 1;
    for(  poly=0; poly<q*ipow(q-1,numin) && col<ncol; poly++  ){
/*      if(  !row  ){
	printf("  p: ");for(i=0;i<numin+1;i++)printf("%d ",coef[i]);printf("\n");
      }*/
      elt = gf->plus[square][gf->times[x[0]][coef[0]]];
      elt = gf->plus[elt][x[monic]];
      for(  i=1; i<numin+1; i++  )
	elt = gf->plus[elt][gf->times[coef[i]][x[indx[i-1]]]];
      A[row][col++] = elt;
      for( i=numin+1-1; i>=0; i--  ){
	coef[i] = (coef[i]+1) % q;
	if(  coef[i]  )
	  break;
	else
	  if( i>0 )coef[i] = 1;
      }
    }
    for(  i=1; i<akn; i++  ){
      s[i] = (s[i]+1) % 2;
      if(  s[i]  )break;
    }
  }

  for(  i=akn-1; i>=0; i--  ){
    x[i] = (x[i]+1) % q;
    if(  x[i]  )break;
  }
}

if(  p !=2  )                    /* Constants kay,b,c,k for odd p */
  akodd(  gf,&kay,b,c,k );
else                             /* Constants kay,b,c,k for even p */
  akeven( gf,&kay,b,c,k );

for(  i=0; i<akn; i++  )
  x[i] = 0;
for(  row=ipow(q,akn); row<2*ipow(q,akn); row++  ){  /* Second q^akn rows */
  col = 0;
  s[0] = 1;          
  for(  i=1; i<akn; i++  )               /* first subset */
    s[i] = 0;           /* nonempty subsets of x indices */
  for(  sub=1; sub<ipow(2,akn) && col<ncol; sub++  ){
    monic = -1;
    numin = 0;
    for(  i=0; i<akn; i++  )
      if( s[i]  )
	if(  monic == -1  )
	  monic = i;
	else
	  indx[numin++] = i;
    for(  i=0; i<numin; i++  )
      coef[i] = 1;
    for(  poly=0; poly<ipow(q-1,numin) && col<ncol; poly++  ){
      elt = x[monic];
      if(  numin && s[0] )
	elt = gf->plus[elt][b[coef[0]]];
      for(  i=0; i<numin; i++  )
	elt = gf->plus[elt][gf->times[coef[i]][x[indx[i]]]];
      A[row][col++] = elt;
      for( i=numin-1; i>=0; i--  ){
	coef[i] = (coef[i]+1) % q;
	if(  coef[i]  )
	  break;
	else
	  coef[i] = 1;
      }
    }
    for(  i=0; i<akn; i++  ){
      s[i] = (s[i]+1) % 2;
      if(  s[i]  )break;
    }      
  }

  ksquare = gf->times[kay][gf->times[x[0]][x[0]]];

  s[1] = 1;          
  for(  i=2; i<akn; i++  )               /* first subset */
    s[i] = 0;           /* nonempty subsets of x indices */
  for(  sub=1; sub<ipow(2,akn-1) && col<ncol; sub++  ){
    monic = -1;
    numin = 0;
    for(  i=1; i<akn; i++  )
      if( s[i]  )
	if(  monic == -1  )
	  monic = i;
	else
	  indx[numin++] = i;
    coef[0] = 0;
    for(  i=1; i<numin+1; i++  )
      coef[i] = 1;
    for(  poly=0; poly<q*ipow(q-1,numin) && col<ncol; poly++  ){
      elt = gf->plus[ksquare][gf->times[x[0]][k[coef[0]]]];
      elt = gf->plus[elt][x[monic]];
      elt = gf->plus[elt][c[coef[0]]];
      for(  i=1; i<numin+1; i++  )
	elt = gf->plus[elt][gf->times[coef[i]][x[indx[i-1]]]];
      A[row][col++] = elt;
      for( i=numin+1-1; i>=0; i--  ){
	coef[i] = (coef[i]+1) % q;
	if(  coef[i]  )
	  break;
	else
	  coef[i] = i>0 ? 1 : 0;
      }
    }
    for(  i=1; i<akn; i++  ){
      s[i] = (s[i]+1) % 2;
      if(  s[i]  )break;
    }
  }

  for(  i=akn-1; i>=0; i--  ){
    x[i] = (x[i]+1) % q;
    if(  x[i]  )break;
  }
}

return 1;
}





