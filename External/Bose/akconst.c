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

/*  Find constants for Addelman Kempthorne designs
  when q is even. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "galois.h"

int  *ivector(), **imatrix();


/*  EVEN    EVEN    EVEN    EVEN    EVEN    EVEN    EVEN    EVEN  */

/*  For even q, only q= 2 or 4 are available.  The prescription
given in Addelman and Kempthorne (1961) does not appear to work.
Commented out code below attempts to implement that prescription.
It seemed to be impossible to find a constant b[1],c[1] pair.
*/


int akeven( gf, kay, b, c, k )
struct GF *gf;
int        *kay, *b, *c, *k;
{
/*
int   **posb, **posc;
int   j, temp, square;
int   v1, v2, x1, x2, z1, z2, ib, ic;
int   *vals1, *vals2, done;
*/
int   i,q,xtn;

q = gf->q;
*kay=1;

if(  q==2  )
  b[1] = c[1] = k[1] = 1;
if(  q==4  ){
  b[1] = c[1] = 2;
  b[2] = c[2] = 1;
  b[3] = c[3] = 3;
  k[1] = 1; k[2] = 2; k[3] = 3;
}

for(  i=1; i<q; i++  )
  k[i]=i;

if(  q>4  ){
  fprintf(stderr,"Addelman Kempthorne designs not yet available for\n");
  fprintf(stderr,"even q >4.");
  exit(1);

/*
  fprintf(stderr,"Using experimental values in akconst.\n");
  xtn = 1;
  for(  i=0; i<gf->n; i++  )
    xtn = gf->times[xtn][2];
  xtn = gf->times[xtn][xtn];
  for(  i=1; i<q; i++  ){
    b[i] = gf->times[xtn][gf->inv[i]];
    c[i] = gf->times[xtn][gf->times[i][i]];
  }
*/
}
return 1;
}

/*
posb  = imatrix(0,q-1,0,q-1);
posc  = imatrix(0,q-1,0,q-1);
vals1 = ivector(0,q-1);
vals2 = ivector(0,q-1);

for(  i=1; i<q; i++  )
for(  j=1; j<q; j++  )
  posb[i][j]=posc[i][j] = 1;
*/


/* Rule out some b */
/*
for(  i =1; i <q;  i++  )       
for(  x1=0; x1<q; x1++  ){
  temp  = gf->times[x1][gf->inv[i]];
  square= gf->times[x1][x1];
  temp  = gf->plus[square][gf->neg[temp]];
  posb[i][temp] = 0;*/        /* b[i] != temp */
/*}*/

/* Rule out some c */
/*
for(  i =1; i <q;  i++  )       
for(  x1=0; x1<q; x1++  ){
  temp = gf->times[x1][i];
  square=gf->times[x1][x1];
  temp = gf->plus[square][temp];
  posc[i][temp] = 0;*/          /* c[i] != temp */
/*}*/

/*
printf("Possible b\n\n");
OA_put( posb,q,q,2 );
printf("Possible c\n\n");
OA_put( posc,q,q,2 );
*/

/* Search for b[i],c[i]  */
/*
for(  i=1; i<q; i++  ){       
  done = 0;
  for(  ib=1; ib<q && !done; ib++  ){
    if(  !posb[i][ib]  )continue;
    for(  ic=1; ic<q && !done; ic++  ){
      if(  !posc[i][ic]  )continue;
      for(  j=0; j<q; j++  )
	vals1[j] = vals2[j] = 0;
      for(  x1=0; x1<q; x1++  ){
	square = gf->times[x1][x1];
	x2 = gf->plus[x1][ib];
	x2 = gf->times[x2][gf->inv[i]];
	x2 = gf->neg[x2];
	v1 = gf->times[k[i]][x1];
	v1 = gf->plus[v1][square];
	v1 = gf->plus[v1][x2];
	v1 = gf->plus[v1][ic];
	vals1[ v1 ] = 1;

	x2 = gf->times[x1][gf->inv[i]];
	x2 = gf->neg[x2];
	v2 = gf->times[k[i]][x1];
	v2 = gf->plus[v2][square];
	v2 = gf->plus[v2][x2];
	vals2[ v2 ] =1;
      }
      for(  j=0; j<q; j++  )
	done += (vals2[j] + vals1[j] == 1);
      if(  done<q  )
	done=0;
      else{
	b[i] = ib;
	c[i] = ic;
      }
    }
  }
  if(  !done  )
    printf("No b[%d] c[%d] combination works.\n",i,i);
  else
    printf("b[%d] = %d and c[%d] = %d works\n",i,b[i],i,c[i]);
}
}
*/





/*  ODD    ODD    ODD    ODD    ODD    ODD    ODD    ODD    ODD  */

int akodd( gf, kay, b, c, k )
struct GF *gf;
int        *kay, *b, *c, *k;
{
int   i, q, p, num, den, four;

q = gf->q; p = gf->p;

if(  p!=3  )
  four = 4;
else
  four = 1;

*kay=0;
for(  i=2; i<q; i++  )
  if( gf->root[i] == -1 )*kay=i;
if(  *kay==0  ){
  fprintf(stderr,"Problem: no rootless element in GF(%d).\n",gf->n);
  return 0;
}

for(  i=1; i<q; i++  ){
  num = gf->plus[*kay][p-1];  /* -1 = +(p-1) */
  den = gf->times[*kay][four];
  den = gf->times[den][i];
  b[i]= gf->times[num][gf->inv[den]];
  k[i]= gf->times[*kay][i];
  c[i]= gf->times[i][i];
  c[i]= gf->times[c[i]][num];
  c[i]= gf->times[c[i]][gf->inv[four]];
  /*    printf("i,num,den,b,k,c %3d %3d %3d %3d %3d %3d\n",i,num,den,b[i],k[i],c[i]);*/
}
return 1;
}
