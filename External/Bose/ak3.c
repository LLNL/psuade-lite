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


/*     Addelman-Kempthorne Constructions with n = 3.

       The article is quite vague on this.  Page 1173 states
   "When n>2 the same procedure will yield the desired plans
   if Lemma 5a is used in place of Lemma 5."  Page 1175
   provides the example n=3,q=3 which is OA( 54,25,3,2 ).
   Based on this example it is possible to make an educated
   guess as to how the construction generalizes to n=3.
   The resulting OA's are seen, by brute force to be of
   strength 2 for q=2,3,4,5,7,11.  These OAs are:
         OA(   16,  13,  2, 2 )   
	 OA(   54,  25,  3, 2 )
	 OA(  128,  41,  4, 2 )   
	 OA(  250,  61,  5, 2 )
	 OA(  686, 113,  7, 2 )
	 OA( 1458, 181,  9, 2 )
	 OA( 2662, 265, 11, 2 )
   The one with q=7 required 212709392 comparisons to determine
   that it really is of strength 2.  This took roughly 11.5 minutes
   on a DEC 5000/240 workstation (real and elapsed in this case).
   
   The array with q=11 took 1.12671e+10 comparisons to verify its strength.
   This took roughly 10 1/2 hours.
   Here is the tail end of the output from oacheck2:
   No violation of strength 2 involves column 262.
   No violation of strength 2 involves column 263.
   No violation of strength 2 involves column 264.
   The array is of strength (at least) 2.
   [2]    Done                 aktest | oacheck2 11 2662 265
   37890.0u 31.7s 10:33:11 99% 61+110k 0+0io 0pf+0w


*/

#include <math.h>
#include <stdio.h>

#include "galois.h"

addelkemp3check( q,p,ncol  )
int q,p,ncol;
{

if(  p==2 && q>4 ){
  fprintf(stderr,"This Addelman-Kempthorne OA(2q^3,ncol,q,2) is only\n");
  fprintf(stderr,"available for odd prime powers q and for even prime\n");
  fprintf(stderr,"powers q<=4.\n");
  return 0;
}

if(  q==8  ){  /* Moot */
  fprintf(stderr,"This Addelman-Kempthorne OA(2*8^3,ncol,8,2) is\n");
  fprintf(stderr,"experimental and not yet working.\n");
  return 1;
}

if(  ncol > 2*q*q + 2*q + 1  ){
  fprintf(stderr,"The Addelman-Kempthorne (n=3) construction needs\n");
   fprintf(stderr,"ncol <= 2q^2+2q+1. Can't have ncol = %d with q = %d,\n",ncol,q);
  return 0;
}

return 1;
}



addelkemp3( gf, A, ncol )
/* Implement Addelman and Kempthorne's 1961 A.M.S. method with n=3 */
struct GF *gf;
int **A, ncol;
{
int i1,i2,i3, m1,m2, p,q;
int kay,*b,*c,*k;  /* A&K notation */
int row, col, square, ksquare;

p=gf->p, q=gf->q;

if(  !addelkemp3check( q,p,ncol  )  )return 0;

b = ivector( 0,q-1 );
c = ivector( 0,q-1 );
k = ivector( 0,q-1 );

for(  i1=0; i1<q; i1++  ){           /* First q^3 rows */
  square = gf->times[i1][i1];
  for(  i2=0; i2<q; i2++  )
  for(  i3=0; i3<q; i3++  ){
    row = i3 + q*i2 + q*q*i1;
    col = 0;
    if(  col<ncol  )A[row][col++]=i2;                     /*      y       */
    for(  m1=1; m1<q && col<ncol; m1++  )                 /* x + my       */
      A[row][col++] = gf->plus[i1][gf->times[m1][i2]];
    if(  col<ncol  )A[row][col++]=i3;                     /*           z  */
    for(  m2=1; m2<q && col<ncol; m2++  )                 /* x      + mz  */
      A[row][col++] = gf->plus[i1][gf->times[m2][i3]];
    for(  m2=1; m2<q && col<ncol; m2++  )                 /*      y + mz  */
      A[row][col++] = gf->plus[i2][gf->times[m2][i3]];
    for(  m1=1; m1<q && col<ncol; m1++  )                 /* x + my + nz  */
    for(  m2=1; m2<q && col<ncol; m2++  )
      A[row][col++] = 
        gf->plus[i1][gf->plus[gf->times[m1][i2]][gf->times[m2][i3]]];
    for(  m1=0; m1<q && col<ncol; m1++  )                 /* x^2 + mx + y */
      A[row][col++] = gf->plus[square][
                       gf->plus[i2][
                       gf->times[m1][i1]]];

    for(  m1=0; m1<q && col<ncol; m1++  )                 /* x^2 + mx + z */
      A[row][col++] = gf->plus[square][
		       gf->plus[i3][
                       gf->times[m1][i1]]];

    for(  m1=0; m1<q && col<ncol; m1++  )              /* x^2 + mx + y + nz */
    for(  m2=1; m2<q && col<ncol; m2++  )              
      A[row][col++] 
	= gf->plus[square][
            gf->plus[i2][
              gf->plus[ gf->times[m2][i3] ][
                        gf->times[m1][i1] 
                ]
              ]
            ];
    if(  col<ncol  )A[row][col++]=i1;                     /* x            */
  }
}

if(  p!=2  )
  akodd(  gf,&kay,b,c,k ); /* Get kay,b,c,k for odd p  */
else                             
  akeven( gf,&kay,b,c,k ); /* Constants kay,b,c,k for even p */

for(  i1=0; i1<q; i1++  ){           /* Second q^3 rows */
  square = gf->times[i1][i1];
  ksquare = gf->times[kay][square];
  for(  i2=0; i2<q; i2++  )
  for(  i3=0; i3<q; i3++  ){
    row = i3 + q*i2 + q*q*i1 + q*q*q;
    col = 0;
    if(  col<ncol  )A[row][col++]=i2;                /*     y        */

    for(  m1=1; m1<q && col<ncol; m1++  ){           /* x + my + b(m)      */
      A[row][col] = gf->plus[i1][gf->times[m1][i2]];
      A[row][col] = gf->plus[A[row][col]][b[m1]];
      col++;
    }
    if(  col<ncol  )A[row][col++]=i3;                /*           z  */

    for(  m2=1; m2<q && col<ncol; m2++  ){      /* x      + mz + b(m) */
      A[row][col] =   gf->plus[i1][gf->times[m2][i3]];
      A[row][col] = gf->plus[A[row][col]][b[m2]];
      col++;
    }
    for(  m2=1; m2<q && col<ncol; m2++  )            /*      y + mz  */
      A[row][col++] = gf->plus[i2][gf->times[m2][i3]];

    for(  m1=1; m1<q && col<ncol; m1++  )      /* x + my + nz + b(m) */
    for(  m2=1; m2<q && col<ncol; m2++  ){
      A[row][col] = 
        gf->plus[i1][gf->plus[gf->times[m1][i2]][gf->times[m2][i3]]];
      A[row][col] = gf->plus[A[row][col]][b[m1]];
      col++;
    }      

    for(  m1=0; m1<q && col<ncol; m1++  ){      /* kx^2 + k(m)x + y + c(m)*/
      A[row][col] = gf->plus[ksquare][
                       gf->plus[i2][
                       gf->times[k[m1]][i1]]];
      A[row][col] = gf->plus[A[row][col]][c[m1]];
      col++;
    }

    for(  m1=0; m1<q && col<ncol; m1++  ){      /* kx^2 + k(m)x + z + c(m)*/
      A[row][col] = gf->plus[ksquare][
                       gf->plus[i3][
                       gf->times[k[m1]][i1]]];
      A[row][col] = gf->plus[A[row][col]][c[m1]];
      col++;
    }

    for(  m1=0; m1<q && col<ncol; m1++  )    /* kx^2 + k(m)x + y + nz +c(m) */
    for(  m2=1; m2<q && col<ncol; m2++  ){
      A[row][col] 
	= gf->plus[ksquare][
            gf->plus[i2][
              gf->plus[ gf->times[m2][i3] ][
                        gf->times[k[m1]][i1] 
                ]
              ]
            ];
      A[row][col] = gf->plus[ A[row][col] ][ c[m1] ];
      col++;
    }
    if(  col<ncol  )A[row][col++]=i1;                /* x            */
  }
}
return 1;
}




