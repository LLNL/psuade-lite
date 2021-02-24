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


/*     Construction of specific Galois fields.  See galois.c 
    for manipulation of generic Galois fields.

       To add a new Galois field GF(q)=GF(p^n) it is necessary 
    to find the polynomial representation of x^n, and use it to
    include another *xtn--- vector like those used below.

*/

#include <math.h>
#include <stdio.h>

#include "galois.h"


/*  Glossary:

       xtn2t2              Polynomial representation of x^2 in gf2t2
       xtn2t3              Polynomial representation of x^3 in gf2t3
       . . .               Polynomial representation of p^n in gfptn
       xtn13t2             Polynomial representation of x^2 in gf13t2

       GF_set_fields       Initialize polynomial representations
       GF_fields_are_set   Indicates initialization done
       GF_getfield         Construct and return GF(q) if possible
*/

int *xtn2t2, *xtn2t3, *xtn2t4, *xtn2t5, *xtn2t6, *xtn2t7, *xtn2t8, *xtn2t9;
int *xtn3t2, *xtn3t3, *xtn3t4, *xtn3t5, *xtn3t6;
int *xtn5t2, *xtn5t3, *xtn5t4, *xtn5t5, *xtn5t6;
int *xtn7t2, *xtn7t3, *xtn7t4, *xtn7t5;
int *xtn11t2;
int *xtn13t2;

int *xtnpt1;

int GF_fields_are_set = 0;

GF_set_fields()
{
/* Brute force set up of defining vectors, from Carmichael */

/* Declare x-to-the-power-n vectors, for GFs p-to-the-n */

if(   GF_fields_are_set   )
  fprintf(stderr,"Warning: Fields being re-initialized.  Possible memory waste.\n");

xtnpt1 = ivector(0,0);

xtn2t2 = ivector(0,1);
xtn3t2 = ivector(0,1);
xtn5t2 = ivector(0,1);
xtn7t2 = ivector(0,1);
xtn11t2 = ivector(0,1);
xtn13t2 = ivector(0,1);

xtn2t3 = ivector(0,2);
xtn3t3 = ivector(0,2);
xtn5t3 = ivector(0,2);
xtn7t3 = ivector(0,2);

xtn2t4 = ivector(0,3);
xtn3t4 = ivector(0,3);
xtn5t4 = ivector(0,3);
xtn7t4 = ivector(0,3);

xtn2t5 = ivector(0,4);
xtn3t5 = ivector(0,4);
xtn5t5 = ivector(0,4);
xtn7t5 = ivector(0,4);

xtn2t6 = ivector(0,5);
xtn3t6 = ivector(0,5);
xtn5t6 = ivector(0,5);

xtn2t7 = ivector(0,6);
xtn2t8 = ivector(0,7);
xtn2t9 = ivector(0,8);

/* Assign values for vectors with p=2 */

xtn2t2[0] = 1;
xtn2t2[1] = 1;

xtn2t3[0] = 1;
xtn2t3[1] = 1;
xtn2t3[2] = 0;

xtn2t4[0] = 1;
xtn2t4[1] = 1;
xtn2t4[2] = 0;
xtn2t4[3] = 0;

xtn2t5[0] = 1;
xtn2t5[1] = 0;
xtn2t5[2] = 1;
xtn2t5[3] = 0;
xtn2t5[4] = 0;

xtn2t6[0] = 1;
xtn2t6[1] = 1;
xtn2t6[2] = 0;
xtn2t6[3] = 0;
xtn2t6[4] = 0;
xtn2t6[5] = 0;

xtn2t7[0] = 1;
xtn2t7[1] = 1;
xtn2t7[2] = 0;
xtn2t7[3] = 0;
xtn2t7[4] = 0;
xtn2t7[5] = 0;
xtn2t7[6] = 0;

xtn2t8[0] = 1;
xtn2t8[1] = 0;
xtn2t8[2] = 1;
xtn2t8[3] = 1;
xtn2t8[4] = 1;
xtn2t8[5] = 0;
xtn2t8[6] = 0;
xtn2t8[7] = 0;

xtn2t9[0] = 1;
xtn2t9[1] = 1;
xtn2t9[2] = 1;
xtn2t9[3] = 1;
xtn2t9[4] = 1;
xtn2t9[5] = 0;
xtn2t9[6] = 0;
xtn2t9[7] = 0;
xtn2t9[8] = 1;

/* Assign values for vectors with p=3 */

xtn3t2[0] = 1;
xtn3t2[1] = 2;

xtn3t3[0] = 2;
xtn3t3[1] = 1;
xtn3t3[2] = 0;

xtn3t4[0] = 1;
xtn3t4[1] = 1;
xtn3t4[2] = 2;
xtn3t4[3] = 2;

xtn3t5[0] = 2;
xtn3t5[1] = 1;
xtn3t5[2] = 0;
xtn3t5[3] = 0;
xtn3t5[4] = 0;

xtn3t6[0] = 1;
xtn3t6[1] = 1;
xtn3t6[2] = 0;
xtn3t6[3] = 0;
xtn3t6[4] = 0;
xtn3t6[5] = 0;

/* Assign values for vectors with p=5 */

xtn5t2[0] = 2;
xtn5t2[1] = 2;

xtn5t3[0] = 3;
xtn5t3[1] = 2;
xtn5t3[2] = 0;

xtn5t4[0] = 2;
xtn5t4[1] = 1;
xtn5t4[2] = 0;
xtn5t4[3] = 1;

xtn5t5[0] = 2;
xtn5t5[1] = 1;
xtn5t5[2] = 0;
xtn5t5[3] = 0;
xtn5t5[4] = 0;

xtn5t6[0] = 3; /* uses some negs in Carmichael */
xtn5t6[1] = 3;
xtn5t6[2] = 0;
xtn5t6[3] = 1;
xtn5t6[4] = 4;
xtn5t6[5] = 1;

/* Assign values for vectors with p=7 */

xtn7t2[0] = 4; /* uses some negs in Carmichael */
xtn7t2[1] = 1;

xtn7t3[0] = 5; /* uses some negs in Carmichael */
xtn7t3[1] = 1;
xtn7t3[2] = 0;

xtn7t4[0] = 2;
xtn7t4[1] = 2;
xtn7t4[2] = 0;
xtn7t4[3] = 2;

xtn7t5[0] = 3;
xtn7t5[1] = 6;
xtn7t5[2] = 0;
xtn7t5[3] = 0;
xtn7t5[4] = 0;

/* Assign values for vectors with p=11,13 */

xtn11t2[0] = 9; /* uses some negs in Carmichael */
xtn11t2[1] = 4;

xtn13t2[0] = 11;
xtn13t2[1] = 12;

xtnpt1[0] = 0; /* Not used */

GF_fields_are_set = 1;
}


GF_getfield( q, gf )
int q;
struct GF *gf;
{
int *xtn;
int p,n,ispp;

if(  !GF_fields_are_set  )
  GF_set_fields();

if(  q<1  ){      /* Impossible argument */
  fprintf(stderr,"Field must have positive number of elements.\n");
  return 0; }
if(  q==1 ){      /* Pointless  argument */
  fprintf(stderr,"Field with 1 element was requested.\n");
  return 0; }

primepow( q,&p,&n,&ispp  );
if(  !ispp  ){
  fprintf(stderr,"q=%d is not a prime power.\n",q);
  return 0; }

xtn = NULL;

/*   4 8 16 32 64 128 256 512 */
if(  q== ipow(2,2)  )xtn = xtn2t2;
if(  q== ipow(2,3)  )xtn = xtn2t3;
if(  q== ipow(2,4)  )xtn = xtn2t4;
if(  q== ipow(2,5)  )xtn = xtn2t5;
if(  q== ipow(2,6)  )xtn = xtn2t6;
if(  q== ipow(2,7)  )xtn = xtn2t7;
if(  q== ipow(2,8)  )xtn = xtn2t8;
if(  q== ipow(2,9)  )xtn = xtn2t9;

/*   9 27 81 243 729          */
if(  q== ipow(3,2)  )xtn = xtn3t2;
if(  q== ipow(3,3)  )xtn = xtn3t3;
if(  q== ipow(3,4)  )xtn = xtn3t4;
if(  q== ipow(3,5)  )xtn = xtn3t5;
if(  q== ipow(3,6)  )xtn = xtn3t6;

/*   25 125 625 3125 15625    */
if(  q== ipow(5,2)  )xtn = xtn5t2;
if(  q== ipow(5,3)  )xtn = xtn5t3;
if(  q== ipow(5,4)  )xtn = xtn5t4;
if(  q== ipow(5,5)  )xtn = xtn5t5;
if(  q== ipow(5,6)  )xtn = xtn5t6;

/*   49 343 2401 16807        */
if(  q== ipow(7,2)  )xtn = xtn7t2;
if(  q== ipow(7,3)  )xtn = xtn7t3;
if(  q== ipow(7,4)  )xtn = xtn7t4;
if(  q== ipow(7,5)  )xtn = xtn7t5;

/*   121 169                  */
if(  q== ipow(11,2) )xtn = xtn11t2;
if(  q== ipow(13,2) )xtn = xtn13t2;

if(  isprime(q)  )xtn = xtnpt1;  /* Could have tested p=q, or n=1 */

if(  xtn   ){
  if(  GF_ready( gf,p,n,xtn )  )
    return 1;
  else{
    fprintf(stderr,"Construction failed for GF(%d).\n",q);
    return 0;
  }
}
else {
  fprintf(stderr,"GF(%d) = GF(%d^%d) is not included in this program.\n",q,p,n);
  fprintf(stderr,"To add it, consider modifying gfields.c.\n",q);
  return 0;
}
}
