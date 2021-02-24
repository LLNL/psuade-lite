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
#include "oa.h"

/*  

       Check whether the array in standard input is really
   of strength 2.  Use brute force.  For OA( nrow, ncol, q, ? )
   it takes work roughly proportional to 
             ncol^2 * nrow * q^2/2
   to decide if ? >= 2.  The user is warned if this is likely
   to be too much work.

       The program calls exit(0) if the input array has strength
   2.  It calls exit(1) if the array is not of strength 2, or if
   the input is invalid, or if it is impossible to allocate enough
   memory to find out.

       The program exits at the first sign that the array is
   not of strength 2.  This can save lots of work if the problem
   shows up early, but it doesn't give a complete list of the
   array's shortcomings.  Such a list could be very large.

       Note that an array of strength larger than 2 is a fortiori
   of strength 2 and will pass this test.  

*/

main(argc,argv)
int  argc;
char *argv[];
{
int q, nrow, ncol, **A;

OA_parsein( argc,argv, &q, &nrow, &ncol, &A );
if(  OA_str2( q,nrow,ncol,A,2   )  )
  exit(1);
else
  exit(0);
}
