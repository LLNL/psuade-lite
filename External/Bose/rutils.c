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
#include <math.h>
#include <stdlib.h>

double **dmatrix(), *dvector();


/*

 Utilities based on runif

*/

unifperm( pi,q )  
/* 
   In S one just does rank(runif(q)).  Here we want
something like rank(runif(q))-1 since the symbols to
be permuted are 0..q-1.  
*/

int q, *pi;
{
int i;
double *z;

z = dvector( 0,q-1 );
if(  !z  ){
  fprintf(stderr,"Could not allocate memory for random permutation.\n");
  exit(1);
}
runif( z,q );

findranks( q,z,pi );

for(  i=0; i<q; i++  )
  pi[i] -= 1;

free_dvector( z,0,q-1);
}


/* These routines perform sorting and ranking chores. */

int rankcomp( a,b )

double **a, **b;
{
/*printf("Randcomp a = (%g,%g) b=(%g,%g)\n",a[0][0],a[0][1],b[0][0],b[0][1]);*/
if(  a[0][0] < b[0][0]  )return(-1);
if(  a[0][0] > b[0][0]  )return( 1);
return(0);
}

findranks( n,v,r )
int n, *r;
double *v;

{
double **temp;
int    i;

temp = dmatrix(0,n-1,0,1);

if(  !temp  ){
  printf("findranks: could not allocate memory to find ranks.\n");
  exit(1);
}

for( i=0; i<n; i++  ){
  temp[i][0] = v[i];
  temp[i][1] = (double) i;
}
qsort((void *)temp,(size_t) n,sizeof(temp[0]),rankcomp);

for(  i=0; i<n; i++  )
  r[ (int)temp[i][1] ] = i+1; /*Ranks go 1..n */
}


int doubcomp( a,b )
double a,b;
{
if(  a < b  )return(-1);
if(  a > b  )return( 1);
return(0);
}
