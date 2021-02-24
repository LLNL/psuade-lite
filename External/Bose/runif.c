#include <stdio.h>
#include <math.h>




/*
c
c Marsaglia - Zaman universal random number generator.
c
c reinitialization: call seed(is,js,ks,ls), with integer arguments
c from 1 to 168, not all 1.
c
c generate n uniform random numbers and store in x(n): call ranums(x,n).
c


  Transliterated from FORTRAN to C by Art Owen, 4 March 1993.

  Functions:

    mod( a,b )             a mod b
    seedok( is,js,ks,ls )  1 if seeds ok, 0 otherwise
    seed( is,js,ks,ls )    sets seed integers, rejects invalid input
    ranums( z,n )          sets z[0] through z[n-1] 
                           to the next n random uniforms between 0 and 1
    runif( z,n )           same as ranums

  The C and f77 programs shown below find the first 2000 random
  numbers starting at 0 and print them to 10 places.  The maximum
  difference between any two corresponding random numbers was
  5e-11.  Equality held in 1822 of the cases.  The original FORTRAN 
  subroutines are appended in a comment at the end of this file.


C test program:

#include <math.h>
#include <stdio.h>

#define N 2000
static  double z[N];

main()
{
int i;
ranums( z,N );
for(  i=0; i<N; i++  )
  printf("%14.10g\n",z[i]);
}

F77 test program:

      integer n
      real z(2000)

      n = 2000
      call ranums( z,n )
      write(6,10)(z(i),i=1,n)
 10   format(F14.10)

      stop
      end

*/


int mod( a,b )
int  a,b;
{
int ans;
ans = a % b;
if(  ans>=0  )return ans;
else return ans+b;
}


static int      jent=0,i=12,j=34,k=56,l=78,      ip,jp;
static double u[97+1],                          c,cd,cm;

seedok( is,js,ks,ls )   /*  1 iff seed is ok   */
int is,js,ks,ls;
{
if(  is==1  && js==1  && ks==1  && ls==1   )return 0;
if(  is<1   || js<1   || ks<1   || ls<1    )return 0;
if(  is>168 || js>168 || ks>168 || ls>168  )return 0;
return 1;
}


seed( is,js,ks,ls )
int is,js,ks,ls;
{
jent=0;

if(  seedok(is,js,ks,ls)  ){
  i=is;  j=js;  k=ks;  l=ls;
}else{
  fprintf(stderr,"Error: Invalid seed %d %d %d %d\n",is,js,ks,ls);
  fprintf(stderr,"Must be four integers between 1 and 168, and\n");
  fprintf(stderr,"must not all be 1.  Seed not changed.\n");
}
}

runif(x,n)
int     n;
double *x;
{
ranums(x,n);
}

ranums(x,n)
int     n;
double *x;
{
int    ii, jj, m;
double s,t,uni;

if(jent != 0) goto L30;
jent=1;
for(  ii=1; ii<=97; ii++  ){     /* do 20 ii=1,97 */
  s=0.0;
  t=0.5;
  for(  jj=1; jj<=24; jj++  ){   /* do 10 jj=1,24 */
    m=mod(mod(i*j,179)*k,179);
    i=j;
    j=k;
    k=m;
    l=mod(53*l+1,169);
    if (mod(l*m,64) >= 32) s=s+t;
    t=0.5*t;
  }                              /* 10   continue */
  u[ii]=s;
}                                /* 20   continue */

c  =   362436.0/16777216.0;
cd =  7654321.0/16777216.0;
cm = 16777213.0/16777216.0;
ip = 97;
jp = 33;

 L30:  for(  ii=1; ii<=n; ii++  ){ /*  ii do 40 ii=1,n */
   uni=u[ip]-u[jp];
   if (uni < 0.0) uni=uni+1.0;
   u[ip]=uni;
   ip=ip-1;
   if (ip == 0) ip=97;
   jp=jp-1;
   if (jp == 0) jp=97;
   c=c-cd;
   if (c < 0.0) c=c+cm;
   uni=uni-c;
   if (uni < 0.0) uni=uni+1.0;
   x[ii-1]=uni;
 }                                 /* 40   continue */
}

/*  Original FORTRAN file:
c
      subroutine seed (is,js,ks,ls)
c
c Marsaglia - Zaman universal random number generator.
c
c reinitialization: call seed(is,js,ks,ls), with integer arguments
c from 1 to 168, not all 1.
c
c generate n uniform random numbers and store in x(n): call ranums(x,n).
c
      real u(97),x(1)
      save u,jent,i,j,k,l,c,cd,cm,ip,jp
      data jent,i,j,k,l /0,12,34,56,78/
      jent=0
      i=is
      j=js
      k=ks
      l=ls
      return
      entry ranums(x,n)
      if (jent.ne.0) go to 30
      jent=1
      do 20 ii=1,97
      s=0.0
      t=0.5
      do 10 jj=1,24
      m=mod(mod(i*j,179)*k,179)
      i=j
      j=k
      k=m
      l=mod(53*l+1,169)
      if (mod(l*m,64).ge.32) s=s+t
      t=0.5*t
 10   continue
      u(ii)=s
 20   continue
      c=362436.0/16777216.0
      cd=7654321.0/16777216.0
      cm=16777213.0/16777216.0
      ip=97
      jp=33
 30   do 40 ii=1,n
      uni=u(ip)-u(jp)
      if (uni.lt.0.0) uni=uni+1.0
      u(ip)=uni
      ip=ip-1
      if (ip.eq.0) ip=97
      jp=jp-1
      if (jp.eq.0) jp=97
      c=c-cd
      if (c.lt.0.0) c=c+cm
      uni=uni-c
      if (uni.lt.0.0) uni=uni+1.0
      x(ii)=uni
 40   continue
      return
      end
*/
