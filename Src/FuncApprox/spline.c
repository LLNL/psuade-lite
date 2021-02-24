#include <stdio.h>
#include <stdlib.h>

/********************************************************************
 * one dimensional spline preprocessing
 * ****************************************************************/
int spline1p(int leng, double *X, double *Y, double *YY)
{
   int    nn;
   double *U, sig, p, ddata;
   /* ======================================= */
   /* compute second derivatives              */
   /* ======================================= */
   U   = (double *) malloc(leng * sizeof(double));
   YY[0] = YY[leng-1] = U[0] = 0;
   for (nn = 1; nn < leng-1; nn++)
   {
      sig = (X[nn] - X[nn-1]) / (X[nn+1] - X[nn-1]);
      p = sig * YY[nn-1] + 2.0;
      YY[nn] = (sig - 1.0) / p;
      ddata = (Y[nn+1] - Y[nn])/(X[nn+1] - X[nn]) -
              (Y[nn] - Y[nn-1])/(X[nn] - X[nn-1]);
      U[nn] = (6.0 * ddata / (X[nn+1]-X[nn-1]) - sig*U[nn-1]) / p;
   }
   for (nn = leng-2; nn >= 0; nn--)
      YY[nn] = YY[nn] * YY[nn+1] + U[nn];
   free(U);
   return 0;
}

/********************************************************************
 * one dimensional spline interpolation
 * ****************************************************************/
double splint(int leng,double *X,double *Y,double *Ypp,double Xval)
{
   int    klo, khi, k;
   double h, Yval, a, b;

   /* search interval */
   klo = 0;
   khi = leng - 1;
   while (khi - klo > 1)
   {
      k = (khi + klo) /2;
      if (X[k] > Xval) khi = k;
      else             klo = k;
   }
   h = X[khi] - X[klo];
   if (h == 0)
   {
      printf("ERROR: bad input for splint %e\n",h);
      printf("The following array should be in ascending order\n");
      for (k = 0; k < leng; k++)
      {
         printf("%e\n", X[k]);
      }
      exit(1);
   }
   a = (X[khi] - Xval) / h;
   b = (Xval - X[klo]) / h;
   Yval = a * Y[klo] + b * Y[khi] + ((a*a*a-a)*Ypp[klo]+(b*b*b-b)*Ypp[khi])*
       (h*h)/6.0;
   return Yval;
}

/********************************************************************
 * one dimensional spline 
 *******************************************************************/
int spline1d(int nx, double *Xn, double *Zn,
             int mx, double *Xm, double *Zm)
{
   int    ii;
   double *ZP, Xval;
   /* ======================================= */
   /* compute second derivatives              */
   /* ======================================= */
   ZP = (double *) malloc(nx * sizeof(double));
   spline1p(nx, Xn, Zn, ZP);
   for (ii = 0; ii < mx; ii++)
   {
      Xval = Xm[ii];
      Zm[ii] = splint(nx, Xn, Zn, ZP, Xval);
   }
   free(ZP);
   return 0;
}

/********************************************************************
 * two dimensional spline 
 *******************************************************************/
int spline2d(int nx, int ny, double *Xn, double *Yn, double *Zn,
             int mm, double *Xm, double *Ym, double *Zm)
{
   int    offset, ii, xx, yy, nn;
   double *ZT, *ZP, *ZU, Xval, Yval;
   /* ======================================= */
   /* compute second derivatives              */
   /* compute spline one row at a time        */
   /* ======================================= */
   nn = nx;
   if (nn < ny) nn = ny;
   ZT = (double *) malloc(nn * sizeof(double));
   ZP = (double *) malloc(nx * ny * sizeof(double));
   ZU = (double *) malloc(nx * sizeof(double));
   for (xx = 0; xx < nx; xx++)
   {
      offset = xx * ny;
      for (yy = 0; yy < ny; yy++) ZT[yy] = Zn[yy*nx+xx];
      spline1p(ny, Yn, ZT, &ZP[offset]);
   }

   /* ======================================= */
   /* interpolate                             */
   /* ======================================= */
   for (ii = 0; ii < mm; ii++)
   {
      Xval = Xm[ii];
      Yval = Ym[ii];
      for (xx = 0; xx < nx; xx++)
      {
         offset = xx * ny;
         for (yy = 0; yy < ny; yy++) ZT[yy] = Zn[yy*nx+xx];
         ZU[xx] = splint(ny, Yn, ZT, &ZP[offset], Yval);
      }
      spline1p(nx, Xn, ZU, ZT);
      Zm[ii] = splint(nx, Xn, ZU, ZT, Xval);
   }
   free(ZT);
   free(ZP);
   free(ZU);
   return 0;
}

/********************************************************************
 * three dimensional spline 
 *******************************************************************/
int spline3d(int nx, int ny,int nz,double *Xn,double *Yn,double *Zn,
             double *Un,int mm,double *Xm,double *Ym,double *Zm, 
             double *Um)
{
   int    offset, ii, xx, yy, zz;
   double *UT, *UD, *UU, *UV, Xval, Yval, Zval;
   /* ======================================= */
   /* allocate memory                         */
   /* ======================================= */
   UD = (double *) malloc(nx * ny * nz * sizeof(double));
   UT = (double *) malloc((nx + ny + nz) * sizeof(double));
   UU = (double *) malloc(nx * ny * sizeof(double));
   UV = (double *) malloc(nx * sizeof(double));
   /* ======================================= */
   /* compute second derivatives in z-axis    */
   /* results stored in UD                    */
   /* ======================================= */
   for (xx = 0; xx < nx; xx++)
   {
      for (yy = 0; yy < ny; yy++)
      {
         offset = xx * ny * nz + yy * nz;
         for (zz = 0; zz < nz; zz++) UT[zz] = Un[zz*nx*ny+yy*nx+xx];
         spline1p(nz, Zn, UT, &UD[offset]);
      }
   }

   /* ======================================= */
   /* interpolate                             */
   /* ======================================= */
   for (ii = 0; ii < mm; ii++)
   {
      Zval = Zm[ii];
      /* use 2nd derivatives in z-axis to interpolate z */
      /* at each of points in the (x,y) plane => UU     */
      for (xx = 0; xx < nx; xx++)
      {
         for (yy = 0; yy < ny; yy++)
         {
            offset = xx * ny * nz + yy * nz;
            for (zz = 0; zz < nz; zz++) UT[zz] = Un[zz*nx*ny+yy*nx+xx];
            UU[xx*ny+yy] = splint(nz,Zn,UT,&UD[offset],Zval);
         }
         /* With a target Z value, UU has the interpolated */
         /* values at all points in the y-axis at a given  */
         /* x point in the original image. Now interpolate */
         /* at the target Y value. First compute second    */
         /* derivatives ==> UT                             */
         spline1p(ny, Yn, &UU[xx*ny], UT);
         /* now interpolate at the given Y for all points  */
         /* in the x axis ==> UV                           */
         Yval = Ym[ii];
         UV[xx] = splint(ny, Yn, &UU[xx*ny], UT, Yval);
      }
      /* now compute second derivative in the x-axis */
      spline1p(nx, Xn, UV, UT);
      Xval = Xm[ii];
      /* finally interpolate the point at the target Xval */
      Um[ii] = splint(nx, Xn, UV, UT, Xval);
   }
   free(UT);
   free(UD);
   free(UU);
   free(UV);
   return 0;
}

