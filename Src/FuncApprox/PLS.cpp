// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class PLS
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sysdef.h"
#include "Psuade.h"
#include "PLS.h"
#include "PDFManager.h"
#include "PsuadeUtil.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

extern "C"
{
   void dgels_(char *, int *, int *, int *, double *, int *,
               double *, int *, double *, int *, int *);
   void dgetrf_(int *, int *, double *, int *, int *, int *);
   void dgetri_(int *, double *, int *, int *, double*, int *, int *);
}

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
PLS::PLS(int nInputs,int nSamples):FuncApprox(nInputs,nSamples)
{
   faID_ = PSUADE_RS_PLS;
   numTerms_  = 0;
   regCoeffs_ = NULL;
   initialized_ = 0;
   regStdevs_ = NULL;
   fuzzyC_    = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PLS::~PLS()
{
   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   if (regStdevs_ != NULL) delete [] regStdevs_;
   if (fuzzyC_ != NULL)
   {
      for (int ii = 0; ii < numTerms_; ii++) delete [] fuzzyC_[ii];
      delete [] fuzzyC_;
   }
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int PLS::initialize(double *X, double *Y)
{
   if (outputLevel_ >= 0)
   {
      printAsterisks(PL_INFO, 0);
      printf("*     Partial Least Squares Linear Regression Analysis\n");
      printf("* R-squared gives a measure of the goodness of the model.\n");
      printf("* R-squared should be close to 1 if it is a good model.\n");
      printf("* TURN ON rs_expert mode to output orthonormal matrix.\n");
      printDashes(PL_INFO, 0);
   }

   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("PLS initialize ERROR: consult PSUADE developers.\n");
      exit( 1 );
   } 

   if (regCoeffs_ != NULL) delete [] regCoeffs_;
   if (regStdevs_ != NULL) delete [] regStdevs_;
   regCoeffs_ = NULL;
   regStdevs_ = NULL;
   if (fuzzyC_ != NULL)
   {
      for (int ii = 0; ii < nInputs_; ii++) delete [] fuzzyC_[ii];
      delete [] fuzzyC_;
      fuzzyC_ = NULL;
   }
   return analyze(X,Y);
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int PLS::genNDGridData(double *X, double *Y, int *NN, double **XX, 
                       double **YY)
{
   int mm, totPts, status=0;

   if (initialized_ == 0)
   {
      status = initialize(X,Y);
      if (status != 0)
      {
         printf("PLS: ERROR detected in regression analysis.\n");
         (*NN) = 0;
         return -1;
      }
   }

   if ((*NN) == -999) return 0;

   genNDGrid(NN, XX);
   if ((*NN) == 0) return 0;
   totPts = (*NN);

   (*YY) = new double[totPts];
   checkAllocate(*YY, "YY in PLS::genNDGridData");

   (*NN) = totPts;
   for (mm = 0; mm < totPts; mm++)
      (*YY)[mm] = evaluatePoint(&((*XX)[mm*nInputs_]));

   return 0;
}

// ************************************************************************
// Generate 1D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
                       int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, status=0;
   double HX, *Xloc;

   if (initialized_ == 0)
   {
      status = initialize(X,Y);
      if (status != 0)
      {
         printf("PLS: ERROR detected in regression analysis.\n");
         (*NN) = 0;
         return -1;
      }
   }

   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PLS::gen1DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      Xloc[ind1] = HX * mm + lowerBounds_[ind1];
      (*XX)[mm] = Xloc[ind1];
      (*YY)[mm] = evaluatePoint(Xloc);
   }

   delete [] Xloc;
   return 0;
}

// ************************************************************************
// Generate 2D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                       double *settings, int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, ind, status=0;
   double *HX, *Xloc;

   if (initialized_ == 0)
   {
      status = initialize(X,Y);
      if (status != 0)
      {
         printf("PLS: ERROR detected in regression analysis.\n");
         (*NN) = 0;
         return -1;
      }
   }

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 2];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PLS::gen2DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         ind = mm * nPtsPerDim_ + nn;
         Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
         Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
         (*XX)[ind*2]   = Xloc[ind1];
         (*XX)[ind*2+1] = Xloc[ind2];
         (*YY)[ind] = evaluatePoint(Xloc);
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3, 
                       double *settings, int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, pp, ind, status=0;
   double *HX, *Xloc;

   if (initialized_ == 0)
   {
      status = initialize(X,Y);
      if (status != 0)
      {
         printf("PLS: ERROR detected in regression analysis.\n");
         (*NN) = 0;
         return -1;
      }
   }

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 3];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PLS::gen3DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            ind = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
            Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
            Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
            Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
            (*XX)[ind*3]   = Xloc[ind1];
            (*XX)[ind*3+1] = Xloc[ind2];
            (*XX)[ind*3+2] = Xloc[ind3];
            (*YY)[ind] = evaluatePoint(Xloc);
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int PLS::gen4DGridData(double *X, double *Y, int ind1, int ind2,
                              int ind3, int ind4, double *settings, 
                              int *NN, double **XX, double **YY)
{
   int    totPts, mm, nn, pp, qq, ind, status=0;
   double *HX, *Xloc;

   if (initialized_ == 0)
   {
      status = initialize(X,Y);
      if (status != 0)
      {
         printf("PLS: ERROR detected in regression analysis.\n");
         (*NN) = 0;
         return -1;
      }
   }
 
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 4];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PLS::gen4DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
   for (mm = 0; mm < nPtsPerDim_; mm++) 
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            for (qq = 0; qq < nPtsPerDim_; qq++)
            {
               ind = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                       nn*nPtsPerDim_*nPtsPerDim_ + pp*nPtsPerDim_ + qq;
               Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
               Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
               Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
               Xloc[ind4] = HX[3] * qq + lowerBounds_[ind4];
               (*XX)[ind*4]   = Xloc[ind1];
               (*XX)[ind*4+1] = Xloc[ind2];
               (*XX)[ind*4+2] = Xloc[ind3];
               (*XX)[ind*4+3] = Xloc[ind4];
               (*YY)[ind] = evaluatePoint(Xloc);
            }
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double PLS::evaluatePoint(double *X)
{
   int    mm;
   double Xdata, Y=0.0;

   if (regCoeffs_ == NULL)
   {
      printf("PLS ERROR: need to call initialize first.\n");
      return 0.0;
   }
   for (mm = 0; mm < nInputs_; mm++)
   {
      Xdata = (X[mm] - XMeans_[mm]) / XStds_[mm]; 
      Y += regCoeffs_[mm] * Xdata;
   }
   Y = Y * YStd_ + YMean_;
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double PLS::evaluatePoint(int npts, double *X, double *Y)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePoint(&X[kk*nInputs_]);
   return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double PLS::evaluatePointFuzzy(double *X, double &std)
{
   int    mm, cc, nTimes=100;
   double accum, *Ys, coef, mean, stds, Xdata;

   if (regCoeffs_ == NULL)
   {
      printf("PLS ERROR: initialize has not been called.\n");
      exit(1);
   }
 
   Ys = new double[nTimes];
   checkAllocate(Ys, "Ys in PLS::evaluatePointFuzzy");

   mean = 0.0;
   for (cc = 0; cc < nTimes; cc++)
   {
      accum = 0.0;
      for (mm = 0; mm < nInputs_; mm++)
      {
         coef = fuzzyC_[mm][cc];
         Xdata = (X[mm] - XMeans_[mm]) / XStds_[mm]; 
         accum += coef * Xdata;
      }
      accum = accum * YStd_ + YMean_;
      Ys[cc] = accum;
      mean += accum;
   }
   mean /= (double) nTimes;
   stds = 0.0;
   for (cc = 0; cc < nTimes; cc++)
      stds += (Ys[cc] - mean) * (Ys[cc] - mean);
   stds = sqrt(stds / (nTimes - 1));
   delete [] Ys;
   std = stds; 
   return mean;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double PLS::evaluatePointFuzzy(int npts, double *X, double *Y,double *Ystd)
{
   if (regCoeffs_ == NULL)
   {
      printf("PLS ERROR: initialize has not been called.\n");
      exit(1);
   }
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePointFuzzy(&(X[kk*nInputs_]), Ystd[kk]);
   return 0.0;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int PLS::analyze(double *Xin, double *Y)
{
   int    ii, jj, kk, nn, info, status=0;
   double ddata, snorm, *P, *T, *R, *S, *B, *X, *WW, *D, esum, ymax;
   double SSresid, SStotal, R2, *XTX, var, *Bstd;
   
   X = new double[nSamples_*nInputs_];
   //Y = new double[nSamples_];
   initInputScaling(Xin, X, 1);
   //initOutputScaling(Yin, Y);
   ddata = 0.0;
   for (ii = 0; ii < nSamples_; ii++) ddata += Y[ii];
   YMean_ = ddata / (double) nSamples_;
   for (ii = 0; ii < nSamples_; ii++) Y[ii] -= YMean_;
   YStd_ = 1.0;

   S = new double[nInputs_];
   checkAllocate(S, "S in PLS::analyze");
   for (ii = 0; ii < nInputs_; ii++)
   {
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
         ddata += X[jj*nInputs_+ii] * Y[jj];
      S[ii] = ddata;
   }
   for (ii = 0; ii < nSamples_; ii++) Y[ii] += YMean_;

   snorm = 0.0;
   for (jj = 0; jj < nInputs_; jj++) snorm += S[jj] * S[jj]; 
   snorm = sqrt(snorm);
   if (snorm == 0.0) 
   {
      printf("PLS ERROR: null correlation matrix.\n");
      status = 1;
      delete [] S;
      delete [] X;
      return status;
   }
   R = new double[nInputs_*nInputs_];
   T = new double[nSamples_*nInputs_];
   P = new double[nInputs_*nInputs_];
   D = new double[nInputs_];
   checkAllocate(D, "D in PLS::analyze");
   for (ii = 0; ii < nInputs_; ii++)
   {
      for (jj = 0; jj < nInputs_; jj++)
         R[ii*nInputs_+jj] = S[jj] / snorm; 
      for (kk = 0; kk < nSamples_; kk++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nInputs_; jj++)
            ddata += X[jj+kk*nInputs_] * R[ii*nInputs_+jj];
         T[ii*nSamples_+kk] = ddata;
      }
      for (kk = 0; kk < nInputs_; kk++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nSamples_; jj++)
            ddata += X[kk+jj*nInputs_] * T[ii*nSamples_+jj];
         P[ii*nInputs_+kk] = ddata;
      }
      for (jj = 0; jj < ii; jj++)
      {
         ddata = 0.0;
         for (kk = 0; kk < nInputs_; kk++)
            ddata = ddata + P[ii*nInputs_+kk] * P[jj*nInputs_+kk];
         for (kk = 0; kk < nInputs_; kk++)
            P[ii*nInputs_+kk] -= ddata * P[jj*nInputs_+kk];
      }
      ddata = 0.0;
      for (jj = 0; jj < nInputs_; jj++)
         ddata += P[ii*nInputs_+jj] * P[ii*nInputs_+jj]; 
      ddata = sqrt(ddata);
      if (ddata == 0.0) 
      {
         printf("PLS ERROR: rank deficient P matrix.\n");
         status = 1;
         break;
      }
      for (jj = 0; jj < nInputs_; jj++) P[ii*nInputs_+jj] /= ddata; 
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
         ddata += T[ii*nSamples_+jj] * T[ii*nSamples_+jj]; 
      ddata = sqrt(ddata);
      if (ddata == 0.0) 
      {
         printf("PLS ERROR: null T vector.\n");
         status = 1;
         break;
      }
      for (jj = 0; jj < nSamples_; jj++) T[ii*nSamples_+jj] /= ddata; 
      for (jj = 0; jj < nInputs_; jj++) R[ii*nInputs_+jj] /= ddata; 
      D[ii] = ddata;
      ddata = 0.0; 
      for (jj = 0; jj < nInputs_; jj++) 
         ddata += P[ii*nInputs_+jj] * S[jj];
      for (jj = 0; jj < nInputs_; jj++) 
         S[jj] = S[jj] - P[ii*nInputs_+jj] * ddata;

      snorm = 0.0;
      for (jj = 0; jj < nInputs_; jj++) snorm += S[jj] * S[jj]; 
      snorm = sqrt(snorm);
      if (snorm < 1.0e-8) break;
   }
   if (ii < nInputs_) ii++;
   numTerms_ = ii;
   snorm = 0.0;
   for (jj = 0; jj < nSamples_; jj++)
   {
      for (ii = 0; ii < nInputs_; ii++)
      {
         for (kk = 0; kk < numTerms_; kk++)
            ddata = T[kk*nSamples_+jj] * P[kk*nInputs_+ii] * D[kk];
         ddata = ddata - X[jj*nInputs_+ii];
         snorm += ddata * ddata;
      }
   }
   snorm = sqrt(snorm);
   printf("PLS residual norm (||X-TDP^t||) = %e\n", snorm);
#if 0
   FILE *fp = fopen("pls.m", "w");
   fprintf(fp, "%% X = T D P^t + E\n");
   fprintf(fp, "X = [\n");
   for (jj = 0; jj < nSamples_; jj++)
   {
      for (ii = 0; ii < numTerms_; ii++)
         fprintf(fp, "%24.16e ", X[jj*nInputs_+ii]);
      fprintf(fp, "\n");
   }
   fprintf(fp, "];\n");
   fprintf(fp, "T = [\n");
   for (jj = 0; jj < nSamples_; jj++)
   {
      for (ii = 0; ii < numTerms_; ii++)
         fprintf(fp,"%24.16e ", T[ii*nSamples_+jj]);
      fprintf(fp, "\n");
   }
   fprintf(fp, "];\n");
   fprintf(fp, "R = [\n");
   for (jj = 0; jj < nInputs_; jj++)
   {
      for (ii = 0; ii < numTerms_; ii++)
         fprintf(fp,"%24.16e ", R[ii*nInputs_+jj]);
      fprintf(fp, "\n");
   }
   fprintf(fp, "];\n");
   fprintf(fp, "P = [\n");
   for (jj = 0; jj < nInputs_; jj++)
   {
      for (ii = 0; ii < numTerms_; ii++)
         fprintf(fp,"%24.16e ", P[ii*nInputs_+jj]);
      fprintf(fp, "\n");
   }
   fprintf(fp, "];\n");
   fprintf(fp, "D = [\n");
   for (ii = 0; ii < numTerms_; ii++)
      fprintf(fp,"%24.16e ", D[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "D1 = diag(D);\n");
   fprintf(fp, "X1 = T * D1 * P';\n");
   fprintf(fp, "sqrt(sum(sum((X1-X)^.2)))\n");
   fclose(fp);
#endif
   if (status == 1) 
   {
      printf("PLS ERROR: failed.\n");
      status = 1;
      delete [] P;
      delete [] R;
      delete [] T;
      delete [] S;
      delete [] X;
      delete [] D;
      //delete [] Y;
      return status;
   }

   if (outputLevel_ > 1)
   {
      printf("Basis set T (so that X = T D P') is:\n");
      for (jj = 0; jj < nSamples_; jj++)
      {
         for (ii = 0; ii < numTerms_; ii++)
            printf("%12.4e ", T[ii*nSamples_+jj]);
         printf("\n");
      }
      printf("Given new x, transform to reduced space via x*P*D\n");
      printf("where x is 1xn\n");
      printf("      P is nxp  (p << n is the reduced input size)\n");
      printf("      D is pxp\n");
   }
   B = new double[nInputs_];
   checkAllocate(B, "B in PLS::analyze");
   for (ii = 0; ii < numTerms_; ii++)
   {
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
         ddata += T[jj+ii*nSamples_] * Y[jj];
      P[ii] = ddata;
   }
   for (ii = 0; ii < nInputs_; ii++)
   {
      ddata = 0.0;
      for (jj = 0; jj < nInputs_; jj++)
         ddata += R[jj*nInputs_+ii] * P[jj];
      B[ii] = ddata;
   }
   initialized_ = 1;
   
   delete [] P;
   delete [] R;
   delete [] T;
   delete [] S;
   delete [] D;

   WW = new double[nSamples_];
   checkAllocate(WW, "WW in PLS::analyze");
   esum = ymax = 0.0;
   for (ii = 0; ii < nSamples_; ii++)
   {
      WW[ii] = 0.0;
      for (jj = 0; jj < nInputs_; jj++) WW[ii] += X[ii*nInputs_+jj]*B[jj];
      WW[ii] -= Y[ii];
      esum = esum + WW[ii] * WW[ii];
      if (PABS(Y[ii]) > ymax) ymax = PABS(Y[ii]);
   }
   esum /= (double) nSamples_;
   esum = sqrt(esum);
   if (outputLevel_ > 1)
      printf("* PLS: rms of interpolation error = %11.4e (Ymax=%9.2e)\n",
             esum, ymax); 

   computeSS(numTerms_, X, Y, B, SSresid, SStotal);
   R2 = 1.0;
   if (SStotal != 0.0) R2  = 1.0 - SSresid / SStotal;
   if (nSamples_ > numTerms_)
      var = SSresid / (double) (nSamples_ - numTerms_);
   else               var = 0.0;
   if (var < 0)
   { 
      if (PABS(var) > 1.0e-12)
           printf("PLS WARNING: var < 0.\n");
      else var = 0;
   }

   Bstd = new double[nInputs_];
   checkAllocate(Bstd, "Bstd in PLS::analyze");
   computeXTX(nInputs_, X, &XTX);
   computeCoeffVariance(nInputs_, XTX, var, Bstd);
   regCoeffs_ = B;
   regStdevs_ = Bstd;

   PDFManager *pdfman = new PDFManager();
   int    cc, nTimes=100;
   int    *inPDFs = new int[nInputs_];
   double *inMeans = new double[nInputs_];
   double *inStds = new double[nInputs_];
   double *inUppers = new double[nInputs_];
   double *inLowers = new double[nInputs_];
   checkAllocate(inLowers, "inLowers in PLS::analyze");
   for (nn = 0; nn < nInputs_; nn++)
   {
      inPDFs[nn] = PSUADE_PDF_NORMAL;
      inMeans[nn] = regCoeffs_[nn];
      inStds[nn] = regStdevs_[nn];
      inUppers[nn] = inMeans[nn] + 4.0 * inStds[nn];
      inLowers[nn] = inMeans[nn] - 4.0 * inStds[nn];
      if (inUppers[nn] == inLowers[nn])
      {
         if (inUppers[nn] > 0) inUppers[nn] *= (1.0 + 1.0e-14);
         else                  inUppers[nn] *= (1.0 - 1.0e-14);
         if (inLowers[nn] > 0) inLowers[nn] *= (1.0 - 1.0e-14);
         else                  inLowers[nn] *= (1.0 + 1.0e-14);
         if (inUppers[nn] == 0.0) inUppers[nn] = 1e-14;
         if (inLowers[nn] == 0.0) inLowers[nn] = -1e-14;
      }
   }
   pdfman->initialize(nInputs_,inPDFs,inMeans,inStds,covMatrix_,NULL,NULL);
   psVector vLower, vUpper, vOut;
   vLower.load(nInputs_, inLowers);
   vUpper.load(nInputs_, inUppers);
   vOut.setLength(nInputs_*nTimes);
   pdfman->genSample(nTimes, vOut, vLower, vUpper);
   fuzzyC_ = new double*[nInputs_];
   checkAllocate(fuzzyC_, "fuzzyC_ in PLS::analyze");
   for (nn = 0; nn < nInputs_; nn++)
   {
      fuzzyC_[nn] = new double[nTimes];
      checkAllocate(fuzzyC_[nn], "fuzzyC_[nn] in PLS::analyze");
      for (cc = 0; cc < nTimes; cc++)
         fuzzyC_[nn][cc] = vOut[cc*nInputs_+nn];
   }
   delete pdfman;
   delete [] inPDFs;
   delete [] inStds;
   delete [] inMeans;
   delete [] inLowers;
   delete [] inUppers;

   if (outputLevel_ >= 0)
   {
      printRC(nInputs_, B, Bstd, X, Y);
      printf("* PLS R-square = %12.4e (SSresid, SStotal =%10.2e,%10.2e)\n",
             R2, SSresid, SStotal);
      if ((nSamples_ - numTerms_ - 1) > 0)
         printf("* Adjusted R2  = %12.4e\n",
                1.0-(1.0-R2)*((nSamples_-1)/(nSamples_-numTerms_-1)));
      if (SSresid > 0)
         printf("* F-statistics = %12.4e (>6 for fit to be significant)\n",
                (SStotal-SSresid)/SSresid*(nSamples_-numTerms_));
      printEquals(PL_INFO, 0);
   }

   delete [] WW;
   delete [] X;
   //delete [] Y;
   delete [] XTX;
   return 0;
}

// *************************************************************************
// compute SS (sum of squares) statistics
// -------------------------------------------------------------------------
int PLS::computeSS(int N, double *XX, double *Y, double *B, double &SSresid, 
                   double &SStotal)
{
   int    nn, mm;
   double rdata, ymean, SSreg, ddata;

   SSresid = SStotal = SSreg = ymean = 0.0;
   for (mm = 0; mm < nSamples_; mm++) ymean += Y[mm];
   ymean /= (double) nSamples_;
   for (mm = 0; mm < nSamples_; mm++)
   {
      ddata = 0.0;
      for (nn = 0; nn < N; nn++) ddata += (XX[mm*N+nn] * B[nn]);
      ddata += YMean_;
      rdata = Y[mm] - ddata;
      SSresid += rdata * rdata;
      SSreg += (ddata - ymean) * (ddata - ymean);
   }
   for (mm = 0; mm < nSamples_; mm++)
      SStotal += (Y[mm] - ymean) * (Y[mm] - ymean);
   if (outputLevel_ > 0)
   {
      printf("* PLS: SStot  = %16.8e\n", SStotal);
      printf("* PLS: SSres  = %16.8e\n", SSresid);
      printf("* PLS: SSreg  = %16.8e\n", SSreg);
   }
   if (outputLevel_ > 0 && nSamples_ != N)
      printf("* PLS: MSE of residual = %16.8e\n",SSresid/(nSamples_-N));
   return 0;
}

// *************************************************************************
// compute coefficient variances (diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int PLS::computeCoeffVariance(int N,double *XX,double var,double *B)
{
   int    nn, nn2, lwork, iOne=1, info, errCnt=0;
   double *B2, *work, *XT;
   char   trans[1];

   (*trans) = 'N';
   B2 = new double[N];
   XT = new double[N*N];
   lwork = 2 * N * N;
   work  = new double[lwork];
   checkAllocate(work, "work in PLS::computeCoeffVariance");
   for (nn = 0; nn < N; nn++)
   {
      for (nn2 = 0; nn2 < N*N; nn2++) XT[nn2] = XX[nn2];
      for (nn2 = 0; nn2 < N; nn2++) B2[nn2] = 0.0;
      B2[nn] = var;
      dgels_(trans, &N, &N, &iOne, XT, &N, B2, &N, work, &lwork, &info);
      if (info != 0)
         printf("PLS WARNING: dgels returns error %d.\n",info);
      if (B2[nn] < 0) errCnt++;
      if (B2[nn] < 0) B[nn] = sqrt(-B2[nn]);
      else            B[nn] = sqrt(B2[nn]);
   }
   if (errCnt > 0)
   {
      printf("* PLS WARNING: some of the coefficient variances\n");
      printf("*              are < 0. May spell trouble but will\n");
      printf("*              proceed anyway (%d).\n", errCnt);
   }
   delete [] B2;
   delete [] XT;

   int    *ipiv = new int[N+1];
   double *invA = new double[lwork];
   double ddata, ddata2;
   checkAllocate(invA, "invA in PLS::computeCoeffVariance");
   for (nn = 0; nn < N*N; nn++) invA[nn] = XX[nn];
   dgetrf_(&N, &N, invA, &N, ipiv, &info);
   if (info != 0)
      printf("PLS WARNING: dgels returns error %d.\n",info);
   dgetri_(&N, invA, &N, ipiv, work, &lwork, &info);
   covMatrix_.setDim(N,N);
   for (nn = 0; nn < N; nn++)
   {
      for (nn2 = 0; nn2 < N; nn2++)
      {
         ddata = invA[nn*N+nn2] * var;
         covMatrix_.setEntry(nn,nn2,ddata);
      }
   }
   for (nn = 0; nn < N; nn++)
   {
      ddata = covMatrix_.getEntry(nn,nn);
      ddata = sqrt(ddata);
      for (nn2 = 0; nn2 < N; nn2++)
      {
         if (nn != nn2)
         {
            ddata2 = covMatrix_.getEntry(nn,nn2);
            if (ddata != 0) ddata2 /= ddata;
            covMatrix_.setEntry(nn,nn2,ddata2);
         }
      }
   }
   for (nn2 = 0; nn2 < N; nn2++)
   {
      ddata = covMatrix_.getEntry(nn2,nn2);
      ddata = sqrt(ddata);
      for (nn = 0; nn < N; nn++)
      {
         if (nn != nn2)
         {
            ddata2 = covMatrix_.getEntry(nn,nn2);
            if (ddata != 0) ddata2 /= ddata;
            covMatrix_.setEntry(nn,nn2,ddata2);
         }
      }
   }
   ddata = 1.0;
   for (nn = 0; nn < N; nn++) covMatrix_.setEntry(nn,nn,ddata);
   for (nn = 0; nn < N; nn++)
   {
      for (nn2 = 0; nn2 < nn; nn2++)
      {
         ddata  = covMatrix_.getEntry(nn,nn2);
         ddata2 = covMatrix_.getEntry(nn2,nn);
         ddata  = 0.5 * (ddata + ddata2);
         covMatrix_.setEntry(nn,nn2,ddata);
         covMatrix_.setEntry(nn2,nn,ddata);
      }
   }
   errCnt = 0;
   for (nn = 0; nn < N; nn++)
   {
      for (nn2 = 0; nn2 < N; nn2++)
      {
         ddata = covMatrix_.getEntry(nn,nn2);
         if (nn != nn2 && (ddata >=1 || ddata <= -1))
         {
            errCnt++;
            covMatrix_.setEntry(nn,nn2,0.0);
         }
      }
   }
   char inStr[1001];
   if (errCnt > 0)
   {
      printf("PLS WARNING:\n");
      printf("  Correlation matrix has invalid entries (%d out of %d).\n",
             errCnt, N*(N-1));
      printf("  EVALUATION MAY BE INCORRECT.\n");
      printf("  CONTINUE ANYWAY (will set them to zeros)? (y or n)");
      scanf("%s", inStr);
      if (inStr[0] != 'y') exit(1);
      fgets(inStr, 100, stdin);
   }
   delete [] work;
   delete [] ipiv;
   delete [] invA;
   return info;
}

// *************************************************************************
// print statistics
// -------------------------------------------------------------------------
int PLS::printRC(int N,double *B,double *Bvar,double *XX,double *Y)
{
   int    nn1;
   double coef, scaled;

   printEquals(PL_INFO, 0);
   scaled = 0.0;
   for (nn1 = 0; nn1 < nInputs_; nn1++) scaled += PABS(XMeans_[nn1]);
   if (scaled == 0)
   {
      for (nn1 = 0; nn1 < nInputs_; nn1++) 
         if (XStds_[nn1] != 1) break;
      if (nn1 != nInputs_) scaled = 1;
   }
   if (scaled != 0)
   {
      printf("* NOTE: The coefficients below have been scaled.\n");
      printf("*       See above for input scaling information.\n");
   }
   printDashes(PL_INFO, 0);
   printf("*            ");
   printf("  coefficient   std. error   t-value\n");
   printDashes(PL_INFO, 0);
   printf("* Constant   = %12.4e\n", YMean_);
   for (nn1 = 0; nn1 < nInputs_; nn1++)
   {
      if (PABS(Bvar[nn1]) < 1.0e-15) coef = 0.0;
      else                           coef = B[nn1] / Bvar[nn1]; 
      {
         printf("* Input  %3d ", nn1+1);
         printf("= %12.4e %12.4e %12.4e\n", B[nn1], Bvar[nn1], coef);
      }
   }
   printEquals(PL_INFO, 0);
   return 0;
}

// *************************************************************************
// print coefficients 
// -------------------------------------------------------------------------
int PLS::printCoefs(int N, double *B)
{
   int    ii, jj, nn1, **indexTable;
   int    *indices;
   double *trueCoefs, Bmax;

   indexTable = new int*[N];
   for (ii = 0; ii < N; ii++) indexTable[ii] = new int[nInputs_];
   checkAllocate(indexTable[N-1], "indexTable in PLS::printCoefs");
   for (ii = 0; ii < nInputs_; ii++) indexTable[0][ii] = 0; 
   for (ii = 0; ii < nInputs_; ii++) 
   {
      for (jj = 0; jj < nInputs_; jj++) 
         indexTable[ii+1][jj] = 0; 
      indexTable[ii+1][ii] = 1; 
   }
   for (ii = 0; ii < N; ii++)
   {
      for (jj = 0; jj < nInputs_; jj++)
         printf("%d ", indexTable[ii][jj]);
      printf("\n");
   }
   printEquals(PL_INFO, 0);
   printf("*** Note: these coefficients are true coefficients.\n");
   printDashes(PL_INFO, 0);
   printf("*            ");
   printf("  coefficient\n");
   printDashes(PL_INFO, 0);

   trueCoefs = new double[N];
   for (ii = 0; ii < N; ii++) trueCoefs[ii] = 0.0;
   trueCoefs[0] = B[0];
   indices = new int[nInputs_];
   checkAllocate(indices, "indices in PLS::printCoefs");
   for (nn1 = 1; nn1 <= nInputs_; nn1++)
   {
      trueCoefs[nn1] = B[nn1] / XStds_[nn1-1];
      trueCoefs[0]  -= B[nn1] * XMeans_[nn1-1] / XStds_[nn1-1];
   }
   Bmax = trueCoefs[0];
   for (nn1 = 1; nn1 < N; nn1++)
      if (PABS(trueCoefs[nn1]) > Bmax) 
         Bmax = PABS(trueCoefs[nn1]);
   if (Bmax == 0) Bmax = 1;
   for (nn1 = 0; nn1 < N; nn1++)
      if (PABS(trueCoefs[nn1]/Bmax) < 1.0e-8)
         trueCoefs[nn1] = 0;

   for (nn1 = 0; nn1 < nInputs_; nn1++)
   {
      if (trueCoefs[nn1] != 0)
         printf("* Input %3d = %16.8e \n", nn1, trueCoefs[nn1]);
   }
   printDashes(PL_INFO, 0);
   delete [] trueCoefs;
   for (ii = 0; ii < N; ii++) delete [] indexTable[ii];
   delete [] indexTable;
   return 0;
}

// *************************************************************************
// form X^T X 
// -------------------------------------------------------------------------
int PLS::computeXTX(int N, double *X, double **XXOut)
{
   int    nn, nn2, mm;
   double *XX, coef;

   XX = new double[nSamples_*N];
   checkAllocate(XX, "XX in PLS::computeXTX");
   for (nn = 0; nn < N; nn++)
   {
      for (nn2 = nn; nn2 < N; nn2++)
      {
         coef = 0.0;
         for (mm = 0; mm < nSamples_; mm++)
            coef += X[nn+mm*N] * X[nn2+mm*N];
         XX[nn*N+nn2] = coef;
         XX[nn2*N+nn] = coef;
      }
   }
   (*XXOut) = XX;
   return 0;
}


