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
// Functions for the class NPLearning
// AUTHOR : CHARLES TONG
// DATE   : 2013
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "NPLearning.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"

extern "C"
{
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
}

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class NPLearning
// ------------------------------------------------------------------------
NPLearning::NPLearning(int nInputs,int nSamples) : 
                       FuncApprox(nInputs,nSamples)
{
   char pString[501];

   faID_ = PSUADE_RS_NPL;

   depth1_  = 25;
   depth2_  = 6;
   epsilon_ = 1.0e-4;

   if (psRSExpertMode_ == 1)
   {
      sprintf(pString,"Enter the desired number of levels (>10): ");
      depth1_ = getInt(1, 10, pString);
      sprintf(pString,"Enter the desired number of levels (>10): ");
      depth2_ = getInt(1, 10, pString);
      epsilon_ = 0.0;
      sprintf(pString,"Enter the tolerance (> 0, < 1): ");
      if (epsilon_ <= 0 || epsilon_ >= 1)
         epsilon_ = getDouble(pString);
   }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
NPLearning::~NPLearning()
{
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int NPLearning::initialize(double *X, double *Y)
{
   learn(X, Y);
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int NPLearning::genNDGridData(double *X, double *Y, int *N2, double **X2, 
                              double **Y2)
{
   int totPts;

   learn(X, Y);

   if ((*N2) == -999) return 0;
  
   genNDGrid(N2, X2);
   if ((*N2) == 0) return 0;
   totPts = (*N2);

   (*Y2) = new double[totPts];
   evaluatePoint(totPts, *X2, *Y2);

   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int NPLearning::gen1DGridData(double *X, double *Y, int ind1, 
                              double *settings, int *N, double **X2, 
                              double **Y2)
{
   int    ii, ss, totPts;
   double *XT, *XX, *YY, HX;

   learn(X, Y);

   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ss = 0; ss < totPts; ss++) 
   {
      XT[ss*nInputs_+ind1]  = HX * ss + lowerBounds_[ind1];
      XX[ss] = HX * ss + lowerBounds_[ind1];
      YY[ss] = 0.0;
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int NPLearning::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                              double *settings, int *N, double **X2, 
                              double **Y2)
{
   int    ii, ss, jj, index, totPts;
   double *XT, *XX, *YY, *HX;
 
   learn(X, Y);

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[2*totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         XT[index*nInputs_+ind1] = HX[0] * ii + lowerBounds_[ind1];
         XT[index*nInputs_+ind2] = HX[1] * jj + lowerBounds_[ind2];
         XX[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         XX[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int NPLearning::gen3DGridData(double *X, double *Y, int ind1, int ind2, 
                              int ind3, double *settings, int *N, 
                              double **X2, double **Y2)
{
   int    ii, ss, jj, ll, index, totPts;
   double *XT, *XX, *YY, *HX;

   learn(X, Y);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[3*totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++)
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
            XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
            XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
            XX[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            XX[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            XX[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int NPLearning::gen4DGridData(double *X, double *Y, int ind1, int ind2, 
                              int ind3, int ind4, double *settings, 
                              int *N, double **X2, double **Y2)
{
   int    ii, ss, jj, ll, mm, index, totPts;
   double *XT, *XX, *YY, *HX;

   learn(X, Y);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[4*totPts];
   XX = (*X2);
   (*Y2) = new double[totPts];
   YY = (*Y2);
   (*N) = totPts;

   XT = new double[totPts*nInputs_];
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            for (mm = 0; mm < nPtsPerDim_; mm++)
            {
               index = ii*nPtsPerDim_*nPtsPerDim_ * nPtsPerDim_ +
                       jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
               XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
               XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
               XT[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
               XX[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               XX[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               XX[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               XX[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double NPLearning::evaluatePoint(double *X)
{
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double NPLearning::evaluatePoint(int npts, double *X, double *Y)
{
   int ss;
   for (ss = 0; ss < npts; ss++) Y[ss] = 0.0;
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double NPLearning::evaluatePointFuzzy(double *X, double &std)
{
   std = 0.0;
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double NPLearning::evaluatePointFuzzy(int npts, double *X, double *Y,
                                      double *Ystd)
{
   for (int ss = 0; ss < npts; ss++) Y[ss] = Ystd[ss] = 0.0;
   return 0.0;
}

// ************************************************************************
// learn
// ------------------------------------------------------------------------
int NPLearning::learn(double *X, double *Y)
{
   int      count, *flags, ii, jj, iInd, pInd, index, loopFlag=1;
   int      nTrain, nTest, input1, input2;
   double   *XTrain, *YTrain, *XT, *YT, mse;
   double   *XTest, *YTest, *XBS, *YBS, *Bcoefs;

printf("Not implemented yet.\n");
exit(1);
   nTrain = nSamples_ * 3 / 4;
   nTest  = nSamples_ - nTrain;
   if (nTest == 0)
   {
      printf("NPLearning ERROR: sample set too small for this method.\n");
      exit(1);
   }
   XTrain = new double[nTrain*nInputs_];
   XTest  = new double[nTest*nInputs_];
   YTrain = new double[nTrain];
   YTest  = new double[nTest];
   flags  = new int[nSamples_];
   for (ii = 0; ii < nSamples_; ii++) flags[ii] = 0;
   count  = 0;
   for (ii = 0; ii < nTrain; ii++)
   {
      index = PSUADE_rand() % nSamples_;
      for (jj = 0; jj < nInputs_; jj++)
         XTrain[count*nInputs_+jj] = X[index*nInputs_+jj];
      YTrain[count] = Y[index];
      flags[index] = 1;
      count++;
   }
   count = 0;
   for (ii = 0; ii < nSamples_; ii++)
   {
      if (flags[ii] == 0)
      {
         for (jj = 0; jj < nInputs_; jj++)
            XTest[count*nInputs_+jj] = X[ii*nInputs_+jj];
         YTest[count] = Y[ii];
         count++;
      }
   }
   delete [] flags;

   iInd = 1;
   pInd = 1;
   flags = new int[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) flags[ii] = 0;
   XBS = new double[nTrain*2];
   YBS = new double[nTrain];
   XT  = new double[nSamples_*2];
   YT  = new double[nSamples_];
   while (loopFlag)
   {
      for (ii = 0; ii < depth1_; ii++)
      {
         count = 0;
         for (jj = 0; jj < nInputs_; jj++) count += (1 - flags[jj]);
         if (count == 0)
         {
            loopFlag = 0;
            break;
         }
         if (count >= 2)
         {
            input1 = -1;
            while (input1 < 0)
            {
               input1 = PSUADE_rand() % nInputs_;
               if (flags[input1] == 1) input1 = -1;
            } 
            flags[input1] = 1;
            input2 = -1;
            while (input2 < 0)
            {
               input2 = PSUADE_rand() % nInputs_;
               if (flags[input2] == 1) input2 = -1;
            } 
            flags[input2] = 1;

            for (jj = 0; jj < nTrain; jj++)
            {
               index = PSUADE_rand() % nTrain;
               XBS[jj*2]   = XTrain[index*nInputs_+input1];
               XBS[jj*2+1] = XTrain[index*nInputs_+input2];
               YBS[jj] = YTrain[index];
            }

            RegressionAnalysis(nTrain, 2, XBS, YBS, &Bcoefs);

            for (jj = 0; jj < nTest; jj++)
            {
               XT[jj*2]   = XTest[jj*nInputs_+input1];
               XT[jj*2+1] = XTest[jj*nInputs_+input2];
            }
            evaluate(nTest, 2, Bcoefs, XT, YT);
            mse = meanSquaredError(nTest, YT, YTest);

            for (jj = 0; jj < nTest; jj++) YTest[jj] -= YT[jj];
            for (jj = 0; jj < nTrain; jj++)
            {
               XT[jj*2]   = XTrain[jj*nInputs_+input1];
               XT[jj*2+1] = XTrain[jj*nInputs_+input2];
            }
            evaluate(nTrain, 2, Bcoefs, XT, YT);
            for (jj = 0; jj < nTrain; jj++) YTrain[jj] -= YT[jj];
         }
         else
         {
            for (jj = 0; jj < nInputs_; jj++)
            {
               if (flags[jj] == 0)
               {
                  input1 = jj;
                  flags[jj] = 1;
               }
            }

            for (jj = 0; jj < nTrain; jj++)
            {
               index = PSUADE_rand() % nTrain;
               XBS[jj] = XTrain[index*nInputs_+input1];
               YBS[jj] = YTrain[index];
            }

            RegressionAnalysis(nTrain, 1, XBS, YBS, &Bcoefs);

            for (jj = 0; jj < nTest; jj++)
               XT[jj] = XTest[jj*nInputs_+input1];
            evaluate(nTest, 1, Bcoefs, XT, YT);
            mse = meanSquaredError(nTest, YT, YTest);

            for (jj = 0; jj < nTest; jj++) YTest[jj] -= YT[jj];
            for (jj = 0; jj < nTrain; jj++)
               XT[jj] = XTrain[jj*nInputs_+input1];
            evaluate(nTrain, 1, Bcoefs, XT, YT);
            for (jj = 0; jj < nTrain; jj++) YTrain[jj] -= YT[jj];
         }
         delete [] Bcoefs;
         printf("mse = %e\n", mse);
      }
   }

   delete [] XTrain;
   delete [] YTrain;
   delete [] XTest;
   delete [] YTest;
   delete [] XT;
   delete [] YT;
   delete [] XBS;
   delete [] YBS;
   delete [] flags;
   return 0;
}

// ************************************************************************
// perform regression analysis
// ------------------------------------------------------------------------
int NPLearning::RegressionAnalysis(int nSamples, int nInputs, double *Xin, 
                                  double *Y, double **BB)
{
   int    M, N, ii, mm, nn, info, wlen, last, pOrder=3;
   double *B, *XX, *txArray, *AA, *SS, *UU, *VV, *WW;
   char   jobu  = 'A', jobvt = 'A';

   if (nSamples < 10)
   {
      printf("ERROR: sample data not suitable for order = 3.\n");
      exit(1);
   }

   txArray = new double[nSamples];
   for (ii = 0; ii < nInputs; ii++)
   {
      for (mm = 0; mm < nSamples; mm++)
         txArray[mm] = Xin[mm*nInputs+ii];
      sortDbleList(nSamples, txArray);
      last = 1;
      for (mm = 1; mm < nSamples; mm++)
      {
         if (txArray[mm] != txArray[last-1])
         {
            txArray[last] = txArray[mm];
            last++;
         }
      }
      if (pOrder >= last)
      {
         printf("ERROR: sample data not suitable for order = 3.\n");
         exit(1);
      }
   }
   delete [] txArray;

   loadXMatrix(nSamples, nInputs, pOrder, Xin, &XX); 
   M = nSamples_;
   N = 10;

   wlen = 5 * M;
   AA = new double[M*N];
   UU = new double[M*M];
   SS = new double[N];
   VV = new double[M*N];
   WW = new double[wlen];
   for (mm = 0; mm < M; mm++) 
      for (nn = 0; nn < N; nn++) AA[mm+nn*M] = XX[mm+nn*M];
   dgesvd_(&jobu, &jobvt, &M, &N, AA, &M, SS, UU, &M, VV, &N, WW,
           &wlen, &info);
   if (info != 0)
   {
      printf("* Regression ERROR: dgesvd returns a nonzero (%d).\n",info);
      printf("* Regression terminates further processing.\n");
      exit(1);
   }

   for (mm = 0; mm < N; mm++) 
   {
      WW[mm] = 0.0;
      for (nn = 0; nn < M; nn++) WW[mm] += UU[mm*M+nn] * Y[nn]; 
   }
   for (nn = 0; nn < N; nn++) 
   {
      if (SS[nn] > 1.0e-15) WW[nn] /= SS[nn];
      else                  WW[nn] = 0.0;
   }
   B = new double[N];
   for (mm = 0; mm < N; mm++) 
   {
      B[mm] = 0.0;
      for (nn = 0; nn < N; nn++) B[mm] += VV[mm*N+nn] * WW[nn]; 
   }

   (*BB) = B;
   delete [] SS;
   delete [] UU;
   delete [] VV;
   delete [] AA;
   delete [] XX;
   delete [] WW;
   return 0;
}

// *************************************************************************
// load the X matrix
// -------------------------------------------------------------------------
int NPLearning::loadXMatrix(int nSamples, int nInputs, int pOrder,
                            double *X, double **XXOut)
{
   int    M, N=0, mm, nn, nn2, nn3, ind;
   double *XX=NULL;

   (*XXOut) = NULL;
   M = nSamples;
   N = 1;
   N += nInputs;
   N += nInputs * (nInputs + 1) / 2;
   for (nn = 0; nn < nInputs; nn++)
      for (nn2 = nn; nn2 < nInputs; nn2++)
         for (nn3 = nn2; nn3 < nInputs; nn3++) N++;
   XX = new double[M*N];
   for (mm = 0; mm < M; mm++) XX[mm] = 1.0;
   for (mm = 0; mm < M; mm++)
   {
      XX[mm] = 1.0;
      for (nn = 0; nn < nInputs; nn++)
         XX[M*(nn+1)+mm] = X[mm*nInputs+nn];
   }
   ind = nInputs + 1;
   for (nn = 0; nn < nInputs; nn++)
   {
      for (nn2 = nn; nn2 < nInputs; nn2++)
      {
         for (mm = 0; mm < M; mm++)
            XX[M*ind+mm] = X[mm*nInputs+nn] * X[mm*nInputs+nn2];
         ind++;
      }
   }
   for (nn = 0; nn < nInputs; nn++)
   {
      for (nn2 = nn; nn2 < nInputs; nn2++)
      {
         for (nn3 = nn2; nn3 < nInputs; nn3++)
         {
            for (mm = 0; mm < M; mm++)
               XX[M*ind+mm] = X[mm*nInputs+nn] * X[mm*nInputs+nn2] * 
                                X[mm*nInputs+nn3];
            ind++;
         }
      }
   }
   (*XXOut) = XX;
   return N;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
int NPLearning::evaluate(int npts, int nInputs, double *B, double *X, 
                            double *Y)
{
   int    ii, mm, nn, pp, offset;
   double Xdata, Xdata2, Xdata3;

   for (ii = 0; ii < npts; ii++)
   {
      Y[ii] = B[0];
      for (mm = 0; mm < nInputs; mm++)
      {
         Xdata = X[ii*nInputs+mm];
         Y[ii] += B[mm+1] * Xdata;
      }
      offset = nInputs_ + 1;
      for (mm = 0; mm < nInputs_; mm++)
      {
         Xdata = X[ii*nInputs+mm];
         for (nn = mm; nn < nInputs_; nn++)
         {
            Xdata2 = X[ii*nInputs+nn];
            Y[ii] += (B[offset++] * Xdata * Xdata2);
         }
      }
      for (mm = 0; mm < nInputs_; mm++)
      {
         Xdata = X[ii*nInputs+mm];
         for (nn = mm; nn < nInputs_; nn++)
         {
            Xdata2 = X[ii*nInputs+nn];
            for (pp = nn; pp < nInputs_; pp++)
            {
               Xdata3 = X[ii*nInputs+pp];
               Y[ii] += B[offset++] * Xdata * Xdata2 * Xdata3;
            }
         }
      }
   }
   return 0;
}

// ************************************************************************
// get the information about splitting 
// ------------------------------------------------------------------------
double NPLearning::meanSquaredError(int npts, double *Y1, double *Y2)
{
   int    ii;
   double error=0;
   for (ii = 0; ii < npts; ii++) error += pow(Y1[ii] - Y2[ii], 2.0);
   error /= (double) npts;
   return error;
}


