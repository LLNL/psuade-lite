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
// Functions for the class CSRegression
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sysdef.h"
#include "Psuade.h"
#include "CSRegression.h"
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
CSRegression::CSRegression(int nInputs,int nSamples):
                    FuncApprox(nInputs,nSamples)
{
  faID_ = PSUADE_RS_CSREG;
  numTerms_  = 0;
  termOrders_ = NULL;
  regCoeffs_ = NULL;
  regStdevs_ = NULL;
  fuzzyC_    = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
CSRegression::~CSRegression()
{
  cleanUp();
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int CSRegression::initialize(double *X, double *Y)
{
  int status, ii;
 
  cleanUp();
   
  status = analyze(X, Y);
  if (status != 0)
  {
    printf("CSRegression: ERROR detected in regression analysis.\n");
    return -1;
  }
  return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int CSRegression::genNDGridData(double *XIn, double *YIn, int *NOut, 
                                double **XOut, double **YOut)
{
  int mm, totPts;

  if (initialize(XIn,YIn) != 0)
  {
    printf("CSRegression: ERROR detected in regression analysis.\n");
    (*NOut) = 0;
    return -1;
  }

  if ((*NOut) == -999) return 0;

  genNDGrid(NOut, XOut);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

   (*YOut) = new double[totPts];
   checkAllocate(*YOut, "YOut in CSRegression::genNDGrid");
   (*NOut) = totPts;
   for (mm = 0; mm < totPts; mm++)
     (*YOut)[mm] = evaluatePoint(&((*XOut)[mm*nInputs_]));
   return 0;
}

// ************************************************************************
// Generate 1D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int CSRegression::gen1DGridData(double *XIn, double *YIn, int ind1,
                                double *settings, int *NOut, 
                                double **XOut, double **YOut)
{
  int    totPts, mm, nn;
  double HX;
  psVector vecXT;

  if (initialize(XIn,YIn) != 0)
  {
    printf("CSRegression: ERROR detected in regression analysis.\n");
    (*NOut) = 0;
    return -1;
  }

  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in CSRegression::gen1DGrid");
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    vecXT[ind1] = HX * mm + lowerBounds_[ind1];
    (*XOut)[mm] = vecXT[ind1];
    (*YOut)[mm] = evaluatePoint(vecXT.getDVector());
  }
  return 0;
}

// ************************************************************************
// Generate 2D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int CSRegression::gen2DGridData(double *XIn, double *YIn, int ind1,
                                int ind2, double *settings, int *NOut, 
                                double **XOut, double **YOut)
{
  int totPts, mm, nn, ind;
  psVector vecHX, vecXT;

  if (initialize(XIn,YIn) != 0)
  {
    printf("CSRegression: ERROR detected in regression analysis.\n");
    (*NOut) = 0;
    return -1;
  }

  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 2];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in CSRegression::gen2DGrid");
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      ind = mm * nPtsPerDim_ + nn;
      vecXT[ind1] = vecHX[0] * mm + lowerBounds_[ind1];
      vecXT[ind2] = vecHX[1] * nn + lowerBounds_[ind2];
      (*XOut)[ind*2]   = vecXT[ind1];
      (*XOut)[ind*2+1] = vecXT[ind2];
      (*YOut)[ind] = evaluatePoint(vecXT.getDVector());
    }
  }
  return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int CSRegression::gen3DGridData(double *XIn, double *YIn, int ind1,
                                int ind2, int ind3, double *settings, 
                                int *NOut, double **XOut, double **YOut)
{
  int totPts, mm, nn, pp, ind;
  psVector vecHX, vecXT;

  if (initialize(XIn,YIn) != 0)
  {
    printf("CSRegression: ERROR detected in regression analysis.\n");
    (*NOut) = 0;
    return -1;
  }

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 3];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in CSRegression::gen3DGrid");
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        ind = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
        vecXT[ind1] = vecHX[0] * mm + lowerBounds_[ind1];
        vecXT[ind2] = vecHX[1] * nn + lowerBounds_[ind2];
        vecXT[ind3] = vecHX[2] * pp + lowerBounds_[ind3];
        (*XOut)[ind*3]   = vecXT[ind1];
        (*XOut)[ind*3+1] = vecXT[ind2];
        (*XOut)[ind*3+2] = vecXT[ind3];
        (*YOut)[ind] = evaluatePoint(vecXT.getDVector());
      }
    }
  }
  return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int CSRegression::gen4DGridData(double *XIn, double *YIn,int ind1,int ind2,
                                int ind3, int ind4, double *settings, 
                                int *NOut, double **XOut, double **YOut)
{
  int totPts, mm, nn, pp, qq, ind;
  psVector vecHX, vecXT;

  if (initialize(XIn,YIn) != 0)
  {
    printf("CSRegression: ERROR detected in regression analysis.\n");
    (*NOut) = 0;
    return -1;
  }
 
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  vecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 4];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in CSRegression::gen4DGrid");
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
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
          vecXT[ind1] = vecHX[0] * mm + lowerBounds_[ind1];
          vecXT[ind2] = vecHX[1] * nn + lowerBounds_[ind2];
          vecXT[ind3] = vecHX[2] * pp + lowerBounds_[ind3];
          vecXT[ind4] = vecHX[3] * qq + lowerBounds_[ind4];
          (*XOut)[ind*4]   = vecXT[ind1];
          (*XOut)[ind*4+1] = vecXT[ind2];
          (*XOut)[ind*4+2] = vecXT[ind3];
          (*XOut)[ind*4+3] = vecXT[ind4];
          (*YOut)[ind] = evaluatePoint(vecXT.getDVector());
        }
      }
    }
  }
  return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double CSRegression::evaluatePoint(double *X)
{
  int    mm, nn;
  double Xdata, Xdata2, Y;

  if (regCoeffs_ == NULL)
  {
    printf("CSRegression ERROR: need to call initialize first.\n");
    return 0.0;
  }
  Y = regCoeffs_[0];
  for (mm = 1; mm < numTerms_; mm++)
  {
    Xdata = 0.0;
    for (nn = 0; nn < nInputs_; nn++)
    {
      Xdata2 = (X[nn] - XMeans_[nn]) / XStds_[nn]; 
      Xdata += pow(Xdata2, 1.0*termOrders_[mm][nn]);
    }
    Y += regCoeffs_[mm] * Xdata;
  }
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double CSRegression::evaluatePoint(int npts, double *X, double *Y)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double CSRegression::evaluatePointFuzzy(double *X, double &std)
{
  int    mm, nn, cc, nTimes=100;
  double accum, mean, stds, Xdata, Xdata2;
  psVector vecYs;

  if (regCoeffs_ == NULL)
  {
    printf("CSRegression ERROR: initialize has not been called.\n");
    exit(1);
  }
 
  vecYs.setLength(nTimes);

  mean = 0.0;
  for (cc = 0; cc < nTimes; cc++)
  {
    accum = fuzzyC_[0][cc];
    for (mm = 1; mm < numTerms_; mm++)
    {
      Xdata = 0.0;
      for (nn = 0; nn < nInputs_; nn++)
      {
        Xdata2 = (X[nn] - XMeans_[nn]) / XStds_[nn]; 
        Xdata += pow(Xdata2, 1.0*termOrders_[mm][nn]);
      }
      accum += fuzzyC_[mm][cc] * Xdata;
    }
    vecYs[cc] = accum;
    mean += accum;
  }
  mean /= (double) nTimes;
  stds = 0.0;
  for (cc = 0; cc < nTimes; cc++)
    stds += (vecYs[cc] - mean) * (vecYs[cc] - mean);
  stds = sqrt(stds / (nTimes - 1));
  std = stds; 
  mean = mean * YStd_ + YMean_;
  std  = std * YStd_;
  return mean;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double CSRegression::evaluatePointFuzzy(int npts, double *X, double *Y,
                                        double *Ystd)
{
  if (regCoeffs_ == NULL)
  {
    printf("CSRegression ERROR: initialize has not been called.\n");
    exit(1);
  }
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePointFuzzy(&(X[kk*nInputs_]), Ystd[kk]);
  return 0.0;
}

// ************************************************************************
// perform regression analysis
// ------------------------------------------------------------------------
int CSRegression::analyze(double *Xin, double *YIn)
{
  int    M, N, ii, mm, nn, status;
  double *B, *XX, SSresid, SStotal, R2, *XTX, var, *Bstd, *X;
  double enorm, ymax, dtemp;
  psVector vecXT;

  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("CSRegression ERROR: consult PSUADE developers.\n");
    exit( 1 );
  } 
   
  if (outputLevel_ >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printf("*               Compressed Sensing Regression Analysis\n");
    printf("* R-squared gives a measure of the goodness of the model.\n");
    printf("* R-squared should be close to 1 if it is a good model.\n");
    printDashes(PL_INFO, 0);
  }

  vecXT.setLength(nSamples_*nInputs_);
  if (psMasterMode_ == 1) 
  {
    printf("* CSRegression INFO: scaling turned off.\n");
    printf("*              To turn on scaling, use rs_expert mode.\n");
    initInputScaling(Xin, vecXT.getDVector(), 0);
  }
  else initInputScaling(Xin, vecXT.getDVector(), 1);

  status = optimize(vecXT.getDVector(), YIn);
  if (status < 0)
  {
    printf("* CSRegression ERROR: error in compression.\n");
    exit(1);
  }   

  loadXMatrix(vecXT.getDVector(), &XX);
  N = numTerms_;
  M = nSamples_;

  enorm = ymax = 0.0;
  for (mm = 0; mm < M; mm++)
  {
    dtemp = 0.0;
    for (nn = 0; nn < N; nn++) 
      dtemp += XX[mm+nn*nSamples_] * regCoeffs_[nn];
    dtemp -= YIn[mm];
    enorm = enorm + dtemp * dtemp;
    if (PABS(YIn[mm]) > ymax) ymax = PABS(YIn[mm]);
  }
  enorm /= (double) nSamples_;
  enorm = sqrt(enorm);
  if (outputLevel_ > 1)
    printf("* CSRegression: interpolation RMS error = %11.4e (Ymax=%9.2e)\n",
           enorm, ymax); 

  computeSS(N, XX, YIn, regCoeffs_, SSresid, SStotal);
  R2 = 1.0;
  if (SStotal != 0.0) R2  = 1.0 - SSresid / SStotal;
  if (nSamples_ > N) var = SSresid / (double) (nSamples_ - N);
  else               var = 0.0;
  if (var < 0)
  { 
    if (PABS(var) > 1.0e-12)
         printf("CSRegression WARNING: var < 0.\n");
    else var = 0;
  }

  Bstd = new double[N];
  checkAllocate(Bstd, "Bstd in CSRegression::analyze");
  computeXTX(N, XX, &XTX);
  computeCoeffVariance(N, XTX, var, Bstd);
  regStdevs_ = Bstd;

  PDFManager *pdfman = new PDFManager();
  int    cc, nTimes=100;
  int    *inPDFs = new int[N];
  double *inMeans = new double[N];
  double *inStds = new double[N];
  double *inUppers = new double[N];
  double *inLowers = new double[N];
  checkAllocate(inLowers, "inLowers in CSRegression::analyze");
  for (nn = 0; nn < N; nn++)
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
  pdfman->initialize(N,inPDFs,inMeans,inStds,covMatrix_,NULL,NULL);
  psVector vLower, vUpper, vOut;
  vLower.load(N, inLowers);
  vUpper.load(N, inUppers);
  vOut.setLength(N*nTimes);
  pdfman->genSample(nTimes, vOut, vLower, vUpper);
  fuzzyC_ = new double*[N];
  for (nn = 0; nn < N; nn++)
  {
    fuzzyC_[nn] = new double[nTimes];
    for (cc = 0; cc < nTimes; cc++)
      fuzzyC_[nn][cc] = vOut[cc*N+nn];
  }
  checkAllocate(fuzzyC_[N-1], "FuzzyC in CSRegression::analyze");
  delete pdfman;
  delete [] inPDFs;
  delete [] inStds;
  delete [] inMeans;
  delete [] inLowers;
  delete [] inUppers;

  if (outputLevel_ >= 0)
  {
    B = regCoeffs_;
    Bstd = regStdevs_;
    printRC(N, B, Bstd, XX, YIn);
    //printCoefs(N, B);
    printf("* CSRegression R-square = %12.4e ", R2);
    printf("(SSresid,SStotal=%10.2e,%10.2e)\n", SSresid, SStotal);
    if ((M - N - 1) > 0)
      printf("* adjusted   R-square = %12.4e\n",
             1.0 - (1.0 - R2) * ((M - 1) / (M - N - 1)));
  }

  numTerms_  = N;
  delete [] XX;
  delete [] XTX;
  return 0;
}

// *************************************************************************
// load the X matrix
// -------------------------------------------------------------------------
int CSRegression::loadXMatrix(double *X, double **XXOut)
{
  int    M, N=0, mm, nn, ii;
  double *XX=NULL, dtemp, dprod;

  (*XXOut) = NULL;

  M = nSamples_;
  XX = new double[M*N];
  checkAllocate(XX, " CSRegression XX");

  for (mm = 0; mm < M; mm++) XX[mm] = 1.0;
  for (nn = 1; nn < numTerms_; nn++)
  {
    for (mm = 0; mm < M; mm++)
    {
      dprod = 1.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        dtemp = X[mm*nInputs_+ii];
        dprod *= pow(dtemp, 1.0*termOrders_[nn][ii]);
      }
      XX[M*(nn+1)+mm] = dprod;
    }
  }
  (*XXOut) = XX;
  return 0;
}

// *************************************************************************
// form X^T X 
// -------------------------------------------------------------------------
int CSRegression::computeXTX(int N, double *X, double **XXOut)
{
  int    nn, nn2, mm;
  double *XX, coef;

  XX = new double[nSamples_*N];
  checkAllocate(XX, " CSRegression XX");
  for (nn = 0; nn < N; nn++)
  {
    for (nn2 = 0; nn2 < N; nn2++)
    {
      coef = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
        coef += X[nn*nSamples_+mm] * weights_[mm] * X[nn2*nSamples_+mm];
      XX[nn*N+nn2] = coef;
    }
  }
  (*XXOut) = XX;
  return 0;
}

// *************************************************************************
// compute SS (sum of squares) statistics
// -------------------------------------------------------------------------
int CSRegression::computeSS(int N, double *XX, double *Y,
                          double *B, double &SSresid, double &SStotal)
{
  int    nn, mm;
  double rdata, ymean, SSreg, ddata, SSresidCheck;

  SSresid = SSresidCheck = SStotal = SSreg = ymean = 0.0;
  for (mm = 0; mm < nSamples_; mm++) ymean += sqrt(weights_[mm]) * Y[mm];
  ymean /= (double) nSamples_;
  for (mm = 0; mm < nSamples_; mm++)
  {
    ddata = 0.0;
    for (nn = 0; nn < N; nn++) ddata += (XX[mm+nn*nSamples_] * B[nn]);
    rdata = Y[mm] - ddata;
    SSresidCheck += rdata * rdata * weights_[mm];
    SSresid += rdata * Y[mm] * weights_[mm];
    SSreg += (ddata - ymean) * (ddata - ymean);
  }
  for (mm = 0; mm < nSamples_; mm++)
    SStotal += weights_[mm] * (Y[mm] - ymean) * (Y[mm] - ymean);
  if (outputLevel_ > 0)
  {
    printf("* CSRegression: SStot  = %24.16e\n", SStotal);
    printf("* CSRegression: SSreg  = %24.16e\n", SSreg);
    printf("* CSRegression: SSres  = %24.16e\n", SSresid);
    printf("* CSRegression: SSres  = %24.16e (true)\n", SSresidCheck);
  }
  SSresid = SSresidCheck;
  if (outputLevel_ > 0 && nSamples_ != N)
  {
    printf("* CSRegression: eps(Y) = %24.16e\n",
           SSresidCheck/(nSamples_-N));
  }
  return 0;
}

// *************************************************************************
// compute coefficient variances (diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int CSRegression::computeCoeffVariance(int N,double *XX,double var,double *B)
{
  int    nn, nn2, lwork, iOne=1, info, errCnt=0;
  double *B2, *work, *XT;
  char   trans[1];

  (*trans) = 'N';
  B2 = new double[N];
  XT = new double[N*N];
  lwork = 2 * N * N;
  work  = new double[lwork];
  checkAllocate(work, " CSRegression work");

  for (nn = 0; nn < N; nn++)
  {
    for (nn2 = 0; nn2 < N*N; nn2++) XT[nn2] = XX[nn2];
    for (nn2 = 0; nn2 < N; nn2++) B2[nn2] = 0.0;
    B2[nn] = var;
    dgels_(trans, &N, &N, &iOne, XT, &N, B2, &N, work, &lwork, &info);
    if (info != 0)
      printf("CSRegression WARNING: dgels returns error %d.\n",info);
    if (B2[nn] < 0) errCnt++;
    if (B2[nn] < 0) B[nn] = sqrt(-B2[nn]);
    else            B[nn] = sqrt(B2[nn]);
  }
  if (errCnt > 0)
  {
    printf("* CSRegression WARNING: some of the coefficient variances\n");
    printf("*            are < 0. May spell trouble but will\n");
    printf("*            proceed anyway (%d).\n", errCnt);
  }
  delete [] B2;
  delete [] XT;

  int    *ipiv = new int[N+1];
  double *invA = new double[lwork];
  checkAllocate(invA, " CSRegression invA");
  double ddata, ddata2;
  FILE   *fp;
  for (nn = 0; nn < N*N; nn++) invA[nn] = XX[nn];
  dgetrf_(&N, &N, invA, &N, ipiv, &info);
  if (info != 0)
    printf("CSRegression WARNING: dgels returns error %d.\n",info);
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
    printf("CSRegression WARNING:\n");
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
int CSRegression::printRC(int N,double *B,double *Bvar,double *XX,double *Y)
{
  int mm, ii;
  printEquals(PL_INFO, 0);
  printf("*** Note: these coefficients may not be true coefficients due\n");
  printf("***       to sample matrix scaling (i.e., they may be scaled).\n");
  printDashes(PL_INFO, 0);
  printf("* ");
  for (ii = 0; ii < nInputs_; ii++) printf("   ");
  printf("  coefficient   std. error   t-value\n");
  printDashes(PL_INFO, 0);
  for (mm = 0; ii < numTerms_; ii++)
  {
    for (ii = 0; ii < nInputs_; ii++)
      printf(" %d ", termOrders_[mm][ii]);
    printf("= %12.4e %12.4e %12.4e\n", B[mm], Bvar[mm], B[ii]/Bvar[mm]);
  }
  printDashes(PL_INFO, 0);
  return 0;
}

// ************************************************************************
// clean up
// ------------------------------------------------------------------------
void CSRegression::cleanUp()
{
  int ii;
  if (termOrders_ != NULL) 
  {
    for (ii = 1; ii < numTerms_; ii++) delete [] termOrders_[ii];
    delete [] termOrders_;
  }
  if (regCoeffs_ != NULL) delete [] regCoeffs_;
  if (regStdevs_ != NULL) delete [] regStdevs_;
  if (fuzzyC_ != NULL)
  {
    for (ii = 0; ii < numTerms_; ii++) delete [] fuzzyC_[ii];
    delete [] fuzzyC_;
  }
}

// *************************************************************************
// optimize 
// -------------------------------------------------------------------------
int CSRegression::optimize(double *X, double *Y)
{
  printf("* CSRegression ERROR: CSRegression not implemented yet.\n");
  return -1;
}

