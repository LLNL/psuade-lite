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
// Functions for the class GradLegendreRegression
// AUTHOR : CHARLES TONG
// DATE   : 2013
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GradLegendreRegression.h"
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PsuadeData.h"
#include "pData.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
GradLegendreRegression::GradLegendreRegression(int nInputs,int nSamples):
                                               FuncApprox(nInputs,nSamples)
{
  faID_ = PSUADE_RS_REGRGL;
  pOrder_ = -1;
  numPerms_ = 0;
  pcePerms_ = NULL;
  derivSample_ = NULL;
  normalizeFlag_ = 1;
  nOutputs_ = 1;
  printAsterisks(PL_INFO, 0);
  printf("* GradLegendreRegression constructor\n");
  printDashes(PL_INFO, 0);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GradLegendreRegression::~GradLegendreRegression()
{
  int ii;
  if (pcePerms_ != NULL)
  {
    for (ii = 0; ii < numPerms_; ii++) 
      if (pcePerms_[ii] != NULL) delete [] pcePerms_[ii];
    delete [] pcePerms_;
  }
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int GradLegendreRegression::initialize(double *X, double *Y)
{
  int  ii, status;
  char pString[101];
  FILE *fp;

  if (pcePerms_ != NULL)
  {
    for (ii = 0; ii < numPerms_; ii++) 
      if (pcePerms_[ii] != NULL) delete [] pcePerms_[ii];
    delete [] pcePerms_;
    pcePerms_ = NULL;
  }

  if (derivSample_ == NULL)
  {
    printf("Users need to provide a sample file (PSUADE format) that\n");
    printf("has derivative information (where the sample outputs are\n");
    printf("derivatives of Y with respect to each input).\n");
    printf("Note: The number of inputs/outputs in this file should match\n");
    printf("      that of the sample used for creating response surface.\n");
    sprintf(pString,"Enter name of derivative file (PSUADE format): ");
    getString(pString, gradFile_);
    ii = strlen(gradFile_);
    gradFile_[ii-1] = '\0';
    fp = fopen(gradFile_, "r");
    if (fp == NULL)
    {
      printf("GradLegendreRegr ERROR: cannot open derivative file (%s)\n",
             gradFile_);
      exit(1);
    }
  }

  status = analyze(X, Y);
  if (status != 0)
  {
    printf("GradLegendreRegression::initialize - ERROR detected.\n");
    return -1;
  }
  return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int GradLegendreRegression::genNDGridData(double *X, double *Y, int *NOut,
                                          double **XOut, double **YOut)
{
  int totPts, ss, status;
  status = analyze(X, Y);
  if (status != 0)
  {
    printf("GradLegendreRegression::genNDGridData - ERROR detected.\n");
    (*NOut) = 0;
    return -1;
  }

  if ((*NOut) == -999) return 0;

  genNDGrid(NOut, XOut);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in GradLegendre::genNDGridData");

  for (ss = 0; ss < totPts; ss++)
    (*YOut)[ss] = evaluatePoint(&((*XOut)[ss*nInputs_]));
  return 0;
}

// ************************************************************************
// Generate 1D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int GradLegendreRegression::gen1DGridData(double *X, double *Y, int ind1,
                                          double *settings, int *NOut, 
                                          double **XOut, double **YOut)
{
  int    totPts, mm, nn;
  double HX;
  psVector vecXT;

  (*NOut) = -999;
  genNDGridData(X, Y, NOut, NULL, NULL);

  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts];
  (*YOut) = new double[totPts];
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
// Generate 2D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int GradLegendreRegression::gen2DGridData(double *XIn,double *YIn,int ind1,
                                     int ind2, double *settings, int *NOut, 
                                     double **XOut, double **YOut)
{
  int totPts, mm, nn, index;
  psVector vecHX, vecXT;

  (*NOut) = -999;
  genNDGridData(XIn, YIn, NOut, NULL, NULL);

  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 2];
  (*YOut) = new double[totPts];
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      index = mm * nPtsPerDim_ + nn;
      vecXT[ind1] = vecHX[0] * mm + lowerBounds_[ind1];
      vecXT[ind2] = vecHX[1] * nn + lowerBounds_[ind2];
      (*XOut)[index*2]   = vecXT[ind1];
      (*XOut)[index*2+1] = vecXT[ind2];
      (*YOut)[index] = evaluatePoint(vecXT.getDVector());
    }
  }

  return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int GradLegendreRegression::gen3DGridData(double *XIn,double *YIn,int ind1,
                                   int ind2, int ind3, double *settings, 
                                   int *NOut, double **XOut, double **YOut)
{
  int totPts, mm, nn, pp, index;
  psVector vecHX, vecXT;

  (*NOut) = -999;
  genNDGridData(XIn, YIn, NOut, NULL, NULL);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 3];
  (*YOut) = new double[totPts];
  vecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) vecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        index = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
        vecXT[ind1] = vecHX[0] * mm + lowerBounds_[ind1];
        vecXT[ind2] = vecHX[1] * nn + lowerBounds_[ind2];
        vecXT[ind3] = vecHX[2] * pp + lowerBounds_[ind3];
        (*XOut)[index*3]   = vecXT[ind1];
        (*XOut)[index*3+1] = vecXT[ind2];
        (*XOut)[index*3+2] = vecXT[ind3];
        (*YOut)[index] = evaluatePoint(vecXT.getDVector());
      }
    }
  }

  return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int GradLegendreRegression::gen4DGridData(double *XIn,double *YIn,int ind1, 
                                          int ind2, int ind3, int ind4, 
                                          double *settings, int *NOut, 
                                          double **XOut, double **YOut)
{
  int totPts, mm, nn, pp, qq, index;
  psVector vecHX, vecXT;

  (*NOut) = -999;
  genNDGridData(XIn, YIn, NOut, NULL, NULL);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  vecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 4];
  (*YOut) = new double[totPts];
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
          index = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                  nn*nPtsPerDim_*nPtsPerDim_ + pp * nPtsPerDim_ + qq;
          vecXT[ind1] = vecHX[0] * mm + lowerBounds_[ind1];
          vecXT[ind2] = vecHX[1] * nn + lowerBounds_[ind2];
          vecXT[ind3] = vecHX[2] * pp + lowerBounds_[ind3];
          vecXT[ind4] = vecHX[3] * qq + lowerBounds_[ind4];
          (*XOut)[index*4]   = vecXT[ind1];
          (*XOut)[index*4+1] = vecXT[ind2];
          (*XOut)[index*4+2] = vecXT[ind3];
          (*XOut)[index*4+3] = vecXT[ind4];
          (*YOut)[index] = evaluatePoint(vecXT.getDVector());
        }
      }
    }
  }

  return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double GradLegendreRegression::evaluatePoint(double *X)
{
  int    ii, nn;
  double Y, multiplier, **LTable, normalX;

  if (regCoeffs_.length() <= 0)
  {
    printf("GradLegendre ERROR: initialize has not been called.\n");
    exit(1);
  }

  LTable = new double*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) LTable[ii] = new double[pOrder_+1];
  checkAllocate(LTable[nInputs_-1], "LTable in GradLegendre::evaluatePoint");

  Y = 0.0;
  for (nn = 0; nn < numPerms_; nn++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (normalizeFlag_ == 0)
      {
        EvalLegendrePolynomials(X[ii], LTable[ii]);
      }
      else
      {
        normalX = X[ii] - lowerBounds_[ii];
        normalX /= (upperBounds_[ii] - lowerBounds_[ii]);
        normalX = normalX * 2.0 - 1.0;
        EvalLegendrePolynomials(normalX, LTable[ii]);
      }
    }
    multiplier = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
      multiplier *= LTable[ii][pcePerms_[nn][ii]];
    Y += regCoeffs_[nn]* multiplier;
  }

  for (ii = 0; ii < nInputs_; ii++) delete [] LTable[ii];
  delete [] LTable;
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double GradLegendreRegression::evaluatePoint(int npts, double *X, double *Y)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double GradLegendreRegression::evaluatePointFuzzy(double *X, double &std)
{
  int      ii, nn;
  double   Y, multiplier, **LTable, normalX, stdev, dtmp;
  psVector VecXs;

  if (regCoeffs_.length() <= 0)
  {
    printf("GradLegendre ERROR: initialize has not been called.\n");
    exit(1);
  }

  LTable = new double*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) LTable[ii] = new double[pOrder_+1];
  checkAllocate(LTable[nInputs_-1], "LTable in GradLegendre::evaluatePoint");

  VecXs.setLength(numPerms_);
  Y = 0.0;
  for (nn = 0; nn < numPerms_; nn++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (normalizeFlag_ == 0)
      {
        EvalLegendrePolynomials(X[ii], LTable[ii]);
      }
      else
      {
        normalX = X[ii] - lowerBounds_[ii];
        normalX /= (upperBounds_[ii] - lowerBounds_[ii]);
        normalX = normalX * 2.0 - 1.0;
        EvalLegendrePolynomials(normalX, LTable[ii]);
      }
    }
    multiplier = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
      multiplier *= LTable[ii][pcePerms_[nn][ii]];
    Y += regCoeffs_[nn]* multiplier;
    VecXs[nn] = multiplier;
  }

  stdev = 0.0;
  for (ii = 0; ii < numPerms_; ii++)
  {
    dtmp = 0.0;
    for (nn = 0; nn < numPerms_; nn++)
      dtmp += invCovMat_.getEntry(ii,nn) * VecXs[nn];
    stdev += dtmp * VecXs[ii];
  }
  std = sqrt(stdev);

  for (ii = 0; ii < nInputs_; ii++) delete [] LTable[ii];
  delete [] LTable;
  return Y;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double GradLegendreRegression::evaluatePointFuzzy(int npts, double *X,
                                                  double *Y, double *Ystd)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePointFuzzy(&X[kk*nInputs_], Ystd[kk]);
  return 0.0;
}

// ************************************************************************
// perform regression analysis
// ------------------------------------------------------------------------
int GradLegendreRegression::analyze(double *Xin, double *Y)
{
  psVector VecX, VecY;
  VecX.load(nSamples_*nInputs_, Xin);
  VecY.load(nSamples_, Y);
  return analyze(VecX, VecY);
}

// ************************************************************************
// perform regression analysis
// ------------------------------------------------------------------------
int GradLegendreRegression::analyze(psVector VecX, psVector VecY)
{
  int    N, M, ii, mm, nn, wlen, info, NRevised, status, nInps;
  double SSreg, SStotal, R2, var, esum, ymax, *UU, *VV, *arrayXX;
  char   pString[100];
  FILE   *fp;
  pData    pPtr;
  psMatrix eigMatT, MatXX, MatA, MatU, MatV;
  psVector eigVals, VecYY, VecS, VecW, VecB;
  PsuadeData *ioPtr;

  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("GradLegendreRegression::analyze ERROR - invalid arguments.\n");
    exit(1);
  } 
  GenPermutations();
   
  if (outputLevel_ >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printf("*                GradLegendre Regression Analysis\n");
    printf("* R-square gives a measure of the goodness of the model.\n");
    printf("* R-square should be close to 1 if it is a good model.\n");
    printf("* Turn on rs_expert mode to output regression matrix.\n");
    printDashes(PL_INFO, 0);
    printf("* Each input will be scaled to [-1, 1].\n");
    printf("* With this, statistics such as mean, variances, and\n");
    printf("* conditional variances are readily available.\n");
    printEquals(PL_INFO, 0);
  }

  if (derivSample_ == NULL)
  {
    status = ioPtr->readPsuadeFile(gradFile_);
    if (status != 0)
    {
      printf("GradLegendreReg ERROR: problem reading derivative file.\n");
      exit(1);
    }
    ioPtr->getParameter("input_ninputs", pPtr);
    int nInps = pPtr.intData_;
    if (nInps != nInputs_)
    {
      printf("GradLegendreReg ERROR: nInput mismatch in derivative file.\n");
      exit(1);
    }
    ioPtr->getParameter("output_noutputs", pPtr);
    if (pPtr.intData_ != nInputs_)
    {
      printf("GradLegendreReg ERROR: invalid nOutputs in derivative file.\n");
      printf("                       Should be equal to %d\n", nInputs_);
      exit(1);
    }
    ioPtr->getParameter("method_nsamples", pPtr);
    if (pPtr.intData_ != nSamples_)
    {
      printf("GradLegendreReg ERROR: invalid sample size in derivative file\n");
      printf("                        Should be equal to %d\n", nSamples_);
      exit(1);
    }
    ioPtr->getParameter("output_noutputs", pPtr);
    nn = pPtr.intData_;
    if (nn != nInputs_)
    {
      printf("GradLegendreReg ERROR - derivative file should have %d outputs\n",
             nInputs_);
      printf("                               You only have %d outputs.\n",nn);
      exit(1);
    }
    ioPtr->getParameter("output_sample", pPtr);
    delete ioPtr;
  }

  double *derivatives;
  if (derivSample_ == NULL) derivatives = pPtr.dbleArray_;
  else                      derivatives = derivSample_;
  M = nSamples_ * (nInputs_ + 1);
  VecYY.setLength(M);
  for (mm = 0; mm < nSamples_; mm++)
  {
    VecYY[mm*(nInputs_+1)] = VecY[mm];
    for (nn = 0; nn < nInputs_; nn++)
      VecYY[mm*(nInputs_+1)+nn+1] = derivatives[mm*nInputs_+nn];
  }

  N = loadXMatrix(VecX, MatXX);
  psVector VecA;
  VecA.setLength(M*N);
  arrayXX = MatXX.getMatrix1D();
  for (mm = 0; mm < M; mm++)
    for (nn = 0; nn < N; nn++) VecA[mm+nn*M] = arrayXX[mm+nn*M];
  MatA.load(M, N, VecA.getDVector());

  if (psRSExpertMode_ == 1 && outputLevel_ > 0)
  {
    fp = fopen("Grad_Legendre_regression_matrix.m", "w");
    if(fp == NULL)
    {
      printf("fopen returned NULL in file %s line %d, exiting\n",
             __FILE__, __LINE__);
      exit(1);
    }
    fprintf(fp, "%% the sample matrix where svd is computed\n");
    fprintf(fp, "%% the last column is the right hand side\n");
    fprintf(fp, "%% B is the vector of coefficients\n");
    fprintf(fp, "AA = [\n");
    for (mm = 0; mm < M; mm++)
    {
      for (nn = 0; nn < N; nn++) 
        fprintf(fp, "%16.6e ", VecA[mm+nn*M]);
      fprintf(fp, "%16.6e \n",VecYY[mm]);
    }
    fprintf(fp, "];\n");
    fprintf(fp, "A = AA(:,1:%d);\n", N);
    fprintf(fp, "B = AA(:,%d);\n", N+1);
    fclose(fp);
    printf("Regression matrix written to Grad_Legendre_regression_matrix.m\n");
  }

  if (outputLevel_ > 3) printf("Running SVD ...\n");
  info = MatA.computeSVD(MatU, VecS, MatV);
  if (outputLevel_ > 3)
    printf("SVD completed: status = %d (should be 0).\n",info);

  if (info != 0)
  {
    printf("* GradLegendreRegression Info: dgesvd returns a nonzero (%d).\n",
           info);
    printf("* GradLegendreRegression terminates further processing.\n");
    printf("* Consult PSUADE developers for advice.\n");
    return -1;
  }

  if (VecS[0] == 0.0) NRevised = 0;
  else
  {
    NRevised = N;
    for (nn = 1; nn < N; nn++)
      if (VecS[nn-1] > 0 && VecS[nn]/VecS[nn-1] < 1.0e-8) NRevised--;
  }
  if (NRevised < N)
  {
    printf("* GradLegendreRegression ERROR: true rank of sample = %d\n",
           NRevised);
    printf("*                               need %d\n", N);
    printf("*                               Try lower order.\n");
    return -1;
  }
  if (psRSExpertMode_ == 1)
  {
    printf("GradLegendreReg: singular values for the Vandermonde matrix\n");
    printf("The VERY small ones may cause poor numerical accuracy,\n");
    printf("but not keeping them may ruin the approximation power.\n");
    printf("So, select them judiciously.\n");
    for (nn = 0; nn < N; nn++)
      printf("Singular value %5d = %e\n", nn+1, VecS[nn]);
    sprintf(pString, "How many to keep (1 - %d) ? ", N);
    NRevised = getInt(1,N,pString);
    for (nn = NRevised; nn < N; nn++) VecS[nn] = 0.0;
  }
  else
  {
    NRevised = N;
    for (nn = 1; nn < N; nn++)
    {
      if (VecS[nn-1] > 0 && VecS[nn]/VecS[nn-1] < 1.0e-8)
      {
        VecS[nn] = 0.0;
        NRevised--;
      }
    }
  }

  VecW.setLength(M+N);
  UU = MatU.getMatrix1D();
  for (mm = 0; mm < NRevised; mm++)
  {
    VecW[mm] = 0.0;
    for (nn = 0; nn < M; nn++) VecW[mm] += UU[mm*M+nn] * VecYY[nn];
  }
  for (nn = 0; nn < NRevised; nn++) VecW[nn] /= VecS[nn];
  for (nn = NRevised; nn < N; nn++) VecW[nn] = 0.0;
  VecB.setLength(N);
  VV = MatV.getMatrix1D();
  for (mm = 0; mm < N; mm++)
  {
    VecB[mm] = 0.0;
    for (nn = 0; nn < NRevised; nn++) VecB[mm] += VV[nn+mm*N] * VecW[nn];
  }

  eigMatT.load(N, N, VV);
  eigVals.load(N, VecS.getDVector());
  for (nn = 0; nn < N; nn++) eigVals[nn] = pow(eigVals[nn], 2.0);

  fp = NULL;
  if (psMasterMode_ == 1)
  {
    fp = fopen("grad_legendre_regression_err_file.m", "w");
    if(fp == NULL)
    {
      printf("fopen returned NULL in file %s line %d, exiting\n",
             __FILE__, __LINE__);
      exit(1);
    }
    fprintf(fp, "%% This file contains errors of each data point.\n");
  }

  esum = ymax = 0.0;
  for (mm = 0; mm < M; mm++)
  {
    VecW[mm] = 0.0;
    for (nn = 0; nn < N; nn++)
      VecW[mm] = VecW[mm] + arrayXX[mm+nn*M] * VecB[nn];
    VecW[mm] = VecW[mm] - VecYY[mm];
    esum = esum + VecW[mm] * VecW[mm];
    if (fp != NULL) 
      fprintf(fp, "%6d %24.16e\n",mm+1,VecW[mm]);
    if (PABS(VecYY[mm]) > ymax) ymax = PABS(VecYY[mm]);
  }
  esum /= (double) M;
  esum = sqrt(esum);
  printf("* GradLegendreRegression: LS mean error = %10.3e (max=%10.3e)\n",
         esum, ymax); 

  if (fp != NULL)
  {
    fclose(fp);
    printf("FILE grad_legendre_regression_err_file contains data errors\n");
  }

  status = computeSS(MatXX, VecYY, VecB, SSreg, SStotal);
  
  if (status == 0)
  {
    if (SStotal == 0) R2 = 1.0;
    else              R2  = SSreg / SStotal;
    if (M > N) var = (SStotal - SSreg) / (double) (M - N);
    else       var = 0.0;
    if (var < 0)
    {
      if (PABS(var) > 1.0e-12)
      {
        printf("GradLegendreRegression WARNING: variance %e < 0.\n",var);
        printf("                       Temporary absolutize var.\n");
        var = PABS(var);
      }
      else var = 0;
    }
  }
  else
  {
    if (SStotal == 0) R2 = 1.0;
    else              R2  = 1.0 - SSreg / SStotal;
    if (M > N) var = SSreg / (double) (M - N);
    if (var < 0)
    {
      if (PABS(var) > 1.0e-12)
           printf("GradLegendreRegression WARNING: var %e < 0.\n", var);
      else var = 0;
    }
  }
  regCoeffs_.load(VecB);

  computeCoeffVariance(eigMatT, eigVals, var);
  psVector VecBstd;
  VecBstd.setLength(N);
  for (ii = 0; ii < N; ii++)
    VecBstd[ii] = sqrt(invCovMat_.getEntry(ii,ii));

  if (outputLevel_ >= 0)
  {
    printRC(VecB, VecBstd, MatXX, VecYY);
    printf("* Regression R-square = %10.3e\n", R2);
    if (M-N-1 > 0)
      printf("* adjusted   R-square = %10.3e\n",
             1.0 - (1.0 - R2) * ((M - 1) / (M - N - 1)));
  }
  printAsterisks(PL_INFO, 0);

  fp = NULL;
  if (psRSCodeGen_ == 1) fp = fopen("psuade_rs.info", "w");
  if (fp != NULL)
  {
    fprintf(fp,"/* ***********************************************/\n");
    fprintf(fp,"/* GradLegendre regression interpolator from PSUADE. */\n");
    fprintf(fp,"/* ==============================================*/\n");
    fprintf(fp,"/* This file contains information for interpolation\n");
    fprintf(fp,"   using response surface. Follow the steps below:\n");
    fprintf(fp,"   1. move this file to *.c file (e.g. main.c)\n");
    fprintf(fp,"   2. Compile main.c (cc -o main main.c -lm) \n");
    fprintf(fp,"   3. run: main input output\n");
    fprintf(fp,"          where input has the number of inputs and\n");
    fprintf(fp,"          the input values\n");
    fprintf(fp,"*/\n");
    fprintf(fp,"/* ==============================================*/\n");
    fprintf(fp,"#include <math.h>\n");
    fprintf(fp,"#include <stdlib.h>\n");
    fprintf(fp,"#include <stdio.h>\n");
    fprintf(fp,"int interpolate(int,double *,double *,double *);\n");
    fprintf(fp,"main(int argc, char **argv) {\n");
    fprintf(fp,"  int    i, iOne=1, nInps;\n");
    fprintf(fp,"  double X[%d], Y, S;\n",nInputs_);
    fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
    fprintf(fp,"  if (argc != 3) {\n");
    fprintf(fp,"     printf(\"ERROR: not enough argument.\\n\");\n");
    fprintf(fp,"     exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  fIn = fopen(argv[1], \"r\");\n");
    fprintf(fp,"  if (fIn == NULL) {\n");
    fprintf(fp,"     printf(\"ERROR: cannot open input file.\\n\");\n");
    fprintf(fp,"     exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  fscanf(fIn, \"%%d\", &nInps);\n");
    fprintf(fp,"  if (nInps != %d) {\n", nInputs_);
    fprintf(fp,"    printf(\"ERROR - wrong nInputs.\\n\");\n");
    fprintf(fp,"    exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  for (i=0; i<%d; i++) fscanf(fIn, \"%%lg\", &X[i]);\n",
            nInputs_);
    fprintf(fp,"  fclose(fIn);\n");
    fprintf(fp,"  interpolate(iOne, X, &Y, &S);\n");
    fprintf(fp,"  printf(\"Y = %%e\\n\", Y);\n");
    fprintf(fp,"  printf(\"S = %%e\\n\", S);\n");
    fprintf(fp,"  fOut = fopen(argv[2], \"w\");\n");
    fprintf(fp,"  if (fOut == NULL) {\n");
    fprintf(fp,"     printf(\"ERROR: cannot open output file.\\n\");\n");
    fprintf(fp,"     exit(1);\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  fprintf(fOut,\" %%e\\n\", Y);\n");
    fprintf(fp,"  fclose(fOut);\n");
    fprintf(fp,"}\n");
    fprintf(fp,"/* ==============================================*/\n");
    fprintf(fp,"/* GradLegendre regression interpolation function*/\n");
    fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
    fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
    fprintf(fp," * ... */\n");
    fprintf(fp,"/* ==============================================*/\n");
    fprintf(fp,"static int\n");
    fprintf(fp,"pcePerms[%d][%d] = \n", numPerms_, nInputs_);
    fprintf(fp,"{\n");
    for (mm = 0; mm < numPerms_; mm++)
    {
      fprintf(fp,"  {");
      for (ii = 0; ii < nInputs_-1; ii++)
        fprintf(fp," %d,", pcePerms_[mm][ii]);
      fprintf(fp," %d },\n", pcePerms_[mm][nInputs_-1]);
    }
    fprintf(fp,"};\n");
    fprintf(fp,"static double\n");
    fprintf(fp,"invCovMat[%d][%d] = \n", numPerms_, numPerms_);
    fprintf(fp,"{\n");
    for (mm = 0; mm < numPerms_; mm++)
    {
       fprintf(fp,"  {");
       for (ii = 0; ii < numPerms_-1; ii++)
          fprintf(fp," %24.16e,", invCovMat_.getEntry(mm,ii));
       fprintf(fp," %24.16e },\n", invCovMat_.getEntry(mm,numPerms_-1));
    }
    fprintf(fp,"};\n");
    fprintf(fp,"static double\n");
    fprintf(fp,"regCoefs[%d] = \n", numPerms_);
    fprintf(fp,"{\n");
    for (mm = 0; mm < numPerms_; mm++)
      fprintf(fp," %24.16e,", regCoeffs_[mm]);
    fprintf(fp,"};\n");
    fprintf(fp,"/* ==============================================*/\n");
    fprintf(fp,"int interpolate(int npts,double *X,double *Y,double *S){\n");
    fprintf(fp,"  int    ii, kk, ss, nn;\n");
    fprintf(fp,"  double *x, y, **LTable, normX, mult;\n");
    fprintf(fp,"  double mean, std, *x2, dtmp;\n");
    fprintf(fp,"  LTable = (double **) malloc(%d * sizeof(double*));\n",
               nInputs_);
    fprintf(fp,"  for (ii = 0; ii < %d; ii++)\n", nInputs_);
    fprintf(fp,"    LTable[ii] = (double *) malloc((%d+1)*sizeof(double));\n",
            pOrder_);
    fprintf(fp,"  x2 = (double *) malloc(%d * sizeof(double));\n",numPerms_);
    fprintf(fp,"  for (ss = 0; ss < npts; ss++) {\n");
    fprintf(fp,"    x = &X[ss * %d];\n", nInputs_);
    fprintf(fp,"    y = 0.0;\n");
    fprintf(fp,"    for (nn = 0; nn < %d; nn++) {\n", numPerms_);
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (normalizeFlag_ == 0)
      {
        fprintf(fp,"      normX = X[%d] - %24.16e;\n",ii, XMeans_[ii]);
        fprintf(fp,"      normX /= %24.16e;\n", XStds_[ii]);
        fprintf(fp,"      EvalLegendrePolynomials(normX,LTable[%d]);\n",ii);
      }
      else
      {
        fprintf(fp,"      normX = X[%d] - %24.16e;\n",
                ii, lowerBounds_[ii]);
        fprintf(fp,"      normX /= (%24.16e - %24.16e);\n",
                upperBounds_[ii], lowerBounds_[ii]);
        fprintf(fp,"      normX = normX * 2.0 - 1.0;\n");
        fprintf(fp,"      EvalLegendrePolynomials(normX,LTable[%d]);\n",
                ii);
      }
    }
    fprintf(fp,"      mult = 1.0;\n");
    for (ii = 0; ii < nInputs_; ii++)
    fprintf(fp,"      mult *= LTable[%d][pcePerms[nn][%d]];\n",ii,ii);
    fprintf(fp,"      y += regCoefs[nn] * mult;\n");
    fprintf(fp,"      x2[nn] = mult;\n");
    fprintf(fp,"    }\n");
    fprintf(fp,"    Y[ss] = y * %e + %e;\n", YStd_, YMean_);
    fprintf(fp,"    std = 0.0;\n");
    fprintf(fp,"    for (ii = 0; ii < %d; ii++) {\n",numPerms_);
    fprintf(fp,"      dtmp = 0.0;\n");
    fprintf(fp,"      for (kk = 0; kk < %d; kk++)\n",numPerms_);
    fprintf(fp,"        dtmp += invCovMat[ii][kk] * x2[kk];\n");
    fprintf(fp,"      std += dtmp * x2[ii];\n");
    fprintf(fp,"    }\n");
    fprintf(fp,"    std = sqrt(std);\n");
    fprintf(fp,"    S[ss] = std;\n");
    fprintf(fp,"  }\n");
    for (ii = 0; ii < nInputs_; ii++)
       fprintf(fp,"  free(LTable[%d]);\n", ii);
    fprintf(fp,"  free(LTable);\n");
    fprintf(fp,"  free(x2);\n");
    fprintf(fp,"  return 0;\n");
    fprintf(fp,"}\n");
    fprintf(fp,"int EvalLegendrePolynomials(double X, double *LTable) {\n");
    fprintf(fp,"  int    ii;\n");
    fprintf(fp,"  LTable[0] = 1.0;\n");
    fprintf(fp,"  if (%d >= 1) {\n", pOrder_);
    fprintf(fp,"     LTable[1] = X;\n");
    fprintf(fp,"     for (ii = 2; ii <= %d; ii++)\n", pOrder_);
    fprintf(fp,"        LTable[ii] = ((2 * ii - 1) * X * LTable[ii-1] -\n");
    fprintf(fp,"                      (ii - 1) * LTable[ii-2]) / ii;\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  return 0;\n");
    fprintf(fp,"}\n");
    fprintf(fp,"/* ==============================================*/\n");
    fclose(fp);
    printf("FILE psuade_rs.info contains information about\n");
    printf("     the Legendre polynomial.\n");
  }
  fp = NULL;
  if (psRSCodeGen_ == 1) fp = fopen("psuade_rs.py", "w");
  if (fp != NULL)
  {
    fwriteRSPythonHeader(fp);
    fprintf(fp,"#==================================================\n");
    fprintf(fp,"# GradLegendre Regression interpolation\n");
    fprintf(fp,"#==================================================\n");
    fwriteRSPythonCommon(fp);
    fprintf(fp,"pcePerms = [\n");
    for (mm = 0; mm < numPerms_; mm++)
    {
      fprintf(fp," [ %d", pcePerms_[mm][0]);
      for (ii = 1; ii < nInputs_; ii++)
        fprintf(fp,", %d", pcePerms_[mm][ii]);
      fprintf(fp," ],\n");
    }
    fprintf(fp,"]\n");
    fprintf(fp,"invCovMat = [\n");
    for (mm = 0; mm < numPerms_; mm++)
    {
      fprintf(fp," [ %24.16e", invCovMat_.getEntry(mm,0));
      for (nn = 1; nn < numPerms_; nn++)
        fprintf(fp,", %24.16e", invCovMat_.getEntry(mm,nn));
      fprintf(fp," ],\n");
    }
    fprintf(fp,"]\n");
    fprintf(fp,"regCoefs = [\n");
    for (mm = 0; mm < numPerms_-1; mm++)
      fprintf(fp," %24.16e,\n", regCoeffs_[mm]);
    fprintf(fp," %24.16e ]\n", regCoeffs_[numPerms_-1]);
    fprintf(fp,"###################################################\n");
    fprintf(fp,"def EvalLegendrePolynomials(X) :\n");
    fprintf(fp,"  LTable = %d * [0.0]\n", pOrder_+1);
    fprintf(fp,"  LTable[0] = 1.0;\n");
    fprintf(fp,"  if (%d >= 1) :\n", pOrder_);
    fprintf(fp,"    LTable[1] = X;\n");
    fprintf(fp,"    for ii in range(%d) : \n", pOrder_-1);
    fprintf(fp,"      jj = ii + 2\n");
    fprintf(fp,"      LTable[jj] = ((2 * jj - 1) * X * LTable[jj-1] -\n");
    fprintf(fp,"                    (jj - 1) * LTable[jj-2]) / jj;\n");
    fprintf(fp,"  return LTable\n");
    fprintf(fp,"###################################################\n");
    fprintf(fp,"# Regression interpolation function  \n");
    fprintf(fp,"# X[0], X[1],   .. X[m-1]   - first point\n");
    fprintf(fp,"# X[m], X[m+1], .. X[2*m-1] - second point\n");
    fprintf(fp,"# ... \n");
    fprintf(fp,"#==================================================\n");
    fprintf(fp,"def interpolate(X): \n");
    fprintf(fp,"  nSamp = int(len(X) / %d + 1.0e-8)\n",nInputs_);
    fprintf(fp,"  Xt = %d * [0.0]\n", nInputs_);
    fprintf(fp,"  Xs = %d * [0.0]\n", numPerms_);
    fprintf(fp,"  Ys = 2 * nSamp * [0.0]\n");
    fprintf(fp,"  for ss in range(nSamp): \n");
    fprintf(fp,"    for ii in range(%d) : \n", nInputs_);
    fprintf(fp,"      Xt[ii] = X[ss*%d+ii]\n",nInputs_);
    fprintf(fp,"    Y = 0.0\n");
    fprintf(fp,"    for nn in range(%d): \n", numPerms_);
    fprintf(fp,"      mult = 1.0\n");
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (normalizeFlag_ == 0)
      {
        fprintf(fp,"      x2 = Xt[%d] - %24.16e\n",ii, XMeans_[ii]);
        fprintf(fp,"      x2 = x2 / %24.16e\n", XStds_[ii]);
        fprintf(fp,"      LTable = EvalLegendrePolynomials(x2)\n");
      }
      else
      {
        fprintf(fp,"      x2 = Xt[%d] - %24.16e\n",ii,lowerBounds_[ii]);
        fprintf(fp,"      x2 = x2 / (%24.16e - %24.16e)\n",
                upperBounds_[ii], lowerBounds_[ii]);
        fprintf(fp,"      x2 = x2 * 2.0 - 1.0\n");
        fprintf(fp,"      LTable = EvalLegendrePolynomials(x2)\n");
      }
      fprintf(fp,"      mult *= LTable[pcePerms[nn][%d]]\n",ii);
    }
    fprintf(fp,"      Xs[nn] = mult\n");
    fprintf(fp,"      Y = Y + regCoefs[nn] * mult\n");
    fprintf(fp,"    Ys[2*ss] = Y * %e + %e\n",YStd_, YMean_);
    fprintf(fp,"    std = 0.0\n");
    fprintf(fp,"    for jj in range(%d): \n", numPerms_);
    fprintf(fp,"      dtmp = 0.0\n");
    fprintf(fp,"      for kk in range(%d): \n", numPerms_);
    fprintf(fp,"        dtmp = dtmp + invCovMat[jj][kk] * Xs[kk]\n");
    fprintf(fp,"      std = std + dtmp * Xs[jj]\n");
    fprintf(fp,"    Ys[ss*2+1] = math.sqrt(std)\n");
    fprintf(fp,"  return Ys\n");
    fprintf(fp,"###################################################\n");
    fprintf(fp,"# main program\n");
    fprintf(fp,"#==================================================\n");
    fprintf(fp,"infileName  = sys.argv[1]\n");
    fprintf(fp,"outfileName = sys.argv[2]\n");
    fprintf(fp,"inputs = getInputData(infileName)\n");
    fprintf(fp,"outputs = interpolate(inputs)\n");
    fprintf(fp,"genOutputFile(outfileName, outputs)\n");
    fprintf(fp,"###################################################\n");
    printf("FILE psuade_rs.py contains the final Legendre polynomial\n");
    printf("     functional form.\n");
    fclose(fp);
  }
  return 0;
}

// *************************************************************************
// load the X matrix
// -------------------------------------------------------------------------
int GradLegendreRegression::loadXMatrix(psVector VecX, psMatrix &MatXX)
{
  int    M, N=0, ss, ii, kk, nn;
  double multiplier, **LTable, **DTable, normalX, range;
  psVector VecXX;

  for (ii = 0; ii < nInputs_; ii++)
  {
    multiplier = upperBounds_[ii] - lowerBounds_[ii];
    if (multiplier == 0.0)
    {
      printf("GradLegendreRegression ERROR: inputs not normalized.\n");
      exit(1);
    }
  }

  M = nSamples_ * (nInputs_ + 1);
  N = numPerms_;
  VecXX.setLength(M*N);
  LTable = new double*[nInputs_];
  DTable = new double*[nInputs_];
  checkAllocate(DTable, "DTable in GradLegendre::loadXMatrix");
  for (ii = 0; ii < nInputs_; ii++)
  {
    LTable[ii] = new double[pOrder_+1];
    DTable[ii] = new double[pOrder_+1];
  }

  for (ss = 0; ss < nSamples_; ss++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (normalizeFlag_ == 0) normalX = VecX[ss*nInputs_+ii]; 
      else
      {
        normalX = VecX[ss*nInputs_+ii] - lowerBounds_[ii];
        normalX /= (upperBounds_[ii] - lowerBounds_[ii]);
        normalX = normalX * 2.0 - 1.0;
      }
      EvalLegendrePolynomials(normalX, LTable[ii]);
      EvalLegendrePolynomialsDerivative(normalX, DTable[ii]);
      if (normalizeFlag_ == 1)
      {
        range = 2.0 / (upperBounds_[ii] - lowerBounds_[ii]);
        for (kk = 0; kk < pOrder_+1; kk++) DTable[ii][kk] *= range;
      }
    }
    for (nn = 0; nn < numPerms_; nn++)
    {
      multiplier = 1.0;
      for (ii = 0; ii < nInputs_; ii++)
        multiplier *= LTable[ii][pcePerms_[nn][ii]];
      VecXX[M*nn+ss*(nInputs_+1)] = multiplier;
    }
    for (kk = 0; kk < nInputs_; kk++)
    {
      for (nn = 0; nn < numPerms_; nn++)
      {
        multiplier = 1.0;
        for (ii = 0; ii < nInputs_; ii++)
        {
          if (kk == ii) multiplier *= DTable[ii][pcePerms_[nn][ii]];
          else          multiplier *= LTable[ii][pcePerms_[nn][ii]];
        }
        VecXX[M*nn+ss*(nInputs_+1)+kk+1] = multiplier;
      }
    }
  }
  MatXX.setFormat(PS_MAT1D);
  MatXX.load(M, N, VecXX.getDVector());

  for (ii = 0; ii < nInputs_; ii++)
  {
    delete [] LTable[ii];
    delete [] DTable[ii];
  }
  delete [] LTable;
  delete [] DTable;
  return N;
}

// *************************************************************************
// compute SS (sum of squares) statistics
// -------------------------------------------------------------------------
int GradLegendreRegression::computeSS(psMatrix MatXX, psVector VecY,
                                      psVector VecB, double &SSreg, 
                                      double &SStotal)
{
  int    nn, mm, M, N;
  double ymean, SSresid, SSresidCheck, ddata, rdata, *arrayXX;
                                                                               
  N = VecB.length();
  arrayXX = MatXX.getMatrix1D();
  M = nSamples_ * (nInputs_ + 1);
  SSresid = SSresidCheck = SStotal = ymean = SSreg = 0.0;
  for (mm = 0; mm < M; mm++) ymean += VecY[mm];
  ymean /= (double) M;
  for (mm = 0; mm < M; mm++)
  {
    ddata = 0.0;
    for (nn = 0; nn < N; nn++) ddata += (arrayXX[mm+nn*M] * VecB[nn]);
    rdata = VecY[mm] - ddata;
    SSresid += VecY[mm] * rdata;
    SSresidCheck += rdata * rdata;
    SSreg += (ddata - ymean) * (ddata - ymean);
  }
  for (mm = 0; mm < M; mm++)
    SStotal += (VecY[mm] - ymean) * (VecY[mm] - ymean);
  printf("* GradLegendreRegression: SStot = %e\n", SStotal);
  printf("* GradLegendreRegression: SSreg = %e\n", SSreg);
  printf("* GradLegendreRegression: SSres = %e\n", SSresid);
  printf("* GradLegendreRegression: SSres = %e (true)\n", SSresidCheck);
  SSresid = SSresidCheck;
  if (SStotal < SSreg)
  {
    SSreg = SSresid;
    return 1;
  }
  return 0;
}

// *************************************************************************
// compute coefficient variance ((diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int GradLegendreRegression::computeCoeffVariance(psMatrix &eigMatT,
                               psVector &eigVals, double var)
{
  int      ii, jj, nRows;
  double   invEig, dtmp;
  psMatrix tMat;

  nRows = eigMatT.nrows();
  tMat.setDim(nRows, nRows);

  for (ii = 0; ii < nRows; ii++)
  {
    invEig = eigVals[ii];
    if (invEig != 0.0) invEig = 1.0 / invEig;
    for (jj = 0; jj < nRows; jj++)
    {
      dtmp = invEig * eigMatT.getEntry(ii,jj) * var;
      tMat.setEntry(jj, ii, dtmp);
    }
  }
  eigMatT.matmult(tMat, invCovMat_);
  return 0;
}

// *************************************************************************
// print statistics
// -------------------------------------------------------------------------
int GradLegendreRegression::printRC(psVector VecB, psVector VecBstd,
                                    psMatrix, psVector VecY)
{
  int    ii, jj, kk, maxTerms, flag;
  double coef, ddata, variance;
  FILE   *fp;

  maxTerms = 0;
  for (ii = 0; ii < numPerms_; ii++) 
    if (pcePerms_[ii][0] > maxTerms) maxTerms = pcePerms_[ii][0];

  printEquals(PL_INFO, 0);
  if (normalizeFlag_ == 1)
    printf("* Note: coefficients below are for normalized input ranges.\n");
  printDashes(PL_INFO, 0);
  printf("*      ");
  for (ii = 0; ii < nInputs_; ii++) printf("     ");
  printf("               coefficient  std. error  t-value\n");

  for (ii = 0; ii < numPerms_; ii++)
  {
    if (PABS(VecBstd[ii]) < 1.0e-15) coef = 0.0;
    else                             coef = VecB[ii] / VecBstd[ii]; 
    {
      printf("* Input orders: ");
      for (jj = 0; jj < nInputs_; jj++)
        printf(" %4d ", pcePerms_[ii][jj]);
      printf("= %11.3e %11.3e %11.3e\n", VecB[ii], VecBstd[ii], coef);
    }
  }
  flag = 1;
  for (ii = 0; ii < nInputs_; ii++)
    if (upperBounds_[ii] != 1.0) flag = 0;
  for (ii = 0; ii < nInputs_; ii++)
    if (lowerBounds_[ii] != -1.0) flag = 0;
  if (normalizeFlag_ == 1 || flag == 1)
  {
    printDashes(PL_INFO, 0);
    printf("* Mean     = %12.4e\n", VecB[0]);
    coef = 0.0;
    for (jj = 1; jj < numPerms_; jj++) 
    {
      ddata = VecB[jj];
      for (kk = 0; kk < nInputs_; kk++)
        ddata /= sqrt(1.0+pcePerms_[jj][kk]*2); 
      coef = coef + ddata * ddata;
    }
    printf("* Variance = %12.4e\n", coef);
    variance = coef;
    fp = fopen("matlabgradlegendre.m", "w");
    fwritePlotCLF(fp);
    fprintf(fp, "A = [\n");
    for (ii = 0; ii < nInputs_; ii++)
    {
      coef = 0.0;
      for (jj = 1; jj < numPerms_; jj++)
      {
        flag = 1;
        for (kk = 0; kk < nInputs_; kk++)
          if (kk != ii && pcePerms_[jj][kk] != 0) flag = 0;
        if (flag == 1)
        {
          ddata = VecB[jj];
          for (kk = 0; kk < nInputs_; kk++)
            ddata /= sqrt(1.0+pcePerms_[jj][kk]*2); 
          coef = coef + ddata * ddata;
        }
      }
      fprintf(fp, "%e\n", coef/variance);
      printf("* Conditional variance %4d = %12.4e\n", ii+1, coef);
    }
    fprintf(fp, "];\n");
    fprintf(fp, "bar(A, 0.8);\n");
    fwritePlotAxes(fp);
    fwritePlotTitle(fp, "Legendre VCE Rankings");
    fwritePlotXLabel(fp, "Input parameters");
    fwritePlotYLabel(fp, "Rank Metric");
    fclose(fp);
    printf("Legendre VCE ranking is now in matlabgradlegendre.m.\n");
  }
  printEquals(PL_INFO, 0);
  return 0;
}

// *************************************************************************
// generate all combinations of a multivariate Legendre expansion
// This code is a direct translation from Burkardt's matlab code)
// -------------------------------------------------------------------------
int GradLegendreRegression::GenPermutations()
{
  int  ii, kk, orderTmp, rvTmp, M;
  char pString[500];

  M = nSamples_ * (nInputs_ + 1);
  if (pOrder_ <= 0)
  {
    numPerms_ = 0;
    pOrder_ = 0;
    while (numPerms_ < M)
    { 
      pOrder_++;
      numPerms_ = 1;
      if (nInputs_ < pOrder_)
      {
        for (ii = nInputs_+pOrder_; ii > nInputs_; ii--)
          numPerms_ = numPerms_ * ii / (nInputs_+pOrder_-ii+1);
      }
      else
      {
        for (ii = nInputs_+pOrder_; ii > pOrder_; ii--)
          numPerms_ = numPerms_ * ii / (nInputs_+pOrder_-ii+1);
      }
    }
    if (numPerms_ > M) pOrder_--;
    printf("* Legendre polynomial maximum order = %d\n", pOrder_);
    sprintf(pString, "Desired order (>=1 and <= %d) ? ", pOrder_);
    pOrder_ = getInt(1, pOrder_, pString);
  }
  numPerms_ = 1;
  if (nInputs_ < pOrder_)
  {
    for (ii = nInputs_+pOrder_; ii > nInputs_; ii--)
      numPerms_ = numPerms_ * ii / (nInputs_+pOrder_-ii+1);
  }
  else
  {
    for (ii = nInputs_+pOrder_; ii > pOrder_; ii--)
      numPerms_ = numPerms_ * ii / (nInputs_+pOrder_-ii+1);
  }
  printf("* GradLegendreRegression: order of polynomials   = %d\n", pOrder_);
  printf("* GradLegendreRegression: number of permutations = %d\n",numPerms_);
   
  pcePerms_ = new int*[numPerms_];
  for (ii = 0; ii < numPerms_; ii++) pcePerms_[ii] = new int[nInputs_];
  checkAllocate(pcePerms_[numPerms_-1],
                "pcePerms in GradLegendre::genPermutations");

  numPerms_ = 0;
  for (kk = 0; kk <= pOrder_; kk++)
  {
    orderTmp = kk;
    rvTmp = 0;
    pcePerms_[numPerms_][0] = orderTmp;
    for (ii = 1; ii < nInputs_; ii++) pcePerms_[numPerms_][ii] = 0;
    while (pcePerms_[numPerms_][nInputs_-1] != kk)
    {
      numPerms_++;
      for (ii = 0; ii < nInputs_; ii++)
        pcePerms_[numPerms_][ii] = pcePerms_[numPerms_-1][ii];
      if (orderTmp > 1) rvTmp = 1;
      else              rvTmp++;
      pcePerms_[numPerms_][rvTmp-1] = 0;
      orderTmp = pcePerms_[numPerms_-1][rvTmp-1];
      pcePerms_[numPerms_][0] = orderTmp - 1;
      pcePerms_[numPerms_][rvTmp] = pcePerms_[numPerms_-1][rvTmp] + 1;
    }
    numPerms_++;
  }
  return 0;
}

// *************************************************************************
// Purpose: evaluate 1D Legendre polynomials (semi-normalized)
// -------------------------------------------------------------------------
int GradLegendreRegression::EvalLegendrePolynomials(double X, double *LTable)
{
  int    ii;
  LTable[0] = 1.0;
  if (pOrder_ >= 1)
  {
    LTable[1] = X;
    for (ii = 2; ii <= pOrder_; ii++)
      LTable[ii] = ((2 * ii - 1) * X * LTable[ii-1] -
                    (ii - 1) * LTable[ii-2]) / ii;
  }
  return 0;
}

// *************************************************************************
// Purpose: evaluate 1D Legendre polynomials derivative (semi-normalized)
// -------------------------------------------------------------------------
int GradLegendreRegression::EvalLegendrePolynomialsDerivative(double X, 
                                                              double *DTable)
{
  int    ii;
  DTable[0] = 0.0;
  if (pOrder_ >= 1) DTable[1] = 1.0;
  if (pOrder_ >= 2)
  {
    DTable[2] = 3.0 * X;
    for (ii = 3; ii <= pOrder_; ii++)
      DTable[ii] = ((2 * ii - 1) * X * DTable[ii-1] -
                    ii * DTable[ii-2]) / (ii - 1.0);
  }
  return 0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double GradLegendreRegression::setParams(int targc, char **targv)
{
  if (targc == 1)
  {
    pOrder_ = *(int *) targv[0];
    if (pOrder_ <= 0)
    {
      printf("LegendreRegression setParams: pOrder not valid.\n");
      exit(1);
    }
    printf("LegendreRegression setParams: pOrder set to %d.\n",pOrder_);
  }
  else if (targc == 2 && !strcmp(targv[0], "nOutputs"))
  {
    nOutputs_ = *(int *) targv[1];
    if (nOutputs_ <= 0)
    {
      printf("LegendreRegression setParams: nOutputs not valid.\n");
      exit(1);
    }
  }
  else if (targc == 2 && !strcmp(targv[0], "deriv_sample"))
  {
    derivSample_ = (double *) targv[1];
  }
  else
  {
    printf("LegendreRegression setParams ERROR: keyword not recognized.\n");
    exit(1);
  }
  return 0.0;
}

