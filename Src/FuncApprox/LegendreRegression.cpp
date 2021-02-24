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
// Functions for the class LegendreRegression
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "LegendreRegression.h"
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PsuadeConfig.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
LegendreRegression::LegendreRegression(int nInputs,int nSamples):
                                       FuncApprox(nInputs,nSamples)
{
  int  ii;
  char line[101], *cString, winput1[500], winput2[500];

  faID_ = PSUADE_RS_REGRL;
  pOrder_ = -1;
  numPerms_ = 0;
  pcePerms_ = NULL;
  normalizeFlag_ = 0;

  printAsterisks(PL_INFO, 0);
  printf("*                Legendre Regression Analysis\n");
  printf("* R-square gives a measure of the goodness of the model.\n");
  printf("* R-square should be close to 1 if it is a good model.\n");
  printf("* Turn on rs_expert mode to output regression matrix.\n");
  printf("* Set print level to 5 to output regression error splot.\n");
  printDashes(PL_INFO, 0);
  printf("* Turn on rs_expert mode to scale the inputs to [-1, 1].\n");
  printf("* With this, statistics such as mean, variances, and\n");
  printf("* conditional variances are readily available.\n");
  printf("* Otherwise, default is: scale the inputs to [0,1].\n");
  printEquals(PL_INFO, 0);

  if (psRSExpertMode_ != 1 && psConfig_ != NULL)
  {
    cString = psConfig_->getParameter("normalize_inputs");
    if (cString != NULL) normalizeFlag_ = 1;
    cString = psConfig_->getParameter("Legendre_order");
    if (cString != NULL)
    {
      sscanf(cString, "%s %s %d", winput1, winput2, &ii);
      if (ii < 0)
      {
        printf("Legendre INFO: polynomial order %d not valid.\n",ii);
        printf("               polynomial order unchanged at %d.\n",
               pOrder_);
      }
      else
      {
        pOrder_ = ii;
        printf("Legendre INFO: polynomial order set to %d (config)\n",
               pOrder_);
      }
    }
  }
  if (psRSExpertMode_ == 1 && psInteractive_ == 1)
  {
    printf("Normalize the input parameters to [-1, 1]? (y - yes) ");
    fgets(line, 100, stdin);
    if (line[0] == 'y') normalizeFlag_ = 1;
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
LegendreRegression::~LegendreRegression()
{
  int ii;
  if (pcePerms_ != NULL)
  {
    for (ii = 0; ii < numPerms_; ii++) 
      if (pcePerms_[ii] != NULL) delete [] pcePerms_[ii];
    delete [] pcePerms_;
    pcePerms_ = NULL;
  }
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int LegendreRegression::initialize(double *X, double *Y)
{
  int      ii, status;
  psVector VecX2;

  if (pcePerms_ != NULL)
  {
    for (ii = 0; ii < numPerms_; ii++) 
      if (pcePerms_[ii] != NULL) delete [] pcePerms_[ii];
    delete [] pcePerms_;
    pcePerms_ = NULL;
  }

  if (normalizeFlag_ == 0)
  {
    VecX2.setLength(nInputs_*nSamples_);
    initInputScaling(X, VecX2.getDVector(), 1);
    status = analyze(VecX2.getDVector(), Y);
  }
  else status = analyze(X, Y);
  return status; 
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int LegendreRegression::genNDGridData(double *XIn, double *YIn, int *NOut,
                                      double **XOut, double **YOut)
{
  int totPts, ss, status;

  if (initialize(XIn, YIn) != 0)
  {
    printf("LegendreRegression::genNDGridData - ERROR detected.\n");
    (*NOut) = 0;
    return -1;
  }

  if ((*NOut) == -999) return 0;

  genNDGrid(NOut, XOut);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in LegendreRegression::genNDGridData");
  for (ss = 0; ss < totPts; ss++)
    (*YOut)[ss] = evaluatePoint(&((*XOut)[ss*nInputs_]));

  return 0;
}

// ************************************************************************
// Generate 1D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int LegendreRegression::gen1DGridData(double *XIn, double *YIn, int ind1,
                                      double *settings, int *NOut, 
                                      double **XOut, double **YOut)
{
  int      totPts, mm, nn;
  double   HX;
  psVector VecXT;

  if (initialize(XIn, YIn) != 0)
  {
    printf("LegendreRegression::gen1DGridData - ERROR detected.\n");
    (*NOut) = 0;
    return -1;
  }

  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "(*YOut) in LegendreRegression::gen1DGridData");
  VecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) VecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    VecXT[ind1] = HX * mm + lowerBounds_[ind1];
    (*XOut)[mm] = VecXT[ind1];
    (*YOut)[mm] = evaluatePoint(VecXT.getDVector());
  }
  return 0;
}

// ************************************************************************
// Generate 2D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int LegendreRegression::gen2DGridData(double *XIn, double *YIn, int ind1,
                                      int ind2,double *settings,int *NOut, 
                                      double **XOut, double **YOut)
{
  int      totPts, mm, nn, index;
  psVector VecHX, VecXT;

  if (initialize(XIn, YIn) != 0)
  {
    printf("LegendreRegression::gen2DGridData - ERROR detected.\n");
    (*NOut) = 0;
    return -1;
  }

  totPts = nPtsPerDim_ * nPtsPerDim_;
  VecHX.setLength(2);
  VecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  VecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 2];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "(*YOut) in LegendreRegression::gen2DGridData");
  VecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) VecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      index = mm * nPtsPerDim_ + nn;
      VecXT[ind1] = VecHX[0] * mm + lowerBounds_[ind1];
      VecXT[ind2] = VecHX[1] * nn + lowerBounds_[ind2];
      (*XOut)[index*2]   = VecXT[ind1];
      (*XOut)[index*2+1] = VecXT[ind2];
      (*YOut)[index] = evaluatePoint(VecXT.getDVector());
    }
  }
  return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int LegendreRegression::gen3DGridData(double *XIn, double *YIn, int ind1,
                                     int ind2, int ind3, double *settings, 
                                     int *NOut,double **XOut,double **YOut)
{
  int      totPts, mm, nn, pp, index;
  psVector VecHX, VecXT;

  if (initialize(XIn, YIn) != 0)
  {
    printf("LegendreRegression::gen3DGridData - ERROR detected.\n");
    (*NOut) = 0;
    return -1;
  }

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  VecHX.setLength(3);
  VecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  VecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  VecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 3];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "(*YOut) in LegendreRegression::gen3DGridData");
  VecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) VecXT[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        index = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
        VecXT[ind1] = VecHX[0] * mm + lowerBounds_[ind1];
        VecXT[ind2] = VecHX[1] * nn + lowerBounds_[ind2];
        VecXT[ind3] = VecHX[2] * pp + lowerBounds_[ind3];
        (*XOut)[index*3]   = VecXT[ind1];
        (*XOut)[index*3+1] = VecXT[ind2];
        (*XOut)[index*3+2] = VecXT[ind3];
        (*YOut)[index] = evaluatePoint(VecXT.getDVector());
      }
    }
  }
  return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int LegendreRegression::gen4DGridData(double *XIn, double *YIn, int ind1, 
                                      int ind2, int ind3, int ind4, 
                                      double *settings, int *NOut, 
                                      double **XOut, double **YOut)
{
  int      totPts, mm, nn, pp, qq, index;
  psVector VecHX, VecXT;

  if (initialize(XIn, YIn) != 0)
  {
    printf("LegendreRegression::gen4DGridData - ERROR detected.\n");
    (*NOut) = 0;
    return -1;
  }

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  VecHX.setLength(4);
  VecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  VecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  VecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  VecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  (*NOut) = totPts;
  (*XOut) = new double[totPts * 4];
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "(*YOut) in LegendreRegression::gen4DGridData");
  VecXT.setLength(nInputs_);
  for (nn = 0; nn < nInputs_; nn++) VecXT[nn] = settings[nn]; 
    
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
          VecXT[ind1] = VecHX[0] * mm + lowerBounds_[ind1];
          VecXT[ind2] = VecHX[1] * nn + lowerBounds_[ind2];
          VecXT[ind3] = VecHX[2] * pp + lowerBounds_[ind3];
          VecXT[ind4] = VecHX[3] * qq + lowerBounds_[ind4];
          (*XOut)[index*4]   = VecXT[ind1];
          (*XOut)[index*4+1] = VecXT[ind2];
          (*XOut)[index*4+2] = VecXT[ind3];
          (*XOut)[index*4+3] = VecXT[ind4];
          (*YOut)[index] = evaluatePoint(VecXT.getDVector());
        }
      }
    }
  }
  return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double LegendreRegression::evaluatePoint(double *X)
{
  int    ii, nn;
  double Y, multiplier, **LTable, normalX;

  if (regCoeffs_.length() <= 0)
  {
     printf("LegendreRegression ERROR: initialize has not been called.\n");
     exit(1);
  }

  LTable = new double*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) LTable[ii] = new double[pOrder_+1];
  checkAllocate(LTable[nInputs_-1],"LTable in LegendreRegr::evaluatePt");

  Y = 0.0;
  for (nn = 0; nn < numPerms_; nn++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (normalizeFlag_ == 0)
      {
        normalX = (X[ii] - XMeans_[ii]) / XStds_[ii];
        EvalLegendrePolynomials(normalX, LTable[ii]);
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
    Y += regCoeffs_[nn] * multiplier;
  }
  Y = Y * YStd_ + YMean_;

  for (ii = 0; ii < nInputs_; ii++) delete [] LTable[ii];
  delete [] LTable;
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double LegendreRegression::evaluatePoint(int npts, double *X, double *Y)
{
  int kk;
  for (kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double LegendreRegression::evaluatePointFuzzy(double *X, double &std)
{
  int      iOne=1, ii, nn;
  double   Y, stdev, dtmp, multiplier, **LTable, normalX, *Xs;
  psVector VecXs;

  if (regCoeffs_.length() <= 0)
  {
     printf("LegendreRegression ERROR: initialize has not been called.\n");
     exit(1);
  }

  LTable = new double*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) LTable[ii] = new double[pOrder_+1];
  checkAllocate(LTable[nInputs_-1],"LTable in LegendreRegr::evaluatePtFuzzy");

  VecXs.setLength(numPerms_);
  Y = 0.0;
  for (nn = 0; nn < numPerms_; nn++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (normalizeFlag_ == 0)
      {
        normalX = (X[ii] - XMeans_[ii]) / XStds_[ii];
        EvalLegendrePolynomials(normalX, LTable[ii]);
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
    Y += regCoeffs_[nn] * multiplier;
    VecXs[nn] = multiplier;
  }
  Y = Y * YStd_ + YMean_;

  stdev = 0.0;
  for (ii = 0; ii < numPerms_; ii++)
  {
    dtmp = 0.0;
    for (nn = 0; nn < numPerms_; nn++)
      dtmp += invCovMat_.getEntry(ii,nn) * VecXs[nn];
    stdev += dtmp * VecXs[ii];
  }
  std = sqrt(stdev) * YStd_;

  for (ii = 0; ii < nInputs_; ii++) delete [] LTable[ii];
  delete [] LTable;
  return Y;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double LegendreRegression::evaluatePointFuzzy(int npts,double *X,double *Y,
                                              double *Ystd)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePointFuzzy(&X[kk*nInputs_], Ystd[kk]);
  return 0.0;
}

// ************************************************************************
// perform regression analysis
// ------------------------------------------------------------------------
int LegendreRegression::analyze(double *Xin, double *Y)
{
  psVector VecX, VecY;
  VecX.load(nSamples_*nInputs_, Xin);
  VecY.load(nSamples_, Y);
  return analyze(VecX, VecY);
}

// ************************************************************************
// perform regression analysis
// ------------------------------------------------------------------------
int LegendreRegression::analyze(psVector VecX, psVector VecY)
{
  int    N, M, ii, mm, nn, info, NRevised;
  double *X, *Y, SSresid, SStotal, R2, var;
  double esum, ymax, *arrayA, *arrayXX, *SS, *UU, *VV;
  char   pString[100], response[1000];
  FILE   *fp;
  psMatrix eigMatT, MatXX;
  psVector eigVals;

  if (nSamples_ <= nInputs_)
  {
    printf("LegendreRegression::analyze ERROR - sample size too small.\n");
    return -1;
  }
  X = VecX.getDVector();
  Y = VecY.getDVector();

  GenPermutations();
   
  N = loadXMatrix(VecX, MatXX);
  M = nSamples_;

  psMatrix MatA;
  psVector VecA;
  VecA.setLength(M*N);
  arrayA = VecA.getDVector();
  arrayXX = MatXX.getMatrix1D();
  for (mm = 0; mm < M; mm++)
    for (nn = 0; nn < N; nn++)
      arrayA[mm+nn*M] = sqrt(weights_[mm]) * arrayXX[mm+nn*M];
  MatA.load(M, N, arrayA);

  if (psMasterMode_ == 1)
  {
    printf("You have the option to store the regression matrix (that\n");
    printf("is, the matrix A in Ax=b) in a matlab file for inspection.\n");
    sprintf(pString, "Store regression matrix? (y or n) ");
    getString(pString, response);
    if (response[0] == 'y')
    {
      fp = fopen("legendre_regression_matrix.m", "w");
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
          fprintf(fp, "%16.6e ", arrayA[mm+nn*M]);
        fprintf(fp, "%16.6e \n",Y[mm]);
      }
      fprintf(fp, "];\n");
      fprintf(fp, "A = AA(:,1:%d);\n", N);
      fprintf(fp, "Y = AA(:,%d);\n", N+1);
      fprintf(fp, "B = A\\Y\n");
      fclose(fp);
      printf("Regression matrix written to legendre_regression_matrix.m\n");
    }
  }

  psMatrix Umat, Vmat;
  psVector Svec;
  if (outputLevel_ > 3) printf("Running SVD ...\n");
  info = MatA.computeSVD(Umat, Svec, Vmat);
  if (outputLevel_ > 3) 
    printf("SVD completed: status = %d (should be 0).\n",info);

  if (info != 0)
  {
    printf("* LegendreRegression Info: dgesvd returns a nonzero (%d).\n",
           info);
    printf("* LegendreRegression terminates further processing.\n");
    return -1;
  }

  SS = Svec.getDVector();
  mm = 0;
  for (nn = 0; nn < N; nn++) if (SS[nn] < 0) mm++;
  if (mm > 0)
  {
    printf("* LegendreRegression WARNING: some of the singular values\n"); 
    printf("*            are < 0. May spell trouble but will\n");
    printf("*            proceed anyway (%d).\n",mm);
  }
  if (SS[0] == 0.0) NRevised = 0;
  else
  {
    NRevised = N;
    for (nn = 1; nn < N; nn++)
      if (SS[nn-1] > 0 && SS[nn]/SS[nn-1] < 1.0e-8) NRevised--;
  }
  if (NRevised < N)
  {
    printf("* LegendreRegression ERROR: \n");
    printf("*         true rank of sample matrix = %d (need %d)\n",
           NRevised, N);
    printf("*         Try lower order polynomials.\n");
    return -1;
  }
  if (psMasterMode_ == 1)
  {
    printf("* LegendreRegression: matrix singular values \n");
    printf("* The VERY small ones may cause poor numerical accuracy,\n");
    printf("* but not keeping them may ruin the approximation power.\n");
    printf("* So, select them judiciously.\n");
    for (nn = 0; nn < N; nn++)
      printf("* Singular value %5d = %e\n", nn+1, SS[nn]);
    sprintf(pString, "How many to keep (1 - %d, 0 - all) ? ", N);
    NRevised = getInt(0,N,pString);
    if (NRevised == 0) NRevised = N;
    for (nn = NRevised; nn < N; nn++) SS[nn] = 0.0;
  }
  else
  {
    NRevised = N;
    for (nn = 1; nn < N; nn++)
    {
      if (SS[nn-1] > 0 && SS[nn]/SS[nn-1] < 1.0e-8)
      {
        SS[nn] = 0.0;
        NRevised--;
      }
    }
    if (NRevised < N) 
      printf("LegendreRegression: %d singular values taken out.\n",
             N-NRevised);
  }

  psVector VecW, VecB;
  VecW.setLength(M+N);
  UU = Umat.getMatrix1D();
  for (mm = 0; mm < NRevised; mm++)
  {
    VecW[mm] = 0.0;
    for (nn = 0; nn < M; nn++)
      VecW[mm] += UU[nn+mm*M] * sqrt(weights_[nn]) * Y[nn];
  }
  for (nn = 0; nn < NRevised; nn++) VecW[nn] /= SS[nn];
  for (nn = NRevised; nn < N; nn++) VecW[nn] = 0.0;
  VecB.setLength(N);
  VV = Vmat.getMatrix1D();
  for (mm = 0; mm < N; mm++)
  {
    VecB[mm] = 0.0;
    for (nn = 0; nn < NRevised; nn++) VecB[mm] += VV[nn+mm*N] * VecW[nn];
  }

  eigMatT.load(N, N, VV);
  eigVals.load(N, SS);
  for (nn = 0; nn < N; nn++) eigVals[nn] = pow(eigVals[nn], 2.0);

  if (psMasterMode_ == 1)
  {
    fp = fopen("regression_error_file.m", "w");
    if(fp == NULL)
    {
      printf("fopen returned NULL in file %s line %d, exiting\n", 
             __FILE__, __LINE__);
      exit(1);
    }
    fprintf(fp, "%% This file contains errors of each data point.\n");
    printf("Regression_error_file.m file is to be created to show\n");
    printf("interpolation errors of each data point.\n");
  }
  else fp = NULL;

  esum = ymax = 0.0;
  for (mm = 0; mm < nSamples_; mm++)
  {
    VecW[mm] = 0.0;
    for (nn = 0; nn < N; nn++)
      VecW[mm] = VecW[mm] + arrayXX[mm+nn*nSamples_] * VecB[nn];
    VecW[mm] = VecW[mm] - Y[mm];
    esum = esum + VecW[mm] * VecW[mm] * weights_[mm];
    if (fp != NULL) 
      fprintf(fp, "%6d %24.16e\n",mm+1,VecW[mm]*sqrt(weights_[mm]));
    if (PABS(Y[mm]) > ymax) ymax = PABS(Y[mm]);
  }
  esum /= (double) nSamples_;
  esum = sqrt(esum);
  printf("* LegendreR:: interpolation rms error = %10.3e (max=%10.3e)\n",
         esum, ymax); 

  if (fp != NULL)
  {
    fclose(fp);
    printf("Now regression_error_file.m file contains data errors.\n");
  }

  computeSS(MatXX, VecY, VecB, SSresid, SStotal);
  if (SStotal == 0) R2 = 1.0;
  else              R2 = 1.0 - SSresid / SStotal;
  if (nSamples_ > N) var = SSresid / (double) (nSamples_ - N);
  else               var = 0.0;
  if (var < 0)
  { 
    if (PABS(var) > 1.0e-12)
    {
      printf("LegendreRegression WARNING: variance < 0.\n");
      printf("    Temporarily absolutize var (may have problems).\n");
      var = PABS(var);
    }
    else var = 0.0;
  }
  regCoeffs_.load(VecB);

  computeCoeffVariance(eigMatT, eigVals, var);
  psVector VecBstd;
  VecBstd.setLength(N);
  for (ii = 0; ii < N; ii++)
    VecBstd[ii] = sqrt(invCovMat_.getEntry(ii,ii));

  if (outputLevel_ >= 0)
  {
    printRC(VecB, VecBstd, MatXX, VecY);
    printf("* Regression R-square = %10.3e\n", R2);
    if (M-N-1 > 0)
      printf("* adjusted   R-square = %10.3e\n",
             1.0 - (1.0 - R2) * ((M - 1) / (M - N - 1)));
    if (outputLevel_ > 1) printSRC(VecX, VecB, SStotal);
  }
  printAsterisks(PL_INFO, 0);
 
  fp = NULL;
  if (psRSCodeGen_ == 1) fp = fopen("psuade_rs.info", "w");
  if (fp != NULL)
  {
    fprintf(fp,"/* ***********************************************/\n");
    fprintf(fp,"/* Legendre regression interpolator from PSUADE. */\n");
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
    fprintf(fp,"/* Legendre regression interpolation function    */\n");
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
    fprintf(fp,"  double std, *x2, dtmp;\n");
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
    fprintf(fp,"# Legendre Regression interpolation\n");
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
int LegendreRegression::loadXMatrix(psVector VecX, psMatrix &MatXX)
{
  int    M, N, ss, ii, nn;
  double multiplier, **LTable, normalX;
  psVector VecXX;

  if (normalizeFlag_ == 1)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      multiplier = upperBounds_[ii] - lowerBounds_[ii];
      if (multiplier == 0.0)
      {
        normalizeFlag_ = 0;
        printf("Legendre INFO: inputs not normalized - bounds not set.\n");
        break;
      }
    }
  }
  M = nSamples_;
  N = numPerms_;
  VecXX.setLength(M*N);
  LTable = new double*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) LTable[ii] = new double[pOrder_+1];
  checkAllocate(LTable[nInputs_-1],"LTable in LegendreRegr::loadXMatrix");
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
    }
    for (nn = 0; nn < numPerms_; nn++)
    {
      multiplier = 1.0;
      for (ii = 0; ii < nInputs_; ii++)
        multiplier *= LTable[ii][pcePerms_[nn][ii]];
      VecXX[nSamples_*nn+ss] = multiplier;
    }
  }
  MatXX.setFormat(PS_MAT1D);
  MatXX.load(M, N, VecXX.getDVector());
  for (ii = 0; ii < nInputs_; ii++) delete [] LTable[ii];
  delete [] LTable;
  return N;
}

// *************************************************************************
// compute SS (sum of squares) statistics
// -------------------------------------------------------------------------
int LegendreRegression::computeSS(psMatrix MatXX, psVector VecY,
                            psVector VecB, double &SSresid, double &SStotal)
{
  int    N, nn, mm;
  double *B, *Y, *arrayXX, ymean, SSresidCheck, SSreg, ddata, rdata;
                                                                                
  N = VecB.length();
  B = VecB.getDVector();
  Y = VecY.getDVector();
  arrayXX = MatXX.getMatrix1D();

  SSresid = SSresidCheck = SStotal = ymean = SSreg = 0.0;
  for (mm = 0; mm < nSamples_; mm++)
    ymean += (sqrt(weights_[mm]) * Y[mm]);
  ymean /= (double) nSamples_;
  for (mm = 0; mm < nSamples_; mm++)
  {
    ddata = 0.0;
    for (nn = 0; nn < N; nn++) ddata += (arrayXX[mm+nn*nSamples_] * B[nn]);
    rdata = Y[mm] - ddata;
    SSresid += rdata * Y[mm] * weights_[mm];
    SSresidCheck += rdata * rdata * weights_[mm];
    SSreg += (ddata - ymean) * (ddata - ymean);
  }
  for (mm = 0; mm < nSamples_; mm++)
    SStotal += weights_[mm] * (Y[mm] - ymean) * (Y[mm] - ymean);
  if (outputLevel_ > 0)
  {
    printf("* LegendreRegression: SStot  = %24.16e\n", SStotal);
    printf("* LegendreRegression: SSreg  = %24.16e\n", SSreg);
    printf("* LegendreRegression: SSres  = %24.16e\n", SSresid);
    printf("* LegendreRegression: SSres  = %24.16e (true)\n", SSresidCheck);
  }
  SSresid = SSresidCheck;
  if (outputLevel_ > 0 && nSamples_ != N)
  {
    printf("* LegendreRegression: eps(Y) = %24.16e\n", 
           SSresidCheck/(nSamples_-N));
  }
  return 0;
}

// *************************************************************************
// compute coefficient variance ((diagonal of sigma^2 (X' X)^(-1))
// -------------------------------------------------------------------------
int LegendreRegression::computeCoeffVariance(psMatrix &eigMatT, 
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
int LegendreRegression::printRC(psVector VecB, psVector VecBstd, 
                                psMatrix MatXX, psVector VecY)
{
  int    ii, jj, kk, maxTerms, flag, N;
  double coef, ddata, variance, *arrayXX;
  FILE   *fp;

  maxTerms = 0;
  for (ii = 0; ii < numPerms_; ii++) 
    if (pcePerms_[ii][0] > maxTerms) maxTerms = pcePerms_[ii][0];

  printEquals(PL_INFO, 0);
  if (normalizeFlag_ == 1)
    printf("* Note: the coefficients below have been normalized\n");
  printDashes(PL_INFO, 0);
  printf("*      ");
  for (ii = 0; ii < nInputs_; ii++) printf("     ");
  printf("               coefficient  std. error  t-value\n");

  N = VecB.length();
  arrayXX = MatXX.getMatrix1D();
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
    fp = fopen("matlablegendre.m", "w");
    fwriteHold(fp,0);
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
    printf("Legendre VCE ranking is now in matlablegendre.m.\n");
  }
  printEquals(PL_INFO, 0);
  return 0;
}

// *************************************************************************
// print standardized regression coefficients
// -------------------------------------------------------------------------
int LegendreRegression::printSRC(psVector VecX,psVector VecB,double SStotal)
{
  int      nn, mm, ii;
  double   denom, xmean, coef, Bmax, coef1;
  psVector VecB2;

  printEquals(PL_INFO, 0);
  printf("* Standardized Regression Coefficients (SRC)\n");
  printf("* When R-square is acceptable (order assumption holds), the\n");
  printf("* absolute values of SRCs provide variable importance.\n"); 
  printDashes(PL_INFO, 0);
  printf("* based on nSamples = %d\n", nSamples_);

  VecB2.setLength(nSamples_);
  denom = sqrt(SStotal / (double) (nSamples_ - 1));
  Bmax  = 0.0;
  for (nn = 0; nn < numPerms_; nn++)
  {
    coef = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
    {
      xmean = 0.0;
      for (mm = 0; mm < nSamples_; mm++) xmean += VecX[mm*nInputs_+ii];
      xmean /= (double) nSamples_;
      coef1 = 0.0;
      for (mm = 0; mm < nSamples_; mm++)
        coef1 += (VecX[mm*nInputs_+ii]-xmean)*(VecX[mm*nInputs_+ii]-xmean);
      coef1 = sqrt(coef1 / (double) (nSamples_ - 1));
      coef *= coef1;
    }
    VecB2[nn] = VecB[nn] * coef / denom;
    if (PABS(VecB2[nn]) > Bmax) Bmax = PABS(VecB2[nn]);
  }
  for (nn = 0; nn < numPerms_; nn++)
  {
    if (PABS(VecB2[nn]) > 1.0e-12 * Bmax)
    {
      printf("* Input orders: ");
      for (ii = 0; ii < nInputs_; ii++)
        printf(" %2d ",pcePerms_[nn][ii]);
      printf("= %12.4e \n", VecB2[nn]);
    }
  }
  printAsterisks(PL_INFO, 0);
  return 0;
}

// *************************************************************************
// generate all combinations of a multivariate Legendre expansion
// This code is a direct translation from Burkardt's matlab code)
// -------------------------------------------------------------------------
int LegendreRegression::GenPermutations()
{
  int  ii, kk, orderTmp, rvTmp, setFlag=0;
  char pString[500], *cString, winput1[500], winput2[500];

  if (pOrder_ < 0)
  {
    numPerms_ = 0;
    pOrder_ = 0;
    while (numPerms_ < nSamples_)
    { 
      pOrder_++;
      numPerms_ = 1;
      if (nInputs_ > pOrder_)
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
    if (numPerms_ > nSamples_) pOrder_--;
    printf("* Legendre polynomial maximum order = %d\n", pOrder_);
    if (psConfig_ != NULL)
    {
      cString = psConfig_->getParameter("Legendre_order");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d",winput1,winput2,&orderTmp);
        if (orderTmp >= 0 && orderTmp <= pOrder_)
        {
          printf("LegendreRegression: order from config file = %d\n",
                 orderTmp);
          pOrder_ = orderTmp;
          setFlag = 1;
        }
        else
        {
          printf("LegendreRegression ERROR: order from config file ");
          printf("is not valid (%d).\n",orderTmp);
        }
      }
    }
    if (setFlag == 0)
    {
      sprintf(pString, "Desired order (>=1 and <= %d) ? ", pOrder_);
      pOrder_ = getInt(1, pOrder_, pString);
    }
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
  printf("* LegendreRegression: order of polynomials   = %d\n", pOrder_);
  printf("* LegendreRegression: number of permutations = %d\n",numPerms_);
  
  pcePerms_ = new int*[numPerms_];
  for (ii = 0; ii < numPerms_; ii++) pcePerms_[ii] = new int[nInputs_];
  checkAllocate(pcePerms_[numPerms_-1],
                "pcePerms in LegendreRegression::genPermutations");

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
// Purpose: evaluate 1D Legendre polynomials (normalized)
// -------------------------------------------------------------------------
int LegendreRegression::EvalLegendrePolynomials(double X, double *LTable)
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

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double LegendreRegression::setParams(int targc, char **targv)
{
  pOrder_ = *(int *) targv[0];
  if (pOrder_ <= 0)
  {
    pOrder_ = -1;
    printf("LegendreRegression setParams: pOrder not valid.\n");
  }
  else printf("LegendreRegression setParams: pOrder set to %d.\n", pOrder_);
  return 0.0;
}

