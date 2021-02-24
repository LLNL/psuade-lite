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
// Functions for the class SparseGridRegression
// AUTHOR : CHARLES TONG
// DATE   : 2011
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SGRegression.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Globals.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
SparseGridRegression::SparseGridRegression(int nInputs,int nSamples):
                                           FuncApprox(nInputs,nSamples)
{
  int    ii, kk;
  double ddata, ddata2;
  FILE   *fp;

  faID_     = PSUADE_RS_REGSG;
  pcePerms_ = NULL;
  numPerms_ = 0;

  if (outputLevel_ >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printf("*       Sparse Grid Regression Analysis\n");
    printDashes(PL_INFO, 0);
    printf("* Note: This function looks for a ps_sparse_grid_info file\n");
    printf("*       for regression information.\n");
    printf("The file should be in the format:\n");
    printf(" line 1 : <nSamples> <nInputs> <pOrder>\n");
    printf(" line 2 : 1 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf(" line 3 : 2 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf("          ...\n");
    printEquals(PL_INFO, 0);
  }
  fp = fopen("ps_sparse_grid_info", "r");
  if (fp == NULL)
  {
    printf("SparseGridRegression ERROR: ps_sparse_grid_info file not found.\n");
    printf("This file is used to specify information needed for \n");
    printf("sparse grid regression. The file has the format:\n");
    printf(" line 1 : <nSamples> <nInputs> <pOrder>\n");
    printf(" line 2 : 1 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf(" line 3 : 2 <input1Val> <input2Val> .. <inputnVal> <Weight>\n");
    printf("          ...\n");
    exit(1);
  }
  fscanf(fp, "%d %d %d", &nSamples_, &nInputs_, &pOrder_);
  printf("SparseGridRegression INFO: polynomial order = %d.\n", pOrder_);
  if (nSamples != nSamples_ || nInputs != nInputs_)
  {
    printf("SparseGridRegression ERROR: nSamples or nInputs does not match.\n");
    printf("                            SparseGrid is rigid in its sample.\n");
    fclose(fp);
    return;
  }
  sampleInputs_.setLength(nSamples_*nInputs_);
  sampleWeights_.setLength(nSamples_);
  for (ii = 0; ii < nSamples_; ii++)
  {
    for (kk = 0; kk < nInputs_; kk++)
    {
      fscanf(fp, "%lg", &ddata);
      sampleInputs_[ii*nInputs_+kk] = ddata;
    }
    fscanf(fp, "%lg", &ddata);
    sampleWeights_[ii] = ddata;
  }
  lBounds_.setLength(nInputs_);
  uBounds_.setLength(nInputs_);
  for (kk = 0; kk < nInputs_; kk++)
  {
    fscanf(fp, "%lg %lg", &ddata, &ddata2);
    lBounds_[kk] = ddata;
    uBounds_[kk] = ddata2;
    if (lBounds_[kk] >= uBounds_[kk])
    {
      printf("SparseGridRegression ERROR: invalid input bounds.\n");
      printf("       lbound, ubound (input %d) = %e %e\n", kk+1, 
             lBounds_[kk], uBounds_[kk]);
      exit(1);
    }
  }
  fclose(fp); 
  GenPermutations();
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SparseGridRegression::~SparseGridRegression()
{
  if (pcePerms_ != NULL)
  {
    for (int ii = 0; ii < numPerms_; ii++) 
      if (pcePerms_[ii] != NULL) delete [] pcePerms_[ii];
    delete [] pcePerms_;
  }
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int SparseGridRegression::initialize(double *X, double *Y)
{
  int totPts, ss;

  if (sampleInputs_.length() == 0)
  {
    printf("SparseGridRegression::initialize ERROR - invalid sample.\n");
    return -1;
  }
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("SparseGridRegression::initialize ERROR - invalid argument.\n");
    exit(1);
  } 
   
  analyze(X, Y);
  if (psRSCodeGen_ == 1) 
  {
    printf("SparseGridRegression INFO: response surface stand-alone ");
    printf("code not available.\n");
  }
  return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int SparseGridRegression::genNDGridData(double *X, double *Y, int *N2,
                                        double **X2, double **Y2)
{
  int totPts, ss;

  if (sampleInputs_.length() == 0)
  {
    printf("SparseGridRegression::genNDGridData ERROR - invalid sample.\n");
    return -1;
  }
  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("SparseGridRegression::genNDGridData ERROR - invalid argument.\n");
    exit(1);
  } 
   
  analyze(X, Y);

  if ((*N2) == -999) return 0;

  genNDGrid(N2, X2);
  if ((*N2) == 0) return 0;
  totPts = (*N2);

  (*Y2) = new double[totPts];
  checkAllocate(*Y2, "Y2 in SGRegression::genNDGridData");
  for (ss = 0; ss < totPts; ss++)
    (*Y2)[ss] = evaluatePoint(&((*X2)[ss*nInputs_]));

  return 0;
}

// ************************************************************************
// Generate 1D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int SparseGridRegression::gen1DGridData(double *X, double *Y, int ind1,
                                        double *settings, int *NN, 
                                        double **XX, double **YY)
{
  int    totPts, mm, nn;
  double HX, *Xloc;

  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*NN) = totPts;
  (*XX) = new double[totPts];
  (*YY) = new double[totPts];
  Xloc  = new double[nInputs_];
  checkAllocate(Xloc, "Xloc in SGRegression::gen1DGridData");
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
// Generate 2D mesh results (set all other inputs to nominal values)
// ------------------------------------------------------------------------
int SparseGridRegression::gen2DGridData(double *X, double *Y, int ind1,
                                      int ind2, double *settings, int *NN, 
                                      double **XX, double **YY)
{
  int    totPts, mm, nn, index;
  double *HX, *Xloc;

  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

  totPts = nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[2];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  (*NN) = totPts;
  (*XX) = new double[totPts * 2];
  (*YY) = new double[totPts];
  Xloc  = new double[nInputs_];
  checkAllocate(Xloc, "Xloc in SGRegression::gen2DGridData");
  for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      index = mm * nPtsPerDim_ + nn;
      Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
      Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
      (*XX)[index*2]   = Xloc[ind1];
      (*XX)[index*2+1] = Xloc[ind2];
      (*YY)[index] = evaluatePoint(Xloc);
    }
  }

  delete [] Xloc;
  delete [] HX;
  return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values) 
// ------------------------------------------------------------------------
int SparseGridRegression::gen3DGridData(double *X, double *Y, int ind1,
                                  int ind2, int ind3, double *settings, 
                                  int *NN, double **XX, double **YY)
{
  int    totPts, mm, nn, pp, index;
  double *HX, *Xloc;

  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  HX    = new double[3];
  HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  (*NN) = totPts;
  (*XX) = new double[totPts * 3];
  (*YY) = new double[totPts];
  Xloc  = new double[nInputs_];
  checkAllocate(Xloc, "Xloc in SGRegression::gen3DGridData");
  for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
  for (mm = 0; mm < nPtsPerDim_; mm++) 
  {
    for (nn = 0; nn < nPtsPerDim_; nn++)
    {
      for (pp = 0; pp < nPtsPerDim_; pp++)
      {
        index = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
        Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
        Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
        Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
        (*XX)[index*3]   = Xloc[ind1];
        (*XX)[index*3+1] = Xloc[ind2];
        (*XX)[index*3+2] = Xloc[ind3];
        (*YY)[index] = evaluatePoint(Xloc);
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
int SparseGridRegression::gen4DGridData(double *X, double *Y, int ind1, 
                             int ind2, int ind3, int ind4, double *settings, 
                             int *NN, double **XX, double **YY)
{
  int    totPts, mm, nn, pp, qq, index;
  double *HX, *Xloc;

  (*NN) = -999;
  genNDGridData(X, Y, NN, NULL, NULL);

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
  checkAllocate(Xloc, "Xloc in SGRegression::gen4DGridData");
  for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
    
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
          Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
          Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
          Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
          Xloc[ind4] = HX[3] * qq + lowerBounds_[ind4];
          (*XX)[index*4]   = Xloc[ind1];
          (*XX)[index*4+1] = Xloc[ind2];
          (*XX)[index*4+2] = Xloc[ind3];
          (*XX)[index*4+3] = Xloc[ind4];
          (*YY)[index] = evaluatePoint(Xloc);
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
double SparseGridRegression::evaluatePoint(double *X)
{
  int    ii, kk;
  double Y, multiplier, **LTable, ddata;

  if (regCoefs_.length() == 0) return 0.0;
  Y = 0.0;
  LTable = new double*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) LTable[ii] = new double[pOrder_+1];
  checkAllocate(LTable[nInputs_-1],"LTable in SGRegression::evaluatePoint");
  Y = 0.0;
  for (kk = 0; kk < numPerms_; kk++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = (X[ii] - lBounds_[ii])/(uBounds_[ii] - lBounds_[ii])*2-1;
      EvalLegendrePolynomials(ddata, LTable[ii]);
    }
    multiplier = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
      multiplier *= LTable[ii][pcePerms_[kk][ii]];
    Y += regCoefs_[kk] * multiplier;
  }
  for (ii = 0; ii < nInputs_; ii++) delete [] LTable[ii];
  delete [] LTable;
  return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double SparseGridRegression::evaluatePoint(int npts, double *X, double *Y)
{
  for (int kk = 0; kk < npts; kk++)
    Y[kk] = evaluatePoint(&X[kk*nInputs_]);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point and also its standard deviation
// ------------------------------------------------------------------------
double SparseGridRegression::evaluatePointFuzzy(double *X, double &std)
{
  printf("SparseGridRegression INFO: not implemented yet.\n");
  std = 0.0;
  return 0.0;
}

// ************************************************************************
// Evaluate a number of points and also their standard deviations
// ------------------------------------------------------------------------
double SparseGridRegression::evaluatePointFuzzy(int npts, double *X, 
                                                double *Y, double *Ystd)
{
   printf("SparseGridRegression INFO: not implemented yet.\n");
   for (int kk = 0; kk < npts; kk++) Y[kk] = 0.0;
   return 0.0;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
int SparseGridRegression::analyze(double *X, double *Y)
{
  int    ii, jj, kk, ll;
  double **LTables, *coefs, ddata, ddata2, wt;

  if (nInputs_ <= 0 || nSamples_ <= 0)
  {
    printf("SparseGridRegression::analyze ERROR - invalid arguments.\n");
    exit(1);
  } 
  if (lBounds_.length() == 0)
  {
    printf("SparseGridRegression::analyze ERROR - input bounds not set.\n");
    exit(1);
  } 
   
  for (ii = 0; ii < nSamples_; ii++)
  {
    for (kk = 0; kk < nInputs_; kk++)
    {
      ddata = sampleInputs_[ii*nInputs_+kk];
      ddata = ddata * (uBounds_[kk] - lBounds_[kk]) + lBounds_[kk];
      if (PABS(ddata - X[ii*nInputs_+kk]) > 1.0e-13)
      {
        printf("SparseGridRegression::analyze ERROR - sample mismatch\n");
        exit(1);
      }
    }
  }

  regCoefs_.setLength(numPerms_);
  LTables = new double*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) LTables[ii] = new double[pOrder_+1];
  coefs = new double[numPerms_];
  checkAllocate(coefs, "coefs in SGRegression::analyze");
  for (ii = 0; ii < numPerms_; ii++) coefs[ii] = 0.0;

  for (ii = 0; ii < nSamples_; ii++)
  {
    for (jj = 0; jj < nInputs_; jj++)
    {
      ddata = 2.0 * sampleInputs_[ii*nInputs_+jj] - 1.0;
      EvalLegendrePolynomials(ddata, LTables[jj]);
    }
    for (kk = 0; kk < numPerms_; kk++)
    {
      ddata = 1.0;
      for (ll = 0; ll < nInputs_; ll++)
        ddata *= LTables[ll][pcePerms_[kk][ll]];
      ddata2 = ddata * ddata;
      wt = sampleWeights_[ii];
      regCoefs_[kk] += (ddata * Y[ii] * wt);
      coefs[kk] += (ddata2 * wt);
    }
  } 
  if (outputLevel_ >= 0)
  {
    printf("Legendre polynomial functional forms: \n");
    printf("X normalized to Z in [-1 1] (e.g. X in [0,1] -> Z=2X-1)\n");
    printf("P_0(Z) = 1\n");
    printf("P_1(Z) = Z\n");
    printf("P_{n+1} = 1/(n+1) {(2n + 1) Z P_n(Z) + n P_{n-1}(Z)}\n");
    printEquals(PL_INFO, 0);
  }
  for (kk = 0; kk < numPerms_; kk++)
  {
    if (coefs[kk] == 0.0) 
         printf("ERROR in SparseGridRegression: divide by 0.\n");
    else regCoefs_[kk] /= coefs[kk];
    if (outputLevel_ >= 0)
    {
      printf("Legendre polynomial (");
      for (jj = 0; jj < nInputs_; jj++)
        printf("%d ", pcePerms_[kk][jj]);
      ddata = regCoefs_[kk];
      if (PABS(ddata) < 1.0e-13) ddata = 0.0;
      printf(") coefficient = %e\n", ddata);
    }
  }
  if (outputLevel_ >= 0) printAsterisks(PL_INFO, 0);

  delete [] coefs;
  for (ii = 0; ii < nInputs_; ii++) delete [] LTables[ii];
  delete [] LTables;
  return 0;
}

// *************************************************************************
// generate all combinations of a multivariate Legendre expansion
// This code is a direct translation from Burkardt's matlab code)
// -------------------------------------------------------------------------
int SparseGridRegression::GenPermutations()
{
  int  ii, kk, orderTmp, rvTmp;

  numPerms_ = computeNumPCEPermutations(nInputs_, pOrder_);

  pcePerms_ = new int*[numPerms_];
  for (ii = 0; ii < numPerms_; ii++) pcePerms_[ii] = new int[nInputs_];
  checkAllocate(pcePerms_[numPerms_-1], 
                "pcePerms_ in SGRegression::GenPermutations");

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
int SparseGridRegression::EvalLegendrePolynomials(double X, double *LTable)
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

