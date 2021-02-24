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
// Functions for the class FuncApprox  
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "Psuade.h"
#include "FuncApprox.h"
#include "Mars.h"
#include "Earth.h"
#include "GP1.h"
#include "GP2.h"
#include "GP3.h"
#include "MGP3.h"
#include "TBGP.h"
#include "SVM.h"
#include "SelectiveRegression.h"
#include "UserRegression.h"
#include "LegendreRegression.h"
#include "Regression.h"
#include "Ann.h"
#include "PWLinear.h"
#include "MarsBagg.h"
#include "SumOfTrees.h"
#include "SGRegression.h"
#include "GradLegendreRegression.h"
#include "Kriging.h"
#include "NPLearning.h"
#include "Splines.h"
#include "KNN.h"
#include "RBF.h"
#include "RBFBagg.h"
#include "MRBF.h"
#include "Acosso.h"
#include "BSSAnova.h"
#include "PsuadeRegression.h"
#include "PLS.h"
#include "MMars.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PDFBase.h"
#include "pData.h"
#include "PsuadeData.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FuncApprox::FuncApprox(int nInputs, int nSamples)
{
  if (nSamples <= 0 || nInputs <= 0) 
  {
    printf("FuncApprox::FuncApprox ERROR - invalid inputs. \n");
    printf("            nSamples = %d\n", nSamples);
    printf("            nInputs  = %d\n", nInputs);
    exit(-1);
  }
  outputLevel_ = 0;
  nSamples_    = nSamples;
  nInputs_     = nInputs;
  nPtsPerDim_  = 10;
  lowerBounds_ = new double[nInputs_];
  upperBounds_ = new double[nInputs_];
  for (int ii = 0 ; ii < nInputs_; ii++)
    lowerBounds_[ii] = upperBounds_[ii] = 0.0;
  weights_ = new double[nSamples_];;
  for (int jj = 0 ; jj < nSamples_; jj++) weights_[jj] = 1.0;
  XMeans_ = new double[nInputs_];
  XStds_ = new double[nInputs_];
  checkAllocate(XStds_, "XStds_ in FuncApprox::constructor");

  for (int ii = 0; ii < nInputs_; ii++)
  {
    XMeans_[ii] = 0.0;
    XStds_[ii] = 1.0;
  }
  YMean_ = 0.0;
  YStd_ = 1.0;
}

// ************************************************************************
// Copy constructor by Bill Oliver 
// ------------------------------------------------------------------------
FuncApprox::FuncApprox(const FuncApprox & fa)
{
  outputLevel_ = fa.outputLevel_;
  nSamples_ = fa.nSamples_;
  nInputs_ = fa.nInputs_;
  nPtsPerDim_ = fa.nPtsPerDim_;
  faID_ = fa.faID_;
  lowerBounds_ = new double[nInputs_];
  upperBounds_ = new double[nInputs_];
  weights_     = new double[nSamples_];
  checkAllocate(weights_, "weights_ in FuncApprox::constructor");
 
  for(int ii = 0; ii < nInputs_; ii++)
  {
    lowerBounds_[ii] = fa.lowerBounds_[ii];
    upperBounds_[ii] = fa.upperBounds_[ii];
  }
  for(int ii = 0; ii < nSamples_; ii++) weights_[ii] = fa.weights_[ii];
}

// ************************************************************************
// operator= by Bill Oliver 
// ------------------------------------------------------------------------
FuncApprox & FuncApprox::operator=(const FuncApprox & fa)
{
  if(this == &fa)  return *this;
  // free lowerBounds_, upperBounds_, and weights_
  delete [] lowerBounds_;
  delete [] upperBounds_;
  delete [] weights_;
  
  outputLevel_ = fa.outputLevel_;
  nSamples_ = fa.nSamples_;
  nInputs_ = fa.nInputs_;
  nPtsPerDim_ = fa.nPtsPerDim_;
  faID_ = fa.faID_;
  lowerBounds_ = new double[nInputs_];
  upperBounds_ = new double[nInputs_];
  weights_     = new double[nSamples_];
  checkAllocate(weights_, "weights_ in FuncApprox::operator=");

  for (int ii = 0; ii < nInputs_; ii++)
  {
    lowerBounds_[ii] = fa.lowerBounds_[ii];
    upperBounds_[ii] = fa.upperBounds_[ii];
  }
  for(int ii = 0; ii < nSamples_; ii++) weights_[ii] = fa.weights_[ii];
  return *this;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FuncApprox::~FuncApprox()
{
  if (lowerBounds_ != NULL) delete [] lowerBounds_;
  if (upperBounds_ != NULL) delete [] upperBounds_;
  if (weights_     != NULL) delete [] weights_;
  if (XMeans_      != NULL) delete [] XMeans_;
  if (XStds_       != NULL) delete [] XStds_;
}

// ************************************************************************
// get function approximator identifier
// ------------------------------------------------------------------------
int FuncApprox::getID()
{
  return faID_;
}

// ************************************************************************
// Set print level
// ------------------------------------------------------------------------
int FuncApprox::setOutputLevel(int level)
{
  if (level >= -1) outputLevel_ = level;
  return 0;
}

// ************************************************************************
// Set bounds for object class FuncApprox
// ------------------------------------------------------------------------
int FuncApprox::setBounds( double *lower, double *upper )
{
  for (int ii=0 ; ii<nInputs_; ii++) 
  {
    lowerBounds_[ii] = lower[ii];
    upperBounds_[ii] = upper[ii];
  }
  return 0;
}

// ************************************************************************
// load output weights
// ------------------------------------------------------------------------
int FuncApprox::loadWeights(int n, double *wgts)
{
  if (n != nSamples_)
  {
    printf("FuncApprox::loadWeights ERROR : invalid length %d.\n",n);
    exit(1);
  }
  if (weights_ != NULL) delete [] weights_;
  weights_ = NULL;
  for (int ii = 0 ; ii < n; ii++) 
  {
    if (wgts[ii] < 0.0)
    {
      printf("FuncApprox::loadWeights WARNING : weight < 0 - not used.\n");
      return 0;
    }
  }
  weights_ = new double[nSamples_];
  checkAllocate(weights_, "weights_ in FuncApprox::loadWeights");
  for (int jj = 0 ; jj < n; jj++) weights_[jj] = wgts[jj];
  return 0;
}

// ************************************************************************
// set number of points to generate in each dimension
// ------------------------------------------------------------------------
void FuncApprox::setNPtsPerDim(int npoints)
{
  if (npoints > 0) nPtsPerDim_ = npoints;
}

// ************************************************************************
// get number of points to generate in each dimension
// ------------------------------------------------------------------------
int FuncApprox::getNPtsPerDim()
{
  return nPtsPerDim_;
}

// ************************************************************************
// generate N dimensional data
// ------------------------------------------------------------------------
int FuncApprox::genNDGridData(double*,double*,int*,double**,double**) 
{
  return -1;
}

// ************************************************************************
// generate 1 dimensional data
// ------------------------------------------------------------------------
int FuncApprox::gen1DGridData(double*,double*,int,double*,int*,double**,
                              double**) 
{
  return -1;
}

// ************************************************************************
// generate 2 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen2DGridData(double*,double*,int,int,double*,int*,
                              double**,double**) 
{
  return -1;
}

// ************************************************************************
// generate 3 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen3DGridData(double*,double*,int,int,int,double*,int*, 
                              double**,double**) 
{
  return -1;
}

// ************************************************************************
// generate 4 dimensional surface data
// ------------------------------------------------------------------------
int FuncApprox::gen4DGridData(double*,double*,int,int,int,int,double*, 
                              int*,double**,double**) 
{
  return -1;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double FuncApprox::evaluatePointFuzzy(double *, double &std)
{
  std = 0.0;
  return 0.0;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double FuncApprox::evaluatePointFuzzy(int npts, double *, double *Y,
                                      double *Ystd)
{
  for (int ii = 0; ii < npts; ii++)
  {
    Y[ii]  = 0.0;
    Ystd[ii]  = 0.0;
  }
  return 0.0;
}

// ************************************************************************
// Set parameters
// ------------------------------------------------------------------------
double FuncApprox::setParams(int, char **)
{
  return -1.0;
}

// ************************************************************************
// initialize scaling for inputs
// ------------------------------------------------------------------------
int FuncApprox::initInputScaling(double *XIn, double *XOut, int flag)
{
  int    ii, jj;
  double ddata;
  char   pString[500], response[100];
                                                                                
  if (XMeans_ != NULL) delete [] XMeans_;
  if (XStds_  != NULL) delete [] XStds_;
  XMeans_ = new double[nInputs_];
  XStds_ = new double[nInputs_];
  checkAllocate(XStds_, "XStds_ in FuncApprox::initInputScaling");
  for (ii = 0; ii < nInputs_; ii++)
  {
    XMeans_[ii] = 0.0;
    XStds_[ii] = 1.0;
  }
  for (ii = 0; ii < nInputs_*nSamples_; ii++) XOut[ii] = XIn[ii];
  if (flag == 1)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++) ddata += XIn[jj*nInputs_+ii];
      XMeans_[ii] = ddata / (double) nSamples_;
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
        ddata += pow(XIn[jj*nInputs_+ii] - XMeans_[ii], 2.0);
      XStds_[ii] = sqrt(ddata / (double) (nSamples_ - 1));
      if (XStds_[ii] == 0.0) XStds_[ii] = 1.0;
      for (jj = 0; jj < nSamples_; jj++)
        XOut[jj*nInputs_+ii] = (XIn[jj*nInputs_+ii]-XMeans_[ii])/
                                XStds_[ii];
      if (outputLevel_ > 3)
        printf("Input %d scaling info : mean, std = %e %e\n",ii+1,
                XMeans_[ii], XStds_[ii]);
    }
  }
  else if (psRSExpertMode_ == 1)
  {
    sprintf(pString, "Scale the sample matrix ? (y or n) ");
    getString(pString, response);
    if (response[0] == 'y')
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = 0.0;
        for (jj = 0; jj < nSamples_; jj++) ddata += XIn[jj*nInputs_+ii];
        XMeans_[ii] = ddata / (double) nSamples_;
        ddata = 0.0;
        for (jj = 0; jj < nSamples_; jj++)
          ddata += pow(XIn[jj*nInputs_+ii] - XMeans_[ii], 2.0);
        XStds_[ii] = sqrt(ddata / (double) (nSamples_ - 1));
        if (XStds_[ii] == 0.0) XStds_[ii] = 1.0;
        for (jj = 0; jj < nSamples_; jj++)
          XOut[jj*nInputs_+ii] = (XIn[jj*nInputs_+ii]-XMeans_[ii])/
                                  XStds_[ii];
        printf("Input %d scaling info : mean, std = %e %e\n",ii+1,
               XMeans_[ii], XStds_[ii]);
      }
    }
  }
  return 0;
}

// ************************************************************************
// initialize scaling for output
// ------------------------------------------------------------------------
int FuncApprox::initOutputScaling(double *YIn, double *YOut)
{
  int    ii;
  double ddata;
                                                                                
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++) ddata += YIn[ii];
  YMean_ = ddata / (double) nSamples_;
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++)
     ddata += pow(YIn[ii] - YMean_, 2.0);
  YStd_ = sqrt(ddata / (double) (nSamples_ - 1));
  if (YStd_ == 0.0) YStd_ = 1.0;
  for (ii = 0; ii < nSamples_; ii++)
    YOut[ii] = (YIn[ii] - YMean_) / YStd_;
  return 0;
}

// ************************************************************************
// generate m-dimensional grid data
// ------------------------------------------------------------------------
int FuncApprox::genNDGrid(int *nPts, double **XOut) 
{
  int    ii, mm, totPts;
  psVector vecHX, vecXT;

  if (nInputs_ > 21)
  {
    printf("FuncApprox genNDGrid INFO: nInputs > 21 not supported.\n");
    (*nPts) = 0;
    (*XOut) = NULL;
    return 0;
  }

  if (nInputs_ == 21 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 20 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 19 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 18 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 17 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 16 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 15 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 14 && nPtsPerDim_ >    2) nPtsPerDim_ =  2;
  if (nInputs_ == 13 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
  if (nInputs_ == 12 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
  if (nInputs_ == 11 && nPtsPerDim_ >    3) nPtsPerDim_ =  3;
  if (nInputs_ == 10 && nPtsPerDim_ >    4) nPtsPerDim_ =  4;
  if (nInputs_ ==  9 && nPtsPerDim_ >    5) nPtsPerDim_ =  5;
  if (nInputs_ ==  8 && nPtsPerDim_ >    6) nPtsPerDim_ =  6;
  if (nInputs_ ==  7 && nPtsPerDim_ >    8) nPtsPerDim_ =  8;
  if (nInputs_ ==  6 && nPtsPerDim_ >   10) nPtsPerDim_ = 10;
  if (nInputs_ ==  5 && nPtsPerDim_ >   16) nPtsPerDim_ = 16;
  if (nInputs_ ==  4 && nPtsPerDim_ >   32) nPtsPerDim_ = 32;
  if (nInputs_ ==  3 && nPtsPerDim_ >   64) nPtsPerDim_ = 64;
  if (nInputs_ ==  2 && nPtsPerDim_ > 1024) nPtsPerDim_ = 1024;
  if (nInputs_ ==  1 && nPtsPerDim_ > 8192) nPtsPerDim_ = 8192;
  totPts = nPtsPerDim_;
  for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
  vecHX.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    vecHX[ii] = (upperBounds_[ii] - lowerBounds_[ii]) /
                 (double) (nPtsPerDim_ - 1);

  (*XOut) = new double[nInputs_ * totPts];
  (*nPts) = totPts;
  vecXT.setLength(nInputs_);

  for (ii = 0; ii < nInputs_; ii++) vecXT[ii] = lowerBounds_[ii];

  for (mm = 0; mm < totPts; mm++)
  {
    for (ii = 0; ii < nInputs_; ii++ ) (*XOut)[mm*nInputs_+ii] = vecXT[ii];
    for (ii = 0; ii < nInputs_; ii++ )
    {
      vecXT[ii] += vecHX[ii];
      if (vecXT[ii] < upperBounds_[ii] ||
           PABS(vecXT[ii] - upperBounds_[ii]) < 1.0E-7) break;
      else vecXT[ii] = lowerBounds_[ii];
    }
  }
  return 0;
}

// ************************************************************************
// ************************************************************************
// ************************************************************************
// friend function (print current function approximator)
// ------------------------------------------------------------------------
extern "C" 
int getFAType(char *pString)
{
  int faType;

  faType = getInt(0, PSUADE_NUM_RS-1, pString);
#ifndef HAVE_MARS
  if (faType == PSUADE_RS_MARS)  faType = -1;
  if (faType == PSUADE_RS_MARSB) faType = -1;
  if (faType == PSUADE_RS_MMARS) faType = -1;
#endif
#ifndef HAVE_SNNS
  if (faType == 5) faType = -1;
#endif
#ifndef HAVE_TPROS
  if (faType == PSUADE_RS_GP1) faType = PSUADE_RS_GP2;
#endif
#ifndef HAVE_SVM
  if (faType == PSUADE_RS_SVM) faType = -1;
#endif
#ifndef HAVE_EARTH
  if (faType == PSUADE_RS_EARTH) faType = -1;
#endif
#ifndef HAVE_TGP
  if (faType == PSUADE_RS_TGP) faType = -1;
#endif
  return faType;
}

// ************************************************************************
// friend function (print current function approximator)
// ------------------------------------------------------------------------
extern "C" 
void printThisFA(int faType)
{
  switch( faType )
  {
#ifdef HAVE_MARS
    case PSUADE_RS_MARS: 
         printf("MARS model\n"); 
         break;
#else
    case PSUADE_RS_MARS: 
         printf("MARS model (not installed)\n"); 
         break;
#endif
    case PSUADE_RS_REGR1: 
         printf("Linear regression model\n"); 
         break;
    case PSUADE_RS_REGR2: 
         printf("Quadratic regression model\n"); 
         break;
    case PSUADE_RS_REGR3: 
         printf("Cubic regression model\n"); 
         break;
    case PSUADE_RS_REGR4: 
         printf("Quartic regression model\n"); 
         break;
#ifdef HAVE_SNNS
    case PSUADE_RS_ANN: 
         printf("Artificial neural network model\n"); 
         break;
#else
    case PSUADE_RS_ANN: 
         printf("Artificial neural network model (not installed)\n"); 
         break;
#endif
    case PSUADE_RS_REGRS: 
         printf("User-defined regression model\n"); 
         break;
#ifdef HAVE_TPROS
    case PSUADE_RS_GP1: 
         printf("Gaussian process (MacKay) model\n"); 
         break;
#else
    case PSUADE_RS_GP1: 
         printf("Gaussian process (MacKay) model (not installed)\n");
         break;
#endif
    case PSUADE_RS_GP2: 
         printf("Gaussian process (Tong) model\n"); 
         break;
#ifdef HAVE_SVM
    case PSUADE_RS_SVM: 
         printf("SVM-light (Joachims) model\n"); 
         break;
#else
    case PSUADE_RS_SVM: 
         printf("SVM-light (Joachims) model (not installed)\n");
         break;
#endif
    case PSUADE_RS_REGRGL: 
         printf("Derivative-based Legendre polynomial regression\n"); 
         break;
#ifdef HAVE_TGP
    case PSUADE_RS_TGP: 
         printf("Tree-based Gaussian Process\n"); 
         break;
#else
    case PSUADE_RS_TGP: 
         printf("Tree-based Gaussian Process (not installed)\n");
         break;
#endif
#ifdef HAVE_MARS
    case PSUADE_RS_MARSB: 
         printf("MARS with bagging model\n");
         break;
#else
    case PSUADE_RS_MARSB: 
         printf("MARS with bagging model (not installed)\n");
         break;
#endif
#ifdef HAVE_EARTH
    case PSUADE_RS_EARTH: 
         printf("Earth model\n"); 
         break;
#else
    case PSUADE_RS_EARTH: 
         printf("Earth model (not installed)\n"); 
         break;
#endif
    case PSUADE_RS_SOTS: 
         printf("Sum-of-trees model\n"); 
         break;
    case PSUADE_RS_REGRL: 
         printf("Legendre polynomial regression\n"); 
         break;
    case PSUADE_RS_REGRU: 
         printf("User-defined (nonpolynomial) regression\n"); 
         break;
    case PSUADE_RS_REGSG: 
         printf("Sparse Grid polynomial regression\n"); 
         break;
    case PSUADE_RS_KR: 
         printf("Kriging\n"); 
         break;
    case PSUADE_RS_SPLINES: 
         printf("Splines on regular grid (1D, 2D, or 3D only)\n");
         break;
    case PSUADE_RS_KNN: 
         printf("K-nearest neighbor\n"); 
         break;
    case PSUADE_RS_RBF: 
         printf("Radial Basis Function\n"); 
         break;
    case PSUADE_RS_RBFB: 
         printf("Radial Basis Function with bagging\n"); 
         break;
    case PSUADE_RS_MRBF: 
         printf("Multi-Radial Basis Function\n"); 
         break;
    case PSUADE_RS_MGP2: 
         printf("Multi-Gaussian process (Tong)\n"); 
         break;
#ifdef HAVE_MARS
    case PSUADE_RS_MMARS: 
         printf("Multiple MARS model\n");
         break;
#else
    case PSUADE_RS_MMARS: 
         printf("Multiple MARS model (MARS not installed)\n");
         break;
#endif
  }
}

// ************************************************************************
// friend function (print function approximator information)
// ------------------------------------------------------------------------
extern "C" 
int writeFAInfo(int level)
{
  printDashes(PL_INFO, 0);
  printf("Available response surface tools: \n");
  printDashes(PL_INFO, 0);
  if (level > 0)
  {
   printf("Expert advices: \n");
#ifdef HAVE_MARS
   printf(" MARS - may have accuracy problem near domain boundary. Use\n"); 
   printf("   this option if sample size is sufficiently large (>100).\n"); 
#endif
   printf(" LINEAR, QUADRATIC, CUBIC, QUARTIC - good for small sample\n");
   printf("   sizes; and when the function is sufficiently smooth.\n");
   printf("   For higher than fourth order, use LEGENDRE (option 15) with\n");
   printf("   response surface expert mode turned on to select order.\n");
#ifdef HAVE_SNNS
   printf(" ANN - supported but currently not maintained.\n");
#endif
   printf(" SELECTIVE POLYNOMIAL REGRESSION - for high order polynomials\n");
   printf("   but your sample size is too small. (So you select certain\n");
   printf("   terms, provided you know which ones.)\n");

   printf(" GAUSSIAN PROCESS - may  encounter non-definite covariance\n");
   printf("   matrix problem. However, if no such problem appears, it\n");
   printf("   especially good for small samples (a few to a few tens).\n");
   printf("   GP is relatively slow, so for sample sizes of more than a\n");
   printf("   few hundred, be patiet. For nonsmooth functions, try TGP.\n");
#ifdef HAVE_SVM
   printf(" SVM - provides 3 options: turn on rs_expert to select. Also,\n");
   printf("   use svmfind to search for best settings.\n");
#endif
   printf(" BOOTSTRAPPED MARS - intended to be used with adaptive sample\n");
   printf("   refinement, which adds more sample points near the boundary\n");
   printf("   of the parameter space.\n");
   printf(" SUM-OF-TREES REGRESSION - usually gives non-smooth response\n");
   printf("   responses. It is provided here for completeness, but is not\n");
   printf("   generally recommended.\n");
   printf(" SPARSE GRID REGRESSION - has to use sparse grid designs. Also,\n");
   printf("   you cannot use cross validation on sparse grid regression.\n");
   printf("   Use a test set (rstest) to validate your response surface.\n");
   printf(" KRIGING - This is another form of the Gaussian process which\n");
   printf("   uses deterministic optimization to compute hyperparameters.\n");
   printf("   This method is good for up to about 2000 sample points;\n");
   printf("   otherwise it may be computationally expensive.\n");
   printf(" SPLINES - currently only supports 1D, 2D, or 3D. This method\n");
   printf("   works only with full factorial designs. Also, you cannot use\n");
   printf("   cross validation with this method.\n");
   printf("   Use a test set (rstest) to validate your response surface.\n");
   printf(" K-NEAREST NEIGHBOR - for large data set when data points are\n");
   printf("   relatively close to one another, this may be useful.\n");
   printf(" RBF  - for small to medium data set (too expensive otherwise).\n");
   printf(" MRBF - for larger data set (> 2000).\n");
   printf(" MGP2 - (Gaussian process) for larger data set (> 2000).\n");
   printf(" MMARS - (multiple MARS) for even larger data set (> 20000).\n");
  }
#ifdef HAVE_MARS
  printDashes(PL_INFO, 0);
  printf("0. MARS \n");
#endif
  printf("1. Linear regression \n");
  printf("2. Quadratic regression \n");
  printf("3. Cubic regression \n");
  printf("4. Quartic regression \n");
#ifdef HAVE_SNNS
  printf("5. Artificial neural network \n");
#endif
  printf("6. Selective polynomial regression \n");
#ifdef HAVE_TPROS
  printf("7. Gaussian process (MacKay)\n");
#endif
  printf("8. Gaussian process (Tong)\n");
#ifdef HAVE_SVM
  printf("9. SVM-light (Joachims)\n");
#endif
  printf("10. Derivative-based Legendre polynomial regression\n");
#ifdef HAVE_TGP
  printf("11. Tree-based Gaussian Process\n");
#endif
#ifdef HAVE_MARS
  printf("12. MARS with bootstrap aggregating (bagging)\n");
#endif
#ifdef HAVE_EARTH
  printf("13. Earth (another MARS)\n");
#endif
  printf("14. Sum-of-trees model\n");
  printf("15. Legendre polynomial regression\n");
  printf("16. User-defined (nonpolynomial) regression\n");
  printf("17. Sparse Grid polynomial regression\n"); 
  printf("18. Kriging\n"); 
  printf("19. Splines on regular grid (1D, 2D, or 3D only)\n");
  printf("20. K nearest neighbors \n");
  printf("21. Radial Basis Function\n");
  printf("22. Acosso (by Storlie, LANL. Need R to run)\n");
  printf("23. BSSAnova (by Storlie, LANL. Need R to run)\n");
  printf("24. Radial Basis Function with bagging\n");
  printf("25. Partial Least Squares Linear Regression (PLS)\n");
  printf("26. Multi-Radial Basis Function (for large samples)\n");
  printf("27. Multi-Gaussian process (Tong, for large samples)\n");
  printf("28. Multi-MARS (for large samples)\n");
  return PSUADE_NUM_RS;
}

// ************************************************************************
// friend function (create a function approximator from a few parameters)
// ------------------------------------------------------------------------
extern "C" 
FuncApprox *genFA(int faType, int nInputs, int outLevel, int nSamples)
{
  int        rsType, nsize;
  FuncApprox *faPtr=NULL;
  char       *params[1], winput[10000], *strPtr, equal[100];

  if (faType >= 0) rsType = faType;
  else
  {
    rsType = -1;
    while (rsType < 0 || rsType >= PSUADE_NUM_RS)
    {
      writeFAInfo(outLevel);
      sprintf(winput, "Please enter your choice ? ");
      rsType = getInt(0, PSUADE_NUM_RS, winput);
    }
  }

  if      (rsType == PSUADE_RS_MARS) faPtr = new Mars(nInputs, nSamples);
  else if (rsType == PSUADE_RS_ANN)  faPtr = new Ann(nInputs, nSamples);
  else if (rsType == PSUADE_RS_REGRS)
          faPtr = new SelectiveRegression(nInputs, nSamples);
  else if (rsType == PSUADE_RS_GP1) faPtr = new GP1(nInputs, nSamples);
  else if (rsType == PSUADE_RS_GP2) 
  {
    faPtr = NULL;
    if (psConfig_ != NULL)
    {
      strPtr = psConfig_->getParameter("RS_no_multi_domain");
      if (strPtr != NULL) faPtr = new RBF(nInputs, nSamples);
    }
    if (faPtr == NULL)
    {
      nsize = 2000; 
      if (psConfig_ != NULL)
      {
        strPtr = psConfig_->getParameter("MGP_max_samples_per_group");
        if (strPtr != NULL)
        {
          sscanf(strPtr, "%s %s %d", winput, equal, &nsize);
          if (nsize < 100) nsize = 1000;
        }
      }
      if (nSamples > nsize)
      {
        printf("nSamples > %d ==> switch to MGP3\n", nsize);
        printf("NOTE: to remain in GP3, set RS_no_multi_domain\n");
        faPtr = new MGP3(nInputs, nSamples);
      }
      else faPtr = new GP3(nInputs, nSamples);
    }
  }
  else if (rsType == PSUADE_RS_SVM) faPtr = new SVM(nInputs, nSamples);
  else if (rsType == PSUADE_RS_REGRGL)
          faPtr = new GradLegendreRegression(nInputs, nSamples);
  else if (rsType == PSUADE_RS_TGP)   faPtr = new TGP(nInputs, nSamples);
  else if (rsType == PSUADE_RS_MARSB) faPtr = new MarsBagg(nInputs, nSamples);
  else if (rsType == PSUADE_RS_EARTH) faPtr = new Earth(nInputs, nSamples);
  else if (rsType == PSUADE_RS_SOTS) faPtr = new SumOfTrees(nInputs,nSamples);
  else if (rsType == PSUADE_RS_REGRL)
          faPtr = new LegendreRegression(nInputs,nSamples);
  else if (rsType == PSUADE_RS_REGRU)
          faPtr = new UserRegression(nInputs, nSamples);
  else if (rsType == PSUADE_RS_REGSG)
          faPtr = new SparseGridRegression(nInputs, nSamples);
  else if (rsType == PSUADE_RS_KR)
          faPtr = new Kriging(nInputs, nSamples);
  else if (rsType == PSUADE_RS_SPLINES)
  {
    if (nInputs > 3)
    {
      printf("genFA ERROR: Splines does not support nInputs > 3.\n");
      exit(1);
    }
    faPtr = new Splines(nInputs, nSamples);
  }
  else if (rsType == PSUADE_RS_KNN) faPtr = new KNN(nInputs, nSamples);
  else if (rsType == PSUADE_RS_RBF)
  {
    faPtr = NULL;
    if (psConfig_ != NULL)
    {
      strPtr = psConfig_->getParameter("RS_no_multi_domain");
      if (strPtr != NULL) faPtr = new RBF(nInputs, nSamples);
    }
    if (faPtr == NULL)
    {
      nsize = 5000; 
      if (psConfig_ != NULL)
      {
        strPtr = psConfig_->getParameter("MRBF_max_samples_per_group");
        if (strPtr != NULL)
        {
          sscanf(strPtr, "%s %s %d", winput, equal, &nsize);
          if (nsize < 100) nsize = 2000;
        }
      }
      if (nSamples > nsize)
      {
        printf("nSamples > %d ==> switch to MRBF\n", nsize);
        printf("NOTE: to remain in RBF, set RS_no_multi_domain\n");
        faPtr = new MRBF(nInputs, nSamples);
      }
      else faPtr = new RBF(nInputs, nSamples);
    }
  }
  else if (rsType == PSUADE_RS_ACOSSO)
       faPtr = new Acosso(nInputs, nSamples);
  else if (rsType == PSUADE_RS_BSSANOVA)
       faPtr = new BSSAnova(nInputs, nSamples);
  else if (rsType == PSUADE_RS_RBFB)
       faPtr = new RBFBagg(nInputs, nSamples);
  else if (rsType == PSUADE_RS_PLS)
       faPtr = new PLS(nInputs, nSamples);
  else if (rsType == PSUADE_RS_MRBF)
       faPtr = new MRBF(nInputs, nSamples);
  else if (rsType == PSUADE_RS_MGP2)
       faPtr = new MGP3(nInputs, nSamples);
  else if (rsType == PSUADE_RS_MMARS)
       faPtr = new MMars(nInputs, nSamples);
  else if (faType == PSUADE_RS_LOCAL)
       faPtr = new PsuadeRegression(nInputs, nSamples);
  else if (rsType == PSUADE_RS_NPL)
       faPtr = new NPLearning(nInputs, nSamples);
  else
  {
    printf("INFO: rstype has been set to default = regression.\n");
    faPtr = new Regression(nInputs, nSamples);
    params[0] = (char *) &rsType;
    faPtr->setParams(1, params);
  }
  return faPtr;
}

// ************************************************************************
// friend function (create a function approximator from a data file)
// ------------------------------------------------------------------------
extern "C" 
FuncApprox *genFAInteractive(PsuadeData *psuadeIO, int flag)
{
  int        faType, nInputs, nSamples, nOutputs, wgtID, ii, nPtsPerDim;
  int        totPts, outputID, printLevel, nsize;
  double     *wghts, *Y;
  FuncApprox *faPtr;
  char       *params[3], winput[5001], equal[100], *strPtr;
  pData      pPtr, pInputs, pOutputs, pStates, pLower, pUpper;

  if (psuadeIO == NULL)
  {
    printf("ERROR: PsuadeData does not exist.\n");
    return NULL;
  }

  if ((flag & 1) == 0)
  {
    assert(psuadeIO->getParameter("ana_rstype", pPtr) == 0);
    faType = pPtr.intData_;
    if (faType < 0 || faType >= PSUADE_NUM_RS)
    {
      printf("createFA : faType (%d) not valid.\n", faType);
      exit(1);
    }
  }
  else
  {
    faType = -1;
    while (faType < 0 || faType >= PSUADE_NUM_RS)
    {
      writeFAInfo(-1);
      sprintf(winput, "Please enter your choice ? ");
      faType = getInt(0, PSUADE_NUM_RS-1, winput);
    }
  }

  psuadeIO->getParameter("input_ninputs", pPtr);
  nInputs = pPtr.intData_;
  psuadeIO->getParameter("output_noutputs", pPtr);
  nOutputs = pPtr.intData_;
  psuadeIO->getParameter("method_nsamples", pPtr);
  nSamples = pPtr.intData_;
  if (nSamples > psFAMaxDataPts_)
  {
    printf("PSUADE WARNING: For nSamples > %d,\n", psFAMaxDataPts_);
    printf("  it can be extremely expensive to create\n");
    printf("  response surfaces.\n");
    printf("  Please consult PSUADE developers before moving on.\n");
    exit(1);
  }
  psuadeIO->getParameter("ana_outputid", pPtr);
  outputID = pPtr.intData_;
  psuadeIO->getParameter("ana_diagnostics", pPtr);
  printLevel = pPtr.intData_;
  psuadeIO->getParameter("ana_regressionwgtid", pPtr);
  wgtID = pPtr.intData_;
  psuadeIO->getParameter("output_sample", pOutputs);
  psuadeIO->getParameter("input_sample", pInputs);
  psuadeIO->getParameter("output_states", pStates);
  psuadeIO->getParameter("input_lbounds", pLower);
  psuadeIO->getParameter("input_ubounds", pUpper);

  if      (faType == PSUADE_RS_MARS) faPtr = new Mars(nInputs, nSamples);
  else if (faType == PSUADE_RS_ANN)  faPtr = new Ann(nInputs, nSamples);
  else if (faType == PSUADE_RS_REGRS)
          faPtr = new SelectiveRegression(nInputs, nSamples);
  else if (faType == PSUADE_RS_GP1)   faPtr = new GP1(nInputs, nSamples);
  else if (faType == PSUADE_RS_GP2)
  {
    faPtr = NULL;
    if (psConfig_ != NULL)
    {
      strPtr = psConfig_->getParameter("RS_no_multi_domain");
      if (strPtr != NULL) faPtr = new RBF(nInputs, nSamples);
    }
    if (faPtr == NULL)
    {
      nsize = 2000; 
      if (psConfig_ != NULL)
      {
        strPtr = psConfig_->getParameter("MGP_max_samples_per_group");
        if (strPtr != NULL)
        {
          sscanf(strPtr, "%s %s %d", winput, equal, &nsize);
          if (nsize < 100) nsize = 1000;
        }
      }
      if (nSamples > nsize)
      {
        printf("nSamples > %d ==> switch to MGP3\n", nsize);
        faPtr = new MGP3(nInputs, nSamples);
      }
      else faPtr = new GP3(nInputs, nSamples);
    }
  }
  else if (faType == PSUADE_RS_SVM)   faPtr = new SVM(nInputs, nSamples);
  else if (faType == PSUADE_RS_REGRGL)
          faPtr = new GradLegendreRegression(nInputs, nSamples);
  else if (faType == PSUADE_RS_TGP)   faPtr = new TGP(nInputs, nSamples);
  else if (faType == PSUADE_RS_MARSB) faPtr = new MarsBagg(nInputs, nSamples);
  else if (faType == PSUADE_RS_EARTH) faPtr = new Earth(nInputs, nSamples);
  else if (faType == PSUADE_RS_SOTS) faPtr = new SumOfTrees(nInputs,nSamples);
  else if (faType == PSUADE_RS_REGRL)
  {
    faPtr = new LegendreRegression(nInputs,nSamples);
    psuadeIO->getParameter("ana_poly_order", pPtr);
    ii = pPtr.intData_;
    if (ii > 0)
    {
      params[0] = (char *) &ii;
      faPtr->setParams(1, params);
    }
  }
  else if (faType == PSUADE_RS_REGRU)
          faPtr = new UserRegression(nInputs, nSamples);
  else if (faType == PSUADE_RS_REGSG)
          faPtr = new SparseGridRegression(nInputs, nSamples);
  else if (faType == PSUADE_RS_KR)
          faPtr = new Kriging(nInputs, nSamples);
  else if (faType == PSUADE_RS_SPLINES)
  {
    if (nInputs > 3)
    {
      printf("genFA ERROR: Splines does not support nInputs > 3.\n");
      exit(1);
    }
    faPtr = new Splines(nInputs, nSamples);
  }
  else if (faType == PSUADE_RS_KNN) faPtr = new KNN(nInputs, nSamples);
  else if (faType == PSUADE_RS_RBF)
  {
    faPtr = NULL;
    if (psConfig_ != NULL)
    {
      strPtr = psConfig_->getParameter("RS_no_multi_domain");
      if (strPtr != NULL) faPtr = new RBF(nInputs, nSamples);
    }
    if (faPtr == NULL)
    {
      nsize = 5000; 
      if (psConfig_ != NULL)
      {
        strPtr = psConfig_->getParameter("MRBF_max_samples_per_group");
        if (strPtr != NULL)
        {
          sscanf(strPtr, "%s %s %d", winput, equal, &nsize);
          if (nsize < 100) nsize = 2000;
        }
      }
      if (nSamples > nsize)
      {
        printf("nSamples > %d ==> switch to MRBF\n", nsize);
        faPtr = new MRBF(nInputs, nSamples);
      }
      else faPtr = new RBF(nInputs, nSamples);
    }
  }
  else if (faType == PSUADE_RS_ACOSSO)
       faPtr = new Acosso(nInputs, nSamples);
  else if (faType == PSUADE_RS_BSSANOVA)
       faPtr = new BSSAnova(nInputs, nSamples);
  else if (faType == PSUADE_RS_RBFB)
       faPtr = new RBFBagg(nInputs, nSamples);
  else if (faType == PSUADE_RS_PLS)
       faPtr = new PLS(nInputs, nSamples);
  else if (faType == PSUADE_RS_MRBF)
       faPtr = new MRBF(nInputs, nSamples);
  else if (faType == PSUADE_RS_MMARS)
       faPtr = new MMars(nInputs, nSamples);
  else if (faType == PSUADE_RS_LOCAL)
       faPtr = new PsuadeRegression(nInputs, nSamples);
  else if (faType == PSUADE_RS_NPL)
       faPtr = new NPLearning(nInputs, nSamples);
  else
  {
    faPtr = new Regression(nInputs, nSamples);
    params[0] = (char *) &faType;
    faPtr->setParams(1, params);
  }
  faPtr->setBounds(pLower.dbleArray_, pUpper.dbleArray_);
  faPtr->setOutputLevel(printLevel);

  nPtsPerDim = 256;
  totPts = 1000001;
  while (totPts > 1000000)
  {
    nPtsPerDim = nPtsPerDim / 2;
    totPts = 1;
    for (ii = 0; ii < nInputs; ii++)
    {
      totPts *= nPtsPerDim;
      if (totPts > 1000000) break;
    }
  }
  faPtr->setNPtsPerDim(nPtsPerDim);

  if (wgtID >= 0 && wgtID < nOutputs)
  {
    wghts = new double[nSamples];
    for (ii = 0; ii < nSamples; ii++)
      wghts[ii] = pOutputs.dbleArray_[nOutputs*ii+wgtID];
    faPtr->loadWeights(nSamples, wghts);
    delete [] wghts;
  }

  if (flag & 2)
  {
    Y = new double[nSamples];
    checkAllocate(Y, "Y in FuncApprox::genFAInteractive");
    for (ii = 0; ii < nSamples; ii++)
       Y[ii] = pOutputs.dbleArray_[nOutputs*ii+outputID];
    faPtr->initialize(pInputs.dbleArray_, Y);
    delete [] Y;
  }
  return faPtr;
}

// ************************************************************************
// friend function (create a function approximator given a file name)
// perform PDF transformation
// check invalid sample points
// RS type from file 
// ------------------------------------------------------------------------
extern "C" 
FuncApprox *genFAFromFile(char *fname, int outputID)
{
  int        ii, nInputs, nOutputs, nSamples, status, *sampleStates;
  double     *sampleInputs, *sampleOutputs;
  PsuadeData *psuadeIO = new PsuadeData();
  FuncApprox *faPtr;
  pData      pPtr, pInpData, pOutData, pStates;

  for (ii = strlen(fname)-1; ii >= 0; ii--) if (fname[ii] == '/') break;
  status = psuadeIO->readPsuadeFile(fname);
  if (status != 0)
  {
    delete psuadeIO;
    return NULL;
  }

  psuadeIO->getParameter("input_ninputs", pPtr);
  nInputs = pPtr.intData_;
  psuadeIO->getParameter("output_noutputs", pPtr);
  nOutputs = pPtr.intData_;
  psuadeIO->getParameter("method_nsamples", pPtr);
  nSamples = pPtr.intData_;
  psuadeIO->getParameter("input_sample", pInpData);
  sampleInputs = pInpData.dbleArray_;
  psuadeIO->getParameter("output_sample", pOutData);
  sampleOutputs = pOutData.dbleArray_;
  psuadeIO->getParameter("output_states", pStates);
  sampleStates = pStates.intArray_;

  for (ii = 0; ii < nSamples; ii++)
  {
    if (sampleStates[ii] != 1 || 
        sampleOutputs[nOutputs*ii+outputID] == PSUADE_UNDEFINED)
    {
      printf("FuncApprox::genRSModel ERROR - invalid output.\n");
      printf("  Advice: check your sample data file to see if there is\n");
      printf("          any sample output having the value 9.999e+34,\n");
      printf("          any sample status flag not equal to 1.\n");
      exit(1);
    }
  }

  PDFTransform(psuadeIO, nSamples, nInputs, sampleInputs);
  psuadeIO->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                               sampleInputs,NULL,NULL,NULL,NULL,NULL);
  psuadeIO->getParameter("ana_outputid", pPtr);
  ii = pPtr.intData_;
  psuadeIO->updateAnalysisSection(-1,-1,-1,-1,outputID,-1);
  faPtr = genFAInteractive(psuadeIO, 2);
  psuadeIO->updateAnalysisSection(-1,-1,-1,-1,ii,-1);
  delete psuadeIO;
  return faPtr;
}

