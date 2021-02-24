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
// Functions for the class GP3
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#ifdef WINDOWS
#include <windows.h>
#undef ERROR  //Windows already has a macro defined as ERROR
#undef IS_ERROR
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "GP3.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "Sampling.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// ************************************************************************
// external functions
// ------------------------------------------------------------------------
extern "C" {
  void dpotrf_(char *, int *, double *, int *, int *);
  void dpotrs_(char *, int *, int *, double *, int *, double *,int *,int *);
  void newuoa_(int *,int *,double *,double *,double *,int *,int *,double*);
#ifdef HAVE_LBFGS
#include "../../External/L-BFGS-B-C/src/lbfgsb.h"
#endif
}
GP3 *GP3Obj=NULL;
psVector GP3_OptX;
double   GP3_OptY=1e35;

// ************************************************************************
// resident function perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" 
{
#endif
  void *gp3newuoaevalfunc_(int *nInps, double *XValues, double *YValue)
  {
    int status;
    psVector VecG;
    if (GP3Obj == NULL)
    {
      printf("GP3 ERROR: no GP3 object in function evalution\n");
      exit(1);
    }
    psVector VecP;
    VecP.load(*nInps, XValues);
    (*YValue) = GP3Obj->computeGradients(VecP, VecG, status);
    if ((*YValue) < GP3_OptY)
    {
      GP3_OptY = *YValue;
      GP3_OptX = VecP;
    }
    return NULL;
  }    
#ifdef __cplusplus
}
#endif

// ************************************************************************
// Constructor for object class GP3
// ------------------------------------------------------------------------
GP3::GP3(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  optLinTerm_ = 0;
  faID_ = PSUADE_RS_GP2;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GP3::~GP3()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int GP3::initialize(double *XIn, double *YIn)
{
  XDataN_.setLength(nSamples_*nInputs_);
  initInputScaling(XIn, XDataN_.getDVector(), 1);
  YData_.setLength(nSamples_);
  initOutputScaling(YIn, YData_.getDVector());
   
  if (outputLevel_ > 1) printf("GP3 training begins....\n");
  train();
  if (outputLevel_ > 1) printf("GP3 training completed.\n");
  if (psRSCodeGen_ == 1) genCode();
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int GP3::genNDGridData(double *XIn, double *YIn, int *NOut, double **XOut, 
                      double **YOut)
{
  int    totPts;
  double *XX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  genNDGrid(NOut, &XX);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  (*YOut) = new double[totPts];
  evaluatePoint(totPts, XX, (*YOut));
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*XOut) = XX;
  return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int GP3::gen1DGridData(double *XIn, double *YIn, int ind1,double *settings, 
                       int *NOut, double **XOut, double **YOut)
{
  int    totPts, ii, kk;
  double HX;
  psVector VecXX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[totPts];
  VecXX.setLength(totPts*nInputs_);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (kk = 0; kk < nInputs_; kk++) 
      VecXX[ii*nInputs_+kk] = settings[kk]; 
    VecXX[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
    (*XOut)[ii] = HX * ii + lowerBounds_[ind1];
  }
   
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, VecXX.getDVector(), *YOut);
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int GP3::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                  double *settings,int *NOut,double **XOut,double **YOut)
{
  int totPts, ii, jj, kk, index;
  psVector VecXX, VecHX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_;
  VecHX.setLength(2);
  VecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  VecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  VecXX.setLength(totPts*nInputs_);
  (*XOut) = new double[2*totPts];
  checkAllocate(*XOut, "XOut in GP3::gen2DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      index = ii * nPtsPerDim_ + jj;
      for (kk = 0; kk < nInputs_; kk++) 
        VecXX[index*nInputs_+kk] = settings[kk]; 
      VecXX[index*nInputs_+ind1]  = VecHX[0] * ii + lowerBounds_[ind1];
      VecXX[index*nInputs_+ind2]  = VecHX[1] * jj + lowerBounds_[ind2];
      (*XOut)[index*2]   = VecHX[0] * ii + lowerBounds_[ind1];
      (*XOut)[index*2+1] = VecHX[1] * jj + lowerBounds_[ind2];
    }
  }
    
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, VecXX.getDVector(), (*YOut));
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int GP3::gen3DGridData(double *XIn, double *YIn,int ind1,int ind2,int ind3,
                  double *settings,int *NOut,double **XOut,double **YOut)
{
  int totPts, ii, jj, kk, ll, index;
  psVector VecXX, VecHX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  VecHX.setLength(3);
  VecHX[0] = (upperBounds_[ind1]-lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  VecHX[1] = (upperBounds_[ind2]-lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  VecHX[2] = (upperBounds_[ind3]-lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  VecXX.setLength(totPts*nInputs_);
  (*XOut) = new double[3*totPts];
  checkAllocate(*XOut, "XOut in GP3::gen3DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      for (ll = 0; ll < nPtsPerDim_; ll++) 
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        for (kk = 0; kk < nInputs_; kk++) 
          VecXX[index*nInputs_+kk] = settings[kk]; 
        VecXX[index*nInputs_+ind1]  = VecHX[0] * ii + lowerBounds_[ind1];
        VecXX[index*nInputs_+ind2]  = VecHX[1] * jj + lowerBounds_[ind2];
        VecXX[index*nInputs_+ind3]  = VecHX[2] * ll + lowerBounds_[ind3];
        (*XOut)[index*3]   = VecHX[0] * ii + lowerBounds_[ind1];
        (*XOut)[index*3+1] = VecHX[1] * jj + lowerBounds_[ind2];
        (*XOut)[index*3+2] = VecHX[2] * ll + lowerBounds_[ind3];
      }
    }
  }
    
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, VecXX.getDVector(), (*YOut));
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int GP3::gen4DGridData(double *XIn,double *YIn,int ind1,int ind2, int ind3,
                       int ind4, double *settings, int *NOut,double **XOut, 
                       double **YOut)
{
  int totPts, ii, jj, kk, ll, mm, index;
  psVector VecXX, VecHX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
 
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  VecHX.setLength(4);
  VecHX[0] = (upperBounds_[ind1]-lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  VecHX[1] = (upperBounds_[ind2]-lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  VecHX[2] = (upperBounds_[ind3]-lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  VecHX[3] = (upperBounds_[ind4]-lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  VecXX.setLength(totPts*nInputs_);
  (*XOut) = new double[4*totPts];
  checkAllocate(*XOut, "XOut in GP3::gen4DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      for (ll = 0; ll < nPtsPerDim_; ll++) 
      {
        for (mm = 0; mm < nPtsPerDim_; mm++) 
        {
          index = ii*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                  jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
          for (kk = 0; kk < nInputs_; kk++) 
            VecXX[index*nInputs_+kk] = settings[kk]; 
          VecXX[index*nInputs_+ind1]  = VecHX[0] * ii + lowerBounds_[ind1];
          VecXX[index*nInputs_+ind2]  = VecHX[1] * jj + lowerBounds_[ind2];
          VecXX[index*nInputs_+ind3]  = VecHX[2] * ll + lowerBounds_[ind3];
          VecXX[index*nInputs_+ind4]  = VecHX[3] * mm + lowerBounds_[ind4];
          (*XOut)[index*4]   = VecHX[0] * ii + lowerBounds_[ind1];
          (*XOut)[index*4+1] = VecHX[1] * jj + lowerBounds_[ind2];
          (*XOut)[index*4+2] = VecHX[2] * ll + lowerBounds_[ind3];
          (*XOut)[index*4+3] = VecHX[3] * mm + lowerBounds_[ind4];
        }
      }
    }
  }
    
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("GP3 interpolation begins....\n");
  evaluatePoint(totPts, VecXX.getDVector(), (*YOut));
  if (outputLevel_ >= 1) printf("GP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double GP3::evaluatePoint(double *X)
{
  int    ii, iOne=1;
  double Y=0.0;
  psVector vecXT;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  vecXT.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    vecXT[ii] = (X[ii] - XMeans_[ii]) / XStds_[ii];
  interpolate(iOne, vecXT.getDVector(), &Y, NULL);
  Y = Y * YStd_ + YMean_;
  return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double GP3::evaluatePoint(int npts, double *X, double *Y)
{
  int      ii, jj;
  psVector VecXX;

  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  VecXX.setLength(npts*nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    for (jj = 0; jj < npts; jj++)
      VecXX[jj*nInputs_+ii] = (X[jj*nInputs_+ii]-XMeans_[ii])/XStds_[ii];
  interpolate(npts, VecXX.getDVector(), Y, NULL);
  for (jj = 0; jj < npts; jj++) Y[jj] = Y[jj] * YStd_ + YMean_;
  return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double GP3::evaluatePointFuzzy(double *X, double &std)
{
  int      ii, iOne=1;
  double   Y=0.0;
  psVector VecXX;

  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  VecXX.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    VecXX[ii] = (X[ii] - XMeans_[ii]) / XStds_[ii];
  interpolate(iOne, VecXX.getDVector(), &Y, &std);
  Y = Y * YStd_ + YMean_;
  if (std < 0) printf("GP3 ERROR: variance (%e) < 0\n",std);
  else         std = sqrt(std) * YStd_;
  return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double GP3::evaluatePointFuzzy(int npts,double *X, double *Y,double *Ystds)
{
  int      ii, jj;
  psVector VecXX;

  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  VecXX.setLength(npts*nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    for (jj = 0; jj < npts; jj++)
      VecXX[jj*nInputs_+ii] = (X[jj*nInputs_+ii]-XMeans_[ii])/XStds_[ii];
  interpolate(npts, VecXX.getDVector(), Y, Ystds);
  for (int ii = 0; ii < npts; ii++)
  {
    Y[ii] = Y[ii] * YStd_ + YMean_;
    if (Ystds[ii] < 0) printf("GP3 ERROR: variance (%e) < 0\n", Ystds[ii]);
    else               Ystds[ii] = sqrt(Ystds[ii]) * YStd_;
  }
  return 0.0;
}

// ************************************************************************
// train to get lengthScales and sigF2 and sigN2
// ------------------------------------------------------------------------
int GP3::train()
{
  int     ii, jj, kk, its, count, status, newnInps, nhypers;
  int     nSamp=1000, fail, nHist, maxHist=10;
  double  minVal=1e35, dmax, dmin, dtemp, dmean, dstd, FValue;
  double  sqrt2pi=sqrt(2*3.1415928);
  Sampling *sampler;
  psVector YT, VecLB, VecUB, VecPVals, VecGrads, VecW, VecXS, VecYS;
  psVector VecHist;

  nhypers = nInputs_ + 4;
  if (optLinTerm_) nhypers += nInputs_;

  hyperparameters_.setLength(nhypers);
  for (ii = 0; ii < nInputs_; ii++)
  {
    dmax = -PSUADE_UNDEFINED;
    dmin =  PSUADE_UNDEFINED;
    for (kk = 0; kk < nSamples_; kk++)
    {
      if (XDataN_[kk*nInputs_+ii] > dmax) dmax = XDataN_[kk*nInputs_+ii];
      if (XDataN_[kk*nInputs_+ii] < dmin) 
        dmin = XDataN_[kk*nInputs_+ii];
    }
    hyperparameters_[ii] = 0.25 * (dmax - dmin);
    if (hyperparameters_[ii] == 0)
    {
      printf("GP3 ERROR: Input %d is a constant.\n", ii+1);
      printf("           Prune this input first (use idelete).\n");
      exit(1);
    }
    hyperparameters_[ii] = log(hyperparameters_[ii]);
  }
  dtemp = 0.0;
  for (kk = 0; kk < nSamples_; kk++) 
    if (PABS(YData_[kk]) > dtemp) dtemp = PABS(YData_[kk]);
  hyperparameters_[nInputs_] = dtemp;
    
  hyperparameters_[nInputs_+1] = log(1e-12);
  hyperparameters_[nInputs_+2] = 0;
  hyperparameters_[nInputs_+3] = log(1e-10);
  count = nInputs_ + 4;
  if (optLinTerm_) 
  {
    for (ii = 0; ii < nInputs_; ii++) 
      hyperparameters_[count++] = 0.0;
  }

  int  optimizeFlag=1;
  char pString[1000], winput[1000], *inStr;
  if (psConfig_ != NULL)
  {
    inStr = psConfig_->getParameter("GP3_optimize2");
    if (inStr != NULL) optimizeFlag = 2;
  }
  if (optimizeFlag != 2 && psRSExpertMode_ == 1 && isScreenDumpModeOn() == 1)
  {
    sprintf(pString,"Use user-provided hyperparameters? (y or n) ");
    getString(pString, winput);
    if (winput[0] == 'y')
    {
      for (ii = 0; ii < hyperparameters_.length(); ii++)
      {
        sprintf(pString,"Enter hyperparameter %d : ",ii+1);
        hyperparameters_[ii] = getDouble(pString);
      }
      optimizeFlag = 0;
    }
    else
    {
      sprintf(pString,"Select optimizer (1 or 2, default = 1): ");
      optimizeFlag = getInt(1, 2, pString);
      if (psConfig_ != NULL && optimizeFlag == 2)
        psConfig_->putParameter("GP3_optimize2");
    }
  }

  status = computeDistances();
  if (status != 0)
  {
    printf("GP3 INFO: since there are repeated sample points.\n");
    printf("    GP will continue without optimizing the length\n");
    printf("    scales but instead set them to 1.\n");
    printf("    This may not give a good quality response surface.\n");
    printf("    So if you want to prune the sample and do it again,\n");
    printf("    terminate now. I am giving you 30 seconds to decide.\n");
    optimizeFlag = 0;
    for (ii = 0; ii < nInputs_; ii++) hyperparameters_[ii] = 1.0;
    for (ii = nInputs_; ii < hyperparameters_.length(); ii++)
      hyperparameters_[ii] = 0.0;
#ifdef WINDOWS
    Sleep(20000);
#else
    sleep(20);
#endif
    printf("    10 more seconds.\n");
#ifdef WINDOWS
    Sleep(10000);
#else
    sleep(10);
#endif
  }

  XLows_.setLength(hyperparameters_.length());
  XHighs_.setLength(hyperparameters_.length());
  for (ii = 0; ii < hyperparameters_.length(); ii++) 
  {
    XLows_[ii]  = +PSUADE_UNDEFINED;
    XHighs_[ii] = -PSUADE_UNDEFINED;
  }

  psIVector VecIS;
  if (optimizeFlag == 1)
  {
    nhypers = hyperparameters_.length();
    VecPVals.setLength(nhypers);
    VecLB.setLength(nhypers);
    VecUB.setLength(nhypers);
    for (ii = 0; ii < nInputs_; ii++) 
    {
      VecLB[ii] = -4.0;
      VecUB[ii] =  4.0;
    }
    VecLB[nInputs_] = -3.0;
    VecUB[nInputs_] = 10.0;
    VecLB[nInputs_+1] = -28;
    VecUB[nInputs_+1] = - 3;
    VecLB[nInputs_+2] = -0.5;
    VecUB[nInputs_+2] =  0.5;
    VecLB[nInputs_+3] = -24;
    VecUB[nInputs_+3] = - 2;
    count = nInputs_ + 4;
    if (optLinTerm_) 
    {
      for (ii = 0; ii < nInputs_; ii++) 
      {
        VecLB[count] = - 5*hyperparameters_[nInputs_]/hyperparameters_[ii];
        VecUB[count++] = + 5*hyperparameters_[nInputs_]/hyperparameters_[ii];
      }
    }
    if (outputLevel_ > 0) 
    {
      for (ii = 0; ii < nhypers; ii++)
        printf("   Hyperparameter bounds %3d = %12.4e %12.4e\n",ii+1,
               VecLB[ii],VecUB[ii]);
    }
    nSamp = 2;
    if (nhypers > 51)
      sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
    else
      sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    sampler->setInputBounds(nhypers,VecLB.getDVector(),VecUB.getDVector());
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nSamp, 1, 0);
    sampler->initialize(0);
    VecIS.setLength(nSamp);
    VecXS.setLength(nSamp*nhypers);
    VecYS.setLength(nSamp);
    sampler->getSamples(nSamp,nhypers,1,VecXS.getDVector(),
                        VecYS.getDVector(), VecIS.getIVector());
    delete sampler;
    for (ii = 0; ii < nhypers; ii++) VecXS[ii] = hyperparameters_[ii];
 
    integer nInps, iprint=0, itask, *task=&itask, lsave[4], isave[44];
    integer *iwork, nCorr=5, *nbds, csave[60], nLBFGS=nSamp;
    double  factr, pgtol, dsave[29];

    nInps = nhypers;
    VecGrads.setLength(nInps);
    nbds = new integer[nInps];
    for (ii = 0; ii < nInps; ii++) nbds[ii] = 2;
    factr = 1e7;
    pgtol = 1e-8;
    kk = (2 * nCorr + 5) * nInps + 11 * nCorr * nCorr + 8 * nCorr;
    VecW.setLength(kk);
    iwork = new integer[3*nInps];
    its   = 0;
    VecHist.setLength(maxHist);

    for (ii = 0; ii < nLBFGS; ii++)
    {
      if (outputLevel_ > 0 && nLBFGS > 1 && ii > 0) 
        printf("GP3 LBFGS sample %d (%ld)\n",ii+1, nLBFGS);
      for (jj = 0; jj < nInps; jj++) VecPVals[jj] = VecXS[ii*nInps+jj];
      if (outputLevel_ > 0 && ii > 0) 
      {
        printf("   GP3 initial hyperparameter values:\n");
        for (ii = 0; ii < hyperparameters_.length(); ii++) 
          printf("   Hyperparameter %2d = %e\n",ii+1,
                 VecPVals[ii]);
      }
      nHist = 0;
      fail = 0;
      its = 0;
      *task = (integer) START;
      while (1 && its < 1000)
      {
        its++;
        setulb(&nInps,&nCorr,VecPVals.getDVector(),VecLB.getDVector(),
               VecUB.getDVector(),nbds,&FValue,VecGrads.getDVector(),
               &factr,&pgtol,VecW.getDVector(),iwork,task,&iprint,
               csave,lsave,isave,dsave);
        if (IS_FG(*task))
        {
          //if (outputLevel_ > 2) 
          //  for (kk = 0; kk < nInps; kk++)
          //    printf("     Hyperparameter %5d = %e\n",kk+1,VecPVals[kk]);

#if 1
          FValue = computeGradients(VecPVals, VecGrads, status);
          if (FValue == 1e35) 
          {
            fail = 1;
            break;
          }
          dtemp = 0.0;
          for (kk = 0; kk < nInps; kk++)
            dtemp += pow(VecGrads[kk], 2.0);
          dtemp = sqrt(dtemp);
#else
          status = 0;
          FValue = computeLikelihood(XDataN_, YData_, VecPVals);
          if (FValue == 1e35) 
          {
            fail = 1;
            break;
          }
          dtemp = 0.0;
          for (kk = 0; kk < nInps; kk++)
          {
            VecPVals[kk] += 1.0e-11;
            VecGrads[kk] = computeLikelihood(XDataN_,YData_,VecPVals);
            VecGrads[kk] = 1e11 * (VecGrads[kk] - FValue);
            VecPVals[kk] -= 1.0e-11;
            dtemp += pow(VecGrads[kk], 2.0);
          }
          dtemp = sqrt(dtemp);
#endif
          if (outputLevel_ > 1) 
            printf("   LBFGS Current FValue = %e (its=%d, |grad|=%e)\n",
                   FValue,its,dtemp);
          if (outputLevel_ > 2) 
            for (kk = 0; kk < nInps; kk++)
              printf("         Gradient %5d = %e\n",kk+1,VecGrads[kk]);
          if (FValue < minVal && status == 0)
          {
            for (jj = 0; jj < nhypers; jj++) 
              hyperparameters_[jj] = VecPVals[jj];
            minVal = FValue;
          }
          if (nHist < maxHist && FValue != 1e35) VecHist[nHist++] = FValue;
          else if (FValue != 1e35)
          {
            for (kk = 0; kk < maxHist-1; kk++) VecHist[kk] = VecHist[kk+1];
            VecHist[maxHist-1] = FValue;
          }
          if (nHist >= maxHist && FValue != 1e35)
          {
            dmean = 0.0;
            for (kk = 0; kk < maxHist; kk++) dmean += VecHist[kk];
            dmean /= (double) maxHist;
            dstd = 0.0;
            for (kk = 0; kk < maxHist; kk++) 
              dstd += pow(VecHist[kk]-dmean,2.0);
            dstd = sqrt(dstd/(maxHist-1.0));
            if (dmean == 0) dmean = 1;
            if (PABS(dstd/dmean) < 1e-3)
            {
              *task = STOP_ITER;
              if (outputLevel_ > 1) 
              {
                printf("INFO: PSUADE issues a stop (converged))\n");
                printf("INFO: current iteration   = %d\n", its);
                printf("INFO: current best FValue = %e (%e)\n",minVal,
                       PABS(dstd/dmean));
                printf("INFO: current grad norm   = %e\n",dtemp);
              }
            }
            else
            {
              if (outputLevel_ > 1) 
                printf("INFO: convergence check: %e < 1e-3?\n",
                       PABS(dstd/dmean));
            }
          }
          if (isave[33] >= 300) 
          {
            *task = STOP_ITER;
            if (outputLevel_ > 1) 
            {
              printf("INFO: PSUADE issues a stop (> 300 iterations)\n");
              printf("INFO: current best FValue = %e\n", minVal);
            }
          }
        }
        else if (*task == NEW_X)
        {
          if (isave[33] >= 300) 
          {
            *task = STOP_ITER;
            if (outputLevel_ > 1) 
              printf("INFO: PSUADE issues a stop (> 300 iterations)\n");
          }
        }
        else
        {
          if (outputLevel_ > 1) 
          {
            printf("INFO: LBFGS issues a stop (its = %d)\n",its);
            printf("INFO: current best FValue = %e\n", minVal);
          }
          break;
        }
      }
      if (outputLevel_ > 1) 
      {
        printf("   GP3 final hyperparameter values:\n");
        for (jj = 0; jj < hyperparameters_.length(); jj++) 
          printf("   Hyperparameter %2d = %e\n",jj+1,
                 hyperparameters_[jj]);
        printf("   Sample %5d: Fvalue (final) = %e\n", ii+1, minVal);
      }
      if (ii == (nLBFGS-1) && (fail == 1)) nLBFGS++;
    }
    delete [] iwork;
    delete [] nbds;

    if (outputLevel_ > 0) 
    {
      printf("GP3 very final hyperparameter values:\n");
      for (ii = 0; ii < hyperparameters_.length(); ii++) 
        printf("   Hyperparameter %2d = %24.16e\n",ii+1,
               hyperparameters_[ii]);
    }
  }
  else if (optimizeFlag == 2)
  {
    printf("GP3: you are choosing optimization option 2.\n");
    printf("     This isn't working as good as option 1.\n");
    nhypers = hyperparameters_.length();
    VecPVals.setLength(nhypers);
    VecLB.setLength(nhypers);
    VecUB.setLength(nhypers);
    for (ii = 0; ii < nInputs_; ii++) 
    {
      VecLB[ii] = -4.0;
      VecUB[ii] =  4.0;
    }
    VecLB[nInputs_] = -3.0;
    VecUB[nInputs_] = 10.0;
    VecLB[nInputs_+1] = -28;
    VecUB[nInputs_+1] = - 3;
    VecLB[nInputs_+2] = -0.5;
    VecUB[nInputs_+2] =  0.5;
    VecLB[nInputs_+3] = -24;
    VecUB[nInputs_+3] = - 2;
    count = nInputs_ + 4;
    if (optLinTerm_) 
    {
      for (ii = 0; ii < nInputs_; ii++) 
      {
        VecLB[count] = - 5*hyperparameters_[nInputs_]/hyperparameters_[ii];
        VecUB[count++] = + 5*hyperparameters_[nInputs_]/hyperparameters_[ii];
      }
    }
    if (outputLevel_ > 0) 
    {
      for (ii = 0; ii < nhypers; ii++)
        printf("   Hyperparameter bounds %3d = %12.4e %12.4e\n",ii+1,
               VecLB[ii],VecUB[ii]);
    }
    nSamp = 5;
    if (nhypers > 51)
         sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
    else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    sampler->setInputBounds(nhypers,VecLB.getDVector(),VecUB.getDVector());
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nSamp, 1, 0);
    sampler->initialize(0);
    VecIS.setLength(nSamp);
    VecXS.setLength(nSamp*nhypers);
    VecYS.setLength(nSamp);
    sampler->getSamples(nSamp,nhypers,1,VecXS.getDVector(),
                        VecYS.getDVector(), VecIS.getIVector());
    delete sampler;
    for (ii = 0; ii < nhypers; ii++) VecXS[ii] = hyperparameters_[ii];
 
    int      nUOA=nSamp, nPts, pLevel=5555, maxfun=1000;
    double   tol=1e-5, rhobeg, rhoend;

    nPts = (nhypers + 1) * (nhypers + 2) / 2;
    jj = (nPts+13) * (nPts+nhypers) + 3*nhypers*(nhypers+3)/2;
    kk = (nPts+5)*(nPts+nhypers)+3*nInputs_*(nhypers+5)/2+1;
    if (jj > kk) VecW.setLength(jj);
    else         VecW.setLength(jj);
    rhobeg = VecUB[0] - VecLB[0];
    for (ii = 1; ii < nhypers; ii++)
    {
      dtemp = VecUB[ii] - VecLB[ii];
      if (dtemp < rhobeg) rhobeg = dtemp;
    }
    rhobeg *= 0.5;
    rhoend = rhobeg * tol;

    its = 0;
    GP3_OptY = 1e35;
    if (GP3Obj != NULL) delete GP3Obj;
    GP3Obj = new GP3(nInputs_, nSamples_);
    GP3Obj->XDataN_ = XDataN_;
    GP3Obj->YData_ = YData_;
    GP3Obj->XDistances_ = XDistances_;

    if (outputLevel_ > 0) printf("GP3 NEWUOA optimization (%d):\n",nUOA);
    for (ii = 0; ii < nUOA; ii++)
    {
      if (outputLevel_ > 0 && nUOA > 1 && ii > 0) 
        printf("GP3 WLS sample %d (%d)\n",ii+1, nUOA);
      for (jj = 0; jj < nhypers; jj++) VecPVals[jj] = VecXS[ii*nhypers+jj];
      if (outputLevel_ > 0) 
      {
        printf("   GP3 initial hyperparameter values:\n");
        for (jj = 0; jj < VecPVals.length(); jj++) 
          printf("   Hyperparameter %2d = %e\n",jj+1,VecPVals[jj]);
      }
      newuoa_(&nhypers,&nPts,VecPVals.getDVector(),&rhobeg,&rhoend,&pLevel,
              &maxfun, VecW.getDVector());
      hyperparameters_ = GP3_OptX;
      //THIS NEEDS DEBUGGED
      if (outputLevel_ > 0) 
      {
        printf("   GP3 final hyperparameter values:\n");
        for (jj = 0; jj < VecPVals.length(); jj++) 
          printf("   Hyperparameter %2d = %e\n",jj+1,
                 VecPVals[jj]);
        printf("   Best objective function so far = %e\n",GP3_OptY);
      }
    }
    delete GP3Obj;
    GP3Obj = NULL;
    if (outputLevel_ > 0) 
    {
      printf("GP3 very final hyperparameter values:\n");
      for (ii = 0; ii < hyperparameters_.length(); ii++) 
        printf("   Hyperparameter %2d = %24.16e\n",ii+1,
               hyperparameters_[ii]);
      printf("   Best objective function = %e\n",GP3_OptY);
    }
  }

  constructCMatrix(CMatrix_, hyperparameters_);
  CMatrix_.LUDecompose();
  CInvY_.setLength(nSamples_);
  YT.setLength(nSamples_);
  for (jj = 0; jj < nSamples_; jj++) 
    YT[jj] = YData_[jj] - hyperparameters_[nInputs_+2];
  CMatrix_.LUSolve(YT, CInvY_);

  //Cannot do because inverse seems to have problems 3/2017
  //double likelihood = 0.0;
  //for (jj = 0; jj < nSamples_; jj++) likelihood += CInvY_[jj] * YT[jj];
  //likelihood *= 0.5;
  //for (jj = 0; jj < nSamples_; jj++) 
  //{
  //  dtemp = CMatrix.getEntry(jj, jj);
  //  likelihood += log(dtemp*sqrt2pi+1e-20);
  //}
  //if (outputLevel_ > 0) 
  //  printf("GP3 Best likelihood = %e\n", likelihood);
  if (psGMMode_ == 1)
  { 
    psMatrix CMatrix;
    constructCMatrix(CMatrix, hyperparameters_);
    FILE *fp = fopen("Cmat.m", "w");
    fprintf(fp, "A = [\n");
    for (ii = 0; ii < nSamples_; ii++) 
    {
      for (jj = 0; jj < nSamples_; jj++) 
        fprintf(fp, "%16.8e ", CMatrix.getEntry(ii,jj));
      fprintf(fp, "\n");
    }
    fprintf(fp,"];\n");
    fclose(fp);
    printf("*** Final Cmatrix given in Cmat.m\n");
  }
  if (outputLevel_ > 1 && optimizeFlag == 1)
  {
    printf("Optimization summary: \n");
    for (ii = 0; ii < nInputs_+4; ii++) 
      printf("Input %2d search bounds = %12.4e %12.4e\n",ii+1,
             XLows_[ii], XHighs_[ii]);
  }
  return 0;
}

// ************************************************************************
// compute log likelihood value (called directly from external)
// ------------------------------------------------------------------------
double GP3::computeLikelihood(psVector VecXD,psVector VecYN,psVector VecP)
{
  int    ii, jj, kk, status, count;
  double likelihood, sqrt2pi=sqrt(2*3.1415928), ddata, fudge, *arrayC;
  double dist;
  psMatrix MatC, MatCI;
  psVector VecY, VecCIY;

  if (outputLevel_ > 3) 
  {
    printf("computeLikelihood:\n");
    for (jj = 0; jj < VecP.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, VecP[jj]);
  }
  MatC.setDim(nSamples_, nSamples_);
  arrayC = MatC.getMatrix1D();
  count = 0;
  fudge = sqrt(nSamples_);
  for (jj = 0; jj < nSamples_; jj++)
  {
    arrayC[jj*nSamples_+jj] = exp(VecP[nInputs_]) +
                fudge*exp(VecP[nInputs_+1]) + exp(VecP[nInputs_+3]);
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
        dist += pow(VecXD[count*nInputs_+ii],2.0)/
                exp(2.0*VecP[ii]);
      dist *= 0.5;
      dist = exp(VecP[nInputs_]) * exp(-dist);
      if (dist < 1.0e-50) dist = 0;
      arrayC[jj*nSamples_+kk] = dist + exp(VecP[nInputs_+1]);
      arrayC[kk*nSamples_+jj] = dist + exp(VecP[nInputs_+1]);
      count++;
    }
  }

  status = MatC.computeInverse(MatCI);
  if (status != 0)
  {
    printf("GP3 ERROR (1): failed in matrix inversion (%d).\n",status);
    for (jj = 0; jj < VecP.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, VecP[jj]);
    return 1e35;
  }
  VecY.setLength(nSamples_);
  for (jj = 0; jj < nSamples_; jj++) 
    VecY[jj] = VecYN[jj] - VecP[nInputs_+2];
  VecCIY.setLength(nSamples_);
  MatCI.matvec(VecY, VecCIY, 0);

  int errCount=0;
  psMatrix MatEVecs;
  psVector VecEVals;
  MatC.eigenSolve(MatEVecs, VecEVals, 1);
  likelihood = 0.0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    ddata = VecEVals[jj];
    if (ddata <= 0) errCount++;
    if (ddata != 0)
      likelihood += log(sqrt(PABS(ddata))*sqrt2pi+1e-50);
  }
  if (errCount > 0)
  {
    if (outputLevel_ > 0)
      printf("GP3 WARNING: Covariance Matrix not positive definite\n");
  }
  for (jj = 0; jj < nSamples_; jj++)
    likelihood += 0.5 * (VecCIY[jj] * VecY[jj]);
  if (likelihood < 0) likelihood = - likelihood;
  if (outputLevel_ > 3) 
    printf("   computeLikelihood: Likelihood = %e\n", likelihood);
  return likelihood;
}

// ************************************************************************
// compute gradients of log likelihood 
// ------------------------------------------------------------------------
double GP3::computeGradients(psVector VecHypers, psVector &VecGrads, 
                             int &retstat)
{
  int    jj, kk, ii, ll, status, count;
  double ddata, *TMat,likelihood,sqrt2pi=sqrt(2*3.1415928),dist;
  double fudge;
  psMatrix CMatrix, CInverse, CPrime, CT;
  psVector YVec, CInvY, YT;

  if (outputLevel_ > 2) 
  {
    printf("computeGradients:\n");
    for (jj = 0; jj < VecHypers.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, VecHypers[jj]);
  }
  constructCMatrix(CMatrix, VecHypers);

  if (XHighs_.length() == VecHypers.length())
  {
    for (ii = 0; ii < VecHypers.length(); ii++)
    {
      if (VecHypers[ii] > XHighs_[ii]) XHighs_[ii] = VecHypers[ii]; 
      if (VecHypers[ii] < XLows_[ii])  XLows_[ii]  = VecHypers[ii]; 
    }
  }

  VecGrads.setLength(VecHypers.length());
  status = CMatrix.computeInverse(CInverse);
  if (status != 0)
  {
    printf("GP3 ERROR (1): failed in matrix inversion (%d).\n",status);
    for (jj = 0; jj < VecHypers.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, VecHypers[jj]);
    for (jj = 0; jj < VecHypers.length(); jj++) VecGrads[jj] = 100;
    return 1e35;
  }
  YVec.setLength(nSamples_);
  for (jj = 0; jj < nSamples_; jj++) 
    YVec[jj] = YData_[jj] - VecHypers[nInputs_+2];
  CInvY.setLength(nSamples_);
  CInverse.matvec(YVec, CInvY, 0);

  YT.setLength(nSamples_);
  CPrime.setDim(nSamples_, nSamples_);
  TMat = CPrime.getMatrix1D();
  for (ll = 0; ll < nInputs_; ll++)
  {
    count = 0;
    for (jj = 0; jj < nSamples_; jj++)
    {
      TMat[jj*nSamples_+jj] = 0;
      for (kk = jj+1; kk < nSamples_; kk++)
      {
        dist = 0.0;
        for (ii = 0; ii < nInputs_; ii++)
          dist += pow(XDistances_[count*nInputs_+ii],2.0)/
                  exp(2.0*VecHypers[ii]);
        dist *= 0.5;
        dist = exp(VecHypers[nInputs_]) * exp(-dist);
        ddata = pow(XDistances_[count*nInputs_+ll]/
                exp(2.0*VecHypers[ll]),2.0);
        ddata = ddata * exp(2.0*VecHypers[ll]);
        dist *= ddata;
        TMat[jj*nSamples_+kk] = dist;
        TMat[kk*nSamples_+jj] = dist;
        count++;
      }
    }
    CPrime.matvec(CInvY, YT, 0);
    ddata = 0.0;
    for (jj = 0; jj < nSamples_; jj++) ddata += YT[jj] * CInvY[jj];
    VecGrads[ll] = - ddata;
    ddata = 0.0;
    for (ii = 0; ii < nSamples_; ii++)
      for (jj = 0; jj < nSamples_; jj++) 
        ddata += CPrime.getEntry(ii, jj) * CInverse.getEntry(ii, jj);
    VecGrads[ll] += ddata;
    VecGrads[ll] *= 0.5;
  }

  count = 0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    TMat[jj*nSamples_+jj] = 1;
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
        dist += pow(XDistances_[count*nInputs_+ii],2.0)/
                exp(2.0*VecHypers[ii]);
      dist = exp(-0.5 * dist);
      if (dist < 1.0e-50) dist = 0;
      TMat[jj*nSamples_+kk] = dist;
      TMat[kk*nSamples_+jj] = dist;
      count++;
    }
  }
  CPrime.matvec(CInvY, YT, 0);
  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += YT[jj] * CInvY[jj];
  VecGrads[nInputs_] = - ddata;
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++) 
    for (jj = 0; jj < nSamples_; jj++) 
      ddata += CInverse.getEntry(ii, jj) * CPrime.getEntry(ii, jj);
  VecGrads[nInputs_] += ddata;
  VecGrads[nInputs_] *= 0.5 * exp(VecHypers[nInputs_]);

  fudge = sqrt(1.0 * nSamples_);
  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += CInvY[jj];
  VecGrads[nInputs_+1] = - ddata * ddata - (fudge - 1.0) * ddata;
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++)
    for (jj = 0; jj < nSamples_; jj++) 
      ddata += CInverse.getEntry(ii, jj);
  VecGrads[nInputs_+1] += ddata;
  ddata = exp(VecHypers[nInputs_+1]);
  VecGrads[nInputs_+1] *= 0.5 * ddata;

  //*/ compute: - 1' * C' * C^{-1} (Y - mu)
  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += CInvY[jj];
  VecGrads[nInputs_+2] = - ddata;

  ddata = 0.0;
  for (jj = 0; jj < nSamples_; jj++) ddata += CInvY[jj] * CInvY[jj];
  VecGrads[nInputs_+3] = - ddata;
  ddata = 0.0;
  for (ii = 0; ii < nSamples_; ii++) 
    ddata += CInverse.getEntry(ii, ii);
  VecGrads[nInputs_+3] += ddata;
  VecGrads[nInputs_+3] *= 0.5 * exp(VecHypers[nInputs_+3]);

  int errCount=0;
  CMatrix.eigenSolve(CT, YT, 1);
  likelihood = 0.0;
  for (jj = 0; jj < nSamples_; jj++) 
  {
    ddata = YT[jj];
    if (ddata <= 0) errCount++;
    if (ddata != 0)
      likelihood += log(sqrt(PABS(ddata))*sqrt2pi+1e-50);
  }
  retstat = 0;
  if (errCount > 0)
  {
    retstat = 1;
    if (outputLevel_ > 0)
    {
      printf("GP3 WARNING: CMatrix non-positive definite (%d)\n",
             errCount);
      if (outputLevel_ > 2)
      {
        for (jj = 0; jj < nInputs_; jj++)
          printf("   hyperparameter %2d = %e\n",jj+1,VecHypers[jj]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+1,VecHypers[nInputs_]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+2,VecHypers[nInputs_+1]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+3,VecHypers[nInputs_+2]);
        printf("   hyperparameter %2d = %e\n",
               nInputs_+4,VecHypers[nInputs_+3]);
      }
    }
    if (psGMMode_ == 1)
    { 
      printf("INFO: non-PD covariance matrix, terminate in GM mode.\n");
      printf("      diagnostics information is stored in errmat.m\n");
      FILE *fp = fopen("errmat.m", "w");
      fprintf(fp, "n = %d;\n", nSamples_);
      fprintf(fp, "f = %e;\n", sqrt(1.0*nSamples_));
      fprintf(fp, "A = [\n");
      for (ii = 0; ii < nSamples_; ii++) 
      {
        for (jj = 0; jj < nSamples_; jj++) 
          fprintf(fp, "%16.8e ", CMatrix.getEntry(ii,jj));
        fprintf(fp, "\n");
      }
      fprintf(fp, "];\n");
      fprintf(fp, "L = [\n");
      for (ii = 0; ii < nSamples_*(nSamples_-1)/2*nInputs_; ii++) 
        fprintf(fp, "%24.16e\n", XDistances_[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "H = [\n");
      for (jj = 0; jj < nInputs_; jj++)
        fprintf(fp, "%e\n", exp(2.0*VecHypers[jj]));
      fprintf(fp, "];\n");
      fprintf(fp, "%e\n", exp(VecHypers[nInputs_]));
      fprintf(fp, "%e\n", exp(VecHypers[nInputs_+1]));
      fprintf(fp, "%e\n", VecHypers[nInputs_+2]);
      fprintf(fp, "%e\n", exp(VecHypers[nInputs_+3]));
      fprintf(fp, "count = 0;\n");
      fprintf(fp, "for ii = 1 : n\n");
      fprintf(fp, "  B(ii,ii) = H(4) + ll * H(5) + H(7);\n");
      fprintf(fp, "  for jj = ii+1 : n \n");
      fprintf(fp, "    dist = 0;\n");
      fprintf(fp, "    for kk = 1 : 3\n");
      fprintf(fp, "      dist = dist + L(count*3+kk)^2 / H(kk);\n");
      fprintf(fp, "    end;\n");
      fprintf(fp, "    ddata = H(4) * exp(-0.5*dist);\n");
      fprintf(fp, "    B(ii,jj) = ddata + H(5);\n");
      fprintf(fp, "    B(jj,ii) = ddata + H(5);\n");
      fprintf(fp, "    count = count + 1;\n");
      fprintf(fp, "  end;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "Bmin = min(eig(B))\n");
      fclose(fp);
      exit(1);
    }
  }
  for (jj = 0; jj < nSamples_; jj++) 
    likelihood += 0.5 * (CInvY[jj] * YVec[jj]);
  if (outputLevel_ > 2) 
    printf("   computeGradients: Likelihood = %e\n", likelihood);
  return likelihood;
}

// ************************************************************************
// compute pairwise distances 
// ------------------------------------------------------------------------
int GP3::computeDistances()
{
  int    ii, jj, kk, count, error=0;
  double *LDists, dist;

  LDists = new double[(nSamples_*(nSamples_-1)/2)*nInputs_];
  checkAllocate(LDists, "LDists in GP3::computeDistances");
  count = 0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        LDists[count*nInputs_+ii] = XDataN_[jj*nInputs_+ii] -
                                    XDataN_[kk*nInputs_+ii];
        if (LDists[count*nInputs_+ii] < 0)
           LDists[count*nInputs_+ii] = - LDists[count*nInputs_+ii];
        dist += pow(LDists[count*nInputs_+ii], 2.0);
      }
      if (dist == 0.0)
      {
        printf("GP3 ERROR: repeated sample points.\n");
        printf("           Prune repeated points and re-run.\n");
        printf("Sample %d : (with sample point %d)\n", kk+1, jj+1);
        for (ii = 0; ii < nInputs_; ii++)
          printf("   Input %d : %e\n",ii+1,
                 XDataN_[kk*nInputs_+ii]*XStds_[ii] +XMeans_[ii]);
        error = 1;
      }
      count++;
    }
  }
  XDistances_.load(count*nInputs_, LDists);
  return error;
}

// ************************************************************************
// interpolation 
// ------------------------------------------------------------------------
int GP3::interpolate(int npts, double *XX, double *Y, double *Ystds)
{
  int    ii, kk, nn, iOne=1, status, offset;
  double dist, expn, ddata;
  psVector YT, YT2;

  YT.setLength(nSamples_);
  YT2.setLength(nSamples_);
  Y[0] = 0;
  offset = nInputs_ + 4;
  if (optLinTerm_) offset += nInputs_;
  for (nn = 0; nn < npts; nn++)
  {
    for (kk = 0; kk < nSamples_; kk++)
    {
      expn = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        dist = XDataN_[kk*nInputs_+ii] - XX[nn*nInputs_+ii];
        expn += dist * dist / exp(2.0*hyperparameters_[ii]);
      }
      expn *= 0.5;
      YT[kk] = exp(hyperparameters_[nInputs_]) * exp(-expn);
    }
    Y[nn] = hyperparameters_[nInputs_+2];
    for (kk = 0; kk < nSamples_; kk++)
      Y[nn] += YT[kk] * CInvY_[kk];
    if (Ystds != NULL)
    {
      Ystds[nn] = 0.0;
      ddata = exp(hyperparameters_[nInputs_]) + 
              exp(hyperparameters_[nInputs_+1]) +
              exp(hyperparameters_[nInputs_+3]);
      CMatrix_.LUSolve(YT, YT2);
      for (kk = 0; kk < nSamples_; kk++)
        ddata -= YT[kk] * YT2[kk];
      if (ddata < 0) ddata = - ddata;
      Ystds[nn] = ddata;
    }
  }
  return 0;
}

// ************************************************************************
// construct C matrix
// ------------------------------------------------------------------------
void GP3::constructCMatrix(psMatrix &CMatrix, psVector VecP)
{
  int    ii, jj, kk, count;
  double dist, ddata, *TMat, fudge;

  if (outputLevel_ > 3) 
  {
    printf("Construct CMatrix:\n");
    for (jj = 0; jj < VecP.length(); jj++)
      printf("  hyperparameter %3d = %e\n", jj+1, VecP[jj]);
  }
  CMatrix.setDim(nSamples_, nSamples_);
  TMat = CMatrix.getMatrix1D();
  fudge = sqrt(nSamples_);
  count = 0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    TMat[jj*nSamples_+jj] = exp(VecP[nInputs_]) + 
                fudge*exp(VecP[nInputs_+1]) + exp(VecP[nInputs_+3]);
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
        dist += pow(XDistances_[count*nInputs_+ii],2.0)/
                exp(2.0*VecP[ii]);
      dist *= 0.5;
      dist = exp(VecP[nInputs_]) * exp(-dist);
      if (dist < 1.0e-50) dist = 0;
      TMat[jj*nSamples_+kk] = dist + exp(VecP[nInputs_+1]);
      TMat[kk*nSamples_+jj] = dist + exp(VecP[nInputs_+1]);
      count++;
    }
  }
  return;
}

// ************************************************************************
// generate C code
// ------------------------------------------------------------------------
void GP3::genCode()
{
  int ii, jj;
  FILE *fp = fopen("psuade_rs.info", "w");
  if (psRSCodeGen_ == 0) return;
  if (fp == NULL) return;
  fprintf(fp,"This file contains information to re-construct GP\n");
  fprintf(fp,"response surface offline. Follow the steps below:\n");
  fprintf(fp,"1. Search for the keywords 'SPLIT HERE' in this file.\n");
  fprintf(fp,"2. Store the lines below keywords into main.c\n");
  fprintf(fp,"3. Add the to-be-evaluated points inside main() in main.c\n");
  fprintf(fp,"   (or, embed all except main() in your own program).\n");
  fprintf(fp,"4. Compile main.c (cc -o main main.c -lm) and run\n");
  fprintf(fp,"   (Do not change the psuade_rs.info file.\n");
  fprintf(fp,"Note: if std dev is desired, uncomment the corresponding\n");
  fprintf(fp,"      section and compile with: cc main.c -llapack -lm\n");
  fprintf(fp,"\n");
  fprintf(fp,"PSUADE_BEGIN\n");
  fprintf(fp, "%d %d\n", nSamples_, nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    fprintf(fp, "%24.16e %24.16e\n", XMeans_[ii], XStds_[ii]);
  fprintf(fp, "%24.16e %24.16e\n", YMean_, YStd_);
  fprintf(fp, "%d\n", hyperparameters_.length());
  for (jj = 0; jj < hyperparameters_.length(); jj++)
    fprintf(fp, "%24.16e \n", hyperparameters_[jj]);
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, "%24.16e ", XDataN_[jj*nInputs_+ii]);
    fprintf(fp, "\n");
  }
  for (jj = 0; jj < nSamples_; jj++)
    fprintf(fp, "%24.16e\n", CInvY_[jj]);
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nSamples_; ii++)
      fprintf(fp, "%24.16e ", CMatrix_.getEntry(jj,ii));
    fprintf(fp, "\n");
  }
  for (jj = 0; jj < nSamples_; jj++)
    fprintf(fp, "%d\n", CMatrix_.pivots_[jj]);
  fprintf(fp,"==================== SPLIT HERE =====================\n");
  fprintf(fp,"/* ***********************************************/ \n");
  fprintf(fp,"/* GP interpolator from PSUADE.                  */ \n");
  fprintf(fp,"/* To estimate prediction uncertainty uncomment, */ \n");
  fprintf(fp,"/* dgetrs and the corresponding code segment.    */ \n");
  fprintf(fp,"/* ==============================================*/ \n");
  fprintf(fp,"#include <math.h>\n");
  fprintf(fp,"#include <stdlib.h>\n");
  fprintf(fp,"#include <stdio.h>\n");
  fprintf(fp,"/*void dgetrs_(char*,int*,int*,double*,int*,int*,\n");
  fprintf(fp,"               double*,int*,int*);*/\n");
  fprintf(fp,"int initialize();\n");
  fprintf(fp,"int finalize();\n");
  fprintf(fp,"int interpolate(int,double*,double*,double*);\n");
  fprintf(fp,"main(int argc, char **argv) {\n");
  fprintf(fp,"  int    i, iOne=1, nInps;\n");
  fprintf(fp,"  double X[%d], Y, Std;\n",nInputs_);
  fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
  fprintf(fp,"  if (argc < 2) {\n");
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
  fprintf(fp,"  initialize();\n");
  fprintf(fp,"  interpolate(iOne, X, &Y, &Std);\n");
  fprintf(fp,"  printf(\"Y = %%e (stdev = %%e)\\n\", Y, Std);\n");
  fprintf(fp,"  finalize();\n");
  fprintf(fp,"  if (argc == 3) {\n");
  fprintf(fp,"    fOut = fopen(argv[2], \"w\");\n");
  fprintf(fp,"    if (fOut == NULL) {\n");
  fprintf(fp,"      printf(\"ERROR: cannot open output file.\\n\");\n");
  fprintf(fp,"      exit(1);\n");
  fprintf(fp,"    }\n");
  fprintf(fp,"    fprintf(fOut,\" %%e\\n\", Y);\n");
  fprintf(fp,"  fclose(fOut);}\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"/* Regression interpolation function         */\n");
  fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
  fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
  fprintf(fp," * ... */\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int    nInps, nSamples, nParams, *pivs=NULL;\n");
  fprintf(fp,"double *XMeans=NULL,*XStds=NULL,YMean,YStd,*Thetas=NULL;\n");
  fprintf(fp,"double *CMat=NULL,*CInvY=NULL, *XN=NULL; \n");
  fprintf(fp,"int interpolate(int npts,double *X,double *Y,\n");
  fprintf(fp,"                double *YStds) {\n");
  fprintf(fp,"  int    ss, ii, jj, kk;\n");
  fprintf(fp,"  double expn, *xt, *yt, *zt, dist, ddata;\n");
  fprintf(fp,"  xt = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  yt = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  zt = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  for (ss = 0; ss < npts; ss++) {\n");
  fprintf(fp,"    for (ii = 0; ii < nInps; ii++)\n");
  fprintf(fp,"      xt[ii] = (X[ss*nInps+ii]-XMeans[ii])/XStds[ii];\n");
  fprintf(fp,"    for (kk = 0; kk < nSamples; kk++) {\n");
  fprintf(fp,"      expn = 0.0;\n");
  fprintf(fp,"      for (ii = 0; ii < nInps; ii++) {\n");
  fprintf(fp,"        dist = XN[kk*nInps+ii] - xt[ii];\n");
  fprintf(fp,"        expn += dist*dist/exp(2.0*Thetas[ii]);}\n");
  fprintf(fp,"      expn *= 0.5;\n");
  fprintf(fp,"      yt[kk] = exp(Thetas[nInps]) * exp(-expn);}\n");
  fprintf(fp,"    Y[ss] = Thetas[nInps+2];\n");
  fprintf(fp,"    for (kk = 0; kk < nSamples; kk++)\n");
  fprintf(fp,"      Y[ss] += yt[kk] * CInvY[kk];\n");
  fprintf(fp,"    Y[ss] = Y[ss] * YStd + YMean;\n");
  fprintf(fp,"    YStds[ss] = 0.0;\n");
  fprintf(fp,"    /* ==== if need to compute std dev. =====\n");
  fprintf(fp,"    ddata = exp(Thetas[nInps])+exp(Thetas[nInps+1])+");
  fprintf(fp,"exp(Thetas[nInps+3]);\n");
  fprintf(fp,"    LUSolve(yt, zt);\n");
  fprintf(fp,"    for (kk = 0; kk < nSamples; kk++)\n");
  fprintf(fp,"      ddata -= yt[kk] * zt[kk];\n");
  fprintf(fp,"    YStd);\n");
  fprintf(fp,"  fscanf(fp, \"%%d\", &nParams);\n");
  fprintf(fp,"  for (ii = 0; ii < nParams; ii++) {\n");
  fprintf(fp,"    fscanf(fp, \"%%lg\", &ddata);\n");
  fprintf(fp,"    Thetas[ii] = ddata;\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  XN = (double *) malloc(nSamples*nInps*sizeof(double));\n");
  fprintf(fp,"  for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"    for (ii = 0; ii < nInps; ii++) \n");
  fprintf(fp,"      fscanf(fp, \"%%lg\", &XN[jj*nInps+ii]);\n");
  fprintf(fp,"  CInvY = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"    fscanf(fp, \"%%lg\", &CInvY[ii]);\n");
  fprintf(fp,"  CMat=(double*) malloc(nSamples*nSamples*sizeof(double));\n");
  fprintf(fp,"  for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"    for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"      fscanf(fp, \"%%lg \", &CMat[jj+ii*nSamples]);\n");
  fprintf(fp,"  pivs = (int *) malloc(nSamples*sizeof(int));\n");
  fprintf(fp,"  for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"    fscanf(fp, \"%%d\", &pivs[ii]);\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int finalize() {\n");
  fprintf(fp,"  if (XMeans != NULL) free(XMeans);\n");
  fprintf(fp,"  if (XStds  != NULL) free(XStds);\n");
  fprintf(fp,"  if (Thetas != NULL) free(Thetas);\n");
  fprintf(fp,"  if (CInvY  != NULL) free(CInvY);\n");
  fprintf(fp,"  if (CMat   != NULL) free(CMat);\n");
  fprintf(fp,"  if (XN     != NULL) free(XN);\n");
  fprintf(fp,"  if (pivs   != NULL) free(pivs);\n");
  fprintf(fp,"}\n");
  fclose(fp);
  printf("FILE psuade_rs.info contains the GP interpolator.\n");
}

// ************************************************************************
// get hyperparameters 
// ------------------------------------------------------------------------
void GP3::getHyperparameters(psVector &inhyper)
{
  inhyper = hyperparameters_;
}


