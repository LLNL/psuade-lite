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
// Functions for the class MGP3
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "GP3.h"
#include "MGP3.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"
#include "Psuade.h"
#include "Sampling.h"
#include "MainEffectAnalyzer.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class MGP3
// ------------------------------------------------------------------------
MGP3::MGP3(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  int    ii, idata;
  double ddata;
  char   *strPtr, winput[1000], equal[100], pString[1000];
  faID_ = PSUADE_RS_MGP2;
  nPartitions_ = 0;

  boxes_ = NULL;
  partSize_ = 500;

  // display banner and additional information
  if (isScreenDumpModeOn() == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,
       "*   Multi-Gaussian Process Function (MGP3) Analysis\n");
    printOutTS(PL_INFO,"* Set printlevel to 1-4 to see details.\n");
    printOutTS(PL_INFO,"* Turn on rs_expert mode to make changes.\n");
    printEquals(PL_INFO, 0);
  }

  // vecXd contains the amount of overlaps between partitions
  vecXd_.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) vecXd_[ii] = 0.05;
  if (psRSExpertMode_ == 1 && psInteractive_ == 1)
  {
    printf("You can improve smoothness across partitions by allowing\n");
    printf("overlaps. The recommended overlap is 0.1 (or 10%%).\n");
    sprintf(pString, "Enter the degree of overlap (0 - 0.4) : ");
    ddata = getDouble(pString);
    if (ddata < 0 || ddata > 0.4)
    {
      ddata = 0.05;
      printf("ERROR: Degree of overlap should be > 0 and <= 0.4.\n");
      printf("INFO:  Degree of overlap set to default = 0.05\n");
    }
    for (ii = 0; ii < nInputs_; ii++) vecXd_[ii] = ddata;
    printf("You can decide the sample size of each partition.\n");
    printf("Larger sample size per partition will take more setup time.\n");
    printf("The default is 1000 (will have more if there is overlap).\n");
    sprintf(pString, "Enter the partition sample size (500 - 10000) : ");
    partSize_ = getInt(500, 20000, pString);
  }

  // user-adjustable parameters
  if (psConfig_ != NULL)
  {
    strPtr = psConfig_->getParameter("MGP_max_samples_per_group");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %d", winput, equal, &partSize_);
      if (partSize_ < 200) 
      {
        printf("MGP3 INFO: config parameter setting not done.\n");
        printf("     max_samples_per_group %d too small.\n",partSize_);
        printf("     max_samples_per_group should be >= 200.\n");
        partSize_ = 200;
      }
    }
  }
  nPartitions_ = (nSamples + partSize_/2) / partSize_;
  ddata = log(1.0*nPartitions_) / log(2.0);
  idata = (int) ddata;
  //if (idata > nInputs_) idata = nInputs_;
  nPartitions_ = 1 << idata;
  printf("MGP3: number of partitions = %d (psize=%d)\n", nPartitions_,
         partSize_);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MGP3::~MGP3()
{
  if (boxes_ != NULL)
  {
    for (int ii = 0; ii < nPartitions_; ii++)
    {
      if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
      if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
      if (boxes_[ii]->rsPtr_   != NULL) delete boxes_[ii]->rsPtr_;
    }
    delete [] boxes_;
  }
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int MGP3::initialize(double *XIn, double *YIn)
{
  int    ii, jj, ss, incr, nSubs, index, samCnt;
  double diff, var, ddata;
  psVector  vecVces, vecXX, vecYY;
  psIVector vecI;

  if (boxes_ != NULL)
  {
    for (ii = 0; ii < nPartitions_; ii++)
    {
      if (boxes_[ii]->lBounds_ != NULL) delete [] boxes_[ii]->lBounds_;
      if (boxes_[ii]->uBounds_ != NULL) delete [] boxes_[ii]->uBounds_;
      if (boxes_[ii]->rsPtr_   != NULL) delete boxes_[ii]->rsPtr_;
    }
    delete [] boxes_;
  }
  boxes_ = NULL;

  if (lowerBounds_ == NULL)
  {
    printOutTS(PL_ERROR,
         "MGP3 initialize ERROR: sample bounds not set yet.\n");
    return -1;
  }

  XData_.setLength(nSamples_*nInputs_);
  for (ii = 0; ii < nSamples_*nInputs_; ii++) XData_[ii] = XIn[ii];
  YData_.setLength(nSamples_);
  for (ii = 0; ii < nSamples_; ii++) YData_[ii] = YIn[ii];
   

  MainEffectAnalyzer *me = new MainEffectAnalyzer();
  int rsExpertTmp = psRSExpertMode_;
  psRSExpertMode_ = 0;
  turnPrintTSOff();
  vecVces.setLength(nInputs_);
  me->computeVCECrude(nInputs_,nSamples_,XIn,YIn,lowerBounds_,
                      upperBounds_, var, vecVces.getDVector());
  delete me;
  psRSExpertMode_ = rsExpertTmp;
  vecI.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) vecI[ii] = ii;
  ddata = 0;
  for (ii = 0; ii < nInputs_; ii++) ddata += vecVces[ii];
  for (ii = 0; ii < nInputs_; ii++) vecVces[ii] /= ddata;
  sortDbleList2a(nInputs_, vecVces.getDVector(), vecI.getIVector());
  if (outputLevel_ > 1)
  {
    for (ii = 0; ii < nInputs_; ii++)
      printf("VCE %d = %e\n", vecI[ii], vecVces[ii]);
  }

  nSubs = int(log(1.0*nPartitions_ + 1.0e-9) / log(2.0));
  boxes_ = new MGP3_Box*[nPartitions_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    boxes_[ii] = new MGP3_Box();
    boxes_[ii]->lBounds_ = new double[nInputs_];
    boxes_[ii]->uBounds_ = new double[nInputs_];
    boxes_[ii]->rsPtr_ = NULL;
    for (jj = 0; jj < nInputs_; jj++)
    {
      boxes_[ii]->lBounds_[jj] = lowerBounds_[jj];
      boxes_[ii]->uBounds_[jj] = upperBounds_[jj];
    }
  }
  for (ii = 0; ii < nSubs; ii++)
  {
    index = vecI[nInputs_-1];
    if (outputLevel_ > 0)
    {
      for (jj = 0; jj < nInputs_; jj++)
        printf("vce %3d = %e\n", vecI[jj]+1, vecVces[jj]);
      printf("Selected input for binary partitioning %d = %d\n",ii+1, 
             index+1);
    }
    incr = 1 << (nSubs - ii - 1);
    for (jj = 0; jj < nPartitions_; jj++)
    {
      if (((jj / incr) % 2) == 0)
      {
        boxes_[jj]->uBounds_[index] = 0.5 *
          (boxes_[jj]->uBounds_[index] + boxes_[jj]->lBounds_[index]);
      }
      else
      {
        boxes_[jj]->lBounds_[index] = 0.5 *
          (boxes_[jj]->uBounds_[index] + boxes_[jj]->lBounds_[index]);
      }
    }
    vecVces[nInputs_-1] *= 0.25;
    sortDbleList2a(nInputs_, vecVces.getDVector(), vecI.getIVector());
  }
  if (outputLevel_ > 3)
  {
    for (jj = 0; jj < nPartitions_; jj++)
    {
      printf("Partition %d:\n", jj);
        for (ii = 0; ii < nInputs_; ii++)
          printf("Input %2d = %12.4e %12.4e\n",ii+1,
            boxes_[jj]->lBounds_[ii], boxes_[jj]->uBounds_[ii]);
    }
  }
  double dcheck1=0, dcheck2=0;
  for (jj = 0; jj < nPartitions_; jj++)
  {
    ddata = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
      ddata *= (boxes_[jj]->uBounds_[ii]-boxes_[jj]->lBounds_[ii]);
    dcheck2 += ddata;
  }
  dcheck1 = 1;
  for (ii = 0; ii < nInputs_; ii++)
    dcheck1 *= (upperBounds_[ii] - lowerBounds_[ii]);
  if (outputLevel_ > 2)
    printf("MGP3: Partition Coverage check: %e (sum) ?= %e (orig)\n",
           dcheck1,dcheck2);

  if (outputLevel_ > 1) printf("MGP3 training begins....\n");
  int total=0;
  vecXX.setLength(nSamples_*nInputs_);
  vecYY.setLength(nSamples_);
  rsExpertTmp = psRSExpertMode_;
  psRSExpertMode_ = 0;
  for (ii = 0; ii < nPartitions_; ii++)
  {
    double *lbs = boxes_[ii]->lBounds_;
    double *ubs = boxes_[ii]->uBounds_;
    samCnt = 0;
    for (ss = 0; ss < nSamples_; ss++)
    {
      for (jj = 0; jj < nInputs_; jj++)
      {
        diff = vecXd_[jj] * (ubs[jj] - lbs[jj]);
        ddata = XData_[ss*nInputs_+jj];
        if (ddata < lbs[jj]-diff || ddata > ubs[jj]+diff) break;
      }
      if (jj == nInputs_)
      {
        for (jj = 0; jj < nInputs_; jj++)
           vecXX[samCnt*nInputs_+jj] = XData_[ss*nInputs_+jj];
        vecYY[samCnt] = YData_[ss];
        samCnt++;
      }
    }
    if (outputLevel_ > 0)
    {
      printf("Partition %d has %d sample points (ovlap=%e).\n",ii+1,
             samCnt, vecXd_[0]);
    }
    if (samCnt == 0)
    {
      printf("MGP3 INFO: some partition has no sample points.\n");
      boxes_[ii]->rsPtr_ = NULL;
    }
    else
    {
      boxes_[ii]->rsPtr_ = new GP3(nInputs_, samCnt);
      boxes_[ii]->rsPtr_->setOutputLevel(0);
      boxes_[ii]->rsPtr_->initialize(vecXX.getDVector(), 
                                     vecYY.getDVector());
    }
    boxes_[ii]->nSamples_ = samCnt;
    total += samCnt;
  }
  psRSExpertMode_ = rsExpertTmp;
  if (outputLevel_ >= 0)
  {
    printf("Sample size = %d\n",nSamples_);
    printf("Total sample sizes from all partitions = %d\n",total);
    printf("INFO: Total from all partitions may be larger than original\n");
    printf("      sample due to overlap. If the total is too large so\n");
    printf("      that it is close to the original size, partitioning\n");
    printf("      is not worthwhile -> you may want to reduce overlap.\n");
  }
  turnPrintTSOn();

  if (outputLevel_ > 1) printf("MGP3 training completed.\n");
  if (psRSCodeGen_ == 1) 
    printf("MGP3 INFO: response surface stand-alone code not available.\n");
  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MGP3::genNDGridData(double *XIn, double *YIn, int *NOut, double **XOut, 
                        double **YOut)
{
  int    totPts;
  double *XX, *YY;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  genNDGrid(NOut, XOut);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  (*YOut) = new double[totPts];
  evaluatePoint(totPts, *XOut, *YOut);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int MGP3::gen1DGridData(double *XIn, double *YIn, int ind1,
                        double *settings, int *NOut, double **XOut, 
                        double **YOut)
{
  int    totPts, ii, kk;
  double HX;
  psVector vecXT;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[totPts];
  vecXT.setLength(totPts*nInputs_);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
     for (kk = 0; kk < nInputs_; kk++) 
        vecXT[ii*nInputs_+kk] = settings[kk]; 
     vecXT[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
     (*XOut)[ii] = HX * ii + lowerBounds_[ind1];
  }
   
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, vecXT.getDVector(), *YOut);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int MGP3::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                        double *settings, int *NOut, double **XOut, 
                        double **YOut)
{
  int totPts, ii, jj, kk, index;
  psVector vecXT, vecHX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2])/(nPtsPerDim_ - 1); 

  vecXT.setLength(totPts*nInputs_);
  (*XOut) = new double[2*totPts];
  checkAllocate(*XOut, "XOut in MGP3::gen2DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      index = ii * nPtsPerDim_ + jj;
      for (kk = 0; kk < nInputs_; kk++) 
        vecXT[index*nInputs_+kk] = settings[kk]; 
      vecXT[index*nInputs_+ind1] = vecHX[0] * ii + lowerBounds_[ind1];
      vecXT[index*nInputs_+ind2] = vecHX[1] * jj + lowerBounds_[ind2];
      (*XOut)[index*2]   = vecHX[0] * ii + lowerBounds_[ind1];
      (*XOut)[index*2+1] = vecHX[1] * jj + lowerBounds_[ind2];
    }
  }
    
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, vecXT.getDVector(), *YOut);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int MGP3::gen3DGridData(double *XIn,double *YIn,int ind1,int ind2,int ind3,
                        double *settings, int *NOut, double **XOut, 
                        double **YOut)
{
  int totPts, ii, jj, kk, ll, index;
  psVector vecXT, vecHX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_-1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_-1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_-1); 

  vecXT.setLength(totPts*nInputs_);
  (*XOut) = new double[3*totPts];
  checkAllocate(*XOut, "XOut in MGP3::gen3DGrid");
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      for (ll = 0; ll < nPtsPerDim_; ll++) 
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        for (kk = 0; kk < nInputs_; kk++) 
          vecXT[index*nInputs_+kk] = settings[kk]; 
        vecXT[index*nInputs_+ind1] = vecHX[0] * ii + lowerBounds_[ind1];
        vecXT[index*nInputs_+ind2] = vecHX[1] * jj + lowerBounds_[ind2];
        vecXT[index*nInputs_+ind3] = vecHX[2] * ll + lowerBounds_[ind3];
        (*XOut)[index*3]   = vecHX[0] * ii + lowerBounds_[ind1];
        (*XOut)[index*3+1] = vecHX[1] * jj + lowerBounds_[ind2];
        (*XOut)[index*3+2] = vecHX[2] * ll + lowerBounds_[ind3];
      }
    }
  }
    
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, vecXT.getDVector(), *YOut);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int MGP3::gen4DGridData(double *XIn,double *YIn,int ind1,int ind2,int ind3,
                        int ind4,double *settings,int *NOut,double **XOut, 
                        double **YOut)
{
  int totPts, ii, jj, kk, ll, mm, index;
  psVector vecXT, vecHX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
 
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_-1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_-1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_-1); 
  vecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_-1); 

  vecXT.setLength(totPts*nInputs_);
  (*XOut) = new double[4*totPts];
  checkAllocate(*XOut, "XOut in MGP3::gen4DGrid");
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
            vecXT[index*nInputs_+kk] = settings[kk]; 
          vecXT[index*nInputs_+ind1] = vecHX[0] * ii + lowerBounds_[ind1];
          vecXT[index*nInputs_+ind2] = vecHX[1] * jj + lowerBounds_[ind2];
          vecXT[index*nInputs_+ind3] = vecHX[2] * ll + lowerBounds_[ind3];
          vecXT[index*nInputs_+ind4] = vecHX[3] * mm + lowerBounds_[ind4];
          (*XOut)[index*4]   = vecHX[0] * ii + lowerBounds_[ind1];
          (*XOut)[index*4+1] = vecHX[1] * jj + lowerBounds_[ind2];
          (*XOut)[index*4+2] = vecHX[2] * ll + lowerBounds_[ind3];
          (*XOut)[index*4+3] = vecHX[3] * mm + lowerBounds_[ind4];
        }
      }
    }
  }
    
  (*YOut) = new double[totPts];
  if (outputLevel_ >= 1) printf("MGP3 interpolation begins....\n");
  evaluatePoint(totPts, vecXT.getDVector(), *YOut);
  if (outputLevel_ >= 1) printf("MGP3 interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double MGP3::evaluatePoint(double *X)
{
  int    iOne=1;
  double Y=0.0;
  interpolate(iOne, X, &Y, NULL);
  return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double MGP3::evaluatePoint(int npts, double *X, double *Y)
{
  interpolate(npts, X, Y, NULL);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double MGP3::evaluatePointFuzzy(double *X, double &std)
{
  int    iOne=1;
  double Y=0.0;
  interpolate(iOne, X, &Y, &std);
  return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double MGP3::evaluatePointFuzzy(int npts,double *X, double *Y,double *Ystds)
{
  interpolate(npts, X, Y, Ystds);
  return 0.0;
}

// ************************************************************************
// interpolation 
// ------------------------------------------------------------------------
int MGP3::interpolate(int npts, double *X, double *Y, double *Ystds)
{
  int    ii, pp, ss, count, highFlag, hasStds=0;
  double Yt, ddata, diff, *lbs, *ubs;

  if (Ystds != NULL) hasStds = 1;
  for (ss = 0; ss < npts; ss++)
  {
    count = 0;
    Yt = 0.0;
    for (pp = 0; pp < nPartitions_; pp++)
    {
      lbs = boxes_[pp]->lBounds_;
      ubs = boxes_[pp]->uBounds_;
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (ubs[ii] == upperBounds_[ii]) highFlag = 1;
        else                             highFlag = 0;
        ddata = X[ss*nInputs_+ii];
        diff = 0.05 * (ubs[ii] - lbs[ii]);
        if (highFlag == 0)
        {
          if (ddata < (lbs[ii]-diff) || ddata >= (ubs[ii]+diff)) break;
        }
        else
        {
          if (ddata < (lbs[ii]-diff) || ddata > (ubs[ii]+diff)) break;
        }
      }
      if (ii == nInputs_)
      {
        if (hasStds) 
          Yt += boxes_[pp]->rsPtr_->evaluatePointFuzzy(&X[ss*nInputs_],
                                                       Ystds[ss]);
        else
          Yt += boxes_[pp]->rsPtr_->evaluatePoint(&X[ss*nInputs_]);
        count++;
      }
    }
    if (count == 0)
    {
      printf("MGP3 evaluate ERROR: sample point outside range.\n");
      for (ii = 0; ii < nInputs_; ii++)
        printf("Input %d = %e (in [%e, %e]?)\n",ii+1,
           X[ss*nInputs_+ii],lowerBounds_[ii],upperBounds_[ii]);
      return 0.0;
    }
    Y[ss] = Yt / (double) count;
  }
  return 0;
}

