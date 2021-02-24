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
// Functions for the class MMars (for large data set)
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "MMars.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"
#include "Psuade.h"
#include "MainEffectAnalyzer.h"

// ************************************************************************
// Constructor for object class Multi-domain Mars
// ------------------------------------------------------------------------
MMars::MMars(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  int    idata, ii;
  double ddata;
  char   pString[501], *strPtr, equal[100], winput[5000];

  boxes_ = NULL;
  partSize_ = 6000;

  // display banner and additonal information
  if (isScreenDumpModeOn() == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,
       "*   Multi-Multivariate Regression Function (MMARS) Analysis\n");
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
    sprintf(pString, "Enter the partition sample size (1000 - 10000) : ");
    partSize_ = getInt(500, 20000, pString);
  }

  if (psConfig_ != NULL)
  {
    strPtr = psConfig_->getParameter("MMARS_max_samples_per_group");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %d", winput, equal, &idata);
      if (idata >= 1000) partSize_ = idata;
      else 
      {
        printf("MMars INFO: config parameter setting not done.\n");
        printf("            max_samples_per_group %d too small.\n",idata);
        printf("            max_samples_per_group should be >= 1000.\n");
      }
    }
  }
  nPartitions_ = (nSamples + partSize_/2) / partSize_;
  ddata = log(1.0*nPartitions_) / log(2.0);
  idata = (int) ddata;
  //if (idata > nInputs_) idata = nInputs_;
  nPartitions_ = 1 << idata;
  printf("MMars: number of partitions = %d\n", nPartitions_);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MMars::~MMars()
{
  if (boxes_ != NULL)
  {
    for (int ii = 0; ii < nPartitions_; ii++) 
    {
      if (boxes_[ii] != NULL)
        if (boxes_[ii]->marsPtr_ != NULL) delete boxes_[ii]->marsPtr_;
    }
    delete [] boxes_;
  }
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int MMars::initialize(double *X, double *Y)
{
  int    ii, jj, ss, nSubs, index, incr, samCnt;
  double range, var=0,ddata, diff;
  char   pString[10000], winput[10000];
  psVector  vecVces, vecXT, vecYT;
  psIVector vecInds;

  if (boxes_ != NULL)
  {
    for (ii = 0; ii < nPartitions_; ii++) 
    {
      if (boxes_[ii] != NULL)
        if (boxes_[ii]->marsPtr_ != NULL) delete boxes_[ii]->marsPtr_;
    }
    delete [] boxes_;
  }
  boxes_ = NULL;

  if (lowerBounds_ == NULL)
  {
    printOutTS(PL_ERROR,
         "MMars initialize ERROR: sample bounds not set yet.\n");
    return -1;
  }

  MainEffectAnalyzer *me = new MainEffectAnalyzer();
  int rsExpertTmp = psAnaExpertMode_;
  psAnaExpertMode_ = 0;
  turnPrintTSOff();
  vecVces.setLength(nInputs_);
  me->computeVCECrude(nInputs_,nSamples_,X,Y,lowerBounds_,upperBounds_, 
                      var, vecVces.getDVector());
  delete me;
  psAnaExpertMode_ = rsExpertTmp;
  vecInds.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) vecInds[ii] = ii;
  ddata = 0;
  for (ii = 0; ii < nInputs_; ii++) ddata += vecVces[ii];
  for (ii = 0; ii < nInputs_; ii++) vecVces[ii] /= ddata;
  sortDbleList2a(nInputs_, vecVces.getDVector(), vecInds.getIVector());
  if (outputLevel_ > 1) 
  {
    for (ii = 0; ii < nInputs_; ii++) 
      printf("VCE %4d = %12.4e\n", vecInds[ii], vecVces[ii]);
  }

  int nBasis, maxVarPerBasis, normalizeY;
  if (psRSExpertMode_ == 1 && isScreenDumpModeOn() == 1)
  {
    printOutTS(PL_INFO,"MMars: Current number of basis functions = 100\n");
    nBasis = nSamples_;
    if (nSamples_ > 10)
    {
      sprintf(pString,"Enter the number of basis functions (>=10, <= %d): ",
              nSamples_);
      nBasis = getInt(10, nSamples_, pString);
    }
    maxVarPerBasis = 8;
    if (nInputs_ < maxVarPerBasis) maxVarPerBasis = nInputs_;
    printOutTS(PL_INFO,
         "MMars: Current degree of interactions    = %d\n",maxVarPerBasis);
    sprintf(pString, "Enter the degree of interactions (<=%d) : ",nInputs_);
    maxVarPerBasis = getInt(1, nInputs_, pString);
    sprintf(pString, "Mars: normalize output? (y or n) ");
    getString(pString, winput);
    normalizeY = 0;
    if (winput[0] == 'y') normalizeY = 1;
    if (psConfig_ != NULL)
    {
      sprintf(pString, "MARS_num_basis = %d", nBasis);
      psConfig_->putParameter(pString);
      sprintf(pString, "MARS_interaction = %d", maxVarPerBasis);
      psConfig_->putParameter(pString);
      if (normalizeY == 1)
      {
        sprintf(pString, "normalize_outputs");
        psConfig_->putParameter(pString);
      }
    }
    rsExpertTmp = psRSExpertMode_;
    psRSExpertMode_ = 0;
  }

  nSubs = int(log(1.0*nPartitions_ + 1.0e-9) / log(2.0));
  boxes_ = new MMars_Box*[nPartitions_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    boxes_[ii] = new MMars_Box();
    boxes_[ii]->vecLBs_.setLength(nInputs_);
    boxes_[ii]->vecUBs_.setLength(nInputs_);
    boxes_[ii]->marsPtr_ = NULL;
    for (jj = 0; jj < nInputs_; jj++)
    {
      boxes_[ii]->vecLBs_[jj] = lowerBounds_[jj];
      boxes_[ii]->vecUBs_[jj] = upperBounds_[jj];
    }
  }

  for (ii = 0; ii < nSubs; ii++)
  {
    index = vecInds[nInputs_-1]; 
    if (outputLevel_ > 0) 
    {
      printf("Selected input for partition %4dd = %4dd, VCE = %12.4e\n",
             ii+1, index+1, vecVces[nInputs_-1]);
    }
    incr = 1 << (nSubs - ii - 1);
    for (jj = 0; jj < nPartitions_; jj++)
    {
      if (((jj / incr) % 2) == 0)
      { 
        boxes_[jj]->vecUBs_[index] = 0.5 * 
          (boxes_[jj]->vecUBs_[index] + boxes_[jj]->vecLBs_[index]);
      }
      else
      { 
        boxes_[jj]->vecLBs_[index] = 0.5 * 
          (boxes_[jj]->vecUBs_[index] + boxes_[jj]->vecLBs_[index]);
      }
    }
    vecVces[nInputs_-1] *= 0.25;
    sortDbleList2a(nInputs_, vecVces.getDVector(), vecInds.getIVector());
  }
  if (outputLevel_ > 3) 
  {
    for (jj = 0; jj < nPartitions_; jj++)
    {
      printf("Partition %d:\n", jj);
        for (ii = 0; ii < nInputs_; ii++)
          printf("Input %2d = %12.4e %12.4e\n",ii+1,
            boxes_[jj]->vecLBs_[ii], boxes_[jj]->vecUBs_[ii]);
    }
  }
  double dcheck1=0, dcheck2=0;
  for (jj = 0; jj < nPartitions_; jj++)
  {
    ddata = 1.0;
    for (ii = 0; ii < nInputs_; ii++)
      ddata *= (boxes_[jj]->vecUBs_[ii]-boxes_[jj]->vecLBs_[ii]);
    dcheck2 += ddata;
  }
  dcheck1 = 1;
  for (ii = 0; ii < nInputs_; ii++)
    dcheck1 *= (upperBounds_[ii] - lowerBounds_[ii]);
  printf("MMars: Partition coverage check: %e (sum) ?= %e (orig)\n", 
         dcheck1, dcheck2);

  if (outputLevel_ > 1) printf("MMars training begins....\n");
  int total=0;
  double *lbs, *ubs;
  vecXT.setLength(nSamples_*nInputs_);
  vecYT.setLength(nSamples_);
  for (ii = 0; ii < nPartitions_; ii++)
  {
    lbs = boxes_[ii]->vecLBs_.getDVector();
    ubs = boxes_[ii]->vecUBs_.getDVector();
    samCnt = 0;
    for (ss = 0; ss < nSamples_; ss++)
    {
      for (jj = 0; jj < nInputs_; jj++)
      {
        diff = vecXd_[jj] * (ubs[jj] - lbs[jj]);
        ddata = X[ss*nInputs_+jj];
        if (ddata < lbs[jj]-diff || ddata > ubs[jj]+diff) break;
      } 
      if (jj == nInputs_)
      {
        for (jj = 0; jj < nInputs_; jj++)
          vecXT[samCnt*nInputs_+jj] = X[ss*nInputs_+jj];
        vecYT[samCnt] = Y[ss];
        samCnt++;
      }
    }
    if (outputLevel_ >= 0) 
      printf("Partition %d has %d sample points.\n",ii+1,samCnt);
    if (samCnt == 0)
    {
      printf("MMars INFO: some partition has no sample points.\n");
      boxes_[ii]->marsPtr_ = NULL;
    }
    else
    {
      boxes_[ii]->marsPtr_ = new Mars(nInputs_, samCnt);
      boxes_[ii]->marsPtr_->initialize(vecXT.getDVector(),
                                       vecYT.getDVector());
    }
    boxes_[ii]->nSamples_ = samCnt;
    total += samCnt;
  }
  if (outputLevel_ > 0) 
  {
    printf("Original sample size = %d\n",nSamples_);
    printf("Total sample sizes from all partitions = %d\n",total);
    printf("INFO: Total from all partitions may be larger than original\n");
    printf("      sample due to overlap. If the total is too large so\n");
    printf("      that it is close to the original size, partitioning\n");
    printf("      is not worthwhile -> you may want to reduce overlap.\n");
  }
  if (rsExpertTmp == 1 && isScreenDumpModeOn() == 1)
  {
    if (psConfig_ != NULL)
    {
      psConfig_->removeParameter("MARS_num_basis");
      psConfig_->removeParameter("MARS_interaction");
      psConfig_->removeParameter("normalize_outputs");
    }
    psRSExpertMode_ = rsExpertTmp;
  }
  turnPrintTSOn();

  if (outputLevel_ > 1) printf("MMars training completed.\n");
  if (psRSCodeGen_ == 1)
    printf("MMars INFO: response surface stand-alone code not available.\n");

  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MMars::genNDGridData(double *XIn,double *YIn,int *NOut,double **XOut,
                         double **YOut)
{
  int totPts;

  initialize(XIn,YIn);

  if ((*NOut) == -999) return 0;
 
  genNDGrid(NOut, XOut);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in MMars::genNDGridData");
  evaluatePoint(totPts, *XOut, *YOut);

  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MMars::gen1DGridData(double *XIn,double *YIn,int ind1,double *settings, 
                         int *NOut, double **XOut, double **YOut)
{
  int    ii, ss, totPts;
  double HX;
  psVector vecXT;

  initialize(XIn,YIn);

  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[totPts];
  (*YOut) = new double[totPts];
  (*NOut) = totPts;

  vecXT.setLength(totPts*nInputs_);
  for (ss = 0; ss < totPts; ss++) 
    for (ii = 0; ii < nInputs_; ii++) vecXT[ss*nInputs_+ii] = settings[ii]; 
   
  for (ss = 0; ss < totPts; ss++) 
  {
    vecXT[ss*nInputs_+ind1]  = HX * ss + lowerBounds_[ind1];
    (*XOut)[ss] = HX * ss + lowerBounds_[ind1];
    (*YOut)[ss] = 0.0;
  }

  evaluatePoint(totPts, vecXT.getDVector(), *YOut);

  return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MMars::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                         double *settings, int *NOut, double **XOut, 
                         double **YOut)
{
  int ii, ss, jj, index, totPts;
  psVector vecXT, vecHX;
 
  initialize(XIn,YIn);

  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2])/(nPtsPerDim_ - 1); 

  (*XOut) = new double[2*totPts];
  (*YOut) = new double[totPts];
  (*NOut) = totPts;

  vecXT.setLength(totPts*nInputs_);
  for (ss = 0; ss < totPts; ss++) 
    for (ii = 0; ii < nInputs_; ii++) vecXT[ss*nInputs_+ii] = settings[ii]; 
    
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      index = ii * nPtsPerDim_ + jj;
      vecXT[index*nInputs_+ind1] = vecHX[0] * ii + lowerBounds_[ind1];
      vecXT[index*nInputs_+ind2] = vecHX[1] * jj + lowerBounds_[ind2];
      (*XOut)[index*2]   = vecHX[0] * ii + lowerBounds_[ind1];
      (*XOut)[index*2+1] = vecHX[1] * jj + lowerBounds_[ind2];
    }
  }

  evaluatePoint(totPts, vecXT.getDVector(), *YOut);

  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int MMars::gen3DGridData(double *XIn, double *YIn, int ind1, int ind2,
                         int ind3, double *settings, int *NOut, 
                         double **XOut, double **YOut)
{
  int ii, ss, jj, ll, index, totPts;
  psVector vecXT, vecHX;

  initialize(XIn,YIn);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2])/(nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3])/(nPtsPerDim_ - 1); 

  (*XOut) = new double[3*totPts];
  (*YOut) = new double[totPts];
  (*NOut) = totPts;

  vecXT.setLength(totPts*nInputs_);
  for (ss = 0; ss < totPts; ss++)
    for (ii = 0; ii < nInputs_; ii++) vecXT[ss*nInputs_+ii] = settings[ii];

  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++)
    {
      for (ll = 0; ll < nPtsPerDim_; ll++)
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        vecXT[index*nInputs_+ind1]  = vecHX[0] * ii + lowerBounds_[ind1];
        vecXT[index*nInputs_+ind2]  = vecHX[1] * jj + lowerBounds_[ind2];
        vecXT[index*nInputs_+ind3]  = vecHX[2] * ll + lowerBounds_[ind3];
        (*XOut)[index*3]   = vecHX[0] * ii + lowerBounds_[ind1];
        (*XOut)[index*3+1] = vecHX[1] * jj + lowerBounds_[ind2];
        (*XOut)[index*3+2] = vecHX[2] * ll + lowerBounds_[ind3];
      }
    }
  }

  evaluatePoint(totPts, vecXT.getDVector(), *YOut);

  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int MMars::gen4DGridData(double *XIn,double *YIn, int ind1, int ind2, 
                         int ind3, int ind4, double *settings, int *NOut, 
                         double **XOut, double **YOut)
{
  int ii, ss, jj, ll, mm, index, totPts;
  psVector vecXT, vecHX;

  initialize(XIn,YIn);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2])/(nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3])/(nPtsPerDim_ - 1); 
  vecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4])/(nPtsPerDim_ - 1); 

  (*XOut) = new double[4*totPts];
  (*YOut) = new double[totPts];
  (*NOut) = totPts;

  vecXT.setLength(totPts*nInputs_);
  for (ss = 0; ss < totPts; ss++) 
    for (ii = 0; ii < nInputs_; ii++) vecXT[ss*nInputs_+ii] = settings[ii]; 
    
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

  evaluatePoint(totPts, vecXT.getDVector(), *YOut);

  return 0;
}

// ************************************************************************
// Evaluate a point
// ------------------------------------------------------------------------
double MMars::evaluatePoint(double *X)
{
  int    iOne=1;
  double Y;
  evaluatePoint(1, X, &Y);
  return Y;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double MMars::evaluatePoint(int nPts, double *X, double *Y)
{
  int    ss, pp, ii, nSamp, count, highFlag;
  double diff, Yt, ddata, *lbs, *ubs;

  for (ss = 0; ss < nPts; ss++) 
  {
    count = 0;
    Yt = 0.0;
    for (pp = 0; pp < nPartitions_; pp++)
    {
      lbs = boxes_[pp]->vecLBs_.getDVector();
      ubs = boxes_[pp]->vecUBs_.getDVector();
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (ubs[ii] == upperBounds_[ii]) highFlag = 1;
        else                             highFlag = 0;
        ddata = X[ss*nInputs_+ii];
        diff = vecXd_[ii] * (ubs[ii] - lbs[ii]);
        if (highFlag == 0)
        {
          if (ddata < (lbs[ii]-diff) || ddata >= (ubs[ii]+diff)) break;
        }
        else 
        {
          if (ddata < (lbs[ii]-diff) || ddata > (ubs[ii]+diff)) break;
        }
      } 
      if (ii == nInputs_ && boxes_[pp]->nSamples_ > 0)
      {
        if (boxes_[pp]->marsPtr_ != NULL)
        {
          Yt += boxes_[pp]->marsPtr_->evaluatePoint(&X[ss*nInputs_]);
          count++;
        }
      } 
    }
    if (count == 0)
    {
      printf("MMars evaluate WARNING: sample point not in any partition.\n");
      printf("INFO: this may happen during cross validation.\n");
      printf("INFO: prediction - average of all partitions (maybe wrong).\n");
      Yt = 0;
      for (pp = 0; pp < nPartitions_; pp++)
      {
        if (boxes_[pp]->marsPtr_ != NULL)
        {
          Yt += boxes_[pp]->marsPtr_->evaluatePoint(&X[ss*nInputs_]);
          count++;
        }
      }
    }
    Y[ss] = Yt / (double) count;
  }
  return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double MMars::evaluatePointFuzzy(double *X, double &std)
{
  int    iOne=1;
  double Y=0.0;
  evaluatePoint(iOne, X, &Y);
  std = 0.0;
  return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double MMars::evaluatePointFuzzy(int npts,double *X,double *Y,double *Ystd)
{
  evaluatePoint(npts, X, Y);
  for (int ss = 0; ss < npts; ss++) Ystd[ss] = 0.0;
  return 0.0;
}

