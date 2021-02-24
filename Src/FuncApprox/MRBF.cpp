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
// Functions for the class RBF
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "MRBF.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "MainEffectAnalyzer.h"
#include "PrintingTS.h"

extern "C" {
#if 0
   void dgetrf_(int *, int *, double *, int *, int *, int *);
   void dgetrs_(char *,int *,int*,double*,int*,int*,double*,int*,int*);
#endif
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
}

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#define PS_RBF1

// ************************************************************************
// Constructor for object class RBF
// ------------------------------------------------------------------------
MRBF::MRBF(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  int    idata, ii;
  double ddata;
  char   pString[501], *strPtr, equal[100], winput[5000];

  type_ = 0;

  svdThresh_ = 1e-15;
  gaussScale_ = 1;
  boxes_ = NULL;
  partSize_ = 500;

  // display banner and additonal information
  if (isScreenDumpModeOn() == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,
         "*   Multiple Radial Basis Function (RBF) Analysis\n");
    printOutTS(PL_INFO,"* Set printlevel to 1-4 to see RBF details.\n");
    printOutTS(PL_INFO,"* Default kernel    = multi-quadratic \n");
    printOutTS(PL_INFO,
         "* Default threshold = 1.0e-15 (for SVD truncation)\n");
    printOutTS(PL_INFO,"* Turn on rs_expert mode to make changes.\n");
    printEquals(PL_INFO, 0);
  }
  
  vecXd_.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) vecXd_[ii] = 0.05;
  if (psRSExpertMode_ == 1 && psInteractive_ == 1)
  {
    printf("In the following you have the option to select the kernel. \n");
    printf("0. multi-quadratic\n");
    printf("1. inverse multi-quadratic\n");
    printf("2. Gaussian\n");
    printf("3. thin plate spline\n");
    sprintf(pString,"Enter your choice (0 - 3) : ");
    type_ = getInt(0, 3, pString);
    if (type_ == 2)
    {
      sprintf(pString,
         "Enter scaling factor for Gaussian kernel (default=1) : ");
      gaussScale_ = getDouble(pString);
    }
    printOutTS(PL_INFO,
         "The RBF matrix to be constructed may be near-singular.\n");
    printOutTS(PL_INFO,
         "Currently, singular values < max(svd)*1e-15 are truncated.\n");
    printOutTS(PL_INFO,
         "You have the option to change this threshold (1e-15).\n");
    printOutTS(PL_INFO,
         "NOTE: truncating singular values may lead to erroneous results.\n");
    sprintf(pString, "Enter new threshold for SVD (> 0 but << 1) : ");
    svdThresh_ = getDouble(pString);
    printf("You can improve smoothness across partitions by allowing\n");
    printf("overlaps. The recommended overlap is 0.1 (or 10%%).\n");
    sprintf(pString, "Enter the degree of overlap (0 - 0.4) : ");
    ddata = getDouble(pString);
    if (ddata < 0 || ddata > 0.4)
    {
      ddata = 0.1;
      printf("ERROR: Degree of overlap should be > 0 and <= 0.4.\n");
      printf("INFO:  Degree of overlap set to default = 0.05\n");
    }
    for (ii = 0; ii < nInputs_; ii++) vecXd_[ii] = ddata;
    printf("You can decide the sample size of each partition.\n");
    printf("Larger sample size per partition will take more setup time.\n");
    printf("The default is 500 (will have more if there is overlap).\n");
    sprintf(pString, "Enter the partition sample size (500 - 2000) : ");
    partSize_ = getInt(500, 20000, pString);
  }

  if (psConfig_ != NULL)
  {
    strPtr = psConfig_->getParameter("MRBF_max_samples_per_group");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %d", winput, equal, &idata);
      if (idata >= 100) partSize_ = idata;
      else
      {
        printf("MRBF INFO: config parameter setting not done.\n");
        printf("           max_samples_per_group %d too small.\n",idata);
        printf("           max_samples_per_group should be >= 100.\n");
      }
    }
  }
  nPartitions_ = (nSamples + partSize_/2) / partSize_;
  ddata = log(1.0*nPartitions_) / log(2.0);
  idata = (int) ddata;
  //if (idata > nInputs_) idata = nInputs_;
  nPartitions_ = 1 << idata;
  printf("MRBF: number of partitions = %d\n", nPartitions_);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MRBF::~MRBF()
{
  if (boxes_ != NULL) delete [] boxes_;
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int MRBF::initialize(double *X, double *Y)
{
  int    ii, jj, ss, nSubs, index, incr, samCnt;
  double range, var=0,ddata, diff;
  psVector  vecVces, vecXT, vecYT;
  psIVector vecInds;

  if (boxes_ != NULL) delete [] boxes_;
  boxes_ = NULL;

  if (lowerBounds_ == NULL)
  {
    printOutTS(PL_ERROR,
         "MRBF initialize ERROR: sample bounds not set yet.\n");
    return -1;
  }

  vecXN_.setLength(nSamples_*nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
  {
    range = 1.0 / (upperBounds_[ii] - lowerBounds_[ii]);
    for (ss = 0; ss < nSamples_; ss++)
      vecXN_[ss*nInputs_+ii] = 
         (X[ss*nInputs_+ii] - lowerBounds_[ii]) * range;
  }
  vecYN_.setLength(nSamples_);
  initOutputScaling(Y, vecYN_.getDVector());
  for (ii = 0; ii < nSamples_; ii++) vecYN_[ii] = Y[ii] - YMean_;

  MainEffectAnalyzer *me = new MainEffectAnalyzer();
  turnPrintTSOff();
  vecVces.setLength(nInputs_);
  me->computeVCECrude(nInputs_,nSamples_,X,Y,lowerBounds_,upperBounds_, 
                      var, vecVces.getDVector());
  delete me;
  vecInds.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) vecInds[ii] = ii;
  sortDbleList2a(nInputs_, vecVces.getDVector(), vecInds.getIVector());
  if (outputLevel_ > 1) 
  {
    for (ii = 0; ii < nInputs_; ii++) 
      printf("VCE %d = %e\n", vecInds[ii], vecVces[ii]);
  }

  nSubs = int(log(1.0*nPartitions_ + 1.0e-9) / log(2.0));
  boxes_ = new MRBF_Box*[nPartitions_];
  for (ii = 0; ii < nPartitions_; ii++)
  {
    boxes_[ii] = new MRBF_Box();
    boxes_[ii]->vecLBs_.setLength(nInputs_);
    boxes_[ii]->vecUBs_.setLength(nInputs_);
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
      for (jj = 0; jj < nInputs_; jj++)
        printf("vce %3d = %e\n", vecInds[jj]+1, vecVces[jj]);
      printf("Selected input for partition = %d\n", index+1);
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
  printf("MRBF: Partition coverage check: %e (sum) ?= %e (orig)\n", 
         dcheck1, dcheck2);

  if (outputLevel_ > 1) printf("MRBF training begins....\n");
  int total=0;
  double *lbs, *ubs, *coefs;
  vecYT.setLength(nSamples_);
  for (ii = 0; ii < nPartitions_; ii++)
  {
    lbs = boxes_[ii]->vecLBs_.getDVector();
    ubs = boxes_[ii]->vecUBs_.getDVector();
    vecXT.setLength(nSamples_*nInputs_);
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
           vecXT[samCnt*nInputs_+jj] = vecXN_[ss*nInputs_+jj];
        vecYT[samCnt] = vecYN_[ss];
        samCnt++;
      }
    }
    if (outputLevel_ > 0) 
      printf("Partition %d has %d sample points.\n",ii+1,samCnt);
    if (samCnt == 0)
    {
      printf("MRBF INFO: some partition has no sample points.\n");
    }
    else
    {
      initialize(samCnt, vecXT.getDVector(),vecYT.getDVector(),&coefs);
      boxes_[ii]->vecCoeffs_.load(samCnt+1, coefs);
    }
    boxes_[ii]->vecX_ = vecXT;
    boxes_[ii]->nSamples_ = samCnt;
    total += samCnt;
  }
  if (outputLevel_ >= 0) 
  {
    printf("Original sample size = %d\n",nSamples_);
    printf("Total sample sizes from all partitions = %d\n",total);
    printf("INFO: Total from all partitions may be larger than original\n"); 
    printf("      sample due to overlap. If the total is too large so\n");
    printf("      that it is close to the original size, partitioning\n"); 
    printf("      is not worthwhile -> you may want to reduce overlap.\n");
  }
  turnPrintTSOn();

  if (outputLevel_ > 1) printf("MRBF training completed.\n");
  if (psRSCodeGen_ == 1)
    printf("MRBF INFO: response surface stand-alone code not available.\n");

  return 0;
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int MRBF::initialize(int nSamples, double *X, double *Y, double **coefs)
{
  int    ii, kk, ss, ss2, nSamp1;
  double ddata;
  psVector vecDmat;

#ifdef PS_RBF1
  nSamp1 = nSamples + 1;
#else
  nSamp1 = nSamples;
#endif
  vecDmat.setLength(nSamp1*nSamp1);
  switch(type_) 
  {
    case 0: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = multi-quadratic\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        vecDmat[ss*nSamp1+ss] = 1.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          vecDmat[ss*nSamp1+ss2] = 
                   vecDmat[ss2*nSamp1+ss] = sqrt(ddata+1.0);
        }
      }
      break;

    case 1: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = inverse multi-quadratic\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        vecDmat[ss*nSamp1+ss] = 1.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          vecDmat[ss*nSamp1+ss2] = 
                vecDmat[ss2*nSamp1+ss] = 1.0/sqrt(ddata+1.0);
        }
      }
      break;

    case 2: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = Gaussian\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        vecDmat[ss*nSamp1+ss] = 1.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          vecDmat[ss*nSamp1+ss2] = 
            vecDmat[ss2*nSamp1+ss] = exp(-gaussScale_*ddata/2.0);
        }
      }
      break;

    case 3: 
      if (outputLevel_ > 0) 
        printOutTS(PL_INFO,"Kernel = thin plate spline\n");
      for (ss = 0; ss < nSamples; ss++)
      {
        vecDmat[ss*nSamp1+ss] = 0.0; 
        for (ss2 = ss+1; ss2 < nSamples; ss2++)
        {
          ddata = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
            ddata += pow((X[ss*nInputs_+ii]-X[ss2*nInputs_+ii]),2.0);
          vecDmat[ss*nSamp1+ss2] = 
             vecDmat[ss2*nSamp1+ss] = (ddata+1.0)*log(sqrt(ddata+1.0));
        }
      }
      break;
  }
#ifdef PS_RBF1
  for (ss = 0; ss < nSamples; ss++)
    vecDmat[ss*nSamp1+nSamples] = vecDmat[nSamples*nSamp1+ss] = 1.0; 
  vecDmat[nSamples*nSamp1+nSamples] = 0.0;
#endif

  int info, cnt=0;
  psMatrix MatA, MatU, MatV;
  psVector vecS, vecW;
  if (outputLevel_ > 3) printf("Running SVD ...\n");
  MatA.load(nSamp1, nSamp1, vecDmat.getDVector());
  info = MatA.computeSVD(MatU, vecS, MatV);
  if (outputLevel_ > 3)
    printf("SVD completed: status = %d (should be 0).\n",info);

  if (info != 0) 
  {
    printOutTS(PL_WARN,"RBF ERROR: dgesvd returns error %d.\n",info);
    return -1;
  }

  double *regCoeffs = new double[nSamp1];
  checkAllocate(regCoeffs, "regCoeffs in MRBF::initialize");
  for (ii = 0; ii < nSamples; ii++) regCoeffs[ii] = Y[ii];

#ifdef PS_RBF1
  regCoeffs[nSamples] = 0.0;
#endif
  vecW.setLength(nSamp1);
  double *UU = MatU.getMatrix1D();
  double *VV = MatV.getMatrix1D();
  for (ss = 1; ss < nSamp1; ss++)
  {
    if (vecS[ss]/vecS[0] < svdThresh_)
    {
      vecS[ss] = 0;
      cnt++;
    }
  }
  if (cnt > 0 && psInteractive_ == 1 && outputLevel_ > 0) 
  {
    printOutTS(PL_WARN,
         "WARNING: RBF matrix is near-singular. Small singular values\n");
    printOutTS(PL_WARN,
         "         (%d out of %d) are truncated.\n",cnt,nSamp1);
    printOutTS(PL_WARN,"         Approximation may be inaccurate.\n");
  }
  for (ss = 0; ss < nSamp1; ss++)
  {
    vecW[ss] = 0.0;
    for (ss2 = 0; ss2 < nSamp1; ss2++)
      vecW[ss] += UU[ss*nSamp1+ss2] * regCoeffs[ss2];
  }
  for (ss = 0; ss < nSamp1; ss++) 
  {
    if (vecS[ss] != 0) vecW[ss] /= vecS[ss];
    else               vecW[ss] = 0;
  }
  for (ss = 0; ss < nSamp1; ss++)
  {
    regCoeffs[ss] = 0.0;
    for (ss2 = 0; ss2 < nSamp1; ss2++) 
      regCoeffs[ss] += VV[ss*nSamp1+ss2] * vecW[ss2];
  }
  (*coefs) = regCoeffs;
  return 0;
}
  
// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MRBF::genNDGridData(double *XIn,double *YIn,int *NOut,double **XOut,
                        double **YOut)
{
   int totPts;

   initialize(XIn,YIn);

   if ((*NOut) == -999) return 0;
  
   genNDGrid(NOut, XOut);
   if ((*NOut) == 0) return 0;
   totPts = (*NOut);

   (*YOut) = new double[totPts];
   checkAllocate(*YOut, "YOut in MRBF::genNDGridData");
   evaluatePoint(totPts, *XOut, *YOut);

   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int MRBF::gen1DGridData(double *XIn,double *YIn,int ind1,double *settings, 
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
      for (ii = 0; ii < nInputs_; ii++) 
        vecXT[ss*nInputs_+ii] = settings[ii]; 
    
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
int MRBF::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                        double *settings, int *NOut, double **XOut, 
                        double **YOut)
{
   int    ii, ss, jj, index, totPts;
   psVector vecXT, vecHX;
 
   initialize(XIn,YIn);

   totPts = nPtsPerDim_ * nPtsPerDim_;
   vecHX.setLength(2);
   vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_-1); 
   vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_-1); 

   (*XOut) = new double[2*totPts];
   (*YOut) = new double[totPts];
   (*NOut) = totPts;

   vecXT.setLength(totPts*nInputs_);
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) 
        vecXT[ss*nInputs_+ii] = settings[ii]; 
    
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
int MRBF::gen3DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                        int ind3, double *settings, int *NOut, 
                        double **XOut, double **YOut)
{
   int    ii, ss, jj, ll, index, totPts;
   double *XT, *XX, *YY, *HX;
   psVector vecXT, vecHX;

   initialize(XIn,YIn);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   vecHX.setLength(3);
   vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_-1); 
   vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_-1); 
   vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_-1); 

   (*XOut) = new double[3*totPts];
   (*YOut) = new double[totPts];
   (*NOut) = totPts;

   vecXT.setLength(totPts*nInputs_);
   for (ss = 0; ss < totPts; ss++)
      for (ii = 0; ii < nInputs_; ii++) 
        vecXT[ss*nInputs_+ii] = settings[ii];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            vecXT[index*nInputs_+ind1] = vecHX[0] * ii + lowerBounds_[ind1];
            vecXT[index*nInputs_+ind2] = vecHX[1] * jj + lowerBounds_[ind2];
            vecXT[index*nInputs_+ind3] = vecHX[2] * ll + lowerBounds_[ind3];
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
int MRBF::gen4DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                        int ind3, int ind4, double *settings, int *NOut, 
                        double **XOut, double **YOut)
{
  int    ii, ss, jj, ll, mm, index, totPts;
  psVector vecXT, vecHX;

  initialize(XIn,YIn);

  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_-1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_-1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_-1); 
  vecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_-1); 

  (*XOut) = new double[4*totPts];
  (*YOut) = new double[totPts];
  (*NOut) = totPts;

  vecXT.setLength(totPts*nInputs_);
  for (ss = 0; ss < totPts; ss++) 
    for (ii = 0; ii < nInputs_; ii++) 
      vecXT[ss*nInputs_+ii] = settings[ii]; 
    
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
          vecXT[index*nInputs_+ind1] = vecHX[0]*ii+lowerBounds_[ind1];
          vecXT[index*nInputs_+ind2] = vecHX[1]*jj+lowerBounds_[ind2];
          vecXT[index*nInputs_+ind3] = vecHX[2]*ll+lowerBounds_[ind3];
          vecXT[index*nInputs_+ind4] = vecHX[3]*mm+lowerBounds_[ind4];
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
double MRBF::evaluatePoint(double *X)
{
   int    iOne=1;
   double Y;
   evaluatePoint(1, X, &Y);
   return Y;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double MRBF::evaluatePoint(int nPts, double *X, double *Y)
{
  int    ss, ss2, ii, pp, nSamp, count, highFlag;
  double dist, diff, Yt, ddata, *lbs, *ubs, *XP;
  psVector vecRanges;

  vecRanges.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
     vecRanges[ii] = 1.0 / (upperBounds_[ii] - lowerBounds_[ii]);
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
        nSamp = boxes_[pp]->nSamples_;
        XP = boxes_[pp]->vecX_.getDVector();
        for (ss2 = 0; ss2 < nSamp; ss2++) 
        {
          dist = 0.0;
          for (ii = 0; ii < nInputs_; ii++) 
          {
            ddata = X[ss*nInputs_+ii];
            ddata = (ddata - lowerBounds_[ii]) * vecRanges[ii];
            ddata -= XP[ss2*nInputs_+ii];
            dist += ddata * ddata;
          }
          switch (type_)
          {
            case 0: dist = sqrt(dist + 1.0); break;
            case 1: dist = 1.0/sqrt(dist + 1.0); break;
            case 2: dist = exp(-0.5*dist*gaussScale_); break;
            case 3: dist = (dist+1)*log(sqrt(dist+1)); break;
          }
          Yt += dist * boxes_[pp]->vecCoeffs_[ss2];
        }
        Yt += boxes_[pp]->vecCoeffs_[nSamp];
        count++;
      } 
    }
    if (count == 0)
    {
      printf("MRBF evaluate ERROR: sample point outside range.\n");
      printf("INFO: this may happen during cross validation.\n"); 
      for (ii = 0; ii < nInputs_; ii++)
        printf("Input %d = %e (in [%e, %e]?)\n",ii+1,
           X[ss*nInputs_+ii],lowerBounds_[ii],upperBounds_[ii]);
      printf("Sample prediction will be set to Ymean: may be wrong.\n");
      Yt = 0;
      count = 1;
    }
    Yt /= (double) count;
    Y[ss] = Yt + YMean_;
  }
  return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double MRBF::evaluatePointFuzzy(double *X, double &std)
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
double MRBF::evaluatePointFuzzy(int npts, double *X, double *Y, double *Ystd)
{
   evaluatePoint(npts, X, Y);
   for (int ss = 0; ss < npts; ss++) Ystd[ss] = 0.0;
   return 0.0;
}

