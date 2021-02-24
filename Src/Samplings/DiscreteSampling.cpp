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
// Functions for the discrete sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "DiscreteSampling.h"

#define PABS(x) ((x) > 0 ? x : -(x))
// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
DiscreteSampling::DiscreteSampling() : Sampling()
{
  samplingID_   = PSUADE_SAMP_DISCRETE;
}

// ************************************************************************
// copy constructor by Bill Oliver
// ------------------------------------------------------------------------
DiscreteSampling::DiscreteSampling(const DiscreteSampling & ds) : Sampling()
{
  printLevel_ = ds.printLevel_;
  samplingID_ = ds.samplingID_;
  nSamples_ = ds.nSamples_;
  nInputs_ = ds.nInputs_;
  nOutputs_ = ds.nOutputs_;
  randomize_ = ds.randomize_;
  nReplications_ = ds.nReplications_;
  vecLBs_ = ds.vecLBs_;
  vecUBs_ = ds.vecUBs_;
  vecSamInps_ = ds.vecSamInps_;
  vecSamOuts_ = ds.vecSamOuts_;
  vecSamStas_ = ds.vecSamStas_;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
DiscreteSampling::~DiscreteSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int DiscreteSampling::initialize(int initLevel)
{
  int    inLeng, ii, jj, kk, ll, sampleID, idata, cnt, num;
  double dsum, ddata;
  char   lineIn[1001], filename[1001], cword[1001];
  FILE   *fp=NULL;

  maxLevels_ = 100;

  printf("To use discrete sampling, you will have to provide a\n");
  printf("parameter file. The parameter file should be of the\n");
  printf("following format:\n");
  printf("PSUADE_BEGIN\n");
  printf("<number of inputs> <number of sample points>\n");
  printf("1 numLevels(n) L1 P1 ... Ln Pn <for input 1>\n");
  printf("2 numLevels(n) L1 P1 ... Ln Pn <for input 2>\n");
  printf("...\n");
  printf("PSUADE_END\n");
  printf("Note: numLevels - number of discrete values for the input.\n");
  printf("      L1        - value of the first level\n");
  printf("      P1        - probability of the first level\n");
  printf("      Ln        - value of the last level\n");
  printf("      Pn        - probability of the last level\n");
  printf("Note: the sum of P's for an input should be 1.\n");
  printf("Enter the name of the parameter file: ");
  scanf("%s", filename);
  fgets(lineIn, 500, stdin);
  inLeng = strlen(filename);
  if (inLeng < 500)
  {
    filename[inLeng] = '\0';
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
      printf("DiscreteSampling ERROR: cannot open file %s.\n",filename);
      return -1;
    }
  }
  else
  {
    printf("DiscreteSampling ERROR: file name too long.\n");
    return -1;
  }
  fgets(lineIn, 1000, stdin);
  sscanf(lineIn, "%s", cword);
  if (!strcmp(cword, "PSUADE_BEGIN"))
  {
    scanf("%d %d", &nInputs_, &nSamples_);
    if (nInputs_ <= 0)
    {
      printf("DiscreteSampling ERROR : nInputs <= 0.\n");
      fclose(fp);
      return -1;
    }
    if (nSamples_ <= 0)
    {
      printf("DiscreteSampling ERROR : nSamples <= 0.\n");
      fclose(fp);
      return -1;
    }
    vecInpValCnts_.setLength(nInputs_);
    vecInpVals_.setLength(nInputs_*maxLevels_);
    vecInpProbs_.setLength(nInputs_*maxLevels_);
    for (ii = 0; ii < nInputs_; ii++)
    {
      scanf("%d", &kk);
      if (kk != ii+1)
      {
        printf("DiscreteSampling ERROR: invalid input index %d (!= %d).\n",
               kk,ii+1);
        fclose(fp);
        return -1;
      }
      scanf("%d", &num);
      if (num <= 0 || num > maxLevels_)
      {
        printf("DiscreteSampling ERROR: invalid numLevels %d (input %d).\n",
               num,ii+1);
        if (num > maxLevels_) 
          printf("   Maximum number of levels = %d\n",maxLevels_);
        fclose(fp);
        return -1;
      }
      vecInpValCnts_[ii] = num;
      dsum = 0.0;
      if (printLevel_ > 1) printf("Input %4d: \n", ii+1);
      for (jj = 0; jj < num; jj++)
      {
        scanf("%d %lg",&idata, &ddata);
        vecInpVals_[ii*maxLevels_+jj] = idata;
        vecInpProbs_[ii*maxLevels_+jj] = ddata;
        if (ddata <= 0.0)
        {
          fclose(fp);
          printf("DiscreteSampling ERROR: probability should be > 0.\n");
          return -1;
        }
        if (printLevel_ > 1)
          printf("     Level = %8d, Probability = %e\n", 
             vecInpVals_[ii*maxLevels_+jj],vecInpProbs_[ii*maxLevels_+jj]);
        dsum += vecInpProbs_[ii*maxLevels_+jj];
      } 
      if (dsum != 1.0)
      {
        fclose(fp);
        printf("DiscreteSampling ERROR: sum of probability should be 1.\n");
        return -1;
      }
    }
  }
  else
  {
    printf("DiscreteSampling ERROR: PSUADE_BEGIN not found.\n");
    fclose(fp);
    return -1;
  }
  fgets(lineIn, 1000, stdin);
  fgets(lineIn, 1000, stdin);
  sscanf(lineIn, "%s", cword);
  if (strcmp(cword, "PSUADE_END"))
  {
    printf("DiscreteSampling ERROR: PSUADE_END not found.\n");
    fclose(fp);
    return -1;
  }
  fclose(fp);

  vecSamInps_.setLength(nSamples_*nInputs_);
  vecSamOuts_.setLength(nSamples_*nOutputs_);
  vecSamStas_.setLength(nSamples_);
  psIVector vecPT;
  vecPT.setLength(nInputs_*101);
  for (ii = 0; ii < nInputs_; ii++)
  {
    cnt = 0;
    kk = vecInpValCnts_[ii];
    for (jj = 0; jj < kk; jj++)
    {
      num = (int) ((vecInpProbs_[ii*maxLevels_+jj] + 1.0e-12) * 100);
      for (ll = 0; ll < num; ll++)
        vecPT[ii*101+cnt+ll] = vecInpVals_[ii*maxLevels_+jj];
      cnt += num;
    }
  }
  for (sampleID = 0; sampleID < nSamples_; sampleID++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      kk = PSUADE_rand() % 100;
      vecSamInps_[sampleID*nInputs_+ii] = (double) vecPT[ii*101+kk];
    }
  }

  if (printLevel_ > 4)
  {
    printf("DiscreteSampling::initialize: nSamples = %d\n", nSamples_);
    printf("DiscreteSampling::initialize: nInputs  = %d\n", nInputs_);
  }
  return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int DiscreteSampling::refine(int refineRatio,int randomize,double thresh,
                              int nSamples, double *sampleErrors)
{
  (void) refineRatio;
  (void) randomize;
  (void) thresh;
  (void) nSamples;
  (void) sampleErrors;
  printf("DiscreteSampling refine: not implemented yet.\n");
  return -1;
} 

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
DiscreteSampling& DiscreteSampling::operator=(const DiscreteSampling & ds)
{
  if(this == & ds) return *this;
  printLevel_ = ds.printLevel_;
  samplingID_ = ds.samplingID_;
  nSamples_ = ds.nSamples_;
  nInputs_ = ds.nInputs_;
  nOutputs_ = ds.nOutputs_;
  randomize_ = ds.randomize_;
  nReplications_ = ds.nReplications_;
  vecLBs_ = ds.vecLBs_;
  vecUBs_ = ds.vecUBs_;
  vecSamInps_ = ds.vecSamInps_;
  vecSamOuts_ = ds.vecSamOuts_;
  vecSamStas_ = ds.vecSamStas_;
  return *this;
}

