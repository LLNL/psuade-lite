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
// Functions for the Box-Behnken sampling class 
// This sample design is compatible with quadratic regression
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>

#include "sysdef.h"
#include "PsuadeUtil.h"
#include "BoxBehnkenSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
BoxBehnkenSampling::BoxBehnkenSampling() : Sampling()
{
  samplingID_ = PSUADE_SAMP_BBD;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BoxBehnkenSampling::~BoxBehnkenSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int BoxBehnkenSampling::initialize(int initLevel)
{
   int inputID, sampleID, inputID2;

  if (nSamples_ == 0)
  {
     printf("BoxBehnkenSampling::initialize ERROR - nSamples = 0.\n");
     exit(1);
  }
  if (nInputs_ == 0)
  {
    printf("BoxBehnkenSampling::initialize ERROR - input not set up.\n");
    exit(1);
  }
  if (nInputs_ < 2)
  {
    printf("BoxBehnkenSampling::initialize ERROR - nInputs < 2.\n");
    printf("                    nInputs should be [2,7].\n");
    exit(1);
  }
  if (nInputs_ > 7)
  {
    printf("BoxBehnkenSampling::initialize ERROR - nInputs > 7.\n");
    printf("                    nInputs should be [2,7].\n");
    exit(1);
  }
   
  if (initLevel == 1) return 0;
  if      (nInputs_ == 2) nSamples_ = 5;
  else if (nInputs_ == 3) nSamples_ = 13;
  else if (nInputs_ == 4) nSamples_ = 25;
  else if (nInputs_ == 5) nSamples_ = 41;
  else if (nInputs_ == 6) nSamples_ = 48;
  else if (nInputs_ == 7) nSamples_ = 57;
  allocSampleData();

  if (printLevel_ > 4)
  {
    printf("BoxBehnkenSampling::initialize: nSamples = %d\n", nSamples_);
    printf("BoxBehnkenSampling::initialize: nInputs  = %d\n", nInputs_);
    printf("BoxBehnkenSampling::initialize: nOutputs = %d\n", nOutputs_);
    for (inputID = 0; inputID < nInputs_; inputID++)
      printf("    BoxBehnkenSampling input %3d = [%e %e]\n", inputID+1,
             vecLBs_[inputID], vecUBs_[inputID]);
  }

  for (sampleID = 0; sampleID < nSamples_; sampleID++)
    for (inputID = 0; inputID < nInputs_; inputID++)
      vecSamInps_[sampleID*nInputs_+inputID] = 
                0.5 * (vecLBs_[inputID] + vecUBs_[inputID]);
  if (nInputs_ >= 2 && nInputs_ <= 5)
  {
    sampleID = 0;
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      for (inputID2 = inputID+1; inputID2 < nInputs_; inputID2++)
      {
        vecSamInps_[sampleID*nInputs_+inputID]  = vecLBs_[inputID];
        vecSamInps_[sampleID*nInputs_+inputID2] = vecLBs_[inputID2];
        vecSamInps_[(sampleID+1)*nInputs_+inputID]  = vecLBs_[inputID];
        vecSamInps_[(sampleID+1)*nInputs_+inputID2] = vecUBs_[inputID2];
        vecSamInps_[(sampleID+2)*nInputs_+inputID]  = vecUBs_[inputID];
        vecSamInps_[(sampleID+2)*nInputs_+inputID2] = vecLBs_[inputID2];
        vecSamInps_[(sampleID+3)*nInputs_+inputID]  = vecUBs_[inputID];
        vecSamInps_[(sampleID+3)*nInputs_+inputID2] = vecUBs_[inputID2];
        sampleID += 4;
      }
    }
  }
  if (nInputs_ == 6)
  {
    sampleID = 0;
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      vecSamInps_[sampleID*nInputs_+inputID]      = vecLBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[sampleID*nInputs_+inputID2]     = vecLBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[sampleID*nInputs_+inputID2]     = vecLBs_[inputID2];
      vecSamInps_[(sampleID+1)*nInputs_+inputID]  = vecLBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[(sampleID+1)*nInputs_+inputID2] = vecLBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[(sampleID+1)*nInputs_+inputID2] = vecUBs_[inputID2];
      vecSamInps_[(sampleID+2)*nInputs_+inputID]  = vecLBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[(sampleID+2)*nInputs_+inputID2] = vecUBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[(sampleID+2)*nInputs_+inputID2] = vecLBs_[inputID2];
      vecSamInps_[(sampleID+3)*nInputs_+inputID]  = vecLBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[(sampleID+3)*nInputs_+inputID2] = vecUBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[(sampleID+3)*nInputs_+inputID2] = vecUBs_[inputID2];

      vecSamInps_[(sampleID+4)*nInputs_+inputID]  = vecUBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[(sampleID+4)*nInputs_+inputID2] = vecLBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[(sampleID+4)*nInputs_+inputID2] = vecLBs_[inputID2];
      vecSamInps_[(sampleID+5)*nInputs_+inputID]  = vecUBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[(sampleID+5)*nInputs_+inputID2] = vecLBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[(sampleID+5)*nInputs_+inputID2] = vecUBs_[inputID2];
      vecSamInps_[(sampleID+6)*nInputs_+inputID]  = vecUBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[(sampleID+6)*nInputs_+inputID2] = vecUBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[(sampleID+6)*nInputs_+inputID2] = vecLBs_[inputID2];
      vecSamInps_[(sampleID+7)*nInputs_+inputID]  = vecUBs_[inputID];
      inputID2 = (inputID + 1) % nInputs_;
      vecSamInps_[(sampleID+7)*nInputs_+inputID2] = vecUBs_[inputID2];
      inputID2 = (inputID + 3) % nInputs_;
      vecSamInps_[(sampleID+7)*nInputs_+inputID2] = vecUBs_[inputID2];
      sampleID += 8;
    }
  }
  if (nInputs_ == 7)
  {
    sampleID = 0;
    vecSamInps_[(sampleID)*nInputs_+3]   = vecLBs_[3];
    vecSamInps_[(sampleID)*nInputs_+4]   = vecLBs_[4];
    vecSamInps_[(sampleID)*nInputs_+5]   = vecLBs_[5];
    vecSamInps_[(sampleID+1)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+1)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+1)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+2)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+2)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+2)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+3)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+3)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+3)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+4)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+4)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+4)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+5)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+5)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+5)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+6)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+6)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+6)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+7)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+7)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+7)*nInputs_+5] = vecUBs_[5];
    sampleID += 8;
    vecSamInps_[(sampleID)*nInputs_+0]   = vecLBs_[0];
    vecSamInps_[(sampleID)*nInputs_+5]   = vecLBs_[5];
    vecSamInps_[(sampleID)*nInputs_+6]   = vecLBs_[6];
    vecSamInps_[(sampleID+1)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+1)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+1)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+2)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+2)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+2)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+3)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+3)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+3)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+4)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+4)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+4)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+5)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+5)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+5)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+6)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+6)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+6)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+7)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+7)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+7)*nInputs_+6] = vecUBs_[6];
    sampleID += 8;
    vecSamInps_[(sampleID)*nInputs_+1]   = vecLBs_[1];
    vecSamInps_[(sampleID)*nInputs_+4]   = vecLBs_[4];
    vecSamInps_[(sampleID)*nInputs_+6]   = vecLBs_[6];
    vecSamInps_[(sampleID+1)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+1)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+1)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+2)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+2)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+2)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+3)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+3)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+3)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+4)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+4)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+4)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+5)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+5)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+5)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+6)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+6)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+6)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+7)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+7)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+7)*nInputs_+6] = vecUBs_[6];
    sampleID += 8;
    vecSamInps_[(sampleID)*nInputs_+0]   = vecLBs_[0];
    vecSamInps_[(sampleID)*nInputs_+1]   = vecLBs_[1];
    vecSamInps_[(sampleID)*nInputs_+3]   = vecLBs_[3];
    vecSamInps_[(sampleID+1)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+1)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+1)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+2)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+2)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+2)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+3)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+3)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+3)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+4)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+4)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+4)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+5)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+5)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+5)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+6)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+6)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+6)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+7)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+7)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+7)*nInputs_+3] = vecUBs_[3];
    sampleID += 8;
    vecSamInps_[(sampleID)*nInputs_+2]   = vecLBs_[2];
    vecSamInps_[(sampleID)*nInputs_+3]   = vecLBs_[3];
    vecSamInps_[(sampleID)*nInputs_+6]   = vecLBs_[6];
    vecSamInps_[(sampleID+1)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+1)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+1)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+2)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+2)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+2)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+3)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+3)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+3)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+4)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+4)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+4)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+5)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+5)*nInputs_+3] = vecLBs_[3];
    vecSamInps_[(sampleID+5)*nInputs_+6] = vecUBs_[6];
    vecSamInps_[(sampleID+6)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+6)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+6)*nInputs_+6] = vecLBs_[6];
    vecSamInps_[(sampleID+7)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+7)*nInputs_+3] = vecUBs_[3];
    vecSamInps_[(sampleID+7)*nInputs_+6] = vecUBs_[6];
    sampleID += 8;
    vecSamInps_[(sampleID)*nInputs_+0]   = vecLBs_[0];
    vecSamInps_[(sampleID)*nInputs_+2]   = vecLBs_[2];
    vecSamInps_[(sampleID)*nInputs_+4]   = vecLBs_[4];
    vecSamInps_[(sampleID+1)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+1)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+1)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+2)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+2)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+2)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+3)*nInputs_+0] = vecLBs_[0];
    vecSamInps_[(sampleID+3)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+3)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+4)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+4)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+4)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+5)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+5)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+5)*nInputs_+4] = vecUBs_[4];
    vecSamInps_[(sampleID+6)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+6)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+6)*nInputs_+4] = vecLBs_[4];
    vecSamInps_[(sampleID+7)*nInputs_+0] = vecUBs_[0];
    vecSamInps_[(sampleID+7)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+7)*nInputs_+4] = vecUBs_[4];
    sampleID += 8;
    vecSamInps_[(sampleID)*nInputs_+1]   = vecLBs_[1];
    vecSamInps_[(sampleID)*nInputs_+2]   = vecLBs_[2];
    vecSamInps_[(sampleID)*nInputs_+5]   = vecLBs_[5];
    vecSamInps_[(sampleID+1)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+1)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+1)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+2)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+2)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+2)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+3)*nInputs_+1] = vecLBs_[1];
    vecSamInps_[(sampleID+3)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+3)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+4)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+4)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+4)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+5)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+5)*nInputs_+2] = vecLBs_[2];
    vecSamInps_[(sampleID+5)*nInputs_+5] = vecUBs_[5];
    vecSamInps_[(sampleID+6)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+6)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+6)*nInputs_+5] = vecLBs_[5];
    vecSamInps_[(sampleID+7)*nInputs_+1] = vecUBs_[1];
    vecSamInps_[(sampleID+7)*nInputs_+2] = vecUBs_[2];
    vecSamInps_[(sampleID+7)*nInputs_+5] = vecUBs_[5];
  }
  return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int BoxBehnkenSampling::refine(int, int, double, int, double *)
{
  printf("BoxBehnkenSampling::refine ERROR - not available.\n");
  exit(1);
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
BoxBehnkenSampling& BoxBehnkenSampling::operator=(const BoxBehnkenSampling &)
{
  printf("BoxBehnkenSampling operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

