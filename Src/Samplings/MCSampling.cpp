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
// Functions for the Monte Carlo sampling class
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "MCSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
MCSampling::MCSampling() : Sampling()
{
  samplingID_ = PSUADE_SAMP_MC;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
MCSampling::~MCSampling()
{
}

// ************************************************************************
// initialization 
// ------------------------------------------------------------------------
int MCSampling::initialize(int initLevel)
{
  int    ii, kk;
  double range;

  if (nSamples_ == 0)
  {
    printf("MCSampling::initialize ERROR - nSamples = 0.\n");
    exit(1);
  }
  if (nInputs_ == 0)
  {
    printf("MCSampling::initialize ERROR - input not set up.\n");
    exit(1);
  }

  if (printLevel_ > 4)
  {
    printf("MCSampling::initialize: nSamples = %d\n", nSamples_);
    printf("MCSampling::initialize: nInputs  = %d\n", nInputs_);
    printf("MCSampling::initialize: nOutputs = %d\n", nOutputs_);
    for (ii = 0; ii < nInputs_; ii++)
      printf("    MCSampling input %3d = [%e %e]\n", ii+1,
             vecLBs_[ii], vecUBs_[ii]);
  }

  if (initLevel != 0) return 0;

  allocSampleData();
  for (ii = 0; ii < nInputs_; ii++)
  { 
    range = vecUBs_[ii] - vecLBs_[ii];
    for (kk = 0; kk < nSamples_; kk++) 
      vecSamInps_[kk*nInputs_+ii] = PSUADE_drand() * range + vecLBs_[ii];
  }
  return 0;
}

// ************************************************************************
// refine 
// ------------------------------------------------------------------------
int MCSampling::refine(int refineRatio, int randomize, double thresh,
                       int nSamples, double *sampleErrors)
{
  int    ii, inputID, oldNumSamples, nLevels;
  double range;
  psVector  vecSamInpsOld, vecSamOutsOld;
  psIVector vecSamStasOld;

  (void) randomize;
  (void) thresh;
  (void) nSamples;
  (void) sampleErrors;

  nLevels       = refineRatio;
  oldNumSamples = nSamples_;
  vecSamInpsOld = vecSamInps_;
  vecSamOutsOld = vecSamOuts_;
  vecSamStasOld = vecSamStas_;

  nSamples_ *= nLevels;
  vecSamInps_.setLength(nSamples_*nInputs_);
  vecSamOuts_.setLength(nSamples_*nOutputs_);
  vecSamStas_.setLength(nSamples_);

  if (printLevel_ > 4)
  {
    printf("MCSampling::refine: nSamples = %d\n", nSamples_);
    printf("MCSampling::refine: nInputs  = %d\n", nInputs_);
    printf("MCSampling::refine: nOutputs = %d\n", nOutputs_);
    for (inputID = 0; inputID < nInputs_; inputID++)
      printf("    MCSampling input %3d = [%e %e]\n", inputID+1,
             vecLBs_[inputID], vecUBs_[inputID]);
  }

  for (ii = 0; ii < oldNumSamples*nInputs_; ii++) 
    vecSamInps_[ii] = vecSamInpsOld[ii];
  for (ii = 0; ii < oldNumSamples*nOutputs_; ii++) 
    vecSamOuts_[ii] = vecSamOutsOld[ii];
  for (ii = 0; ii < oldNumSamples; ii++) 
    vecSamStas_[ii] = vecSamStasOld[ii];

  for (inputID = 0; inputID < nInputs_; inputID++) 
  {
    range = vecUBs_[inputID] - vecLBs_[inputID];
    for (ii = oldNumSamples; ii < nSamples_; ii++)
      vecSamInps_[ii*nInputs_+inputID] = PSUADE_drand() * range + 
                                         vecLBs_[inputID];
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MCSampling& MCSampling::operator=(const MCSampling &)
{
  printf("MCSampling operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

