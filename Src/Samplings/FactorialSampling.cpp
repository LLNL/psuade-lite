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
// Functions for the factorial sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdio.h>

#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FactorialSampling.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FactorialSampling::FactorialSampling() : Sampling()
{
  samplingID_ = PSUADE_SAMP_FACT;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
FactorialSampling::~FactorialSampling()
{
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int FactorialSampling::initialize(int initLevel)
{
  int    inputID, sampleID, curSymbol, nsym;
  double scale, ddata, dpower;

  if (nSamples_ == 0)
  {
    printf("FactorialSampling::initialize ERROR - nSamples = 0.\n");
    exit(1);
  }
  if (nInputs_ == 0)
  {
    printf("FactorialSampling::initialize ERROR - input not set up.\n");
    exit(1);
  }

  if (vecInpSettings_ != NULL)
  {
    if (printLevel_ > 4)
       printf("FactorialSampling::initialize: inputSettings used.\n");
    if (printLevel_ > 4 && vecSymTable_.length() > 0)
       printf("FactorialSampling::initialize: symbol table overwritten.\n");
    maxNumSettings_ = 0;
    for (inputID = 0; inputID < nInputs_; inputID++) 
      if (vecInpNumSettings_[inputID] > maxNumSettings_)
        maxNumSettings_ = vecInpNumSettings_[inputID];
    
    vecSymTable_.setLength(nInputs_);
    for (inputID = 0; inputID < nInputs_; inputID++) 
      vecSymTable_[inputID] = vecInpNumSettings_[inputID];
    nSamples_ = 1;
    for (inputID = 0; inputID < nInputs_; inputID++) 
      nSamples_ *= vecSymTable_[inputID];
    if (nSamples_ == 0)
    {
      printf("FactorialSampling::initialize ERROR - ");
      printf("incomplete inputSettings.\n");
      exit(1);
    }
  } 
  else if (vecSymTable_.length() != nInputs_)
  {
    ddata  = (double) nSamples_;
    dpower = 1.0 / (double) nInputs_;
    ddata  = pow(ddata, dpower+1.0E-12);
    nsym   = (int) ddata;
    vecSymTable_.setLength(nInputs_);
    for (inputID = 0; inputID < nInputs_; inputID++) 
      vecSymTable_[inputID] = nsym;
    nSamples_ = 1;
    for (inputID = 0; inputID < nInputs_; inputID++) 
      nSamples_ *= vecSymTable_[inputID];
  }
  if (initLevel == 1) return 0;

  if (printLevel_ > 4)
  {
    printf("FactorialSampling::initialize: nSamples = %d\n", nSamples_);
    printf("FactorialSampling::initialize: nInputs  = %d\n", nInputs_);
    printf("FactorialSampling::initialize: nOutputs = %d\n", nOutputs_);
    if (randomize_ != 0)
         printf("FactorialSampling::initialize: randomize on\n");
    else printf("FactorialSampling::initialize: randomize off\n");
    for (inputID = 0; inputID < nInputs_; inputID++)
      printf("    FactorialSampling input %3d = [%e %e]\n", inputID+1,
             vecLBs_[inputID], vecUBs_[inputID]);
    if (vecSymTable_.length() == nInputs_)
      for (inputID = 0; inputID < nInputs_; inputID++) 
            printf("    FactorialSampling symbol table %2d = %d\n",
                   inputID+1, vecSymTable_[inputID]);
  }

  allocSampleData();
  if (nSamples_ == 1)
  {
    for (inputID = 0; inputID < nInputs_; inputID++) 
    {
      vecSamInps_[0*nInputs_+inputID] = 0.5 * (vecLBs_[inputID] +
                                               vecUBs_[inputID]);
    }
    return 0;
  }
  psVector  vecRanges;
  psIVector vecCntMax, vecCnts, vecSyms;
  vecRanges.setLength(nInputs_);
  for (inputID = 0; inputID < nInputs_; inputID++)
    vecRanges[inputID] = vecUBs_[inputID] - vecLBs_[inputID];
  vecCntMax.setLength(nInputs_);
  vecCnts.setLength(nInputs_);
  vecSyms.setLength(nInputs_);
  vecCntMax[0] = 1;
  for (inputID = 1; inputID < nInputs_; inputID++)
    vecCntMax[inputID] = vecCntMax[inputID-1] * vecSymTable_[inputID-1];
  for (inputID = 0; inputID < nInputs_; inputID++)
    vecCnts[inputID] = vecSyms[inputID] = 0;

  if (randomize_ != 0)
  {
    for (sampleID = 0; sampleID < nSamples_; sampleID++)
    {
      for (inputID = 0; inputID < nInputs_; inputID++) 
      {
        curSymbol = vecSyms[inputID];
        vecCnts[inputID]++;
        if (vecCnts[inputID] == vecCntMax[inputID])
        {
          vecCnts[inputID] = 0;
          vecSyms[inputID]++;
          if (vecSyms[inputID] >= vecSymTable_[inputID]) 
            vecSyms[inputID] = 0;
        }
        scale = vecRanges[inputID] / ((double) (vecSymTable_[inputID]));
        vecSamInps_[sampleID*nInputs_+inputID] = 
                        curSymbol * scale + vecLBs_[inputID];
        vecSamInps_[sampleID*nInputs_+inputID] += (PSUADE_drand() * scale);
      }
    }
  }
  else
  {
    for (sampleID = 0; sampleID < nSamples_; sampleID++)
    {
      for (inputID = 0; inputID < nInputs_; inputID++) 
      {
        curSymbol = vecSyms[inputID];
        vecCnts[inputID]++;
        if (vecCnts[inputID] == vecCntMax[inputID])
        {
          vecCnts[inputID] = 0;
          vecSyms[inputID]++;
          if (vecSyms[inputID] >= vecSymTable_[inputID]) 
            vecSyms[inputID] = 0;
        }
        scale = vecRanges[inputID] / ((double) (vecSymTable_[inputID] - 1));
        if (vecInpNumSettings_.length() > 0)
           vecSamInps_[sampleID*nInputs_+inputID] =
                         vecInpSettings_[inputID][curSymbol]; 
        else
           vecSamInps_[sampleID*nInputs_+inputID] = 
                         curSymbol * scale + vecLBs_[inputID];
      }
    }
  }
  return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int FactorialSampling::refine(int refineRatio, int randomize, double thresh,
                              int nSamples, double *sampleErrors)
{
  int    newNSamples, inputID, *useArray, *iCounts, *iCntMax, *iSymbol;
  int    sampleID, index, curSymbol, *newSampleStates, outputID, nLevels;
  double **newSampleMatrix, stepSize, ddata;
  double *newSampleOutput;

  (void) randomize;
  (void) thresh;
  (void) nSamples;
  (void) sampleErrors;

  nLevels = refineRatio;
  if (nLevels != 2)
  {
    printf("FactorialSampling::refine WARNING - nLevels set to 2.\n");
    nLevels = 2;
  }

  newNSamples = 1;
  if (randomize_ != 0)
  {
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      vecSymTable_[inputID] *= 2;
      newNSamples *= vecSymTable_[inputID];
    }
  }
  else
  {
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      vecSymTable_[inputID] = vecSymTable_[inputID] * 2 - 1;
      newNSamples *= vecSymTable_[inputID];
    }
  }

  psVector vecRanges;
  vecRanges.setLength(nInputs_);
  for (inputID = 0; inputID < nInputs_; inputID++)
    vecRanges[inputID] = vecUBs_[inputID] - vecLBs_[inputID];
  psIVector vecUseArray;
  vecUseArray.setLength(newNSamples);
  for (sampleID = 0; sampleID < newNSamples; sampleID++)
    vecUseArray[sampleID] = 0;
  psIVector vecCntMax, vecCnts, vecSyms;
  vecCntMax.setLength(nInputs_);
  vecCnts.setLength(nInputs_);
  vecSyms.setLength(nInputs_);
  vecCntMax[0] = 1;
  for (inputID = 1; inputID < nInputs_; inputID++)
    vecCntMax[inputID] = vecCntMax[inputID-1] * vecSymTable_[inputID-1];
  for (inputID = 0; inputID < nInputs_; inputID++)
    vecCnts[inputID] = vecSyms[inputID] = 0;

  for (sampleID = 0; sampleID < nSamples_; sampleID++)
  {
    index = 0;
    for (inputID = 0; inputID < nInputs_; inputID++) 
    {
      if (randomize_ != 0)
         stepSize = ((double) (vecSymTable_[inputID]))/vecRanges[inputID];
      else
         stepSize = ((double) (vecSymTable_[inputID]-1))/vecRanges[inputID];
      ddata = vecSamInps_[sampleID*nInputs_+inputID] - vecLBs_[inputID];
      curSymbol = (int) (ddata * stepSize * 1.00000000001);
      index += curSymbol * vecCntMax[inputID];
    } 
    vecUseArray[index] = sampleID + 1;
    //printf("Fact: use %d = %d\n", index, useArray[index]);
  }

  psVector  vecSamInps2, vecSamOuts2;
  psIVector vecSamStas2;
  
  vecSamInps2 = vecSamInps_;
  vecSamOuts2 = vecSamOuts2;
  vecSamStas2 = vecSamStas2;
  vecSamInps_.setLength(newNSamples*nInputs_);
  vecSamOuts_.setLength(newNSamples*nOutputs_);
  vecSamStas_.setLength(newNSamples);

  for (sampleID = 0; sampleID < newNSamples*nOutputs_; sampleID++)
    vecSamOuts_[sampleID] = PSUADE_UNDEFINED;
  for (sampleID = 0; sampleID < newNSamples; sampleID++)
    vecSamStas_[sampleID] = 0;
   
  for (sampleID = 0; sampleID < newNSamples; sampleID++)
  {
    index = -1;
    if (vecUseArray[sampleID] > 0) index = vecUseArray[sampleID] - 1;

    if (index >= 0)
    {
      for (inputID = 0; inputID < nInputs_; inputID++) 
        vecSamInps_[sampleID*nInputs_+inputID] = 
                        vecSamInps2[index*nInputs_+inputID]; 
      for (outputID = 0; outputID < nOutputs_; outputID++) 
        vecSamOuts_[sampleID*nOutputs_+outputID] = 
                        vecSamOuts2[index*nOutputs_+outputID]; 
      vecSamStas_[sampleID] = vecSamStas2[index];
    }
    for (inputID = 0; inputID < nInputs_; inputID++) 
    {
      curSymbol = vecSyms[inputID];
      vecCnts[inputID]++;
      if (vecCnts[inputID] == vecCntMax[inputID])
      {
        vecCnts[inputID] = 0;
        vecSyms[inputID]++;
        if (vecSyms[inputID] >= vecSymTable_[inputID]) 
          vecSyms[inputID] = 0;
      }
      if (index < 0)
      {
        if (randomize_ != 0)
        {
          stepSize = vecRanges[inputID] / ((double) (vecSymTable_[inputID]));
          vecSamInps_[sampleID*nInputs_+inputID] = vecLBs_[inputID] +
                                (curSymbol + PSUADE_drand()) * stepSize;
        }
        else
        {
          stepSize = vecRanges[inputID] / ((double) (vecSymTable_[inputID]-1));
          vecSamInps_[sampleID*nInputs_+inputID] = 
                           curSymbol * stepSize + vecLBs_[inputID];
        }
      }
    }
  }
  nSamples_ = newNSamples;

  if (printLevel_ > 4)
  {
    printf("FactorialSampling::refine: nSamples = %d\n", nSamples_);
    printf("FactorialSampling::refine: nInputs  = %d\n", nInputs_);
    printf("FactorialSampling::refine: nOutputs = %d\n", nOutputs_);
    if (randomize_ != 0)
         printf("FactorialSampling::refine: randomize on\n");
    else printf("FactorialSampling::refine: randomize off\n");
    for (inputID = 0; inputID < nInputs_; inputID++)
      printf("    FactorialSampling input %3d = [%e %e]\n", inputID+1,
             vecLBs_[inputID], vecUBs_[inputID]);
    if (vecSymTable_.length() == nInputs_)
      for (inputID = 0; inputID < nInputs_; inputID++) 
        printf("    FactorialSampling symbol table %2d = %d\n",
               inputID+1, vecSymTable_[inputID]);
  }
  return 0;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int FactorialSampling::setInputParams(int nInputs, int *counts, 
                                      double **settings, int *symtable)
{
  int inputID, sID;

  if (nInputs_ != 0 && nInputs != nInputs_)
  {
    printf("FactorialSampling::setInputParams - nInputs mismatch.\n");
    exit(1);
  }
  nInputs_ = nInputs;
  if (symtable != NULL)
  {
    vecSymTable_.setLength(nInputs);
    for (inputID = 0; inputID < nInputs_; inputID++)
      vecSymTable_[inputID] = symtable[inputID];
  }
  if (counts != NULL)
  {
    maxNumSettings_ = 0;
    for (inputID = 0; inputID < nInputs_; inputID++) 
      if (counts[inputID] > maxNumSettings_)
        maxNumSettings_ = counts[inputID];
    vecInpNumSettings_.setLength(nInputs_);
    vecInpSettings_ = new psVector[nInputs_];
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      vecInpSettings_[inputID].setLength(counts[inputID]);
      vecInpNumSettings_[inputID] = counts[inputID];
      for (sID = 0; sID < counts[inputID]; sID++)
        vecInpSettings_[inputID][sID] = settings[inputID][sID];
    }
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FactorialSampling& FactorialSampling::operator=(const FactorialSampling &)
{
  printf("FactorialSampling operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

