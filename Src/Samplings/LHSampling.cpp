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
// Functions for the Latin hypercube class 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "LHSampling.h"
#define PABS(x) ((x) > 0 ? (x) : -(x))

// ************************************************************************
// external functions
// ************************************************************************
extern "C" 
{
  void OA_strength(int q,int nrow,int ncol,int** A,int *str,int verbose);
}

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
LHSampling::LHSampling() : Sampling()
{
  samplingID_ = PSUADE_SAMP_LHS;
  trueRandom_ = 0;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
LHSampling::~LHSampling()
{
}

// ************************************************************************
// initialize the sampling data
// ------------------------------------------------------------------------
int LHSampling::initialize(int initLevel)
{
  int    ss, jj, ii, repID, nReps, strength, index;
  int    ir, maxMinDist, minDist, dist, ss2, ntimes;
  double scale, ddata;
  char   pString[1000];
  psVector vecRanges;

  if (nInputs_ == 0)
  {
    printf("LHSampling::initialize ERROR - input not set up.\n");
    exit(1);
  }
  if (nSamples_ == 0)
  {
    printf("LHSampling::initialize ERROR - nSamples = 0.\n");
    exit(1);
  }

  vecSamInps_.clean();
  vecSamOuts_.clean();
  vecSamStas_.clean();
  nReps = nReplications_;
  trueRandom_ = 0;
  if (randomize_ & 2) trueRandom_ = 1;
  if (nSamples_ / nReps * nReps != nSamples_) 
  {
    printf("LHSampling : nSamples must be multiples of replications.\n");
    exit(1);
  }
  nSymbols_ = nSamples_ / nReps;
  if (initLevel != 0) return 0;

  if (printLevel_ > 4)
  {
    printf("LHSampling::initialize: nSamples = %d\n", nSamples_);
    printf("LHSampling::initialize: nInputs  = %d\n", nInputs_);
    printf("LHSampling::initialize: nOutputs = %d\n", nOutputs_);
    if (randomize_ != 0)
         printf("LHSampling::initialize: randomize on\n");
    else printf("LHSampling::initialize: randomize off\n");
    if (trueRandom_ != 0)
         printf("LHSampling::initialize: more randomize on\n");
    else printf("LHSampling::initialize: more randomize off\n");
    for (ii = 0; ii < nInputs_; ii++)
      printf("    LHSampling input %3d = [%e %e]\n", ii+1,
             vecLBs_[ii], vecUBs_[ii]);
    if (vecInpNumSettings_.length() > 0)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (vecInpNumSettings_[ii] != nSymbols_)
        {
          printf("LHSampling ERROR: inputSetting not compabible.\n");
          exit(1);
        }
      }
    }
  }

  allocSampleData();
  if (nSymbols_ == 1)
  {
    for (ss = 0; ss < nSamples_; ss++)
    {
      if (trueRandom_ != 0)
      {
        for (ii = 0; ii < nInputs_; ii++) 
        {
          ddata = vecUBs_[ii] - vecLBs_[ii];
          vecSamInps_[ss*nInputs_+ii] = PSUADE_drand() * ddata +
                                     vecLBs_[ii];
        }
      }
      else
      {
        for (ii = 0; ii < nInputs_; ii++)
          vecSamInps_[ss*nInputs_+ii] = 0.5*(vecUBs_[ii] + vecLBs_[ii]);
      }
    }
    return 0;
  }

  psIVector vecPerm, vecIStore, vecIT1, vecIT2;
  vecPerm.setLength(nSamples_*nInputs_);
  vecIStore.setLength(nSamples_*nInputs_);
  for (ss = 0; ss < nSamples_; ss+=nSymbols_) 
    for (jj = 0; jj < nSymbols_; jj++) 
      for (ii = 0; ii < nInputs_; ii++) vecPerm[(ss+jj)*nInputs_+ii] = jj;

  vecIT1.setLength(nSymbols_);
  vecIT2.setLength(nSymbols_);
  maxMinDist = 0;
  if      (nSamples_ > 10000) ntimes = 1;
  else if (nSamples_ > 9000)  ntimes = 2;
  else if (nSamples_ > 4000)  ntimes = 5;
  else if (nSamples_ > 1000)  ntimes = 10;
  else                        ntimes = 50;
   
  if (psSamExpertMode_ == 1 && isScreenDumpModeOn() == 1)
  {
    printf("LHSampling uses repeated generation and checking to\n");
    printf("  find the sub-optimal sample (using maxi-min, that is,\n");
    printf("  select the sample with the largest minimum distance\n");
    printf("  between sample points.\n");
    sprintf(pString,
      "LHSampling: number of times to search for maxi-min (1-1000): \n");
    ntimes = getInt(1, 10000000, pString);
  }

  for (ir = 0; ir < ntimes; ir++)
  {
    if (printLevel_ >= 5 && (ir % (ntimes/10+1) == 0)) 
      printf("LHSampling::initialize - max-min cycle %d (out of %d).\n",
             ir+1, ntimes); 
    for (ss = 0; ss < nSamples_; ss += nSymbols_) 
    {
      for (ii = 0; ii < nInputs_; ii++) 
      { 
        generateRandomIvector(nSymbols_, vecIT1.getIVector());
        for (jj = 0; jj < nSymbols_; jj++) 
           vecIT2[jj] = vecPerm[(ss+vecIT1[jj])*nInputs_+ii];
        for (jj = 0; jj < nSymbols_; jj++) 
           vecPerm[(ss+jj)*nInputs_+ii] = vecIT2[jj];
      }
    }
    if (ntimes == 1)
    {
      for (ss = 0; ss < nSamples_; ss++) 
        for (ii = 0; ii < nInputs_; ii++) 
          vecIStore[ss*nInputs_+ii] = vecPerm[ss*nInputs_+ii];
      break;
    }
    minDist = nInputs_ * nSymbols_;
    for (ss = 0; ss < nSamples_; ss++) 
    {
      for (ss2 = ss+1; ss2 < nSamples_; ss2++) 
      {
        dist = 0;
        for (ii = 0; ii < nInputs_; ii++) 
          dist += PABS(vecPerm[ss*nInputs_+ii] - vecPerm[ss2*nInputs_+ii]);
        if (dist > 0 && dist < minDist) minDist = dist;
      }
    }
    if (minDist > maxMinDist) 
    {
      for (ss = 0; ss < nSamples_; ss++) 
        for (ii = 0; ii < nInputs_; ii++) 
          vecIStore[ss*nInputs_+ii] = vecPerm[ss*nInputs_+ii];
      maxMinDist = minDist; 
      if (printLevel_ >= 1)
      {
        printf("LHSampling::current max-min distance    = %d (sample %d)\n",
               maxMinDist, ir+1);
        printf("            each input has max distance = %d\n",
               nSamples_-1);
        printf("            number of inputs            = %d\n",nInputs_);
      }
    }
  }
  vecPerm = vecIStore;

  if (nSamples_ < 1000)
  {
    int *intPtr = vecPerm.getIVector();
    int **int2Ptr = new int*[nSymbols_];
    if (printLevel_ >= 1 && ntimes > 1) 
      printf("LHSampling::initialize - checking.\n");
    for (repID = 0; repID < nReps; repID++)
    {
      for (ss = 0; ss < nSymbols_; ss++)
        int2Ptr[ss] = &(intPtr[(repID*nSymbols_+ss)*nInputs_]);
      OA_strength(nSymbols_, nSymbols_, nInputs_, int2Ptr, &strength, 0);
      if (strength < 1)
      {
        printf("LHSampling : failed strength test (%d,%d).\n",
               strength,repID);
        printf("   ==> Please consult PSUADE developers.\n");
        exit(1);
      }
    }
  }

  vecRanges.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) 
    vecRanges[ii] = vecUBs_[ii] - vecLBs_[ii];

  if (trueRandom_ != 0)
  {
    scale = 1.0 / ((double) nSymbols_);
    for (ss = 0; ss < nSamples_; ss++) 
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        index = vecPerm[ss*nInputs_+ii];
        ddata = (PSUADE_drand() + index) * scale;
        vecSamInps_[ss*nInputs_+ii] = ddata * vecRanges[ii] + vecLBs_[ii];
      }
    }
  } 

  else if (randomize_ != 0)
  {
    psVector vecPerturb;
    vecPerturb.setLength(nInputs_*nSymbols_);
    for (ii = 0; ii < nInputs_; ii++)
    {
      for (jj = 0; jj < nSymbols_; jj++)
        vecPerturb[ii*nSymbols_+jj] = PSUADE_drand() - 0.5;
    }
    scale = 1.0 / ((double) nSymbols_);
    for (ss = 0; ss < nSamples_; ss++)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        index = vecPerm[ss*nInputs_+ii];
        ddata = (vecPerturb[ii*nSymbols_+index] + index + 0.5) * scale;
        vecSamInps_[ss*nInputs_+ii] = ddata * vecRanges[ii] + vecLBs_[ii];
      }
    }
  }

  else 
  {
    scale = 1.0 / ((double) (nSymbols_ - 1));
    for (ss = 0; ss < nSamples_; ss++) 
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        index = vecPerm[ss*nInputs_+ii];
        if (vecInpNumSettings_.length() > 0 &&
            index < vecInpNumSettings_[ii])
        {
          vecSamInps_[ss*nInputs_+ii] = vecInpSettings_[ii][index];
        }
        else
        {
          ddata = scale * index;
          vecSamInps_[ss*nInputs_+ii] = ddata*vecRanges[ii]+vecLBs_[ii];
        }
      }
    }
  }
  for (ss = 0; ss < nSamples_*nOutputs_; ss++)
    vecSamOuts_[ss] = PSUADE_UNDEFINED;
  return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int LHSampling::refine(int refineRatio, int randomize, double thresh,
                       int nSamples, double *sampleErrors)
{
  int    ss, ii, jj, newNSymbols, sampleOffset, index, nReps, outputID;
  int    newSampleOffset, repID, strength, binCount, addNSymbols, nLevels;
  double scale, ddata;
  psIVector vecIT1, vecIT2, vecBinFlags, vecPerm, vecNewSamStas;
  psVector  vecPerturb, vecNewSamInps, vecNewSamOuts, vecBounds, vecRanges;

  (void) randomize;
  (void) thresh;
  (void) nSamples;
  (void) sampleErrors;

  nLevels = refineRatio;
  vecRanges.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) 
    vecRanges[ii] = vecUBs_[ii] - vecLBs_[ii];
  nReps = nSamples_ / nSymbols_;

  if (randomize_ != 0 || trueRandom_ != 0)
  {
    vecNewSamInps.setLength(nSamples_*nLevels*nInputs_);
    vecNewSamOuts.setLength(nSamples_*nLevels*nOutputs_);
    vecNewSamStas.setLength(nSamples_*nLevels);
    for (ss = 0; ss < nSamples_*nLevels; ss++)
    {
      vecNewSamStas[ss] = 0;
      for (outputID = 0; outputID < nOutputs_; outputID++)
        vecNewSamOuts[ss*nOutputs_+outputID] = PSUADE_UNDEFINED;
    }

    newNSymbols = nSymbols_ * nLevels;
    vecBounds.setLength(nInputs_*(newNSymbols+1));
    for (ii = 0; ii < nInputs_; ii++) 
    {
      for (jj = 0; jj <= newNSymbols; jj++) 
        vecBounds[ii*(newNSymbols+1)+jj] = vecRanges[ii] / newNSymbols * 
                             jj + vecLBs_[ii];
    }

    vecPerturb.setLength(nInputs_*nSymbols_*nLevels);
    for (ii = 0; ii < nInputs_; ii++)
    {
      for (jj = 0; jj < nSymbols_*nLevels; jj++)
        vecPerturb[ii*nSymbols_*nLevels+jj] = PSUADE_drand();
    }

    for (ss = 0; ss < nSamples_; ss++)
    {
      for (ii = 0; ii < nInputs_; ii++) 
      {
        ddata = vecSamInps_[ss*nInputs_+ii]; 
        if (ddata == vecBounds[ii*(newNSymbols+1)+newNSymbols]) 
             jj = newNSymbols;
        else
        {
          for (jj = 1; jj <= newNSymbols; jj++) 
            if (ddata < vecBounds[ii*(newNSymbols+1)+jj]) break;
        }
        jj--;
        if (jj >= newNSymbols)
        {
          printf("LHSampling::refine ERROR (3) - %d %d %e\n",jj,
                 newNSymbols,ddata);
          exit(1);
        }
        ddata -= vecLBs_[ii];
        ddata  = ddata / vecRanges[ii] * (double) newNSymbols;
        ddata -= (double) jj;
        vecPerturb[ii*nSymbols_*nLevels+jj] = ddata;
      }
    }

    vecPerm.setLength(nSymbols_*nInputs_);
    vecIT1.setLength(nSymbols_);
    vecIT2.setLength(nSymbols_);
    vecBinFlags.setLength(newNSymbols*nInputs_);

    sampleOffset = 0;
    newSampleOffset = 0;
    for (repID = 0; repID < nReps; repID++)
    {
      for (jj = 0; jj < nSymbols_; jj++) 
      {
        for (ii = 0; ii < nInputs_; ii++)
          vecNewSamInps[(newSampleOffset+jj)*nInputs_+ii] = 
                     vecSamInps_[(sampleOffset+jj)*nInputs_+ii];
        for (outputID = 0; outputID < nOutputs_; outputID++)
          vecNewSamOuts[(newSampleOffset+jj)*nOutputs_+outputID] = 
              vecSamOuts_[(sampleOffset+jj)*nOutputs_+outputID];
        vecNewSamStas[newSampleOffset+jj] = vecSamStas_[sampleOffset+jj]; 
      }
      newSampleOffset += nSymbols_;

      for (ii = 0; ii < nInputs_; ii++) 
        for (jj = 0; jj < newNSymbols; jj++) 
          vecBinFlags[ii*newNSymbols+jj] = 0;

      for (ss = sampleOffset; ss < sampleOffset+nSymbols_; ss++) 
      {
        for (ii = 0; ii < nInputs_; ii++) 
        {
          ddata = vecSamInps_[ss*nInputs_+ii]; 
          if (ddata == vecBounds[ii*(newNSymbols+1)+newNSymbols]) 
            jj = newNSymbols;
          else
          {
            for (jj = 1; jj < (newNSymbols+1); jj++) 
              if (ddata < vecBounds[ii*(newNSymbols+1)+jj]) break;
          }
          vecBinFlags[ii*newNSymbols+jj-1] = -1;
        }
      }

      for (ii = 0; ii < nInputs_; ii++)
      {
        for (jj = 0; jj < newNSymbols; jj++) 
          if (vecBinFlags[ii*newNSymbols+jj] == 0)
            vecBinFlags[ii*newNSymbols+jj] = jj;
        binCount = 0;
        for (jj = 0; jj < newNSymbols; jj++) 
        {
          if (vecBinFlags[ii*newNSymbols+jj] >= 0)
          {
            vecBinFlags[ii*newNSymbols+binCount] = 
                     vecBinFlags[ii*newNSymbols+jj];
            binCount++;
          }
        }
      }

      for (jj = 0; jj < binCount; jj++) 
        for (ii = 0; ii < nInputs_; ii++) 
          vecPerm[jj*nInputs_+ii] = vecBinFlags[ii*newNSymbols+jj];
      for (ii = 0; ii < nInputs_; ii++) 
      {
        generateRandomIvector(binCount, vecIT1.getIVector());
        for (jj = 0; jj < binCount; jj++)
          vecIT2[jj] = vecPerm[vecIT1[jj]*nInputs_+ii];
        for (jj = 0; jj < binCount; jj++)
          vecPerm[jj*nInputs_+ii] = vecIT2[jj];
      }

      if (trueRandom_ != 0)
      {
        scale = 1.0 / ((double) (nSymbols_*nLevels));
        for (ss = 0; ss < binCount; ss++) 
        {
          for (ii = 0; ii < nInputs_; ii++)
          {
            index =  vecPerm[ss*nInputs_+ii];
            ddata = (PSUADE_drand() + index) * scale;
            vecNewSamInps[(newSampleOffset+ss)*nInputs_+ii] = 
                      ddata * vecRanges[ii] + vecLBs_[ii];
          }
        }
      }
      else
      {
        scale = 1.0 / ((double) (nSymbols_*nLevels));
        for (ss = 0; ss < binCount; ss++) 
        {
          for (ii = 0; ii < nInputs_; ii++)
          {
            index = vecPerm[ss*nInputs_+ii];
            ddata = (vecPerturb[ii*nSymbols_*nLevels+index]+index)*scale;
            vecNewSamInps[(newSampleOffset+ss)*nInputs_+ii] = 
                      ddata * vecRanges[ii] + vecLBs_[ii];
          }
        }
      }
      newSampleOffset += binCount;
      sampleOffset += nSymbols_;
    }

    nSamples_ = nSamples_ * nLevels;
    nSymbols_ = nSymbols_ * nLevels;
    vecSamInps_ = vecNewSamInps;
    vecSamOuts_ = vecNewSamOuts;
    vecSamStas_ = vecNewSamStas;
  }
  else
  {

    newNSymbols = (nSymbols_ - 1) * nLevels + 1;
    addNSymbols = newNSymbols - nSymbols_;
    vecNewSamInps.setLength(newNSymbols*nReps*nInputs_);
    vecNewSamOuts.setLength(newNSymbols*nReps*nOutputs_);
    vecNewSamStas.setLength(newNSymbols*nReps);
    for (ss = 0; ss < newNSymbols*nReps; ss++)
    {
      vecNewSamStas[ss] = 0;
      for (outputID = 0; outputID < nOutputs_; outputID++)
        vecNewSamOuts[ss*nOutputs_+outputID] = PSUADE_UNDEFINED;
    }

    vecPerm.setLength(nSymbols_*nLevels*nInputs_);;
    vecIT1.setLength(nSymbols_);
    vecIT2.setLength(nSymbols_);

    sampleOffset = 0;
    newSampleOffset = 0;
    for (repID = 0; repID < nReps; repID++)
    {
      for (jj = 0; jj < nSymbols_; jj++) 
      {
        for (ii = 0; ii < nInputs_; ii++)
          vecNewSamInps[(newSampleOffset+jj)*nInputs_+ii] = 
                        vecSamInps_[(sampleOffset+jj)*nInputs_+ii];
        for (outputID = 0; outputID < nOutputs_; outputID++)
          vecNewSamOuts[(newSampleOffset+jj)*nOutputs_+outputID] = 
               vecSamOuts_[(sampleOffset+jj)*nOutputs_+outputID];
        vecNewSamStas[newSampleOffset+jj] = 
                      vecSamStas_[sampleOffset+jj]; 
      }
      newSampleOffset += nSymbols_;

      for (ii = 0; ii < nInputs_; ii++)
      {
        addNSymbols = 0;
        for (jj = 0; jj < newNSymbols; jj++) 
        {
          if ((jj % nLevels) != 0)
          {
            vecPerm[addNSymbols*nInputs_+ii] = jj;
            addNSymbols++;
          }
        }
        generateRandomIvector(addNSymbols, vecIT1.getIVector());
        for (jj = 0; jj < addNSymbols; jj++)
          vecIT2[jj] = vecPerm[vecIT1[jj]*nInputs_+ii];
        for (jj = 0; jj < addNSymbols; jj++)
          vecPerm[jj*nInputs_+ii] = vecIT2[jj];
      }

      scale = 1.0 / ((double) (newNSymbols - 1));
      for (jj = 0; jj < addNSymbols; jj++) 
      {
        for (ii = 0; ii < nInputs_; ii++)
        {
          ddata = (double) vecPerm[jj*nInputs_+ii];
          ddata = ddata * scale;
          vecNewSamInps[(newSampleOffset+jj)*nInputs_+ii] = 
                     ddata * vecRanges[ii] + vecLBs_[ii];
        }
      }
      newSampleOffset += addNSymbols;
      sampleOffset += nSymbols_;
    }

    nSymbols_ = newNSymbols;
    nSamples_ = nSymbols_ * nReps;
    vecSamInps_ = vecNewSamInps;
    vecSamOuts_ = vecNewSamOuts;
    vecSamStas_ = vecNewSamStas;
  }

  if (printLevel_ > 4)
  {
    printf("LHSampling refine: nSamples = %d\n", nSamples_);
    printf("LHSampling refine: nInputs  = %d\n", nInputs_);
    printf("LHSampling refine: nOutputs = %d\n", nOutputs_);
    if (randomize_ != 0)
         printf("LHSampling refine: randomize on\n");
    else printf("LHSampling refine: randomize off\n");
    if (trueRandom_ != 0)
         printf("LHSampling refine: more randomize on\n");
    else printf("LHSampling refine: more randomize off\n");
    for (ii = 0; ii < nInputs_; ii++)
      printf("    LHSampling input %3d = [%e %e]\n", ii+1,
             vecLBs_[ii], vecUBs_[ii]);
    if (vecInpNumSettings_.length() > 0 || vecSymTable_.length() > 0)
      printf("LHSampling refine: diable input settings, symbol table.\n");
  }

  if (nSamples_ < 1000)
  {
    vecBounds.setLength(nInputs_*(nSymbols_+1));
    for (ii = 0; ii < nInputs_; ii++) 
    {
      for (jj = 0; jj <= nSymbols_; jj++) 
        vecBounds[ii*(nSymbols_+1)+jj] = vecRanges[ii] / nSymbols_ * 
                             jj + vecLBs_[ii];
    }
    int **permMatrix = new int*[nSamples_/nReps];
    for (ss = 0; ss < nSamples_/nReps; ss++)
      permMatrix[ss] = new int[nInputs_];

    sampleOffset = 0;
    for (repID = 0; repID < nReps; repID++)
    {
      for (ss = 0; ss < nSamples_/nReps; ss++) 
      {
        for (ii = 0; ii < nInputs_; ii++) 
        {
          ddata = vecSamInps_[(sampleOffset+ss)*nInputs_+ii]; 
          if (ddata == vecBounds[ii*(nSymbols_+1)+nSymbols_]) 
            jj = nSymbols_;
          else
          {
            for (jj = 1; jj < (nSymbols_+1); jj++) 
              if (ddata < vecBounds[ii*(nSymbols_+1)+jj]) break;
          }
          permMatrix[ss][ii] = jj-1;
        }
      }
      OA_strength(nSamples_/nReps, nSymbols_, nInputs_,
                  permMatrix, &strength, 0);
      if (strength != 1)
        printf("LHS refine ERROR : replication %d : OA_strength = %d\n",
               repID, strength);
      sampleOffset += (nSamples_/nReps);
    }
    for (ss = 0; ss < nSamples_/nReps; ss++) delete [] permMatrix[ss];
    delete [] permMatrix;
  }
  return 0;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int LHSampling::setInputParams(int nInputs, int *counts, double **settings,
                               int *symtable)
{
  int ii, inputCnt, ss;

  if (nInputs_ != 0 && nInputs != nInputs_) 
  { 
    printf("LHSampling::setInputParams - nInputs mismatch.\n");
    exit(1);
  }
  nInputs_ = nInputs;
  if (symtable != NULL)
  {
    vecSymTable_.setLength(nInputs);
    for (ii = 0; ii < nInputs_; ii++) vecSymTable_[ii] = symtable[ii];
  }
  if (counts != NULL)
  {
    inputCnt = 0;
    maxNumSettings_ = 0;
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (counts[ii] != 0 && counts[ii] != nSymbols_)
      {
        printf("LHSampling::setInputParams - counts mismatch.\n");
        printf("            count data = %d %d\n", nSymbols_,
               counts[ii]);
        exit(1);
      }
      else if (counts[ii] == nSymbols_) inputCnt++;
      if (counts[ii] > maxNumSettings_) maxNumSettings_ = counts[ii];
    }
    if (inputCnt > 0)
    {
      vecInpSettings_ = new psVector[nInputs_];
      vecInpNumSettings_.setLength(nInputs_);
      for (ii = 0; ii < nInputs_; ii++)
      {
        if (counts[ii] == nSymbols_)
        {
          vecInpSettings_[ii].setLength(counts[ii]);
          vecInpNumSettings_[ii] = counts[ii];
          for (ss = 0; ss < counts[ii]; ss++)
            vecInpSettings_[ii][ss] = settings[ii][ss];
        }
        else vecInpNumSettings_[ii] = 0;
      } 
    }
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
LHSampling& LHSampling::operator=(const LHSampling &)
{
  printf("LHSampling operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

