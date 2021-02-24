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
// Functions for the Metis sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <unistd.h>

#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FuncApprox.h"
#include "Mars.h"
#include "MetisSampling.h"

#define PABS(x) ((x) > 0 ? x : -(x))

#ifdef HAVE_METIS
extern "C" 
{
void METIS_PartGraphRecursive(int *, int *, int *, int *, int *,
                              int *, int *, int *, int *, int *, int *);
}
#endif

int modeSet_ = 0;
int pruneSet_ = 0;

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
MetisSampling::MetisSampling() : Sampling()
{
  samplingID_ = PSUADE_SAMP_METIS;
  refineType_ = 0;
  refineSize_ = 1000000;
  n1d_        = -1;
  nAggrs_     = 0;
  vecAggrLabels_ = NULL;
  changeInfoName_ = 0;
}

// ************************************************************************
// copy constructor added by Bill Oliver
// ------------------------------------------------------------------------
MetisSampling::MetisSampling(const MetisSampling & ms) : Sampling()
{
  nSamples_ = ms.nSamples_;
  nInputs_  = ms.nInputs_;
  nOutputs_ = ms.nOutputs_;
  refineType_ = ms.refineType_;
  refineSize_ = ms.refineSize_;
  n1d_ = ms.n1d_;
  nAggrs_ = ms.nAggrs_;
  graphN_ = ms.graphN_;
  if (ms.vecAggrLabels_ != NULL)
  {
    vecAggrLabels_ = new psIVector[nAggrs_];
    for (int ii = 0; ii < nAggrs_; ii++)
      vecAggrLabels_[ii] = ms.vecAggrLabels_[ii]; 
  }
  vecAggrCnts_ = ms.vecAggrCnts_;
  printLevel_ = ms.printLevel_;
  samplingID_ = ms.samplingID_;
  randomize_ = ms.randomize_;
  nReplications_ = ms.nReplications_;
  vecLBs_ = ms.vecLBs_;
  vecUBs_ = ms.vecUBs_;
  vecSamInps_ = ms.vecSamInps_;
  vecSamOuts_ = ms.vecSamOuts_;
  vecSamStas_ = ms.vecSamStas_;
  vecGraphI_ = ms.vecGraphI_;
  vecGraphJ_ = ms.vecGraphJ_;
  vecCellsOccupied_ = ms.vecCellsOccupied_;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MetisSampling::~MetisSampling()
{
  if (vecAggrLabels_ != NULL)
  {
    for (int ii = 0; ii < nAggrs_; ii++) vecAggrLabels_[ii].clean();
    delete [] vecAggrLabels_;
    vecAggrLabels_ = NULL;
  }
}

// ************************************************************************
// initialize the sampler data
// ------------------------------------------------------------------------
int MetisSampling::initialize(int initLevel)
{
  int    inputID, ii, jj, kk, itmp, jtmp, nnz, sampleID, randFlag;
  int    options[10], index, count, saveFlag=1;
#ifdef HAVE_METIS
  int    wgtflag=0, numflag=0, edgeCut=0;
#endif
  double dtmp, expand;
  char   response[1001], pString[1001], filename[1001];
  FILE   *fp;
  psIVector vecIncrs;
  psVector  vecRanges, vecLBs, vecUBs;
 
  if( nSamples_ <= 0)
  {
    printf("nSamples_ in file %s line %d is <= 0\n", __FILE__, __LINE__);
    exit(1);
  }
  if (nInputs_ == 0)
  {
    printf("MetisSampling::initialize ERROR - input not set up.\n");
    exit(1);
  }
  randFlag  = randomize_;

  if (nInputs_ > 22)
  {
    printf("MetisSampling ERROR: nInputs > 22 currently not supported.\n");
    exit(1);
  }
  if (nInputs_ == 1 ) n1d_ = nSamples_*10;
  if (nInputs_ == 2 ) n1d_ = 2048;
  if (nInputs_ == 3 ) n1d_ = 150;
  if (nInputs_ == 4 ) n1d_ = 44;
  if (nInputs_ == 5 ) n1d_ = 20;
  if (nInputs_ == 6 ) n1d_ = 12;
  if (nInputs_ == 7 ) n1d_ = 8;
  if (nInputs_ == 8 ) n1d_ = 6;
  if (nInputs_ == 9 ) n1d_ = 5;
  if (nInputs_ == 10) n1d_ = 4;
  if (nInputs_ == 11) n1d_ = 4;
  if (nInputs_ == 12) n1d_ = 3;
  if (nInputs_ == 13) n1d_ = 3;
  if (nInputs_ >= 14) n1d_ = 3;
  if (nInputs_ >= 15) n1d_ = 2;

  vecIncrs.setLength(nInputs_+1);
  graphN_ = 1;
  vecIncrs[0] = graphN_;
  for (inputID = 1; inputID <= nInputs_; inputID++)
  {
    graphN_ *= n1d_;
    vecIncrs[inputID] = graphN_;
  }
  if (nSamples_ > 2*graphN_)
  {
    printf("MetisSampling ERROR : nSamples %d too large.\n",nSamples_);
    exit(1);
  }

  vecGraphI_.setLength(graphN_+1);
  vecGraphJ_.setLength(graphN_*nInputs_*2+1);
  nnz = 0;
  vecGraphI_[0] = nnz;
  for (ii = 0; ii < graphN_; ii++)
  {
    itmp = ii;
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      jtmp = itmp % n1d_;
      itmp = itmp / n1d_;
      if (jtmp > 0     ) vecGraphJ_[nnz++] = ii - vecIncrs[inputID];
      if (jtmp < n1d_-1) vecGraphJ_[nnz++] = ii + vecIncrs[inputID];
    }
    vecGraphI_[ii+1] = nnz;
  }
  vecCellsOccupied_.setLength(graphN_);

  if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "r");
  else                      fp = fopen("psuadeMetisInfo.tmp", "r");
  if (fp != NULL)
  {
    if (vecAggrLabels_ != NULL)
    {
      for (int ii = 0; ii < nAggrs_; ii++) vecAggrLabels_[ii].clean();
      delete [] vecAggrLabels_;
      vecAggrLabels_ = NULL;
    }
    printf("INFO: psuadeMetisInfo file found. Reading it in ...\n");
    fscanf(fp, "%d %d %d", &jj, &itmp, &jtmp);
    if (itmp != nSamples_ || jtmp != nInputs_ || jj != nSamples_)
    {
      fclose(fp);
      printf("MetisSampling INFO: a partition file is found but\n");
      printf("      the data is not consistent with this setup\n");
      printf("      (The file name is psuadeMetisInfo).\n");
      if (itmp != nSamples_ || jj != nSamples_)
        printf("      nSamples : %d != %d.\n", nSamples_, itmp);
      if (jtmp != nInputs_)
        printf("      nInputs  : %d != %d.\n", nInputs_, jtmp);
      sprintf(pString,"Would you like to provide another file? (y or n) ");
      getString(pString, response);
      if (response[0] == 'y')
      {
        sprintf(pString,"Enter partition file name : ");
        getString(pString, filename);
        fp = fopen(filename, "r");
        if (fp == NULL)
        {
          printf("MetisSampling ERROR: partition file not found.\n");
          exit(1);
        }
      }
      else
      {
        printf("MetisSampling INFO: delete psuadeMetisInfo file and.\n");
        printf("                    re-launch.\n");
        exit(1);
      }
      fscanf(fp, "%d %d %d", &jj, &itmp, &jtmp);
      if (itmp != nSamples_ || jtmp != nInputs_)
      {
        printf("MetisSampling ERROR: partition file not valid.\n");
        exit(1);
      }
      printf("Metis INFO: partition file has %d subdomains\n",jj);
      printf("            and %d sample points.\n",itmp);
      printf("Metis INFO: reconstructing the partitioning.\n");
    }
    printf("      Incoming nSamples : %d.\n", itmp);
    printf("      Incoming nInputs  : %d.\n", jtmp);
    nAggrs_ = jj;
    vecAggrCnts_.setLength(nAggrs_);
    vecAggrLabels_ = new psIVector[nAggrs_];
    for (ii = 0; ii < nAggrs_; ii++)
    {
      fscanf(fp, "%d", &count);
      if (printLevel_ > 4) 
        printf("Metis read: aggr %8d, size = %d\n", ii+1, count);
      if (count > 0)
      {
        vecAggrCnts_[ii] = count;
        vecAggrLabels_[ii].setLength(count);
      }
      else vecAggrCnts_[ii] = count = 0;
      for (jj = 0; jj < count; jj++)
      {
        fscanf(fp, "%d", &kk);
        vecAggrLabels_[ii][jj] = kk;
        if (vecAggrLabels_[ii][jj] < 0 || vecAggrLabels_[ii][jj] >= graphN_)
        {
          printf("Metis ERROR: psuadeMetisInfo file has invalid info.\n");
          exit(1);
        }
        vecCellsOccupied_[vecAggrLabels_[ii][jj]] = ii;
      }
    }
    fscanf(fp, "%d", &count);
    for (ii = 0; ii < count; ii++)
    {
      fscanf(fp, "%d %d", &jj, &kk);
      if (jj < 0 || jj >= graphN_)
      {
        printf("Metis ERROR: psuadeMetisInfo file has invalid info (2).\n");
        exit(1);
      }
      vecCellsOccupied_[jj] = - kk - 1;
    }
    fclose(fp);
    saveFlag = 0;
  }

  if (vecAggrCnts_.length() == 0)
  {
    options[0] = 0;
    if (printLevel_ > 1)
      printf("MetisSampling:: calling domain partitioner.\n");
#ifdef HAVE_METIS
    METIS_PartGraphRecursive(&graphN_, vecGraphI_.getIVector(), 
          vecGraphJ_.getIVector(), NULL, NULL, &wgtflag,&numflag,
          &nSamples_,options,&edgeCut,vecCellsOccupied_.getIVector());
#else
    printf("MetisSampling ERROR : METIS not installed.\n");
    exit(1);
#endif
    if (printLevel_ > 1)
      printf("MetisSampling:: subdomains created.\n");

    nAggrs_ = nSamples_;
    vecAggrCnts_.setLength(nAggrs_);
    for (ii = 0; ii < nAggrs_; ii++) vecAggrCnts_[ii] = 0;
    for (ii = 0; ii < graphN_; ii++)
    {
      if (vecCellsOccupied_[ii] < 0 || vecCellsOccupied_[ii] >= nAggrs_)
      {
        printf("MetisSampling INTERNAL ERROR.\n");
        exit(1);
      }
      vecAggrCnts_[vecCellsOccupied_[ii]]++;
    }  
    //**  belong to it 
    vecAggrLabels_ = new psIVector[nAggrs_];
    for (ii = 0; ii < nAggrs_; ii++)
    {
      if (vecAggrCnts_[ii] <= 0)
      {
        printf("MetisSampling INTERNAL ERROR (2).\n");
        exit(1);
      }
      vecAggrLabels_[ii].setLength(vecAggrCnts_[ii]);
      vecAggrCnts_[ii] = 0;
    }
    for (ii = 0; ii < graphN_; ii++)
    {
      index = vecCellsOccupied_[ii];
      if (index < 0 || index >= nAggrs_)
      {
        printf("MetisSampling INTERNAL ERROR (3).\n");
        exit(1);
      }
      vecAggrLabels_[index][vecAggrCnts_[index]] = ii;  
      vecAggrCnts_[index]++;
    }
  }

  if (saveFlag == 1)
  {
    if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "w");
    else                      fp = fopen("psuadeMetisInfo.tmp", "w");
    if (fp != NULL)
    {
      fprintf(fp, "%d %d %d\n", nAggrs_, nSamples_, nInputs_);
      for (ii = 0; ii < nAggrs_; ii++)
      {
        fprintf(fp, "%d\n", vecAggrCnts_[ii]);
        for (jj = 0; jj < vecAggrCnts_[ii]; jj++)
        {
          kk = vecAggrLabels_[ii][jj];
          fprintf(fp, "%d ", kk);
          if (jj != 0 && jj % 10 == 0) fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
      }
      jj = 0;
      for (ii = 0; ii < graphN_; ii++) if (vecCellsOccupied_[ii] < 0) jj++;
      fprintf(fp, "%d\n", jj);
      for (ii = 0; ii < graphN_; ii++)
        if (vecCellsOccupied_[ii] < 0) 
          fprintf(fp, "%d %d\n",ii,-(vecCellsOccupied_[ii]+1));
      fclose(fp);
    }
  }
  if (initLevel != 0) return 0;

  char *inStr, winput1[500], winput2[500];
  expand = 0.0;
  if (psSamExpertMode_ == 1)
  {
    printf("For the Metis sampling, you can divert more data points\n");
    printf("to the surface of the parameter space by sampling from\n");
    printf("the expanded parameter space folowed by projecting the\n");
    printf("sample points back to the original parameter space.\n");
    printf("In order to do so, an expansion ratio needs to be given.\n");
    printf("An expansion ratio of 0.0 means no expansion.\n");
    printf("Usually the expansion ratio should be no more than 0.1-0.2.\n");
    sprintf(pString, "Enter an expansion ratio (default=0): ");
    expand = -0.1;
    while (expand < 0.0) expand = getDouble(pString);
  }
  else if (psConfig_ != NULL)
  {
    inStr = psConfig_->getParameter("METIS_expand_ratio");
    if (inStr != NULL)
    {
      sscanf(inStr, "%s %s %lg\n", winput1, winput2, &expand);
      if (winput2[0] != '=')
      {
        printf("METIS read config file syntax error : %s\n", inStr);
        expand = 0.0;
      }
      else if (expand < 0.0 || expand >= 1.0)
      {
        printf("METIS read config file read expand error : %e\n",expand);
        expand = 0.0;
      }
      else printf("METIS : expansion ratio set to %e\n", expand);
    }
  }

  allocSampleData();
  vecRanges.setLength(nInputs_);
  vecLBs.setLength(nInputs_);
  vecUBs.setLength(nInputs_);
  for (inputID = 0; inputID < nInputs_; inputID++)
  {
    vecRanges[inputID] = (vecUBs_[inputID]-vecLBs_[inputID])*0.5*expand;
    vecLBs[inputID] = vecLBs_[inputID] - vecRanges[inputID];
    vecUBs[inputID] = vecUBs_[inputID] + vecRanges[inputID];
  }
  for (inputID = 0; inputID < nInputs_; inputID++)
    vecRanges[inputID] = vecUBs[inputID] - vecLBs[inputID];

  for (sampleID = 0; sampleID < nSamples_; sampleID++)
  {
    if ((randFlag & 1) == 1)
         index = (int) (PSUADE_drand() * vecAggrCnts_[sampleID]);
    else index = vecAggrCnts_[sampleID] / 2;
    if (index == vecAggrCnts_[sampleID]) index--;
    index = vecAggrLabels_[sampleID][index];
    vecCellsOccupied_[index] = -(vecCellsOccupied_[index] + 1);
    itmp = index;
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      jtmp = itmp % n1d_;
      itmp = itmp / n1d_;
      if ((randFlag & 1) == 1)
           dtmp = ((double) (jtmp + 0.999*PSUADE_drand())) / (double) n1d_;
      else dtmp = ((double) (jtmp + 0.5)) / (double) n1d_;
      vecSamInps_[sampleID*nInputs_+inputID] = dtmp * vecRanges[inputID] +
                                               vecLBs[inputID];
      if (vecSamInps_[sampleID*nInputs_+inputID] < vecLBs_[inputID])
        vecSamInps_[sampleID*nInputs_+inputID] = vecLBs_[inputID];
      if (vecSamInps_[sampleID*nInputs_+inputID] > vecUBs_[inputID])
        vecSamInps_[sampleID*nInputs_+inputID] = vecUBs_[inputID];
    }
  }

  if (printLevel_ > 4)
  {
    printf("MetisSampling::initialize: nSamples = %d\n", nSamples_);
    printf("MetisSampling::initialize: nInputs  = %d\n", nInputs_);
    printf("MetisSampling::initialize: nOutputs = %d\n", nOutputs_);
    if (randomize_ != 0)
         printf("MetisSampling::initialize: randomize on\n");
    else printf("MetisSampling::initialize: randomize off\n");
    for (inputID = 0; inputID < nInputs_; inputID++)
       printf("    MetisSampling input %3d = [%e %e]\n", inputID+1,
              vecLBs_[inputID], vecUBs_[inputID]);
  }
  return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int MetisSampling::refine(int nLevels, int randFlag, double threshold,
                          int nSamples, double *sampleErrors)
{
  int    inputID, ii, jj, count, itmp, jtmp, localN, count1, iOne=1;
  int    maxN, rowInd, colInd, localNNZ, index, count0, kk, options[10];
#ifdef HAVE_METIS
  int    wgtflag=0, numflag=0, edgeCut=0, itwo=2;
#endif
  int    currNAggr, newCell, oldNumSamples, mode;
  int    status, outputID, ss, maxNNZ, marsCnt, aggr1, splitCount;
  int    splitSuccess, maxCellSize, minCellSize, avgCellSize, cellCnt=0;
  int    aggr2, rowInd2, colInd2, localN2, ii2, jj2, limit;
  double ddmax, diffCnt, ddmin, ddmean, thresh, ddsd, dtmp;
  char   pString[500];
  FILE   *fp;
  FuncApprox *faPtr;
  psIVector vecSubLabels, vecLocalIA, vecLocalJA, vecTmpCnts, vecTags; 
  psIVector vecOldSamStas, vecLabels, vecRefineFlags, vecNode2Aggr;
  psIVector vecTags2, vecSubLabels2, *vecTmpLabels;
  psVector  vecRanges, vecOldSamInps, vecOldSamOuts, vecLBs, vecUBs;
  psVector  vecDiffs, vecSortList1, vecSortList2, vecMarsIn, vecMarsOut;

  if (nSamples != 0 && nSamples != nAggrs_)
  {
    printf("MetisSampling::refine ERROR - sampleErrors length mismatch.\n");
    printf("               This may be due to the use of wrong\n");
    printf("               psuadeMetisInfo file.\n");
    exit(1);
  }
  if (vecCellsOccupied_.length() == 0 || vecAggrLabels_ == NULL)
  {
    printf("MetisSampling::refine ERROR - need to call initialize first.\n");
    exit(1);
  }
  if (printLevel_ > 0)
  {
    printf("MetisSampling::refine(1): nSamples = %d\n", nSamples_);
    printf("MetisSampling::refine(1): nInputs  = %d\n", nInputs_);
    printf("MetisSampling::refine(1): nOutputs = %d\n", nOutputs_);
  }

  vecRanges.setLength(nInputs_);
  vecLBs.setLength(nInputs_);
  vecUBs.setLength(nInputs_);
  for (inputID = 0; inputID < nInputs_; inputID++)
  {
    vecRanges[inputID]  = (vecUBs_[inputID] - vecLBs_[inputID])*0.0;
    vecLBs[inputID] = vecLBs_[inputID] - vecRanges[inputID];
    vecUBs[inputID] = vecUBs_[inputID] + vecRanges[inputID];
  }
  for (inputID = 0; inputID < nInputs_; inputID++)
  {
    vecRanges[inputID] = vecUBs[inputID] - vecLBs[inputID];
    if (vecRanges[inputID] <= 0.0)
    {
      printf("MetisSampling::refine ERROR - lbound/ubound mismatch.\n");
      exit(1);
    }
  }

  count = 0;
  for (ii = 0; ii < graphN_; ii++) 
    if (vecCellsOccupied_[ii] < 0) count++;

  if (count == 0)
  {
    for (ss = 0; ss < nSamples_; ss++)
    {
      itmp = 0;
      for (inputID = nInputs_-1; inputID >= 0; inputID--)
      {
        itmp = itmp * n1d_;
        dtmp = vecSamInps_[ss*nInputs_+inputID];
        dtmp = (dtmp - vecLBs_[inputID]) / vecRanges[inputID];
        if (dtmp == 1.0) jtmp = n1d_ - 1;
        else             jtmp = (int) (dtmp * n1d_);
        itmp += jtmp;
      }
      if (itmp < 0 || itmp >= graphN_)
      {
        printf("MetisSampling::refine INTERNAL ERROR.\n");
        printf("               Consult PSUADE developer.\n");
      }
      vecCellsOccupied_[itmp] = -(vecCellsOccupied_[itmp] + 1); 
    }
  }
  else if (count != nSamples_)
  {
    printf("MetisSampling::refine ERROR - in re-initialize.\n");
    printf("               Consult PSUADE developer.\n");
    printf("               Sample size mismatch: %d vs %d\n",
           count, nSamples_);
    printf("                between psuadeMetisInfo and current data set.\n");
    exit(1);
  }
   
  count = 0;
  for (ii = 0; ii < graphN_; ii++) 
    if (vecCellsOccupied_[ii] < 0) count -= (vecCellsOccupied_[ii] + 1);
  itmp = (nSamples_ - 1) * nSamples_ / 2;
  if (nSamples_ < 40000 && count != itmp)
  {
    printf("MetisSampling::refine ERROR - METIS not used in initialize.\n");
    printf("               so Metis cannot be used in refine (%d,%d).\n",
           count, itmp);
    printf("Note: It can be due to psuadeMetisInfo file being modified.\n");
    exit(1);
  }
  else if (nSamples_ >= 40000)
  {
    printf("MetisSampling::refine INFO - less error checking (N>40000).\n");
  }

  vecTmpCnts = vecAggrCnts_;
  vecTmpLabels = vecAggrLabels_;
  vecAggrCnts_.setLength(2*nAggrs_);
  vecAggrLabels_ = new psIVector[2*nAggrs_];
  for (ii = 0; ii < nAggrs_; ii++)
  {
    vecAggrCnts_[ii] = vecTmpCnts[ii];
    vecAggrLabels_[ii].setLength(vecAggrCnts_[ii]);
    for (jj = 0; jj < vecAggrCnts_[ii]; jj++) 
      vecAggrLabels_[ii][jj] = vecTmpLabels[ii][jj];
    vecTmpLabels[ii].clean();
  }
  delete [] vecTmpLabels;
  for (ii = nAggrs_; ii < 2*nAggrs_; ii++) 
  {
    vecAggrCnts_[ii] = 0;
    vecAggrLabels_[ii].clean();
  }

  maxN = 0;
  for (ss = 0; ss < nAggrs_; ss++)
    if (vecAggrCnts_[ss] > maxN) maxN = vecAggrCnts_[ss];
  vecLabels.setLength(maxN);
  vecLocalIA.setLength(2*maxN+1);
  maxNNZ = maxN * (2 * nInputs_ + 1);
  vecLocalJA.setLength(maxNNZ);

  if (refineType_ == 1 && sampleErrors == NULL)
  {
    printf("MetisSampling::refine ERROR - error based but no error given.\n");
    exit(1);
  }
  if (refineType_ == 1 && sampleErrors != NULL)
  {
    printf("MetisSampling::refine - maximum number of new points = %d.\n",
           refineSize_);
    printf("Note: error scaled by the cell size.\n");
    vecNode2Aggr.setLength(graphN_);
    for (ii = 0; ii < graphN_; ii++) vecNode2Aggr[ii] = -1;
    for (ii = 0; ii < nAggrs_; ii++)
    {
      localN = vecAggrCnts_[ii];
      vecSubLabels = vecAggrLabels_[ii];
      for (jj = 0; jj < localN; jj++)
      {
        rowInd = vecSubLabels[jj];
        if (rowInd < 0 || rowInd >= graphN_)
        {
          printf("MetisSampling::refine ERROR (index out of bound.)\n");
          printf("               rowInd = %d (should be in [0,%d])\n",
                 rowInd, graphN_-1);
          printf("               Aggregate = %d (loc=%d)\n", ii, jj);
          printf("               Please consult PSUADE developers.\n");
          exit(1);
        }
        vecNode2Aggr[rowInd] = ii;
      }
    }
    for (ii = 0; ii < graphN_; ii++)
    {
      if (vecNode2Aggr[ii] == -1)
      {
        printf("MetisSampling::refine ERROR - node2Aggr not correct.\n");
        exit(1);
      }
    }

    vecSortList1.setLength(nAggrs_);
    vecSortList2.setLength(nAggrs_);
    for (ii = 0; ii < nAggrs_; ii++)
    {
      localN = vecAggrCnts_[ii];
      vecSortList1[ii] = PABS(sampleErrors[ii]);
      vecSortList1[ii] *= pow(1.0*localN,1.0/nInputs_);
    }
    for (ii = 0; ii < nAggrs_; ii++) vecSortList2[ii] = (double) ii;
    sortDbleList2(nAggrs_,vecSortList1.getDVector(),
                  vecSortList2.getDVector());
    vecRefineFlags.setLength(nAggrs_);
    for (ii = nAggrs_-1; ii >= 0; ii--) vecRefineFlags[ii] = 0;
    cellCnt = 0;

    for (ii = nAggrs_-1; ii >= 0; ii--)
    {
      index = (int) vecSortList2[ii];
      localN = vecAggrCnts_[index];
      if (sampleErrors[index] != 0.0) 
      {
        //printf("MetisSampling::refine - chosen aggregate, error = %e\n",
        //        sampleErrors[index]);
        printf("MetisSampling: selected refinement: error %13.5e, leng %d\n",
               sampleErrors[index], localN);
        if (vecRefineFlags[index] == 0)
        {
          cellCnt++;
          vecRefineFlags[index] = 1;
          if (cellCnt >= refineSize_) break;
        }
#if 0
        vecSubLabels = vecAggrLabels_[index];
        for (jj = 0; jj < localN; jj++)
        {
          rowInd = vecSubLabels[jj];
          for (kk = vecGraphI_[rowInd]; kk < vecGraphI_[rowInd+1]; kk++)
          {
            colInd = vecGraphJ_[kk];
            aggr1 = vecNode2Aggr[colInd];
            if (vecRefineFlags[aggr1] == 0)
            {
              cellCnt++;
              vecRefines[aggr1] = 1;
              if (cellCnt >= refineSize_) break;
            }
          }
          if (cellCnt >= refineSize_) break;
        }
#endif
        if (cellCnt >= refineSize_) break;
      }
      if (cellCnt >= refineSize_) break;
    }
  }

  if (refineType_ == 2)
  {
    printf("MetisSampling::refine - maximum number of new points = %d.\n",
           refineSize_);
    vecTags.setLength(nAggrs_);
    for (ii = 0; ii < nAggrs_; ii++) vecTags[ii] = 0;
    vecNode2Aggr.setLength(graphN_);
    for (ii = 0; ii < nAggrs_; ii++)
    {
      localN = vecAggrCnts_[ii];
      vecSubLabels = vecAggrLabels_[ii];
      for (jj = 0; jj < localN; jj++)
      {
        rowInd = vecSubLabels[jj];
        if (rowInd < 0 || rowInd >= graphN_)
        {
          printf("MetisSampling::refine ERROR (index out of bound.)\n");
          printf("               rowInd = %d (should be in [0,%d])\n",
                 rowInd, graphN_-1);
          printf("               Aggregate = %d (loc=%d)\n", ii, jj);
          printf("               Please consult PSUADE developers.\n");
          printf("               next loc index = %d\n",vecSubLabels[jj+1]);
          exit(1);
        }
        vecNode2Aggr[rowInd] = ii;
      }
    }

    vecDiffs.setLength(nAggrs_);
    for (ss = 0; ss < nAggrs_; ss++)
    {
      vecDiffs[ss] = 0;
      localN = vecAggrCnts_[ss];
      vecSubLabels = vecAggrLabels_[ss];

      ddmean = vecSamOuts_[ss];
      count0 = 1;
      for (ii = 0; ii < localN; ii++)
      {
        rowInd = vecSubLabels[ii];
        for (jj = vecGraphI_[rowInd]; jj < vecGraphI_[rowInd+1]; jj++)
        {
          colInd = vecGraphJ_[jj];
          aggr1 = vecNode2Aggr[colInd];
          if (aggr1 != ss && vecTags[aggr1] == 0)
          {
            vecTags[aggr1] = 1;
            ddmean += vecSamOuts_[aggr1];
            count0++;
            localN2 = vecAggrCnts_[aggr1];
            vecSubLabels2 = vecAggrLabels_[aggr1];
            for (ii2 = 0; ii2 < localN2; ii2++)
            {
              rowInd2 = vecSubLabels2[ii2];
              for (jj2 = vecGraphI_[rowInd2];jj2 < vecGraphI_[rowInd2+1];jj2++)
              {
                colInd2 = vecGraphJ_[jj2];
                aggr2 = vecNode2Aggr[colInd2];
                if (aggr2 != ss && vecTags[aggr2] == 0)
                {
                  vecTags[aggr2] = 1;
                  ddmean += vecSamOuts_[aggr2];
                  count0++;
                }
              }
            }
          }
        }
      }
      ddmean /= count0;
      ddsd = 0.0;
      if (count0 > 1)
      {
        ddsd += (vecSamOuts_[ss] - ddmean) * (vecSamOuts_[ss] - ddmean);
        for (ii = 0; ii < localN; ii++)
        {
          rowInd = vecSubLabels[ii];
          for (jj = vecGraphI_[rowInd]; jj < vecGraphI_[rowInd+1]; jj++)
          {
            colInd = vecGraphJ_[jj];
            aggr1 = vecNode2Aggr[colInd];
            if (aggr1 != ss && vecTags[aggr1] == 1)
            {
              dtmp = vecSamOuts_[aggr1];
              ddsd += (dtmp - ddmean) * (dtmp - ddmean);
              vecTags[aggr1] = 0;
              localN2 = vecAggrCnts_[aggr1];
              vecSubLabels2 = vecAggrLabels_[aggr1];
              for (ii2 = 0; ii2 < localN2; ii2++)
              {
                rowInd2 = vecSubLabels2[ii2];
                for (jj2=vecGraphI_[rowInd2];jj2<vecGraphI_[rowInd2+1];jj2++)
                {
                  colInd2 = vecGraphJ_[jj2];
                  aggr2 = vecNode2Aggr[colInd2];
                  if (aggr2 != ss && vecTags[aggr2] == 1)
                  {
                    vecTags[aggr2] = 0;
                    dtmp = vecSamOuts_[aggr2];
                    ddsd += (dtmp - ddmean) * (dtmp - ddmean);
                  }
                }
              }
            }
          }
        }
      }
      if (count0 > 0) vecDiffs[ss] = sqrt(ddsd) / count0;

      if (printLevel_ > 4 && nInputs_ == 2)
      {
        for (ii = 0; ii < nAggrs_; ii++) vecTags[ii] = 0;
        printf(" ..... Sample %6d inputs : %e %e = %e\n", ss+1, 
               vecSamInps_[ss*nInputs_],vecSamInps_[ss*nInputs_+1],
               vecSamOuts_[ss]);
        for (ii = 0; ii < localN; ii++)
        {
          rowInd = vecSubLabels[ii];
          for (jj = vecGraphI_[rowInd]; jj < vecGraphI_[rowInd+1]; jj++)
          {
            colInd = vecGraphJ_[jj];
            aggr1 = vecNode2Aggr[colInd];
            if (aggr1 != ss && vecTags[aggr1] == 0)
            {
              printf(" .....          neighbor : %e %e = %e\n", 
                 vecSamInps_[aggr1*nInputs_],vecSamInps_[aggr1*nInputs_+1],
                 vecSamOuts_[aggr1]);
              vecTags[aggr1] = 1;
            }
          }
        }
      }
    }

    diffCnt = 0.0;
    for (ii = 0; ii < nAggrs_; ii++) diffCnt += vecDiffs[ii];

    vecRefineFlags.setLength(nAggrs_);
    for (ii = 0; ii < nAggrs_; ii++) vecRefineFlags[ii] = 0;
    maxCellSize = cellCnt = 0;
    avgCellSize = 0;
    minCellSize = 1000000000;
    if (diffCnt == 0 || diffCnt == nAggrs_)
    { 
      for (ii = 0; ii < nAggrs_; ii++) 
      {
        vecRefineFlags[ii] = 0;
        if (localN > 1)
        {
          vecRefineFlags[ii] = 1;
          localN = vecAggrCnts_[ii];
          maxCellSize = (localN > maxCellSize) ? localN : maxCellSize;
          minCellSize = (localN < maxCellSize) ? localN : minCellSize;
          avgCellSize += localN;
          cellCnt++;
        }
      }
    }
    else
    {
      mode = 2;
      if (psSamExpertMode_ == 1 && modeSet_ == 0)
      {
        modeSet_ = 1;
      }
      limit = 0;
      if (mode == 2)
      {
        for (ss = 0; ss < nAggrs_; ss++) vecDiffs[ss] *= vecAggrCnts_[ss];
      }
      else if (mode == 3 || mode == 4)
      {
        if (mode == 4)
        {
          limit = nInputs_ + 1 + nInputs_ * (nInputs_ + 1) / 2;
        }
        vecMarsIn.setLength(nAggrs_*nInputs_);
        vecMarsOut.setLength(nAggrs_);
        vecTags2.setLength(nAggrs_);
        for (ss = 0; ss < nAggrs_; ss++)
        {
          for (ii = 0; ii < nAggrs_; ii++) vecTags2[ii] = 0;
          vecTags2[ss] = 1;
          marsCnt = 0;
          localN = vecAggrCnts_[ss];
          vecSubLabels = vecAggrLabels_[ss];
          for (ii = 0; ii < localN; ii++)
          {
            rowInd = vecSubLabels[ii];
            for (jj = vecGraphI_[rowInd]; jj < vecGraphI_[rowInd+1]; jj++)
            {
              colInd = vecGraphJ_[jj];
              aggr1 = vecNode2Aggr[colInd];
              if (vecTags2[aggr1] == 0)
              {
                for (kk = 0; kk < nInputs_; kk++)
                  vecMarsIn[marsCnt*nInputs_+kk]=vecSamInps_[aggr1*nInputs_+kk];
                vecMarsOut[marsCnt++] = vecSamOuts_[aggr1];
                vecTags2[aggr1] = 1;
              }
            }
          }
          if (marsCnt < limit)
          {
            for (ii = 0; ii < localN; ii++)
            {
              rowInd = vecSubLabels[ii];
              for (jj = vecGraphI_[rowInd]; jj < vecGraphI_[rowInd+1]; jj++)
              {
                colInd = vecGraphJ_[jj];
                aggr1 = vecNode2Aggr[colInd];
                localN2 = vecAggrCnts_[aggr1];
                vecSubLabels2 = vecAggrLabels_[aggr1];
                for (ii2 = 0; ii2 < localN2; ii2++)
                {
                  rowInd2 = vecSubLabels2[ii2];
                  for (jj2=vecGraphI_[rowInd2];jj2<vecGraphI_[rowInd2+1];jj2++)
                  {
                    colInd2 = vecGraphJ_[jj2];
                    aggr2 = vecNode2Aggr[colInd2];
                  }
                  if (vecTags2[aggr2] == 0)
                  {
                    for (kk = 0; kk < nInputs_; kk++)
                      vecMarsIn[marsCnt*nInputs_+kk] = 
                                     vecSamInps_[aggr2*nInputs_+kk];
                    vecMarsOut[marsCnt++] = vecSamOuts_[aggr2];
                    vecTags2[aggr2] = 1;
                  }
                }
              }
            }
          }
          if (printLevel_ > 2)
            printf("Metis refine aggr %d (%d): sample size = %d\n",
                   ss+1,nAggrs_,marsCnt);
          if (mode == 3)
               faPtr = genFA(0, nInputs_, iOne, marsCnt);
          else faPtr = genFA(2, nInputs_, iOne, marsCnt);
          faPtr->setBounds(vecLBs.getDVector(), vecUBs.getDVector());
          faPtr->setOutputLevel(1);
          faPtr->initialize(vecMarsIn.getDVector(), vecMarsOut.getDVector());
          double *dPtr = vecSamInps_.getDVector();
          vecDiffs[ss] = faPtr->evaluatePoint(&(dPtr[ss*nInputs_]));
          vecDiffs[ss] = PABS(vecDiffs[ss] - vecSamOuts_[ss]);
          delete faPtr;
        }
        dtmp = 0.0;
        for (ss = 0; ss < nAggrs_; ss++) dtmp += vecDiffs[ss];
        printf("Average interpolation error = %e\n", dtmp/nAggrs_);
      }

      thresh = 0.0;
      ddmax = -PSUADE_UNDEFINED;
      ddmin =  PSUADE_UNDEFINED;
      for (ii = 0; ii < nAggrs_; ii++)
      {
        if (vecDiffs[ii] > ddmax) ddmax = vecDiffs[ii];
        if (vecDiffs[ii] < ddmin) ddmin = vecDiffs[ii];
      }
      thresh = 0.0;
      if (psSamExpertMode_ == 1 && pruneSet_ == 0)
      {
        printf("Maximum discrepancy between aggregates = %e\n",ddmax);
        printf("Minimum discrepancy between aggregates = %e\n",ddmin);
        fp = fopen("arsmPDF.m","w");
        // add check for NULL by Bill Oliver
        if(fp != NULL)
        {
          fprintf(fp,"A = [\n");
          for (ii = 0; ii < nAggrs_; ii++)
            fprintf(fp,"%e\n",vecDiffs[ii]);
          fprintf(fp,"];\n");
          fprintf(fp,"hist(A,10)\n");
          fclose(fp);
        }
        printf("A PDF for error has been given to you in arsmPDF.m.\n");
        printf("Default threshold for pruning = 0.\n");
        sprintf(pString,"Enter thresh to prune refinement candidates : ");
        thresh = getDouble(pString);
        if (thresh < 0.0)
        {
          printf("WARNING: threshold < 0 not allowed. Reset to 0.\n");
          thresh = 0.0;
        }
        pruneSet_ = 1;
      }
      for (ii = 0; ii < nAggrs_; ii++)
        if (vecDiffs[ii] < thresh) vecDiffs[ii] = 0.0;

      vecSortList2.setLength(nAggrs_);
      for (ii = 0; ii < nAggrs_; ii++) vecSortList2[ii] = (double) ii;
      sortDbleList2(nAggrs_,vecDiffs.getDVector(),vecSortList2.getDVector());
      cellCnt = 0;
      ddmax = -PSUADE_UNDEFINED;
      for (ii = nAggrs_-1; ii >= 0; ii--) vecRefineFlags[ii] = 0;
      for (ii = nAggrs_-1; ii >= 0; ii--)
      {
        jj = (int) vecSortList2[ii];
        localN = vecAggrCnts_[jj];
        if (localN > 1 && vecDiffs[ii] > 0.0)
        {
          vecRefineFlags[jj] = 1;
          maxCellSize = (localN > maxCellSize) ? localN : maxCellSize;
          minCellSize = (localN < maxCellSize) ? localN : minCellSize;
          avgCellSize += localN;
          if (vecDiffs[ii] > ddmax) ddmax = vecDiffs[ii];
          cellCnt++;
        }
        if (cellCnt >= refineSize_) break;
      }
      //printf("MetisSampling:: maximum neighbor discrepancy = %e\n",ddmax);
    }
    printf("MetisSampling::refine - number of new sample points = %d.\n",
           cellCnt);
  }

  options[0] = 0;
  currNAggr = nAggrs_;
  splitCount = splitSuccess = 0;
  for (ss = 0; ss < nAggrs_; ss++)
  {
    localN = vecAggrCnts_[ss];
    if (localN > maxN) 
    {
      printf("MetisSampling INTERNAL ERROR (1): consult PSUADE developers\n");
      exit(1);
    }
    if (vecRefineFlags.length() == 0 || vecRefineFlags[ss] == 1) splitCount++;
    if (localN == 1 && (vecRefineFlags.length() == 0 || 
        vecRefineFlags[ss] == 1))
      printf("MetisSampling::refine INFO- cannot split cell (too small).\n");
    if (localN > 1 && 
        (vecRefineFlags.length() == 0 || vecRefineFlags[ss] == 1))
    {
      if (printLevel_ > 4)
        printf("MetisSampling::refine - split sample %d\n",ss+1);
      splitSuccess++;
      localNNZ = 0;
      vecSubLabels = vecAggrLabels_[ss];
      vecLocalIA[0] = localNNZ;
      for (ii = 0; ii < localN; ii++)
      {
        rowInd = vecSubLabels[ii];
        for (jj = vecGraphI_[rowInd]; jj < vecGraphI_[rowInd+1]; jj++)
        {
          colInd = vecGraphJ_[jj];
          status = binarySearchInt(colInd,vecSubLabels.getIVector(),localN);
          if (status >= 0) vecLocalJA[localNNZ++] = status;
        }
        if (ii >= maxN)
        {
          printf("MetisSampling INTERNAL ERROR (2): consult developers\n");
          exit(1);
        }
        vecLocalIA[ii+1] = localNNZ;
      }
      if (localNNZ > maxNNZ) 
      {
        printf("MetisSampling INTERNAL ERROR (3): consult developers\n");
        exit(1);
      }
 
#ifdef HAVE_METIS
      METIS_PartGraphRecursive(&localN, vecLocalIA.getIVector(), 
           vecLocalJA.getIVector(), NULL, NULL, &wgtflag,&numflag,&itwo,
           options,&edgeCut,vecSubLabels.getIVector());
#else
      printf("MetisSampling ERROR : METIS not installed.\n");
      exit(1);
#endif

      count0 = 0;
      for (ii = 0; ii < localN; ii++) if (vecSubLabels[ii] == 0) count0++;
      count1 = localN - count0;
 
      for (ii = 0; ii < localN; ii++) 
      {
        index = vecAggrLabels_[ss][ii];
        if (vecCellsOccupied_[index] < 0) break;
      }
      newCell = 0;
      if (vecLabels[ii] == 0) newCell = 1;
 
      if (newCell == 0)
      {
        vecAggrCnts_[currNAggr] = count0;
        vecAggrLabels_[currNAggr].setLength(count0); 
        count = 0;
        for (ii = 0; ii < localN; ii++) 
        {
          if (vecLabels[ii] == 0) 
          {
            vecAggrLabels_[currNAggr][count++] = vecAggrLabels_[ss][ii];
            vecCellsOccupied_[vecAggrLabels_[ss][ii]] = currNAggr;
          }
        }
        if (count != count0)
        {
          printf("MetisSampling INTERNAL ERROR (4): consult developers\n");
          exit(1);
        }
        vecAggrCnts_[ss] = count1; 
        count = 0;
        for (ii = 0; ii < localN; ii++) 
          if (vecLabels[ii] == 1) 
            vecAggrLabels_[ss][count++] = vecAggrLabels_[ss][ii];
        if (count != count1)
        {
          printf("MetisSampling INTERNAL ERROR (5): consult developers\n");
          exit(1);
        }
      }
      else
      {
        vecAggrCnts_[currNAggr] = count1;
        vecAggrLabels_[currNAggr].setLength(count1); 
        count = 0;
        for (ii = 0; ii < localN; ii++) 
        {
          if (vecSubLabels[ii] == 1) 
          {
            vecAggrLabels_[currNAggr][count++] = vecAggrLabels_[ss][ii];
            vecCellsOccupied_[vecAggrLabels_[ss][ii]] = currNAggr;
          }
        }
        if (count != count1)
        {
          printf("MetisSampling INTERNAL ERROR (6): consult developers\n");
          exit(1);
        }
        vecAggrCnts_[ss] = count0; 
        count = 0;
        for (ii = 0; ii < localN; ii++) 
          if (vecLabels[ii] == 0) 
            vecAggrLabels_[ss][count++] = vecAggrLabels_[ss][ii];
        if (count != count0)
        {
          printf("MetisSampling INTERNAL ERROR (7): consult developers\n");
          exit(1);
        }
      }
      currNAggr++;
    }
  }
  if (printLevel_ > 4 && refineType_ == 2 && maxCellSize != 0)
    printf("MetisSampling:: maximum resolution = %d \n", maxCellSize);
  if (printLevel_ > 4 && refineType_ == 2 && minCellSize != 1000000000)
    printf("MetisSampling:: minimum resolution = %d \n", minCellSize);
  if (printLevel_ > 4 && refineType_ == 2 && cellCnt != 0)
    printf("MetisSampling:: average resolution = %e (%d)\n", 
           1.0*avgCellSize/cellCnt, cellCnt);
           
  if (printLevel_ > 4 && splitSuccess != splitCount)
    printf("MetisSampling:: number of successful splits = %d (out of %d)\n",
           splitSuccess, splitCount);

  oldNumSamples = nSamples_;
  vecOldSamInps = vecSamInps_;
  vecOldSamOuts = vecSamOuts_;
  vecOldSamStas = vecSamStas_;
  nSamples_ = currNAggr;
  allocSampleData();
  for (ss = 0; ss < oldNumSamples; ss++)
  {
    for (inputID = 0; inputID < nInputs_; inputID++)
      vecSamInps_[ss*nInputs_+inputID] = vecOldSamInps[ss*nInputs_+inputID];
    for (outputID = 0; outputID < nOutputs_; outputID++)
      vecSamOuts_[ss*nOutputs_+outputID] =
                 vecOldSamOuts[ss*nOutputs_+outputID];
    vecSamStas_[ss] = vecOldSamStas[ss];
  }

  for (ss = oldNumSamples; ss < currNAggr; ss++)
  {
    if ((randFlag & 1) == 1)
         index = (int) (PSUADE_drand() * vecAggrCnts_[ss]);
    else index = vecAggrCnts_[ss] / 2;
    if (index == vecAggrCnts_[ss]) index--;
    index = vecAggrLabels_[ss][index];
    vecCellsOccupied_[index] = -(vecCellsOccupied_[index] + 1);
    itmp = index;
    for (inputID = 0; inputID < nInputs_; inputID++)
    {
      jtmp = itmp % n1d_;
      itmp = itmp / n1d_;
      if ((randFlag & 1) == 1)
      dtmp = ((double) (jtmp + 0.999*PSUADE_drand())) / (double) n1d_;
      else dtmp = ((double) (jtmp + 0.5)) / (double) n1d_;
      vecSamInps_[ss*nInputs_+inputID] = dtmp * vecRanges[inputID] +
                                         vecLBs[inputID];
      if (vecSamInps_[ss*nInputs_+inputID] < vecLBs_[inputID])
        vecSamInps_[ss*nInputs_+inputID] = vecLBs_[inputID];
      if (vecSamInps_[ss*nInputs_+inputID] > vecUBs_[inputID])
        vecSamInps_[ss*nInputs_+inputID] = vecUBs_[inputID];
    }
  }
  nAggrs_ = currNAggr;

  if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "w");
  else                      fp = fopen("psuadeMetisInfo.tmp", "w");
  if (fp != NULL)
  {
    fprintf(fp, "%d %d %d\n", nAggrs_, nSamples_, nInputs_);
    for (ii = 0; ii < nAggrs_; ii++)
    {
      fprintf(fp, "%d\n", vecAggrCnts_[ii]);
      for (jj = 0; jj < vecAggrCnts_[ii]; jj++)
      {
        fprintf(fp, "%d ", vecAggrLabels_[ii][jj]);
        if (jj % 10 == 0) fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
    }
    jj = 0;
    for (ii = 0; ii < graphN_; ii++) if (vecCellsOccupied_[ii] < 0) jj++;
    fprintf(fp, "%d\n", jj);
    for (ii = 0; ii < graphN_; ii++)
      if (vecCellsOccupied_[ii] < 0) 
        fprintf(fp, "%7d %7d\n",ii,-(vecCellsOccupied_[ii]+1));
    fclose(fp);
  }

  if (printLevel_ > 0)
  {
    printEquals(PL_INFO, 0);
    printf("MetisSampling::refine(2): nSamples = %d\n", nSamples_);
    printf("MetisSampling::refine(2): nInputs  = %d\n", nInputs_);
    printf("MetisSampling::refine(2): nOutputs = %d\n", nOutputs_);
    if (randomize_ != 0)
         printf("MetisSampling::refine: randomize on\n");
    else printf("MetisSampling::refine: randomize off\n");
    printEquals(PL_INFO, 0);
  }
  return 0;
}

// ************************************************************************
// set internal scheme
// ------------------------------------------------------------------------
int MetisSampling::setParam(char * sparam)
{
  int  ii, curVol,count, sID, localN, sID2, rowInd;
  char winput[501];
  FILE *fp;
  psIVector vecSubLabels, vecGridFlags;

  sscanf(sparam, "%s", winput);
  if (!strcmp(winput, "reset"))
  {
    if (changeInfoName_ == 0) fp = fopen("psuadeMetisInfo", "r");
    else                      fp = fopen("psuadeMetisInfo.tmp", "r");
    if (fp != NULL)
    {
      fclose(fp);
      if (changeInfoName_ == 0) unlink("psuadeMetisInfo");
      else                      unlink("psuadeMetisInfo.tmp");
    }
  }
  else if (!strcmp(winput, "changeInfoName"))
  {
    changeInfoName_ = 1;
  }
  else if (!strcmp(winput, "setUniformRefinement"))
  {
    refineType_ = 0;
  }
  else if (!strcmp(winput, "setAdaptiveRefinementBasedOnErrors"))
  {
    refineType_ = 1;
  }
  else if (!strcmp(winput, "setAdaptiveRefinementBasedOnOutputs"))
  {
    refineType_ = 2;
  }
  else if (!strcmp(winput, "calVolumes"))
  {
    curVol = count = 0;
    for (ii = 0; ii < nAggrs_; ii++)
    {
      if (vecSamOuts_[ii] == 1)
      {
        curVol += vecAggrCnts_[ii];
        count++;
      }
    }
    return curVol;
  }
  else if (!strcmp(winput, "totalVolumes"))
  {
    curVol = 0;
    for (ii = 0; ii < nAggrs_; ii++) curVol += vecAggrCnts_[ii];
    return curVol;
  }
  else if (!strcmp(winput, "setRefineSize"))
  {
    sscanf(sparam, "%s %d", winput, &refineSize_);
  }
  else if (!strcmp(winput, "genMeshPlot") && nInputs_ == 2)
  {
    vecGridFlags.setLength(graphN_);
    for (sID = 0; sID < nAggrs_; sID++)
    {
      localN = vecAggrCnts_[sID];
      vecSubLabels = vecAggrLabels_[sID];
      for (ii = 0; ii < localN; ii++)
      {
        rowInd = vecSubLabels[ii];
        if (vecSamOuts_[sID] == 1) vecGridFlags[rowInd] = 1;
        else                       vecGridFlags[rowInd] = 0;
      }
    }
    fp = fopen("metisMeshPlot.m", "w");
    if (fp != NULL)
    {
      fprintf(fp, "meshA = [ \n");
      for (sID = 0; sID < n1d_; sID++)
      {
        for (sID2 = 0; sID2 < n1d_; sID2++)
          fprintf(fp, "%d ", vecGridFlags[sID*n1d_+sID2]);
        fprintf(fp, "\n");
      }
      fprintf(fp, "];\n");
      fprintf(fp, "contour(meshA)\n");
      fclose(fp);
    }
    return 0;
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MetisSampling& MetisSampling::operator=(const MetisSampling & ms)
{
  if (this == &ms) return *this;
  nSamples_ = ms.nSamples_;
  nInputs_ = ms.nInputs_;
  nOutputs_ = ms.nOutputs_;
  refineType_ = ms.refineType_;
  refineSize_ = ms.refineSize_;
  n1d_ = ms.n1d_;
  nAggrs_ = ms.nAggrs_;
  graphN_ = ms.graphN_;
  if (ms.vecAggrLabels_ != NULL)
  {
    vecAggrLabels_ = new psIVector[nAggrs_];
    for (int ii = 0; ii < nAggrs_; ii++)
      vecAggrLabels_[ii] = ms.vecAggrLabels_[ii]; 
  }
  vecAggrCnts_ = ms.vecAggrCnts_;
  printLevel_ = ms.printLevel_;
  samplingID_ = ms.samplingID_;
  randomize_ = ms.randomize_;
  nReplications_ = ms.nReplications_;
  vecLBs_ = ms.vecLBs_;
  vecUBs_ = ms.vecUBs_;
  vecSamInps_ = ms.vecSamInps_;
  vecSamOuts_ = ms.vecSamOuts_;
  vecSamStas_ = ms.vecSamStas_;
  vecGraphI_ = ms.vecGraphI_;
  vecGraphJ_ = ms.vecGraphJ_;
  vecCellsOccupied_ = ms.vecCellsOccupied_;
  return (*this);
}

