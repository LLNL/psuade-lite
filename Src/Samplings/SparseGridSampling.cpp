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
// Functions for the sparse grid sampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <sstream>
#include "sysdef.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "Vector.h"
#include "SparseGridSampling.h"

#define HUQ_nLevels_ 8
#define HUQ_nTerms_  4
static double 
HUQ_GL[HUQ_nLevels_][HUQ_nTerms_] =
{
  {0.5, 0.0, 0.0, 0.0},
  {7.8867513459481287e-1, 0.0, 0.0, 0.0},
  {5.0e-01, 8.8729833462074170e-1, 0.0, 0.0},
  {6.6999052179242813e-1, 9.3056815579702623e-1, 0.0, 0.0},
  {0.5, 7.6923465505284150e-1, 9.5308992296933193e-1, 0.0},
  {6.1930959304159849e-1, 8.3060469323313235e-1, 9.6623475710157603e-1, 0},
  {0.5, 7.0292257568869854e-1, 8.7076559279969723e-1, 9.7455395617137919e-1},
  {5.9171732124782495e-1, 7.6276620495816450e-1, 8.9833323870681348e-1, 9.8014492824876809e-1}

};
static double HUQ_GLW[HUQ_nLevels_][HUQ_nTerms_] =
{
  {1.0, 0.0, 0.0, 0.0},
  {0.5, 0.0, 0.0, 0.0},
  {4.4444444444444570e-1, 2.7777777777777712e-1, 0, 0},
  {3.2607257743127516e-1, 1.7392742256872484e-1, 0, 0},
  {2.8444444444444655e-1, 2.3931433524968501e-1, 1.1846344252809174e-1, 0},
  {2.3395696728634746e-1, 1.8038078652407072e-1, 8.5662246189581834e-2, 0},
  {2.089795918367362e-1, 1.909150252525609e-1, 1.3985269574463935e-1, 6.4742483084431701e-2},
  {1.8134189168918213e-1, 1.5685332293894469e-1, 1.1119051722668793e-1, 5.0614268145185180e-2}
};
static int HUQ_GLn[HUQ_nLevels_] =
{
  1, 1, 2, 2, 3, 3, 4, 4
};

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
SparseGridSampling::SparseGridSampling() : Sampling()
{
  samplingID_ = PSUADE_SAMP_SG;
  pOrder_ = 2; /* default order = 2 */
}

// ************************************************************************
// copy constructor added by Oliver
// ------------------------------------------------------------------------
SparseGridSampling::SparseGridSampling(const SparseGridSampling &gs) : 
                                 Sampling()
{
  pOrder_ = gs.pOrder_;
  nSamples_ = gs.nSamples_;
  vecSamWeights_ = gs.vecSamWeights_;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SparseGridSampling::~SparseGridSampling()
{
}

//*************************************************************************
//* initialize the sampling data
//*------------------------------------------------------------------------
int SparseGridSampling::initialize(int initLevel)
{
  int    ii, jj, kk, ll, minQ, maxQ, bQ, index, size = 0;
  int    nPerms, **pcePerms, total, newLeng, nVecs, numNew;
  double val;
  FILE   *fp;
  psVector  *Vnodes, Vweights, vecRanges, vecNewW;
  psIVector vecCounts, vecMidx, vecFlags;
  psMatrix  matDT, matNewN;
  psIMatrix matPCEPerms;

  if (nSamples_ == 0)
  {
    printf("SparseGridSampling::initialize ERROR - nSamples = 0.\n");
    exit(1);
  }
  if (nInputs_ == 0)
  {
    printf("SparseGridSampling::initialize ERROR - input not set up.\n");
    exit(1);
  }

  if (initLevel != 0) return 0;

  if (psSamExpertMode_ == 1 && isScreenDumpModeOn() == 1)
  {
    pOrder_ = 0;
    while (pOrder_ < 2 || pOrder_ > 4)
    {
      printf("SparseGridSampling: enter polynomial order (2 - 4): ");
      scanf("%d", &pOrder_);
    }
  }
  else
  {
    pOrder_ = 2;
    printf("SparseGridSampling: polynomial order has been set to 2.\n");
  } 

  minQ = pOrder_ + 1 - nInputs_;
  if (minQ < 0) minQ = 0;
  maxQ = pOrder_;
  vecMidx.setLength(nInputs_);
  Vnodes = new psVector[nInputs_];
  nVecs = 0;
  for (ii = minQ; ii <= maxQ; ii++)
  {
    val = (pow(-1.0, 1.0*(maxQ-ii))) * 
           nChooseK(nInputs_-1,nInputs_+ii-(pOrder_+1));
    bQ = (int) val;
    nPerms = GenSequence(nInputs_, nInputs_+ii, matPCEPerms);
    pcePerms = matPCEPerms.getMatrix2D();
    vecCounts.setLength(nPerms);
    total = 0;
    for (jj = 0; jj < nPerms; jj++)
    {
      vecCounts[jj] = 1;
      for (kk = 0; kk < nInputs_; kk++)
      {
        index = pcePerms[jj][kk];
        vecCounts[jj] *= HUQ_GLn[index];
      } 
      total += vecCounts[jj];
    }
    for (kk = 0; kk < nInputs_; kk++) Vnodes[kk].addElements(total, NULL);
    Vweights.addElements(total,NULL);
    for (jj = 0; jj < nPerms; jj++)
    {
      for (kk = 0; kk < nInputs_; kk++)
      {
        index = pcePerms[jj][kk];
        vecMidx[kk] = index;
      }

      size = KronProd(nInputs_, vecMidx, matNewN, vecNewW);
      double **newn = matNewN.getMatrix2D();
      if (size <= 0)
      {
        printf("size variable is <= 0 in file %s line %d\n",
               __FILE__, __LINE__);
        exit(1);
      }
      for (ll = nVecs; ll < nVecs+vecCounts[jj]; ll++)
      {
        for (kk = 0; kk < nInputs_; kk++)
          Vnodes[kk][ll] = newn[ll-nVecs][kk];
        Vweights[ll] = bQ * vecNewW[ll-nVecs];
      } 
      nVecs += vecCounts[jj];
    }
    /* need to prune repeated ones */
    matDT.setFormat(2);
    matDT.setDim(nVecs, nInputs_+1);
    for (jj = 0; jj < nVecs; jj++) 
    {
      for (kk = 0; kk < nInputs_; kk++) 
        matDT.setEntry(jj, kk, Vnodes[kk][jj]);
      matDT.setEntry(jj, nInputs_, Vweights[jj]);
    }
    newLeng = sortAndDelete(nVecs, nInputs_+1, matDT.getMatrix2D());
    for (jj = 0; jj < newLeng; jj++) 
    {
      for (kk = 0; kk < nInputs_; kk++) 
        Vnodes[kk][jj] = matDT.getEntry(jj, kk);
      Vweights[jj] = matDT.getEntry(jj, nInputs_);
    }
    nVecs = newLeng;
  }

  val = HUQ_GL[0][0];
  for (ii = 0; ii < nInputs_; ii++)
  {
    vecFlags.setLength(nVecs);
    for (jj = 0; jj < nVecs; jj++) vecFlags[jj] = -1;
    numNew = 0;
    for (jj = 0; jj < nVecs; jj++) 
    {
      if (Vnodes[ii][jj] != val)
      {
        vecFlags[numNew] = jj;
        numNew++;
      }
    }
    if (numNew > 0)
    {
      for (jj = 0; jj < nInputs_; jj++) 
        Vnodes[jj].addElements(numNew,NULL);
      Vweights.addElements(numNew,NULL);
      for (jj = 0; jj < nInputs_; jj++)
      {
        numNew = 0;
        for (kk = 0; kk < nVecs; kk++) 
        {
          if (vecFlags[kk] != -1)
          {
            Vnodes[jj][nVecs+numNew] = Vnodes[jj][vecFlags[kk]];
            Vweights[nVecs+numNew] = Vweights[vecFlags[kk]];
            numNew++;
          }
        }
      }
      for (kk = nVecs; kk < nVecs+numNew; kk++) 
        Vnodes[ii][kk] = 2 * val - Vnodes[ii][kk];
      nVecs += numNew;
    }
  }

  /* need to prune repeated ones */
  matDT.setFormat(2);
  matDT.setDim(nVecs, nInputs_+1);
  for (jj = 0; jj < nVecs; jj++) 
  {
    for (kk = 0; kk < nInputs_; kk++) 
      matDT.setEntry(jj, kk, Vnodes[kk][jj]);
    matDT.setEntry(jj, nInputs_, Vweights[jj]);
  }
  newLeng = sortAndDelete(nVecs, nInputs_+1, matDT.getMatrix2D());
  for (jj = 0; jj < newLeng; jj++) 
  {
    for (kk = 0; kk < nInputs_; kk++) 
      Vnodes[kk][jj] = matDT.getEntry(jj,kk);
    Vweights[jj] = matDT.getEntry(jj,nInputs_);
  }
  nVecs = newLeng;
  val = 0.0;
  for (jj = 0; jj < nVecs; jj++) val += Vweights[jj];
  for (jj = 0; jj < nVecs; jj++) Vweights[jj] /= val;

  nSamples_ = nVecs;
  allocSampleData();
  vecRanges.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++) 
    vecRanges[ii] = vecUBs_[ii] - vecLBs_[ii];
  vecSamWeights_.setLength(nSamples_);
  for (kk = 0; kk < nSamples_; kk++)
  {
    for (ii = 0; ii < nInputs_; ii++) 
      vecSamInps_[kk*nInputs_+ii] = 
                 Vnodes[ii][kk] * vecRanges[ii] + vecLBs_[ii];
  }
   
  fp = fopen("ps_sparse_grid_info", "w");
  fprintf(fp, "%d %d %d\n", nSamples_, nInputs_, pOrder_);
  for (kk = 0; kk < nSamples_; kk++)
  {
    for (ii = 0; ii < nInputs_; ii++) 
      fprintf(fp, "%24.16e ", Vnodes[ii][kk]);
    vecSamWeights_[kk] = Vweights[kk]; 
    fprintf(fp, "%24.16e ", vecSamWeights_[kk]);
    fprintf(fp, "\n");
  }
  for (ii = 0; ii < nInputs_; ii++) 
    fprintf(fp, "%24.16e %24.16e\n", vecLBs_[ii], vecUBs_[ii]);
  fclose(fp);
  printf("Sparse grid data has been stored to ps_sparse_grid_info.\n");
  printf("You need this file to build sparse grid response surfaces.\n");

  delete [] Vnodes;
  return 0;
}

// ************************************************************************
// refine the sample space
// ------------------------------------------------------------------------
int SparseGridSampling::refine(int, int, double , int, double *)
{
  printf("SparseGridSampling::refine ERROR - not available.\n");
  exit(1);
  return 0;
}

// ************************************************************************
// Kronecker product 
// ------------------------------------------------------------------------
int SparseGridSampling::KronProd(int nn, psIVector &vecMidx,
                                 psMatrix &matNewN, psVector &vecNewW)
{
  int total, ii, index, cnt;
  psIVector vecCounts;

  /* total number of new rows */
  total = 1;
  for (ii = 0; ii < nn; ii++)
  {
    index = vecMidx[ii];
    total *= HUQ_GLn[index];
  }

  /* create storage */
  matNewN.setFormat(2);
  matNewN.setDim(total, nn);
  vecNewW.setLength(total);
  for (ii = 0; ii < total; ii++) vecNewW[ii] = 1.0;
  vecCounts.setLength(nn);
  for (ii = 0; ii < nn; ii++) vecCounts[ii] = 0;

  /* create data */
  double **newn = matNewN.getMatrix2D();
  cnt  = 0;
  while (cnt < total)
  {
    for (ii = 0; ii < nn; ii++) 
    {
      newn[cnt][ii] = HUQ_GL[vecMidx[ii]][vecCounts[ii]];
      vecNewW[cnt] *= HUQ_GLW[vecMidx[ii]][vecCounts[ii]];
    }
    vecCounts[nn-1]++;
    ii = nn - 1; 
    while (vecCounts[ii] >= vecMidx[ii] && ii > 0)
    {
      vecCounts[ii-1]++;
      vecCounts[ii] = 0;
      ii--;
    }
    cnt++;
  }
  /* This function is used in a loop so we need to return total */
  return total;
}

// ************************************************************************
// generate sequence
// ------------------------------------------------------------------------
int SparseGridSampling::GenSequence(int nn, int rsum, psIMatrix &matPerms)
{
  int nPerms, cnt, ii, jj, idata, idata2;
  psIVector vecFlags;

  vecFlags.setLength(nn);
  for (ii = 0; ii < nn; ii++) vecFlags[ii] = 0;
  idata = rsum - nn;
  vecFlags[0] = idata;
  if      ((rsum-nn) == 0) nPerms = 1;
  else if ((rsum-nn) == 1) nPerms = nn;
  else                     
  {
    nPerms = computeNumPCEPermutations(nn, idata) - 
             computeNumPCEPermutations(nn, idata-1);
  }
  if(nPerms <= 0)
  {
    printf("Problem in SparseGridSampling::GenSequence with calculation\n");
    printf("of nPerms (%d).\n", nPerms);
    exit(1);
  }

  matPerms.setFormat(2);
  matPerms.setDim(nPerms, nn);
  int **outPerms = matPerms.getMatrix2D();
  for (ii = 0; ii < nn; ii++) outPerms[0][ii] = vecFlags[ii];
  idata2 = 0;
  cnt = 1;
  while (vecFlags[nn-1] < idata)
  {
    if (idata2 == nn-1)
    {
      for (jj = idata2-1; jj >= 0; jj--)
      {
        idata2 = jj;
        if (vecFlags[jj] != 0) break;
      }
    }
    vecFlags[idata2]--;
    idata2++;
    vecFlags[idata2] = idata;
    for (jj = 0; jj < idata2; jj++) vecFlags[idata2] -= vecFlags[jj];
    if (idata2 < nn)
      for (jj = idata2+1; jj < nn; jj++) vecFlags[jj] = 0;
    for (jj = 0; jj < nn; jj++) outPerms[cnt][jj] = vecFlags[jj];
    cnt++;
  }
  return nPerms;
}

// ************************************************************************
// n choose k
// ------------------------------------------------------------------------
int SparseGridSampling::nChooseK(int n, int k)
{
  int ii, idata=1;
  for (ii = n; ii > n-k; ii--) idata = idata * ii / (n - ii + 1);
  return idata;
}

// ************************************************************************
// set input settings
// ------------------------------------------------------------------------
int SparseGridSampling::setParam(char *sparam)
{
  char winput[1001];
  sscanf(sparam, "%s", winput);
  if (!strcmp(winput, "pOrder"))
  {
    sscanf(sparam, "%s %d", winput, &pOrder_);
    if (pOrder_ < 2) pOrder_ = 2;
    if (pOrder_ > 4) pOrder_ = 4;
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
SparseGridSampling& SparseGridSampling::operator=(const SparseGridSampling &gs)
{
  // Bill Oliver modified the operator= to work
  if (this == &gs) return *this;
  pOrder_ = gs.pOrder_;
  nSamples_ = gs.nSamples_;
  vecSamWeights_ = gs.vecSamWeights_;
  return (*this);
}

