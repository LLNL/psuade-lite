// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team. 
// All rights reserved.
//
// Please see the COPYRIGHT and LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class DeltaAnalyzer (Delta test)
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG/Michael Snow
// DATE   : 2009
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "DeltaAnalyzer.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
DeltaAnalyzer::DeltaAnalyzer(): Analyzer(),nBins_(0),nInputs_(0),nConfig_(0),
                   minDeltas_(0), deltaBins_(0), dOrder_(0), ranks_(0)
{
   setName("DELTATEST");
   mode_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
DeltaAnalyzer::~DeltaAnalyzer()
{
   if (minDeltas_) delete []minDeltas_;
   if (deltaBins_)
   {
      for (int jj = 0; jj < nBins_; jj++) delete [] deltaBins_[jj];
      delete [] deltaBins_;
   }
   if (dOrder_) delete []dOrder_;
   if (ranks_) delete []ranks_;
}

// ************************************************************************
// perform delta test
// ------------------------------------------------------------------------
double DeltaAnalyzer::analyze(aData &adata)
{
   int    printLevel, nSamples, nInputs, nOutputs, outputID, ss, ss2;
   int    *inputBins, *auxBins, nSelected=0, *minIndices, nIndex=3, info;
   int    ii, jj, kk, ll, **deltaBins, nBins=1000, *iPtr, converged=0;
   int    place, count, uniqueFlag=0, reverseCnt=0, *ranks;
   int    nConfig=20, iter=100, stagnate=100;
   double *X, *Y, distance, delta, minDist, *minDeltas, *iLowerB, *iUpperB;
   double dtemp, temperature=.0001, oldDelta=PSUADE_UNDEFINED, minDelta;
   double bestDelta=PSUADE_UNDEFINED, *dOrder, *rangesInv2, *distPairs;
   double *YY, alpha=0.98, r, ddata, accum, auxMin, deltaSave;
   char   pString[500];
   FILE   *fp;

   printLevel = adata.printLevel_;
   nSamples   = adata.nSamples_;
   X          = adata.sampleInputs_;
   YY         = adata.sampleOutputs_;
   nInputs    = adata.nInputs_;
   nInputs_   = nInputs;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   outputID   = adata.outputID_;
   iLowerB    = adata.iLowerB_;
   iUpperB    = adata.iUpperB_;
   nBins_     = nBins;
   nConfig_   = nConfig;

   if (adata.inputPDFs_ != NULL)
   {
      count = 0;
      for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
      if (count > 0)
      {
         printOutTS(PL_WARN, 
              "DeltaTest INFO: some inputs have non-uniform PDFs, but\n");
         printOutTS(PL_WARN,
              "          they are not relevant in this analysis.\n");
      }
   }

   if (nSamples <= 1)
   {
      printOutTS(PL_ERROR, 
           "DeltaTest INFO: test not meaningful for nSamples <= 1.\n");
      return PSUADE_UNDEFINED;
   }
   if (X == NULL || YY == NULL)
   {
      printOutTS(PL_ERROR, "DeltaTest ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   }
   info = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (YY[nOutputs*ii+outputID] == PSUADE_UNDEFINED) info = 1;
   if (info == 1)
   {
      printOutTS(PL_ERROR, 
           "DeltaTest ERROR: Some outputs are undefined.\n");
      printOutTS(PL_ERROR, 
           "                 Prune the undefined's first.\n");
      return PSUADE_UNDEFINED;
   }

   if (minDeltas_) delete [] minDeltas_;
   if (deltaBins_)
   {
      for (jj = 0; jj < nBins_; jj++) delete [] deltaBins_[jj];
      delete [] deltaBins_;
   }
   if (dOrder_) delete []dOrder_;
   if (ranks_)  delete []ranks_;
   ranks_ = NULL;
   dOrder_ = NULL;
   deltaBins_ = NULL;
   minDeltas_ = NULL;
 
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"DeltaTest for variable selection\n");
   printOutTS(PL_INFO,
        "This test has the characteristics that the more important\n");
   printOutTS(PL_INFO,
        "a parameter is relative to the others, the smaller the \n");
   printOutTS(PL_INFO,
        "subset is at the end of the test (sharp zoom into the most\n");
   printOutTS(PL_INFO,"important subset).\n");
   printOutTS(PL_INFO,
        "Thus, the purpose of this test is to identify a subset of\n");
   printOutTS(PL_INFO,"important parameters.\n");
   printOutTS(PL_INFO,
        "Note: If both nInputs and nSamples are large, this test\n");
   printOutTS(PL_INFO,
        "      may take a long time to run. So, be patient.)\n");
   printEquals(PL_INFO, 0);
   auxBins = new int[nInputs];
   inputBins = new int[nInputs];
   ranks = new int[nInputs];
   ranks_ = new int[nInputs_];
   dOrder = new double[nInputs];
   dOrder_ = new double[nInputs_];
   rangesInv2 = new double[nInputs];
   distPairs = new double[nSamples*(nSamples-1)/2];
   Y = new double[nSamples];
   checkAllocate(Y, "Y in DeltaTest::analyze");
   for (ii = 0; ii < nSamples; ii++) Y[ii] = YY[ii*nOutputs+outputID];

   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO,
           "DeltaTest Option: set the number of neighbors K.\n");
      printOutTS(PL_INFO, 
           "The larger K is, the larger the distinguishing power is.\n");
      sprintf(pString, "What is K (>= 1, <= 20, default=3)? ");
      nIndex = getInt(1, 20, pString);
      sprintf(pString,"How many inputs to select FOR SURE? (0 if not sure) ");
      nSelected = getInt(0, nInputs-1, pString);
      for (ii = 0; ii < nInputs; ii++) auxBins[ii] = 0;
      for (ii = 0; ii < nSelected; ii++)
      {
         sprintf(pString,"Enter the %d-th input to be selected : ", ii+1);
         kk = getInt(1, nInputs, pString);
         auxBins[kk-1] = 1;
      }
      sprintf(pString,"How many iterations for optimization? (> 100) ");
      iter = getInt(1, 100000, pString);
      //sprintf(pString,"How many configurations for ranking? (1 if not sure) ");
      //nConfig = getInt(1, nBins, pString);
      printEquals(PL_INFO, 0);
   }

   for (ii = 0; ii < nInputs; ii++) 
   {
      if ((iUpperB[ii] - iLowerB[ii]) > 0)
         rangesInv2[ii] = 1.0 / (iUpperB[ii] - iLowerB[ii]) /
                                (iUpperB[ii] - iLowerB[ii]);
      else
      {
         printOutTS(PL_ERROR, "DeltaTest ERROR: problem with input range.\n");
         exit(1);
      }
   }
   for (ii = 0; ii < nInputs; ii++) inputBins[ii] = 0;
   deltaBins = new int*[nBins];
   deltaBins_ = new int*[nBins_];
   for (ii = 0; ii < nBins; ii++)
   {
      deltaBins[ii] = new int[nInputs];
      deltaBins_[ii] = new int[nInputs];
      for (jj = 0; jj < nInputs; jj++) deltaBins[ii][jj] = 0;
   }
   minDeltas_ = new double[nBins_];
   minIndices = new int[nIndex];
   minDeltas  = new double[nBins];
   checkAllocate(minDeltas, "minDeltas in DeltaTest::analyze");
   for (ii = 0; ii < nBins; ii++) minDeltas[ii] = PSUADE_UNDEFINED;

   if (nSelected == 0)
      for (ii = 0; ii < nInputs;ii++) inputBins[ii]=PSUADE_rand()%2;
   else
      for (ii = 0; ii < nInputs;ii++) inputBins[ii]=auxBins[ii];

   for (ss = 1; ss < nSamples; ss++)
   {
      for (ss2 = 0; ss2 < ss; ss2++)
      {
         distance = 0.0;
         for (ii = 0; ii < nInputs; ii++)
         {
            if (inputBins[ii] == 1)
            {
               dtemp = X[ss*nInputs+ii] - X[ss2*nInputs+ii];
               distance += dtemp * dtemp * rangesInv2[ii];
            }
         }
         distPairs[ss*(ss-1)/2+ss2] = distance;
      }
   }

   delta = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      ddata = 0.0;
      for (jj = 0; jj < nIndex; jj++)
      {
         minDist = PSUADE_UNDEFINED;
         minIndices[jj] = -1;
         for (ss2 = 0; ss2 < ss; ss2++)
         {
            kk = ss * (ss - 1) / 2 + ss2;
            if (distPairs[kk] < minDist)
            {
               for (ll = 0; ll < jj; ll++)
                  if (ss2 == minIndices[ll]) break;
               if (jj == 0 || ll == jj)
               {
                  minDist = distPairs[kk];
                  minIndices[jj] = ss2;
               }
            }
         }
         for (ss2 = ss+1; ss2 < nSamples; ss2++)
         {
            kk = ss2 * (ss2 - 1) / 2 + ss;
            if (distPairs[kk] < minDist)
            {
               for (ll = 0; ll < jj; ll++)
                  if (ss2 == minIndices[ll]) break;
               if (jj == 0 || ll == jj)
               {
                  minDist = distPairs[kk];
                  minIndices[jj] = ss2;
               }
            }
         }
         if (minIndices[jj] == -1)
         {
            printOutTS(PL_ERROR, "DeltaTest ERROR (1).\n");
            exit(1);
         }
         ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
      }
      delta += ddata / (double) nIndex;
   }
   printOutTS(PL_INFO,"Current best solution for output %d:\n",outputID+1);
   printOutTS(PL_INFO,
        "To stop the search, create a psuade_stop file in local directory.\n");
   printDashes(PL_INFO, 0);
   delta /= (2.0 * nSamples);
   for (ii = 0; ii < nInputs; ii++) printOutTS(PL_INFO, "%d ", inputBins[ii]);
   printOutTS(PL_INFO, " = %e\n", delta);
 
   count = 1;
   auxMin = - PSUADE_UNDEFINED;
   while (count <= 3*iter*nInputs)
   {
      fflush(stdout);
      count++;

      if (reverseCnt >= 4*nInputs)
      {	
         temperature*=nInputs*nInputs;
         for (ii = 0;ii <= nInputs/5;ii++) inputBins[PSUADE_rand()%nInputs] ^=1;
         for (ss = 1; ss < nSamples; ss++)
         {
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               distance = 0.0;
               for (ii = 0; ii < nInputs; ii++)
               {
                  if (inputBins[ii] == 1)
                  {
                     dtemp = X[ss*nInputs+ii] - X[ss2*nInputs+ii];
                     distance += dtemp * dtemp * rangesInv2[ii];
                  }
               }
               distPairs[ss*(ss-1)/2+ss2] = distance;
            }
         }
         reverseCnt = 0;
         place = PSUADE_rand()%(nInputs);
      }
      else 
      {
         if (reverseCnt >= 3*nInputs)
         {
            //printOutTS(PL_WARN, "suspected local minima, checking");
            place = reverseCnt - 3 * nInputs;
         }
         else 
         {
            place = PSUADE_rand()%(nInputs);
         }
      }
      temperature *= alpha;
      inputBins[place] ^= 1;

      delta = 0.0;
      for (ss = 1; ss < nSamples; ss++)
      {
         for (ss2 = 0; ss2 < ss; ss2++)
         {
            kk = ss * (ss - 1) / 2 + ss2;
            dtemp = X[ss*nInputs+place] - X[ss2*nInputs+place];
            if (inputBins[place] == 1)
                 distPairs[kk] += dtemp * dtemp * rangesInv2[place];
            else distPairs[kk] -= dtemp * dtemp * rangesInv2[place];
         }
      }

      delta = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nIndex; jj++)
         {
            minDist = PSUADE_UNDEFINED;
            minIndices[jj] = -1;
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               kk = ss * (ss - 1) / 2 + ss2;
               if (distPairs[kk] < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distPairs[kk];
                     minIndices[jj] = ss2;
                  }
               }
            }
            for (ss2 = ss+1; ss2 < nSamples; ss2++)
            {
               kk = ss2 * (ss2 - 1) / 2 + ss;
               if (distPairs[kk] < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distPairs[kk];
                     minIndices[jj] = ss2;
                  }
               }
            }
            if (minIndices[jj] == -1)
            {
               printOutTS(PL_ERROR, "DeltaTest ERROR (1).\n");
               exit(1);
            }
            ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
         }
         delta += ddata / (double) nIndex;
      }
      delta /= (2.0 * nSamples);
      if ((count % (3*nInputs) == 0))
      {
         for (ii = 0; ii < nInputs; ii++) 
            printOutTS(PL_INFO, "%d ", deltaBins[nBins-1][ii]);
         printOutTS(PL_INFO,
              " = %e (%d of %d)\n", bestDelta, count/(3*nInputs), iter);
      }

      if (delta < minDeltas[0])
      {
         uniqueFlag = 1;
         for (ii = 0; ii < nBins; ii++)
         {
            if (minDeltas[ii] != PSUADE_UNDEFINED)
            {
               for (jj = 0; jj < nInputs; jj++)
                 if (inputBins[jj] != deltaBins[ii][jj]) break;
               if (jj == nInputs) {uniqueFlag = 0; break;}
            }
         }
         if (uniqueFlag == 1)
         {
            minDeltas[0] = delta;
            for (ii = 0; ii < nInputs; ii++) deltaBins[0][ii] = inputBins[ii];
            for (ii = 1; ii < nBins; ii++)
            {
               if (minDeltas[ii] > minDeltas[ii-1])
               {
                  dtemp = minDeltas[ii];
                  minDeltas[ii] = minDeltas[ii-1];
                  minDeltas[ii-1] = dtemp;
                  iPtr = deltaBins[ii];
                  deltaBins[ii] = deltaBins[ii-1];
                  deltaBins[ii-1] = iPtr;
               }
            }
         }
      }

      if (minDeltas[nBins-1] == auxMin) converged++;
      else
      {
         converged = 0;
         auxMin = minDeltas[nBins-1];
      }
      if (converged > stagnate*3*nInputs)
      {
         printOutTS(PL_INFO, "DeltaTest: stagnate for %d iterations, ", 
                    stagnate);
         printOutTS(PL_INFO, "considered converged.\n");
         break;
      }

      if (delta >= oldDelta) 
      {
         r = PSUADE_rand()%100000;
         r /= 100000;

         if (r>=exp(-.1*(delta-oldDelta)/(temperature)))
         {
            inputBins[place] ^=1;
            reverseCnt++;
            for (ss = 1; ss < nSamples; ss++)
            {
               for (ss2 = 0; ss2 < ss; ss2++)
               {
                  kk = ss * (ss - 1) / 2 + ss2;
                  dtemp = X[ss*nInputs+place] - X[ss2*nInputs+place];
                  if (inputBins[place] == 1)
                       distPairs[kk] += dtemp * dtemp * rangesInv2[place];
                  else distPairs[kk] -= dtemp * dtemp * rangesInv2[place];
               }
            }
         }
         else
         {
            oldDelta = delta;
            reverseCnt = 0;
         }
      }
      else 
      {
         oldDelta = delta;
         reverseCnt = 0;
      }
      if (oldDelta <= bestDelta) bestDelta = oldDelta;
      fp = fopen("psuade_stop","r");
      if (fp != NULL)
      {
         printOutTS(PL_INFO, "psuade_stop file found ==> terminate.\n");
         printOutTS(PL_INFO, "To restart, delete psuade_stop first.\n");
         fclose(fp);
         break;
      }
   }

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, 
        "Final Selections (based on %d neighbors) = \n", nIndex);

   //save minDeltas and deltaBins
   for (ii=0; ii < nBins_; ii++)
   {
	   minDeltas_[ii] = minDeltas[ii];
	   for (kk = 0; kk < nInputs_; kk++)
		   deltaBins_[ii][kk] = deltaBins[ii][kk];
   }

   for (kk = 0; kk < 10; kk++)
   {
      if (minDeltas[nBins-kk-1] < 0.99 * PSUADE_UNDEFINED)
      {
         printOutTS(PL_INFO, "Rank %2d => ", kk+1);
         for (ii = 0; ii < nInputs; ii++) 
            printOutTS(PL_INFO, "%d ", deltaBins[nBins-kk-1][ii]);
         printOutTS(PL_INFO, ": delta = %11.4e\n", minDeltas[nBins-kk-1]);
      }
   }
   printDashes(PL_INFO, 0);
   count = 0;
   for (ii = 0; ii < nInputs; ii++)
   {
      ddata = 0;
      accum = 0.0;
      for (kk = 0; kk < nConfig; kk++)
      {
         if (minDeltas[nBins-kk-1] != PSUADE_UNDEFINED)
         {
            ddata += (minDeltas[nBins-1]*deltaBins[nBins-kk-1][ii]/
                      minDeltas[nBins-kk-1]);
            accum += (minDeltas[nBins-1]/minDeltas[nBins-kk-1]);
         }
      }
      ranks[ii] = (int) (ddata / accum * 100);
   }

   if (psPlotTool_ == 1) fp = fopen("scilabdelta.sci", "w");
   else                  fp = fopen("matlabdelta.m", "w");
   if (fp == NULL)
   {
      printOutTS(PL_INFO, "Delta test ERROR: cannot open graphics files.\n");
      printOutTS(PL_INFO, "                  ==> graphics not generated.\n");
   }
   else
   {
      fwritePlotCLF(fp);
      fprintf(fp, "A = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%e\n", 0.01 * ranks[ii]);
      fprintf(fp, "];\n");
      fprintf(fp, "bar(A, 0.8);\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Delta Test Rankings");
      fwritePlotXLabel(fp, "Input parameters");
      fwritePlotYLabel(fp, "Delta Metric (normalized)");
      fclose(fp);
      if (psPlotTool_ == 1) 
           printOutTS(PL_INFO,
              "Delta test ranking is now in scilabdelta.sci.\n");
      else printOutTS(PL_INFO,"Delta test ranking is now in matlabdelta.m.\n");
   }

   for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
   sortIntList2a(nInputs, ranks, dOrder);
   printOutTS(PL_INFO,
        "Order of importance (based on %d best configurations):\n",
        nConfig);
   for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, "(D)Rank %4d : input %4d (score = %d )\n", ii+1, 
                 (int) dOrder[nInputs-ii-1]+1, ranks[nInputs-ii-1]);
   printAsterisks(PL_INFO, 0);

   //save dOrder and ranks
   for (ii = 0; ii < nInputs_; ii++)
   {
	   dOrder_[ii] = dOrder[ii];
	   ranks_[ii] = ranks[ii];
   }
   printOutTS(PL_INFO, 
        "Final test using the most important parameters incrementally:\n");
   printDashes(PL_INFO, 0);
   for (ii = 0; ii < nInputs; ii++) inputBins[ii] = 0;
   for (ii = 1; ii >= 0; ii--)
   {
      inputBins[(int) dOrder[nInputs-ii-1]] = 1;
      delta = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nIndex; jj++)
         {
            minDist = PSUADE_UNDEFINED;
            minIndices[jj] = -1;
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            for (ss2 = ss+1; ss2 < nSamples; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
         }
         delta += ddata / (double) nIndex;
      }
      delta /= (2.0 * nSamples);
      minDeltas[ii] = delta;
   }
   deltaSave = minDeltas[1] - minDeltas[0];
   inputBins[(int) dOrder[nInputs-2]] = 0;
   inputBins[(int) dOrder[nInputs-1]] = 0;
   minDelta = 1.0e35;
   for (ii = 0; ii < nInputs; ii++)
   {
      inputBins[(int) dOrder[nInputs-ii-1]] = 1;
      delta = 0.0;
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = 0.0;
         for (jj = 0; jj < nIndex; jj++)
         {
            minDist = PSUADE_UNDEFINED;
            minIndices[jj] = -1;
            for (ss2 = 0; ss2 < ss; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            for (ss2 = ss+1; ss2 < nSamples; ss2++)
            {
               distance = 0.0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (inputBins[kk] == 1)
                  {
                     dtemp = X[ss*nInputs+kk] - X[ss2*nInputs+kk];
                     distance += dtemp * dtemp * rangesInv2[kk];
                  }
               }
               if (distance < minDist)
               {
                  for (ll = 0; ll < jj; ll++)
                     if (ss2 == minIndices[ll]) break;
                  if (jj == 0 || ll == jj)
                  {
                     minDist = distance;
                     minIndices[jj] = ss2;
                  }
               }
            }
            ddata += pow(Y[ss] - Y[minIndices[jj]], 2.0);
         }
         delta += ddata / (double) nIndex;
      }
      delta /= (2.0 * nSamples);
      if (ii == 0)
      {
         deltaSave += delta;
         for (kk = 0; kk < nInputs; kk++) printOutTS(PL_INFO, "0 ");
         printOutTS(PL_INFO, " = %e\n", deltaSave);
      }
      for (kk = 0; kk < nInputs; kk++) 
         printOutTS(PL_INFO, "%d ", inputBins[kk]);
      printOutTS(PL_INFO, " = %e\n", delta);
      minDeltas[ii] = delta;
      if (delta < minDelta) minDelta = delta;
   }
   printAsterisks(PL_INFO, 0);

   delete [] auxBins;
   delete [] inputBins;
   for (ii = 0; ii < nBins; ii++) delete [] deltaBins[ii];
   delete [] deltaBins;
   delete [] minDeltas;
   delete [] ranks;
   delete [] dOrder;
   delete [] rangesInv2;
   delete [] distPairs;
   delete [] minIndices;
   delete [] Y;
   return minDelta;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int DeltaAnalyzer::setParams(int argc, char **argv)
{
   char *request = (char *) argv[0];
   Analyzer::setParams(argc, argv);
   if (!strcmp(request, "gdelta")) mode_ = 1;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
DeltaAnalyzer& DeltaAnalyzer::operator=(const DeltaAnalyzer &)
{
   printOutTS(PL_ERROR,"DeltaTest operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int DeltaAnalyzer::get_mode()
{
   return mode_;
}
int DeltaAnalyzer::get_nBins()
{
   return nBins_;
}
int DeltaAnalyzer::get_nInputs()
{
   return nInputs_;
}
int DeltaAnalyzer::get_nConfig()
{
   return nConfig_;
}
double *DeltaAnalyzer::get_minDeltas()
{
   double* retVal = NULL;
   if (minDeltas_)
   {
      retVal = new double[nBins_];
      checkAllocate(retVal, "retVal in DeltaTest::get_minDeltas");
      std::copy(minDeltas_, minDeltas_+nBins_+1, retVal);
   }
   return retVal;
}
int **DeltaAnalyzer::get_deltaBins()
{
   int** retVal = NULL;
   if (deltaBins_)
   {
      retVal = new int*[nBins_];
      checkAllocate(retVal, "retVal in DeltaTest::get_deltaBins");
      for (int i=0; i<nBins_; i++)
      {
    	  retVal[i] = new int[nInputs_];
          checkAllocate(retVal[i],"retVal in DeltaTest::get_deltaBins");
    	  std::copy(deltaBins_[i], deltaBins_[i]+nInputs_, retVal[i]);
      }
   }
   return retVal;
}
double *DeltaAnalyzer::get_dOrder()
{
   double* retVal = NULL;
   if (dOrder_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in DeltaTest::get_dOrder");
      std::copy(dOrder_, dOrder_+nInputs_+1, retVal);
   }
   return retVal;
}
int *DeltaAnalyzer::get_ranks()
{
   int* retVal = NULL;
   if (ranks_)
   {
      retVal = new int[nInputs_];
      checkAllocate(retVal, "retVal in DeltaTest::get_ranks");
      std::copy(ranks_, ranks_+nInputs_+1, retVal);
   }
   return retVal;
}

