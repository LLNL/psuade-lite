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
// Functions for the class EtaAnalyzer (Eta test)
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONGSnow
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
#include "EtaAnalyzer.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
EtaAnalyzer::EtaAnalyzer() : Analyzer(), nInputs_(0), dOrder_(0), dRanks_(0),
                             inputViolations_(0), totalViolation_(0)
{
   setName("DELTATEST");
   mode_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
EtaAnalyzer::~EtaAnalyzer()
{
   if (dOrder_) delete [] dOrder_;
   if (dRanks_) delete [] dRanks_;
   if (inputViolations_) delete [] inputViolations_;
}

// ************************************************************************
// perform test
// ------------------------------------------------------------------------
double EtaAnalyzer::analyze(aData &adata)
{
   int    printLevel, nSamples, nInputs, nOutputs, outputID, ss, ss2, ii2;
   int    *inputBins, ii, jj, kk, ii3, info, *indices, ind, nNeighs=3, imax;
   int    ***minDistIndices, *inputViolations, count;
   double *X, *Y, dist, *iLowerB, *iUpperB, *rangesInv2, *dRanks, dmax;
   double ***minDistMap, mean, stdev, dtemp, *minNeighs, *means, *dOrder;
   double ddata;
   char   pString[500];
   FILE   *fp;

   printLevel = adata.printLevel_;
   nSamples   = adata.nSamples_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   nInputs    = adata.nInputs_;
   nInputs_   = nInputs;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   outputID   = adata.outputID_;
   iLowerB    = adata.iLowerB_;
   iUpperB    = adata.iUpperB_;
   if (adata.inputPDFs_ != NULL)
   {
      count = 0;
      for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
      if (count > 0)
      {
         printOutTS(PL_WARN, 
              "EtaTest INFO: some inputs have non-uniform PDFs, but\n");
         printOutTS(PL_WARN, 
              "              they are not relevant in this analysis.\n");
      }
   }

   if (nSamples <= 1)
   {
      printOutTS(PL_ERROR, "EtaTest INFO: not meaningful to do this");
      printOutTS(PL_ERROR, "                    test when nSamples < 2.\n");
      return PSUADE_UNDEFINED;
   }
   if (X == NULL || Y == NULL)
   {
      printOutTS(PL_ERROR, "EtaTest ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   }
   info = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y[nOutputs*ii+outputID] == PSUADE_UNDEFINED) info = 1;
   if (info == 1)
   {
      printOutTS(PL_ERROR, "EtaTest: Some outputs are undefined.\n");
      printOutTS(PL_ERROR, "             Prune the undefined's first.\n");
      return PSUADE_UNDEFINED;
   }
   rangesInv2 = new double[nInputs];
   checkAllocate(rangesInv2, "rangesInv2 in etaTest::analyze");

   for (ii = 0; ii < nInputs; ii++) 
   {
      if ((iUpperB[ii] - iLowerB[ii]) > 0)
         rangesInv2[ii] = 1.0 / (iUpperB[ii] - iLowerB[ii]) /
                                (iUpperB[ii] - iLowerB[ii]);
      else
      {
         printOutTS(PL_ERROR, "EtaTest ERROR: problem with input range.\n");
         delete [] rangesInv2;
         return PSUADE_UNDEFINED;
      }
   }
 
   if (dOrder_) delete [] dOrder_;
   if (dRanks_) delete [] dRanks_;
   if (inputViolations_) delete [] inputViolations_;
   dOrder_ = NULL;
   dRanks_ = NULL;
   inputViolations_ = NULL;

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, "EtaTest for variable selection\n");
   printDashes(PL_INFO, 0);
   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO, "EtaTest Option: to use k > 1 neighbors.\n");
      printOutTS(PL_INFO, 
           "The larger k is, the larger the distinguishing power is.\n");
      sprintf(pString, "k neighbors analysis. What is k (>= 1, <= 20)? ");
      nNeighs = getInt(1, 20, pString);
   }
   if (nNeighs > nSamples/2)
   {
      nNeighs = nSamples / 2; 
      printOutTS(PL_INFO, 
           "EtaTest INFO: number of neighbors reset to be %d.\n",nNeighs);
   }  
   printOutTS(PL_INFO,"EtaTest INFO: number of neighbors = %d.\n",nNeighs);
   printEquals(PL_INFO, 0);

   inputBins = new int[nInputs];
   minDistMap = new double**[nSamples];
   checkAllocate(minDistMap, "minDistMap in etaTest::analyze");
   for (ii = 0; ii < nSamples; ii++)
   {
      minDistMap[ii] = new double*[nInputs];
      for (jj = 0; jj < nInputs; jj++)
      {
         minDistMap[ii][jj] = new double[nNeighs];
         for (ss = 0; ss < nNeighs; ss++)
            minDistMap[ii][jj][ss] = 0.0;
      }
   }
   minDistIndices = new int**[nSamples];
   for (ii = 0; ii < nSamples; ii++)
   {
      minDistIndices[ii] = new int*[nInputs];
      for (jj = 0; jj < nInputs; jj++)
      {
         minDistIndices[ii][jj] = new int[nNeighs];
         for (ss = 0; ss < nNeighs; ss++)
            minDistIndices[ii][jj][ss] = -1;
      }
   }
   minNeighs = new double[nNeighs];
   dRanks = new double[nInputs];
   dRanks_ = new double[nInputs_];
   means = new double[nInputs];
   indices = new int[nNeighs];
   inputViolations = new int[nInputs];
   inputViolations_ = new int[nInputs_];
   dOrder = new double[nInputs];
   dOrder_ = new double[nInputs_];
   checkAllocate(dOrder_, "dOrder in etaTest::analyze");

   srand(time(NULL));  
   for (ii = 0; ii < nInputs; ii++)
   {
      inputBins[ii] = 1;
      inputViolations[ii] = 0;
   }
   for (ss = 0; ss < nSamples; ss++)
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         inputBins[ii] = 0;
         for (jj = 0; jj < nNeighs; jj++)
         {
            minNeighs[jj] = PSUADE_UNDEFINED;
            indices[jj] = -1;
         }
         for (ss2 = 0; ss2 < nSamples; ss2++)
         {
            if (ss != ss2)
            {
               dist = 0.0;
               for (jj = 0; jj < nInputs; jj++)
               {
                  if (inputBins[jj] == 1)
                  {
                     dtemp = X[ss*nInputs+jj] - X[ss2*nInputs+jj];
                     dist += dtemp * dtemp * rangesInv2[jj];
                  }
               }
               if (dist > 0.0)
               {
                  if (dist < minNeighs[0])
                  {
                     minNeighs[0] = dist; 
                     indices[0] = ss2; 
                     for (jj = 1; jj < nNeighs; jj++)
                     {
                        if (minNeighs[jj] > minNeighs[jj-1])
                        {
                           dmax = minNeighs[jj];
                           minNeighs[jj] = minNeighs[jj-1];
                           minNeighs[jj-1] = dmax;
                           imax = indices[jj];
                           indices[jj] = indices[jj-1];
                           indices[jj-1] = imax;
                        }
                     }
                  }
               }
            }
         }
         for (jj = 0; jj < nNeighs; jj++) if (indices[jj] == -1) break;
         if (jj != nNeighs)
         {
            printOutTS(PL_ERROR, "EtaTest ERROR: cannot find neighbor.\n");
            printOutTS(PL_ERROR, 
                 "    Sample %d: nNeigh = %d (%d)\n",ss+1,jj,nNeighs);
            exit(1);
         }
         for (jj = 0; jj < nNeighs; jj++)
         {
            ind = indices[jj];
            minDistIndices[ss][ii][jj] = ind;
            minDistMap[ss][ii][jj] = 
               PABS((Y[ind*nOutputs+outputID]-Y[ss*nOutputs+outputID]));
            ddata = PABS((X[ind*nInputs]-X[ss*nInputs]))/
                         (iUpperB[0]-iLowerB[0]);
            if (ii == 0) ddata = 0.0;
            ii3 = 0;
            for (kk = 1; kk < nInputs; kk++)
            {
               dtemp = PABS((X[ind*nInputs+kk]-X[ss*nInputs+kk]))/
                            (iUpperB[kk]-iLowerB[kk]);
               if (kk != ii && dtemp > ddata)
               {
                  ddata = dtemp;
                  ii3 = kk;
               }
            }
            dtemp = PABS((X[ind*nInputs+ii]-X[ss*nInputs+ii]))/
                         (iUpperB[ii]-iLowerB[ii]);
            ddata /= dtemp;
            if (ddata < 1) inputViolations[ii3]++;
            if (printLevel > 3)
            {
               printOutTS(PL_INFO, 
                    "EtaTest: sample %5d, input %3d (%3d), dist = %e\n",
                    ss+1, ii3+1, ii+1, ddata);
            }
         }
         inputBins[ii] = 1;
      }
   }

   for (ii = 0; ii < nInputs; ii++) dRanks[ii] = 0;
   for (jj = 0; jj < nNeighs; jj++)
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         mean = stdev = 0.0;
         for (ss = 0; ss < nSamples; ss++) mean += minDistMap[ss][ii][jj];
         mean /= (double) (nSamples);
         for (ss = 0; ss < nSamples; ss++) 
            stdev += pow(minDistMap[ss][ii][jj]-mean,2.0);
         stdev = sqrt(stdev / nSamples);
         means[ii] = mean;
      }
      for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
      sortDbleList2(nInputs, means, dOrder);
      for (ii = 0; ii < nInputs; ii++)
      {
         kk = (int) dOrder[ii];
         dRanks[kk] += ii;
      } 
   }
   dist = dRanks[0];
   for (ii = 1; ii < nInputs; ii++) 
      if (dRanks[ii] > dist) dist = dRanks[ii];
   if (dist > 0.0) 
      for (ii = 0; ii < nInputs; ii++) dRanks[ii] /= dist;
   for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
   sortDbleList2(nInputs, dRanks, dOrder);
   printOutTS(PL_INFO, 
      "Order of importance based on %d nearest neighbors single effect:\n",
      nNeighs);
   for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, 
           "(E1)Rank %4d: input %4d (score = %d)\n", ii+1, (int)
           dOrder[nInputs-ii-1]+1, (int) (dRanks[nInputs-ii-1]*100));
   printEquals(PL_INFO, 0);
   kk = 0;
   for (ii = 0; ii < nInputs; ii++)
   {
      printOutTS(PL_INFO, "(E1) input %4d violation score = %5.2f\n", 
           ii+1, 100.0*inputViolations[ii]/nNeighs/nSamples);
      kk += inputViolations[ii];
   }
   totalViolation_ = 100.0*kk/nSamples/nInputs/nNeighs;
   printOutTS(PL_INFO, "(E1) Violation score = %5.2f\n",totalViolation_) ;
   printAsterisks(PL_INFO, 0);

   //save results
   for (ii = 0; ii < nInputs_; ii++)
   {
      dOrder_[ii] = dOrder[ii];
      dRanks_[ii] = dRanks[ii];
      inputViolations_[ii] = inputViolations[ii];
   }
   if (psAnaExpertMode_ == 1)
   {
      printf("Perform second order analysis ? (y or n) ");
      scanf("%s", pString);
   }
   else pString[0] = 'n';

   if (pString[0] == 'y')
   {
      for (ii = 0; ii < nInputs; ii++) inputBins[ii] = 1;
      for (ss = 0; ss < nSamples; ss++)
      {
         for (ii = 0; ii < nInputs; ii++)
         {
            inputBins[ii] = 0;
            for (ii2 = ii+1; ii2 < nInputs; ii2++)
            {
               inputBins[ii2] = 0;
               for (jj = 0; jj < nNeighs; jj++)
               {
                  minNeighs[jj] = PSUADE_UNDEFINED;
                  indices[jj] = -1;
               }
               for (ss2 = 0; ss2 < nSamples; ss2++)
               {
                  if (ss != ss2)
                  {
                     dist = 0.0;
                     for (jj = 0; jj < nInputs; jj++)
                     {
                        if (inputBins[jj] == 1)
                        {
                           dtemp = X[ss*nInputs+jj] - X[ss2*nInputs+jj];
                           dist += dtemp * dtemp * rangesInv2[jj];
                        }
                     }
                     if (dist > 0.0)
                     {
                        if (dist < minNeighs[0])
                        {
                           minNeighs[0] = dist; 
                           indices[0] = ss2; 
                           for (jj = 1; jj < nNeighs; jj++)
                           {
                              if (minNeighs[jj] > minNeighs[jj-1])
                              {
                                 dmax = minNeighs[jj];
                                 minNeighs[jj] = minNeighs[jj-1];
                                 minNeighs[jj-1] = dmax;
                                 imax = indices[jj];
                                 indices[jj] = indices[jj-1];
                                 indices[jj-1] = imax;
                              }
                           }
                        }
                     }
                  }
               }
               for (jj = 0; jj < nNeighs; jj++) if (indices[jj] == -1) break;
               if (jj != nNeighs)
               {
                  printOutTS(PL_ERROR, "EtaTest ERROR: cannot find neighbor.\n");
                  printOutTS(PL_ERROR, 
                       "    Sample %d: nNeigh = %d (%d)\n",ss+1,jj,nNeighs);
                  exit(1);
               }
               for (jj = 0; jj < nNeighs; jj++)
               {
                  ind = indices[jj];
                  minDistMap[ss][ii][jj] += 
                     PABS((Y[ind*nOutputs+outputID]-Y[ss*nOutputs+outputID]));
                  minDistMap[ss][ii2][jj] += 
                     PABS((Y[ind*nOutputs+outputID]-Y[ss*nOutputs+outputID]));
               }
               inputBins[ii2] = 1;
            }
            inputBins[ii] = 1;
         }
      }

      for (ii = 0; ii < nInputs; ii++) dRanks[ii] = 0;
      for (jj = 0; jj < nNeighs; jj++)
      {
         for (ii = 0; ii < nInputs; ii++)
         {
            mean = stdev = 0.0;
            for (ss = 0; ss < nSamples; ss++) mean += minDistMap[ss][ii][jj];
            mean /= (double) (nSamples);
            for (ss = 0; ss < nSamples; ss++) 
               stdev += pow(minDistMap[ss][ii][jj]-mean,2.0);
            stdev = sqrt(stdev / nSamples);
            means[ii] = mean;
         }
         for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
         sortDbleList2(nInputs, means, dOrder);
         for (ii = 0; ii < nInputs; ii++)
         {
            kk = (int) dOrder[ii];
            dRanks[kk] += ii;
         } 
      }
      dist = dRanks[0];
      for (ii = 1; ii < nInputs; ii++) 
         if (dRanks[ii] > dist) dist = dRanks[ii];
      if (dist > 0.0) 
         for (ii = 0; ii < nInputs; ii++) dRanks[ii] /= dist;
      for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
      sortDbleList2(nInputs, dRanks, dOrder);
      printOutTS(PL_INFO, 
         "Order of importance based on %d nearest neighbors twin effect:\n",
         nNeighs);
      for (ii = 0; ii < nInputs; ii++)
         printOutTS(PL_INFO, "(E2)Rank %4d: input %4d (score = %d)\n", ii+1, 
              (int) dOrder[nInputs-ii-1]+1,(int) (dRanks[nInputs-ii-1]*100));
      printAsterisks(PL_INFO, 0);
   }

   if (psAnaExpertMode_ == 1)
   {
      printf("Perform third order analysis ? (y or n) ");
      scanf("%s", pString);
   }
   else pString[0] = 'n';
   if (pString[0] == 'y')
   {
      for (ii = 0; ii < nInputs; ii++) inputBins[ii] = 1;
      for (ss = 0; ss < nSamples; ss++)
      {
         for (ii = 0; ii < nInputs; ii++)
         {
            inputBins[ii] = 0;
            for (ii2 = ii+1; ii2 < nInputs; ii2++)
            {
               inputBins[ii2] = 0;
               for (ii3 = ii2+1; ii3 < nInputs; ii3++)
               {
                  inputBins[ii3] = 0;
                  for (jj = 0; jj < nNeighs; jj++)
                  {
                     minNeighs[jj] = PSUADE_UNDEFINED;
                     indices[jj] = -1;
                  }
                  for (ss2 = 0; ss2 < nSamples; ss2++)
                  {
                     if (ss != ss2)
                     {
                        dist = 0.0;
                        for (jj = 0; jj < nInputs; jj++)
                        {
                           if (inputBins[jj] == 1)
                           {
                              dtemp = X[ss*nInputs+jj]-X[ss2*nInputs+jj];
                              dist += dtemp * dtemp * rangesInv2[jj];
                           }
                        }
                        if (dist > 0.0)
                        {
                           if (dist < minNeighs[0])
                           {
                              minNeighs[0] = dist; 
                              indices[0] = ss2; 
                              for (jj = 1; jj < nNeighs; jj++)
                              {
                                 if (minNeighs[jj] > minNeighs[jj-1])
                                 {
                                    dmax = minNeighs[jj];
                                    minNeighs[jj] = minNeighs[jj-1];
                                    minNeighs[jj-1] = dmax;
                                    imax = indices[jj];
                                    indices[jj] = indices[jj-1];
                                    indices[jj-1] = imax;
                                 }
                              }
                           }
                        }
                     }
                  }
                  for (jj = 0; jj < nNeighs; jj++) 
                     if (indices[jj] == -1) break;
                  if (jj != nNeighs)
                  {
                     printOutTS(PL_ERROR, 
                          "EtaTest ERROR: cannot find neighbor.\n");
                     printOutTS(PL_ERROR, 
                          "   Sample %d: nNeigh = %d (%d)\n",ss+1,jj,nNeighs);
                     exit(1);
                  }
                  for (jj = 0; jj < nNeighs; jj++)
                  {
                     ind = indices[jj];
                     minDistMap[ss][ii][jj] += 
                     PABS((Y[ind*nOutputs+outputID]-Y[ss*nOutputs+outputID]));
                     minDistMap[ss][ii2][jj] += 
                     PABS((Y[ind*nOutputs+outputID]-Y[ss*nOutputs+outputID]));
                     minDistMap[ss][ii3][jj] += 
                     PABS((Y[ind*nOutputs+outputID]-Y[ss*nOutputs+outputID]));
                  }
                  inputBins[ii3] = 1;
               }
               inputBins[ii2] = 1;
            }
            inputBins[ii] = 1;
         }
      }

      for (ii = 0; ii < nInputs; ii++) dRanks[ii] = 0;
      for (jj = 0; jj < nNeighs; jj++)
      {
         for (ii = 0; ii < nInputs; ii++)
         {
            mean = stdev = 0.0;
            for (ss = 0; ss < nSamples; ss++) mean += minDistMap[ss][ii][jj];
            mean /= (double) (nSamples);
            for (ss = 0; ss < nSamples; ss++) 
               stdev += pow(minDistMap[ss][ii][jj]-mean,2.0);
            stdev = sqrt(stdev / nSamples);
            means[ii] = mean;
         }
         for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
         sortDbleList2(nInputs, means, dOrder);
         for (ii = 0; ii < nInputs; ii++)
         {
            kk = (int) dOrder[ii];
            dRanks[kk] += ii;
         } 
      }
      dist = dRanks[0];
      for (ii = 1; ii < nInputs; ii++) 
         if (dRanks[ii] > dist) dist = dRanks[ii];
      if (dist > 0.0) 
         for (ii = 0; ii < nInputs; ii++) dRanks[ii] /= dist;
      for (ii = 0; ii < nInputs; ii++) dOrder[ii] = 1.0 * ii;
      sortDbleList2(nInputs, dRanks, dOrder);
      printOutTS(PL_INFO, 
         "Order of importance based on %d nearest neighbors 3-way effect:\n",
         nNeighs);
      for (ii = 0; ii < nInputs; ii++)
         printOutTS(PL_INFO, "(E3)Rank %4d: input %4d (score = %d)\n", ii+1, 
            (int) dOrder[nInputs-ii-1]+1, (int) (dRanks[nInputs-ii-1]*100));
      printAsterisks(PL_INFO, 0);
   }

   if (psAnaExpertMode_ == 1)
   {
      fp = fopen("psEtaData.m", "w");
      if (fp != NULL)
      {
         fprintf(fp, "clf\n");
         fprintf(fp, "n = %d;\n", nInputs);
         fprintf(fp, "X = 1:n;\n");
         for (jj = 0; jj < nNeighs; jj++) 
         {
            fprintf(fp, "A%d = [\n", jj+1);
            for (ss = 0; ss < nSamples; ss++) 
            {
               for (ii = 0; ii < nInputs; ii++)
                  fprintf(fp, "%13.5e ", minDistMap[ss][ii][jj]);
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            if (jj % 5 == 0)
               fprintf(fp, "plot(X,A%d,'bx')\n",jj+1);
            else if (jj % 5 == 1)
               fprintf(fp, "plot(X,A%d,'rx')\n",jj+1);
            else if (jj % 5 == 2)
               fprintf(fp, "plot(X,A%d,'kx')\n",jj+1);
            else if (jj % 5 == 3)
               fprintf(fp, "plot(X,A%d,'gx')\n",jj+1);
            else if (jj % 5 == 4)
               fprintf(fp, "plot(X,A%d,'cx')\n",jj+1);
            fprintf(fp, "hold on;\n");
            fprintf(fp, "pause\n");
         }
         fclose(fp);
         printOutTS(PL_INFO, 
              "Use the psEtaData.m file to examine the distributions.\n");
      }
   }
   printAsterisks(PL_INFO, 0);
 
   delete [] inputBins;
   delete [] minNeighs;
   for (ii = 0; ii < nSamples; ii++)
   {
      for (jj = 0; jj < nInputs; jj++)
         delete [] minDistMap[ii][jj];
      delete [] minDistMap[ii];
   }
   delete [] minDistMap;
   for (ii = 0; ii < nSamples; ii++)
   {
      for (jj = 0; jj < nInputs; jj++)
         delete [] minDistIndices[ii][jj];
      delete [] minDistIndices[ii];
   }
   delete [] minDistIndices;
   delete [] indices;
   delete [] rangesInv2;
   delete [] means;
   delete [] dRanks;
   delete [] dOrder;
   delete [] inputViolations;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
EtaAnalyzer& EtaAnalyzer::operator=(const EtaAnalyzer &)
{
   printOutTS(PL_ERROR, "EtaTest operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int EtaAnalyzer::get_mode()
{
   return mode_;
}

int EtaAnalyzer::get_nInputs()
{
   return nInputs_;
}
double *EtaAnalyzer::get_dOrder()
{
   double* retVal = NULL;
   if (dOrder_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in etaTest::get_dOrder");
      std::copy(dOrder_, dOrder_+nInputs_+1, retVal);
   }
   return retVal;
}
double *EtaAnalyzer::get_dRanks()
{
   double* retVal = NULL;
   if (dRanks_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in etaTest::get_dRanks");
      std::copy(dRanks_, dRanks_+nInputs_+1, retVal);
   }
   return retVal;
}
double *EtaAnalyzer::get_inputViolations()
{
   double* retVal = NULL;
   if (inputViolations_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in etaTest::get_inputViolations");
      std::copy(inputViolations_, inputViolations_+nInputs_+1, retVal);
   }
   return retVal;
}
double EtaAnalyzer::get_totalViolation()
{
   return totalViolation_;
}

