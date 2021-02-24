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
// Function for the class GradStatAnalyzer  
// Functions for the class FASTAnalyzer
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GradStatAnalyzer.h"
#include "Sampling.h"
#include "PDFBase.h"
#include "PDFNormal.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
GradStatAnalyzer::GradStatAnalyzer() : Analyzer()
{
   setName("GRADSTAT");
   sampler_     = NULL;
   inputPDFs_   = NULL;
   inputMeans_  = NULL;
   inputStdevs_ = NULL;
   threshold_   = 0.0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GradStatAnalyzer::~GradStatAnalyzer()
{
   if (sampler_ != NULL) delete sampler_;
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double GradStatAnalyzer::analyze(aData &adata)
{
   int     nInputs, nOutputs, nSamples, outputID, nLevels, *levelSeps;
   int     whichOutput, totalNSamples, prevNSamples, ss;
   int     ss2, inputID, minIndex, nLHSamples, iOne=1, iZero=0;
   int     *LHSampleStates, nReps, randomize;
   double  curThresh, distance, minDistance, interpolatedOutput, xtemp;
   double  *LHSampleInputs, *LHSampleOutputs, sampleMean, sampleStdDev;
   double  *localData1, *localData2, totalVol, ytemp;
   double  *X, *Y, *inputLBounds, *inputUBounds;
   PDFBase **normalPDFs;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   inputLBounds = adata.iLowerB_;
   inputUBounds = adata.iUpperB_;
   X        = adata.sampleInputs_;
   Y        = adata.sampleOutputs_;
   outputID = adata.outputID_;
   nLevels  = adata.currRefineLevel_;
   levelSeps = adata.refineSeparators_;
   threshold_   = adata.analysisThreshold_;
   sampler_     = (Sampling *) adata.sampler_;
   inputPDFs_   = adata.inputPDFs_;
   inputMeans_  = adata.inputMeans_;
   inputStdevs_ = adata.inputStdevs_;
   if (inputPDFs_ != NULL)
   {
      printOutTS(PL_INFO,
           "GradStatAnalyzer INFO: non-uniform distributions\n");
      printOutTS(PL_INFO,
           "    detected. Analysis will use these distributions.\n");
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR,"GradStatAnalyzer ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR,"    nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR,"    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;
   if ((whichOutput + nInputs) > nOutputs)
   {
      printOutTS(PL_ERROR,"GradStatAnalyzer ERROR: not enough outputs.\n");
      printOutTS(PL_ERROR,
           "                        nOutpus should be nInputs+1.\n");
      return PSUADE_UNDEFINED;
   } 
   
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,
      "            Gradient-based Uncertainty Analysis\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,
      " This analysis is for models that produce derviative\n");
   printOutTS(PL_INFO,
      " information in addition to typical outputs. Thus, the\n");
   printOutTS(PL_INFO,
      " number of outputs is expected to be nInputs+1. The\n");
   printOutTS(PL_INFO,
      " outputs and derivatives are used to construct a respone\n");
   printOutTS(PL_INFO,
      " surface which will be probed to compute basic statistics\n");
   printOutTS(PL_INFO," based on the given distributions.\n");
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO,
      "GradStatAnalyzer: non-Gradient-based analysis\n");
   computeMeanVariance(nInputs,nOutputs,nSamples,Y,&sampleMean,
                       &sampleStdDev,whichOutput);
   printDashes(PL_INFO, 0);

   totalNSamples = levelSeps[nLevels];
   if (nLevels > 0)
   {
      prevNSamples  = levelSeps[nLevels-1];
      curThresh = 0.0;
      totalVol = 0.0;
      for (ss = prevNSamples; ss < totalNSamples; ss++)
      {
         minDistance = 1.0e30;
         /* choose the point with the least of the max delta y */
         for (ss2 = 0; ss2 < prevNSamples; ss2++)
         {
            distance = 0.0;
            minIndex = 0;
            for (inputID = 0; inputID < nInputs; inputID++)
            {
               xtemp = X[ss*nInputs+inputID]-X[ss2*nInputs+inputID];
               xtemp *= Y[ss2*nOutputs+whichOutput+inputID+1];
               if (PABS(xtemp) > distance) distance = PABS(xtemp);
            }
            if (distance < minDistance) 
            {
               minIndex = ss2;
               minDistance = distance;
            }
         }
         interpolatedOutput = Y[minIndex*nOutputs+whichOutput];
         for (inputID = 0; inputID < nInputs; inputID++)
         {
            xtemp = X[ss*nInputs+inputID]-X[minIndex*nInputs+inputID];
            interpolatedOutput += 
                xtemp*Y[minIndex*nOutputs+whichOutput+inputID+1];
         }
         ytemp = Y[ss*nOutputs+whichOutput] - interpolatedOutput;

         curThresh += PABS(ytemp);
         totalVol += PABS(Y[ss*nOutputs+whichOutput]);
      }
      //curThresh = PABS(curThresh);
      //curThresh /= sampleMean;
      curThresh /= (double) (totalNSamples - prevNSamples);
      totalVol /= (double) (totalNSamples - prevNSamples);
      curThresh /= totalVol;
      printOutTS(PL_INFO,  "GradStatAnalyzer: convergence check %e < %e ?\n",
             curThresh, threshold_);
   }
   else curThresh = 1.0;


   randomize = 1;
   nReps = 1;
   if (sampler_ != NULL) delete sampler_;
   sampler_ = SamplingCreateFromID(PSUADE_SAMP_LHS);
   sampler_->setInputBounds(nInputs, inputLBounds, inputUBounds);
   sampler_->setOutputParams(nOutputs);
   sampler_->setSamplingParams(nSamples, nReps, randomize);
   sampler_->initialize(0);
   LHSampleInputs = new double[nSamples*nInputs];
   LHSampleOutputs = new double[nSamples*nOutputs];
   LHSampleStates = new int[nSamples];
   checkAllocate(LHSampleStates, "LHSampleStates in GradStat::analyze");
   sampler_->getSamples(nSamples, nInputs, nOutputs, LHSampleInputs,
                        LHSampleOutputs, LHSampleStates);
   sampler_->loadSamples(nSamples, nInputs, nOutputs, LHSampleInputs,
                         LHSampleOutputs, LHSampleStates);
   delete [] LHSampleInputs;
   delete [] LHSampleOutputs;
   delete [] LHSampleStates;

   nLHSamples = nSamples;
   while (nLHSamples < 8000)
   {
      sampler_->refine(2, 1, 1.0, 0, NULL);
      nLHSamples = sampler_->getNumSamples(); 
   }
   if (nLHSamples == nSamples)
   {
      printOutTS(PL_INFO,  "GradStatAnalyzer: nSamples >= 8000 - stop.\n");
      curThresh = 0.0;
   }
   LHSampleInputs = new double[nLHSamples*nInputs];
   LHSampleOutputs = new double[nLHSamples*nOutputs];
   LHSampleStates = new int[nLHSamples];
   checkAllocate(LHSampleStates, "LHSampleStates(2) in GradStat::analyze");
   sampler_->getSamples(nLHSamples, nInputs, nOutputs, LHSampleInputs,
                        LHSampleOutputs, LHSampleStates);

   localData1 = new double[nLHSamples];
   localData2 = new double[nLHSamples];
   checkAllocate(localData2, "localData2 in GradStat::analyze");
   normalPDFs = new PDFBase*[nInputs];
   for (inputID = 0; inputID < nInputs; inputID++)
   {
      if (inputPDFs_ != NULL && inputPDFs_[inputID] == 1)
         normalPDFs[inputID] = (PDFBase *) new PDFNormal(inputMeans_[inputID],
                                               inputStdevs_[inputID]);
      else normalPDFs[inputID] = NULL;
   }
   for (inputID = 0; inputID < nInputs; inputID++)
   {
      if (normalPDFs[inputID] != NULL)
      {
         for (ss = 0; ss < nLHSamples; ss++)
            localData1[ss] = LHSampleInputs[ss*nInputs+inputID];
         normalPDFs[inputID]->invCDF(nLHSamples,localData1,localData2,
                          inputLBounds[inputID],inputUBounds[inputID]);
         for (ss = 0; ss < nLHSamples; ss++)
            LHSampleInputs[ss*nInputs+inputID] = localData2[ss];
      }
   }
   delete [] localData1;
   delete [] localData2;
   for (inputID = 0; inputID < nInputs; inputID++)
      if (normalPDFs[inputID] != NULL) delete normalPDFs[inputID];
   delete [] normalPDFs;

   for (ss = 0; ss < nLHSamples; ss++)
   {
      minDistance = 1.0e30;
      for (ss2 = 0; ss2 < totalNSamples; ss2++)
      {
         distance = 0.0;
         for (inputID = 0; inputID < nInputs; inputID++)
         {
            xtemp = LHSampleInputs[ss*nInputs+inputID] -
                    X[ss2*nInputs+inputID];
            xtemp *= Y[ss2*nOutputs+whichOutput+inputID+1];
            if (PABS(xtemp) > distance) distance = PABS(xtemp);
         }
         if (distance < minDistance) 
         {
            minIndex = ss2;
            minDistance = distance;
         }
      }
      interpolatedOutput = Y[minIndex*nOutputs+whichOutput];
      for (inputID = 0; inputID < nInputs; inputID++)
      {
         xtemp = LHSampleInputs[ss*nInputs+inputID] - 
                 X[minIndex*nInputs+inputID];
         interpolatedOutput += xtemp*Y[minIndex*nOutputs+whichOutput+inputID+1];
      }
      LHSampleOutputs[ss] = interpolatedOutput;
   }
   printOutTS(PL_INFO,  "GradStatAnalyzer: Gradient-based analysis\n");
   computeMeanVariance(nInputs,iOne,nLHSamples,LHSampleOutputs,&sampleMean,
                       &sampleStdDev,iZero);
   printAsterisks(PL_INFO, 0);
   delete [] LHSampleInputs;
   delete [] LHSampleOutputs;
   delete [] LHSampleStates;
   return curThresh;
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int GradStatAnalyzer::computeMeanVariance(int nInputs, int nOutputs, 
                        int nSamples, double *Y, double *aMean, 
                        double *stdDev, int outputID)
{
   int    ss;
   double mean, variance;

   mean = 0.0;
   for (ss = 0; ss < nSamples; ss++) mean += Y[nOutputs*ss+outputID];
   mean /= (double) nSamples;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      variance += ((Y[nOutputs*ss+outputID] - mean) *
                   (Y[nOutputs*ss+outputID] - mean));
   }
   variance /= (double) nSamples;
   (*aMean) = mean;
   (*stdDev) = sqrt(variance);
   printOutTS(PL_INFO,  "GradStatAnalyzer: mean     = %24.16e\n", mean);
   printOutTS(PL_INFO,  "GradStatAnalyzer: std dev  = %24.16e\n", (*stdDev));
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
GradStatAnalyzer& GradStatAnalyzer::operator=(const GradStatAnalyzer &)
{
   printOutTS(PL_ERROR,
        "GradStatAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

