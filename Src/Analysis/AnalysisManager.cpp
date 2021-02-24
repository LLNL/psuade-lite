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
// Functions for AnalysisManager 
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <assert.h>
#include "pData.h"
#include "sysdef.h"
#include "AnalysisManager.h"
#include "AnovaAnalyzer.h"
#include "GradStatAnalyzer.h"
#include "MOATAnalyzer.h"
#include "MainEffectAnalyzer.h"
#include "TwoParamAnalyzer.h"
#include "SobolAnalyzer.h"
#include "RSFuncApproxAnalyzer.h"
#include "MomentAnalyzer.h"
#include "CorrelationAnalyzer.h"
#include "IntegrationAnalyzer.h"
#include "FASTAnalyzer.h"
#include "FFAnalyzer.h"
#include "LSAnalyzer.h"
#include "PCAnalyzer.h"
#include "OneSigmaAnalyzer.h"
#include "FORMAnalyzer.h"
#include "RSMSobol1Analyzer.h"
#include "RSMSobol2Analyzer.h"
#include "RSMSobolTSIAnalyzer.h"
#include "BootstrapAnalyzer.h"
#include "RSMSobolGAnalyzer.h"
#include "OneSampleAnalyzer.h"
#include "TwoSampleAnalyzer.h"
#include "MCMCAnalyzer.h"
#include "DeltaAnalyzer.h"
#include "EtaAnalyzer.h"
#include "GowerAnalyzer.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------ 
AnalysisManager::AnalysisManager()
{
   numAnalyzers_ = 27;
   analyzers_ = new Analyzer*[numAnalyzers_];;
   for (int ii = 0; ii < numAnalyzers_; ii++) analyzers_[ii] = NULL;
   sampler_ = NULL;
   analysisSampleErrors_ = NULL;
   logXsformFlags_ = NULL;

#ifdef HAVE_PYTHON
   AnalysisDataList = PyList_New(0);
#endif
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------ 
AnalysisManager::~AnalysisManager()
{ 
   for (int ii = 0; ii < numAnalyzers_; ii++)
      if (analyzers_[ii] != NULL) delete analyzers_[ii];
   if (analyzers_ != NULL) delete [] analyzers_;
   if (analysisSampleErrors_ != NULL) delete [] analysisSampleErrors_;
   if (logXsformFlags_ != NULL) delete [] logXsformFlags_;

#ifdef HAVE_PYTHON
   Py_DECREF(AnalysisDataList);
#endif
}

// ************************************************************************
// prepare for a new analysis run 
// ------------------------------------------------------------------------ 
int AnalysisManager::clearAnalyzers()
{ 
   for (int ii = 0; ii < numAnalyzers_; ii++)
   {
      if (analyzers_[ii] != NULL)
      {
        delete analyzers_[ii];
        analyzers_[ii] = NULL;
      }
   }
   return 0;
}

// ************************************************************************
// set sampler
// ------------------------------------------------------------------------ 
int AnalysisManager::setSampler(Sampling *sampler)
{ 
   sampler_ = sampler;
   return 0;
}

// ************************************************************************
// set up various analysis
// ------------------------------------------------------------------------
int AnalysisManager::setup(PsuadeData *psuadeIO)
{
   int   anaMethod, faType;
   pData pPtr;

   if (psuadeIO == NULL)
   {
      printOutTS(PL_ERROR,
                 "AnalysisManager setup ERROR: missing psuadeIO.\n");
      exit(1);
   }
   psuadeIO->getParameter("ana_method", pPtr);
   anaMethod = pPtr.intData_;
   psuadeIO->getParameter("ana_rstype", pPtr);
   faType = pPtr.intData_;
   setup(anaMethod, faType);
   return 0;
}

// ************************************************************************
// set up analysis (different parameter set)
// ------------------------------------------------------------------------
int AnalysisManager::setup(int anaMethod, int faType)
{
   int  targc;
   char *targv[3], request[100];

   if ((anaMethod & PSUADE_ANA_MOMENT) != 0 && (analyzers_[0] == NULL))
      analyzers_[0] = new MomentAnalyzer();
   if ((anaMethod & PSUADE_ANA_GLSA) != 0 && (analyzers_[1] == NULL))
      analyzers_[1] = new GradStatAnalyzer();
   if ((anaMethod & PSUADE_ANA_CORRELATION) != 0 && (analyzers_[2] == NULL))
      analyzers_[2] = new CorrelationAnalyzer();
   if ((anaMethod & PSUADE_ANA_ME) != 0 && (analyzers_[3] == NULL))
      analyzers_[3] = new MainEffectAnalyzer();
   if ((anaMethod & PSUADE_ANA_IE) != 0 && (analyzers_[4] == NULL))
      analyzers_[4] = new TwoParamAnalyzer();
   if ((anaMethod & PSUADE_ANA_MOAT) != 0 && (analyzers_[5] == NULL))
      analyzers_[5] = new MOATAnalyzer();
   if ((anaMethod & PSUADE_ANA_SOBOL) != 0 && (analyzers_[6] == NULL))
      analyzers_[6] = new SobolAnalyzer();
   if ((anaMethod & PSUADE_ANA_ANOVA) != 0 && (analyzers_[7] == NULL))
      analyzers_[7] = new AnovaAnalyzer();
   if ((anaMethod & PSUADE_ANA_RSFA) != 0 && (analyzers_[8] == NULL))
   {
      analyzers_[8] = new RSFuncApproxAnalyzer();
      targc = 2;
      strcpy(request, "rstype");
      targv[0] = (char *) request;
      targv[1] = (char *) &faType;
      analyzers_[8]->setParams(targc, targv);
   }
   if ((anaMethod & PSUADE_ANA_INTEGRATION) != 0 && (analyzers_[9] == NULL))
      analyzers_[9] = new IntegrationAnalyzer();
   if ((anaMethod & PSUADE_ANA_FAST) != 0 && (analyzers_[10] == NULL))
      analyzers_[10] = new FASTAnalyzer();
   if ((anaMethod & PSUADE_ANA_FF) != 0 && (analyzers_[11] == NULL))
      analyzers_[11] = new FFAnalyzer();
   if ((anaMethod & PSUADE_ANA_PCA) != 0 && (analyzers_[12] == NULL))
      analyzers_[12] = new PCAnalyzer();
   if ((anaMethod & PSUADE_ANA_ONESIGMA) != 0 && (analyzers_[13] == NULL))
   {
      analyzers_[13] = new OneSigmaAnalyzer();
      targc = 2;
      strcpy(request, "rstype");
      targv[0] = (char *) request;
      targv[1] = (char *) &faType;
      analyzers_[13]->setParams(targc, targv);
   }
   if ((anaMethod & PSUADE_ANA_FORM) != 0 && (analyzers_[14] == NULL))
   {
      analyzers_[14] = new FORMAnalyzer();
      targc = 2;
      strcpy(request, "rstype");
      targv[0] = (char *) request;
      targv[1] = (char *) &faType;
      analyzers_[14]->setParams(targc, targv);
   }
   if ((anaMethod & PSUADE_ANA_RSSOBOL1) != 0 && (analyzers_[15] == NULL))
      analyzers_[15] = new RSMSobol1Analyzer();
   if ((anaMethod & PSUADE_ANA_RSSOBOL2) != 0 && (analyzers_[16] == NULL))
      analyzers_[16] = new RSMSobol2Analyzer();
   if ((anaMethod & PSUADE_ANA_RSSOBOLTSI) != 0 && (analyzers_[17] == NULL))
      analyzers_[17] = new RSMSobolTSIAnalyzer();
   if ((anaMethod & PSUADE_ANA_BSTRAP) != 0 && (analyzers_[18] == NULL))
      analyzers_[18] = new BootstrapAnalyzer();
   if ((anaMethod & PSUADE_ANA_RSSOBOLG) != 0 && (analyzers_[19] == NULL))
      analyzers_[19] = new RSMSobolGAnalyzer();
   if ((anaMethod & PSUADE_ANA_1SAMPLE) != 0 && (analyzers_[20] == NULL))
      analyzers_[20] = new OneSampleAnalyzer();
   if ((anaMethod & PSUADE_ANA_2SAMPLE) != 0 && (analyzers_[21] == NULL))
      analyzers_[21] = new TwoSampleAnalyzer();
   if ((anaMethod & PSUADE_ANA_MCMC) != 0 && (analyzers_[22] == NULL))
      analyzers_[22] = new MCMCAnalyzer();
   if ((anaMethod & PSUADE_ANA_DTEST) != 0 && (analyzers_[23] == NULL))
      analyzers_[23] = new DeltaAnalyzer();
   if ((anaMethod & PSUADE_ANA_GOWER) != 0 && (analyzers_[24] == NULL))
      analyzers_[24] = new GowerAnalyzer();
   if ((anaMethod & PSUADE_ANA_ETEST) != 0 && (analyzers_[25] == NULL))
      analyzers_[25] = new EtaAnalyzer();
   if ((anaMethod & PSUADE_ANA_LSA) != 0 && (analyzers_[26] == NULL))
      analyzers_[26] = new LSAnalyzer();
   return 0;
}

// ************************************************************************
// set log transform flags (use logarithms of the inputs/outputs)
// ------------------------------------------------------------------------
int AnalysisManager::loadLogXsformFlags(int n, int *flags)
{
   if (n <= 0)
   {
      printOutTS(PL_ERROR,
           "AnalysisManager::loadLogXsformFlags ERROR: n <= 0.\n");
      exit(1);
   }
   if (logXsformFlags_ != NULL) delete [] logXsformFlags_;
   logXsformFlags_ = new int[n];
   for (int ii = 0; ii < n; ii++) logXsformFlags_[ii] = flags[ii]; 
   return 0;
}

// ************************************************************************
// analyze
// ------------------------------------------------------------------------
int AnalysisManager::analyze(PsuadeData *psuadeIO, int nLevels, 
                             int *levelSeps, int analysisOutputID)
{
  int    anaMethod, nInputs, nOutputs, nSamples, samplingMethod, nReps;
  int    refineFlag=1, jj, analysisTransform, wgtID, *xsforms, nActive;
  int    ii, outputLevel, *states, errCnt, *auxPDFs, pdfFlag, onlyMCMC;
  double *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs;
  double *auxMeans, *auxStds, dmin;
  double analysisThreshold, analysisData, *tempX, *tempY;
  char   **names, inStr[1000], lineIn[10000];
  pData  pPtr, pLowerB, pUpperB, pInpData, pOutData, pPDFs, pMeans, pStds;
  pData  pStates;
  aData  aPtr;

  nActive = 0;
  for (ii = 0; ii < numAnalyzers_; ii++) 
     if (analyzers_[ii] != NULL) nActive++;
  if (nActive == 0) return 0;
  onlyMCMC = 0;
  if (nActive == 1 && analyzers_[22] != NULL) onlyMCMC = 1;

  if (psuadeIO == NULL && onlyMCMC == 0)
  {
     printOutTS(PL_ERROR, "AnalysisManager ERROR: no DataIO.\n");
     return -1;
  }
  psuadeIO->getParameter("ana_method", pPtr);
  anaMethod = pPtr.intData_;
  psuadeIO->getParameter("input_ninputs", pPtr);
  nInputs = pPtr.intData_;
  psuadeIO->getParameter("output_noutputs", pPtr);
  nOutputs = pPtr.intData_;
  psuadeIO->getParameter("method_nsamples", pPtr);
  nSamples = pPtr.intData_;
  psuadeIO->getParameter("method_sampling", pPtr);
  samplingMethod = pPtr.intData_;
  psuadeIO->getParameter("method_nreplications", pPtr);
  nReps = pPtr.intData_;
  psuadeIO->getParameter("input_lbounds", pLowerB);
  iLowerB = pLowerB.dbleArray_;
  psuadeIO->getParameter("input_ubounds", pUpperB);
  iUpperB = pUpperB.dbleArray_;
  psuadeIO->getParameter("input_sample", pInpData);
  tempX = pInpData.dbleArray_;
  psuadeIO->getParameter("output_sample", pOutData);
  tempY = pOutData.dbleArray_;
  psuadeIO->getParameter("output_states", pStates);
  states = pStates.intArray_;
  psuadeIO->getParameter("ana_threshold", pPtr);
  analysisThreshold = pPtr.dbleData_;
  psuadeIO->getParameter("ana_transform", pPtr);
  analysisTransform = pPtr.intData_;
  if (tempX == NULL)
  {
    printOutTS(PL_WARN,
         "AnalysisManager WARNING: missing sample input data.\n");
    nSamples = 0;
  }
  if (tempY == NULL)
  {
    printOutTS(PL_WARN,
         "AnalysisManager WARNING: missing sample output data.\n");
    nSamples = 0;
  }
  if (analysisOutputID < 0)
  {
    psuadeIO->getParameter("ana_outputid", pPtr);
    analysisOutputID = pPtr.intData_;
    if (analysisOutputID < 0)
    {
      printOutTS(PL_ERROR,"AnalysisManager ERROR: output ID <= 0.\n");
      return -1;
    }
  }
  psuadeIO->getParameter("ana_regressionwgtid", pPtr);
  wgtID = pPtr.intData_;
  psuadeIO->getParameter("ana_diagnostics", pPtr);
  outputLevel = pPtr.intData_;

  if (logXsformFlags_ != NULL)
  {
    xsforms = new int[nInputs];
    aPtr.inputXsforms_ = xsforms;
    for (ii = 0; ii < nInputs; ii++)
      xsforms[ii] = logXsformFlags_[0] & 1;
    analysisTransform = logXsformFlags_[0] | logXsformFlags_[1];
  }
  else   
  {
    xsforms = new int[nInputs];
    aPtr.inputXsforms_ = xsforms;
    for (ii = 0; ii < nInputs; ii++) xsforms[ii] = analysisTransform & 1;
  }
  if (analysisTransform == 0)
  {
    printOutTS(PL_INFO,
         "No transformation (e.g. log) on sample inputs nor outputs.\n");
    sampleInputs = tempX;
    sampleOutputs = tempY;
  }
  else
  {
    if ((analysisTransform & 1) && tempX != NULL)
    {
      for (ii = 0; ii < nInputs*nSamples; ii++)
      {
        if (tempX[ii] <= 0.0)
        {
          analysisTransform &= 2;
          for (jj = 0; jj < nInputs; jj++) xsforms[jj] = 0;
          printOutTS(PL_INFO,
            "Some inputs are < 0 => TURN OFF INPUT TRANSFORMATION.\n");
          printf("Okay to continue with no input transformation (y or n)? ");
          scanf("%s", inStr);
          fgets(lineIn, 1000, stdin);
          if (inStr[0] != 'y') return -1;
          break;
        }
      }
    }
    if ((analysisTransform & 2) && tempY != NULL)
    {
      for (ii = 0; ii < nSamples*nOutputs; ii++)
      {
        if (tempY[ii] <= 0.0)
        {
          analysisTransform &= 1;
          printOutTS(PL_ERROR, 
            "SOME OUTPUT ARE < 0 => TURN OFF OUTPUT TRANSFORMATION.\n");
          printf("E.G. Sample %d output %d is negative.\n",
                 ii/nOutputs+1, ii%nOutputs+1);
          printf("Okay to continue with output no transformation (y or n)? ");
          scanf("%s", inStr);
          fgets(lineIn, 1000, stdin);
          if (inStr[0] != 'y') return -1;
          break;
        }
      }
    }
    if (analysisTransform == 0)
    {
      sampleInputs = tempX;
      sampleOutputs = tempY;
    }
    else
    {
      sampleInputs = new double[nSamples*nInputs];
      sampleOutputs = new double[nSamples*nOutputs];
      if (analysisTransform & 1) 
        printOutTS(PL_INFO, "Log transformation on inputs.\n");
      for (ii = 0; ii < nInputs; ii++)
      {
        if (xsforms[ii] == 1)
        {
          for (jj = 0; jj < nSamples; jj++)
            sampleInputs[jj*nInputs+ii] = log(tempX[jj*nInputs+ii]);
          if (iLowerB[ii] > 0) iLowerB[ii] = log(iLowerB[ii]);
          else
          {
            dmin = PSUADE_UNDEFINED;
            for (jj = 0; jj < nSamples; jj++)
              if (tempX[jj*nInputs+ii] < dmin) dmin = tempX[jj*nInputs+ii];
            iLowerB[ii] = log(dmin);
            iUpperB[ii] = log(iUpperB[ii]);
          }
        }
        else
        {
          for (jj = 0; jj < nSamples; jj++)
            sampleInputs[jj*nInputs+ii] = tempX[jj*nInputs+ii];
        }
      }
      if (analysisTransform & 2) 
        printOutTS(PL_INFO, "Log transformation on outputs.\n");
      if ((analysisTransform & 2) && tempY != NULL)
        for (ii = 0; ii < nSamples*nOutputs; ii++)
          sampleOutputs[ii] = log(tempY[ii]);
      else if (tempY != NULL)
        for (ii = 0; ii < nSamples*nOutputs; ii++)
          sampleOutputs[ii] = tempY[ii];
    }
  }
  errCnt = 0;
  if (sampleOutputs != NULL)
  {
    for (ii = 0; ii < nSamples*nOutputs; ii++)
      if (sampleOutputs[ii] == PSUADE_UNDEFINED) errCnt++;
  }
  if (errCnt > 0)
  {
    printOutTS(PL_WARN, 
         "AnalysisManager ERROR: undefined data found (%d).\n",
         errCnt);
    return -1;
  }

  aPtr.ioPtr_ = psuadeIO;
  aPtr.nSamples_ = nSamples;
  aPtr.nInputs_ = nInputs;
  aPtr.nOutputs_ = nOutputs;
  aPtr.outputID_ = analysisOutputID;
  aPtr.iLowerB_ = iLowerB;
  aPtr.iUpperB_ = iUpperB;
  aPtr.sampleInputs_ = sampleInputs;
  aPtr.sampleOutputs_ = sampleOutputs;
  aPtr.sampleStates_ = states;
  aPtr.printLevel_ = outputLevel;

  psuadeIO->getParameter("input_pdfs", pPDFs);
  pdfFlag = 0;
  for (ii = 0; ii < nInputs; ii++) pdfFlag += pPDFs.intArray_[ii];
  auxPDFs  = new int[nInputs];
  auxMeans = new double[nInputs];
  auxStds  = new double[nInputs];
  aPtr.inputPDFs_   = auxPDFs;
  aPtr.inputMeans_  = auxMeans;
  aPtr.inputStdevs_ = auxStds;
  if (pdfFlag > 0)
  {
    psuadeIO->getParameter("input_pdfs", pPDFs);
    psuadeIO->getParameter("input_means", pMeans);
    psuadeIO->getParameter("input_stdevs", pStds);
    for (ii = 0; ii < nInputs; ii++)
    {
      auxPDFs[ii]  = pPDFs.intArray_[ii];
      auxMeans[ii] = pMeans.dbleArray_[ii];
      auxStds[ii]  = pStds.dbleArray_[ii];
    }
  }
  else
  {
    for (ii = 0; ii < nInputs; ii++)
    {
      auxMeans[ii] = auxStds[ii] = 0.0;
      auxPDFs[ii] = 0;
    }
  } 

  if (analyzers_[0] != NULL)
  {
    aPtr.currRefineLevel_ = 0;
    aPtr.refineSeparators_ = NULL;
    aPtr.nSubSamples_ = nSamples/nReps;
    analysisData = analyzers_[0]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[1] != NULL)
  {
    aPtr.sampler_ = (void *) sampler_;
    aPtr.currRefineLevel_ = nLevels;
    aPtr.refineSeparators_ = levelSeps;
    aPtr.nSubSamples_ = nSamples/nReps;
    aPtr.analysisThreshold_ = analysisThreshold;
    analysisData = analyzers_[1]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager : analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[2] != NULL)
  {
    aPtr.currRefineLevel_ = nLevels;
    aPtr.refineSeparators_ = levelSeps;
    aPtr.nSubSamples_ = nSamples/nReps;
    analysisData = analyzers_[2]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[3] != NULL)
  {
    aPtr.currRefineLevel_ = 0;
    aPtr.refineSeparators_ = NULL;
    aPtr.nSubSamples_ = nSamples / nReps;
    analysisData = analyzers_[3]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[4] != NULL)
  {
    aPtr.currRefineLevel_ = 0;
    aPtr.refineSeparators_ = NULL;
    aPtr.nSubSamples_ = nSamples / nReps;
    analysisData = analyzers_[4]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
           analysisData, analysisThreshold);
  }

  if (analyzers_[5] != NULL)
  {
    aPtr.currRefineLevel_ = nLevels;
    aPtr.refineSeparators_ = levelSeps;
    aPtr.nSubSamples_ = 0;
    analysisData = analyzers_[5]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[6] != NULL)
  {
    aPtr.currRefineLevel_ = nLevels;
    aPtr.refineSeparators_ = levelSeps;
    aPtr.nSubSamples_ = 0;
    analysisData = analyzers_[6]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
           analysisData, analysisThreshold);
  }

  if (analyzers_[7] != NULL)
  {
    aPtr.currRefineLevel_ = nLevels;
    aPtr.refineSeparators_ = levelSeps;
    aPtr.nSubSamples_ = 0;
    analysisData = analyzers_[7]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
           analysisData, analysisThreshold);
  }

  if (analyzers_[8] != NULL)
  {
    aPtr.currRefineLevel_ = nLevels;
    aPtr.refineSeparators_ = levelSeps;
    aPtr.nSubSamples_ = 0;
    aPtr.regWgtID_ = wgtID;
    aPtr.cvFlag_ = 0;
    if (((refineFlag == 1) && (outputLevel >= 3)) || (outputLevel >= 4))
       aPtr.cvFlag_ = 1;
    if (outputLevel >= 5) aPtr.cvFlag_ = 2;
    analysisData = analyzers_[8]->analyze(aPtr);
    if (analysisSampleErrors_ != NULL) delete [] analysisSampleErrors_;
    analysisSampleErrors_ = aPtr.sampleErrors_;
    aPtr.sampleErrors_ = NULL;
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
           analysisData, analysisThreshold);
  }

  if (analyzers_[9] != NULL)
  {
    aPtr.currRefineLevel_ = nLevels;
    aPtr.refineSeparators_ = levelSeps;
    aPtr.nSubSamples_ = 0;
    analysisData = analyzers_[9]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[10] != NULL)
  {
    analysisData = analyzers_[10]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[11] != NULL)
  {
    aPtr.samplingMethod_ = samplingMethod;
    analysisData = analyzers_[11]->analyze(aPtr);
    aPtr.samplingMethod_ = -1;
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[12] != NULL)
  {
    analysisData = analyzers_[12]->analyze(aPtr);
    if (analysisData < analysisThreshold) refineFlag = 0;
    if (outputLevel > 0 && analysisData > 0.0)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
           analysisData, analysisThreshold);
    if (aPtr.nOutputs_ < nOutputs)
    {
      names = new char*[nOutputs];
      for (ii = 0; ii < aPtr.nOutputs_; ii++)
      {
        names[ii] = new char[20];
        sprintf(names[ii], "PC%d", ii+1);
      }
    }
    else names = NULL;
    psuadeIO->updateOutputSection(nSamples, aPtr.nOutputs_,
                             aPtr.sampleOutputs_, states, names);
    if (aPtr.nOutputs_ < nOutputs)
    {
      for (ii = 0; ii < aPtr.nOutputs_; ii++)
        delete [] names[ii];
      delete [] names;
    } 
    psuadeIO->writePsuadeFile(NULL,0);
  }

  if (analyzers_[13] != NULL)
  {
    analysisData = analyzers_[13]->analyze(aPtr);
    if (analysisData < analysisThreshold) refineFlag = 0;
    if (outputLevel > 0 && analysisData > 0.0)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
    psuadeIO->updateOutputSection(nSamples, aPtr.nOutputs_,
                             aPtr.sampleOutputs_, states, NULL);
    if (analysisSampleErrors_ != NULL) delete [] analysisSampleErrors_;
    analysisSampleErrors_ = aPtr.sampleErrors_;
    aPtr.sampleErrors_ = NULL;
  }

  if (analyzers_[14] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[14]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[15] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[15]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[16] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[16]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[17] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[17]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[18] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[18]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[19] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[19]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[20] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[20]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[21] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[21]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
       printOutTS(PL_INFO, 
            "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
            analysisData, analysisThreshold);
  }

  if (analyzers_[22] != NULL)
  {
    aPtr.analysisThreshold_ = analysisThreshold; 
    analysisData = analyzers_[22]->analyze(aPtr);
    if (analysisData < analysisThreshold || 
        analysisData == PSUADE_UNDEFINED) refineFlag = 0;
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
           analysisData, analysisThreshold);
  }
  
  if (analyzers_[23] != NULL)
  {
    analysisData = analyzers_[23]->analyze(aPtr);
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis metric = %8.2e\n",
           analysisData);
  }

  if (analyzers_[24] != NULL)
  {
    analysisData = analyzers_[24]->analyze(aPtr);
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis metric = %8.2e\n",
           analysisData);
  }

  if (analyzers_[25] != NULL)
  {
    analysisData = analyzers_[25]->analyze(aPtr);
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis metric = %8.2e\n",
           analysisData);
  }

  if (analyzers_[26] != NULL)
  {
    analysisData = analyzers_[26]->analyze(aPtr);
    if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
      printOutTS(PL_INFO, 
           "AnalysisManager: analysis metric = %8.2e\n",
           analysisData);
  }
#ifdef HAVE_PYTHON
  for (ii = 0; ii < numAnalyzers_; ii++)
    if (analyzers_[ii] != NULL)
      PyList_Append(AnalysisDataList, 
                    analyzers_[ii]->AnalysisDataDict);
#endif

  if (analysisTransform != 0)
  {
    if (sampleInputs  != NULL) delete [] sampleInputs;
    if (sampleOutputs != NULL) delete [] sampleOutputs;
  }
  delete [] xsforms;
  aPtr.sampleInputs_ = NULL;
  aPtr.sampleOutputs_ = NULL;
  aPtr.sampleStates_ = NULL;
  aPtr.iLowerB_ = NULL;
  aPtr.iUpperB_ = NULL;
  aPtr.inputXsforms_ = NULL;
  aPtr.inputPDFs_ = NULL;
  aPtr.inputMeans_ = NULL;
  aPtr.inputStdevs_ = NULL;
  delete [] auxPDFs;
  delete [] auxMeans;
  delete [] auxStds;
  return refineFlag;
}

// ************************************************************************
// analyze (simple, uniform distribution only)
// ------------------------------------------------------------------------
int AnalysisManager::analyze(int anaMethod, int nSamples, psVector &vLower, 
                             psVector &vUpper, psVector &vInputs, 
                             psVector &vOutputs, int outputLevel)
{
   int    nInputs, ii, *states, refineFlag=1;
   double *iLowerB, *iUpperB, *sampleInputs, *sampleOutputs;
   double analysisThreshold, analysisData;
   aData  aPtr;

   nInputs = vLower.length();
   ii = vInputs.length() / nInputs;
   assert(ii == nSamples);
   assert(nInputs == vUpper.length());
   assert(vOutputs.length() == nSamples);
   sampleInputs = vInputs.getDVector();
   sampleOutputs = vOutputs.getDVector();
   iLowerB = vLower.getDVector(); 
   iUpperB = vUpper.getDVector(); 
   states = new int[nSamples];
   for (ii = 0; ii < nSamples; ii++) states[ii] = 1;
   analysisThreshold = 1.0;

   aPtr.nSamples_ = nSamples;
   aPtr.nInputs_ = nInputs;
   aPtr.nOutputs_ = 1;
   aPtr.outputID_ = 0;
   aPtr.iLowerB_ = iLowerB;
   aPtr.iUpperB_ = iUpperB;
   aPtr.sampleInputs_ = sampleInputs;
   aPtr.sampleOutputs_ = sampleOutputs;
   aPtr.sampleStates_ = states;
   aPtr.printLevel_ = outputLevel;
   aPtr.currRefineLevel_ = 0;
   aPtr.refineSeparators_ = NULL;
   aPtr.nSubSamples_ = nSamples;

   for (ii = 0; ii < numAnalyzers_; ii++) 
   {
      if (analyzers_[ii] != NULL) 
      {
         analysisData = analyzers_[0]->analyze(aPtr);
         if (analysisData < analysisThreshold || 
             analysisData == PSUADE_UNDEFINED) refineFlag = 0;
         if (outputLevel > 0 && analysisData != PSUADE_UNDEFINED)
            printOutTS(PL_INFO, 
                 "AnalysisManager: analysis error = %8.2e <? %8.2e\n",
                 analysisData, analysisThreshold);
      }
   }

   aPtr.sampleInputs_ = NULL;
   aPtr.sampleOutputs_ = NULL;
   aPtr.sampleStates_ = NULL;
   aPtr.iLowerB_ = NULL;
   aPtr.iUpperB_ = NULL;
   aPtr.inputXsforms_ = NULL;
   delete [] states;
   return refineFlag;
}

// ************************************************************************
// analyze (one and two sample test only)
// ------------------------------------------------------------------------
int AnalysisManager::analyze(int anaMethod)
{
   aData aPtr;

   if ((anaMethod & PSUADE_ANA_1SAMPLE) != 0 && (analyzers_[20] != NULL))
        analyzers_[20]->analyze(aPtr);
   else if ((anaMethod & PSUADE_ANA_2SAMPLE) != 0 && 
            (analyzers_[21] != NULL))
        analyzers_[21]->analyze(aPtr);
   else
   {
      printOutTS(PL_INFO, 
           "AnalysisManager ERROR: this analyze call is not valid.\n");
      return 1;
   }
   return 0;
}

// ************************************************************************
// get sample errors (errors for individual sample point)
// ------------------------------------------------------------------------ 
double *AnalysisManager::getSampleErrors()
{ 
   return analysisSampleErrors_;
}

// ************************************************************************
// special request (send specialized parameters to individual methods)
// ------------------------------------------------------------------------ 
int AnalysisManager::specialRequest(int anaMethod, int narg, char **argv)
{ 
   if (((anaMethod & PSUADE_ANA_MOMENT) != 0) && (analyzers_[0] != NULL))
      analyzers_[0]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_GLSA) != 0) && (analyzers_[1] != NULL))
      analyzers_[1]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_CORRELATION) != 0) && 
        (analyzers_[2] != NULL))
      analyzers_[2]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_ME) != 0) && (analyzers_[3] != NULL))
      analyzers_[3]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_IE) != 0) && (analyzers_[4] != NULL))
      analyzers_[4]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_MOAT) != 0) && (analyzers_[5] != NULL))
      analyzers_[5]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_SOBOL) != 0) && (analyzers_[6] != NULL))
      analyzers_[6]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_ANOVA) != 0) && (analyzers_[7] != NULL))
      analyzers_[7]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_RSFA) != 0) && (analyzers_[8] != NULL))
      analyzers_[8]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_INTEGRATION) != 0) && 
        (analyzers_[9] != NULL))
      analyzers_[9]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_FAST) != 0) && (analyzers_[10] != NULL))
      analyzers_[10]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_FF) != 0) && (analyzers_[11] != NULL))
      analyzers_[11]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_PCA) != 0) && (analyzers_[12] != NULL))
      analyzers_[12]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_ONESIGMA) != 0) && 
        (analyzers_[13] != NULL))
      analyzers_[13]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_FORM) != 0) && (analyzers_[14] != NULL))
      analyzers_[14]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_RSSOBOL1) != 0) && 
        (analyzers_[15] != NULL))
      analyzers_[15]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_RSSOBOL2) != 0) && 
        (analyzers_[16] != NULL))
      analyzers_[16]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_RSSOBOLTSI) != 0) && 
        (analyzers_[17] != NULL))
      analyzers_[17]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_BSTRAP) != 0) && (analyzers_[18] != NULL))
      analyzers_[18]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_RSSOBOLG) != 0) && 
        (analyzers_[19] != NULL))
      analyzers_[19]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_1SAMPLE) != 0) && (analyzers_[20] != NULL))
      analyzers_[20]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_2SAMPLE) != 0) && (analyzers_[21] != NULL))
      analyzers_[21]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_MCMC) != 0) && (analyzers_[22] != NULL))
      analyzers_[22]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_DTEST) != 0) && (analyzers_[23] != NULL))
      analyzers_[23]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_GOWER) != 0) && (analyzers_[24] != NULL))
      analyzers_[24]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_ETEST) != 0) && (analyzers_[25] != NULL))
      analyzers_[25]->setParams(narg, argv);
   if (((anaMethod & PSUADE_ANA_LSA) != 0) && (analyzers_[26] != NULL))
      analyzers_[26]->setParams(narg, argv);
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
AnalysisManager& AnalysisManager::operator=(const AnalysisManager &)
{
   printOutTS(PL_ERROR, 
        "AnalysisManager operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// return the MOAT analyzer (this function is for library call to MOAT)
// ------------------------------------------------------------------------
MOATAnalyzer *AnalysisManager::getMOATAnalyzer()
{
   return reinterpret_cast<MOATAnalyzer *>(analyzers_[5]);
}

