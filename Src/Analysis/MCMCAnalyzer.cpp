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
// Functions for the class MCMCAnalyzer
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
using namespace std;

#include "MCMCAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "pData.h"
#include "PDFManager.h"
#include "PDFBase.h"
#include "PDFNormal.h"
#include "PDFLogNormal.h"
#include "PDFTriangle.h"
#include "PDFBeta.h"
#include "PDFWeibull.h"
#include "PDFGamma.h"
#include "Psuade.h"
#include "Sampling.h"
#include "RSConstraints.h"
#include "PrintingTS.h"
#include "TwoSampleAnalyzer.h"
#include "Matrix.h"
#include "psStrings.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))
#define PS_INTERP 0

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MCMCAnalyzer::MCMCAnalyzer(): Analyzer(), nInputs_(0), nOutputs_(0), 
                              means_(0), sigmas_(0), mostLikelyInput_(0), 
                              mostLikelyOutput_(0)
{
  setName("MCMC");
  mode_ = 0;   // RS-based mode (mode = 1 : simulator-based mode)
  bfmode_ = 1; // brute force mode (default = 1, brute force)
  scheme_ = 0; // 0 - Gibbs or brute force, 1 - MH (not available yet)
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MCMCAnalyzer::~MCMCAnalyzer()
{
  if (means_) delete[] means_;
  if (sigmas_) delete[] sigmas_;
  if (mostLikelyInput_) delete[] mostLikelyInput_;
  if (mostLikelyOutput_) delete[] mostLikelyOutput_;
}

// ************************************************************************
// perform MCMC analysis 
// ------------------------------------------------------------------------
double MCMCAnalyzer::analyze(aData &adata)
{
  int    ii, ii2, jj, kk, cnt, count, nbins, iOne=1, iZero=0, nPlots;
  int    maxPts=257, faType, status, genPosteriors=0, rsErrFlag=0;
  int    maxSamples,burnInSamples,modelFormFlag=0,numChains=3, freq=1;
  int    nChainGood=0, masterOption=0, masterCount=1, combFlag=0;
  double *discFuncConstantMeans=NULL, *discFuncConstantStds=NULL;
  double psrfThreshold=1.05, dstatus, ddata;
  char   lineIn[1001], charString[1001], *cString;
  FILE   *fp=NULL;
  pData      pPtr, qData, pOutputs;
  FuncApprox **faPtrs=NULL, **faPtrs1=NULL;
  PDFBase    **inputPDFs;
  FunctionInterface *funcIO=NULL;
  PsuadeData *dataPtr=NULL, *ioOutFile=NULL;
  RSConstraints *constrPtr;
  psVector  vecRSValues, vecSamOut1;
  psIVector vecRSIndices, vecPlotIndices;

  int    printLevel = adata.printLevel_;
  int    nInputs    = adata.nInputs_;
  int    nOutputs   = adata.nOutputs_;
  int    nSamples   = adata.nSamples_;
  double *sampIns    = adata.sampleInputs_;
  double *samOuts    = adata.sampleOutputs_;
  int    *samStats   = adata.sampleStates_;
  double *xLower     = adata.iLowerB_;
  double *xUpper     = adata.iUpperB_;
  int    *pdfTypes = adata.inputPDFs_;
  double *pdfMeans = adata.inputMeans_;
  double *pdfStdvs = adata.inputStdevs_;
  nInputs_ = nInputs;
  nOutputs_ = nOutputs;
  dataPtr = adata.ioPtr_;
  if (dataPtr != NULL) dataPtr->getParameter("input_names", qData);

  // error checking
  if (nInputs <= 0 || nOutputs <= 0)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: in nInputs/nOutputs.\n");
    printOutTS(PL_ERROR,"    nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
    return PSUADE_UNDEFINED;
  }
  if (nSamples <= 0)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: in nSamples.\n");
    printOutTS(PL_ERROR,"    nSamples = %d\n", nSamples);
    printOutTS(PL_ERROR,"NOTE: even in simulator mode, a small sample\n");
    printOutTS(PL_ERROR,"      should be provided initially.\n");
    return PSUADE_UNDEFINED;
  }
  if (mode_ == 1 && dataPtr == NULL)
  {
    printOutTS(PL_ERROR,
         "MCMC ERROR: Direct simulation mode is requested but no\n");
    printOutTS(PL_ERROR,
         "     information on the simulator is given (missing\n");
    printOutTS(PL_ERROR,
         "     PSUADE data file that points to the simulator).\n");
    return PSUADE_UNDEFINED;
  }
  if (nOutputs > 1 && mode_ == 1 && dataPtr != NULL)
  {
    printOutTS(PL_ERROR,
         "MCMC ERROR: Direct simulation mode does not support\n");
    printOutTS(PL_ERROR,"     nOutputs > 1 at present.\n");
    printOutTS(PL_ERROR,
         "     Only response surface mode supports nOutputs > 1.\n");
    printOutTS(PL_ERROR,
         "     Suggestion: for simple cases with no design parameters\n");
    printOutTS(PL_ERROR,
         "            and no discrepancy modeling, nOutputs > 1\n");
    printOutTS(PL_ERROR,
         "            can be reduced to one output by computing the\n");
    printOutTS(PL_ERROR,
         "            likelihood function yourself and use zero mean\n");
    printOutTS(PL_ERROR,
         "            std dev = 1 for your one output).\n");
    return PSUADE_UNDEFINED;
  }
  status = 0;
  for (ii = 0; ii < nSamples*nOutputs; ii++)
    if (samOuts[ii] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: Some outputs are undefined\n");
    printOutTS(PL_ERROR,"     (>= 1e35). Prune the undefined\n");
    printOutTS(PL_ERROR,"     sample points and re-run.\n");
    return PSUADE_UNDEFINED;
  }
  if (samStats != NULL)
  {
    status = 0;
    for (ii = 0; ii < nSamples; ii++) status += samStats[ii];
    if (status != nSamples)
    {
      printOutTS(PL_WARN,"MCMC WARNING: Some outputs may not be ready.\n");
      printOutTS(PL_WARN,"     (Since their states are not 1's. Will\n");
      printOutTS(PL_WARN,"      continue anyway because outputs are not\n");
      printOutTS(PL_WARN,"      undefined.\n");
    }
  }
 
  if (psConfig_ != NULL)
  {
    cString = psConfig_->getParameter("MCMC_brute_force");
    if (cString != NULL) bfmode_ = 1;
    else
    {
      cString = psConfig_->getParameter("MCMC_gibbs");
      if (cString != NULL) bfmode_ = 0;
    }
  }
  if (bfmode_ == 1) return analyze_bf(adata);

  if (pdfTypes != NULL)
  {
    for (ii = 0; ii < nInputs; ii++) 
    {
      if (pdfTypes[ii] == PSUADE_PDF_SAMPLE)
      {
        printOutTS(PL_INFO,"MCMC INFO: user distribution detected.\n");
        printOutTS(PL_INFO,"           Switch to brute force MCMC.\n");
        return analyze_bf(adata);
      }
      else if (pdfTypes[ii] == 1000+PSUADE_PDF_NORMAL)
      {
        printOutTS(PL_INFO,"MCMC INFO: parameter joint PDF detected.\n");
        printOutTS(PL_INFO,"           Switch to brute force MCMC.\n");
        return analyze_bf(adata);
      }
      else if (pdfTypes[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
        printOutTS(PL_INFO,"MCMC INFO: parameter joint PDF detected.\n");
        printOutTS(PL_INFO,"           Switch to brute force MCMC.\n");
        return analyze_bf(adata);
      }
    }
  }
  if (scheme_ == 1) return analyze_mh(adata);

  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO,"*                     MCMC Optimizer\n");
  printEquals(PL_INFO, 0);
  if (printLevel > 0)
  {
    printOutTS(PL_INFO,"TO GAIN ACCESS TO DIFFERENT OPTIONS: TURN ON\n\n");
    printOutTS(PL_INFO," * ana_expert to finetune MCMC parameters, \n");
    printOutTS(PL_INFO,
         "   (e.g. sample size for burn-in can be adjusted).\n");
    printOutTS(PL_INFO,
         " * rs_expert to customize response surface for MCMC,\n");
    printOutTS(PL_INFO," * printlevel 3 to display more diagnostics info.\n");
    printOutTS(PL_INFO," * printlevel 4 to display even more diagnostics.\n");
    printOutTS(PL_INFO," * printlevel >=5 reserved only for expert only.\n");
    printDashes(PL_INFO,0);
    printOutTS(PL_INFO,
         "FEATURES AVAILABLE IN THE CURRENT VERSION OF MCMC:\n");
    printOutTS(PL_INFO," * Support other than uniform prior distributions\n");
    printOutTS(PL_INFO,
         " * Support likelihood functions from multiple outputs\n");
    printOutTS(PL_INFO," * Option to include response surface errors for\n");
    printOutTS(PL_INFO,"   polynomial regressions, bootstrapped MARS, and\n");
    printOutTS(PL_INFO,"   Gaussian process (GP).\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode.\n");
    printOutTS(PL_INFO,
         " * Option to include model form errors in the form of\n");
    printOutTS(PL_INFO,"   discrepancy models.\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode.\n");
    printOutTS(PL_INFO,
         " * Option to set some inputs as design parameters\n");
    printOutTS(PL_INFO,
         "   - to be specified in the observation data spec file\n");
    printOutTS(PL_INFO,
         " * Option to disable some parameters (set to default)\n");
    printOutTS(PL_INFO,
         "   - in case these parameters are not to be calibrated\n");
    printOutTS(PL_INFO,
         "   - use rs_index_file in PSUADE's ANALYSIS section\n");
    printOutTS(PL_INFO,"   - not available with discrepancy modeling\n");
    printOutTS(PL_INFO," * Option to generate a posterior sample\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode\n");
    printOutTS(PL_INFO,
         " * MCMC can be terminated gracefully by creating a file\n");
    printOutTS(PL_INFO,
         "   named 'psuade_stop' in the same directory during the\n");
    printOutTS(PL_INFO,"   run (in case it takes too long).\n");
    printOutTS(PL_INFO,
         " * For multi-modal posteriors, a large number of chains\n");
    printOutTS(PL_INFO,
         "   may be needed. The number of chains can be adjusted\n");
    printOutTS(PL_INFO,"   in ana_expert mode.\n");
    printOutTS(PL_INFO,
         " * In GM mode, you have a few options to choose from:\n");
    printOutTS(PL_INFO,
         "  1. track proposal distribution at each MCMC iteration\n");
    printOutTS(PL_INFO,
         "     (MCMCDistTrack.m will be created at each iteration to\n");
    printOutTS(PL_INFO,
         "      give a snapshot of the current proposal distribution)\n");
    printOutTS(PL_INFO,
         "      NOTE: iteration will pause until further instructions\n");
    printOutTS(PL_INFO,
         "  2. track posterior distributions after each MCMC cycle\n");
    printOutTS(PL_INFO,
         "     (MCMCChainHistogram.m will be created after each cycle\n");
    printOutTS(PL_INFO,
         "      to give snapshot of the current parameter posterior \n");
    printOutTS(PL_INFO,"      distributions.\n");
    printOutTS(PL_INFO,
         "      NOTE: iteration will pause until further instructions\n");
    printEquals(PL_INFO, 0);
  }
  if (psGMMode_ == 1)
  {
    printf("Please choose from the following option (or 0 if none): \n");
    printf("  1. track proposal distribution at each MCMC iteration\n");
    printf("  2. track posterior distributions after each MCMC cycle\n");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == '1') masterOption = 1;
    if (charString[0] == '2') masterOption = 2;
  }

  cleanUp();

  // get experimental data information from the spec file
  int    dnSamples=0, dnInputs=0, *designParams;
  double *dSamInputs=NULL,*dSamMeans=NULL,*dSamStdevs=NULL;
  dstatus = readSpecFile(nInputs, nOutputs, &dnSamples, &dnInputs, 
                 &designParams, &dSamInputs, &dSamMeans, &dSamStdevs,
                 combFlag, printLevel);
  if (dstatus != 0.0) return PSUADE_UNDEFINED;

  if (psAnaExpertMode_ == 1 && mode_ == 0)
  {
    printOutTS(PL_INFO,
        "*** OPTION TO INCLUDE RESPONSE SURFACE UNCERTAINTIES:\n");
    printOutTS(PL_INFO,
        "\nTo incorporate the response surface errors into the\n");
    printOutTS(PL_INFO,
        "likelihood function, make sure that either Gaussian\n");
    printOutTS(PL_INFO,
        "process/Kriging, polynomial regression, or bootstrapped\n");
    printOutTS(PL_INFO,
        "MARS response surface is selected in the simulation\n");
    printOutTS(PL_INFO,
        "data file. Otherwise, no RS uncertainties will be\n");
    printOutTS(PL_INFO,"included.\n\n");
    printOutTS(PL_INFO,
        "NOTE: if you don't know what this is, just say no.\n");
    printf( "===> Include response surface uncertainties? (y or n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') rsErrFlag = 1;
    printEquals(PL_INFO, 0);
  }

  if (dataPtr != NULL)
  {
    vecRSIndices.setLength(nInputs);
    for (ii = 0; ii < nInputs; ii++) vecRSIndices[ii] = 0;
    vecRSValues.setLength(nInputs);
    status = readIndexFile(dataPtr,nInputs,designParams,
                 vecRSIndices.getIVector(), vecRSValues.getDVector());
    count = 0;
    for (ii = 0; ii < nInputs; ii++) count += vecRSIndices[ii];
    if (status == -1)
    {
      if (designParams != NULL) delete [] designParams;
      if (dSamInputs != NULL) delete [] dSamInputs;
      if (dSamMeans != NULL) delete [] dSamMeans;
      if (dSamStdevs != NULL) delete [] dSamStdevs;
      return PSUADE_UNDEFINED;
    }
    dataPtr->getParameter("ana_rstype", pPtr);
    faType = pPtr.intData_;
  }
  else
  {
    printOutTS(PL_INFO,
         "MCMC INFO: since ioPtr=NULL, assume MARS as response surface.\n");
    faType = PSUADE_RS_MARS;
  }
  printOutTS(PL_INFO,
       "MCMC INFO: CREATING RESPONSE SURFACES FOR ALL OUTPUTS.\n");
  faPtrs = new FuncApprox*[nOutputs];
  vecSamOut1.setLength(nSamples);
  for (ii = 0; ii < nOutputs; ii++)
  {
    faType = -1;
    printOutTS(PL_INFO,
         "MCMC INFO: CREATING RESPONSE SURFACE FOR OUTPUT %d.\n",ii+1);
    faPtrs[ii] = genFA(faType, nInputs, iZero, nSamples);
    faPtrs[ii]->setNPtsPerDim(16);
    faPtrs[ii]->setBounds(xLower, xUpper);
    faPtrs[ii]->setOutputLevel(0);
    for (kk = 0; kk < nSamples; kk++) 
      vecSamOut1[kk] = samOuts[kk*nOutputs+ii];

    status = faPtrs[ii]->initialize(sampIns, vecSamOut1.getDVector());
    if (status != 0)
    {
      printOutTS(PL_ERROR,
           "MCMC ERROR: Unable to create response surface.\n");
      printOutTS(PL_ERROR,"            Consult PSUADE developers.\n");
      for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
      delete [] faPtrs;
      if (designParams != NULL) delete [] designParams;
      if (dSamInputs != NULL) delete [] dSamInputs;
      if (dSamMeans != NULL) delete [] dSamMeans;
      if (dSamStdevs != NULL) delete [] dSamStdevs;
      return PSUADE_UNDEFINED;
    }
  }

  //    ==> burnInSamples, maxSamples, nbins, plotIndices, nPlots
  maxSamples = 10000;
  burnInSamples = maxSamples / 2;
  nbins = 20;
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"*** CURRENT SETTINGS OF MCMC PARAMETERS: \n\n");
  printOutTS(PL_INFO,"MCMC Burn-in sample size      (default) = %d\n", 
             burnInSamples);
  printOutTS(PL_INFO,"MCMC sample increment         (default) = %d\n", 
             maxSamples);
  printOutTS(PL_INFO,"MCMC no. of bins in histogram (default) = %d\n",nbins);
  printOutTS(PL_INFO,
     "NOTE: sample increment - sample size to run before convergence check\n");
  printOutTS(PL_INFO,
     "NOTE: histogram nBins  - define granularity of histogram bar graph\n");
  printOutTS(PL_INFO, 
     "Turn on ana_expert mode to change these default settings.\n\n");
  if (psAnaExpertMode_ == 1)
  {
    sprintf(charString,"Enter sample increment (100 - 50000): ");
    maxSamples = getInt(100, 200000, charString);
    burnInSamples = maxSamples / 2;
    sprintf(charString,"Enter the number of histogram bins (5 - 25) : ");
    nbins = getInt(5, 50, charString);
  }
  if (psAnaExpertMode_ == 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
         "*** OPTION TO CREATE POSTERIOR PLOTS FOR SELECTED INPUTS:\n\n");
    printOutTS(PL_INFO,
         "MCMC will create MATLAB files for the posterior distributions.\n");
    printOutTS(PL_INFO,
         "You can choose to generate posterior plots for all inputs, or \n");
    printOutTS(PL_INFO,
         "just a selected few (in case there are too many inputs).\n");
    printf("Select inputs for which posterior plots are to be generated.\n");
    sprintf(charString,"Enter input number (-1 for all, 0 to terminate) : ");
    kk = 1;
    vecPlotIndices.setLength(nInputs);
    nPlots = 0;
    while (kk != 0 || nPlots < 1)
    {
      kk = getInt(-1, nInputs, charString);
      if (kk == -1)
      {
        for (ii = 0; ii < nInputs; ii++)
        {
          if (vecRSIndices.length() == 0 || vecRSIndices[ii] >= 0)
            if (designParams == NULL || designParams[ii] == 0) 
              vecPlotIndices[nPlots++] = ii;
        }
        break;
      }
      if (kk != 0)
      {
        if ((vecRSIndices.length() > 0) && (vecRSIndices[kk-1] < 0)) 
          printOutTS(PL_ERROR,
               "Input %d has been fixed by the rs index file (no plot).\n",
               kk+1);
        else if (designParams != NULL && designParams[kk-1] == 1)
          printOutTS(PL_ERROR,
               "Input %d is a design parameter (no plot)\n",kk);
        else 
          vecPlotIndices[nPlots++] = kk - 1;
      }
      if (kk == 0 && nPlots == 0)
        printOutTS(PL_ERROR,
             "You need to set at least 1 input for plotting posteriors.\n");
    }
    if (nPlots > 1) sortIntList(nPlots, vecPlotIndices.getIVector());
  }
  else
  {
    vecPlotIndices.setLength(nInputs);
    nPlots = 0;
    for (ii = 0; ii < nInputs; ii++) 
      if (vecRSIndices.length() == 0 || vecRSIndices[ii] >= 0)
        if (designParams == NULL || designParams[ii] == 0) 
          vecPlotIndices[nPlots++] = ii;
  }
  printOutTS(PL_INFO, 
       "MCMC Plot summary: input number to be plotted are (%d):\n",nPlots);
  for (ii = 0; ii < nPlots; ii++)
    printOutTS(PL_INFO, "   Input %4d\n", vecPlotIndices[ii]+1);

  // option to add discrepancy function and a posterior sample
  // ==> modelFormFlag, genPosteriors
  if (psAnaExpertMode_ == 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,"*** OPTION TO ADD A DISCREPANCY FUNCTION:\n\n");
    printOutTS(PL_INFO,
         "To use this feature, first make sure that the observation\n");
    printOutTS(PL_INFO,
         "data file specified earlier has design parameters specified\n");
    printOutTS(PL_INFO,
         "since the discrepancy function is to be a function of these\n");
    printOutTS(PL_INFO,
         "design parameters (if not, a constant discrepancy function\n");
    printOutTS(PL_INFO,"is to be created).\n");
    printOutTS(PL_INFO,
         "NOTE: if you don't know what this is, just say NO.\n");
    printf("===> Add discrepancy function ? (y or n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') modelFormFlag = 1;
    if (modelFormFlag == 1 && dnInputs == 0)
    {
      printOutTS(PL_INFO,
        "NOTE: No design inputs ==> discrepancy will be a constant function\n");
    }
    if (modelFormFlag == 1 && dnSamples == 1)
    {
      printOutTS(PL_INFO,
        "NOTE: 1 experiment ==> discrepancy will be a constant function\n");
    }

    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
       "*** OPTION TO CREATE A SAMPLE FROM THE POSTERIOR DISTRIBUTIONS:\n\n");
    printOutTS(PL_INFO,
         "In addition to generating the posterior distributions, you can\n");
    printOutTS(PL_INFO,
         "also draw a sample from these posteriors. The posterior sample\n");
    printOutTS(PL_INFO,
         "can be used as prior sample for another simulator/emulator.\n");
    printOutTS(PL_INFO,"NOTE: if you don't know what this is, just say no.\n");
    printf("==> Create posterior samples for the input parameters? (y/n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') genPosteriors = 1;
  } 
  printEquals(PL_INFO, 0);

  // setup input PDF, if there is any
  if (printLevel > 2) 
    printOutTS(PL_INFO,"*** INFORMATION ON PARAMETER PRIOR DISTRIBUTIONS\n");
  inputPDFs = new PDFBase*[nInputs];
  for (ii = 0; ii < nInputs; ii++)
  {
    inputPDFs[ii] = NULL;
    if (( vecRSIndices.length() == 0 || vecRSIndices[ii] >= 0) &&
        (designParams == NULL || designParams[ii] == 0)) 
    {
      if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_NORMAL)
      {
        inputPDFs[ii] = (PDFBase *) new PDFNormal(pdfMeans[ii], 
                                                  pdfStdvs[ii]);
        if (printLevel > 2) 
           printOutTS(PL_INFO,
                "Parameter %3d has normal prior distribution (%e,%e)\n",
                ii+1, pdfMeans[ii], pdfStdvs[ii]);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_LOGNORMAL)
      {
        inputPDFs[ii] = (PDFBase *) new PDFLogNormal(pdfMeans[ii],
                                                     pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has lognormal prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_TRIANGLE)
      {
        inputPDFs[ii] = (PDFBase *) new PDFTriangle(pdfMeans[ii],
                                                    pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has triangle prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_BETA)
      {
        inputPDFs[ii] = (PDFBase *) new PDFBeta(pdfMeans[ii],
                                                pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has beta prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_WEIBULL)
      {
        inputPDFs[ii] = (PDFBase *) new PDFWeibull(pdfMeans[ii],
                                                   pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has Weibull prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_GAMMA)
      {
        inputPDFs[ii] = (PDFBase *) new PDFGamma(pdfMeans[ii],
                                                 pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
             "Parameter %3d has gamma prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == 1000+PSUADE_PDF_NORMAL)
      {
        inputPDFs[ii] = NULL;
        printOutTS(PL_INFO,
           "Parameter %3d: multi-parameter normal distribution.\n",ii+1);
        printOutTS(PL_INFO,"               curently not supported.\n");
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        delete [] inputPDFs;
        return -1.0;
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
        inputPDFs[ii] = NULL;
        printOutTS(PL_INFO,
           "Parameter %3d: multi-parameter lognormal distribution.\n",ii+1);
        printOutTS(PL_INFO,"               curently not supported.\n");
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        delete [] inputPDFs;
        return -1.0;
      }
      else if (pdfTypes == NULL || pdfTypes[ii] == PSUADE_PDF_UNIFORM)
      {
        inputPDFs[ii] = NULL;
        if (printLevel > 2)
          printOutTS(PL_INFO,
             "Parameter %3d has uniform prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_SAMPLE)
      {
        inputPDFs[ii] = NULL;
        printOutTS(PL_INFO,
           "Parameter %3d: user-provided distribution currently not\n",ii+1);
        printOutTS(PL_INFO,"               supported.\n");
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        delete [] inputPDFs;
        return -1.0;
      }
    }
  }
  if (printLevel > 2) printEquals(PL_INFO, 0);

  //    ==> funcIO, freq (how often to do simulation and emulation)
  maxPts = nbins * 5;
  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO,"*** SETTING PROPOSAL DISTRIBUTION RESOLUTION\n");
    printOutTS(PL_INFO,
         "Since MCMC uses many function evaluations to construct\n");
    printOutTS(PL_INFO,
         "the proposal distributions, you have the option to set\n");
    printOutTS(PL_INFO,
         "how many points are used to construct it in order to\n");
    printOutTS(PL_INFO,"keep the inference cost reasonable.\n");
    printf("Sample size to construct proposal distribution. Default is %d.\n",
           maxPts);
    sprintf(charString,"Enter new sample size (%d - %d): ",nbins*3,nbins*10);
    maxPts = getInt(nbins*3, nbins*10, charString);
    maxPts = maxPts / nbins * nbins;
    printOutTS(PL_INFO,"Proposal distribution sample size = %d.\n", maxPts);
  }
  if (mode_ == 1 && dataPtr != NULL) 
  {
    funcIO = createFunctionInterface(dataPtr);
    if (nSamples == 0)
      printOutTS(PL_INFO, "MCMC: DIRECT SIMULATION has been set up.\n");
    else
      printOutTS(PL_INFO,
         "MCMC: DIRECT SIMULATION PLUS RESPONSE SURFACE have been set up.\n");
    printOutTS(PL_INFO,
      "MCMC INFO: make sure simulation output is in the right form\n");
    printOutTS(PL_INFO,
      "     which should not have been translated (mean) nor scaled (std)\n");
    printOutTS(PL_INFO,
      "     unless you compute the error measure yourself, in which case\n");
    printOutTS(PL_INFO,
      "     you should have use nOutputs=1, mean=0, and std dev=1 in the\n");
    printOutTS(PL_INFO,"     spec file.\n");
    if (psAnaExpertMode_ == 1 && nSamples > 0)
    {
      printOutTS(PL_INFO,
        "Since MCMC uses many function evaluations, you have the option\n");
      printOutTS(PL_INFO,
        "to set how frequent the simulator is invoked (response surface\n");
      printOutTS(PL_INFO,
        "is used otherwise) in constructing the proposal distribution\n");
      printOutTS(PL_INFO,
        "at each MCMC iteration. A frequency of f means that each MCMC\n");
      printOutTS(PL_INFO,
        "step uses f simulator runs. These f simulator runs will then be\n");
      printOutTS(PL_INFO,
        "supplemented with evaluations from the given response surface.\n");
      printOutTS(PL_INFO,
        "The default is f=10 (if you do not know what this is, enter 10).\n");
      sprintf(charString,
        "Max. number of simulator runs per MCMC step (1 - %d, default=10)? ",
              maxPts);
      kk = getInt(1, maxPts, charString);
      freq = maxPts * nInputs / kk;
      if (freq * kk != maxPts) freq++;
    }
    else if (nSamples > 0) freq = maxPts / 10;
    else                   freq = 1;
    if (nSamples > 0)
      printOutTS(PL_INFO,
         "Frequency of invoking the simulator has been set to %d\n", freq);
  }
  if (funcIO == NULL && faPtrs == NULL)
  {
    printOutTS(PL_ERROR,
         "MCMC ERROR: missing simulator and sample data - cannot proceed.\n");
    return PSUADE_UNDEFINED;
  }

  double *discOutputs=NULL;
  if (modelFormFlag == 1)
  {
    int    *ExpSamStates, ind, ExpNSamples, dfaType, dnPerDim=16;
    double *dOneSample, expdata, simdata;
    double *ExpSamInputs, *tSamInputs, *settings;

    ExpNSamples   = dnSamples;
    ExpSamInputs  = new double[ExpNSamples*nInputs];
    discOutputs   = new double[ExpNSamples*nOutputs];
    ExpSamStates  = new int[ExpNSamples];
    discFuncConstantMeans = new double[nOutputs];
    discFuncConstantStds  = new double[nOutputs];
    for (ii2 = 0; ii2 < nOutputs; ii2++)
      discFuncConstantMeans[ii2] = discFuncConstantStds[ii2] = 
                                   PSUADE_UNDEFINED;

    printOutTS(PL_INFO,
         "*** SELECT RESPONSE SURFACE TYPE FOR DISCREPANCY FUNCTION:\n");
    dfaType = -1;
    while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
    {
      writeFAInfo(-1);
      sprintf(charString, "===> Enter your choice : ");
      dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
    }

    settings = new double[nInputs];
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      //if (inputPDFs == NULL || 
      //    (inputPDFs != NULL && inputPDFs[ii2] == NULL))
      //   settings[ii2] = 0.5*(xLower[ii2] + xUpper[ii2]);
      //else
      //   settings[ii2] = inputPDFs[ii2]->getMean();
      settings[ii2] = 0.5*(xLower[ii2] + xUpper[ii2]);
      if (vecRSIndices.length() > 0 && vecRSIndices[ii2] < 0)
        settings[ii2] = vecRSValues[ii2];
    }

    faPtrs1 = new FuncApprox*[nOutputs];
    dOneSample = new double[nInputs];
    tSamInputs = new double[ExpNSamples*nInputs];
    int        askFlag = 0, *states=NULL;
    double     *tLowers = new double[nInputs];
    double     *tUppers = new double[nInputs];
    PsuadeData *dataPtr = new PsuadeData(); 
    char       **iNames;
    for (ii = 0; ii < nOutputs; ii++)
    {
      for (kk = 0; kk < dnSamples; kk++)
      {
        cnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          if (designParams != NULL && designParams[ii2] == 1)
          {
            dOneSample[ii2] = dSamInputs[kk*dnInputs+cnt];
            cnt++;
          }
          else dOneSample[ii2] = settings[ii2];
        }

        simdata = 0.0;
        if (psAnaExpertMode_ == 1 && askFlag == 0)
        {
          printOutTS(PL_INFO,
               "To create discrepancy functions, the calibration\n");
          printOutTS(PL_INFO,
               "parameters need to be set to some nominal values.\n");
          printOutTS(PL_INFO,
               "You can choose the nominal values, or it will be\n");
          printOutTS(PL_INFO,
               "set to the input means or mid points of the ranges.\n");
          printf( "Set nomininal values yourself ? (y or n) ");
          scanf("%s", charString);
          fgets(lineIn,1000,stdin);
          if (charString[0] == 'y')
          {
            for (ii2 = 0; ii2 < nInputs; ii2++)
            {
              if ((vecRSIndices.length() == 0 || vecRSIndices[ii2] >= 0) &&
                  (designParams == NULL || designParams[ii2] == 0))
              {
                printOutTS(PL_INFO,
                     "Input %d has lower and upper bounds = %e %e\n",
                     ii2+1, xLower[ii2], xUpper[ii2]);
                sprintf(charString,"Nominal value for input %d : ",
                        ii2+1);
                dOneSample[ii2] = getDouble(charString);
                settings[ii2]   = dOneSample[ii2];
              }
            }
          }
          else
          {
            for (ii2 = 0; ii2 < nInputs; ii2++)
            {
              if ((vecRSIndices.length() == 0 || vecRSIndices[ii2] >= 0) &&
                  (designParams == NULL || designParams[ii2] == 0))
              {
                printOutTS(PL_INFO,
                     "Nominal value for input %d = %e\n",
                     ii2+1, dOneSample[ii2]);
              }
            }
          }
          askFlag = 1;
        }

        if (funcIO != NULL)
          funcIO->evaluate(kk+1,nInputs,dOneSample,1,&simdata,0);
        else simdata = faPtrs[ii]->evaluatePoint(dOneSample);

        expdata = dSamMeans[kk*nOutputs+ii];

        if (printLevel >= 4)
        {
          printOutTS(PL_INFO,
               "Experiment %4d (out of %d) : ",kk+1,dnSamples);
          for (ii2 = 0; ii2 < nInputs; ii2++)
            printOutTS(PL_INFO,
                 "Input %7d = %12.4e ",ii2+1,dOneSample[ii2]);
          printOutTS(PL_INFO,
               "simuation, experimental data = %12.4e %12.4e\n",
               simdata, expdata);
        }

        discOutputs[ii*dnSamples+kk] = expdata - simdata;
      }

      if (dnInputs > 0) 
      {
        for (kk = 0; kk < ExpNSamples*dnInputs; kk++)
          tSamInputs[kk] = dSamInputs[kk];
        iNames = new char*[dnInputs];
        cnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          if (designParams[ii2] == 1)
          {
            iNames[cnt] = new char[100];
            tLowers[cnt] = xLower[ii2];
            tUppers[cnt] = xUpper[ii2];
            if (qData.strArray_ == NULL)
                 sprintf(iNames[cnt], "X%d", ii2+1);
            else strcpy(iNames[cnt], qData.strArray_[ii2]);
            cnt++;
          }
        }
        dataPtr->updateInputSection(ExpNSamples,dnInputs,NULL,tLowers,
                       tUppers,dSamInputs, iNames, NULL,NULL,NULL,NULL);
        for (ii2 = 0; ii2 < dnInputs; ii2++) delete [] iNames[ii2];
        delete [] iNames;
      }
      else
      {
        iNames = new char*[1];
        iNames[0] = new char[100];
        sprintf(iNames[0], "X0");
        for (ii2 = 0; ii2 < ExpNSamples; ii2++) tSamInputs[ii2] = 0.5;
        tLowers[0] = 0.0;
        tUppers[0] = 1.0;
        dataPtr->updateInputSection(ExpNSamples,iOne,NULL,tLowers,tUppers,
                             tSamInputs, iNames, NULL,NULL,NULL,NULL);
        delete [] iNames[0];
        delete [] iNames;
      }

      states = new int[ExpNSamples];
      for (kk = 0; kk < ExpNSamples; kk++) states[kk] = 1;
      iNames = new char*[1];
      iNames[0] = new char[100];
      sprintf(iNames[0], "Y%d", ii+1);
      dataPtr->updateOutputSection(ExpNSamples,iOne,
                     &discOutputs[dnSamples*ii],states,iNames);
      delete [] states;
      delete [] iNames[0];
      delete [] iNames;
      dataPtr->updateMethodSection(PSUADE_SAMP_MC, ExpNSamples, 1, -1, -1);
      sprintf(charString, "psDiscrepancyModel%d", ii+1);
      dataPtr->writePsuadeFile(charString, 0);

      printOutTS(PL_INFO,
          "Creating discrepancy response surface for output %d\n",ii+1);
      faPtrs1[ii] = NULL;
      if (dnInputs > 0 && dnSamples > 1)
      {
        faPtrs1[ii] = genFA(dfaType,dnInputs,iOne,ExpNSamples);
        if (faPtrs1[ii] == NULL)
        {
          printOutTS(PL_ERROR,
               "MCMC ERROR: cannot create discrepancy func for output %d.\n",
               ii+1);
          return -1.0;
        }
      }
      if (faPtrs1[ii] != NULL)
      {
        faPtrs1[ii]->setNPtsPerDim(dnPerDim);
        faPtrs1[ii]->setBounds(xLower, xUpper);
        faPtrs1[ii]->setOutputLevel(0);
      }
      else
      {
        discFuncConstantMeans[ii] = 0.0;
        for (kk = 0; kk < ExpNSamples; kk++)
          discFuncConstantMeans[ii] += discOutputs[ii*dnSamples+kk];
        discFuncConstantMeans[ii] /= (double) ExpNSamples;
        discFuncConstantStds[ii] = 0.0;
        for (kk = 0; kk < ExpNSamples; kk++)
          discFuncConstantStds[ii] += pow(discOutputs[ii*dnSamples+kk]-
                                          discFuncConstantMeans[ii],2.0);
        discFuncConstantStds[ii] = 
                sqrt(discFuncConstantStds[ii]/ExpNSamples);
      }
    }
    delete dataPtr;
    delete [] ExpSamInputs;
    delete [] dOneSample;
    delete [] settings;
    delete [] tSamInputs;
    delete [] tLowers;
    delete [] tUppers;
  }

  //    set up constraint filters, if any
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"MCMC INFO: creating constraints, if there is any.\n");
  printOutTS(PL_INFO,
       "     Constraints remove infeasible regions from the priors.\n");
  printOutTS(PL_INFO,
       "     Constraints can be specified by RS constraint files.\n");
  constrPtr = new RSConstraints();
  constrPtr->genConstraints(dataPtr);
  printEquals(PL_INFO, 0);

  // get information on how many chains to use and threshold for 
  // convergence check
  if (psAnaExpertMode_ == 1)
  {
    sprintf(charString, "How many MCMC chains? (2-20, default=3) : ");
    numChains = getInt(2,20,charString);
    sprintf(charString, "PSRF threshold? (1.0 - 1.2, default = 1.05) : ");
    psrfThreshold = getDouble(charString);
    if (psrfThreshold < 1.0 || psrfThreshold > 1.2)
    {
      printOutTS(PL_INFO,
           "MCMC : invalid PSRF threshold ==> reset to 1.05.\n");
      psrfThreshold = 1.05;
    }
  }

  //    set up for MCMC iterations
  int    *Ivec, **bins, ****bins2, globalIts, countTrack, dcnt;
  int    mcmcFail=0, sumBins, index2, nFail;
  int    ii3, jj2, kk2, index, length, iChain, chainCnt, mcmcIts;
  int    maxGlobalIts=20, chainCntSave, *chainStatus;
  double *XRange=NULL, *XGuess=NULL, *XDist=NULL, *XDesignS, *YDesignS;
  double *YDesignStds=NULL, *XGuessS=NULL, *YGuessS=NULL, *YGuessStds=NULL;
  double Xtemp, Ytemp, Ytemp2, *Xmax, Ymax, *s2Vec, *SDist;
  double ***XChains=NULL, stdev, stdev2, ddata2, WStat, BStat;
  double *chainMeans=NULL, *chainStdevs=NULL, *psrfs=NULL;
  TwoSampleAnalyzer *s2Analyzer=NULL;
  XRange  = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) XRange[ii] = xUpper[ii] - xLower[ii]; 
  XDist   = new double[maxPts+1];
  SDist   = new double[maxPts+1];
  XGuess  = new double[nInputs];
  XGuessS = new double[dnSamples*nInputs*(maxPts+1)];
  YGuessS = new double[dnSamples*nOutputs*(maxPts+1)];
  YGuessStds = new double[dnSamples*nOutputs*(maxPts+1)];
  XDesignS = new double[dnSamples*nInputs*(maxPts+1)];
  YDesignS = new double[dnSamples*nOutputs*(maxPts+1)];
  YDesignStds = new double[dnSamples*nOutputs*(maxPts+1)];
  means_ = new double[nInputs_];
  sigmas_ = new double[nInputs_];
  for (ii = 0; ii < nInputs; ii++) means_[ii] = sigmas_[ii] = 0.0;
  Xmax = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) Xmax[ii] = 0;
  mostLikelyInput_ = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) mostLikelyInput_[ii] = 0;
  mostLikelyOutput_ = new double[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++) mostLikelyOutput_[ii] = 0;
  Ymax = -PSUADE_UNDEFINED;
  Ivec = new int[nInputs];
  Ivec[nInputs-1] = -1;
  XChains = new double**[numChains];
  for (ii = 0; ii < numChains; ii++)
  {
    XChains[ii] = new double*[maxGlobalIts*maxSamples];
    for (jj = 0; jj < maxGlobalIts*maxSamples; jj++)
      XChains[ii][jj] = new double[nInputs+1];
  }
  chainMeans = new double[numChains];
  chainStdevs = new double[numChains];
  chainStatus  = new int[numChains];
  checkAllocate(chainStatus, "chainStatus in MCMC::analyze");
  for (ii = 0; ii < numChains; ii++) chainMeans[ii] = chainStdevs[ii] = 0.0;
  for (ii = 0; ii < numChains; ii++) chainStatus[ii] = 0;
  psrfs = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) psrfs[ii] = 0.0;
  bins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    bins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) bins[ii][jj] = 0;
  }
  bins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    bins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      bins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        bins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          bins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }
  s2Vec = new double[maxGlobalIts*maxSamples];
  checkAllocate(s2Vec, "s2Vec in MCMC::analyze");
  if (printLevel > 3) s2Analyzer = new TwoSampleAnalyzer();

  Sampling *sampler;
  if (nInputs > 50) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
  else              sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
  sampler->setInputBounds(nInputs, xLower, xUpper);
  sampler->setOutputParams(1);
  sampler->setSamplingParams(numChains, 1, 1);
  sampler->initialize(0);
  double *mcmcSeeds = new double[numChains*nInputs];
  double *tmpOuts = new double[numChains];
  int    *tmpStates = new int[numChains];
  checkAllocate(tmpStates, "tmpStates in MCMC::analyze");
  sampler->getSamples(numChains,nInputs,1,mcmcSeeds,tmpOuts,tmpStates);
  delete [] tmpOuts;
  delete [] tmpStates;
  delete sampler;
  for (iChain = 0; iChain < numChains; iChain++)
  {
    for (ii = 0; ii < nInputs; ii++)
    {
      ddata = mcmcSeeds[iChain*nInputs+ii];
      ddata = (ddata - xLower[ii]) / XRange[ii];
      mcmcSeeds[iChain*nInputs+ii] = ddata;
    }
  }
#if PS_INTERP == 1
  double b12, b11, b21, b22, det, aa, bb, cc, xd1, xd2;
  b12 = 1.0/maxPts; b11 = b12 * b12;
  b22 = 2.0/maxPts; b21 = b22 * b22;
  det = 1.0 / (b11 * b22 - b12 * b21);
#endif
#if PS_INTERP == 2
  FuncApprox *faDist;
  double *ZDist = new double[maxPts+1];
  for (jj = 0; jj <= maxPts; jj++) ZDist[jj] = 1.0 * jj / maxPts;
  faType = PSUADE_RS_MARS;
  faDist = genFA(faType, iOne, iZero, maxPts+1);
  faDist->setNPtsPerDim(16);
  double lo=0.0, hi=1.0;
  faDist->setBounds(&lo, &hi);
  faDist->setOutputLevel(-1);
  double xtrial, xtol=1.0e-4, xbeg, xend;
#endif
   
  // run the Gibbs algorithm
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO, "MCMC begins ... \n");
  fflush(stdout);
  fp = NULL;
  globalIts = chainCnt = 0;
  while (globalIts < maxGlobalIts)
  {
    for (iChain = 0; iChain < numChains; iChain++)
    {
      printOutTS(PL_INFO,"MCMC : Chain %d, iteration = %d\n",iChain+1,
                 globalIts+1);
      if (iChain == 0) chainCntSave = chainCnt;
      else             chainCnt     = chainCntSave;
      if (chainCnt == 0)
      {
        for (ii = 0; ii < nInputs; ii++)
        {
          if (designParams == NULL || designParams[ii] == 0)
               XGuess[ii] = mcmcSeeds[iChain*nInputs+ii];
          else XGuess[ii] = 0.5;
          XGuess[ii] = XGuess[ii]*(xUpper[ii]-xLower[ii])+xLower[ii];
          XChains[iChain][0][ii] = XGuess[ii];
        }
        if (vecRSIndices.length() > 0)
        {
          for (ii = 0; ii < nInputs; ii++)
            if (vecRSIndices[ii] < 0) XGuess[ii] = vecRSValues[ii];
        }
        XChains[iChain][0][nInputs] = -1;
      }
      else
      {
        for (ii = 0; ii < nInputs; ii++)
          XGuess[ii] = XChains[iChain][chainCnt-1][ii];
      }
      if (printLevel >= 0)
      {
        printOutTS(PL_INFO,"       Chain %d current initial guess : \n",
                   iChain+1);
        for (ii = 0; ii < nInputs; ii++)
          printOutTS(PL_INFO,"          Input %4d = %e\n",
                     ii+1,XGuess[ii]);
      }
       
      mcmcIts = countTrack = 0;
      while (mcmcIts < maxSamples)
      {
        count = (mcmcIts+1) / (maxSamples/10);
        if (count != countTrack)
        { 
          countTrack++;
          printOutTS(PL_INFO, "%3.0f%% ",10.0*countTrack );
          fflush(stdout);
        }
        jj = Ivec[nInputs-1];
        generateRandomIvector(nInputs, Ivec);
        if (Ivec[0] == jj && nInputs > 1)
        {
          Ivec[0] = Ivec[nInputs-1];
          Ivec[nInputs-1] = jj;
        }
        for (kk = 0; kk < nInputs; kk++)
        {
          ii = Ivec[kk];
          if ((vecRSIndices.length() == 0 ||
              (vecRSIndices.length() > 0 && vecRSIndices[ii] >=0)) && 
              (designParams == NULL || designParams[ii] == 0))
          {
            Xtemp = XGuess[ii];
            fp = NULL;
            if (masterOption == 1 && masterCount == 1)
            {
              fp = fopen("MCMCDistTrack.m", "w");
              if (fp == NULL)
              {
                printOutTS(PL_WARN,
                  "MCMC ERROR: cannot write to MCMCDistTrack.m file.\n");
              }
              else
              {
                fprintf(fp,
                  "%% This file contains individual terms in the\n");
                fprintf(fp,
                  "%% exponent of the likelihood function, i.e.\n");
                fprintf(fp,"%% in each row of first column : \n");
                fprintf(fp,"%% S = 1/(p*n) sum_{k=1}^p sum_{i=1)^n ");
                fprintf(fp,"(Y_ki - m_ki)^2/sd_ki^2\n");
                fprintf(fp,"%% next columns: log-likelihood terms\n");
                fprintf(fp,"nOuts = %d;\n", nOutputs);
                fprintf(fp,"nObs  = %d;\n", dnSamples);
                fprintf(fp,"nPts  = %d;\n", maxPts+1);
                fprintf(fp,"A = [\n");
              }
            }
            if (nOutputs == 1)
            {
              cnt = 0;
              for (jj = 0; jj <= maxPts; jj++)
              {
                XGuess[ii] = xLower[ii]+jj*XRange[ii]/maxPts;
                   
                cnt++;
                if (cnt >= freq) cnt = 0; 
                index = jj * dnSamples;
                for (kk2 = 0; kk2 < dnSamples; kk2++) 
                {
                  dcnt = 0;
                  for (ii2 = 0; ii2 < nInputs; ii2++) 
                  {
                    XGuessS[(index+kk2)*nInputs+ii2] = XGuess[ii2];
                    if (designParams != NULL && 
                        designParams[ii2] == 1)
                    {
                      XGuessS[(index+kk2)*nInputs+ii2] = 
                                      dSamInputs[kk2*dnInputs+dcnt];
                      XDesignS[(index+kk2)*dnInputs+dcnt] = 
                                      dSamInputs[kk2*dnInputs+dcnt];
                      dcnt++;
                    }
                  }
                }
                for (kk2 = 0; kk2 < dnSamples; kk2++) 
                {
                  YDesignS[index+kk2] = YDesignStds[index+kk2] = 0.0;
                  YGuessS[index+kk2] = YGuessStds[index+kk2] = 0.0;
                }
              }

              if (rsErrFlag == 1)
              {
                if (funcIO == NULL || cnt > 0)
                  faPtrs[0]->evaluatePointFuzzy((maxPts+1)*dnSamples,
                                           XGuessS,YGuessS,YGuessStds);
                else
                {
                  for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                  {
                    funcIO->evaluate(kk2+1,nInputs,
                                     &XGuessS[kk2*nInputs],1,
                                     &YGuessS[kk2],0);
                    YGuessStds[kk2] = 0.0;
                  }
                }
                if (faPtrs1 != NULL && faPtrs1[0] != NULL)
                {
                  for (ii3 = 0; ii3 <= maxPts; ii3++)
                  {
                    for (kk2 = 0; kk2 < dnSamples; kk2++)
                    {
                      index = ii3 * dnSamples + kk2;
                      YDesignS[index] = discOutputs[kk2];
                    } 
                  }
                }
                else if (discFuncConstantMeans != NULL &&
                         discFuncConstantMeans[0] != PSUADE_UNDEFINED)
                {
                  for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                    YDesignS[kk2] = discFuncConstantMeans[0];
                }
              }
              else
              {
                if (funcIO == NULL || cnt > 0)
                         faPtrs[0]->evaluatePoint((maxPts+1)*dnSamples,
                                                  XGuessS,YGuessS);
                else
                {
                  for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                    funcIO->evaluate(kk2+1,nInputs,&XGuessS[kk2*nInputs],
                                     1,&YGuessS[kk2],0);
                }
                if (faPtrs1 != NULL && faPtrs1[0] != NULL)
                {
                  for (ii3 = 0; ii3 <= maxPts; ii3++)
                  {
                    for (kk2 = 0; kk2 < dnSamples; kk2++)
                    {
                      index = ii3 * dnSamples + kk2;
                      YDesignS[index] = discOutputs[kk2];
                    } 
                  }
                }
                else if (discFuncConstantMeans != NULL &&
                         discFuncConstantMeans[0] != PSUADE_UNDEFINED)
                {
                  for (kk2 = 0; kk2 < (maxPts+1)*dnSamples; kk2++) 
                    YDesignS[kk2] = discFuncConstantMeans[0];
                }
              }

              for (jj = 0; jj <= maxPts; jj++)
              {
                index = jj * dnSamples;
                for (kk2 = 0; kk2 < dnSamples; kk2++) 
                  YGuessS[index+kk2] += YDesignS[index+kk2];

                XDist[jj] = 0.0;
                for (kk2 = 0; kk2 < dnSamples; kk2++)
                {
                  Ytemp = YGuessS[index+kk2];
                  stdev = YGuessStds[index+kk2];
                  stdev2 = YDesignStds[index+kk2];
                  Ytemp2 = pow((Ytemp-dSamMeans[kk2]),2.0) /
                                (pow(dSamStdevs[kk2],2.0)+
                                 stdev*stdev+stdev2*stdev2);
                  XDist[jj] += Ytemp2;
                }
                XDist[jj] = XDist[jj] / dnSamples;
                if (masterOption == 1 && fp != NULL) 
                {
                  fprintf(fp, "%e ", XDist[jj]);
                  for (kk2 = 0; kk2 < dnSamples; kk2++)
                  {
                    Ytemp = YGuessS[index+kk2];
                    stdev = YGuessStds[index+kk2];
                    stdev2 = YDesignStds[index+kk2];
                    Ytemp2 = pow((Ytemp-dSamMeans[kk2]),2.0) /
                               (pow(dSamStdevs[kk2],2.0)+
                               stdev*stdev+stdev2*stdev2);
                    fprintf(fp, "%e ", Ytemp2);
                  }
                  fprintf(fp, "\n");
                }
              }
            }
            else
            {
              for (jj = 0; jj <= maxPts; jj++)
              {
                XGuess[ii] = xLower[ii]+jj*XRange[ii]/maxPts;
                index = jj * dnSamples;
                for (kk2 = 0; kk2 < dnSamples; kk2++) 
                {
                  dcnt = 0;
                  for (ii2 = 0; ii2 < nInputs; ii2++) 
                  {
                    XGuessS[(index+kk2)*nInputs+ii2] = XGuess[ii2];
                    if (designParams != NULL && designParams[ii2] == 1)
                    {
                      XGuessS[(index+kk2)*nInputs+ii2] = 
                            dSamInputs[kk2*dnInputs+dcnt];
                      XDesignS[(index+kk2)*dnInputs+dcnt] = 
                            dSamInputs[kk2*dnInputs+dcnt];
                      dcnt++;
                    }
                  }
                }
                for (ii2 = 0; ii2 < dnSamples*nOutputs; ii2++) 
                {
                  YDesignS[index*nOutputs+ii2] = 
                                 YDesignStds[index*nOutputs+ii2] = 0.0;
                  YGuessS[index*nOutputs+ii2] = 
                                 YGuessStds[index*nOutputs+ii2] = 0.0;
                }
              }
              for (ii2 = 0; ii2 < nOutputs; ii2++) 
              {
                if (rsErrFlag == 1)
                {
                  faPtrs[ii2]->evaluatePointFuzzy((maxPts+1)*dnSamples,
                                 XGuessS,&YGuessS[ii2*dnSamples*(maxPts+1)],
                                 &YGuessStds[ii2*dnSamples*(maxPts+1)]);
                  if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                  {
                    for (ii3 = 0; ii3 <= maxPts; ii3++)
                    {
                      for (kk2 = 0; kk2 < dnSamples; kk2++)
                      {
                        index = ii2*dnSamples*(maxPts+1)+
                                       ii3*dnSamples+kk2;
                        YDesignS[index] = 
                                    discOutputs[ii2*dnSamples+kk2];
                      } 
                    }
                  }
                  else if (discFuncConstantMeans != NULL &&
                           discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
                  {
                    for (kk2 = 0; kk2 < dnSamples*(maxPts+1); kk2++)
                    {
                      YDesignS[ii2*dnSamples*(maxPts+1)+kk2] = 
                                           discFuncConstantMeans[ii2];
                      YDesignStds[ii2*dnSamples*(maxPts+1)+kk2]=0.0;
                    }
                  }
                }
                else
                {
                  faPtrs[ii2]->evaluatePoint((maxPts+1)*dnSamples,
                                XGuessS,&YGuessS[ii2*dnSamples*(maxPts+1)]);
                  if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
                  {
                    for (ii3 = 0; ii3 <= maxPts; ii3++)
                    {
                      for (kk2 = 0; kk2 < dnSamples; kk2++)
                      {
                        index = ii2*dnSamples*(maxPts+1)+ii3*dnSamples+kk2;
                        YDesignS[index] = 
                                        discOutputs[ii2*dnSamples+kk2];
                      } 
                    }
                  }
                  else if (discFuncConstantMeans != NULL &&
                           discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
                  {
                    for (kk2 = 0; kk2 < dnSamples*(maxPts+1); kk2++)
                      YDesignS[ii2*dnSamples*(maxPts+1)+kk2] = 
                                           discFuncConstantMeans[ii2];
                  }
                  for (kk2 = 0; kk2 < dnSamples*(maxPts+1); kk2++)
                    YGuessStds[ii2*dnSamples*(maxPts+1)+kk2] = 
                              YDesignStds[ii2*dnSamples*(maxPts+1)+kk2]=0;
                }
              }
              for (jj = 0; jj <= maxPts; jj++)
              {
                XDist[jj] = 0.0;
                index = jj * dnSamples;
                for (ii2 = 0; ii2 < nOutputs; ii2++) 
                {
                  for (kk2 = 0; kk2 < dnSamples; kk2++)
                  {
                    Ytemp=YGuessS[ii2*dnSamples*(maxPts+1)+index+kk2]+ 
                           YDesignS[ii2*dnSamples*(maxPts+1)+index+kk2];
                    stdev=YGuessStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                    stdev2=YDesignStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                    Ytemp2=pow((Ytemp-dSamMeans[kk2*nOutputs+ii2]),2.0)/
                                 (pow(dSamStdevs[kk2*nOutputs+ii2],2.0) + 
                                  stdev*stdev + stdev2*stdev2);
                    XDist[jj] += Ytemp2;
                  }
                }
                if (combFlag == 1)
                  XDist[jj] = XDist[jj] / (double) dnSamples;
                if (masterOption == 1 && fp != NULL) 
                {
                  fprintf(fp, "%e ", XDist[jj]);
                  for (ii2 = 0; ii2 < nOutputs; ii2++) 
                  {
                    for (kk2 = 0; kk2 < dnSamples; kk2++)
                    {
                      Ytemp=YGuessS[ii2*dnSamples*(maxPts+1)+index+kk2]+ 
                             YDesignS[ii2*dnSamples*(maxPts+1)+index+kk2];
                      stdev = YGuessStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                      stdev2 = YDesignStds[ii2*dnSamples*(maxPts+1)+index+kk2];
                      Ytemp2 = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2]),2.0)/
                               (pow(dSamStdevs[kk2*nOutputs+ii2],2.0)+
                                stdev*stdev + stdev2*stdev2);
                      fprintf(fp, "%e ", Ytemp2);
                    }
                  }
                  fprintf(fp, "\n");
                }
              }
            }

            nFail = 0;
            for (jj = 0; jj <= maxPts; jj++)
            {
              XGuess[ii] = xLower[ii]+jj*XRange[ii]/maxPts;
              Ytemp = constrPtr->evaluate(XGuess,XDist[jj],status);
              if (status == 0)
              {
                XDist[jj] = 0.0;
                nFail++;
              }

              ddata = 1.0;
              if (inputPDFs != NULL)
              {
                for (ii2 = 0; ii2 < nInputs; ii2++)
                {
                  if ((designParams == NULL || designParams[ii2] == 0) && 
                       inputPDFs[ii2] != NULL &&
                      (vecRSIndices.length() == 0 || vecRSIndices[ii2] >= 0))
                  {
                    inputPDFs[ii2]->getPDF(iOne,&XGuess[ii2],&ddata2);
                    ddata *= ddata2;
                  }
                }
              }
              SDist[jj] = ddata; 
            }

            ddata = XDist[0];
            for (jj = 1; jj <= maxPts; jj++) 
            {
              ddata2 = XDist[jj];
              if (ddata2 < ddata) ddata = ddata2;
            }
            XChains[iChain][chainCnt][nInputs] = ddata;
            for (jj = 0; jj <= maxPts; jj++)
            {
              XDist[jj] = SDist[jj] * exp(-0.5 * (XDist[jj] - ddata));

              ddata2 = XDist[jj] * exp(-0.5*ddata);
              if (ddata2 > Ymax)
              {
                Ymax = ddata2;
                for (ii2 = 0; ii2 < nInputs; ii2++) 
                  Xmax[ii2] = XGuess[ii2];
                Xmax[ii] = xLower[ii]+(xUpper[ii]-xLower[ii])*jj/maxPts;
              }
              if (jj > 0) XDist[jj] += XDist[jj-1];
            }
            if (masterOption == 1 && fp != NULL)
            {
              fprintf(fp,"];\n");
              fprintf(fp,"S = [\n");
              for (jj = 0; jj <= maxPts; jj++) 
                fprintf(fp, "%e\n", SDist[jj]);
              fprintf(fp,"];\n");
              fprintf(fp,"P = A(:,1) .* S;\n");
              fprintf(fp,"figure(1)\n");
              if (psPlotTool_ == 1)
                   fprintf(fp,"        set(gca(),\"auto_clear\",\"on\")\n");
              else fprintf(fp,"        hold off\n");
              fprintf(fp,"subplot(1,2,1)\n");
              fprintf(fp,"X = %e : %e : %e;\n",xLower[ii],
                      (xUpper[ii]-xLower[ii])/maxPts-1e-8,xUpper[ii]);
              fprintf(fp,"if length(X) > nPts\n");
              fprintf(fp,"   X = X(1:nPts);\n");
              fprintf(fp,"end;\n");
              fprintf(fp,"plot(X,P,'lineWidth',2.0)\n");
              fprintf(fp,"set(gca,'linewidth',2)\n");
              fprintf(fp,"set(gca,'fontweight','bold')\n");
              fprintf(fp,"set(gca,'fontsize',12)\n");
              fprintf(fp,"ylabel('-2 log likelihood','FontWeight','bold'");
              fprintf(fp,",'FontSize',12)\n");
              fprintf(fp,"title('Input %d','FontWeight','bold'",ii+1);
              fprintf(fp,",'FontSize',12)\n");
              fprintf(fp,"grid on\n");
              fprintf(fp,"box on\n");
              fprintf(fp,"subplot(1,2,2)\n");
              fprintf(fp,"P2 = exp(-0.5*P);\n");
              fprintf(fp,"for ii = 2 : %d\n",maxPts+1);
              fprintf(fp,"   P2(ii) = P2(ii) + P2(ii-1);\n");
              fprintf(fp,"end;\n");
              fprintf(fp,"P2 = P2 / P2(nPts);\n");
              fprintf(fp,"plot(X,P2,'lineWidth',2.0)\n");
              fprintf(fp,"set(gca,'linewidth',2)\n");
              fprintf(fp,"set(gca,'fontweight','bold')\n");
              fprintf(fp,"set(gca,'fontsize',12)\n");
              fprintf(fp,"title('Input %d','FontWeight','bold'",ii+1);
              fprintf(fp,",'FontSize',12)\n");
              fprintf(fp,"ylabel('likelihood CDF','FontWeight','bold'");
              fprintf(fp,",'FontSize',12)\n");
              fprintf(fp,"grid on\n");
              fprintf(fp,"box on\n");
            }

            if (printLevel > 3)
              printOutTS(PL_INFO,"proposal distribution max = %e\n", 
                         XDist[maxPts]);
            if (XDist[maxPts] - XDist[0] > 1.0e-16)
            {
              for (jj = 1; jj <= maxPts; jj++)
                XDist[jj] = (XDist[jj] - XDist[0]) / 
                            (XDist[maxPts]-XDist[0]);
              XDist[0] = 0;
              if (masterOption == 1 && fp != NULL)
              {
                fprintf(fp,"XCDF = [\n");
                for (jj = 0; jj <= maxPts; jj++) 
                  fprintf(fp, "%e\n", XDist[jj]);
                fprintf(fp,"];\n");
                fprintf(fp,"X = %e : %e : %e;\n",xLower[ii],
                          (xUpper[ii]-xLower[ii])/maxPts-1e-8,xUpper[ii]);
                fprintf(fp,"figure(2)\n");
                fprintf(fp,"plot(X,XCDF,'lineWidth',2.0)\n");
                fprintf(fp,"set(gca,'linewidth',2)\n");
                fprintf(fp,"set(gca,'fontweight','bold')\n");
                fprintf(fp,"set(gca,'fontsize',12)\n");
                fprintf(fp,"ylabel('likelihood CDF (check)','FontWeight',");
                fprintf(fp,"'bold','FontSize',12)\n");
                fprintf(fp,"title('Input %d','FontWeight','bold'",ii+1);
                fprintf(fp,",'FontSize',12)\n");
                fprintf(fp,"grid on\n");
                fprintf(fp,"box on\n");
              }
              XGuess[ii] = PSUADE_drand();
#if PS_INTERP == 0
              index = binarySearchDble(XGuess[ii], XDist, maxPts+1);
              if (index < 0) index = - index - 1;
              if      (index == maxPts) ddata = (double) index;
              else if (index == maxPts-1)
              {
                if (PABS(XDist[index]-XDist[maxPts]) > 1.0e-16)
                  ddata = index + (XGuess[ii]-XDist[index]) /
                                  (XDist[maxPts]-XDist[index]);
                else ddata = (double) index;
              }
              else
              {
                if (PABS(XDist[index+1]-XDist[index]) > 1.0e-16)
                  ddata=index+(XGuess[ii]-XDist[index])/
                              (XDist[index+1]-XDist[index]);
                else ddata = (double) index;
              }
              XGuess[ii] = xLower[ii]+ddata*XRange[ii]/maxPts;
#endif
#if PS_INTERP == 1
              index = binarySearchDble(XGuess[ii], XDist, maxPts+1);
              if (index < 0) index = - index - 1;
              if      (index == maxPts) ddata = (double) index;
              else if (index == 0)
              {
                if (PABS(XDist[index]-XDist[index+1]) > 1.0e-16)
                  ddata = index + (XGuess[ii]-XDist[index])/
                                  (XDist[index+1]-XDist[index]);
                else ddata = (double) index;
              }
              else
              {
                xd1 = XDist[index] - XDist[index-1];
                xd2 = XDist[index+1] - XDist[index-1];
                aa  = det * (b22 * xd1 - b12 * xd2);
                bb  = det * (b11 * xd2 - b21 * xd1);
                cc  = XDist[index-1] - XGuess[ii];
                ddata = bb * bb - 4.0 * aa *cc;
                if (ddata < 0 || aa == 0)
                {
                  ddata = (XGuess[ii]-XDist[index])/
                          (XDist[index+1]-XDist[index]);
                  ddata += (double) index;
                }
                else
                {
                  ddata = (-bb + sqrt(ddata)) / (2 * aa) + 
                           (index - 1.0)/maxPts;
                  if (ddata < (index-1.0)/maxPts || 
                      ddata > (index+1.0)/maxPts)
                  {
                    ddata = bb * bb - 4.0 * aa *cc;
                    if (ddata < 0)
                      ddata=(XGuess[ii]-XDist[index])/
                            (XDist[index+1]-XDist[index]);
                    else ddata = (-bb - sqrt(ddata)) / (2 * aa);
                    if (ddata < (index-1.0)/maxPts || 
                        ddata > (index+1.0)/maxPts)
                    {
                      printOutTS(PL_INFO,
                         "MCMC INFO: something is wrong in interpolation.\n");
                      printOutTS(PL_INFO,
                         "           Offending ddata = %e (%e, %e)\n",
                      ddata,(index-1.0)/maxPts,(index+1.0)/maxPts);
                      ddata=(XGuess[ii]-XDist[index])/
                                    (XDist[index+1]-XDist[index]);
                    }
                  }
                  ddata *= (double) maxPts;
                }
              }
              XGuess[ii] = xLower[ii]+ddata*XRange[ii]/maxPts;
#endif
#if PS_INTERP == 2
              status = faDist->initialize(ZDist, XDist);
              xbeg = 0; xend = 1.0;
              while ((xend-xbeg) > xtol)
              {
                xtrial = 0.5 * (xbeg + xend);
                ddata = faDist->evaluatePoint(&xtrial);
                if (PABS(ddata-XGuess[ii]) < xtol) break;
                if (XGuess[ii] > ddata) xbeg = xtrial;
                else                    xend = xtrial;
              }
              XGuess[ii] = xLower[ii]+xtrial*XRange[ii];
#endif
            }
            else 
            {
              XGuess[ii] = Xtemp;
              if (printLevel > 2)
                printOutTS(PL_INFO,
                     "MCMC iteration %7d : no modification in input %d\n",
                     mcmcIts+1,ii+1);
              if (nFail == maxPts+1) 
              {
                printOutTS(PL_ERROR,
                     "ERROR: Constraints have resulted in zero proposal.\n");
                exit(1);
              }
            }
            if (masterOption == 1 && fp != NULL)
            {
              if (psPlotTool_ == 1)
                   fprintf(fp,"        set(gca(),\"auto_clear\",\"off\")\n");
              else fprintf(fp,"        hold on\n");
              fprintf(fp,"XX = %e * ones(2,1);\n", XGuess[ii]);
              fprintf(fp,"YY = [0 1]';\n");
              fprintf(fp,"plot(XX,YY,'r-','linewidth',2)\n");
              if (psPlotTool_ == 1)
                   fprintf(fp,"        set(gca(),\"auto_clear\",\"on\")\n");
              else fprintf(fp,"        hold off\n");
              fclose(fp);
              fp = NULL; 
              for (ii2 = 0; ii2 < nInputs; ii2++)
              {
                printOutTS(PL_INFO,
                   "Next guess point %d = %e ",ii2+1,XGuess[ii2]);
                if (ii2 == ii) printOutTS(PL_INFO," ***\n");
                else           printOutTS(PL_INFO,"\n");
              }
              printf("Now examine MCMCDistTrack.m for diagnostics.\n");
              printf("Next, enter 0 to continue without any more\n");
              printf("stops, enter n (an integer > 1) to skip <n>\n");
              printf("iterations, or enter anything else to go to\n");
              printf("the next iteration.\n");
              scanf("%d", &ii2);
              if (ii2 == 0) masterOption = 0;
              if (ii2 > 1) masterCount = ii2 + 1;
            }
            if (masterOption == 1 && masterCount > 1) masterCount--;

            if (mcmcIts >= maxSamples/2 || globalIts > 0)
            {
              for (ii2 = 0; ii2 < nInputs; ii2++) 
              {
                XChains[iChain][chainCnt][ii2] = XGuess[ii2];
                ddata = (XGuess[ii2] - xLower[ii2]) / XRange[ii2];
                index = (int) (ddata * nbins);
                if (index >= nbins) index = nbins - 1;
                bins[index][ii2]++;
              }
              for (ii2 = 0; ii2 < nInputs; ii2++) 
              {
                ddata = (XGuess[ii2] - xLower[ii2]) / XRange[ii2];
                index = (int) (ddata * nbins);
                if (index >= nbins) index = nbins - 1;
                for (ii3 = 0; ii3 < nInputs; ii3++) 
                {
                  ddata2 = (XGuess[ii3] - xLower[ii3]) / XRange[ii3];
                  index2 = (int) (ddata2 * nbins);
                  if (index2 >= nbins) index2 = nbins - 1;
                  bins2[index][index2][ii2][ii3]++;
                }
              }
              chainCnt++;
            }
            mcmcIts++;
          }
          fp = fopen("psuade_stop", "r");
          if (fp != NULL)
          {
            printOutTS(PL_ERROR,
                 "MCMC INFO: psuade_stop FILE FOUND - TERMINATE MCMC.\n");
                fclose(fp);
            return 0.0;
          }
          if (mcmcIts >= maxSamples) break;
        }
      }
      if (countTrack <= 10) printOutTS(PL_INFO,"100%%\n");
      else                  printOutTS(PL_INFO,"\n");
      if (printLevel >= 0)
      {
        printOutTS(PL_INFO,
             "       Chain %d current final guess : \n",iChain+1);
        for (ii = 0; ii < nInputs; ii++)
          printOutTS(PL_INFO,"          Input %4d = %e\n",ii+1,XGuess[ii]);
      }
    } // ichain loop
 
    globalIts++;
    printOutTS(PL_INFO, "\nIteration %d summary: \n", globalIts);
    mcmcFail = nInputs - dnInputs;
    if (vecRSIndices.length() > 0)
    {
      for (ii = 0; ii < nInputs; ii++)
        if (vecRSIndices[ii] < 0) mcmcFail--; 
    }
    for (ii = 0; ii < nInputs; ii++)
    {
      if ((vecRSIndices.length() == 0 || 
          (vecRSIndices.length() > 0 && vecRSIndices[ii] >=0)) && 
          (designParams == NULL || designParams[ii] == 0))
      {
        if (printLevel > 2) printOutTS(PL_INFO, "Input = %d\n", ii+1);
       
        for (iChain = 0; iChain < numChains; iChain++)
        {
          ddata = 0.0;
          for (jj = 0; jj < chainCnt; jj++) 
            ddata += XChains[iChain][jj][ii];
          ddata /= chainCnt;
          ddata2 = 0.0;
          for (jj = 0; jj < chainCnt; jj++) 
             ddata2 += pow(XChains[iChain][jj][ii]-ddata,2.0);
          ddata2 /= (double) (chainCnt - 1);
          chainMeans[iChain] = ddata;
          chainStdevs[iChain] = ddata2;
          if (globalIts > 2 && chainStdevs[iChain] < 1.0e-20) 
          {
            printOutTS(PL_INFO,
                 "MCMC INFO: chain %d disabled.\n",iChain+1);
            chainStatus[iChain] = 1;
          }
        }
        nChainGood = 0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0) nChainGood++;
        }
        if (nChainGood <= 1)
        {
          printOutTS(PL_ERROR,"MCMC ERROR: too few chains <= 1.\n");
          printOutTS(PL_ERROR,
             "Suggestion: You may want to relax the experimental data\n");
          printOutTS(PL_ERROR,"     uncertainties (make them larger).\n");
          printOutTS(PL_ERROR,
             "     To see if this is the problem, turn on printlevel\n");
          printOutTS(PL_ERROR,
             "     to 3 and run again. If the variance of the chains\n");
          printOutTS(PL_ERROR,
             "     are small, small data uncertainties is probably the\n");
          printOutTS(PL_ERROR,"     problem.\n");
          exit(1);
        }
        WStat = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0) WStat += chainStdevs[iChain];
        }
        WStat /= (double) nChainGood;
        if (WStat < 0) WStat = PSUADE_UNDEFINED;
        if (printLevel > 2) 
          printOutTS(PL_INFO,"  Within  chain variance W = %e\n", WStat);
        ddata = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0) ddata += chainMeans[iChain];
        }
        ddata /= (double) nChainGood;
        BStat = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0)
            BStat += pow(chainMeans[iChain]-ddata,2.0);
        }
        BStat = BStat / (nChainGood - 1.0) * chainCnt;
        if (printLevel > 2) 
          printOutTS(PL_INFO,
               "  Between chain variance B = %e\n", BStat/chainCnt);
        ddata = (1 - 1.0/chainCnt) * WStat + BStat / chainCnt;
        ddata = ddata / WStat * (numChains + 1) / numChains - 
                  (chainCnt - 1.0) / (double) (chainCnt * numChains); 
        if (ddata < 0) ddata2 = PSUADE_UNDEFINED;
        else           ddata2 = sqrt(ddata);
        if (printLevel > 2)
        {
          for (iChain = 0; iChain < numChains; iChain++)
            printOutTS(PL_INFO,"  Chain %d mean, var = %e %e\n",iChain+1,
                    chainMeans[iChain]*XRange[ii]+xLower[ii],
                    chainStdevs[iChain]*XRange[ii]*XRange[ii]);
          printOutTS(PL_INFO,"  Chain length             = %d\n",chainCnt);
          printOutTS(PL_INFO,"  Weighted average of B, W = %e\n", ddata);
        }
        printOutTS(PL_INFO,"  Input %d PSRF = %e\n", ii+1, ddata2);
        psrfs[ii] = ddata2;
        if (ddata2 < psrfThreshold)
        {
          printOutTS(PL_INFO,"MCMC INFO : PSRF < %e ==> converged.\n",
                     psrfThreshold);
          mcmcFail--;
        }
        ddata = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0)
            for (jj = 0; jj < chainCnt; jj++)
               ddata += XChains[iChain][jj][ii];
        }
        ddata /= (double) (nChainGood * chainCnt);
        means_[ii] = ddata;
        ddata2 = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0)
            for (jj = 0; jj < chainCnt; jj++)
              ddata2 += pow(XChains[iChain][jj][ii]-ddata,2.0);
        }
        ddata2 /= (double) (chainCnt*nChainGood-1);
        sigmas_[ii] = sqrt(ddata2);
      }
    }

    if (mcmcFail == 0 && printLevel > 3 && s2Analyzer != NULL)
    {
      if (vecRSIndices.length() > 0)
      {
        for (ii = 0; ii < nInputs; ii++)
          if (vecRSIndices[ii] < 0) mcmcFail--; 
      }
      for (ii = 0; ii < nInputs; ii++)
      {
        if ((vecRSIndices.length() == 0 || 
            (vecRSIndices.length() > 0 && vecRSIndices[ii] >=0)) && 
            (designParams == NULL || designParams[ii] == 0))
        {
          printOutTS(PL_INFO, "Geweke Input = %d\n", ii+1);
     
          cnt = chainCnt / 2;
          ddata2 = (double) numChains;
          for (iChain = 0; iChain < numChains; iChain++)
          {
            for (jj = 0; jj < 2*cnt; jj++) 
              s2Vec[jj] = XChains[iChain][chainCnt-2*cnt+jj][ii];
            ddata = s2Analyzer->TAnalyze(cnt,s2Vec,cnt,&s2Vec[cnt],1);
           }
        }
      }
    }

    for (ii = 0; ii < nInputs; ii++) 
    {
      if ((vecRSIndices.length() == 0 || 
          (vecRSIndices.length() > 0 && vecRSIndices[ii] >=0)) && 
          (designParams == NULL || designParams[ii] == 0))
      {
        printOutTS(PL_INFO,
             "MCMC: input %3d value at peak of likelihood = %e\n",
                   ii+1, Xmax[ii]);
        ddata = means_[ii];
        printOutTS(PL_INFO,"MCMC: input %3d mean    = %e\n", ii+1, ddata);
        ddata = sigmas_[ii];
        printOutTS(PL_INFO,"MCMC: input %3d std dev = %e\n", ii+1, ddata);
        mostLikelyInput_[ii] = Xmax[ii];
      }
    }

    if (masterOption == 2)
    {
      fp = fopen("MCMCChainHistogram.m", "w");
      if (fp != NULL)
      {
        for (ii = 0; ii < nInputs; ii++)
        {
          if ((vecRSIndices.length() == 0 || vecRSIndices[ii] >= 0) &&
              (designParams == NULL || designParams[ii] == 0)) 
          {
            fprintf(fp,"X%d = [\n",ii+1);
            for (jj = 0; jj < chainCnt; jj++)
            {
              for (iChain = 0; iChain < numChains; iChain++)
                 fprintf(fp,"%e ", XChains[iChain][jj][ii]);
              fprintf(fp,"\n");
            }
            fprintf(fp,"];\n");
            fprintf(fp,"[nk, nx] = hist(X%d);\n", ii+1);
            fprintf(fp,"nk = nk / %d;\n", chainCnt);
            fprintf(fp,"bar(nx, nk)\n");
            fprintf(fp,"%%plot(X%d)\n", ii+1);
            sprintf(charString, "Input %d", ii+1);
            fwritePlotTitle(fp, charString);
            fwritePlotXLabel(fp, "Input Value");
            fwritePlotYLabel(fp, "Input Count");
            fwritePlotAxes(fp);
            fprintf(fp,
               "text(0.05,0.9,'colors are for different chains','sc')\n");
            fprintf(fp,
               "disp('Colors in histograms are for different chains.\n");
            fprintf(fp,"disp('Press enter to continue to next input')\n");
            fprintf(fp,"pause\n");
          }
        }
        fclose(fp);
        fp = NULL;
        printOutTS(PL_INFO,
           "MCMC: A file called MCMCChainHistogram.m has been created.\n");
        printOutTS(PL_INFO,
           "      Use Matlab to view the histograms of all chains for\n");
        printOutTS(PL_INFO,"      all inputs to assess convergence.\n");
        sprintf(charString, "Enter 1 to continue or 0 to terminate : ");
        ii = getInt(0, 10, charString);
        if (ii == 0) break;
      }
      else
      {
        printOutTS(PL_INFO,
             "MCMC INFO: cannot create MCMCChainHistogram.m file\n");
      } 
    } 

    if (mcmcFail == 0) break;

    genMatlabFile(nInputs,xLower,xUpper,XRange,nPlots,
             vecPlotIndices.getIVector(),nbins,
             NULL,NULL,bins,bins2,qData,numChains,chainCnt,XChains,
             chainStatus,Xmax,0);

    fp = fopen("psuade_stop", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_stop FILE FOUND - TERMINATE MCMC.\n");
      fclose(fp);
      fp = NULL;
      strcpy(charString, "psuade_stop");
      unlink(charString);
      break;
    }

    fp = fopen("psuade_nogm", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_nogm FILE FOUND. GM mode is now off.\n");
      fclose(fp);
      fp = NULL;
      psGMMode_ = 0;
      masterOption = 0;
      strcpy(charString, "psuade_nogm");
      unlink(charString);
    }
    fp = fopen("psuade_gm", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_gm FILE FOUND. GM mode is now on.\n");
      fclose(fp);
      fp = NULL;
      psGMMode_ = 1;
      printf("Please choose from the following option (or 0 if none): \n");
      printf("  1. track proposal distribution at each MCMC iteration\n");
      printf("  2. track posterior distributions after each MCMC cycle\n");
      scanf("%s", charString);
      fgets(lineIn,1000,stdin);
      if (charString[0] == '1') masterOption = 1;
      if (charString[0] == '2') masterOption = 2;
      strcpy(charString, "psuade_gm");
      unlink(charString);
    }
    fp = fopen("psuade_print", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_print FILE FOUND. Print level is set to 3.\n");
      fclose(fp);
      fp = NULL;
      printLevel = 3;
      strcpy(charString, "psuade_print");
      unlink(charString);
    }
  } // while loop

  if (globalIts >= maxGlobalIts)
  {
    mcmcFail = 0;
    for (ii = 0; ii < nInputs; ii++) 
    {
      if ((vecRSIndices.length() == 0 || 
          (vecRSIndices.length() > 0 && vecRSIndices[ii] >=0)) && 
          (designParams == NULL || designParams[ii] == 0))
         if (psrfs[ii] > psrfThreshold) mcmcFail = 1;
    }
    if (mcmcFail == 1) 
      printOutTS(PL_INFO,
           "MCMC maximum iterations exceeded but no convergence.\n");
  }
  else printOutTS(PL_INFO, "MCMC iterations completed\n");

  for (ii = 0; ii < nInputs; ii++) 
    for (jj = 0; jj < nbins; jj++) bins[jj][ii] = 0;
  for (ii = 0; ii < nInputs; ii++) 
    for (ii2 = 0; ii2 < nInputs; ii2++) 
      for (jj = 0; jj < nbins; jj++)
        for (jj2 = 0; jj2 < nbins; jj2++)
          bins2[jj][jj2][ii][ii2] = 0;
  for (iChain = 0; iChain < numChains; iChain++) 
  { 
    if (chainStatus[iChain] == 0)
    {
      for (jj = 0; jj < chainCnt; jj++) 
      {
        for (ii2 = 0; ii2 < nInputs; ii2++) 
        {
          ddata = XChains[iChain][jj][ii2];
          ddata = (ddata - xLower[ii2]) / (xUpper[ii2] - xLower[ii2]);
          index = (int) (ddata * nbins);
          if (index >= nbins) index = nbins - 1;
          bins[index][ii2]++;
        }
        for (ii2 = 0; ii2 < nInputs; ii2++) 
        {
          ddata = XChains[iChain][jj][ii2];
          ddata = (ddata - xLower[ii2]) / (xUpper[ii2] - xLower[ii2]);
          index = (int) (ddata * nbins);
          if (index >= nbins) index = nbins - 1;
          for (ii3 = 0; ii3 < nInputs; ii3++) 
          {
            ddata2 = XChains[iChain][jj][ii3];
            ddata2 = (ddata2 - xLower[ii3]) / (xUpper[ii3] - xLower[ii3]);
            index2 = (int) (ddata2 * nbins);
            if (index2 >= nbins) index2 = nbins - 1;
            bins2[index][index2][ii2][ii3]++;
          }
        }
      }
    }
  }

  delete [] XDist;
  delete [] SDist;
  delete [] XGuess;
  delete [] XGuessS;
  delete [] YGuessS;
  delete [] YGuessStds;
  delete [] XDesignS;
  delete [] YDesignS;
  delete [] YDesignStds;
 
  dataPtr->getParameter("method_sampling", pPtr);
  int methodSave = pPtr.intData_;
  dataPtr->updateMethodSection(PSUADE_SAMP_MC,-1,-1,-1,-1);
  PDFManager *pdfman = new PDFManager();
  pdfman->initialize(dataPtr);
  psVector vecLB, vecUB, vecOut;
  vecLB.load(nInputs, xLower);
  vecUB.load(nInputs, xUpper);
  int nSamps = 100000;
  vecOut.setLength(nSamps*nInputs);
  pdfman->genSample(nSamps, vecOut, vecLB, vecUB);
  dataPtr->updateMethodSection(methodSave,-1,-1,-1,-1);
  delete pdfman;

  int **pbins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    pbins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) pbins[ii][jj] = 0;
  }
  int ****pbins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    pbins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      pbins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        pbins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          pbins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }
  double ddata1;
  for (ii = 0; ii < nInputs; ii++)
  {
    ddata1 = nbins / (xUpper[ii] - xLower[ii]);
    for (jj = 0; jj < nSamps; jj++)
    {
      ddata = vecOut[jj*nInputs+ii];
      ddata -= xLower[ii];
      ddata *= ddata1;
      kk = (int) ddata;
      if (kk >= nbins) kk = nbins - 1; 
      pbins[kk][ii]++;
    }
  }
  for (ii = 0; ii < nInputs; ii++)
  {
    ddata1 = nbins / (xUpper[ii] - xLower[ii]);
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      ddata2 = nbins / (xUpper[ii2] - xLower[ii2]);
      for (jj = 0; jj < nSamps; jj++)
      {
        ddata = vecOut[jj*nInputs+ii];
        ddata -= xLower[ii];
        ddata *= ddata1;
        kk = (int) ddata;
        if (kk >= nbins) kk = nbins - 1; 
        ddata = vecOut[jj*nInputs+ii2];
        ddata -= xLower[ii2];
        ddata *= ddata2;
        kk2 = (int) ddata;
        if (kk2 >= nbins) kk2 = nbins - 1; 
        pbins2[kk][kk2][ii][ii2]++;
      }
    }
  }
  genMatlabFile(nInputs,xLower,xUpper,XRange,nPlots,
        vecPlotIndices.getIVector(),nbins,
        pbins,pbins2,bins,bins2,qData,numChains,chainCnt,XChains,
        chainStatus, Xmax,0);
  for (ii = 0; ii < nbins; ii++) delete [] pbins[ii];
  delete [] pbins;
  for (jj = 0; jj < nbins; jj++)
  {
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      for (ii = 0; ii < nInputs; ii++) delete [] pbins2[jj][jj2][ii];
      delete [] pbins2[jj][jj2];
    }
    delete [] pbins2[jj];
  }
  delete [] pbins2;

  if (genPosteriors == 1)
  {
    cnt = nChainGood * chainCnt;
    if (cnt > 200000) cnt = 200000;
    cnt /= nChainGood;
    genPostLikelihood(nInputs,xLower,xUpper,XRange,numChains, chainCnt,
                  XChains, chainStatus, cnt, vecRSIndices.getIVector(), 
                  vecRSValues.getDVector(),
                  designParams, dnInputs, dnSamples, dSamInputs, faPtrs,
                  faPtrs1, nOutputs, discOutputs, discFuncConstantMeans,
                  dSamMeans, dSamStdevs);
  }

  if (genPosteriors == 1)
  {
    fp = fopen("MCMCPostSample", "w");
    if (fp != NULL)
    {
      fprintf(fp, "PSUADE_BEGIN\n");
      cnt = nChainGood * chainCnt;
      if (cnt > 50000) cnt = 50000;
      cnt /= nChainGood;
      fprintf(fp, "%d %d\n", cnt*nChainGood,nInputs);
      if (qData.strArray_ != NULL)
      {
        fprintf(fp, "# ");
        for (jj = 0; jj < nInputs; jj++)
          fprintf(fp,"%s ", qData.strArray_[jj]);
        fprintf(fp, "\n");
      }
      ii2 = 0;
      for (iChain = 0; iChain < numChains; iChain++)
      { 
        if (chainStatus[iChain] == 0)
        {
          for (ii = chainCnt-cnt; ii < chainCnt; ii++)
          {
            fprintf(fp, "%d ", ii2+1);
            for (jj = 0; jj < nInputs; jj++)
            {
              if ((vecRSIndices.length() == 0 || vecRSIndices[jj] >= 0) &&
                  (designParams == NULL || designParams[jj] == 0)) 
              {
                ddata = XChains[iChain][ii][jj];
                fprintf(fp, "%e ", ddata);
              }
              else if (vecRSIndices.length() > 0 && vecRSIndices[jj] < 0)
                fprintf(fp, "%e ", vecRSValues[jj]);
              else if (designParams != NULL && designParams[jj] != 0) 
                fprintf(fp, "%e ", 0.5 * (xUpper[jj] + xLower[jj]));
            }
            fprintf(fp, "\n");
            ii2++;
          }
        }
      }
      fprintf(fp, "PSUADE_END\n");
      fprintf(fp, "#N=%d;\n",nChainGood*cnt);
      fprintf(fp, "#m=%d;\n",cnt);
      for (iChain = 0; iChain < numChains; iChain++)
      {
        if (chainStatus[iChain] == 0)
          fprintf(fp, "#A%d = A(%d*m+1:%d*m,:);\n",iChain,iChain,iChain+1);
      }
      fprintf(fp, "#for ii = 2 : %d\n", nInputs+1);
      for (iChain = 0; iChain < numChains; iChain++)
      {
        if (chainStatus[iChain] == 0)
        {
          fprintf(fp, "#subplot(*,*,%d)\n",iChain+1);
          fprintf(fp, "#hist(A%d(:,ii))\n",iChain+1);
        }
        fprintf(fp, "#ii-1\n");
        fprintf(fp, "#pause;\n");
      }
      fprintf(fp, "#end;\n");
      fclose(fp);
    }
    printOutTS(PL_INFO,
         "MCMC: 'MCMCPostSample' file has a posterior sample.\n");
  }

  int    nInps, nOuts, nSams, *states;
  double *allOuts;
  char   **oNames;
  PsuadeData *filePtr1, *filePtr2;
  if (modelFormFlag == 1)
  {
    sprintf(charString, "psDiscrepancyModel1");
    filePtr1 = new PsuadeData();
    status = filePtr1->readPsuadeFile(charString);
    if (status != 0)
    {
      printOutTS(PL_ERROR,
        "MCMC ERROR: cannot read file %s in PSUADE format.\n",charString);
      exit(1);
    } 
  }
  if (modelFormFlag == 1 && status == 0)
  {
    filePtr1->getParameter("input_ninputs", pPtr);
    nInps = pPtr.intData_;
    filePtr1->getParameter("output_noutputs", pPtr);
    nOuts = pPtr.intData_;
    filePtr1->getParameter("method_nsamples", pPtr);
    nSams = pPtr.intData_;
    filePtr1->getParameter("output_sample", pOutputs);
    unlink(charString);
    allOuts = new double[nOutputs * nSams];
    for (jj = 0; jj < nSams; jj++) 
      allOuts[jj*nOutputs] = pOutputs.dbleArray_[jj];
    pOutputs.clean();
    oNames = new char*[nOutputs];
    oNames[0] = new char[100];
    sprintf(oNames[0], "Y1");
    for (ii = 1; ii < nOutputs; ii++)
    {
      filePtr2 = new PsuadeData();
      sprintf(charString, "psDiscrepancyModel%d", ii+1);
      status = filePtr2->readPsuadeFile(charString);
      if (status != 0) break;
      filePtr2->getParameter("input_ninputs", pPtr);
      if (pPtr.intData_ != nInps) break;
      filePtr2->getParameter("output_noutputs", pPtr);
      if (pPtr.intData_ != nOuts) break;
      filePtr2->getParameter("method_nsamples", pPtr);
      if (pPtr.intData_ != nSams) break;
      filePtr2->getParameter("output_sample", pOutputs);
      delete filePtr2;
      unlink(charString);
      for (jj = 0; jj < nSams; jj++) 
        allOuts[jj*nOutputs+ii] = pOutputs.dbleArray_[jj];
      pOutputs.clean();
      oNames[ii] = new char[100];
      sprintf(oNames[ii], "Y%d", ii+1);
    }
    if (nOutputs == 1)
    { 
      sprintf(charString, "psDiscrepancyModel");
      filePtr1->writePsuadeFile(charString, 0);
    }
    else if (ii == nOutputs)
    {
      states = new int[nSams];
      for (jj = 0; jj < nSams; jj++) states[jj] = 1;
      filePtr1->updateOutputSection(nSams,nOutputs,allOuts,states,oNames);
      sprintf(charString, "psDiscrepancyModel");
      filePtr1->writePsuadeFile(charString, 0);
      printOutTS(PL_INFO,
           "MCMC INFO: a sample (inputs/outputs) the discrepancy model\n");
      printOutTS(PL_INFO,"           is now in psDiscrepancyModel.\n");
      delete [] states;
      for (ii = 1; ii < nOutputs; ii++) delete [] oNames[ii];
    }
    else
    {
      printOutTS(PL_INFO,
           "MCMC INFO: unsuccessful creation of discrepancy sample file\n");
    }
    delete [] oNames[0];
    delete [] oNames;
    delete filePtr1;
    delete [] allOuts;
  }
  
  // final clean up
  for (ii = 0; ii < numChains; ii++)
  {
    for (jj = 0; jj < maxGlobalIts*maxSamples; jj++)
      delete [] XChains[ii][jj];
    delete [] XChains[ii];
  }
  delete [] XChains;
  if (inputPDFs != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
      if (inputPDFs[ii] != NULL) delete inputPDFs[ii];
    delete [] inputPDFs;
  }
  delete [] Xmax;
  delete [] XRange;
  delete [] psrfs;
  if (discOutputs != NULL) delete [] discOutputs;
  delete [] s2Vec;
  delete [] mcmcSeeds;
  if (s2Analyzer != NULL) delete s2Analyzer;
  if (discFuncConstantMeans != NULL) delete [] discFuncConstantMeans;
  if (discFuncConstantStds  != NULL) delete [] discFuncConstantStds;
  if (dSamMeans != NULL) delete [] dSamMeans;
  if (dSamStdevs != NULL) delete [] dSamStdevs;
#if 0
  delete [] ZDist;
  delete faDist;
#endif
  for (ii = 0; ii < nbins; ii++) delete [] bins[ii];
  delete [] bins;
  for (jj = 0; jj < nbins; jj++)
  {
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      for (ii = 0; ii < nInputs; ii++) delete [] bins2[jj][jj2][ii];
      delete [] bins2[jj][jj2];
    }
    delete [] bins2[jj];
  }
  delete [] bins2;
  if (dSamInputs != NULL) delete [] dSamInputs;
  delete constrPtr;
  if (faPtrs != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++) 
      if (faPtrs[ii] != NULL) delete faPtrs[ii];
    delete [] faPtrs;
  }
  if (faPtrs1 != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++) 
      if (faPtrs1[ii] != NULL) delete faPtrs1[ii];
    delete [] faPtrs1;
  }
  delete [] Ivec;
  return 0.0;
}

// ************************************************************************
// perform MCMC-like analysis (brute force)
// ------------------------------------------------------------------------
double MCMCAnalyzer::analyze_bf(aData &adata)
{
  int    ii, ii2, jj, jj2, kk, kk2, status, cnt, iOne=1, iZero=0, nInputs;
  int    nOutputs, nSamples, nbins, printLevel, faType, genPosteriors=0;
  int    *pdfFlags, maxSamples, modelFormFlag=0, noPDF, combFlag=0;
  int    nPlots, *plotIndices=NULL, *designParams=NULL;
  int    *rsIndices=NULL, dnSamples=0, dnInputs=0, rsErrFlag=0;
  double *dSamInputs=NULL, *dSamMeans=NULL, *dSamStdevs=NULL, dstatus;
  double *X=NULL, *Y=NULL, *lower=NULL, *upper=NULL, *rsValues=NULL;
  double *discFuncConstantMeans=NULL, *discFuncConstantStds=NULL;
  char   lineIn[1001], charString[1001], *rsFile=NULL;
  FILE   *fp=NULL;
  pData      pPtr, qData, pOutputs;
  FuncApprox **faPtrs=NULL, **faPtrs1=NULL;
  PsuadeData *dataPtr=NULL;
  RSConstraints *constrPtr;

  // display header 
  printLevel  = adata.printLevel_;
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO,"*              Brute Force Optimizer\n");
  printEquals(PL_INFO, 0);
  if (printLevel > 0)
  {
    printOutTS(PL_INFO,"TO GAIN ACCESS TO DIFFERENT OPTIONS: TURN ON\n\n");
    printOutTS(PL_INFO," * ana_expert to finetune MCMC parameters, \n");
    printOutTS(PL_INFO,"   (e.g. sample size for burn-in can be adjusted).\n");
    printOutTS(PL_INFO,
         " * rs_expert to customize response surface for MCMC,\n");
    printOutTS(PL_INFO," * printlevel 3 to display more diagnostics info.\n");
    printOutTS(PL_INFO," * printlevel 4 to display even more diagnostics.\n");
    printOutTS(PL_INFO," * printlevel >=5 reserved only for expert only.\n");
    printDashes(PL_INFO,0);
    printOutTS(PL_INFO,"FEATURES AVAILABLE IN THE CURRENT VERSION OF MCMC:\n");
    printOutTS(PL_INFO,
         " * Support likelihood functions from multiple outputs\n");
    printOutTS(PL_INFO," * Option to include response surface errors for\n");
    printOutTS(PL_INFO,"   polynomial regressions, bootstrapped MARS, and\n");
    printOutTS(PL_INFO,"   boostrapped RBF and Gaussian process.\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode.\n");
    printOutTS(PL_INFO,
         " * Option to include model form errors in the form of\n");
    printOutTS(PL_INFO,"   discrepancy models.\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode.\n");
    printOutTS(PL_INFO," * Option to set some inputs as design parameters\n");
    printOutTS(PL_INFO,
         "   - to be specified in the observation data spec file\n");
    printOutTS(PL_INFO,
         " * Option to disable some parameters (set to default)\n");
    printOutTS(PL_INFO,
         "   - in case these parameters are not to be calibrated\n");
    printOutTS(PL_INFO,
         "   - use rs_index_file in PSUADE's ANALYSIS section\n");
    printOutTS(PL_INFO,"   - not available with discrepancy modeling\n");
    printOutTS(PL_INFO," * Option to generate a posterior sample\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode\n");
    printOutTS(PL_INFO,
         " * MCMC can be terminated gracefully by creating a file\n");
    printOutTS(PL_INFO,
         "   named 'psuade_stop' in the same directory during the\n");
    printOutTS(PL_INFO,"   run (in case it takes too long).\n");
    printEquals(PL_INFO, 0);
  }

  // extract data from aData object (passed in from outside)
  nInputs     = adata.nInputs_;
  nInputs_    = nInputs;
  nOutputs    = adata.nOutputs_;
  nOutputs_   = nOutputs;
  nSamples    = adata.nSamples_;
  X           = adata.sampleInputs_;
  Y           = adata.sampleOutputs_;
  lower       = adata.iLowerB_;
  upper       = adata.iUpperB_;
  pdfFlags    = adata.inputPDFs_;
  noPDF       = 1;
  if (pdfFlags != NULL)
    for (ii = 0; ii < nInputs; ii++) if (pdfFlags[ii] != 0) noPDF = 0;

  dataPtr = adata.ioPtr_;
  if (dataPtr != NULL) dataPtr->getParameter("input_names", qData);

  // clean up
  if (means_) delete[] means_;
  if (sigmas_) delete[] sigmas_;
  if (mostLikelyInput_) delete[] mostLikelyInput_;
  if (mostLikelyOutput_) delete[] mostLikelyOutput_;
  means_ = NULL;
  sigmas_ = NULL;
  mostLikelyInput_ = NULL;
  mostLikelyOutput_ = NULL;

  if (dataPtr != NULL)
  {
    dataPtr->getParameter("ana_rstype", pPtr);
    faType = pPtr.intData_;
    dataPtr->getParameter("ana_rsindexfile", pPtr);
    rsFile = pPtr.strArray_[0];
    if (strcmp(rsFile, "NONE"))
    {
      printOutTS(PL_INFO,
           "A response surface index file has been specified.\n");
      fp = fopen(rsFile, "r");
      if (fp == NULL)
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: rs_index_file %s not found.\n",rsFile);
        return PSUADE_UNDEFINED;
      }
      else
      {
        printOutTS(PL_INFO,"INFO: rs_index_file %s found.\n",rsFile);
        fscanf(fp,"%d", &kk);
        if (kk != nInputs)
        {
          printOutTS(PL_ERROR,
             "MCMC ERROR: invalid nInputs in rs_index_file (%d != %d).\n",
             kk, nInputs);
          printOutTS(PL_ERROR,"  Data format should be: \n");
          printOutTS(PL_ERROR,
              "  line 1: nInputs in rs data (driver) file\n");
          printOutTS(PL_ERROR,
              "  line 2: 1 <1 or 0> <default value if 0, -1 if sample>\n");
          printOutTS(PL_ERROR,
              "  line 3: 2 <2 or 0> <0 if it is == 2>\n");
          printOutTS(PL_ERROR,
              "  line 4: 3 <3 or 0> <default value if 0, -1 if sample>\n");
          printOutTS(PL_ERROR,
              "  line 5: 4 <4 or 0> <0 if it is == 4>\n");
          printOutTS(PL_ERROR,"  ...\n");
          fclose(fp);
          return PSUADE_UNDEFINED;
        }
        rsIndices = new int[nInputs];
        rsValues = new double[nInputs];
        for (ii = 0; ii < nInputs; ii++) rsIndices[ii] = 0;
        for (ii = 0; ii < nInputs; ii++)
        {
          fscanf(fp, "%d", &kk);
          if (kk != ii+1)
          {
            printOutTS(PL_ERROR,
             "MCMC ERROR: 1st index in indexFile = %d (must be %d]).\n",
              kk, ii+1);
            printOutTS(PL_ERROR,"  Data format should be: \n");
            printOutTS(PL_ERROR,
             "  line 1: nInputs in rs data (driver) file\n");
            printOutTS(PL_ERROR,
             "  line 2: 1 <1 or 0> <default value if 0, -1 if sample>\n");
            printOutTS(PL_ERROR,
             "  line 3: 2 <2 or 0> <0 if it is == 2>\n");
            printOutTS(PL_ERROR,
             "  line 4: 3 <3 or 0> <default value if 0, -1 if sample>\n");
            printOutTS(PL_ERROR,
             "  line 5: 4 <4 or 0> <0 if it is == 4>\n");
            printOutTS(PL_ERROR,"  ...\n");
            fclose(fp);
            delete [] rsIndices;
            delete [] rsValues;
            return PSUADE_UNDEFINED;
          }
          fscanf(fp, "%d", &rsIndices[ii]);

          if (rsIndices[ii] == 999)
          {
            printOutTS(PL_ERROR,
                "MCMC_bf INFO: input %3d has a sample ==> switch to bf2.\n",
                ii+1);
            fclose(fp);
            delete [] rsIndices;
            delete [] rsValues;
            return (analyze_bf2(adata));
          }

          if (rsIndices[ii] < 0 || rsIndices[ii] > nInputs)
          {
            printOutTS(PL_ERROR,
                "MCMC_bf ERROR: input %3d = %d invalid\n",ii+1,
                rsIndices[ii]);
            fclose(fp);
            delete [] rsIndices;
            delete [] rsValues;
            return PSUADE_UNDEFINED;
          }
          rsIndices[ii]--;
          fscanf(fp, "%lg", &rsValues[ii]);
        }
        fclose(fp);
      }
    }
  }
  else if (mode_ == 0)
  {
    printOutTS(PL_INFO,
         "MCMC INFO: since ioPtr=NULL, assume MARS as reponse surface.\n");
    faType = PSUADE_RS_MARS;
  }

  // get experimental data information from the spec file
  dstatus = readSpecFile(nInputs, nOutputs, &dnSamples, &dnInputs, 
                    &designParams, &dSamInputs, &dSamMeans, &dSamStdevs,
                    combFlag, printLevel);
  if (dstatus != 0.0)
  {
    printf("MCMC ERROR: fail to read experimental data file.\n");
    return PSUADE_UNDEFINED;
  }

  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO,
        "*** OPTION TO INCLUDE RESPONSE SURFACE UNCERTAINTIES:\n");
    printOutTS(PL_INFO,
        "\nTo incorporate response surface uncertainties into the\n");
    printOutTS(PL_INFO,
        "likelihood function, make sure stochastic response surfaces\n");
    printOutTS(PL_INFO,
        "are used (GP/Kriging, polynomial regression, or bootstrapped\n");
    printOutTS(PL_INFO,
        "methods). Otherwise, no RS uncertainties will be included.\n");
    printOutTS(PL_INFO,
        "NOTE: if you don't know what this is, just say no below.\n");
    printf( "===> Include response surface uncertainties? (y or n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') rsErrFlag = 1;
    printEquals(PL_INFO, 0);
  }

  if (rsIndices != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
    {
      if (rsIndices[ii] == 0)
        printOutTS(PL_INFO,"MCMC INFO: input %3d inactive\n",ii+1);

      if (rsIndices[ii] == -1 && designParams != NULL && 
          designParams[ii] == 1)
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: inactive input %d cannot be design parameter\n",
             ii+1);
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        delete [] rsIndices;
        delete [] rsValues;
        return PSUADE_UNDEFINED;
      }

      if (rsIndices[ii] == 999 && 
          (designParams != NULL && designParams[ii] == 1))
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: inactive input %d cannot be uncertain parameter\n",
             ii+1);
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        delete [] rsIndices;
        delete [] rsValues;
        return PSUADE_UNDEFINED;
      }

      if (rsIndices[ii] < -1 || rsIndices[ii] > nInputs)
      {
        printOutTS(PL_ERROR,
             "MCMC INFO: input %3d = %d invalid\n",ii+1,rsIndices[ii]);
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        delete [] rsIndices;
        delete [] rsValues;
        return PSUADE_UNDEFINED;
      }
    }

    printOutTS(PL_INFO, "Response surface index information: \n");
    for (ii = 0; ii < nInputs; ii++)
    {
      if (rsIndices[ii] == -1)
        printOutTS(PL_INFO, "Input %4d: fixed at default value  = %e\n",
                   ii+1, rsValues[ii]);
      else if (rsIndices[ii] >= 1000)
        printOutTS(PL_INFO, "Input %4d: uncertain, sample index = %4d\n",
                   ii+1, rsIndices[ii]-999);
      else if (designParams != NULL && designParams[ii] == 1)
        printOutTS(PL_INFO, "Input %4d: design parameter\n", ii+1);
      else
        printOutTS(PL_INFO, "Input %4d: calibration parameter\n",ii+1);
    }
  }

  printOutTS(PL_INFO,
       "MCMC INFO: CREATING RESPONSE SURFACES FOR ALL OUTPUTS.\n");
  double *YY=NULL;
  faPtrs = new FuncApprox*[nOutputs];
  if (nSamples > 0)
  {
    YY = new double[nSamples];
    for (ii = 0; ii < nOutputs; ii++)
    {
      faType = -1;
      printOutTS(PL_INFO,
           "MCMC INFO: CREATING RESPONSE SURFACE FOR OUTPUT %d.\n",ii+1);
      faPtrs[ii] = genFA(faType, nInputs, iZero, nSamples);
      faPtrs[ii]->setNPtsPerDim(16);
      faPtrs[ii]->setBounds(lower, upper);
      faPtrs[ii]->setOutputLevel(0);
      for (kk = 0; kk < nSamples; kk++) YY[kk] = Y[kk*nOutputs+ii];

      status = faPtrs[ii]->initialize(X, YY);
      if (status != 0)
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: Unable to create response surface.\n");
        printOutTS(PL_ERROR,"            Consult PSUADE developers.\n");
        for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
        delete [] faPtrs;
        delete [] YY;
        if (rsValues  != NULL) delete [] rsValues;
        if (rsIndices != NULL) delete [] rsIndices;
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        return PSUADE_UNDEFINED;
      }
    }
    delete [] YY;
  }
  else
  {
    int        nInps, nSamps, nOuts;
    double     *XX, *tlower, *tupper;
    PsuadeData *ioOutFile;
    for (ii = 0; ii < nOutputs; ii++)
    {
      faType = -1;
      sprintf(charString,"Enter file name for output %d : ", ii+1);
      getString(charString, lineIn);
      kk = strlen(lineIn);
      lineIn[kk-1] = '\0';
      ioOutFile = new PsuadeData;
      status = ioOutFile->readPsuadeFile(lineIn);
      ioOutFile->getParameter("input_ninputs", pPtr);
      nInps = pPtr.intData_;
      ioOutFile->getParameter("output_noutputs", pPtr);
      nOuts = pPtr.intData_;
      if (nInps != nInputs)
      {
        printOutTS(PL_ERROR,
         "MCMC ERROR: Unable to create response surface for output %d\n",
         ii+1);
        printOutTS(PL_ERROR,"            due to nInputs mismatch.\n");
        for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
        delete [] faPtrs;
        if (rsValues  != NULL) delete [] rsValues;
        if (rsIndices != NULL) delete [] rsIndices;
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        return PSUADE_UNDEFINED;
      }
      if (nOuts != 1)
      {
        printOutTS(PL_ERROR,
         "MCMC ERROR: Unable to create response surface for output %d\n",
         ii+1);
        printOutTS(PL_ERROR,"            due to nOutputs != 1.\n");
        for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
        delete [] faPtrs;
        if (rsValues  != NULL) delete [] rsValues;
        if (rsIndices != NULL) delete [] rsIndices;
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        return PSUADE_UNDEFINED;
      }
      ioOutFile->getParameter("method_nsamples", pPtr);
      nSamps = pPtr.intData_;
      ioOutFile->getParameter("input_sample", pPtr);
      XX = pPtr.dbleArray_;
      pPtr.dbleArray_ = NULL;
      ioOutFile->getParameter("output_sample", pPtr);
      YY = pPtr.dbleArray_;
      pPtr.dbleArray_ = NULL;
      ioOutFile->getParameter("input_lbounds", pPtr);
      tlower = pPtr.dbleArray_;
      pPtr.dbleArray_ = NULL;
      ioOutFile->getParameter("input_ubounds", pPtr);
      tupper = pPtr.dbleArray_;
      pPtr.dbleArray_ = NULL;
      faPtrs[ii] = genFA(faType, nInps, iZero, nSamps);
      faPtrs[ii]->setNPtsPerDim(16);
      faPtrs[ii]->setBounds(tlower, tupper);
      faPtrs[ii]->setOutputLevel(0);
      status = faPtrs[ii]->initialize(XX, YY);
      delete [] YY;
      delete [] XX;
      delete [] tlower;
      delete [] tupper;
      delete ioOutFile;
      if (status != 0)
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: Unable to create response surface.\n");
        printOutTS(PL_ERROR,"            Consult PSUADE developers.\n");
        for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
        delete [] faPtrs;
        if (rsValues  != NULL) delete [] rsValues;
        if (rsIndices != NULL) delete [] rsIndices;
        if (designParams != NULL) delete [] designParams;
        if (dSamInputs != NULL) delete [] dSamInputs;
        if (dSamMeans != NULL) delete [] dSamMeans;
        if (dSamStdevs != NULL) delete [] dSamStdevs;
        return PSUADE_UNDEFINED;
      }
    }
  }

  nbins = 20;
  maxSamples = 500000;
  if (nInputs >= 10) maxSamples = 1000000;
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"*** CURRENT DEFAULT PARAMETER SETTINGS : \n\n");
  printOutTS(PL_INFO,"Inference max sample size = %d\n",maxSamples);
  printOutTS(PL_INFO,"Posterior histogram nbins = %d\n",nbins);
  printOutTS(PL_INFO,
       "NOTE: histogram nBins  - resolution of histogram bar graph\n");
  printOutTS(PL_INFO,
       "Turn on ana_expert mode to change these default settings.\n\n");

  if (psAnaExpertMode_ == 1)
  {
    sprintf(charString,
            "Enter maximum inference sample size (500000 - 5000000): ");
    maxSamples = getInt(1000, 50000000, charString);
    if (maxSamples < 500000) maxSamples = 500000;
    if (nInputs >= 10) maxSamples = 1000000;
    sprintf(charString,"Enter the number of histogram bins (10 - 25) : ");
    nbins = getInt(10, 50, charString);
  }
  if (psAnaExpertMode_ == 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
        "*** OPTION TO CREATE POSTERIOR PLOTS FOR SELECTED INPUTS:\n\n");
    printOutTS(PL_INFO,
        "MCMC will create MATLAB files for the posterior distributions.\n");
    printOutTS(PL_INFO,
        "You can choose to generate posterior plots for all inputs, or \n");
    printOutTS(PL_INFO,
        "just a selected few (in case there are too many inputs).\n");
    printf("Select inputs for which posterior plots are to be generated.\n");
    sprintf(charString,"Enter input number (-1 for all, 0 to terminate) : ");
    kk = 1;
    plotIndices = new int[nInputs];
    nPlots = 0;
    while (kk != 0 || nPlots < 1)
    {
      kk = getInt(-1, nInputs, charString);
      if (kk == -1)
      {
        for (ii = 0; ii < nInputs; ii++)
        {
          if (rsIndices == NULL || rsIndices[ii] >= 0)
            if (designParams == NULL || designParams[ii] == 0)
              plotIndices[nPlots++] = ii;
        }
        break;
      }
      if (kk != 0)
      {
        if (rsIndices != NULL && rsIndices[kk-1] < 0)
          printOutTS(PL_ERROR,
              "Input %d has been fixed by the rs index file (no plot).\n",
              kk+1);
        else if (designParams != NULL && designParams[kk-1] == 1)
          printOutTS(PL_ERROR,
              "Input %d is a design parameter (no plot)\n",kk);
        else
          plotIndices[nPlots++] = kk - 1;
      }
      if (kk == 0 && nPlots == 0)
        printOutTS(PL_ERROR,
            "You need to set at least 1 input for plotting posteriors.\n");
    }
    if (nPlots > 1) sortIntList(nPlots, plotIndices);
  }
  else
  {
    plotIndices = new int[nInputs];
    nPlots = 0;
    for (ii = 0; ii < nInputs; ii++)
      if (rsIndices == NULL || rsIndices[ii] >= 0)
        if (designParams == NULL || designParams[ii] == 0)
          plotIndices[nPlots++] = ii;
  }
  printOutTS(PL_INFO,
       "MCMC Plot summary: input number to be plotted are (%d):\n",nPlots);
  for (ii = 0; ii < nPlots; ii++)
    printOutTS(PL_INFO, "   Input %4d\n", plotIndices[ii]+1);

  // option to add discrepancy function and a posterior sample
  // ==> modelFormFlag, genPosteriors
  if (psAnaExpertMode_ == 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,"*** OPTION TO ADD A DISCREPANCY FUNCTION:\n\n");
    printOutTS(PL_INFO,
         "To use this feature, first make sure that the observation\n");
    printOutTS(PL_INFO,
         "data file specified earlier has design parameters specified\n");
    printOutTS(PL_INFO,
         "since the discrepancy function is to be a function of these\n");
    printOutTS(PL_INFO,
         "design parameters (if not, a constant discrepancy function\n");
    printOutTS(PL_INFO,"is to be created).\n");
    printOutTS(PL_INFO,
         "NOTE: if you don't know what this is, just say NO.\n");
    printf("===> Add discrepancy function ? (y or n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') modelFormFlag = 1;
    if (modelFormFlag == 1 && dnInputs == 0)
    {
      printOutTS(PL_INFO,
       "NOTE: No design inputs ==> discrepancy will be a constant function.\n");
    }
    if (modelFormFlag == 1 && dnSamples == 1)
    {
      printOutTS(PL_INFO,
       "NOTE: 1 experiment ==> discrepancy will be a constant function.\n");
    }

    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
       "*** OPTION TO CREATE A SAMPLE FROM THE POSTERIOR DISTRIBUTIONS:\n\n");
    printOutTS(PL_INFO,
       "In addition to generating the posterior distributions, you can\n");
    printOutTS(PL_INFO,
       "also draw a sample from these posteriors. The posterior sample\n");
    printOutTS(PL_INFO,
       "can be used as prior sample for another simulator/emulator.\n");
    printOutTS(PL_INFO,
       "NOTE: if you don't know what this is, just say no.\n");
    printf("==> Create posterior sample for the input parameters? (y/n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') genPosteriors = 1;
  }
  printEquals(PL_INFO, 0);

  if (modelFormFlag == 1)
  {
    int    *ExpSamStates, ind, ExpNSamples, dfaType, dnPerDim=16;
    double *dOneSample, expdata, simdata;
    double *ExpSamOutputs, *ExpSamInputs, *tSamInputs, *settings;

    ExpNSamples   = dnSamples;
    ExpSamInputs  = new double[ExpNSamples*nInputs];
    ExpSamOutputs = new double[ExpNSamples];
    ExpSamStates  = new int[ExpNSamples];
    discFuncConstantMeans = new double[nOutputs];
    discFuncConstantStds  = new double[nOutputs];
    for (ii2 = 0; ii2 < nOutputs; ii2++)
      discFuncConstantMeans[ii2] = discFuncConstantStds[ii2] = 
                                      PSUADE_UNDEFINED;

    printOutTS(PL_INFO,
         "*** SELECT RESPONSE SURFACE TYPE FOR DISCREPANCY FUNCTION:\n");
    dfaType = -1;
    while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
    {
      writeFAInfo(-1);
      sprintf(charString, "===> Enter your choice : ");
      dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
    }

    settings = new double[nInputs];
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      settings[ii2] = 0.5*(lower[ii2] + upper[ii2]);
      if (rsIndices != NULL && rsIndices[ii2] < 0)
        settings[ii2] = rsValues[ii2];
    }

    faPtrs1 = new FuncApprox*[nOutputs];
    dOneSample = new double[nInputs];
    tSamInputs = new double[ExpNSamples*nInputs];
    int        askFlag = 0, *states=NULL;
    double     *tLowers = new double[nInputs];
    double     *tUppers = new double[nInputs];
    PsuadeData *dPtr = new PsuadeData();
    char       **iNames;
    for (ii = 0; ii < nOutputs; ii++)
    {
      for (kk = 0; kk < dnSamples; kk++)
      {
        cnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          if (designParams != NULL && designParams[ii2] == 1)
          {
            dOneSample[ii2] = dSamInputs[kk*dnInputs+cnt];
            cnt++;
          }
          else dOneSample[ii2] = settings[ii2];
        }

        simdata = 0.0;
        if (psAnaExpertMode_ == 1 && askFlag == 0)
        {
          printOutTS(PL_INFO,
               "To create discrepancy functions, the calibration\n");
          printOutTS(PL_INFO,
               "parameters need to be set to some nominal values.\n");
          printOutTS(PL_INFO,
               "You can choose the nominal values, or it will be\n");
          printOutTS(PL_INFO,
               "set to the input means or mid points of the ranges.\n");
          printf( "Set nomininal values yourself ? (y or n) ");
          scanf("%s", charString);
          fgets(lineIn,1000,stdin);
          if (charString[0] == 'y')
          {
            for (ii2 = 0; ii2 < nInputs; ii2++)
            {
              if ((rsIndices == NULL || rsIndices[ii2] >= 0) &&
                  (designParams == NULL || designParams[ii2] == 0))
              {
                printOutTS(PL_INFO,
                    "Input %d has lower and upper bounds = %e %e\n",
                    ii2+1, lower[ii2], upper[ii2]);
                sprintf(charString, 
                        "Nominal value for input %d : ",ii2+1);
                dOneSample[ii2] = getDouble(charString);
                settings[ii2]   = dOneSample[ii2];
              }
            }
          }
          askFlag = 1;
        }
        simdata = faPtrs[ii]->evaluatePoint(dOneSample);
        expdata = dSamMeans[kk*nOutputs+ii];

        if (printLevel >= 4)
        {
          printOutTS(PL_INFO, 
               "Experiment %4d (out of %d) : ",kk+1,dnSamples);
          for (ii2 = 0; ii2 < nInputs; ii2++)
             printOutTS(PL_INFO,
                  "Input %7d = %12.4e ",ii2+1,dOneSample[ii2]);
          printOutTS(PL_INFO, 
               "simuation, experimental data = %12.4e %12.4e\n",
                     simdata, expdata);
        }
        ExpSamOutputs[kk] = expdata - simdata;
      }

      if (dnInputs > 0)
      {
        for (kk = 0; kk < ExpNSamples*dnInputs; kk++)
          tSamInputs[kk] = dSamInputs[kk];
        iNames = new char*[dnInputs];
        cnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          if (designParams[ii2] == 1)
          {
            iNames[cnt] = new char[100];
            tLowers[cnt] = lower[ii2];
            tUppers[cnt] = upper[ii2];
            if (qData.strArray_ == NULL)
                 sprintf(iNames[cnt], "X%d", ii2+1);
            else strcpy(iNames[cnt], qData.strArray_[ii2]);
            cnt++;
          }
        }
        dPtr->updateInputSection(ExpNSamples,dnInputs,NULL,tLowers,
                    tUppers,dSamInputs, iNames, NULL,NULL,NULL,NULL);
        for (ii2 = 0; ii2 < dnInputs; ii2++) delete [] iNames[ii2];
        delete [] iNames;
      }
      else
      {
        iNames = new char*[1];
        iNames[0] = new char[100];
        sprintf(iNames[0], "X0");
        for (ii2 = 0; ii2 < ExpNSamples; ii2++) tSamInputs[ii2] = 0.5;
        tLowers[0] = 0.0;
        tUppers[0] = 1.0;
        dPtr->updateInputSection(ExpNSamples,iOne,NULL,tLowers,tUppers,
                         tSamInputs, iNames, NULL,NULL,NULL,NULL);
        delete [] iNames[0];
        delete [] iNames;
      }
      states = new int[ExpNSamples];
      for (kk = 0; kk < ExpNSamples; kk++) states[kk] = 1;
      iNames = new char*[1];
      iNames[0] = new char[100];
      sprintf(iNames[0], "Y%d", ii+1);
      dPtr->updateOutputSection(ExpNSamples,iOne,ExpSamOutputs,
                                   states,iNames);
      delete [] states;
      delete [] iNames[0];
      delete [] iNames;
      dPtr->updateMethodSection(PSUADE_SAMP_MC, ExpNSamples, 1, -1, -1);
      sprintf(charString, "psDiscrepancyModel%d", ii+1);
      dPtr->writePsuadeFile(charString, 0);

      printOutTS(PL_INFO,
           "Creating discrepancy response surface for output %d\n",ii+1);
      faPtrs1[ii] = NULL;
      if (dnInputs > 0 && dnSamples > 1)
      {
        faPtrs1[ii] = genFA(dfaType,dnInputs,iOne,ExpNSamples);
        if (faPtrs1[ii] == NULL)
        {
          printOutTS(PL_ERROR,
             "MCMC ERROR: cannot create discrepancy func for output %d.\n",
               ii+1);
          return -1.0;
        }
      }
      if (faPtrs1[ii] != NULL)
      {
        faPtrs1[ii]->setNPtsPerDim(dnPerDim);
        faPtrs1[ii]->setBounds(lower, upper);
        faPtrs1[ii]->setOutputLevel(0);
        faPtrs1[ii]->initialize(tSamInputs,ExpSamOutputs);
      }
      else
      {
        discFuncConstantMeans[ii] = 0.0;
        for (kk = 0; kk < ExpNSamples; kk++)
          discFuncConstantMeans[ii] += ExpSamOutputs[kk];
        discFuncConstantMeans[ii] /= (double) ExpNSamples;
        discFuncConstantStds[ii] = 0.0;
        for (kk = 0; kk < ExpNSamples; kk++)
          discFuncConstantStds[ii] +=
                pow(ExpSamOutputs[kk]-discFuncConstantMeans[ii],2.0);
        discFuncConstantStds[ii] = 
                 sqrt(discFuncConstantStds[ii]/ExpNSamples);
      }
      printOutTS(PL_INFO,
           "Discrepancy response surface for output %d created.\n",ii+1);
    }
    delete [] ExpSamInputs;
    delete [] dOneSample;
    delete [] ExpSamOutputs;
    delete [] settings;
    delete [] tSamInputs;
    delete [] tLowers;
    delete [] tUppers;
    delete dPtr;
  }

  //    set up constraint filters, if any
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"MCMC INFO: creating constraints, if there is any.\n");
  printOutTS(PL_INFO,
       "     Constraints remove infeasible regions from the priors.\n");
  printOutTS(PL_INFO,
       "     Constraints can be specified by RS constraint files.\n");
  constrPtr = new RSConstraints();
  constrPtr->genConstraints(dataPtr);
  printEquals(PL_INFO, 0);

  //  set up for inference
  double *XRange  = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) XRange[ii] = 1.0 / (upper[ii]-lower[ii]);
  double *Xmax = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) Xmax[ii] = 0;
  mostLikelyInput_ = new double[nInputs_];
  means_  = new double[nInputs_];
  sigmas_ = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) mostLikelyInput_[ii] = 0;
  mostLikelyOutput_ = new double[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++) mostLikelyOutput_[ii] = 0;

  int **pbins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    pbins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) pbins[ii][jj] = 0;
  }
  int ****pbins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    pbins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      pbins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        pbins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          pbins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }

  int **bins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    bins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) bins[ii][jj] = 0;
  }
  int ****bins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    bins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      bins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        bins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          bins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }
  double **dbins = new double*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    dbins[ii] = new double[nInputs];
    for (jj = 0; jj < nInputs; jj++) dbins[ii][jj] = 0;
  }
  double ****dbins2 = new double***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    dbins2[jj] = new double**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      dbins2[jj][jj2] = new double*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        dbins2[jj][jj2][ii] = new double[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          dbins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }
  int samInc = 10000;
  if (samInc * dnSamples > 100000) samInc = 100000 / dnSamples; 
  maxSamples = maxSamples / samInc;
  maxSamples = maxSamples * samInc;

  int      methodSave, *SSS;
  double   *XXX, *YYY;
  Sampling *sampler;
  psVector vecLB, vecUB, vecOut; 
  if (noPDF == 1)
  {
    printf("MCMC_BF INFO: no PDF, use uniform for priors.\n");
    sampler = SamplingCreateFromID(PSUADE_SAMP_MC);
    sampler->setInputBounds(nInputs, lower, upper);
    sampler->setOutputParams(1);
    sampler->setSamplingParams(maxSamples, 1, 1);
    sampler->initialize(0);
    vecOut.setLength(maxSamples*nInputs);
    SSS = new int[maxSamples];
    XXX = vecOut.getDVector();
    YYY = new double[maxSamples];
    sampler->getSamples(maxSamples, nInputs, 1, XXX, YYY, SSS);
    delete [] SSS;
    delete [] YYY;
    delete sampler;
  }
  else
  {
    printf("MCMC_BF INFO: has PDF, draw sample from distribution.\n");
    dataPtr->getParameter("method_sampling", pPtr);
    methodSave = pPtr.intData_;
    dataPtr->updateMethodSection(PSUADE_SAMP_MC,-1,-1,-1,-1);
    PDFManager *pdfman = new PDFManager();
    pdfman->initialize(dataPtr);
    vecLB.load(nInputs, lower);
    vecUB.load(nInputs, upper);
    vecOut.setLength(maxSamples*nInputs);
    pdfman->genSample(maxSamples, vecOut, vecLB, vecUB);
    dataPtr->updateMethodSection(methodSave,-1,-1,-1,-1);
    delete pdfman;
  }

  //  perform inference
  psVector vecXSam,vecYSam,vecXDes,vecYDes,vecYSamStd,vecYDesStd;
  psVector vecInferOut, vecTempOut, vecTmpMeans, vecTmpStds;
  vecXSam.setLength(samInc*dnSamples*nInputs);
  vecYSam.setLength(samInc*dnSamples*nOutputs);
  vecYSamStd.setLength(samInc*dnSamples*nOutputs);
  vecXDes.setLength(samInc*dnSamples*nInputs);
  vecYDes.setLength(samInc*dnSamples*nOutputs);
  vecYDesStd.setLength(samInc*dnSamples*nOutputs);
  vecInferOut.setLength(maxSamples);
  vecTempOut.setLength(maxSamples);
  double *XSample = vecXSam.getDVector();
  double *YSample = vecYSam.getDVector();
  double *YSamStd = vecYSamStd.getDVector();
  double *XDesign = vecXDes.getDVector();
  double *YDesign = vecYDes.getDVector();
  double *YDesStd = vecYDesStd.getDVector();
  double *inferenceSamOut = vecInferOut.getDVector();
  double *inferenceTmpOut = vecTempOut.getDVector();
  double *inferenceSamIns = vecOut.getDVector();
  double stdev, stdv2, YT, YT1, YT2;
  int    dcnt, ss = 0, index, nFail, printStep=0, twiceFlag=0;
  int    nInpsActive = 0, passCnt;
  double Ymin = PSUADE_UNDEFINED, dd;
  double Ymax = -PSUADE_UNDEFINED;
  vecTmpMeans.setLength(100*nInputs);
  vecTmpStds.setLength(100*nInputs);
  double *tmpMeans = vecTmpMeans.getDVector();
  double *tmpStds  = vecTmpStds.getDVector();

  fp = NULL;
  for (ii2 = 0; ii2 < nInputs; ii2++)
  {
    if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii2] >=0)) &&
        (designParams == NULL || designParams[ii2] == 0)) nInpsActive++;
  }

  printOutTS(PL_INFO, "MCMC_BF Inference begins ... \n");
  fflush(stdout);
  while (ss < maxSamples)
  {
    cnt = (ss+1) / (maxSamples/20);
    if (cnt != printStep)
    {
      printOutTS(PL_INFO, "%3.0f%% ",5.0*(printStep+1) );
      fflush(stdout);
    }

    for (jj = 0; jj < samInc; jj++)
    {
      index = jj * dnSamples;
      for (kk2 = 0; kk2 < dnSamples; kk2++)
      {
        dcnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          XSample[(index+kk2)*nInputs+ii2] = 
                    inferenceSamIns[(ss+jj)*nInputs+ii2]; 
          if (designParams != NULL && designParams[ii2] == 1)
          {
             XSample[(index+kk2)*nInputs+ii2] =
                                   dSamInputs[kk2*dnInputs+dcnt];
             XDesign[(index+kk2)*dnInputs+dcnt] =
                                   dSamInputs[kk2*dnInputs+dcnt];
             dcnt++;
          }
          if (rsIndices != NULL && rsIndices[ii2] < 0) 
            XSample[(index+kk2)*nInputs+ii2] = rsValues[ii2];
        }
      }
      for (ii2 = 0; ii2 < dnSamples*nOutputs; ii2++)
      {
        YDesign[index*nOutputs+ii2] =
                YDesStd[index*nOutputs+ii2] = 0.0;
        YSample[index*nOutputs+ii2] =
                YSamStd[index*nOutputs+ii2] = 0.0;
      }
    }
    for (ii2 = 0; ii2 < nOutputs; ii2++)
    {
      if (rsErrFlag == 1)
      {
        faPtrs[ii2]->evaluatePointFuzzy(samInc*dnSamples, 
                        XSample,&YSample[ii2*dnSamples*samInc],
                        &YSamStd[ii2*dnSamples*samInc]);
        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
        {
          faPtrs1[ii2]->evaluatePointFuzzy(samInc*dnSamples,
                             XDesign,&YDesign[ii2*dnSamples*samInc],
                             &YDesStd[ii2*dnSamples*samInc]);
        }
        else if (discFuncConstantMeans != NULL &&
                 discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < dnSamples*samInc; kk2++)
          {
            YDesign[ii2*dnSamples*samInc+kk2] =
                             discFuncConstantMeans[ii2];
            YDesStd[ii2*dnSamples*samInc+kk2] = 0.0;
          }
        }
      }
      else
      {
        faPtrs[ii2]->evaluatePoint(samInc*dnSamples,XSample,
                                   &YSample[ii2*dnSamples*samInc]);
        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
        {
          faPtrs1[ii2]->evaluatePoint(dnSamples*samInc,XDesign,
                                    &YDesign[ii2*dnSamples*samInc]);
        }
        else if (discFuncConstantMeans != NULL &&
                 discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < dnSamples*samInc; kk2++)
                  YDesign[ii2*dnSamples*samInc+kk2] =
                                       discFuncConstantMeans[ii2];
        }
        for (kk2 = 0; kk2 < dnSamples*samInc; kk2++)
          YSamStd[ii2*dnSamples*samInc+kk2] =
                       YDesStd[ii2*dnSamples*samInc+kk2] = 0.0;
      }
    }
    for (jj = 0; jj < samInc; jj++)
    {
      inferenceSamOut[ss+jj] = 0.0;
      index = jj * dnSamples;
      for (ii2 = 0; ii2 < nOutputs; ii2++)
      {
        for (kk2 = 0; kk2 < dnSamples; kk2++)
        {
          YT1 = YSample[ii2*dnSamples*samInc+index+kk2] +
                YDesign[ii2*dnSamples*samInc+index+kk2];
          stdev = YSamStd[ii2*dnSamples*samInc+index+kk2];
          stdv2 = YDesStd[ii2*dnSamples*samInc+index+kk2];
          YT2 = pow((YT1-dSamMeans[kk2*nOutputs+ii2]),2.0) /
               (pow(dSamStdevs[kk2*nOutputs+ii2],2.0) +
                stdev*stdev + stdv2*stdv2);
          inferenceSamOut[ss+jj] += YT2;
        }
      }
      //*/ inferenceSamOut[ss+jj] /= (dnSamples*nOutputs);
      if (combFlag == 1)
        inferenceSamOut[ss+jj] /= (double) dnSamples;
    }
    nFail = 0;
    for (jj = 0; jj < samInc; jj++)
    {
      YT1 = constrPtr->evaluate(&inferenceSamIns[(ss+jj)*nInputs],
                                inferenceSamOut[ss+jj], status);
      if (status == 0)
      {
        inferenceSamOut[ss+jj] = 0.0;
        nFail++;
      }
    }
    ss += samInc;

    if (cnt != printStep)
    {
      printStep++;
      for (jj = 0; jj < ss; jj++) 
        inferenceTmpOut[jj] = inferenceSamOut[jj];
      Ymin = PSUADE_UNDEFINED;
      for (jj = 0; jj < ss; jj++) 
      {
        dd = inferenceTmpOut[jj];
        if (dd < Ymin) Ymin = dd;
        if (dd > Ymax) Ymax = dd;
      }
      for (jj = 0; jj < ss; jj++) 
        inferenceTmpOut[jj] = exp(-0.5*(inferenceTmpOut[jj]-Ymin));

      printOutTS(PL_INFO,"\n");
      passCnt = 0;
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
        if ((rsIndices == NULL || 
            (rsIndices != NULL && rsIndices[ii2] >=0)) &&
           (designParams == NULL || designParams[ii2] == 0))
        {
          YT = 0.0;
          for (jj = 0; jj < ss; jj++) YT += inferenceTmpOut[jj];
          YT1 = 0.0;
          for (jj = 0; jj < ss; jj++)
            YT1 += inferenceSamIns[jj*nInputs+ii2] * 
                   inferenceTmpOut[jj] / YT; 
          printOutTS(PL_INFO,"MCMC_BF: input %3d mean    = %e\n", 
                     ii2+1, YT1);
          YT2 = 0.0;
          for (jj = 0; jj < ss; jj++)
            YT2 += pow(inferenceSamIns[jj*nInputs+ii2]-YT1,2.0)*
                       inferenceTmpOut[jj]/YT; 
          printOutTS(PL_INFO,"MCMC: input %3d std dev = %e\n",
                     ii2+1,sqrt(YT2));
          tmpMeans[(printStep-1)+ii2*100] = YT1; 
          tmpStds[(printStep-1)+ii2*100] = sqrt(YT2); 
          if (printStep > 2) 
          {
            jj = checkConvergence(3,&tmpMeans[ii2*100+printStep-3],
                         &tmpStds[ii2*100+printStep-3],ss-samInc);
            if (jj == 1)
            {
              passCnt++;
              printf("MCMC input %3d converged.\n",ii2+1);
            }
          }
        }
      }
      if (passCnt == nInpsActive)
      {
        if      (twiceFlag <= 0) twiceFlag++;
        else if (twiceFlag == 1) maxSamples = ss;
      }
    }
  }

  Ymin = PSUADE_UNDEFINED;
  Ymax = -PSUADE_UNDEFINED;
  for (ss = 0; ss < maxSamples; ss++)
  {
    dd = inferenceSamOut[ss];
    if (dd < Ymin) Ymin = dd;
    if (dd > Ymax) Ymax = dd;
  }
  for (ss = 0; ss < maxSamples; ss++)
    inferenceSamOut[ss] = exp(-0.5*(inferenceSamOut[ss]-Ymin));

  index = -1;
  Ymax = -PSUADE_UNDEFINED;
  for (ss = 0; ss < maxSamples; ss++)
  {
    dd = inferenceSamOut[ss] * exp(-0.5*Ymin);
    if (dd > Ymax)
    {
      Ymax = dd;
      index = ss;
    }
  } 
  if (index >= 0)
  {
    for (ii = 0; ii < nInputs; ii++) 
      Xmax[ii] = inferenceSamIns[index*nInputs+ii];
  }
  printf("Maximum likelihood estimated solution:\n");
  for (ii = 0; ii < nInputs; ii++) 
    printf("Input %3d = %16.8e\n",ii+1,Xmax[ii]);
  printf("Negative log likelihood (unnormalized) = %e\n", -log(Ymax));

  if (Ymax > 0)
  {
    Ymax = Ymax / exp(-0.5*Ymin);
    for (ss = 0; ss < maxSamples; ss++) inferenceSamOut[ss] /= Ymax;
  }

  int    ii3, index2;
  for (ss = 0; ss < maxSamples; ss++)
  {
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      YT1 = (inferenceSamIns[ss*nInputs+ii2] - lower[ii2]) * XRange[ii2];
      index = (int) (YT1 * nbins);
      if (index > nbins)
        printOutTS(PL_ERROR,"MCMC binning error 1 in file %s, line %d.\n",
                   __FILE__, __LINE__);
      if (index < 0)
        printOutTS(PL_ERROR,"MCMC binning error 2 in file %s, line %d.\n",
                   __FILE__, __LINE__);
      if (index >= nbins) index = nbins - 1;
      if (index <  0)     index = 0;
      dbins[index][ii2] += inferenceSamOut[ss];
      pbins[index][ii2]++;
    }
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      YT1 = (inferenceSamIns[ss*nInputs+ii2] - lower[ii2]) * XRange[ii2];
      index = (int) (YT1 * nbins);
      if (index >= nbins) index = nbins - 1;
      if (index < 0)      index = 0;
      for (ii3 = 0; ii3 < nInputs; ii3++)
      {
        YT2 = (inferenceSamIns[ss*nInputs+ii3] - lower[ii3]) * XRange[ii3];
        index2 = (int) (YT2 * nbins);
        if (index2 >= nbins) index2 = nbins - 1;
        if (index2 < 0)      index2 = 0;
        dbins2[index][index2][ii2][ii3] += inferenceSamOut[ss];
        pbins2[index][index2][ii2][ii3]++;
      }
    }
  }

  Ymax = 0;
  for (kk = 0; kk < nbins; kk++)
    for (ii2 = 0; ii2 < nInputs; ii2++) 
      if (dbins[kk][ii2] > Ymax) Ymax = dbins[kk][ii2];
  for (kk = 0; kk < nbins; kk++)
    for (ii2 = 0; ii2 < nInputs; ii2++) 
      dbins[kk][ii2] = dbins[kk][ii2] / Ymax * maxSamples;
  for (kk = 0; kk < nbins; kk++)
    for (ii2 = 0; ii2 < nInputs; ii2++) bins[kk][ii2] = (int) dbins[kk][ii2];
  Ymax = 0;
  for (jj = 0; jj < nbins; jj++)
    for (kk = 0; kk < nbins; kk++)
      for (ii2 = 0; ii2 < nInputs; ii2++)
        for (ii3 = 0; ii3 < nInputs; ii3++)
          if (dbins2[jj][kk][ii2][ii3] > Ymax) 
            Ymax = dbins2[jj][kk][ii2][ii3];
  for (jj = 0; jj < nbins; jj++)
    for (kk = 0; kk < nbins; kk++)
      for (ii2 = 0; ii2 < nInputs; ii2++)
        for (ii3 = 0; ii3 < nInputs; ii3++)
          dbins2[jj][kk][ii2][ii3] = dbins2[jj][kk][ii2][ii3] / Ymax *
                                     maxSamples;
  for (jj = 0; jj < nbins; jj++)
    for (kk = 0; kk < nbins; kk++)
      for (ii2 = 0; ii2 < nInputs; ii2++)
        for (ii3 = 0; ii3 < nInputs; ii3++)
          bins2[jj][kk][ii2][ii3] = (int) dbins2[jj][kk][ii2][ii3];
     
  double psum;
  for (ii = 0; ii < nInputs; ii++)
  {
    if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) &&
        (designParams == NULL || designParams[ii] == 0))
    {
      printOutTS(PL_INFO,
                 "MCMC: input %3d value at peak of likelihood = %e\n",
                 ii+1, Xmax[ii]);
      psum = 0.0;
      for (ss = 0; ss < maxSamples; ss++) psum += inferenceSamOut[ss];
      YT = 0.0;
      for (ss = 0; ss < maxSamples; ss++)
        YT += inferenceSamIns[ss*nInputs+ii] * inferenceSamOut[ss] / psum; 
      means_[ii] = YT;
      printOutTS(PL_INFO,"MCMC: input %3d mean    = %e\n", ii+1, YT);
      YT = 0.0;
      for (ss = 0; ss < maxSamples; ss++)
        YT += pow(inferenceSamIns[ss*nInputs+ii]-means_[ii],2.0)*
                  inferenceSamOut[ss]/psum; 
      sigmas_[ii] = sqrt(YT);
      printOutTS(PL_INFO,"MCMC: input %3d std dev = %e\n",ii+1,sigmas_[ii]);
      mostLikelyInput_[ii] = Xmax[ii];
    }
  }
 
  for (ii = 0; ii < nInputs; ii++) XRange[ii] = 1.0 / XRange[ii];
  genMatlabFile(nInputs,lower,upper,XRange,nPlots,plotIndices,nbins,
        pbins,pbins2,bins,bins2,qData,0,0,NULL,NULL,Xmax,Ymin);

  if (genPosteriors == 1)
  {
    int    maxPostSam=50000;
    double dmax = 0.0;
    for (ss = 0; ss < maxSamples; ss++) 
      if (inferenceSamOut[ss] > dmax) dmax = inferenceSamOut[ss];
    if (dmax == 0)
    {
      printOutTS(PL_ERROR,
           "MCMC_BF: ERROR encountered in posterior sample generation.\n");
      fp = NULL;
    }
    else 
    {
      fp = fopen("MCMCPostSample", "w");
      for (ss = 0; ss < maxSamples; ss++) 
        inferenceSamOut[ss] = inferenceSamOut[ss] / dmax;
    }
    if (fp != NULL)
    {
      fprintf(fp, "PSUADE_BEGIN\n");
      fprintf(fp, "%d %d\n", maxPostSam, nInputs);
      if (qData.strArray_ != NULL)
      {
        fprintf(fp, "# ");
        for (jj = 0; jj < nInputs; jj++)
           fprintf(fp,"%s ", qData.strArray_[jj]);
        fprintf(fp, "\n");
      }
      int    count=0, ichoose, nValid=0;
      double dchoose;
      psIVector vecValid;
      vecValid.setLength(maxSamples);
      for (ss = 0; ss < maxSamples; ss++) 
        if (inferenceSamOut[ss] > 1e-5) vecValid[nValid++] = ss;
      while (count < maxPostSam)
      {
        ichoose = PSUADE_rand() % maxSamples;
        ichoose = vecValid[ichoose];
        dchoose = PSUADE_drand();
        if (inferenceSamOut[ichoose] > dchoose)
        {
          fprintf(fp, "%d ", count+1);
          for (jj = 0; jj < nInputs; jj++)
          {
             if ((rsIndices == NULL || rsIndices[jj] >= 0) &&
                 (designParams == NULL || designParams[jj] == 0))
             {
                fprintf(fp, "%e ", inferenceSamIns[ichoose*nInputs+jj]);
             }
             else if (rsIndices != NULL && rsIndices[jj] < 0)
                fprintf(fp, "%e ", rsValues[jj]);
             else if (designParams != NULL && designParams[jj] != 0)
                fprintf(fp, "%e ", 0.5 * (upper[jj] + lower[jj]));
          }
          fprintf(fp, "\n");
          count++;
        }
      }
      fprintf(fp, "PSUADE_END\n");
      vecXSam.setLength(nInputs*dnSamples);
      vecYSam.setLength(nOutputs*dnSamples);
      vecXDes.setLength(dnInputs*dnSamples);
      vecYDes.setLength(nOutputs*dnSamples);
      YSample = vecYSam.getDVector();
      YDesign = vecYDes.getDVector();
      fprintf(fp, "Optimal parameter values: \n");
      dcnt = 0;
      for (ii2 = 0; ii2 < nInputs; ii2++) 
      {
        if (designParams != NULL && designParams[ii2] == 1) 
        {
          fprintf(fp, "%16.8e\n", dSamInputs[dcnt]);
          dcnt++;
        }
        else if (rsIndices != NULL && rsIndices[ii2] < 0)
          fprintf(fp, "%16.8e\n", rsValues[ii2]);
        else
          fprintf(fp, "%16.8e\n", Xmax[ii2]);
      }
      fprintf(fp,
        "MLE Statistics (prediction, exp data, sd, -loglikelihood)\n");
      for (kk2 = 0; kk2 < dnSamples; kk2++)
      {
        dcnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++) 
        {
          vecXSam[kk2*nInputs+ii2] = Xmax[ii2];
          if (designParams != NULL && designParams[ii2] == 1)
          {
            vecXSam[kk2*nInputs+ii2] = dSamInputs[kk2*dnInputs+dcnt];
            vecXDes[kk2*dnInputs+dcnt] = dSamInputs[kk2*dnInputs+dcnt];
            dcnt++;
          }
          if (rsIndices != NULL && rsIndices[ii2] < 0)
            vecXSam[kk2*nInputs+ii2] = rsValues[ii2];
        }
      }
      for (ii2 = 0; ii2 < dnSamples*nOutputs; ii2++)
        YDesign[ii2] = YSample[ii2] = 0.0;
      for (ii2 = 0; ii2 < nOutputs; ii2++) 
      {
        faPtrs[ii2]->evaluatePoint(dnSamples,vecXSam.getDVector(),
                                   &YSample[ii2*dnSamples]);
        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
        {
          faPtrs1[ii2]->evaluatePoint(dnSamples,vecXDes.getDVector(),
                                      &YDesign[ii2*dnSamples]);
        }
        else if (discFuncConstantMeans != NULL &&
                 discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < dnSamples; kk2++)
            YDesign[ii2*dnSamples+kk2] = discFuncConstantMeans[ii2];
        }
      }
      for (kk2 = 0; kk2 < dnSamples; kk2++)
      {
        for (ii2 = 0; ii2 < nOutputs; ii2++)
        {
          YT1 = YSample[ii2*dnSamples+kk2] + YDesign[ii2*dnSamples+kk2];
          YT2 = dSamMeans[kk2*nOutputs+ii2];
          stdev = dSamStdevs[kk2*nOutputs+ii2];
          stdv2 = 0.5 * pow(YT1 - YT2, 2.0) / (stdev * stdev);
          fprintf(fp,"%4d: %16.8e %16.8e %16.8e %16.8e\n",kk2+1,
                  YT1,YT2,stdev,stdv2);
        }
      }
      fclose(fp);
    }
    printOutTS(PL_INFO,
         "MCMC: 'MCMCPostSample' file has a posterior sample.\n");
  }

  int    nInps, nOuts, nSams, *states;
  double *allOuts;
  char   **oNames;
  PsuadeData *filePtr1, *filePtr2;
  if (modelFormFlag == 1)
  {
    sprintf(charString, "psDiscrepancyModel1");
    filePtr1 = new PsuadeData();
    status = filePtr1->readPsuadeFile(charString);
    if (status != 0)
    {
       printOutTS(PL_ERROR,
            "MCMC ERROR: cannot read file %s in PSUADE format.\n",
            charString);
       exit(1);
    }
  }
  if (modelFormFlag == 1 && status == 0)
  {
    filePtr1->getParameter("input_ninputs", pPtr);
    nInps = pPtr.intData_;
    filePtr1->getParameter("output_noutputs", pPtr);
    nOuts = pPtr.intData_;
    filePtr1->getParameter("method_nsamples", pPtr);
    nSams = pPtr.intData_;
    filePtr1->getParameter("output_sample", pOutputs);
    unlink(charString);
    allOuts = new double[nOutputs * nSams];
    for (jj = 0; jj < nSams; jj++)
      allOuts[jj*nOutputs] = pOutputs.dbleArray_[jj];
    pOutputs.clean();
    oNames = new char*[nOutputs];
    oNames[0] = new char[100];
    sprintf(oNames[0], "Y1");
    for (ii = 1; ii < nOutputs; ii++)
    {
      filePtr2 = new PsuadeData();
      sprintf(charString, "psDiscrepancyModel%d", ii+1);
      status = filePtr2->readPsuadeFile(charString);
      if (status != 0) break;
      filePtr2->getParameter("input_ninputs", pPtr);
      if (pPtr.intData_ != nInps) break;
      filePtr2->getParameter("output_noutputs", pPtr);
      if (pPtr.intData_ != nOuts) break;
      filePtr2->getParameter("method_nsamples", pPtr);
      if (pPtr.intData_ != nSams) break;
      filePtr2->getParameter("output_sample", pOutputs);
      delete filePtr2;
      unlink(charString);
      for (jj = 0; jj < nSams; jj++)
         allOuts[jj*nOutputs+ii] = pOutputs.dbleArray_[jj];
      pOutputs.clean();
      oNames[ii] = new char[100];
      sprintf(oNames[ii], "Y%d", ii+1);
    }
    if (nOutputs == 1)
    {
       sprintf(charString, "psDiscrepancyModel");
       filePtr1->writePsuadeFile(charString, 0);
    }
    else if (ii == nOutputs)
    {
      states = new int[nSams];
      for (jj = 0; jj < nSams; jj++) states[jj] = 1;
      filePtr1->updateOutputSection(nSams,nOutputs,allOuts,states,oNames);
      sprintf(charString, "psDiscrepancyModel");
      filePtr1->writePsuadeFile(charString, 0);
      printOutTS(PL_INFO,
           "MCMC INFO: a sample (inputs/outputs) the discrepancy model\n");
      printOutTS(PL_INFO,"           is now in psDiscrepancyModel.\n");
      delete [] states;
      for (ii = 1; ii < nOutputs; ii++) delete [] oNames[ii];
    }
    else
    {
      printOutTS(PL_INFO,
           "MCMC INFO: unsuccessful creation of discrepancy sample file\n");
    }
    delete [] oNames[0];
    delete [] oNames;
    delete filePtr1;
    delete [] allOuts;
  }

  // clean up
  delete [] Xmax;
  delete [] XRange;
  if (discFuncConstantMeans != NULL) delete [] discFuncConstantMeans;
  if (discFuncConstantStds  != NULL) delete [] discFuncConstantStds;
  if (plotIndices != NULL) delete [] plotIndices;
  if (designParams != NULL) delete [] designParams;
  if (dSamMeans != NULL) delete [] dSamMeans;
  if (dSamStdevs != NULL) delete [] dSamStdevs;
  if (faPtrs != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++)
      if (faPtrs[ii] != NULL) delete faPtrs[ii];
    delete [] faPtrs;
  }
  if (faPtrs1 != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++)
      if (faPtrs1[ii] != NULL) delete faPtrs1[ii];
    delete [] faPtrs1;
  }
  for (ii = 0; ii < nbins; ii++) 
  {
    delete [] bins[ii];
    delete [] dbins[ii];
    delete [] pbins[ii];
  }
  delete [] bins;
  delete [] pbins;
  delete [] dbins;
  for (jj = 0; jj < nbins; jj++)
  {
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      for (ii = 0; ii < nInputs; ii++) delete [] bins2[jj][jj2][ii];
      for (ii = 0; ii < nInputs; ii++) delete [] dbins2[jj][jj2][ii];
      for (ii = 0; ii < nInputs; ii++) delete [] pbins2[jj][jj2][ii];
      delete [] bins2[jj][jj2];
      delete [] pbins2[jj][jj2];
      delete [] dbins2[jj][jj2];
    }
    delete [] bins2[jj];
    delete [] dbins2[jj];
    delete [] pbins2[jj];
  }
  delete [] bins2;
  delete [] dbins2;
  delete [] pbins2;
  if (dSamInputs != NULL) delete [] dSamInputs;
  if (rsIndices != NULL) delete [] rsIndices;
  if (rsValues  != NULL) delete [] rsValues;
  delete constrPtr;
  return 0.0;
}

// ************************************************************************
// perform MCMC-like analysis (brute force with sample uncertain variables)
// ------------------------------------------------------------------------
double MCMCAnalyzer::analyze_bf2(aData &adata)
{
  int    ii, ii2, jj, jj2, kk, kk2, status, cnt, iOne=1, iZero=0;
  int    printLevel, nbins, faType, genPosteriors=0;
  int    maxSamples, modelFormFlag=0, combFlag=0;
  int    nPlots, *plotIndices=NULL;
  int    dnSamples=0, dnInputs=0;
  double dstatus, ddata;
  char   lineIn[1001], charString[1001], *rsFile=NULL;
  FILE   *fp=NULL;
  pData      pPtr, pOutputs;
  FuncApprox **faPtrs=NULL, **faPtrs1=NULL;
  PsuadeData *dataPtr=NULL;
  RSConstraints *constrPtr;

  // clean up
  if (means_) delete[] means_;
  if (sigmas_) delete[] sigmas_;
  if (mostLikelyInput_) delete[] mostLikelyInput_;
  if (mostLikelyOutput_) delete[] mostLikelyOutput_;
  means_ = NULL;
  sigmas_ = NULL;
  mostLikelyInput_ = NULL;
  mostLikelyOutput_ = NULL;

  // display header 
  printLevel  = adata.printLevel_;
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO,"*              Brute Force Optimizer (2a)\n");
  printEquals(PL_INFO, 0);
  displayBanner_bf2(printLevel);

  // extract data from aData object (passed in from outside)
  int    nInputs, nOutputs, nSamples, *pdfFlags, noPDF; 
  double *XIn, *YIn, *lower, *upper;
  nInputs     = adata.nInputs_;
  nInputs_    = nInputs;
  nOutputs    = adata.nOutputs_;
  nOutputs_   = nOutputs;
  nSamples    = adata.nSamples_;
  XIn         = adata.sampleInputs_;
  YIn         = adata.sampleOutputs_;
  lower       = adata.iLowerB_;
  upper       = adata.iUpperB_;
  dataPtr     = adata.ioPtr_;
  pdfFlags    = adata.inputPDFs_;
  noPDF       = 1;
  if (pdfFlags != NULL)
  {
    for (ii = 0; ii < nInputs; ii++) 
      if (pdfFlags[ii] != 0) noPDF = 0;
  }

  int       UParamsNumToUse=1, imax;
  psIVector rsIndices;
  psVector  rsValues;
  psMatrix  UParamsSample;
  faType = PSUADE_RS_MARS;
  if (dataPtr != NULL)
  {
    dataPtr->getParameter("ana_rstype", pPtr);
    faType = pPtr.intData_;

    status = readIndexFile2(dataPtr, rsIndices, rsValues, UParamsSample);
    if (status < 0) return PSUADE_UNDEFINED;

    if (UParamsSample.nrows() > 0)
    {
      printf("A sample for uncertain parameters has been provided.\n");
      printf("The sample size is %d\n", UParamsSample.nrows());
      imax = UParamsSample.nrows();
      if (imax > 10000) imax = 10000;
      sprintf(charString,
         "Enter the sub-sample size to use for inference (1 - %d): ",imax);
      UParamsNumToUse = getInt(1, imax, charString);

      status = 0;
      cnt = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
        if (rsIndices[ii] >= 1000)
        {
          for (jj = 0; jj < UParamsSample.nrows(); jj++)
          {
            ddata = UParamsSample.getEntry(jj, cnt);
            if ((ddata < lower[ii] || ddata > upper[ii]) && status == 0)
            {
              printf("MCMC_BF2 ERROR: fail uncertain sample bound check.\n");
              printf("   Sample number = %d\n", jj+1);
              printf("   Input         = %d\n", ii+1);
              printf("   UParams Input = %d\n", cnt+1);
              printf("   Sample data   = %e\n", ddata);
              printf("   Input Bounds  = %e %e\n", lower[ii], upper[ii]);
              status = 1;
            }
            if (ddata < lower[ii]) lower[ii] = ddata;
            if (ddata > upper[ii]) upper[ii] = ddata;
          }
          cnt++;
        }
      }
    }
  }
  else if (mode_ == 0)
  {
    printOutTS(PL_INFO,
         "MCMC INFO: since ioPtr=NULL, assume MARS as reponse surface.\n");
  }

  // get experimental data information from the spec file
  psMatrix  dSamInputs, dSamMeans, dSamStdvs; 
  psIVector dParams; 
  dstatus = readSpecFile2(nInputs, nOutputs, dParams, dSamInputs,
                          dSamMeans, dSamStdvs, combFlag, printLevel);
  dnSamples = dSamMeans.nrows();
  dnInputs  = dSamInputs.ncols();

  if (dstatus != 0.0)
  {
    printf("MCMC ERROR: fail to read experimental data file.\n");
    return PSUADE_UNDEFINED;
  }
  for (ii = 0; ii < nInputs; ii++)
  {
    if (rsIndices[ii] < 0 && dParams.length() > 0 && dParams[ii] == 1)
    {
      printOutTS(PL_ERROR,
        "MCMC_BF2 ERROR: inactive input %d cannot be design parameter\n",
        ii+1);
      return PSUADE_UNDEFINED;
    }

    if (rsIndices[ii] == 999 && (dParams.length() > 0 && dParams[ii] == 1))
    {
      printOutTS(PL_ERROR,
        "MCMC_BF2 ERROR: inactive input %d cannot be uncertain parameter\n",
        ii+1);
      return PSUADE_UNDEFINED;
    }
  }

  int rsErrFlag=0;
  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO,
        "*** OPTION TO INCLUDE RESPONSE SURFACE UNCERTAINTIES:\n");
    printOutTS(PL_INFO,
        "\nTo incorporate response surface uncertainties into the\n");
    printOutTS(PL_INFO,
        "likelihood function, make sure stochastic response surfaces\n");
    printOutTS(PL_INFO,
        "are used (GP/Kriging, polynomial regression, or bootstrapped\n");
    printOutTS(PL_INFO,
        "methods). Otherwise, no RS uncertainties will be included.\n");
    printOutTS(PL_INFO,
        "NOTE: if you don't know what this is, just say no below.\n");
    printf( "===> Include response surface uncertainties? (y or n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') rsErrFlag = 1;
    printEquals(PL_INFO, 0);
  }

  printOutTS(PL_INFO,
       "MCMC INFO: CREATING RESPONSE SURFACES FOR ALL OUTPUTS.\n");
  if (nSamples <= 0)
  {
    printOutTS(PL_ERROR, "MCMC ERROR: no sample loaded yet.\n");
    return PSUADE_UNDEFINED;
  }
  else
  {
    faPtrs = new FuncApprox*[nOutputs];
    psVector YY;
    YY.setLength(nSamples);
    for (ii = 0; ii < nOutputs; ii++)
    {
      faType = -1;
      printOutTS(PL_INFO,
           "MCMC INFO: CREATING RESPONSE SURFACE FOR OUTPUT %d.\n",ii+1);
      faPtrs[ii] = genFA(faType, nInputs, iZero, nSamples);
      faPtrs[ii]->setNPtsPerDim(16);
      faPtrs[ii]->setBounds(lower, upper);
      faPtrs[ii]->setOutputLevel(0);
      for (kk = 0; kk < nSamples; kk++) YY[kk] = YIn[kk*nOutputs+ii];

      status = faPtrs[ii]->initialize(XIn, YY.getDVector());
      if (status != 0)
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: Unable to create response surface.\n");
        printOutTS(PL_ERROR,"            Consult PSUADE developers.\n");
        for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
        delete [] faPtrs;
        return PSUADE_UNDEFINED;
      }
    }
  }

  nbins = 20;
  maxSamples = 500000;
  if (nInputs >= 10) maxSamples = 1000000;
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"*** CURRENT DEFAULT PARAMETER SETTINGS : \n\n");
  printOutTS(PL_INFO,"Inference max sample size = %d\n",maxSamples);
  printOutTS(PL_INFO,"Posterior histogram nbins = %d\n",nbins);
  printOutTS(PL_INFO,
       "NOTE: histogram nBins  - resolution of histogram bar graph\n");
  printOutTS(PL_INFO,
       "Turn on ana_expert mode to change these default settings.\n\n");

  if (psAnaExpertMode_ == 1)
  {
    sprintf(charString,
            "Enter maximum inference sample size (500000 - 5000000): ");
    maxSamples = getInt(1000, 50000000, charString);
    if (maxSamples < 500000) maxSamples = 500000;
    if (nInputs >= 10 && maxSamples < 1000000)
      maxSamples = 1000000;
    sprintf(charString,"Enter the number of histogram bins (10 - 25) : ");
    nbins = getInt(10, 50, charString);
  }
  if (psAnaExpertMode_ == 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
        "*** OPTION TO CREATE POSTERIOR PLOTS FOR SELECTED INPUTS:\n\n");
    printOutTS(PL_INFO,
        "MCMC will create MATLAB files for the posterior distributions.\n");
    printOutTS(PL_INFO,
        "You can choose to generate posterior plots for all inputs, or \n");
    printOutTS(PL_INFO,
        "just a selected few (in case there are too many inputs).\n");
    printf("Select inputs for which posterior plots are to be generated.\n");
    sprintf(charString,"Enter input number (-1 for all, 0 to terminate) : ");
    kk = 1;
    plotIndices = new int[nInputs];
    nPlots = 0;
    while (kk != 0 || nPlots < 1)
    {
      kk = getInt(-1, nInputs, charString);
      if (kk == -1)
      {
        for (ii = 0; ii < nInputs; ii++)
        {
          if (rsIndices.length() == 0 || 
              (rsIndices[ii] >= 0 && rsIndices[ii] < 1000))
            if (dParams.length() == 0 || dParams[ii] == 0)
              plotIndices[nPlots++] = ii;
        }
        break;
      }
      if (kk != 0)
      {
        if (rsIndices.length() > 0 && rsIndices[kk-1] < 0)
          printOutTS(PL_ERROR,
              "Input %d has been fixed by the rs index file (no plot).\n",
              kk+1);
        else if (dParams.length() > 0 && dParams[kk-1] == 1)
          printOutTS(PL_ERROR,
              "Input %d is a design parameter (no plot)\n",kk);
        else
          plotIndices[nPlots++] = kk - 1;
      }
      if (kk == 0 && nPlots == 0)
        printOutTS(PL_ERROR,
            "You need to set at least 1 input for plotting posteriors.\n");
    }
    if (nPlots > 1) sortIntList(nPlots, plotIndices);
  }
  else
  {
    plotIndices = new int[nInputs];
    nPlots = 0;
    for (ii = 0; ii < nInputs; ii++)
      if (rsIndices.length() == 0 || 
          (rsIndices[ii] >= 0 && rsIndices[ii] < 1000))
        if (dParams.length() == 0 || dParams[ii] == 0)
          plotIndices[nPlots++] = ii;
  }
  printOutTS(PL_INFO,
       "MCMC Plot summary: input number to be plotted are (%d):\n",nPlots);
  for (ii = 0; ii < nPlots; ii++)
    printOutTS(PL_INFO, "   Input %4d\n", plotIndices[ii]+1);

  // option to add discrepancy function and a posterior sample
  // ==> modelFormFlag, genPosteriors
  if (psAnaExpertMode_ == 1 && UParamsSample.nrows() == 0)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,"*** OPTION TO ADD A DISCREPANCY FUNCTION:\n\n");
    printOutTS(PL_INFO,
         "To use this feature, first make sure that the observation\n");
    printOutTS(PL_INFO,
         "data file specified earlier has design parameters specified\n");
    printOutTS(PL_INFO,
         "since the discrepancy function is to be a function of these\n");
    printOutTS(PL_INFO,
         "design parameters (if not, a constant discrepancy function\n");
    printOutTS(PL_INFO,"is to be created).\n");
    printOutTS(PL_INFO,
         "NOTE: if you don't know what this is, just say NO.\n");
    printf("===> Add discrepancy function ? (y or n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') modelFormFlag = 1;
    if (modelFormFlag == 1 && dnInputs == 0)
    {
      printOutTS(PL_INFO,
       "NOTE: No design inputs ==> discrepancy will be a constant function.\n");
    }
    if (modelFormFlag == 1 && dnSamples == 1)
    {
      printOutTS(PL_INFO,
       "NOTE: 1 experiment ==> discrepancy will be a constant function.\n");
    }
    printEquals(PL_INFO, 0);
  }
  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO,
       "*** OPTION TO CREATE A SAMPLE FROM THE POSTERIOR DISTRIBUTIONS:\n\n");
    printOutTS(PL_INFO,
       "In addition to generating the posterior distributions, you can\n");
    printOutTS(PL_INFO,
       "also draw a sample from these posteriors. The posterior sample\n");
    printOutTS(PL_INFO,
       "can be used as prior sample for another simulator/emulator.\n");
    printOutTS(PL_INFO,
       "NOTE: if you don't know what this is, just say no.\n");
    printf("==> Create posterior sample for the input parameters? (y/n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') genPosteriors = 1;
  }
  printEquals(PL_INFO, 0);

  pData qData; 
  if (dataPtr != NULL) dataPtr->getParameter("input_names", qData);
  psVector discFuncConstantMeans, discFuncConstantStds;
  if (modelFormFlag == 1)
  {
    dstatus = createDiscrepancyFunctions(nInputs, nOutputs, lower, upper, 
                  rsIndices, rsValues, dParams, dnInputs, dnSamples, 
                  dSamInputs, dSamMeans, dataPtr, discFuncConstantMeans,
                  discFuncConstantStds, faPtrs, &faPtrs1, printLevel);
    if (dstatus < 0) return PSUADE_UNDEFINED;
  }

  //    set up constraint filters, if any ==> constrPtr
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"MCMC INFO: creating constraints, if there is any.\n");
  printOutTS(PL_INFO,
       "     Constraints remove infeasible regions from the priors.\n");
  printOutTS(PL_INFO,
       "     Constraints can be specified by RS constraint files.\n");
  constrPtr = new RSConstraints();
  constrPtr->genConstraints(dataPtr);
  printEquals(PL_INFO, 0);

  //  set up for inference
  psVector vecXmax, vecRanges;
  vecRanges.setLength(nInputs);
  double *XRange = vecRanges.getDVector();
  for (ii = 0; ii < nInputs; ii++) XRange[ii]=1.0/(upper[ii]-lower[ii]);
  vecXmax.setLength(nInputs);
  double *Xmax = vecXmax.getDVector();
  mostLikelyInput_ = new double[nInputs_];
  means_  = new double[nInputs_];
  sigmas_ = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) mostLikelyInput_[ii] = 0;
  mostLikelyOutput_ = new double[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++) mostLikelyOutput_[ii] = 0;

  int **pbins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    pbins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) pbins[ii][jj] = 0;
  }
  int ****pbins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    pbins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      pbins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        pbins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          pbins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }

  int **bins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    bins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) bins[ii][jj] = 0;
  }
  int ****bins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    bins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      bins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        bins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          bins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }

  double **dbins = new double*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    dbins[ii] = new double[nInputs];
    for (jj = 0; jj < nInputs; jj++) dbins[ii][jj] = 0;
  }
  double ****dbins2 = new double***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    dbins2[jj] = new double**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      dbins2[jj][jj2] = new double*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        dbins2[jj][jj2][ii] = new double[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          dbins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }

  int  samInc = 10000;
  long itmp;
  if (UParamsNumToUse == 1)
  {
    if (samInc * dnSamples > 100000) samInc = 100000 / dnSamples; 
    maxSamples = maxSamples / samInc;
    maxSamples = maxSamples * samInc;
  }
  else
  {
    if (maxSamples < 5000000) maxSamples = 5000000;
    itmp = samInc * dnSamples * UParamsNumToUse;
    if (itmp > 100000) samInc = 100000 / (dnSamples * UParamsNumToUse); 
    if (samInc < 1000) samInc = 1000;
    maxSamples = maxSamples / (samInc * UParamsNumToUse);
    if (maxSamples < 5) maxSamples = 5;
    maxSamples = maxSamples * samInc * UParamsNumToUse;
    while (maxSamples > 100000000)
    {
      maxSamples /= 2;
      samInc /= 2;
    }
    printOutTS(PL_INFO,"MCMC maxSamples         = %d\n",maxSamples);
    printOutTS(PL_INFO,"MCMC sample increment   = %d\n",samInc);
    printOutTS(PL_INFO,"MCMC UParam sample size = %d\n",UParamsNumToUse);
  }

  int      methodSave, *SSS;
  double   *XXX, *YYY;
  Sampling *sampler;
  psVector  vecLB, vecUB, vecOut, vecYT; 
  psIVector vecST;
  if (noPDF == 1)
  {
    printOutTS(PL_INFO,
               "MCMC_BF2 INFO: no PDF, use uniform for priors.\n");
    sampler = SamplingCreateFromID(PSUADE_SAMP_MC);
    sampler->setInputBounds(nInputs, lower, upper);
    sampler->setOutputParams(1);
    cnt = maxSamples / UParamsNumToUse;
    sampler->setSamplingParams(cnt, 1, 1);
    sampler->initialize(0);
    vecOut.setLength(cnt*nInputs);
    XXX = vecOut.getDVector();
    vecYT.setLength(cnt);
    YYY = vecYT.getDVector();
    vecST.setLength(cnt);
    SSS = vecST.getIVector();
    sampler->getSamples(cnt,nInputs,1,XXX,YYY,SSS);
    delete sampler;
  }
  else
  {
    printOutTS(PL_INFO,
               "MCMC_BF2 INFO: has PDF, draw sample from distribution.\n");
    dataPtr->getParameter("method_sampling", pPtr);
    methodSave = pPtr.intData_;
    dataPtr->updateMethodSection(PSUADE_SAMP_MC,-1,-1,-1,-1);
    PDFManager *pdfman = new PDFManager();
    pdfman->initialize(dataPtr);
    vecLB.load(nInputs, lower);
    vecUB.load(nInputs, upper);
    cnt = maxSamples / UParamsNumToUse;
    vecOut.setLength(cnt*nInputs);
    pdfman->genSample(cnt, vecOut, vecLB, vecUB);
    dataPtr->updateMethodSection(methodSave,-1,-1,-1,-1);
    delete pdfman;
  }

  //  allocate storage for inference
  psVector vecXS, vecYS, vecSamStd, vecXDesign, vecYDesign;
  psVector vecYDesStd, vecInfSamInp, vecInfSamOut, vecInfTmpOut;

  vecXS.setLength(samInc*dnSamples*UParamsNumToUse*nInputs);
  double *XSample = vecXS.getDVector();
  vecYS.setLength(samInc*dnSamples*UParamsNumToUse*nOutputs);
  double *YSample = vecYS.getDVector();
  vecSamStd.setLength(samInc*dnSamples*UParamsNumToUse*nOutputs);
  double *YSamStd = vecSamStd.getDVector();
  vecXDesign.setLength(samInc*dnSamples*UParamsNumToUse*nInputs);
  double *XDesign = vecXDesign.getDVector();
  vecYDesign.setLength(samInc*dnSamples*UParamsNumToUse*nOutputs);
  double *YDesign = vecYDesign.getDVector();
  vecYDesStd.setLength(samInc*dnSamples*UParamsNumToUse*nOutputs);
  double *YDesStd = vecYDesStd.getDVector();
  vecInfSamInp.setLength(maxSamples*nInputs);
  double *inferenceSamIns = vecInfSamInp.getDVector();
  vecInfSamOut.setLength(maxSamples);
  double *inferenceSamOut = vecInfSamOut.getDVector();
  vecInfTmpOut.setLength(maxSamples);
  double *inferenceTmpOut = vecInfTmpOut.getDVector();
  psIVector randIVec;
  double stdev, stdv2, YT, YT1, YT2;
  int    dcnt, irand, nInpsActive = 0, passCnt;

  //  duplicate the generated sample UParamsNumToUse-1 times
  for (ii = 0; ii < maxSamples/UParamsNumToUse; ii++)
  {
    for (jj = 0; jj < UParamsNumToUse; jj++)
    {
      cnt = (ii * UParamsNumToUse + jj) * nInputs;
      for (kk = 0; kk < nInputs; kk++)
        inferenceSamIns[cnt+kk] = vecOut[ii*nInputs+kk]; 
    }
  }

  if (UParamsSample.ncols() > 0)
  {
    randIVec.setLength(UParamsNumToUse);
    for (int mm = 0; mm < UParamsNumToUse; mm++)
    {
      if (UParamsNumToUse >= UParamsSample.nrows())
           ii = mm;
      else ii = PSUADE_rand() % UParamsSample.nrows();
      randIVec[mm] = ii;
      printf("Uncertain sample index %d = %d\n", mm+1, randIVec[mm]);
    }
    for (ii = 0; ii < maxSamples; ii+=UParamsNumToUse)
    {
      cnt = 0;
      for (kk = 0; kk < nInputs; kk++)
      {
        if ((rsIndices.length() == 0 || 
            (rsIndices.length() > 0 && rsIndices[kk] >=1000)))
        {
          for (jj = 0; jj < UParamsNumToUse; jj++)
          {
            ii2 = randIVec[jj];     
            ddata = UParamsSample.getEntry(ii2,cnt);
            inferenceSamIns[(ii+jj)*nInputs+kk] = ddata;
          }
          cnt++;
        }
      }
    }
  }

  for (ii2 = 0; ii2 < nInputs; ii2++)
  {
    if ((rsIndices.length() == 0 || 
        (rsIndices.length() > 0 && rsIndices[ii2] >=0 &&
         rsIndices[ii2] < 1000)) &&
        (dParams.length() == 0 || dParams[ii2] == 0)) nInpsActive++;
  }

  //  perform inference
  int    ss, sTotal=0, runSize, index, nFail, printStep=0, thriceFlag=0;
  double Ymin = PSUADE_UNDEFINED, Ymax = -PSUADE_UNDEFINED;
  psVector vecTmpMeans, vecTmpStds;
  vecTmpMeans.setLength(100*nInputs);
  double *tmpMeans = vecTmpMeans.getDVector();
  vecTmpStds.setLength(100*nInputs);
  double *tmpStds = vecTmpStds.getDVector();

  printOutTS(PL_INFO, "MCMC_BF2 Inference begins ... \n");
  fflush(stdout);
  fp = NULL;

  while (sTotal < maxSamples)
  {
    cnt = (sTotal+1) / (maxSamples/20);
    if (cnt != printStep)
    {
      printOutTS(PL_INFO, "%3.0f%% ",5.0*cnt);
      fflush(stdout);
    }

    for (int mm = 0; mm < samInc; mm++)
    {
      dcnt = 0;
      index = mm * dnSamples * UParamsNumToUse;
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
        for (kk2 = 0; kk2 < dnSamples; kk2++)
        {
          for (jj = 0; jj < UParamsNumToUse; jj++)
          {
            XSample[(index+kk2*UParamsNumToUse+jj)*nInputs+ii2] = 
             inferenceSamIns[(sTotal+mm*UParamsNumToUse+jj)*nInputs+ii2]; 
          }
        }
        if (dParams.length() > 0 && dParams[ii2] == 1)
        {
          for (kk2 = 0; kk2 < dnSamples; kk2++)
          {
            for (jj = 0; jj < UParamsNumToUse; jj++)
            {
              XSample[(index+kk2*UParamsNumToUse+jj)*nInputs+ii2] = 
                                     dSamInputs.getEntry(kk2, dcnt);
              XDesign[(index+kk2*UParamsNumToUse+jj)*dnInputs+dcnt] = 
                                     dSamInputs.getEntry(kk2, dcnt);
            }
          }
          dcnt++;
        }
        else if (rsIndices.length() > 0 && rsIndices[ii2] == -1) 
        {
          for (jj = 0; jj < samInc*UParamsNumToUse*dnSamples; jj++)
          {
            XSample[jj*nInputs+ii2] = rsValues[ii2];
          }
        }
      }
    }
    for (jj = 0; jj < samInc*UParamsNumToUse*dnSamples*nOutputs; jj++)
    {
      YDesign[jj] = YDesStd[jj] = 0.0;
      YSample[jj] = YSamStd[jj] = 0.0;
    }
    runSize = samInc * dnSamples * UParamsNumToUse;
    for (ii2 = 0; ii2 < nOutputs; ii2++)
    {
      if (rsErrFlag == 1)
      {
        faPtrs[ii2]->evaluatePointFuzzy(runSize, 
                             XSample,&YSample[ii2*runSize],
                             &YSamStd[ii2*runSize]);
        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
        {
          faPtrs1[ii2]->evaluatePointFuzzy(runSize,
                                XDesign,&YDesign[ii2*runSize],
                                &YDesStd[ii2*runSize]);
        }
        else if (discFuncConstantMeans.length() > 0 &&
                 discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < runSize; kk2++)
          {
            YDesign[ii2*runSize+kk2] = discFuncConstantMeans[ii2];
            YDesStd[ii2*runSize+kk2] = 0.0;
          }
        }
      }
      else
      {
        faPtrs[ii2]->evaluatePoint(runSize,XSample,
                                   &YSample[ii2*runSize]);
        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
        {
          faPtrs1[ii2]->evaluatePoint(runSize,XDesign,
                                    &YDesign[ii2*runSize]);
        }
        else if (discFuncConstantMeans.length() > 0 &&
                 discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < runSize; kk2++)
            YDesign[ii2*runSize+kk2] = discFuncConstantMeans[ii2];
        }
        for (kk2 = 0; kk2 < dnSamples*samInc; kk2++)
          YSamStd[ii2*runSize+kk2] = YDesStd[ii2*runSize+kk2] = 0.0;
      }
    }

    for (ii = 0; ii < samInc; ii++)
    {
      index = ii * UParamsNumToUse;
      for (jj = 0; jj < UParamsNumToUse; jj++)
      {
        inferenceSamOut[sTotal+index+jj] = 0.0;
        for (ii2 = 0; ii2 < nOutputs; ii2++)
        {
          for (kk2 = 0; kk2 < dnSamples; kk2++)
          {
            YT1 = YSample[ii2*runSize+index*dnSamples+kk2*UParamsNumToUse+jj] +
                  YDesign[ii2*runSize+index*dnSamples+kk2*UParamsNumToUse+jj];
            stdev = YSamStd[ii2*runSize+index*dnSamples+kk2*UParamsNumToUse+jj];
            stdv2 = YDesStd[ii2*runSize+index*dnSamples+kk2*UParamsNumToUse+jj];
            YT2 = pow((YT1-dSamMeans.getEntry(kk2,ii2)),2.0) /
                 (pow(dSamStdvs.getEntry(kk2,ii2),2.0) +
                  stdev*stdev + stdv2*stdv2);
            inferenceSamOut[sTotal+index+jj] += YT2;
          }
        }
        //*/ inferenceSamOut[sTotal+jj] /= (dnSamples*nOutputs);
        if (combFlag == 1)
          inferenceSamOut[sTotal+index+jj] /= (double) dnSamples;
      }
    }
    sTotal += (samInc * UParamsNumToUse);

    nFail = 0;
    for (jj = 0; jj < samInc*UParamsNumToUse; jj++)
    {
      dcnt = sTotal - (samInc * UParamsNumToUse) + jj;
      YT1 = constrPtr->evaluate(&inferenceSamIns[dcnt*nInputs],
                                inferenceSamOut[dcnt], status);
      if (status == 0)
      {
        inferenceSamOut[dcnt] = 0.0;
        nFail++;
      }
    }

    if (cnt != printStep)
    {
      printStep++;
      for (jj = 0; jj < sTotal; jj++) 
        inferenceTmpOut[jj] = inferenceSamOut[jj];
      Ymin = PSUADE_UNDEFINED;
      for (jj = 0; jj < sTotal; jj++) 
      {
        ddata = inferenceTmpOut[jj];
        if (ddata < Ymin) Ymin = ddata;
        if (ddata > Ymax) Ymax = ddata;
      }
      for (jj = 0; jj < sTotal; jj++) 
        inferenceTmpOut[jj] = exp(-0.5*(inferenceTmpOut[jj]-Ymin));

      printOutTS(PL_INFO,"\n");
      passCnt = 0;
      printOutTS(PL_INFO,"Convergence Checking (%d) =========>\n",sTotal);
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
        if ((rsIndices.length() == 0 || 
             (rsIndices.length() > 0 && rsIndices[ii2] >=0 &&
              rsIndices[ii2] < 1000)) &&
            (dParams.length() == 0 || dParams[ii2] == 0))
        {
          YT = 0.0;
          for (jj = 0; jj < sTotal; jj++) YT += inferenceTmpOut[jj];
          YT1 = 0.0;
          for (jj = 0; jj < sTotal; jj++)
          {
            YT1 += inferenceSamIns[jj*nInputs+ii2] * 
                   inferenceTmpOut[jj] / YT; 
          }
          printOutTS(PL_INFO,"MCMC_BF2: input %3d mean    = %e\n", 
                     ii2+1, YT1);
          YT2 = 0.0;
          for (jj = 0; jj < sTotal; jj++)
          {
            YT2 += pow(inferenceSamIns[jj*nInputs+ii2]-YT1,2.0)*
                       inferenceTmpOut[jj]/YT; 
          }
          printOutTS(PL_INFO,"MCMC_BF2: input %3d std dev = %e\n",
                     ii2+1,sqrt(YT2));
          tmpMeans[(printStep-1)+ii2*100] = YT1; 
          tmpStds[(printStep-1)+ii2*100] = sqrt(YT2); 
          if (printStep > 2) 
          {
            jj = checkConvergence(3,&tmpMeans[ii2*100+printStep-3],
                     &tmpStds[ii2*100+printStep-3],
                     sTotal-samInc*UParamsNumToUse);
            if (jj == 1)
            {
              passCnt++;
              printf("MCMC_BF2: input %3d converged.\n",ii2+1);
            }
          }
        }
      }
      printOutTS(PL_INFO,"<========= Convergence Checking\n");
      if (passCnt == nInpsActive)
      {
        if (thriceFlag <= 1) 
          thriceFlag++;
        else if (thriceFlag >= 2 && sTotal > 0.2*maxSamples) 
          maxSamples = sTotal;
      }
    }
  }

  Ymin = PSUADE_UNDEFINED;
  Ymax = -PSUADE_UNDEFINED;
  for (ss = 0; ss < maxSamples; ss++)
  {
    ddata = inferenceSamOut[ss];
    if (ddata < Ymin) Ymin = ddata;
    if (ddata > Ymax) Ymax = ddata;
  }
  for (ss = 0; ss < maxSamples; ss++)
    inferenceSamOut[ss] = exp(-0.5*(inferenceSamOut[ss]-Ymin));

  index = -1;
  Ymax = -PSUADE_UNDEFINED;
  for (ss = 0; ss < maxSamples; ss++)
  {
    ddata = inferenceSamOut[ss] * exp(-0.5*Ymin);
    if (ddata > Ymax)
    {
      Ymax = ddata;
      index = ss;
    }
  } 
  if (index >= 0)
  {
    for (ii = 0; ii < nInputs; ii++) 
      Xmax[ii] = inferenceSamIns[index*nInputs+ii];
  }
  printf("Maximum likelihood estimated solution:\n");
  for (ii = 0; ii < nInputs; ii++) 
  {
    if ((rsIndices.length() == 0 || 
        (rsIndices.length() > 0 && rsIndices[ii] < 0)) &&
        (dParams.length() == 0 || dParams[ii] == 0))
      printf("Input %3d = %16.8e\n",ii+1,Xmax[ii]);
  }
  printf("Negative log likelihood (unnormalized) = %e\n", -log(Ymax));

  if (Ymax > 0)
  {
    Ymax = Ymax / exp(-0.5*Ymin);
    for (ss = 0; ss < maxSamples; ss++) inferenceSamOut[ss] /= Ymax;
  }

  int    ii3, index2;
  for (ss = 0; ss < maxSamples; ss++)
  {
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      YT1 = (inferenceSamIns[ss*nInputs+ii2] - lower[ii2]) * XRange[ii2];
      index = (int) (YT1 * nbins);
      if (index > nbins)
      {
        printOutTS(PL_ERROR,"MCMC binning error 1 in file %s, line %d.\n",
                   __FILE__, __LINE__);
        printOutTS(PL_ERROR,"Sample input %d = %e\n", ii2+1, 
                   inferenceSamIns[ss*nInputs+ii2]);
        printOutTS(PL_ERROR,"Sample input %d lower bound = %e\n", ii2+1, 
                   lower[ii2]);
        printOutTS(PL_ERROR,"Sample input %d scaled val  = %e\n",ii2+1,YT1);
        printOutTS(PL_ERROR,"Sample input %d bin number  = %d (<%d?)\n",
                   ii2+1,index,nbins);
      }
      if (index < 0)
      {
        printOutTS(PL_ERROR,"MCMC binning error 2 in file %s, line %d.\n",
                   __FILE__, __LINE__);
        printOutTS(PL_ERROR,"Sample input %d = %e\n", ii2+1, 
                   inferenceSamIns[ss*nInputs+ii2]);
        printOutTS(PL_ERROR,"Sample input %d lower bound = %e\n", ii2+1, 
                   lower[ii2]);
        printOutTS(PL_ERROR,"Sample input %d scaled val  = %e\n",ii2+1,YT1);
        printOutTS(PL_ERROR,"Sample input %d bin number  = %d (>0?)\n",
                   ii2+1,index);
      }
      if (index >= nbins) index = nbins - 1;
      if (index <  0)     index = 0;
      dbins[index][ii2] += inferenceSamOut[ss];
      pbins[index][ii2]++;
    }
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      YT1 = (inferenceSamIns[ss*nInputs+ii2]-lower[ii2])*XRange[ii2];
      index = (int) (YT1 * nbins);
      if (index >= nbins) index = nbins - 1;
      if (index < 0)      index = 0;
      for (ii3 = 0; ii3 < nInputs; ii3++)
      {
        YT2 = (inferenceSamIns[ss*nInputs+ii3]-lower[ii3])*XRange[ii3];
        index2 = (int) (YT2 * nbins);
        if (index2 >= nbins) index2 = nbins - 1;
        if (index2 < 0)      index2 = 0;
        dbins2[index][index2][ii2][ii3] += inferenceSamOut[ss];
        pbins2[index][index2][ii2][ii3]++;
      }
    }
  }

  Ymax = 0;
  for (kk = 0; kk < nbins; kk++)
    for (ii2 = 0; ii2 < nInputs; ii2++) 
      if (dbins[kk][ii2] > Ymax) Ymax = dbins[kk][ii2];
  for (kk = 0; kk < nbins; kk++)
    for (ii2 = 0; ii2 < nInputs; ii2++) 
      dbins[kk][ii2] = dbins[kk][ii2] / Ymax * maxSamples;
  for (kk = 0; kk < nbins; kk++)
    for (ii2 = 0; ii2 < nInputs; ii2++) bins[kk][ii2] = (int) dbins[kk][ii2];
  Ymax = 0;
  for (jj = 0; jj < nbins; jj++)
    for (kk = 0; kk < nbins; kk++)
      for (ii2 = 0; ii2 < nInputs; ii2++)
        for (ii3 = 0; ii3 < nInputs; ii3++)
          if (dbins2[jj][kk][ii2][ii3] > Ymax) 
            Ymax = dbins2[jj][kk][ii2][ii3];
  for (jj = 0; jj < nbins; jj++)
    for (kk = 0; kk < nbins; kk++)
      for (ii2 = 0; ii2 < nInputs; ii2++)
        for (ii3 = 0; ii3 < nInputs; ii3++)
          dbins2[jj][kk][ii2][ii3] = dbins2[jj][kk][ii2][ii3] / Ymax *
                                     maxSamples;
  for (jj = 0; jj < nbins; jj++)
    for (kk = 0; kk < nbins; kk++)
      for (ii2 = 0; ii2 < nInputs; ii2++)
        for (ii3 = 0; ii3 < nInputs; ii3++)
          bins2[jj][kk][ii2][ii3] = (int) dbins2[jj][kk][ii2][ii3];
     
  double psum;
  for (ii = 0; ii < nInputs; ii++)
  {
    if ((rsIndices.length() == 0 || (rsIndices.length() > 0 && 
         rsIndices[ii] >=0 && rsIndices[ii] < 1000)) &&
        (dParams.length() == 0 || dParams[ii] == 0))
    {
      printOutTS(PL_INFO,
                 "MCMC_BF2: input %3d value at peak of likelihood = %e\n",
                 ii+1, Xmax[ii]);
      psum = 0.0;
      for (ss = 0; ss < maxSamples; ss++) psum += inferenceSamOut[ss];
      YT = 0.0;
      for (ss = 0; ss < maxSamples; ss++)
        YT += inferenceSamIns[ss*nInputs+ii]*inferenceSamOut[ss]/psum; 
      means_[ii] = YT;
      printOutTS(PL_INFO,"MCMC_BF2: input %3d mean    = %e\n", ii+1, YT);
      YT = 0.0;
      for (ss = 0; ss < maxSamples; ss++)
        YT += pow(inferenceSamIns[ss*nInputs+ii]-means_[ii],2.0)*
              inferenceSamOut[ss]/psum; 
      sigmas_[ii] = sqrt(YT);
      printOutTS(PL_INFO,"MCMC_BF2: input %3d std dev = %e\n",ii+1,sigmas_[ii]);
      mostLikelyInput_[ii] = Xmax[ii];
    }
  }
 
  for (ii = 0; ii < nInputs; ii++) XRange[ii] = 1.0 / XRange[ii];
  genMatlabFile(nInputs,lower,upper,XRange,nPlots,plotIndices,nbins,
        pbins,pbins2,bins,bins2,qData,0,0,NULL,NULL,Xmax,Ymin);

  if (genPosteriors == 1)
  {
    int    maxPostSam=50000;
    double dmax = 0.0;
    for (ss = 0; ss < maxSamples; ss++) 
      if (inferenceSamOut[ss] > dmax) dmax = inferenceSamOut[ss];
    if (dmax == 0)
    {
      printOutTS(PL_ERROR,
           "MCMC_BF: ERROR encountered in posterior sample generation.\n");
      fp = NULL;
    }
    else 
    {
      fp = fopen("MCMCPostSample", "w");
      for (ss = 0; ss < maxSamples; ss++) 
        inferenceSamOut[ss] = inferenceSamOut[ss] / dmax;
    }
    if (fp != NULL)
    {
      fprintf(fp, "PSUADE_BEGIN\n");
      fprintf(fp, "%d %d\n", maxPostSam, nInputs);
      if (qData.strArray_ != NULL)
      {
        fprintf(fp, "# ");
        for (jj = 0; jj < nInputs; jj++)
           fprintf(fp,"%s ", qData.strArray_[jj]);
        fprintf(fp, "\n");
      }
      int    count=0, ichoose, nValid=0;
      double dchoose;
      psIVector vecValid;
      vecValid.setLength(maxSamples);
      for (ss = 0; ss < maxSamples; ss++) 
        if (inferenceSamOut[ss] > 1e-5) vecValid[nValid++] = ss;
      while (count < maxPostSam)
      {
        ichoose = PSUADE_rand() % nValid;
        ichoose = vecValid[ichoose];
        dchoose = PSUADE_drand();
        if (inferenceSamOut[ichoose] > dchoose)
        {
          fprintf(fp, "%d ", count+1);
          for (jj = 0; jj < nInputs; jj++)
          {
             if ((rsIndices.length() == 0 || rsIndices[jj] >= 0) &&
                 (dParams.length() == 0 || dParams[jj] == 0))
             {
                fprintf(fp, "%e ", 
                    inferenceSamIns[ichoose*nInputs+jj]);
             }
             else if (rsIndices.length() > 0 && rsIndices[jj] == -1)
                fprintf(fp, "%e ", rsValues[jj]);
             else if (dParams.length() > 0 && dParams[jj] != 0)
                fprintf(fp, "%e ", 0.5 * (upper[jj] + lower[jj]));
             else if (rsIndices.length() > 0 && rsIndices[jj] == -2)
                fprintf(fp, "0 ");
          }
          fprintf(fp, "\n");
          count++;
        }
      }
      fprintf(fp, "PSUADE_END\n");
      psVector vecXSam, vecYSam, vecXDes, vecYDes;
      vecXSam.setLength(nInputs*dnSamples);
      vecYSam.setLength(nOutputs*dnSamples);
      vecXDes.setLength(dnInputs*dnSamples);
      vecYDes.setLength(nOutputs*dnSamples);
      YSample = vecYSam.getDVector();
      YDesign = vecYDes.getDVector();
      fprintf(fp, "Optimal parameter values: \n");
      dcnt = 0;
      for (ii2 = 0; ii2 < nInputs; ii2++) 
      {
        if (dParams.length() > 0 && dParams[ii2] == 1) 
        {
          fprintf(fp, "%16.8e\n", dSamInputs.getEntry(0,dcnt));
          dcnt++;
        }
        else if (rsIndices.length() > 0 && rsIndices[ii2] == -1)
          fprintf(fp, "%16.8e\n", rsValues[ii2]);
        else if (rsIndices.length() > 0 && rsIndices[ii2] == -2)
          fprintf(fp, "0.0\n");
        else
          fprintf(fp, "%16.8e\n", Xmax[ii2]);
      }
      fprintf(fp,
        "MLE Statistics (prediction, exp data, sd, -loglikelihood)\n");
      for (kk2 = 0; kk2 < dnSamples; kk2++)
      {
        dcnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++) 
        {
          vecXSam[kk2*nInputs+ii2] = Xmax[ii2];
          if (dParams.length() > 0 && dParams[ii2] == 1)
          {
            vecXSam[kk2*nInputs+ii2] = dSamInputs.getEntry(kk2,dcnt);
            vecXDes[kk2*dnInputs+dcnt] = dSamInputs.getEntry(kk2,dcnt);
            dcnt++;
          }
          if (rsIndices.length() > 0 && rsIndices[ii2] == -1)
            vecXSam[kk2*nInputs+ii2] = rsValues[ii2];
        }
      }
      for (ii2 = 0; ii2 < dnSamples*nOutputs; ii2++)
        vecYDes[ii2] = YSample[ii2] = 0.0;
      for (ii2 = 0; ii2 < nOutputs; ii2++) 
      {
        faPtrs[ii2]->evaluatePoint(dnSamples,vecXSam.getDVector(),
                                   &YSample[ii2*dnSamples]);
        if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
        {
          faPtrs1[ii2]->evaluatePoint(dnSamples,vecXDes.getDVector(),
                                      &YDesign[ii2*dnSamples]);
        }
        else if (discFuncConstantMeans.length() > 0 &&
                 discFuncConstantMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < dnSamples; kk2++)
            vecYDes[ii2*dnSamples+kk2] = discFuncConstantMeans[ii2];
        }
      }
      for (kk2 = 0; kk2 < dnSamples; kk2++)
      {
        for (ii2 = 0; ii2 < nOutputs; ii2++)
        {
          YT1 = YSample[ii2*dnSamples+kk2] + vecYDes[ii2*dnSamples+kk2];
          YT2 = dSamMeans.getEntry(kk2, ii2);
          stdev = dSamStdvs.getEntry(kk2, ii2);
          stdv2 = 0.5 * pow(YT1 - YT2, 2.0) / (stdev * stdev);
          fprintf(fp,"%4d: %16.8e %16.8e %16.8e %16.8e\n",kk2+1,
                  YT1,YT2,stdev,stdv2);
        }
      }
      fclose(fp);
    }
    printOutTS(PL_INFO,
         "MCMC: 'MCMCPostSample' file has a posterior sample.\n");
  }

  int    nInps, nOuts, nSams, *states;
  double *allOuts;
  char   **oNames;
  PsuadeData *filePtr1, *filePtr2;
  if (modelFormFlag == 1)
  {
    sprintf(charString, "psDiscrepancyModel1");
    filePtr1 = new PsuadeData();
    status = filePtr1->readPsuadeFile(charString);
    if (status != 0)
    {
       printOutTS(PL_ERROR,
            "MCMC ERROR: cannot read file %s in PSUADE format.\n",
            charString);
       exit(1);
    }
  }
  if (modelFormFlag == 1 && status == 0)
  {
    filePtr1->getParameter("input_ninputs", pPtr);
    nInps = pPtr.intData_;
    filePtr1->getParameter("output_noutputs", pPtr);
    nOuts = pPtr.intData_;
    filePtr1->getParameter("method_nsamples", pPtr);
    nSams = pPtr.intData_;
    filePtr1->getParameter("output_sample", pOutputs);
    unlink(charString);
    allOuts = new double[nOutputs * nSams];
    for (jj = 0; jj < nSams; jj++)
      allOuts[jj*nOutputs] = pOutputs.dbleArray_[jj];
    pOutputs.clean();
    oNames = new char*[nOutputs];
    oNames[0] = new char[100];
    sprintf(oNames[0], "Y1");
    for (ii = 1; ii < nOutputs; ii++)
    {
      filePtr2 = new PsuadeData();
      sprintf(charString, "psDiscrepancyModel%d", ii+1);
      status = filePtr2->readPsuadeFile(charString);
      if (status != 0) break;
      filePtr2->getParameter("input_ninputs", pPtr);
      if (pPtr.intData_ != nInps) break;
      filePtr2->getParameter("output_noutputs", pPtr);
      if (pPtr.intData_ != nOuts) break;
      filePtr2->getParameter("method_nsamples", pPtr);
      if (pPtr.intData_ != nSams) break;
      filePtr2->getParameter("output_sample", pOutputs);
      delete filePtr2;
      unlink(charString);
      for (jj = 0; jj < nSams; jj++)
         allOuts[jj*nOutputs+ii] = pOutputs.dbleArray_[jj];
      pOutputs.clean();
      oNames[ii] = new char[100];
      sprintf(oNames[ii], "Y%d", ii+1);
    }
    if (nOutputs == 1)
    {
       sprintf(charString, "psDiscrepancyModel");
       filePtr1->writePsuadeFile(charString, 0);
    }
    else if (ii == nOutputs)
    {
      states = new int[nSams];
      for (jj = 0; jj < nSams; jj++) states[jj] = 1;
      filePtr1->updateOutputSection(nSams,nOutputs,allOuts,states,oNames);
      sprintf(charString, "psDiscrepancyModel");
      filePtr1->writePsuadeFile(charString, 0);
      printOutTS(PL_INFO,
           "MCMC INFO: a sample (inputs/outputs) the discrepancy model\n");
      printOutTS(PL_INFO,"           is now in psDiscrepancyModel.\n");
      delete [] states;
      for (ii = 1; ii < nOutputs; ii++) delete [] oNames[ii];
    }
    else
    {
      printOutTS(PL_INFO,
           "MCMC INFO: unsuccessful creation of discrepancy sample file\n");
    }
    delete [] oNames[0];
    delete [] oNames;
    delete filePtr1;
    delete [] allOuts;
  }

  // clean up
  if (plotIndices != NULL) delete [] plotIndices;
  if (faPtrs != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++)
      if (faPtrs[ii] != NULL) delete faPtrs[ii];
    delete [] faPtrs;
  }
  if (faPtrs1 != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++)
      if (faPtrs1[ii] != NULL) delete faPtrs1[ii];
    delete [] faPtrs1;
  }
  for (ii = 0; ii < nbins; ii++) 
  {
    delete [] bins[ii];
    delete [] dbins[ii];
    delete [] pbins[ii];
  }
  delete [] bins;
  delete [] pbins;
  delete [] dbins;
  for (jj = 0; jj < nbins; jj++)
  {
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      for (ii = 0; ii < nInputs; ii++) delete [] bins2[jj][jj2][ii];
      for (ii = 0; ii < nInputs; ii++) delete [] dbins2[jj][jj2][ii];
      for (ii = 0; ii < nInputs; ii++) delete [] pbins2[jj][jj2][ii];
      delete [] bins2[jj][jj2];
      delete [] pbins2[jj][jj2];
      delete [] dbins2[jj][jj2];
    }
    delete [] bins2[jj];
    delete [] dbins2[jj];
    delete [] pbins2[jj];
  }
  delete [] bins2;
  delete [] dbins2;
  delete [] pbins2;
  delete constrPtr;
  return 0.0;
}

// ************************************************************************
// perform MCMC analysis using Metropolis Hasting
// ------------------------------------------------------------------------
double MCMCAnalyzer::analyze_mh(aData &adata)
{
  int    ii, ii2, jj, kk, status, cnt, iOne=1, iZero=0, nInputs, nOutputs;
  int    nSamples, maxPts=257, nbins, printLevel, faType, genPosteriors=0;
  int    maxSamples,burnInSamples,modelFormFlag=0,numChains=3, count;
  int    nPlots, *plotIndices=NULL, *pdfTypes, *designParams=NULL;
  int    *rsIndices=NULL, freq=1, dnSamples=0, dnInputs=0;
  int    nChainGood=0, masterCount=1, combFlag=0;
  double *dSamInputs=NULL, *dSamMeans=NULL, *dSamStdevs=NULL;
  double *sampIns=NULL, *samOuts=NULL, *samOut1=NULL, *xLower=NULL;
  double *xUpper=NULL, *pdfMeans=NULL, *pdfStdvs=NULL, *rsValues=NULL;
  double *discFuncConstantMeans=NULL, *discFuncConstantStds=NULL;
  double psrfThreshold=1.05, dstatus, ddata;
  char   lineIn[1001], charString[1001], *cString;
  FILE   *fp=NULL;
  pData      pPtr, qData, pOutputs;
  FuncApprox **faPtrs=NULL, **faPtrs1=NULL;
  PDFBase    **inputPDFs;
  FunctionInterface *funcIO=NULL;
  PsuadeData *dataPtr=NULL, *ioOutFile=NULL;
  RSConstraints *constrPtr;

  printLevel = adata.printLevel_;
  nInputs_   = nInputs  = adata.nInputs_;
  nOutputs_  = nOutputs = adata.nOutputs_;
  nSamples   = adata.nSamples_;
  sampIns    = adata.sampleInputs_;
  samOuts    = adata.sampleOutputs_;
  xLower     = adata.iLowerB_;
  xUpper     = adata.iUpperB_;
  pdfTypes   = adata.inputPDFs_;
  pdfMeans   = adata.inputMeans_;
  pdfStdvs   = adata.inputStdevs_;
  dataPtr    = adata.ioPtr_;
  if (dataPtr != NULL) dataPtr->getParameter("input_names", qData);

  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO,"*                     MCMC Optimizer (MH)\n");
  printEquals(PL_INFO, 0);
  if (printLevel > 0)
  {
    printOutTS(PL_INFO,"TO GAIN ACCESS TO DIFFERENT OPTIONS: TURN ON\n\n");
    printOutTS(PL_INFO," * ana_expert to finetune MCMC parameters, \n");
    printOutTS(PL_INFO,
         "   (e.g. sample size for burn-in can be adjusted).\n");
    printOutTS(PL_INFO,
         " * rs_expert to customize response surface for MCMC,\n");
    printDashes(PL_INFO,0);
    printOutTS(PL_INFO,
         "FEATURES AVAILABLE IN THE CURRENT VERSION OF MCMC:\n");
    printOutTS(PL_INFO," * Support other than uniform prior distributions\n");
    printOutTS(PL_INFO,
         " * Support likelihood functions from multiple outputs\n");
    printOutTS(PL_INFO,
         " * Option to set some inputs as design parameters\n");
    printOutTS(PL_INFO,
         "   - to be specified in the observation data spec file\n");
    printOutTS(PL_INFO,
         " * Option to disable some parameters (set to default)\n");
    printOutTS(PL_INFO,
         "   - in case these parameters are not to be calibrated\n");
    printOutTS(PL_INFO,
         "   - use rs_index_file in PSUADE's ANALYSIS section\n");
    printOutTS(PL_INFO,"   - not available with discrepancy modeling\n");
    printOutTS(PL_INFO," * Option to generate a posterior sample\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode\n");
    printOutTS(PL_INFO,
         " * MCMC can be terminated gracefully by creating a file\n");
    printOutTS(PL_INFO,
         "   named 'psuade_stop' in the same directory during the\n");
    printOutTS(PL_INFO,"   run (in case it takes too long).\n");
    printOutTS(PL_INFO,
         " * For multi-modal posteriors, a large number of chains\n");
    printOutTS(PL_INFO,
         "   may be needed. The number of chains can be adjusted\n");
    printOutTS(PL_INFO,"   in ana_expert mode.\n");
    printEquals(PL_INFO, 0);
  }

  cleanUp();

  // get experimental data information from the spec file
  dstatus = readSpecFile(nInputs, nOutputs, &dnSamples, &dnInputs, 
                 &designParams, &dSamInputs, &dSamMeans, &dSamStdevs,
                 combFlag, printLevel);
  if (dstatus != 0.0) return PSUADE_UNDEFINED;

  if (dataPtr != NULL)
  {
    rsIndices = new int[nInputs];
    for (ii = 0; ii < nInputs; ii++) rsIndices[ii] = 0;
    rsValues  = new double[nInputs];
    checkAllocate(rsValues, "rsValues in MCMC::analyze");
    status = readIndexFile(dataPtr,nInputs,designParams,rsIndices, 
                           rsValues);
    count = 0;
    for (ii = 0; ii < nInputs; ii++) count += rsIndices[ii];
    if (count == 0)
    {
      delete [] rsIndices;
      delete [] rsValues;
      rsIndices = NULL;
      rsValues = NULL;
    }
    if (status == -1)
    {
      if (rsIndices != NULL) delete [] rsIndices;
      if (rsValues  != NULL) delete [] rsValues;
      return PSUADE_UNDEFINED;
    }
    dataPtr->getParameter("ana_rstype", pPtr);
    faType = pPtr.intData_;
  }
  else
  {
    printOutTS(PL_INFO,
         "MCMC INFO: since ioPtr=NULL, assume MARS as reponse surface.\n");
    faType = PSUADE_RS_MARS;
  }
  printOutTS(PL_INFO,
       "MCMC INFO: CREATING RESPONSE SURFACES FOR ALL OUTPUTS.\n");
  faPtrs = new FuncApprox*[nOutputs];
  samOut1 = new double[nSamples];
  checkAllocate(samOut1, "samOut1 in MCMC::analyze");
  for (ii = 0; ii < nOutputs; ii++)
  {
    faType = -1;
    printOutTS(PL_INFO,
         "MCMC INFO: CREATING RESPONSE SURFACE FOR OUTPUT %d.\n",ii+1);
    faPtrs[ii] = genFA(faType, nInputs, iZero, nSamples);
    faPtrs[ii]->setNPtsPerDim(16);
    faPtrs[ii]->setBounds(xLower, xUpper);
    faPtrs[ii]->setOutputLevel(0);
    for (kk = 0; kk < nSamples; kk++) samOut1[kk] = samOuts[kk*nOutputs+ii];

    status = faPtrs[ii]->initialize(sampIns, samOut1);
    if (status != 0)
    {
      printOutTS(PL_ERROR,
           "MCMC ERROR: Unable to create response surface.\n");
      printOutTS(PL_ERROR,"            Consult PSUADE developers.\n");
      for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
      delete [] faPtrs;
      delete [] samOut1;
      if (rsIndices != NULL) delete [] rsIndices;
      if (rsValues  != NULL) delete [] rsValues;
      return PSUADE_UNDEFINED;
    }
  }
  delete [] samOut1;

  //    ==> burnInSamples, maxSamples, nbins, plotIndices, nPlots
  maxSamples = 10000;
  burnInSamples = maxSamples / 2;
  nbins = 20;
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"*** CURRENT SETTINGS OF MCMC PARAMETERS: \n\n");
  printOutTS(PL_INFO,"MCMC Burn-in sample size      (default) = %d\n", 
             burnInSamples);
  printOutTS(PL_INFO,"MCMC sample increment         (default) = %d\n", 
             maxSamples);
  printOutTS(PL_INFO,"MCMC no. of bins in histogram (default) = %d\n",nbins);
  printOutTS(PL_INFO,
     "NOTE: sample increment - sample size to run before convergence check\n");
  printOutTS(PL_INFO,
     "NOTE: histogram nBins  - define granularity of histogram bar graph\n");
  printOutTS(PL_INFO, 
     "Turn on ana_expert mode to change these default settings.\n\n");
  if (psAnaExpertMode_ == 1)
  {
    sprintf(charString,"Enter sample increment (100 - 50000): ");
    maxSamples = getInt(100, 200000, charString);
    burnInSamples = maxSamples / 2;
    sprintf(charString,"Enter the number of histogram bins (5 - 25) : ");
    nbins = getInt(5, 50, charString);
  }
  if (psAnaExpertMode_ == 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
         "*** OPTION TO CREATE POSTERIOR PLOTS FOR SELECTED INPUTS:\n\n");
    printOutTS(PL_INFO,
         "MCMC will create MATLAB files for the posterior distributions.\n");
    printOutTS(PL_INFO,
         "You can choose to generate posterior plots for all inputs, or \n");
    printOutTS(PL_INFO,
         "just a selected few (in case there are too many inputs).\n");
    printf("Select inputs for which posterior plots are to be generated.\n");
    sprintf(charString,"Enter input number (-1 for all, 0 to terminate) : ");
    kk = 1;
    plotIndices = new int[nInputs];
    nPlots = 0;
    while (kk != 0 || nPlots < 1)
    {
      kk = getInt(-1, nInputs, charString);
      if (kk == -1)
      {
        for (ii = 0; ii < nInputs; ii++)
        {
          if (rsIndices == NULL || rsIndices[ii] >= 0)
            if (designParams == NULL || designParams[ii] == 0) 
              plotIndices[nPlots++] = ii;
        }
        break;
      }
      if (kk != 0)
      {
        if (rsIndices != NULL && rsIndices[kk-1] < 0) 
          printOutTS(PL_ERROR,
               "Input %d has been fixed by the rs index file (no plot).\n",
               kk+1);
        else if (designParams != NULL && designParams[kk-1] == 1)
          printOutTS(PL_ERROR,
               "Input %d is a design parameter (no plot)\n",kk);
        else 
          plotIndices[nPlots++] = kk - 1;
      }
      if (kk == 0 && nPlots == 0)
        printOutTS(PL_ERROR,
             "You need to set at least 1 input for plotting posteriors.\n");
    }
    if (nPlots > 1) sortIntList(nPlots, plotIndices);
  }
  else
  {
    plotIndices = new int[nInputs];
    nPlots = 0;
    for (ii = 0; ii < nInputs; ii++) 
      if (rsIndices == NULL || rsIndices[ii] >= 0)
        if (designParams == NULL || designParams[ii] == 0) 
          plotIndices[nPlots++] = ii;
  }
  printOutTS(PL_INFO, 
       "MCMC Plot summary: input number to be plotted are (%d):\n",nPlots);
  for (ii = 0; ii < nPlots; ii++)
    printOutTS(PL_INFO, "   Input %4d\n", plotIndices[ii]+1);

  // option to add discrepancy function and a posterior sample
  // ==> modelFormFlag, genPosteriors
  if (psAnaExpertMode_ == 1)
  {
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,"*** OPTION TO ADD A DISCREPANCY FUNCTION:\n\n");
    printOutTS(PL_INFO,
         "To use this feature, first make sure that the observation\n");
    printOutTS(PL_INFO,
         "data file specified earlier has design parameters specified\n");
    printOutTS(PL_INFO,
         "since the discrepancy function is to be a function of these\n");
    printOutTS(PL_INFO,
         "design parameters (if not, a constant discrepancy function\n");
    printOutTS(PL_INFO,"is to be created).\n");
    printOutTS(PL_INFO,
         "NOTE: if you don't know what this is, just say NO.\n");
    printf("===> Add discrepancy function ? (y or n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') modelFormFlag = 1;

    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,
       "*** OPTION TO CREATE A SAMPLE FROM THE POSTERIOR DISTRIBUTIONS:\n\n");
    printOutTS(PL_INFO,
         "In addition to generating the posterior distributions, you can\n");
    printOutTS(PL_INFO,
         "also draw a sample from these posteriors. The posterior sample\n");
    printOutTS(PL_INFO,
         "can be used as prior sample for another simulator/emulator.\n");
    printOutTS(PL_INFO,"NOTE: if you don't know what this is, just say no.\n");
    printf("==> Create posterior samples for the input parameters? (y/n) ");
    scanf("%s", charString);
    fgets(lineIn,1000,stdin);
    if (charString[0] == 'y') genPosteriors = 1;
  } 
  printEquals(PL_INFO, 0);

  // setup input PDF, if there is any
  if (printLevel > 2) 
    printOutTS(PL_INFO,"*** INFORMATION ON PARAMETER PRIOR DISTRIBUTIONS\n");
  inputPDFs = new PDFBase*[nInputs];
  for (ii = 0; ii < nInputs; ii++)
  {
    if ((rsIndices == NULL || rsIndices[ii] >= 0) &&
        (designParams == NULL || designParams[ii] == 0)) 
    {
      if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_NORMAL)
      {
        inputPDFs[ii] = (PDFBase *) new PDFNormal(pdfMeans[ii], 
                                                  pdfStdvs[ii]);
        if (printLevel > 2) 
           printOutTS(PL_INFO,
                "Parameter %3d has normal prior distribution (%e,%e)\n",
                ii+1, pdfMeans[ii], pdfStdvs[ii]);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_LOGNORMAL)
      {
        inputPDFs[ii] = (PDFBase *) new PDFLogNormal(pdfMeans[ii],
                                                     pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has lognormal prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_TRIANGLE)
      {
        inputPDFs[ii] = (PDFBase *) new PDFTriangle(pdfMeans[ii],
                                                    pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has triangle prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_BETA)
      {
        inputPDFs[ii] = (PDFBase *) new PDFBeta(pdfMeans[ii],pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has beta prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_WEIBULL)
      {
        inputPDFs[ii] = (PDFBase *) new PDFWeibull(pdfMeans[ii],
                                                   pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has Weibull prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_GAMMA)
      {
        inputPDFs[ii] = (PDFBase *) new PDFGamma(pdfMeans[ii],
                                                 pdfStdvs[ii]);
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has gamma prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == 1000+PSUADE_PDF_NORMAL)
      {
        inputPDFs[ii] = NULL;
        printOutTS(PL_INFO,
             "Parameter %3d: multi-parameter normal distribution.\n",ii+1);
        printOutTS(PL_INFO,"               curently not supported.\n");
        return -1.0;
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == 1000+PSUADE_PDF_LOGNORMAL)
      {
        inputPDFs[ii] = NULL;
        printOutTS(PL_INFO,
           "Parameter %3d: multi-parameter lognormal distribution.\n",ii+1);
        printOutTS(PL_INFO,"               curently not supported.\n");
        return -1.0;
      }
      else if (pdfTypes == NULL || pdfTypes[ii] == PSUADE_PDF_UNIFORM)
      {
        inputPDFs[ii] = NULL;
        if (printLevel > 2)
          printOutTS(PL_INFO,
               "Parameter %3d has uniform prior distribution.\n",ii+1);
      }
      else if (pdfTypes != NULL && pdfTypes[ii] == PSUADE_PDF_SAMPLE)
      {
        inputPDFs[ii] = NULL;
        printOutTS(PL_INFO,
           "Parameter %3d: user-provided distribution currently not\n",ii+1);
        printOutTS(PL_INFO,"               supported.\n");
        return -1.0;
      }
    }
  }
  if (printLevel > 2) printEquals(PL_INFO, 0);

  //    ==> funcIO, freq (how often to do simulation and emulation)
  maxPts = nbins * 5;
  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO,"*** SETTING PROPOSAL DISTRIBUTION RESOLUTION\n");
    printOutTS(PL_INFO,
         "Since MCMC uses many function evaluations to construct\n");
    printOutTS(PL_INFO,
         "the proposal distributions, you have the option to set\n");
    printOutTS(PL_INFO,
         "how many points are used to construct it in order to\n");
    printOutTS(PL_INFO,"keep the inference cost reasonable.\n");
    printf("Sample size to construct proposal distribution. Default is %d.\n",
           maxPts);
    sprintf(charString,"Enter new sample size (%d - %d): ",nbins*3,nbins*10);
    maxPts = getInt(nbins*3, nbins*10, charString);
    maxPts = maxPts / nbins * nbins;
    printOutTS(PL_INFO,"Proposal distribution sample size = %d.\n", maxPts);
  }
  if (mode_ == 1 && dataPtr != NULL) 
  {
    funcIO = createFunctionInterface(dataPtr);
    if (nSamples == 0)
      printOutTS(PL_INFO, "MCMC: DIRECT SIMULATION has been set up.\n");
    else
      printOutTS(PL_INFO,
         "MCMC: DIRECT SIMULATION PLUS RESPONSE SURFACE have been set up.\n");
    printOutTS(PL_INFO,
      "MCMC INFO: make sure simulation output is in the right form\n");
    printOutTS(PL_INFO,
      "     which should not have been translated (mean) nor scaled (std)\n");
    printOutTS(PL_INFO,
      "     unless you compute the error measure yourself, in which case\n");
    printOutTS(PL_INFO,
      "     you should have use nOutputs=1, mean=0, and std dev=1 in the\n");
    printOutTS(PL_INFO,"     spec file.\n");
    if (psAnaExpertMode_ == 1 && nSamples > 0)
    {
      printOutTS(PL_INFO,
        "Since MCMC uses many function evaluations, you have the option\n");
      printOutTS(PL_INFO,
        "to set how frequent the simulator is invoked (response surface\n");
      printOutTS(PL_INFO,
        "is used otherwise) in constructing the proposal distribution\n");
      printOutTS(PL_INFO,
        "at each MCMC iteration. A frequency of f means that each MCMC\n");
      printOutTS(PL_INFO,
        "step uses f simulator runs. These f simulator runs will then be\n");
      printOutTS(PL_INFO,
        "supplemented with evaluations from the given response surface.\n");
      printOutTS(PL_INFO,
        "The default is f=10 (if you do not know what this is, enter 10).\n");
      sprintf(charString,
        "Max. number of simulator runs per MCMC step (1 - %d, default=10)? ",
              maxPts);
      kk = getInt(1, maxPts, charString);
      freq = maxPts * nInputs / kk;
      if (freq * kk != maxPts) freq++;
    }
    else if (nSamples > 0) freq = maxPts / 10;
    else                   freq = 1;
    if (nSamples > 0)
      printOutTS(PL_INFO,
         "Frequency of invoking the simulator has been set to %d\n", freq);
  }
  if (funcIO == NULL && faPtrs == NULL)
  {
    printOutTS(PL_ERROR,
         "MCMC ERROR: missing simulator and sample data - cannot proceed.\n");
    return PSUADE_UNDEFINED;
  }

  double *discOutputs=NULL;
  if (modelFormFlag == 1)
  {
    int    *ExpSamStates, ind, ExpNSamples, dfaType, dnPerDim=16;
    double *dOneSample, expdata, simdata;
    double *ExpSamInputs, *tSamInputs, *settings;

    ExpNSamples   = dnSamples;
    ExpSamInputs  = new double[ExpNSamples*nInputs];
    discOutputs   = new double[ExpNSamples*nOutputs];
    ExpSamStates  = new int[ExpNSamples];
    discFuncConstantMeans = new double[nOutputs];
    discFuncConstantStds  = new double[nOutputs];
    for (ii2 = 0; ii2 < nOutputs; ii2++)
      discFuncConstantMeans[ii2] = discFuncConstantStds[ii2] = 
                                   PSUADE_UNDEFINED;

    printOutTS(PL_INFO,
         "*** SELECT RESPONSE SURFACE TYPE FOR DISCREPANCY FUNCTION:\n");
    dfaType = -1;
    while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
    {
      writeFAInfo(-1);
      sprintf(charString, "===> Enter your choice : ");
      dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
    }

    settings = new double[nInputs];
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      settings[ii2] = 0.5*(xLower[ii2] + xUpper[ii2]);
      if (rsIndices != NULL && rsIndices[ii2] < 0)
        settings[ii2] = rsValues[ii2];
    }

    faPtrs1 = new FuncApprox*[nOutputs];
    dOneSample = new double[nInputs];
    tSamInputs = new double[ExpNSamples*nInputs];
    int        askFlag = 0, *states=NULL;
    double     *tLowers = new double[nInputs];
    double     *tUppers = new double[nInputs];
    PsuadeData *dataPtr = new PsuadeData(); 
    char       **iNames;
    for (ii = 0; ii < nOutputs; ii++)
    {
      for (kk = 0; kk < dnSamples; kk++)
      {
        cnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          if (designParams != NULL && designParams[ii2] == 1)
          {
            dOneSample[ii2] = dSamInputs[kk*dnInputs+cnt];
            cnt++;
          }
          else dOneSample[ii2] = settings[ii2];
        }

        simdata = 0.0;
        if (psAnaExpertMode_ == 1 && askFlag == 0)
        {
          printOutTS(PL_INFO,
               "To create discrepancy functions, the calibration\n");
          printOutTS(PL_INFO,
               "parameters need to be set to some nominal values.\n");
          printOutTS(PL_INFO,
               "You can choose the nominal values, or it will be\n");
          printOutTS(PL_INFO,
               "set to the input means or mid points of the ranges.\n");
          printf( "Set nomininal values yourself ? (y or n) ");
          scanf("%s", charString);
          fgets(lineIn,1000,stdin);
          if (charString[0] == 'y')
          {
            for (ii2 = 0; ii2 < nInputs; ii2++)
            {
              if ((rsIndices == NULL || rsIndices[ii2] >= 0) &&
                  (designParams == NULL || designParams[ii2] == 0))
              {
                printOutTS(PL_INFO,
                     "Input %d has lower and upper bounds = %e %e\n",
                     ii2+1, xLower[ii2], xUpper[ii2]);
                sprintf(charString,"Nominal value for input %d : ",
                        ii2+1);
                dOneSample[ii2] = getDouble(charString);
                settings[ii2]   = dOneSample[ii2];
              }
            }
          }
          else
          {
            for (ii2 = 0; ii2 < nInputs; ii2++)
            {
              if ((rsIndices == NULL || rsIndices[ii2] >= 0) &&
                  (designParams == NULL || designParams[ii2] == 0))
              {
                printOutTS(PL_INFO,
                     "Nominal value for input %d = %e\n",
                     ii2+1, dOneSample[ii2]);
              }
            }
          }
          askFlag = 1;
        }

        if (funcIO != NULL)
          funcIO->evaluate(kk+1,nInputs,dOneSample,1,&simdata,0);
        else simdata = faPtrs[ii]->evaluatePoint(dOneSample);

        expdata = dSamMeans[kk*nOutputs+ii];

        if (printLevel >= 4)
        {
          printOutTS(PL_INFO,
               "Experiment %4d (out of %d) : ",kk+1,dnSamples);
          for (ii2 = 0; ii2 < nInputs; ii2++)
            printOutTS(PL_INFO,
                 "Input %7d = %12.4e ",ii2+1,dOneSample[ii2]);
          printOutTS(PL_INFO,
               "simuation, experimental data = %12.4e %12.4e\n",
               simdata, expdata);
        }

        discOutputs[ii*dnSamples+kk] = expdata - simdata;
      }

      if (dnInputs > 0) 
      {
        for (kk = 0; kk < ExpNSamples*dnInputs; kk++)
          tSamInputs[kk] = dSamInputs[kk];
        iNames = new char*[dnInputs];
        cnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          if (designParams[ii2] == 1)
          {
            iNames[cnt] = new char[100];
            tLowers[cnt] = xLower[ii2];
            tUppers[cnt] = xUpper[ii2];
            if (qData.strArray_ == NULL)
                 sprintf(iNames[cnt], "X%d", ii2+1);
            else strcpy(iNames[cnt], qData.strArray_[ii2]);
            cnt++;
          }
        }
        dataPtr->updateInputSection(ExpNSamples,dnInputs,NULL,tLowers,
                       tUppers,dSamInputs, iNames, NULL,NULL,NULL,NULL);
        for (ii2 = 0; ii2 < dnInputs; ii2++) delete [] iNames[ii2];
        delete [] iNames;
      }
      else
      {
        iNames = new char*[1];
        iNames[0] = new char[100];
        sprintf(iNames[0], "X0");
        for (ii2 = 0; ii2 < ExpNSamples; ii2++) tSamInputs[ii2] = 0.5;
        tLowers[0] = 0.0;
        tUppers[0] = 1.0;
        dataPtr->updateInputSection(ExpNSamples,iOne,NULL,tLowers,tUppers,
                             tSamInputs, iNames, NULL,NULL,NULL,NULL);
        delete [] iNames[0];
        delete [] iNames;
      }

      states = new int[ExpNSamples];
      for (kk = 0; kk < ExpNSamples; kk++) states[kk] = 1;
      iNames = new char*[1];
      iNames[0] = new char[100];
      sprintf(iNames[0], "Y%d", ii+1);
      dataPtr->updateOutputSection(ExpNSamples,iOne,
                     &discOutputs[dnSamples*ii],states,iNames);
      delete [] states;
      delete [] iNames[0];
      delete [] iNames;
      dataPtr->updateMethodSection(PSUADE_SAMP_MC, ExpNSamples, 1, -1, -1);
      sprintf(charString, "psDiscrepancyModel%d", ii+1);
      dataPtr->writePsuadeFile(charString, 0);

      printOutTS(PL_INFO,
          "Creating discrepancy response surface for output %d\n",ii+1);
      faPtrs1[ii] = NULL;
      if (dnInputs > 0 && dnSamples > 1)
      {
        faPtrs1[ii] = genFA(dfaType,dnInputs,iOne,ExpNSamples);
        if (faPtrs1[ii] == NULL)
        {
          printOutTS(PL_ERROR,
               "MCMC ERROR: cannot create discrepancy func for output %d.\n",
               ii+1);
          return -1.0;
        }
      }
      if (faPtrs1[ii] != NULL)
      {
        faPtrs1[ii]->setNPtsPerDim(dnPerDim);
        faPtrs1[ii]->setBounds(xLower, xUpper);
        faPtrs1[ii]->setOutputLevel(0);
      }
      else
      {
        discFuncConstantMeans[ii] = 0.0;
        for (kk = 0; kk < ExpNSamples; kk++)
          discFuncConstantMeans[ii] += discOutputs[ii*dnSamples+kk];
        discFuncConstantMeans[ii] /= (double) ExpNSamples;
        discFuncConstantStds[ii] = 0.0;
        for (kk = 0; kk < ExpNSamples; kk++)
          discFuncConstantStds[ii] += pow(discOutputs[ii*dnSamples+kk]-
                                          discFuncConstantMeans[ii],2.0);
        discFuncConstantStds[ii] = 
                sqrt(discFuncConstantStds[ii]/ExpNSamples);
      }
    }
    delete dataPtr;
    delete [] ExpSamInputs;
    delete [] dOneSample;
    delete [] settings;
    delete [] tSamInputs;
    delete [] tLowers;
    delete [] tUppers;
  }

  //    set up constraint filters, if any
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"MCMC INFO: creating constraints, if there is any.\n");
  printOutTS(PL_INFO,
       "     Constraints remove infeasible regions from the priors.\n");
  printOutTS(PL_INFO,
       "     Constraints can be specified by RS constraint files.\n");
  constrPtr = new RSConstraints();
  constrPtr->genConstraints(dataPtr);
  printEquals(PL_INFO, 0);
  if (psAnaExpertMode_ == 1)
  {
    sprintf(charString, "How many MCMC chains? (1-20, default=3) : ");
    numChains = getInt(1,20,charString);
    sprintf(charString, "PSRF threshold? (1.0 - 1.2, default = 1.05) : ");
    psrfThreshold = 1.0;
    if (numChains > 1)
    {
       psrfThreshold = getDouble(charString);
       if (psrfThreshold < 1.0 || psrfThreshold > 1.2)
       {
         printOutTS(PL_INFO,
              "MCMC : invalid PSRF threshold ==> reset to 1.05.\n");
         psrfThreshold = 1.05;
      }
    }
  }

  //    set up for MCMC iterations
  int    *Ivec, **bins, ****bins2, globalIts, countTrack, dcnt;
  int    mcmcFail=0, sumBins, index2, nFail;
  int    ii3, jj2, kk2, index, length, iChain, chainCnt, mcmcIts;
  int    maxGlobalIts=20, chainCntSave, *chainStatus;
  double *XRange=NULL, *XGuess=NULL, *XDist=NULL, *XDesignS, *YDesignS;
  double *YDesignStds=NULL, *XGuessS=NULL, *YGuessS=NULL, *YGuessStds=NULL;
  double Xtemp, Ytemp, Ytemp2, *Xmax, Ymax, *s2Vec, *SDist;
  double ***XChains=NULL, stdev, stdev2, ddata2, WStat, BStat;
  double *chainMeans=NULL, *chainStdevs=NULL, *psrfs=NULL;
  TwoSampleAnalyzer *s2Analyzer=NULL;
  XRange  = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) XRange[ii] = xUpper[ii] - xLower[ii]; 
  XDist   = new double[maxPts+1];
  SDist   = new double[maxPts+1];
  XGuess  = new double[nInputs];
  XGuessS = new double[dnSamples*nInputs*(maxPts+1)];
  YGuessS = new double[dnSamples*nOutputs*(maxPts+1)];
  YGuessStds = new double[dnSamples*nOutputs*(maxPts+1)];
  XDesignS = new double[dnSamples*nInputs*(maxPts+1)];
  YDesignS = new double[dnSamples*nOutputs*(maxPts+1)];
  YDesignStds = new double[dnSamples*nOutputs*(maxPts+1)];
  means_ = new double[nInputs_];
  sigmas_ = new double[nInputs_];
  for (ii = 0; ii < nInputs; ii++) means_[ii] = sigmas_[ii] = 0.0;
  Xmax = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) Xmax[ii] = 0;
  mostLikelyInput_ = new double[nInputs_];
  for (ii = 0; ii < nInputs_; ii++) mostLikelyInput_[ii] = 0;
  mostLikelyOutput_ = new double[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++) mostLikelyOutput_[ii] = 0;
  Ymax = -PSUADE_UNDEFINED;
  Ivec = new int[nInputs];
  Ivec[nInputs-1] = -1;
  XChains = new double**[numChains];
  for (ii = 0; ii < numChains; ii++)
  {
    XChains[ii] = new double*[maxGlobalIts*maxSamples];
    for (jj = 0; jj < maxGlobalIts*maxSamples; jj++)
      XChains[ii][jj] = new double[nInputs+1];
  }
  chainMeans = new double[numChains];
  chainStdevs = new double[numChains];
  chainStatus  = new int[numChains];
  checkAllocate(chainStatus, "chainStatus in MCMC::analyze");
  for (ii = 0; ii < numChains; ii++) chainMeans[ii] = chainStdevs[ii] = 0.0;
  for (ii = 0; ii < numChains; ii++) chainStatus[ii] = 0;
  psrfs = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) psrfs[ii] = 0.0;
  bins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    bins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) bins[ii][jj] = 0;
  }
  bins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    bins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      bins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        bins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          bins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }
  s2Vec = new double[maxGlobalIts*maxSamples];
  checkAllocate(s2Vec, "s2Vec in MCMC::analyze");
  if (printLevel > 3) s2Analyzer = new TwoSampleAnalyzer();

  Sampling *sampler;
  if (nInputs > 50) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
  else              sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
  sampler->setInputBounds(nInputs, xLower, xUpper);
  sampler->setOutputParams(1);
  sampler->setSamplingParams(numChains, 1, 1);
  sampler->initialize(0);
  double *mcmcSeeds = new double[numChains*nInputs];
  double *tmpOuts = new double[numChains];
  int    *tmpStates = new int[numChains];
  checkAllocate(tmpStates, "tmpStates in MCMC::analyze");
  sampler->getSamples(numChains,nInputs,1,mcmcSeeds,tmpOuts,tmpStates);
  delete [] tmpOuts;
  delete [] tmpStates;
  delete sampler;
  for (iChain = 0; iChain < numChains; iChain++)
  {
    for (ii = 0; ii < nInputs; ii++)
    {
      ddata = mcmcSeeds[iChain*nInputs+ii];
      ddata = (ddata - xLower[ii]) / XRange[ii];
      mcmcSeeds[iChain*nInputs+ii] = ddata;
    }
  }
#if PS_INTERP == 1
  double b12, b11, b21, b22, det, aa, bb, cc, xd1, xd2;
  b12 = 1.0/maxPts; b11 = b12 * b12;
  b22 = 2.0/maxPts; b21 = b22 * b22;
  det = 1.0 / (b11 * b22 - b12 * b21);
#endif
#if PS_INTERP == 2
  FuncApprox *faDist;
  double *ZDist = new double[maxPts+1];
  for (jj = 0; jj <= maxPts; jj++) ZDist[jj] = 1.0 * jj / maxPts;
  faType = PSUADE_RS_MARS;
  faDist = genFA(faType, iOne, iZero, maxPts+1);
  faDist->setNPtsPerDim(16);
  double lo=0.0, hi=1.0;
  faDist->setBounds(&lo, &hi);
  faDist->setOutputLevel(-1);
  double xtrial, xtol=1.0e-4, xbeg, xend;
#endif
   
  // get ready to run the Metropolis Hasting algorithm
  int    nActive, *activePDFTypes, *activeList, simNum;
  double *activeMeans=NULL, *activeStdvs=NULL, *oneSample;
  double *likelyMeans=NULL, *likelyStdvs=NULL, dOne=1.0, *simOuts;
  PDFManager *pdfprior=NULL, *pdflikelihood=NULL;
  psVector vecSam, vecLower, vecUpper;
  psMatrix activeCorMat;
  nActive = 0;
  activeList  = new int[nInputs];
  activePDFTypes = new int[nInputs];
  activeMeans = new double[nInputs];
  activeStdvs = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++)
  {
    if ((rsIndices == NULL || rsIndices[ii] >= 0) &&
        (designParams == NULL || designParams[ii] == 0)) 
    {
      activeList[nActive] = ii;
      if (pdfTypes != NULL) activePDFTypes[nActive] = pdfTypes[ii];
      else                  activePDFTypes[nActive] = 0;
      if (pdfMeans != NULL) activeMeans[nActive] = pdfMeans[ii];
      else                  activeMeans[nActive] = 0;
      if (pdfStdvs != NULL) activeStdvs[nActive] = pdfStdvs[ii];
      else                  activeStdvs[nActive] = 1;
      nActive++;
    }
  }
  vecLower.setLength(nActive);
  vecUpper.setLength(nActive);
  for (ii = 0; ii < nActive; ii++)
  {
     vecLower[ii] = xLower[activeList[ii]];
     vecUpper[ii] = xUpper[activeList[ii]];
  }
  likelyMeans = new double[nInputs];
  likelyStdvs = new double[nInputs];
  activeCorMat.setDim(nActive,nActive);
  for (ii = 0; ii < nActive; ii++) activeCorMat.setEntry(ii,ii,dOne);
  vecSam.setLength(nActive);
  simOuts = new double[nOutputs];
  oneSample = new double[nInputs];
  simNum = 1;

  // run the Metropolis Hasting algorithm
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO, "MCMC begins ... \n");
  fflush(stdout);
  fp = NULL;
  globalIts = chainCnt = 0;
  while (globalIts < maxGlobalIts)
  {
    for (iChain = 0; iChain < numChains; iChain++)
    {
      printOutTS(PL_INFO,"MCMC : Chain %d, iteration = %d\n",iChain+1,
                 globalIts+1);
      if (iChain == 0) chainCntSave = chainCnt;
      else             chainCnt     = chainCntSave;
      if (chainCnt == 0)
      {
        for (ii = 0; ii < nInputs; ii++)
        {
          if (designParams == NULL || designParams[ii] == 0)
               XGuess[ii] = mcmcSeeds[iChain*nInputs+ii];
          else XGuess[ii] = 0.5;
          XGuess[ii] = XGuess[ii]*(xUpper[ii]-xLower[ii])+xLower[ii];
          XChains[iChain][0][ii] = XGuess[ii];
        }
        if (rsIndices != NULL)
        {
          for (ii = 0; ii < nInputs; ii++)
            if (rsIndices[ii] < 0) XGuess[ii] = rsValues[ii];
        }
        XChains[iChain][0][nInputs] = -1;
      }
      else
      {
        for (ii = 0; ii < nInputs; ii++)
          XGuess[ii] = XChains[iChain][chainCnt-1][ii];
      }
      if (printLevel >= 0)
      {
        printOutTS(PL_INFO,"       Chain %d current initial guess : \n",
                   iChain+1);
        for (ii = 0; ii < nInputs; ii++)
          printOutTS(PL_INFO,"          Input %4d = %e\n",
                     ii+1,XGuess[ii]);
      }
       
      mcmcIts = countTrack = 0;
      while (mcmcIts < maxSamples)
      {
        count = (mcmcIts+1) / (maxSamples/10);
        if (count != countTrack)
        { 
          countTrack++;
          printOutTS(PL_INFO, "%3.0f%% ",10.0*countTrack );
          fflush(stdout);
        }
        pdflikelihood = new PDFManager();
        for (ii = 0; ii < nActive; ii++)
          likelyMeans[ii] = XGuess[activeList[ii]];
        pdflikelihood->initialize(nActive,activePDFTypes,likelyMeans,
                                  likelyStdvs,activeCorMat,NULL,NULL);
        pdflikelihood->genSample(iOne, vecSam, vecLower, vecUpper);
        for (ii = 0; ii < nInputs; ii++) oneSample[ii] = XGuess[ii];
        for (ii = 0; ii < nActive; ii++)
          oneSample[activeList[ii]] = vecSam[ii];

        if (funcIO != NULL)
        {
          funcIO->evaluate(simNum,nInputs,oneSample,nOutputs,simOuts,0);
          simNum++;
        }
        else
        {
          for (ii = 0; ii < nOutputs; ii++)
            simOuts[ii] = faPtrs[ii]->evaluatePoint(oneSample);
        }
        /* compute acceptance */
        if (mcmcIts >= maxSamples) break;
      }
      if (countTrack <= 10) printOutTS(PL_INFO,"100%%\n");
      else                  printOutTS(PL_INFO,"\n");
      if (printLevel >= 0)
      {
        printOutTS(PL_INFO,
             "       Chain %d current final guess : \n",iChain+1);
        for (ii = 0; ii < nInputs; ii++)
          printOutTS(PL_INFO,"          Input %4d = %e\n",ii+1,XGuess[ii]);
      }
    }
 
    globalIts++;
    printOutTS(PL_INFO, "\nIteration %d summary: \n", globalIts);
    mcmcFail = nInputs - dnInputs;
    if (rsIndices != NULL)
    {
      for (ii = 0; ii < nInputs; ii++)
        if (rsIndices[ii] < 0) mcmcFail--; 
    }
    for (ii = 0; ii < nInputs; ii++)
    {
      if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
          (designParams == NULL || designParams[ii] == 0))
      {
        if (printLevel > 2) printOutTS(PL_INFO, "Input = %d\n", ii+1);
       
        for (iChain = 0; iChain < numChains; iChain++)
        {
          ddata = 0.0;
          for (jj = 0; jj < chainCnt; jj++) 
            ddata += XChains[iChain][jj][ii];
          ddata /= chainCnt;
          ddata2 = 0.0;
          for (jj = 0; jj < chainCnt; jj++) 
            ddata2 += pow(XChains[iChain][jj][ii]-ddata,2.0);
          ddata2 /= (double) (chainCnt - 1);
          chainMeans[iChain] = ddata;
          chainStdevs[iChain] = ddata2;
          if (globalIts > 2 && chainStdevs[iChain] < 1.0e-20) 
          {
            printOutTS(PL_INFO,
                 "MCMC INFO: chain %d disabled.\n",iChain+1);
            chainStatus[iChain] = 1;
          }
        }
        nChainGood = 0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0) nChainGood++;
        }
        if (nChainGood <= 1)
        {
          printOutTS(PL_ERROR,"MCMC ERROR: too few chains <= 1.\n");
          printOutTS(PL_ERROR,
             "Suggestion: You may want to relax the experimental data\n");
          printOutTS(PL_ERROR,"     uncertainties (make them larger).\n");
          printOutTS(PL_ERROR,
             "     To see if this is the problem, turn on printlevel\n");
          printOutTS(PL_ERROR,
             "     to 3 and run again. If the variance of the chains\n");
          printOutTS(PL_ERROR,
             "     are small, small data uncertainties is probably the\n");
          printOutTS(PL_ERROR,"     problem.\n");
          exit(1);
        }
        WStat = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0) WStat += chainStdevs[iChain];
        }
        WStat /= (double) nChainGood;
        if (WStat < 0) WStat = PSUADE_UNDEFINED;
        if (printLevel > 2) 
          printOutTS(PL_INFO,"  Within  chain variance W = %e\n", WStat);
        ddata = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0) ddata += chainMeans[iChain];
        }
        ddata /= (double) nChainGood;
        BStat = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0)
            BStat += pow(chainMeans[iChain]-ddata,2.0);
        }
        BStat = BStat / (nChainGood - 1.0) * chainCnt;
        if (printLevel > 2) 
          printOutTS(PL_INFO,
               "  Between chain variance B = %e\n", BStat/chainCnt);
        ddata = (1 - 1.0/chainCnt) * WStat + BStat / chainCnt;
        ddata = ddata / WStat * (numChains + 1) / numChains - 
                  (chainCnt - 1.0) / (double) (chainCnt * numChains); 
        if (ddata < 0) ddata2 = PSUADE_UNDEFINED;
        else           ddata2 = sqrt(ddata);
        if (printLevel > 2)
        {
          for (iChain = 0; iChain < numChains; iChain++)
            printOutTS(PL_INFO,"  Chain %d mean, var = %e %e\n",iChain+1,
                    chainMeans[iChain]*XRange[ii]+xLower[ii],
                    chainStdevs[iChain]*XRange[ii]*XRange[ii]);
          printOutTS(PL_INFO,"  Chain length             = %d\n",chainCnt);
          printOutTS(PL_INFO,"  Weighted average of B, W = %e\n", ddata);
        }
        printOutTS(PL_INFO,"  Input %d PSRF = %e\n", ii+1, ddata2);
        psrfs[ii] = ddata2;
        if (ddata2 < psrfThreshold)
        {
          printOutTS(PL_INFO,"MCMC INFO : PSRF < %e ==> converged.\n",
                     psrfThreshold);
          mcmcFail--;
        }
        ddata = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0)
            for (jj = 0; jj < chainCnt; jj++)
               ddata += XChains[iChain][jj][ii];
        }
        ddata /= (double) (nChainGood * chainCnt);
        means_[ii] = ddata;
        ddata2 = 0.0;
        for (iChain = 0; iChain < numChains; iChain++)
        {
          if (chainStatus[iChain] == 0)
            for (jj = 0; jj < chainCnt; jj++)
              ddata2 += pow(XChains[iChain][jj][ii]-ddata,2.0);
        }
        ddata2 /= (double) (chainCnt*nChainGood-1);
        sigmas_[ii] = sqrt(ddata2);
      }
    }
    if (mcmcFail == 0 && printLevel > 3 && s2Analyzer != NULL)
    {
      if (rsIndices != NULL)
      {
        for (ii = 0; ii < nInputs; ii++)
          if (rsIndices[ii] < 0) mcmcFail--; 
      }
      for (ii = 0; ii < nInputs; ii++)
      {
        if ((rsIndices == NULL || 
            (rsIndices != NULL && rsIndices[ii] >=0)) && 
            (designParams == NULL || designParams[ii] == 0))
        {
          printOutTS(PL_INFO, "Geweke Input = %d\n", ii+1);
     
          cnt = chainCnt / 2;
          ddata2 = (double) numChains;
          for (iChain = 0; iChain < numChains; iChain++)
          {
            for (jj = 0; jj < 2*cnt; jj++) 
              s2Vec[jj] = XChains[iChain][chainCnt-2*cnt+jj][ii];
            ddata = s2Analyzer->TAnalyze(cnt,s2Vec,cnt,&s2Vec[cnt],1);
           }
        }
      }
    }

    for (ii = 0; ii < nInputs; ii++) 
    {
      if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
          (designParams == NULL || designParams[ii] == 0))
      {
        printOutTS(PL_INFO,
             "MCMC: input %3d value at peak of likelihood = %e\n",
                   ii+1, Xmax[ii]);
        ddata = means_[ii]*(xUpper[ii]-xLower[ii])+xLower[ii];
        printOutTS(PL_INFO,"MCMC: input %3d mean    = %e\n", ii+1, ddata);
        ddata = sigmas_[ii]*(xUpper[ii]-xLower[ii]);
        printOutTS(PL_INFO,"MCMC: input %3d std dev = %e\n", ii+1, ddata);
        mostLikelyInput_[ii] = Xmax[ii];
      }
    }
    if (mcmcFail == 0) break;

    genMatlabFile(nInputs,xLower,xUpper,XRange,nPlots,plotIndices,nbins,
             NULL,NULL,bins,bins2,qData,numChains,chainCnt,XChains,
             chainStatus,Xmax,0);

    fp = fopen("psuade_stop", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_stop FILE FOUND - TERMINATE MCMC.\n");
      fclose(fp);
      fp = NULL;
      strcpy(charString, "psuade_stop");
      unlink(charString);
      break;
    }
    fp = fopen("psuade_nogm", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_nogm FILE FOUND. GM mode is now off.\n");
      fclose(fp);
      fp = NULL;
      psGMMode_ = 0;
      strcpy(charString, "psuade_nogm");
      unlink(charString);
    }
    fp = fopen("psuade_gm", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_gm FILE FOUND. GM mode is now on.\n");
      fclose(fp);
      fp = NULL;
      psGMMode_ = 1;
      strcpy(charString, "psuade_gm");
      unlink(charString);
    }
    fp = fopen("psuade_print", "r");
    if (fp != NULL)
    {
      printOutTS(PL_INFO,
           "MCMC INFO: psuade_print FILE FOUND. Print level is set to 3.\n");
      fclose(fp);
      fp = NULL;
      printLevel = 3;
      strcpy(charString, "psuade_print");
      unlink(charString);
    }
  }
  if (globalIts >= maxGlobalIts)
  {
    mcmcFail = 0;
    for (ii = 0; ii < nInputs; ii++) 
    {
      if ((rsIndices == NULL || (rsIndices != NULL && rsIndices[ii] >=0)) && 
          (designParams == NULL || designParams[ii] == 0))
         if (psrfs[ii] > psrfThreshold) mcmcFail = 1;
    }
    if (mcmcFail == 1) 
      printOutTS(PL_INFO,
           "MCMC maximum iterations exceeded but no convergence.\n");
  }
  else printOutTS(PL_INFO, "MCMC iterations completed\n");

  for (ii = 0; ii < nInputs; ii++) 
    for (jj = 0; jj < nbins; jj++) bins[jj][ii] = 0;
  for (ii = 0; ii < nInputs; ii++) 
    for (ii2 = 0; ii2 < nInputs; ii2++) 
      for (jj = 0; jj < nbins; jj++)
        for (jj2 = 0; jj2 < nbins; jj2++)
          bins2[jj][jj2][ii][ii2] = 0;
  for (iChain = 0; iChain < numChains; iChain++) 
  { 
    if (chainStatus[iChain] == 0)
    {
      for (jj = 0; jj < chainCnt; jj++) 
      {
        for (ii2 = 0; ii2 < nInputs; ii2++) 
        {
          ddata = XChains[iChain][jj][ii2];
          ddata = (ddata - xLower[ii2]) / (xUpper[ii2] - xLower[ii2]);
          index = (int) (ddata * nbins);
          if (index > nbins)
            printOutTS(PL_ERROR,"MCMC binning error 1 in file %s, line %d.\n",
                       __FILE__, __LINE__);
          if (index < 0)
            printOutTS(PL_ERROR,"MCMC binning error 2 in file %s, line %d.\n",
                       __FILE__, __LINE__);
          if (index >= nbins) index = nbins - 1;
          if (index <  0)     index = 0;
          bins[index][ii2]++;
        }
        for (ii2 = 0; ii2 < nInputs; ii2++) 
        {
          ddata = XChains[iChain][jj][ii2];
          ddata = (ddata - xLower[ii2]) / (xUpper[ii2] - xLower[ii2]);
          index = (int) (ddata * nbins);
          if (index >= nbins) index = nbins - 1;
          if (index <  0)     index = 0;
          for (ii3 = 0; ii3 < nInputs; ii3++) 
          {
            ddata2 = XChains[iChain][jj][ii3];
            ddata2 = (ddata2 - xLower[ii3]) / (xUpper[ii3] - xLower[ii3]);
            index2 = (int) (ddata2 * nbins);
            if (index2 >= nbins) index2 = nbins - 1;
            if (index2 <  0)     index2 = 0;
            bins2[index][index2][ii2][ii3]++;
          }
        }
      }
    }
  }

  delete [] XDist;
  delete [] SDist;
  delete [] XGuess;
  delete [] XGuessS;
  delete [] YGuessS;
  delete [] YGuessStds;
  delete [] XDesignS;
  delete [] YDesignS;
  delete [] YDesignStds;
 
  dataPtr->getParameter("method_sampling", pPtr);
  int methodSave = pPtr.intData_;
  dataPtr->updateMethodSection(PSUADE_SAMP_MC,-1,-1,-1,-1);
  PDFManager *pdfman = new PDFManager();
  pdfman->initialize(dataPtr);
  psVector vecLB, vecUB, vecOut;
  vecLB.load(nInputs, xLower);
  vecUB.load(nInputs, xUpper);
  int nSamps = 100000;
  vecOut.setLength(nSamps*nInputs);
  pdfman->genSample(nSamps, vecOut, vecLB, vecUB);
  dataPtr->updateMethodSection(methodSave,-1,-1,-1,-1);
  delete pdfman;
  int **pbins = new int*[nbins];
  for (ii = 0; ii < nbins; ii++)
  {
    pbins[ii] = new int[nInputs];
    for (jj = 0; jj < nInputs; jj++) pbins[ii][jj] = 0;
  }
  int ****pbins2 = new int***[nbins];
  for (jj = 0; jj < nbins; jj++)
  {
    pbins2[jj] = new int**[nbins];
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      pbins2[jj][jj2] = new int*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
        pbins2[jj][jj2][ii] = new int[nInputs];
        for (ii2 = 0; ii2 < nInputs; ii2++)
          pbins2[jj][jj2][ii][ii2] = 0;
      }
    }
  }
  double ddata1;
  for (ii = 0; ii < nInputs; ii++)
  {
    ddata1 = nbins / (xUpper[ii] - xLower[ii]);
    for (jj = 0; jj < nSamps; jj++)
    {
      ddata = vecOut[jj*nInputs+ii];
      ddata -= xLower[ii];
      ddata *= ddata1;
      kk = (int) ddata;
      if (kk >= nbins) kk = nbins - 1; 
      pbins[kk][ii]++;
    }
  }
  for (ii = 0; ii < nInputs; ii++)
  {
    ddata1 = nbins / (xUpper[ii] - xLower[ii]);
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      ddata2 = nbins / (xUpper[ii2] - xLower[ii2]);
      for (jj = 0; jj < nSamps; jj++)
      {
        ddata = vecOut[jj*nInputs+ii];
        ddata -= xLower[ii];
        ddata *= ddata1;
        kk = (int) ddata;
        if (kk >= nbins) kk = nbins - 1; 
        ddata = vecOut[jj*nInputs+ii2];
        ddata -= xLower[ii2];
        ddata *= ddata2;
        kk2 = (int) ddata;
        if (kk2 >= nbins) kk2 = nbins - 1; 
        pbins2[kk][kk2][ii][ii2]++;
      }
    }
  }
  genMatlabFile(nInputs,xLower,xUpper,XRange,nPlots,plotIndices,nbins,
        NULL,NULL,bins,bins2,qData,numChains,chainCnt,XChains,
        chainStatus, Xmax,0);
  for (ii = 0; ii < nbins; ii++) delete [] pbins[ii];
  delete [] pbins;
  for (jj = 0; jj < nbins; jj++)
  {
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      for (ii = 0; ii < nInputs; ii++) delete [] pbins2[jj][jj2][ii];
      delete [] pbins2[jj][jj2];
    }
    delete [] pbins2[jj];
  }
  delete [] pbins2;

  if (genPosteriors == 1)
  {
    cnt = nChainGood * chainCnt;
    if (cnt > 200000) cnt = 200000;
    cnt /= nChainGood;
    genPostLikelihood(nInputs,xLower,xUpper,XRange,numChains, chainCnt,
                  XChains, chainStatus, cnt, rsIndices, rsValues,
                  designParams, dnInputs, dnSamples, dSamInputs, faPtrs,
                  faPtrs1, nOutputs, discOutputs, discFuncConstantMeans,
                  dSamMeans, dSamStdevs);
  }

  if (genPosteriors == 1)
  {
    fp = fopen("MCMCPostSample", "w");
    if (fp != NULL)
    {
      fprintf(fp, "PSUADE_BEGIN\n");
      cnt = nChainGood * chainCnt;
      if (cnt > 50000) cnt = 50000;
      cnt /= nChainGood;
      fprintf(fp, "%d %d\n", cnt*nChainGood,nInputs);
      if (qData.strArray_ != NULL)
      {
        fprintf(fp, "# ");
        for (jj = 0; jj < nInputs; jj++)
          fprintf(fp,"%s ", qData.strArray_[jj]);
        fprintf(fp, "\n");
      }
      ii2 = 0;
      for (iChain = 0; iChain < numChains; iChain++)
      { 
        if (chainStatus[iChain] == 0)
        {
          for (ii = chainCnt-cnt; ii < chainCnt; ii++)
          {
            fprintf(fp, "%d ", ii2+1);
            for (jj = 0; jj < nInputs; jj++)
            {
              if ((rsIndices == NULL || rsIndices[jj] >= 0) &&
                  (designParams == NULL || designParams[jj] == 0)) 
              {
                ddata = XChains[iChain][ii][jj];
                fprintf(fp, "%e ", ddata);
              }
              else if (rsIndices != NULL && rsIndices[jj] < 0)
                fprintf(fp, "%e ", rsValues[jj]);
              else if (designParams != NULL && designParams[jj] != 0) 
                fprintf(fp, "%e ", 0.5 * (xUpper[jj] + xLower[jj]));
            }
            fprintf(fp, "\n");
            ii2++;
          }
        }
      }
      fprintf(fp, "PSUADE_END\n");
      fprintf(fp, "#N=%d;\n",nChainGood*cnt);
      fprintf(fp, "#m=%d;\n",cnt);
      for (iChain = 0; iChain < numChains; iChain++)
      {
        if (chainStatus[iChain] == 0)
          fprintf(fp, "#A%d = A(%d*m+1:%d*m,:);\n",iChain+1,iChain,iChain+1);
      }
      fprintf(fp, "#for ii = 2 : %d\n", nInputs+1);
      ii2 = sqrt(1.0*nInputs);
      if (ii2 * ii2 < nInputs) ii2++;
      for (iChain = 0; iChain < numChains; iChain++)
      {
        if (chainStatus[iChain] == 0)
        {
          fprintf(fp, "#subplot(%d,%d,%d)\n",ii2,ii2,iChain+1);
          fprintf(fp, "#hist(A%d(:,ii))\n",iChain+1);
        }
        fprintf(fp, "#ii-1\n");
        fprintf(fp, "#pause;\n");
      }
      fprintf(fp, "#end;\n");
      fclose(fp);
    }
    printOutTS(PL_INFO,
         "MCMC: 'MCMCPostSample' file has a posterior sample.\n");
  }

  int    nInps, nOuts, nSams, *states;
  double *allOuts;
  char   **oNames;
  PsuadeData *filePtr1, *filePtr2;
  if (modelFormFlag == 1)
  {
    sprintf(charString, "psDiscrepancyModel1");
    filePtr1 = new PsuadeData();
    status = filePtr1->readPsuadeFile(charString);
    if (status != 0)
    {
      printOutTS(PL_ERROR,
           "MCMC ERROR: cannot read file %s in PSUADE format.\n",charString);
      exit(1);
    } 
  }
  if (modelFormFlag == 1 && status == 0)
  {
    filePtr1->getParameter("input_ninputs", pPtr);
    nInps = pPtr.intData_;
    filePtr1->getParameter("output_noutputs", pPtr);
    nOuts = pPtr.intData_;
    filePtr1->getParameter("method_nsamples", pPtr);
    nSams = pPtr.intData_;
    filePtr1->getParameter("output_sample", pOutputs);
    unlink(charString);
    allOuts = new double[nOutputs * nSams];
    for (jj = 0; jj < nSams; jj++) 
      allOuts[jj*nOutputs] = pOutputs.dbleArray_[jj];
    pOutputs.clean();
    oNames = new char*[nOutputs];
    oNames[0] = new char[100];
    sprintf(oNames[0], "Y1");
    for (ii = 1; ii < nOutputs; ii++)
    {
      filePtr2 = new PsuadeData();
      sprintf(charString, "psDiscrepancyModel%d", ii+1);
      status = filePtr2->readPsuadeFile(charString);
      if (status != 0) break;
      filePtr2->getParameter("input_ninputs", pPtr);
      if (pPtr.intData_ != nInps) break;
      filePtr2->getParameter("output_noutputs", pPtr);
      if (pPtr.intData_ != nOuts) break;
      filePtr2->getParameter("method_nsamples", pPtr);
      if (pPtr.intData_ != nSams) break;
      filePtr2->getParameter("output_sample", pOutputs);
      delete filePtr2;
      unlink(charString);
      for (jj = 0; jj < nSams; jj++) 
        allOuts[jj*nOutputs+ii] = pOutputs.dbleArray_[jj];
      pOutputs.clean();
      oNames[ii] = new char[100];
      sprintf(oNames[ii], "Y%d", ii+1);
    }
    if (nOutputs == 1)
    { 
      sprintf(charString, "psDiscrepancyModel");
      filePtr1->writePsuadeFile(charString, 0);
    }
    else if (ii == nOutputs)
    {
      states = new int[nSams];
      for (jj = 0; jj < nSams; jj++) states[jj] = 1;
      filePtr1->updateOutputSection(nSams,nOutputs,allOuts,states,oNames);
      sprintf(charString, "psDiscrepancyModel");
      filePtr1->writePsuadeFile(charString, 0);
      printOutTS(PL_INFO,
           "MCMC INFO: a sample (inputs/outputs) the discrepancy model\n");
      printOutTS(PL_INFO,"           is now in psDiscrepancyModel.\n");
      delete [] states;
      for (ii = 1; ii < nOutputs; ii++) delete [] oNames[ii];
    }
    else
    {
      printOutTS(PL_INFO,
           "MCMC INFO: unsuccessful creation of discrepancy sample file\n");
    }
    delete [] oNames[0];
    delete [] oNames;
    delete filePtr1;
    delete [] allOuts;
  }
  
  // clean up
  for (ii = 0; ii < numChains; ii++)
  {
    for (jj = 0; jj < maxGlobalIts*maxSamples; jj++)
      delete [] XChains[ii][jj];
    delete [] XChains[ii];
  }
  delete [] XChains;
  if (inputPDFs != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
      if (inputPDFs[ii] != NULL) delete inputPDFs[ii];
    delete [] inputPDFs;
  }
  delete [] Xmax;
  delete [] XRange;
  delete [] psrfs;
  if (discOutputs != NULL) delete [] discOutputs;
  delete [] s2Vec;
  delete [] mcmcSeeds;
  if (s2Analyzer != NULL) delete s2Analyzer;
  if (discFuncConstantMeans != NULL) delete [] discFuncConstantMeans;
  if (discFuncConstantStds  != NULL) delete [] discFuncConstantStds;
  if (plotIndices != NULL) delete [] plotIndices;
  if (dSamMeans != NULL) delete [] dSamMeans;
  if (dSamStdevs != NULL) delete [] dSamStdevs;
#if 0
  delete [] ZDist;
  delete faDist;
#endif
  for (ii = 0; ii < nbins; ii++) delete [] bins[ii];
  delete [] bins;
  for (jj = 0; jj < nbins; jj++)
  {
    for (jj2 = 0; jj2 < nbins; jj2++)
    {
      for (ii = 0; ii < nInputs; ii++) delete [] bins2[jj][jj2][ii];
      delete [] bins2[jj][jj2];
    }
    delete [] bins2[jj];
  }
  delete [] bins2;
  if (dSamInputs != NULL) delete [] dSamInputs;
  if (rsIndices != NULL) delete [] rsIndices;
  if (rsValues  != NULL) delete [] rsValues;
  delete constrPtr;
  if (faPtrs != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++) 
      if (faPtrs[ii] != NULL) delete faPtrs[ii];
    delete [] faPtrs;
  }
  if (faPtrs1 != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++) 
      if (faPtrs1[ii] != NULL) delete faPtrs1[ii];
    delete [] faPtrs1;
  }
  delete [] Ivec;
  return 0.0;
}

// ************************************************************************
// write to matlab file 
// ------------------------------------------------------------------------
double MCMCAnalyzer::genMatlabFile(int nInputs, double *lower, double *upper,
                                   double *XRange, int nPlots, int *plotIndices,
                                   int nbins, int **pbins, int ****pbins2, 
                                   int **bins, int ****bins2, 
                                   pData &qData, int nChains, int chainCnt,
                                   double ***XChains, int *chainStatus,
                                   double *Xmax, double Ymin)
{
  int    kk, kk2, ii2, jj, jj2, sumBins, iChain;
  double ddata, dmean, dstd;
  char   cfname[1001], charString[1001];;
  FILE   *fp;

  if (psPlotTool_ == 1) strcpy(cfname, "scilabmcmc2.sci");
  else                  strcpy(cfname, "matlabmcmc2.m");
  fp = fopen(cfname, "w");
  if (fp == NULL)
  {
    printOutTS(PL_ERROR, "ERROR: cannot open %s file.\n", cfname);
    return 0;
  }
  sprintf(charString,"This file shows posteriors plots");
  fwriteComment(fp, charString);
  sprintf(charString,"ns  - set to 1 for 1-step smoothing of 2D contours");
  fwriteComment(fp, charString);
  sprintf(charString,"ns1 - set to 1 for 1-step smoothing of 1D histgrams");
  fwriteComment(fp, charString);
  fprintf(fp, "ns  = 0;\n");
  fprintf(fp, "ns1 = 0;\n");
  fwritePlotCLF(fp);
  fprintf(fp, "active = [\n");
  for (kk = 0; kk < nInputs; kk++)
  {
    ii2 = binarySearchInt(kk, plotIndices, nPlots);
    if (ii2 < 0) fprintf(fp, "0\n");
    else         fprintf(fp, "1\n");
  }
  fprintf(fp, "];\n");
  fprintf(fp, "L = [\n");
  for (kk = 0; kk < nInputs; kk++) fprintf(fp, "%e ",lower[kk]);
  fprintf(fp, "];\n");
  fprintf(fp, "U = [\n");
  for (kk = 0; kk < nInputs; kk++) fprintf(fp, "%e ",upper[kk]);
  fprintf(fp, "];\n");
  fprintf(fp, "iStr = {\n");
  for (kk = 0; kk < nInputs-1; kk++)
  {
    if (qData.strArray_ != NULL)
         fprintf(fp, "'%s',", qData.strArray_[kk]);
    else fprintf(fp, "'Input %d',", kk+1);
  }
  if (qData.strArray_ != NULL)
       fprintf(fp, "'%s'};\n", qData.strArray_[nInputs-1]);
  else fprintf(fp, "'Input %d'};\n", nInputs);
  fprintf(fp, "X = zeros(%d,%d);\n", nInputs, nbins);
  fprintf(fp, "D = zeros(%d,%d);\n", nInputs, nbins);
  fprintf(fp, "NC = zeros(%d,%d,%d,%d);\n",nInputs,nInputs,nbins,nbins);
  for (kk = 0; kk < nInputs; kk++)
  {
    for (kk2 = 0; kk2 < nInputs; kk2++)
    {
      if (kk == kk2)
      {
        fprintf(fp, "X(%d,:) = [\n", kk+1);
        for (jj = 0; jj < nbins; jj++)
          fprintf(fp, "%e ", XRange[kk]/nbins*(jj+0.5)+lower[kk]);
        fprintf(fp, "];\n");
        sumBins = 0;
        for (jj = 0; jj < nbins; jj++) sumBins += bins[jj][kk];
        if (sumBins == 0) sumBins = 1;
        fprintf(fp, "D(%d,:) = [\n", kk+1);
        for (jj = 0; jj < nbins; jj++)
          fprintf(fp, "%e ", (double) bins[jj][kk]/(double) sumBins);
        fprintf(fp, "];\n");
      }
      else
      {
        fprintf(fp, "NC(%d,%d,:,:) = [\n", kk+1, kk2+1);
        for (jj = 0; jj < nbins; jj++)
        {
          for (jj2 = 0; jj2 < nbins; jj2++)
            fprintf(fp, "%d ", bins2[jj][jj2][kk][kk2]);
          fprintf(fp, "\n");
        }
        fprintf(fp, "]';\n");
      }
    }
  }
  if (pbins != NULL)
    fprintf(fp, "DP = zeros(%d,%d);\n", nInputs, nbins);
  if (pbins2 != NULL)
    fprintf(fp, "NCP = zeros(%d,%d,%d,%d);\n",nInputs,nInputs,nbins,nbins);
  for (kk = 0; kk < nInputs; kk++)
  {
    for (kk2 = 0; kk2 < nInputs; kk2++)
    {
      if (kk == kk2)
      {
        if (pbins != NULL)
        {
          sumBins = 0;
          for (jj = 0; jj < nbins; jj++) sumBins += pbins[jj][kk];
          fprintf(fp, "DP(%d,:) = [\n", kk+1);
          for (jj = 0; jj < nbins; jj++)
            fprintf(fp, "%e ", (double) pbins[jj][kk]/(double) sumBins);
          fprintf(fp, "];\n");
        }
      }
      else
      {
        if (pbins2 != NULL)
        {
          fprintf(fp, "NCP(%d,%d,:,:) = [\n", kk+1, kk2+1);
          for (jj = 0; jj < nbins; jj++)
          {
            for (jj2 = 0; jj2 < nbins; jj2++)
              fprintf(fp, "%d ", pbins2[jj][jj2][kk][kk2]);
            fprintf(fp, "\n");
          }
          fprintf(fp, "]';\n");
        }
      }
    }
  }
  fprintf(fp,"nInps  = length(active);\n");
  fprintf(fp,"nPlots = 0;\n");
  fprintf(fp,"for ii = 1 : nInps\n");
  fprintf(fp,"   if (active(ii) == 1)\n");
  fprintf(fp,"      nPlots = nPlots + 1;\n");
  fprintf(fp,"      active(ii) = nPlots;\n");
  fprintf(fp,"   end;\n");
  fprintf(fp,"end;\n");
  fprintf(fp,"dzero = 0;\n");
  if (Xmax != NULL)
  {
    for (kk = 0; kk < nInputs; kk++)
      fprintf(fp, "XOpt(%d) = %e;\n", kk+1, Xmax[kk]);
  }
  fprintf(fp,"for ii = 1 : nInps\n");
  fprintf(fp,"  for jj = ii : nInps\n");
  fprintf(fp,"    if (active(ii) ~= 0 & active(jj) ~= 0)\n");
  fprintf(fp,"      index = (active(ii)-1) * nPlots + active(jj);\n");
  fprintf(fp,"      subplot(nPlots,nPlots,index)\n");
  fprintf(fp,"      if (ii == jj)\n");
  fprintf(fp,"        n = length(D(ii,:));\n");
  fprintf(fp,"        DN = D(ii,:);\n");
  fprintf(fp,"        for kk = 1 : ns1\n");
  fprintf(fp,"          DN1 = DN;\n");
  fprintf(fp,"          for ll = 2 : n-1\n");
  fprintf(fp,"            DN(ll) = DN(ll) + DN1(ll+1);\n");
  fprintf(fp,"            DN(ll) = DN(ll) + DN1(ll-1);\n");
  fprintf(fp,"            DN(ll) = DN(ll) / 3;\n");
  fprintf(fp,"          end;\n");
  fprintf(fp,"        end;\n");
  fprintf(fp,"        bar(X(ii,:), DN, 1.0);\n");
  if (psPlotTool_ == 1)
       fprintf(fp,"        set(gca(),\"auto_clear\",\"off\")\n");
  else fprintf(fp,"        hold on\n");
  if (Xmax != NULL)
  {
    fprintf(fp,
      "        plot([XOpt(ii) XOpt(ii)],[0 max(DN)],'r-','LineWidth',2)\n");
  }
  if (pbins != NULL)
    fprintf(fp,"        plot(X(ii,:),DP(ii,:),'c-','LineWidth',2);\n");
  fprintf(fp,"        xmin = min(X(ii,:));\n");
  fprintf(fp,"        xmax = max(X(ii,:));\n");
  fprintf(fp,"        xwid = xmax - xmin;\n");
  fprintf(fp,"        xmin = xmin - 0.5 * xwid / %d;\n", nbins);
  fprintf(fp,"        xmax = xmax + 0.5 * xwid / %d;\n", nbins);
  fprintf(fp,"        ymax = max(DN);\n");
  if (psPlotTool_ == 1)
  {
    fprintf(fp,"        e = gce();\n");
    fprintf(fp,"        e.children.thickness = 2;\n");
    fprintf(fp,"        e.children.foreground = 0;\n");
    fprintf(fp,"        e.children.background = 2;\n");
    fprintf(fp,"        a = gca();\n");
    fprintf(fp,"        a.data_bounds=[xmin,0;xmax,ymax];\n");
    fprintf(fp,"        a.x_label.text = iStr(ii);\n");
    fprintf(fp,"        a.x_label.font_size = 3;\n");
    fprintf(fp,"        a.x_label.font_style = 4;\n");
    fprintf(fp,"        a.grid = [1 1];\n");
    fprintf(fp,"        a.y_label.text = iStr(jj);\n");
    fprintf(fp,"        a.y_label.font_size = 3;\n");
    fprintf(fp,"        a.y_label.font_style = 4;\n");
    fprintf(fp,"        a.thickness = 2;\n");
    fprintf(fp,"        a.font_size = 3;\n");
    fprintf(fp,"        a.font_style = 4;\n");
    fprintf(fp,"        a.box = \"on\";\n");
  }
  else
  {
    fprintf(fp,"        axis([xmin xmax 0 ymax])\n");
    fprintf(fp,"        set(gca,'linewidth',2)\n");
    fprintf(fp,"        set(gca,'fontweight','bold')\n");
    fprintf(fp,"        set(gca,'fontsize',12)\n");
    fprintf(fp, 
       "        xlabel(iStr(ii),'FontWeight','bold','FontSize',12)\n");
    fprintf(fp,"        if (ii == 1)\n"); 
    fprintf(fp,
       "        ylabel('Probabilities','FontWeight','bold','FontSize',12)\n");
    fprintf(fp,"        end;\n"); 
    fprintf(fp,"        grid on\n");
    fprintf(fp,"        box on\n");
  }
  fprintf(fp,"      else\n");
  fprintf(fp,"        n = length(X(jj,:));\n");
  fprintf(fp,"        XT = X(jj,:);\n");
  fprintf(fp,"        YT = X(ii,:);\n");
  fprintf(fp,"        HX = (XT(n) - XT(1)) / (n-1);\n");
  fprintf(fp,"        HY = (YT(n) - YT(1)) / (n-1);\n");
  fprintf(fp,"        ZZ = squeeze(NC(ii,jj,:,:));\n");
  fprintf(fp,"        for kk = 1 : ns\n");
  fprintf(fp,"          ZZ1 = ZZ;\n");
  fprintf(fp,"          for ll = 2 : n-1\n");
  fprintf(fp,"            for mm = 2 : n-1\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm+1);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm-1);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm+1);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm-1);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm+1);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm-1);\n");
  fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) / 9;\n");
  fprintf(fp,"            end;\n");
  fprintf(fp,"          end;\n");
  fprintf(fp,"        end;\n");
  fprintf(fp,"        ZZ = ZZ / (sum(sum(ZZ)));\n");
  if (psPlotTool_ == 1)
  {
    fprintf(fp,"        XX = [XT(1):HX:XT(n)];\n");
    fprintf(fp,"        YY = [YT(1):HY:YT(n)];\n");
    fprintf(fp,"        DD = splin2d(XX,YY,ZZ);\n");
    fprintf(fp,"        HX = 0.01 * (XT(n) - XT(1));\n");
    fprintf(fp,"        HY = 0.01 * (YT(n) - YT(1));\n");
    fprintf(fp,"        X2 = [XT(1):HX:XT(n)];\n");
    fprintf(fp,"        Y2 = [YT(1):HY:YT(n)];\n");
    fprintf(fp,"        [XI, YI] = ndgrid(X2, Y2);\n");
    fprintf(fp,"        disp('interpolation')\n");
    fprintf(fp,"        ZI =interp2d(XI, YI, XX, YY, DD, \"natural\");\n");
    fprintf(fp,"        disp('interpolation done')\n");
    fprintf(fp,"        ZB = ZI;\n");
    fprintf(fp,"        nX = length(X2);\n");
    fprintf(fp,"        nY = length(Y2);\n");
    fprintf(fp,"        for kk = 1 : nX\n");
    fprintf(fp,"          for ll = 1 : nY\n");
    fprintf(fp,"            ZI(kk,ll) = ZB(kk,nY-ll+1);\n");
    fprintf(fp,"          end;\n");
    fprintf(fp,"        end;\n");
    fprintf(fp,"        zmax = max(max(ZI));\n");
    fprintf(fp,"        zmin = min(min(ZI)) / zmax;\n");
    fprintf(fp,"        ZI   = ZI / zmax;\n");
    fprintf(fp,"        zmax = 1;\n");
    fprintf(fp,"        Matplot1((ZI-zmin)/(zmax-zmin)*64,[L(jj),L(ii),");
    fprintf(fp,"U(jj),U(ii)]);\n");
    fprintf(fp,"        xset(\"colormap\",jetcolormap(64));\n");
    fprintf(fp,"%%      colorbar(zmin,zmax);\n");
    fprintf(fp,
         "        contour2d(X2,Y2,ZB,5,rect=[L(jj),L(ii),U(jj),U(ii)]);\n");
    fprintf(fp,"        xset(\"fpf\",\" \");\n");
    fprintf(fp,"        a = gca();\n");
    fprintf(fp,"        a.x_label.text = iStr(jj);\n");
    fprintf(fp,"        a.x_label.font_size = 3;\n");
    fprintf(fp,"        a.x_label.font_style = 4;\n");
    fprintf(fp,"        a.y_label.text = iStr(ii);\n");
    fprintf(fp,"        a.y_label.font_size = 3;\n");
    fprintf(fp,"        a.y_label.font_style = 4;\n");
    fwritePlotAxesNoGrid(fp);
  }
  else
  {
    fprintf(fp,"%%      [YY,XX]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));\n");
    fprintf(fp,"%%      HX = 0.01 * (XT(n) - XT(1));\n");
    fprintf(fp,"%%      HY = 0.01 * (YT(n) - YT(1));\n");
    fprintf(fp,"%%      [YI,XI]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));\n");
    fprintf(fp,"%%      ZI=interp2(YY, XX, ZZ, YI, XI, 'spline');\n");
    fprintf(fp,"%%      pcolor(XI,YI,ZI)\n");
    fprintf(fp,"%%      shading interp\n");
    fprintf(fp,"%%      hold on\n");
    fprintf(fp,"%%      contour(XI,YI,ZI,5,'k')\n");
    fprintf(fp,"        imagesc(ZZ')\n");
    fprintf(fp,"        xtick = L(ii):(U(ii)-L(ii))/4:U(ii);\n");
    fprintf(fp,"        set(gca,'XTick',0:n/4:n);\n");
    fprintf(fp,"        set(gca,'XTickLabel', xtick);\n");
    fprintf(fp,"        ytick = L(jj):(U(jj)-L(jj))/4:U(jj);\n");
    fprintf(fp,"        set(gca,'YTick',0:n/4:n);\n");
    fprintf(fp,"        set(gca,'YTickLabel', ytick);\n");
    fprintf(fp,"        set(gca,'YDir', 'normal');\n");
    fprintf(fp,
       "        xlabel(iStr(jj),'FontWeight','bold','FontSize',12)\n");
    fprintf(fp,
       "        ylabel(iStr(ii),'FontWeight','bold','FontSize',12)\n");
    fwritePlotAxesNoGrid(fp);
    fprintf(fp,"%%      colorbar\n");
  }
  if (pbins2 != NULL)
  {
    fprintf(fp,"        index = (active(jj)-1) * nPlots + active(ii);\n");
    fprintf(fp,"        subplot(nPlots,nPlots,index)\n");
    fprintf(fp,"        ZZ = squeeze(NCP(ii,jj,:,:));\n");
    fprintf(fp,"        for kk = 1 : ns\n");
    fprintf(fp,"          ZZ1 = ZZ;\n");
    fprintf(fp,"          for ll = 2 : n-1\n");
    fprintf(fp,"            for mm = 2 : n-1\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm+1);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm-1);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm+1);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm-1);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm+1);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm-1);\n");
    fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) / 9;\n");
    fprintf(fp,"            end;\n");
    fprintf(fp,"          end;\n");
    fprintf(fp,"        end;\n");
    fprintf(fp,"        ZZ = ZZ / (sum(sum(ZZ)));\n");
    if (psPlotTool_ == 1)
    {
      fprintf(fp,"        XX = [XT(1):HX:XT(n)];\n");
      fprintf(fp,"        YY = [YT(1):HY:YT(n)];\n");
      fprintf(fp,"        DD = splin2d(XX,YY,ZZ);\n");
      fprintf(fp,"        HX = 0.01 * (XT(n) - XT(1));\n");
      fprintf(fp,"        HY = 0.01 * (YT(n) - YT(1));\n");
      fprintf(fp,"        X2 = [XT(1):HX:XT(n)];\n");
      fprintf(fp,"        Y2 = [YT(1):HY:YT(n)];\n");
      fprintf(fp,"        [XI, YI] = ndgrid(X2, Y2);\n");
      fprintf(fp,"        disp('interpolation')\n");
      fprintf(fp,"        ZI =interp2d(XI, YI, XX, YY, DD, \"natural\");\n");
      fprintf(fp,"        disp('interpolation done')\n");
      fprintf(fp,"        ZB = ZI;\n");
      fprintf(fp,"        nX = length(X2);\n");
      fprintf(fp,"        nY = length(Y2);\n");
      fprintf(fp,"        for kk = 1 : nX\n");
      fprintf(fp,"          for ll = 1 : nY\n");
      fprintf(fp,"            ZI(kk,ll) = ZB(kk,nY-ll+1);\n");
      fprintf(fp,"          end;\n");
      fprintf(fp,"        end;\n");
      fprintf(fp,"        zmax = max(max(ZI));\n");
      fprintf(fp,"        zmin = min(min(ZI)) / zmax;\n");
      fprintf(fp,"        ZI   = ZI / zmax;\n");
      fprintf(fp,"        zmax = 1;\n");
      fprintf(fp,"        Matplot1((ZI-zmin)/(zmax-zmin)*64,[L(jj),L(ii),");
      fprintf(fp,"U(jj),U(ii)]);\n");
      fprintf(fp,"        xset(\"colormap\",jetcolormap(64));\n");
      fprintf(fp,"        colorbar(zmin,zmax);\n");
      fprintf(fp,
           "        contour2d(X2,Y2,ZB,5,rect=[L(jj),L(ii),U(jj),U(ii)]);\n");
      fprintf(fp,"        xset(\"fpf\",\" \");\n");
      fprintf(fp,"        a = gca();\n");
      fprintf(fp,"        a.x_label.text = iStr(jj);\n");
      fprintf(fp,"        a.x_label.font_size = 3;\n");
      fprintf(fp,"        a.x_label.font_style = 4;\n");
      fprintf(fp,"        a.y_label.text = iStr(ii);\n");
      fprintf(fp,"        a.y_label.font_size = 3;\n");
      fprintf(fp,"        a.y_label.font_style = 4;\n");
      fwritePlotAxesNoGrid(fp);
    }
    else
    {
      fprintf(fp,"        imagesc(ZZ')\n");
      fprintf(fp,"        xtick = L(ii):(U(ii)-L(ii))/4:U(ii);\n");
      fprintf(fp,"        set(gca,'XTick',0:n/4:n);\n");
      fprintf(fp,"        set(gca,'XTickLabel', xtick);\n");
      fprintf(fp,"        ytick = L(jj):(U(jj)-L(jj))/4:U(jj);\n");
      fprintf(fp,"        set(gca,'YTick',0:n/4:n);\n");
      fprintf(fp,"        set(gca,'YTickLabel', ytick);\n");
      fprintf(fp,"        set(gca,'YDir', 'normal');\n");
      fprintf(fp,
         "        xlabel(iStr(jj),'FontWeight','bold','FontSize',12)\n");
      fprintf(fp,
         "        ylabel(iStr(ii),'FontWeight','bold','FontSize',12)\n");
      fprintf(fp,
         "        if (ii == 1 & jj == 2)\n");
      fprintf(fp,
         "          title('Prior','FontWeight','bold','FontSize',12)\n");
      fprintf(fp,"        end;\n");
      fwritePlotAxesNoGrid(fp);
      fprintf(fp,"        colorbar\n");
    }
  }
  fprintf(fp,"      end;\n");
  fprintf(fp,"    end;\n");
  fprintf(fp,"  end;\n");
  fprintf(fp,"end;\n");
  double minData=PSUADE_UNDEFINED;
  if (XChains == NULL) minData = Ymin;
  else
  {
    for (iChain = 0; iChain < nChains; iChain++) 
    { 
      if (chainStatus[iChain] == 0)
      {
        for (jj = 0; jj < chainCnt; jj++) 
        {
          ddata = XChains[iChain][jj][nInputs];
          if (ddata < minData) minData = ddata;
        }
      }
    }
  }
  fprintf(fp,"      subplot(nPlots,nPlots,1)\n");
  if (psPlotTool_ == 0)
  {
    fprintf(fp,"set(gcf,'NextPlot','add');\n");
    fprintf(fp,"axes;\n");
    if (minData < PSUADE_UNDEFINED)
      fprintf(fp,
          "h=title('MCMC Priors/Posteriors, best -log(likelihood)=%e',",
          minData);
    else fprintf(fp,"h=title('MCMC Prior/Posterior Distributions',");
    fprintf(fp,"'fontSize',12,'fontWeight','bold');\n");
    fprintf(fp,"set(gca,'Visible','off');\n");
    fprintf(fp,"set(h,'Visible','on');\n");
  }
  fprintf(fp,"negll = %e;\n", minData);
  if (chainCnt > 0)
  {
    for (kk = 0; kk < nInputs; kk++)
    {
      dmean = dstd = 0.0; 
      kk2 = 0;
      for (iChain = 0; iChain < nChains; iChain++) 
      { 
        if (chainStatus[iChain] == 0)
        {
          for (jj = 0; jj < chainCnt; jj++) 
          {
            dmean += XChains[iChain][jj][kk];
            kk2++;
          }
        }
      }
      if (kk2 > 0) dmean /= (double) kk2;
      for (iChain = 0; iChain < nChains; iChain++) 
      { 
        if (chainStatus[iChain] == 0)
        {
          for (jj = 0; jj < chainCnt; jj++) 
            dstd += pow(XChains[iChain][jj][kk]-dmean,2.0);
        }
      }
      if (kk2 > 0) dstd = sqrt(dstd/(double) kk2);
      fprintf(fp,
       "disp(['Stat for Input %d: mean,std = ',num2str(%e),' ',num2str(%e)])\n",
       kk+1,dmean,dstd);
    }
  }
  fprintf(fp,"disp('Lower diagonal plots: priors')\n");
  fprintf(fp,"disp('Upper diagonal plots: posteriors')\n");
  fclose(fp);
  if (psPlotTool_ == 1)
       printOutTS(PL_INFO, "MCMC: scilabmcmc2.sci file has been created.\n");
  else printOutTS(PL_INFO, "MCMC: matlabmcmc2.m file has been created.\n");
  return 0;
}

// ************************************************************************
// write to another matlab file the mse of the posterior sample with the
// experimental data
// ------------------------------------------------------------------------
int MCMCAnalyzer::genPostLikelihood(int nInputs, double *lower, 
                     double *upper, double *XRange, int numChains, 
                     int chainCnt, double ***XChains, int *chainStatus, 
                     int chainLimit, int *rsIndices, double *rsValues, 
                     int *designParams, int dnInputs, int dnSamples, 
                     double *dSamInputs, FuncApprox **faPtrs, 
                     FuncApprox **faPtrs1, int nOutputs, double *discOutputs,
                     double *discFuncConstantMeans, double *dSamMeans,
                     double *dSamStdevs)
                        
{
  int    iChain, ii, jj, ii2, kk2, dcnt;  
  double ddata,*XGuessS,*XDesignS,*YGuessS,*YDesignS,Ytemp,Ytemp2,ddata2;
  FILE   *fp;

  fp = fopen("matlablpostlikelihood.m", "w");
  if (fp == NULL) return -1;

  XGuessS  = new double[dnSamples * nInputs];
  XDesignS = new double[dnSamples * nInputs];
  YGuessS  = new double[dnSamples * nOutputs];
  YDesignS = new double[dnSamples * nOutputs];
  checkAllocate(YDesignS, "YDesignS in MCMC::genPostLikelihood");
  fprintf(fp, "A = [\n");
  for (iChain = 0; iChain < numChains; iChain++)
  {
    if (chainStatus[iChain] == 0)
    {
      for (ii = chainCnt-chainLimit; ii < chainCnt; ii++)
      {
        for (jj = 0; jj < nInputs; jj++)
        {
          if ((rsIndices == NULL || rsIndices[jj] >= 0) &&
              (designParams == NULL || designParams[jj] == 0))
          {
            ddata = XChains[iChain][ii][jj];
            fprintf(fp, "%e ", ddata);
          }
          else if (rsIndices != NULL && rsIndices[jj] < 0)
            fprintf(fp, "%e ", rsValues[jj]);
          else if (designParams != NULL && designParams[jj] != 0)
            fprintf(fp, "%e ", 0.5 * (upper[jj] + lower[jj]));
        }
        for (kk2 = 0; kk2 < dnSamples; kk2++)
        {
          dcnt = 0;
          for (ii2 = 0; ii2 < nInputs; ii2++)
          {
            XGuessS[kk2*nInputs+ii2] = XChains[iChain][ii][ii2]; 
            if (designParams != NULL && designParams[ii2] == 1)
            {
              XGuessS[kk2*nInputs+ii2] = dSamInputs[kk2*dnInputs+dcnt];
                 XDesignS[kk2*dnInputs+dcnt]=dSamInputs[kk2*dnInputs+dcnt];
              dcnt++;
            }
          }
        }
        for (ii2 = 0; ii2 < nOutputs; ii2++)
        {
          faPtrs[ii2]->evaluatePoint(dnSamples,XGuessS,
                                     &YGuessS[ii2*dnSamples]);
          if (faPtrs1 != NULL && faPtrs1[ii2] != NULL)
          {
            for (kk2 = 0; kk2 < dnSamples; kk2++)
              YDesignS[ii2*dnSamples+kk2] = 
                                     discOutputs[ii2*dnSamples+kk2];
          }
          else if (discFuncConstantMeans != NULL &&
                   discFuncConstantMeans[0] != PSUADE_UNDEFINED)
          {
            for (kk2 = 0; kk2 < dnSamples; kk2++)
              YDesignS[ii2*dnSamples+kk2] = discFuncConstantMeans[ii2];
          }
          else
          {
            for (kk2 = 0; kk2 < dnSamples; kk2++)
              YDesignS[ii2*dnSamples+kk2] = 0.0;
          }
        }
        ddata = ddata2 = 0.0;
        for (ii2 = 0; ii2 < nOutputs; ii2++)
        {
          for (kk2 = 0; kk2 < dnSamples; kk2++)
          {
            Ytemp = YGuessS[ii2*dnSamples+kk2] + 
                    YDesignS[ii2*dnSamples+kk2];
            Ytemp2 = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2]),2.0) /
                    (pow(dSamStdevs[kk2*nOutputs+ii2],2.0));
            ddata += Ytemp2;
            Ytemp2 = pow((Ytemp-dSamMeans[kk2*nOutputs+ii2]),2.0); 
            ddata2 += Ytemp2;
          }
        }
        ddata /= (dnSamples*nOutputs);
        ddata2 /= (dnSamples*nOutputs);
        fprintf(fp, "%e %e\n", ddata, ddata2);
      }
    }
  }
  fprintf(fp,"];\n");
  fprintf(fp,"figure(1)\n");
  fprintf(fp,"Y = A(:,%d);\n", nInputs+1);
  fprintf(fp,"subplot(1,2,1)\n");
  fprintf(fp,"hist(Y, 20);\n");
  fprintf(fp,"set(gca,'linewidth',2)\n");
  fprintf(fp,"set(gca,'fontweight','bold')\n");
  fprintf(fp,"set(gca,'fontsize',12)\n");
  fprintf(fp,"xlabel('Weighted MSE ','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"ylabel('Frequencies','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"grid on\n");
  fprintf(fp,"box on\n");
  fprintf(fp,"subplot(1,2,2)\n");
  fprintf(fp,"plot(Y, 'lineWidth', 2);\n");
  fprintf(fp,"set(gca,'linewidth',2)\n");
  fprintf(fp,"set(gca,'fontweight','bold')\n");
  fprintf(fp,"set(gca,'fontsize',12)\n");
  fprintf(fp,"xlabel('Sample number','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"ylabel('Weighted MSE','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"grid on\n");
  fprintf(fp,"box on\n");
  fprintf(fp,"figure(2)\n");
  fprintf(fp,"Y2 = A(:,%d);\n", nInputs+2);
  fprintf(fp,"subplot(1,2,1)\n");
  fprintf(fp,"hist(Y2, 20);\n");
  fprintf(fp,"set(gca,'linewidth',2)\n");
  fprintf(fp,"set(gca,'fontweight','bold')\n");
  fprintf(fp,"set(gca,'fontsize',12)\n");
  fprintf(fp,"xlabel('MSE','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"ylabel('Frequencies','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"grid on\n");
  fprintf(fp,"box on\n");
  fprintf(fp,"subplot(1,2,2)\n");
  fprintf(fp,"plot(Y2, 'lineWidth', 2);\n");
  fprintf(fp,"set(gca,'linewidth',2)\n");
  fprintf(fp,"set(gca,'fontweight','bold')\n");
  fprintf(fp,"set(gca,'fontsize',12)\n");
  fprintf(fp,"xlabel('Sample number','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"ylabel('Weighted MSE','FontWeight','bold','FontSize',12)\n");
  fprintf(fp,"grid on\n");
  fprintf(fp,"box on\n");
  fclose(fp);
  delete [] XGuessS;
  delete [] XDesignS;
  delete [] YGuessS;
  delete [] YDesignS;
  printOutTS(PL_INFO,"MCMC: matlabpostlikelihood.m file has been created.\n");
  return 0;
}

// ************************************************************************
// read spec (experimental data) file 
// ------------------------------------------------------------------------
double MCMCAnalyzer::readSpecFile(int nInputs, int nOutputs, int *dnSamp, 
                          int *dnInps, int **dParams, double **dSamIns, 
                          double **dMeans, double **dStds, int &combineFlag,
                          int printLevel)
{
  int    ii, jj, kk, cnt, dnSamples, dnInputs, *designParams=NULL;
  double *dSamInputs=NULL, *dSamMeans=NULL, *dSamStdevs=NULL;
  char   lineIn[1001], cfname[1001], cword[101], cword2[5001];
  FILE   *fp=NULL;

  if (printLevel > 0)
  {
    printOutTS(PL_INFO,"*** NEED DATA TO CREATE LIKELIHOOD FUNCTION: \n\n");
    printOutTS(PL_INFO,"MCMC will create a Gaussian likelihood function.\n");
    printOutTS(PL_INFO,"Please provide a data file containing design\n");
    printOutTS(PL_INFO,"parameter values, mean, and std. dev. of the\n");
    printOutTS(PL_INFO,"observation data for each output.\n");
    printOutTS(PL_INFO,"NOTE: Design parameters should be defined in the\n");
    printOutTS(PL_INFO,"   observation data file if the data used in MCMC\n");
    printOutTS(PL_INFO,"   are collected at different design points.\n");
    printOutTS(PL_INFO,"IMPORTANT: IF m DESIGN PARAMETERS ARE SPECIFIED,\n");
    printOutTS(PL_INFO,"   YOU NEED TO SPECIFY WHICH ONES THEY ARE.\n");
    printOutTS(PL_INFO,"   THESE DESIGN PARAMETERS WILL BE EXCLUDED FROM\n");
    printOutTS(PL_INFO,"   BEING IN THE CALIBRATION PARAMETER SET.\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO,"** OBSERVATION DATA FILE FORMAT : (O1 = Output 1, \n");
    printOutTS(PL_INFO,"        M   - no. of design parameters, \n");
    printOutTS(PL_INFO,"        K   - no. of model outputs, \n");
    printOutTS(PL_INFO,"        P   - no. of experiments \n");
    printOutTS(PL_INFO,"        O1m - Output 1 mean\n");
    printOutTS(PL_INFO,"        O1s - Output 1 std. dev.\n");
    printOutTS(PL_INFO,"        OKs - Output K std. dev.\n");
    printOutTS(PL_INFO,"PSUADE_BEGIN\n");
    printOutTS(PL_INFO,"<P> <K> <M> <design parameter identifiers>\n");
    printOutTS(PL_INFO,"1 <design values...> <O1m> <O1s> ... <OKs> \n");
    printOutTS(PL_INFO,"2 <design values...> <O1m> <O1s> ... <OKs> \n");
    printOutTS(PL_INFO,"...\n");
    printOutTS(PL_INFO,"P <design values...> <O1m> <O1s> ... <OK> \n");
    printOutTS(PL_INFO,"PSUADE_END\n");
    printOutTS(PL_INFO,"AVERAGE (Optional command)\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO,"The likelihood function is in the form of:\n");
    printOutTS(PL_INFO,"  C exp(-0.5*S) \n");
    printOutTS(PL_INFO,"where C is the normalization constant and\n");
    printOutTS(PL_INFO,
         "  S=sum_{p=1}^P sum_{k=1)^K (Y_pk-m_pk)^2/sd_pk^2\n");
    printOutTS(PL_INFO,"where K is the number of outputs and m_pk and\n");
    printOutTS(PL_INFO,"  sd_pk are the mean and std. dev. of output k\n");
    printOutTS(PL_INFO,"  of experiment k.\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO,
         "NOTE: Alternately, your simulator (or response surface)\n");
    printOutTS(PL_INFO,
         "   output may be some error measure from comparison of\n");
    printOutTS(PL_INFO,
         "   all model outputs with observation data. In this\n");
    printOutTS(PL_INFO,
         "   case, set nOutputs=1, mean=0 and std. dev.=1 in the\n");
    printOutTS(PL_INFO,
         "   specification file (that is, your simulation output\n");
    printOutTS(PL_INFO,
         "   is S above, and MCMC will compute likelihood as :\n");
    printOutTS(PL_INFO,"   C exp(-0.5 S).\n");
    printOutTS(PL_INFO,
         "   However, if you choose this option and you desire to include\n");
    printOutTS(PL_INFO,
         "   response surface uncertainty, you need to be very careful.\n");
  }
 
  printf("==> Enter name of the spec file for building likelihood function: ");
  scanf("%s", cfname);
  fgets(lineIn, 1000, stdin);
  kk = strlen(cfname);
  if (kk <= 1000)
  {
    cfname[kk] = '\0';
    fp = fopen(cfname, "r");
    if (fp == NULL)
    {
      printOutTS(PL_ERROR,"MCMC ERROR : cannot open spec file %s.\n",cfname);
      return PSUADE_UNDEFINED;
    }
  }
  else
  {
    printOutTS(PL_ERROR,"MCMC ERROR: file name too long.\n");
    return PSUADE_UNDEFINED;
  }
  int lineCnt = 0;
  lineIn[0] = '#';
  combineFlag = 0;
  while (lineIn[0] == '#')
  {
    fgets(lineIn, 2000, fp);
    cword2[0] = 'N';
    sscanf(lineIn, "%s %s", cword, cword2);
    if (!strcmp(cword2, "AVERAGE"))
    {
      combineFlag = 1;
      printf("MCMC INFO: likelihood = 1/m sum_i (exp_i-model_i)^2\n");
    }
    lineCnt++;
  }
  sscanf(lineIn, "%s", cword); 
  if (!strcmp(cword, "PSUADE_BEGIN")) lineCnt++;
  lineCnt--;
  fclose(fp);
  fp = fopen(cfname, "r");
  for (ii = 0; ii < lineCnt; ii++) fgets(lineIn, 2000, fp);
  fscanf(fp, "%d %d %d", &dnSamples, &kk, &dnInputs);
  if (dnSamples <= 0)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: no. of experiments <= 0.\n");
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  printOutTS(PL_INFO,"SPEC FILE: Number of experiments = %d\n",dnSamples);

  if (kk != nOutputs)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: nOutputs in spec file experiment\n");
    printOutTS(PL_ERROR,"     does not match the PSUADE file nOutputs.\n");
    printOutTS(PL_ERROR,"     %d versus %d\n", kk, nOutputs);
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  printOutTS(PL_INFO,"SPEC FILE: Number of outputs = %d\n",nOutputs);

  if (dnInputs < 0)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: number of design variables < 0.\n");
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  if (dnInputs > nInputs)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: number of design variables %d\n",
               dnInputs);
    printOutTS(PL_ERROR,"     cannot be larger than the total number\n");
    printOutTS(PL_ERROR,"     of inputs %d.\n", nInputs);
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  printOutTS(PL_INFO,"SPEC FILE: Number of design parameters = %d\n",dnInputs);
  if (dnInputs > 0)
  {
    designParams = new int[nInputs];
    for (ii = 0; ii < nInputs; ii++) designParams[ii] = 0;
    cnt = 0;
    for (ii = 0; ii < dnInputs; ii++)
    {
      fscanf(fp, "%d", &kk);
      if (kk <= 0 || kk > nInputs)
      {
        printOutTS(PL_ERROR,"MCMC ERROR: invalid design parameter %d.\n",
                   kk);
        printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
        printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
        fclose(fp);
        delete [] designParams;
        return PSUADE_UNDEFINED;
      }
      if (kk <= cnt)
      {
        printOutTS(PL_ERROR,"MCMC ERROR: design parameters should be in\n");
        printOutTS(PL_ERROR,"            ascending order.\n");
        printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
        printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
        fclose(fp);
        delete [] designParams;
        return PSUADE_UNDEFINED;
      }
      designParams[kk-1] = 1;
      printOutTS(PL_INFO,"SPEC FILE: input %d is a design parameter\n", kk);
      cnt = kk;
    }
    dSamInputs = new double[dnSamples*dnInputs];
  }
  dSamMeans = new double[dnSamples*nOutputs];
  dSamStdevs = new double[dnSamples*nOutputs];
  for (ii = 0; ii < dnSamples; ii++)
  {
    fscanf(fp, "%d", &kk);
    if (kk != ii+1)
    {
      printOutTS(PL_ERROR,"MCMC ERROR: invalid experiment number %d\n",kk);
      printOutTS(PL_ERROR,"     at line %d in the spec file.\n",ii+2);
      printOutTS(PL_ERROR,"            (Expecting %d).\n", ii+1);
      printOutTS(PL_ERROR,"==> check line %d\n", ii+3);
      printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
      printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
      fclose(fp);
      if (dSamInputs != NULL) delete [] dSamInputs;
      delete [] dSamMeans;
      delete [] dSamStdevs;
      return PSUADE_UNDEFINED;
    }
    if (printLevel > 0)
      printOutTS(PL_INFO,"Calibration Data Set %d\n", kk);
    for (jj = 0; jj < dnInputs; jj++)
    {
      fscanf(fp, "%lg", &dSamInputs[ii*dnInputs+jj]);
      if (printLevel > 0)
        printOutTS(PL_INFO,"   Design parameter %d = %e\n", jj+1,
                   dSamInputs[ii*dnInputs+jj]);
    }
    for (jj = 0; jj < nOutputs; jj++)
    {
      fscanf(fp, "%lg %lg", &dSamMeans[ii*nOutputs+jj],
             &dSamStdevs[ii*nOutputs+jj]);
      if (printLevel > 0)
         printOutTS(PL_INFO,"      Data mean/stdev = %16.8e %16.8e\n",
                    dSamMeans[ii*nOutputs+jj],dSamStdevs[ii*nOutputs+jj]);
      if (dSamStdevs[ii*nOutputs+jj] < 0.0)
      {
        fclose(fp);
        printOutTS(PL_ERROR,"MCMC ERROR: std dev in spec file <= 0.\n");
        printOutTS(PL_ERROR,"==> check the last entry in line %d\n", ii+3);
        printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
        printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
        if (dSamInputs != NULL) delete [] dSamInputs;
        delete [] dSamMeans;
        delete [] dSamStdevs;
        return PSUADE_UNDEFINED;
      }
    }
  }
  fclose(fp);
  (*dnSamp)  = dnSamples;
  (*dnInps)  = dnInputs;
  (*dParams) = designParams;
  (*dSamIns) = dSamInputs;
  (*dMeans)  = dSamMeans;
  (*dStds)   = dSamStdevs;
  printEquals(PL_INFO, 0);
  return 0;
}

// ************************************************************************
// read spec (experimental data) file - more structured version
// ------------------------------------------------------------------------
double MCMCAnalyzer::readSpecFile2(int nInputs, int nOutputs,  
                          psIVector &dParams, psMatrix &dSamInputs, 
                          psMatrix &dSamMeans, psMatrix &dSamStdvs, 
                          int &combineFlag, int printLevel)
{
  int    ii, jj, kk, cnt, dnSamples, dnInputs;
  double ddata, ddata2;
  char   lineIn[1001], cfname[1001], cword[101], cword2[5001];
  FILE   *fp=NULL;

  if (printLevel > 0)
  {
    printOutTS(PL_INFO,"*** NEED DATA TO CREATE LIKELIHOOD FUNCTION: \n\n");
    printOutTS(PL_INFO,"MCMC will create a Gaussian likelihood function.\n");
    printOutTS(PL_INFO,"Please provide a data file containing design\n");
    printOutTS(PL_INFO,"parameter values, mean, and std. dev. of the\n");
    printOutTS(PL_INFO,"observation data for each output.\n");
    printOutTS(PL_INFO,"NOTE: Design parameters should be defined in the\n");
    printOutTS(PL_INFO,"   observation data file if the data used in MCMC\n");
    printOutTS(PL_INFO,"   are collected at different design points.\n");
    printOutTS(PL_INFO,"IMPORTANT: IF m DESIGN PARAMETERS ARE SPECIFIED,\n");
    printOutTS(PL_INFO,"   YOU NEED TO SPECIFY WHICH ONES THEY ARE.\n");
    printOutTS(PL_INFO,"   THESE DESIGN PARAMETERS WILL BE EXCLUDED FROM\n");
    printOutTS(PL_INFO,"   BEING IN THE CALIBRATION PARAMETER SET.\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO,"** OBSERVATION DATA FILE FORMAT : (O1 = Output 1, \n");
    printOutTS(PL_INFO,"        M   - no. of design parameters, \n");
    printOutTS(PL_INFO,"        K   - no. of model outputs, \n");
    printOutTS(PL_INFO,"        P   - no. of experiments \n");
    printOutTS(PL_INFO,"        O1m - Output 1 mean\n");
    printOutTS(PL_INFO,"        O1s - Output 1 std. dev.\n");
    printOutTS(PL_INFO,"        OKs - Output K std. dev.\n");
    printOutTS(PL_INFO,"PSUADE_BEGIN\n");
    printOutTS(PL_INFO,"<P> <K> <M> <design parameter identifiers>\n");
    printOutTS(PL_INFO,"1 <design values...> <O1m> <O1s> ... <OKs> \n");
    printOutTS(PL_INFO,"2 <design values...> <O1m> <O1s> ... <OKs> \n");
    printOutTS(PL_INFO,"...\n");
    printOutTS(PL_INFO,"P <design values...> <O1m> <O1s> ... <OK> \n");
    printOutTS(PL_INFO,"PSUADE_END\n");
    printOutTS(PL_INFO,"AVERAGE (Optional command)\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO,"The likelihood function is in the form of:\n");
    printOutTS(PL_INFO,"  C exp(-0.5*S) \n");
    printOutTS(PL_INFO,"where C is the normalization constant and\n");
    printOutTS(PL_INFO,
         "  S=sum_{p=1}^P sum_{k=1)^K (Y_pk-m_pk)^2/sd_pk^2\n");
    printOutTS(PL_INFO,"where K is the number of outputs and m_pk and\n");
    printOutTS(PL_INFO,"  sd_pk are the mean and std. dev. of output k\n");
    printOutTS(PL_INFO,"  of experiment k.\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO,
         "NOTE: Alternately, your simulator (or response surface)\n");
    printOutTS(PL_INFO,
         "   output may be some error measure from comparison of\n");
    printOutTS(PL_INFO,
         "   all model outputs with observation data. In this\n");
    printOutTS(PL_INFO,
         "   case, set nOutputs=1, mean=0 and std. dev.=1 in the\n");
    printOutTS(PL_INFO,
         "   specification file (that is, your simulation output\n");
    printOutTS(PL_INFO,
         "   is S above, and MCMC will compute likelihood as :\n");
    printOutTS(PL_INFO,"   C exp(-0.5 S).\n");
    printOutTS(PL_INFO,
         "   However, if you choose this option and you desire to include\n");
    printOutTS(PL_INFO,
         "   response surface uncertainty, you need to be very careful.\n");
  }
 
  printf("==> Enter name of the spec file for building likelihood function: ");
  scanf("%s", cfname);
  fgets(lineIn, 1000, stdin);
  kk = strlen(cfname);
  if (kk <= 1000)
  {
    cfname[kk] = '\0';
    fp = fopen(cfname, "r");
    if (fp == NULL)
    {
      printOutTS(PL_ERROR,"MCMC ERROR : cannot open spec file %s.\n",cfname);
      return PSUADE_UNDEFINED;
    }
  }
  else
  {
    printOutTS(PL_ERROR,"MCMC ERROR: file name too long.\n");
    return PSUADE_UNDEFINED;
  }
  int lineCnt = 0;
  lineIn[0] = '#';
  combineFlag = 0;
  while (lineIn[0] == '#')
  {
    fgets(lineIn, 2000, fp);
    cword2[0] = 'N';
    sscanf(lineIn, "%s %s", cword, cword2);
    if (!strcmp(cword2, "AVERAGE"))
    {
      combineFlag = 1;
      printf("MCMC INFO: likelihood = 1/m sum_i (exp_i-model_i)^2\n");
    }
    lineCnt++;
  }
  sscanf(lineIn, "%s", cword); 
  if (!strcmp(cword, "PSUADE_BEGIN")) lineCnt++;
  lineCnt--;
  fclose(fp);
  fp = fopen(cfname, "r");
  for (ii = 0; ii < lineCnt; ii++) fgets(lineIn, 2000, fp);
  fscanf(fp, "%d %d %d", &dnSamples, &kk, &dnInputs);
  if (dnSamples <= 0)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: no. of experiments <= 0.\n");
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  printOutTS(PL_INFO,"SPEC FILE: Number of experiments = %d\n",dnSamples);

  if (kk != nOutputs)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: nOutputs in spec file experiment\n");
    printOutTS(PL_ERROR,"     does not match the PSUADE file nOutputs.\n");
    printOutTS(PL_ERROR,"     %d versus %d\n", kk, nOutputs);
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  printOutTS(PL_INFO,"SPEC FILE: Number of outputs = %d\n",nOutputs);

  if (dnInputs < 0)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: number of design variables < 0.\n");
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  if (dnInputs > nInputs)
  {
    printOutTS(PL_ERROR,"MCMC ERROR: number of design variables %d\n",
               dnInputs);
    printOutTS(PL_ERROR,"     cannot be larger than the total number\n");
    printOutTS(PL_ERROR,"     of inputs %d.\n", nInputs);
    fclose(fp);
    return PSUADE_UNDEFINED;
  }
  printOutTS(PL_INFO,"SPEC FILE: Number of design parameters = %d\n",dnInputs);
  if (dnInputs > 0)
  {
    dParams.setLength(nInputs);
    cnt = 0;
    for (ii = 0; ii < dnInputs; ii++)
    {
      fscanf(fp, "%d", &kk);
      if (kk <= 0 || kk > nInputs)
      {
        printOutTS(PL_ERROR,"MCMC ERROR: invalid design parameter %d.\n",
                   kk);
        printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
        printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
        fclose(fp);
        return PSUADE_UNDEFINED;
      }
      if (kk <= cnt)
      {
        printOutTS(PL_ERROR,"MCMC ERROR: design parameters should be in\n");
        printOutTS(PL_ERROR,"            ascending order.\n");
        printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
        printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
        fclose(fp);
        return PSUADE_UNDEFINED;
      }
      dParams[kk-1] = 1;
      printOutTS(PL_INFO,"SPEC FILE: input %d is a design parameter\n", kk);
      cnt = kk;
    }
    dSamInputs.setDim(dnSamples, dnInputs);
  }
  dSamMeans.setDim(dnSamples, nOutputs);
  dSamStdvs.setDim(dnSamples, nOutputs);
  for (ii = 0; ii < dnSamples; ii++)
  {
    fscanf(fp, "%d", &kk);
    if (kk != ii+1)
    {
      printOutTS(PL_ERROR,"MCMC ERROR: invalid experiment number %d\n",kk);
      printOutTS(PL_ERROR,"     at line %d in the spec file.\n",ii+2);
      printOutTS(PL_ERROR,"            (Expecting %d).\n", ii+1);
      printOutTS(PL_ERROR,"==> check line %d\n", ii+3);
      printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
      printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
      fclose(fp);
      return PSUADE_UNDEFINED;
    }
    if (printLevel > 0)
      printOutTS(PL_INFO,"Calibration Data Set %d\n", kk);
    for (jj = 0; jj < dnInputs; jj++)
    {
      fscanf(fp, "%lg", &ddata);
      dSamInputs.setEntry(ii, jj, ddata);
      if (printLevel > 0)
        printOutTS(PL_INFO,"   Design parameter %d = %e\n",jj+1,ddata);
    }
    for (jj = 0; jj < nOutputs; jj++)
    {
      fscanf(fp, "%lg %lg", &ddata, &ddata2);
      dSamMeans.setEntry(ii, jj, ddata);
      dSamStdvs.setEntry(ii, jj, ddata2);
      if (printLevel > 0)
         printOutTS(PL_INFO,"      Data mean/stdev = %16.8e %16.8e\n",
                    ddata, ddata2);
      if (ddata2 < 0.0)
      {
        fclose(fp);
        printOutTS(PL_ERROR,"MCMC ERROR: std dev in spec file <= 0.\n");
        printOutTS(PL_ERROR,"==> check the last entry in line %d\n", ii+3);
        printOutTS(PL_ERROR,"   ** Use showformat in command line mode\n");
        printOutTS(PL_ERROR,"   ** to verify your spec file format.\n");
        return PSUADE_UNDEFINED;
      }
    }
  }
  fclose(fp);
  printEquals(PL_INFO, 0);
  return 0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
int MCMCAnalyzer::setParams(int argc, char **argv)
{
  char  *request = (char *) argv[0];
  Analyzer::setParams(argc, argv);
  if      (!strcmp(request, "setsim")) mode_ = 1;
  else if (!strcmp(request, "MCMC_brute_force")) bfmode_ = 1;
  else if (!strcmp(request, "MCMC_gibbs")) bfmode_ = 0;
  else
  {
    printOutTS(PL_ERROR,"MCMCAnalyzer ERROR: setParams - not valid.\n");
    exit(1);
  }
  return 0;
}

// ************************************************************************
// check convergence 
// ------------------------------------------------------------------------
int MCMCAnalyzer::checkConvergence(int num, double *means, double *stds,
                                   int leng)
{
  int    ii;
  double WStat, BStat, ddata, ddata2, thresh=1.02;
  WStat = 0.0;
  for (ii = 0; ii < num; ii++) WStat += stds[ii];
  WStat /= (double) num;
  if (WStat < 0) WStat = PSUADE_UNDEFINED;
  ddata = 0.0;
  for (ii = 0; ii < num; ii++) ddata += means[ii];
  ddata /= (double) num;
  BStat = 0.0;
  for (ii = 0; ii < num; ii++) 
     BStat += pow(means[ii]-ddata,2.0);
  BStat = BStat / (num - 1.0) * leng;
  ddata = (1 - 1.0/leng) * WStat + BStat / leng;
  ddata = ddata / WStat * (num + 1) / num - 
          (leng - 1.0) / (double) (leng * num); 
  if (ddata < 0) ddata2 = PSUADE_UNDEFINED;
  else           ddata2 = sqrt(ddata);
  printf("MCMC convergence check: %e <? %e\n", ddata2, thresh);
  if (ddata2 < thresh) return 1;
  else                 return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MCMCAnalyzer& MCMCAnalyzer::operator=(const MCMCAnalyzer &)
{
   printOutTS(PL_ERROR,
        "MCMCAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// create response surface
// ------------------------------------------------------------------------
int MCMCAnalyzer::readIndexFile(PsuadeData *dataPtr, int nInputs,
                      int *designParams, int *rsIndices, double *rsValues)
{
  int    kk, ii;
  char   *rsFile;
  FILE   *fp;
  pData  pPtr;

  dataPtr->getParameter("ana_rsindexfile", pPtr);
  rsFile = pPtr.strArray_[0];
  if (strcmp(rsFile, "NONE"))
  {
    printOutTS(PL_INFO,"A response surface index file has been specified.\n");
    fp = fopen(rsFile, "r");
    if (fp == NULL)
    {
      printOutTS(PL_ERROR,
           "MCMC ERROR: rs_index_file %s not found.\n",rsFile);
      return -1;
    }
    else
    {
      printOutTS(PL_INFO,"INFO: rs_index_file %s found.\n",rsFile);
      fscanf(fp,"%d", &kk);
      if (kk != nInputs)
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: invalid nInputs in rs_index_file (%d != %d).\n",
             kk, nInputs);
        printOutTS(PL_ERROR,"  Data format should be: \n");
        printOutTS(PL_ERROR,
             "  line 1: nInputs in rs data (driver) file\n");
        printOutTS(PL_ERROR,
             "  line 2: 1 <1 or 0> <default value if first number==0>\n");
        printOutTS(PL_ERROR,
             "  line 3: 2 <2 or 0> <0 if first number != 0>\n");
        printOutTS(PL_ERROR,
             "  line 4: 3 <3 or 0> <default value if first number==0>\n");
        printOutTS(PL_ERROR,
             "  line 5: 4 <4 or 0> <0 if first number != 0>\n");
        printOutTS(PL_ERROR,"  ...\n");
        fclose(fp);
        return -1;
      }
      for (ii = 0; ii < nInputs; ii++)
      {
        rsIndices[ii] = 0;
        fscanf(fp, "%d", &kk);
        if (kk != ii+1)
        {
          printOutTS(PL_ERROR,
               "MCMC ERROR: 1st index in indexFile = %d (must be %d]).\n",
               kk, ii+1);
          printOutTS(PL_ERROR,"  Data format should be: \n");
          printOutTS(PL_ERROR,
               "  line 1: nInputs in rs data (driver) file\n");
          printOutTS(PL_ERROR,
               "  line 2: 1 <1 or 0> <default value if first number==0>\n");
          printOutTS(PL_ERROR,
               "  line 3: 2 <2 or 0> <0 if first number != 0>\n");
          printOutTS(PL_ERROR,
               "  line 4: 3 <3 or 0> <default value if first number==0>\n");
          printOutTS(PL_ERROR,
               "  line 5: 4 <4 or 0> <0 if first number != 0>\n");
          printOutTS(PL_ERROR,"  ...\n");
          fclose(fp);
          return -1;
        }
        fscanf(fp, "%d", &rsIndices[ii]);
        if (rsIndices[ii] == 0)
          printOutTS(PL_INFO,"MCMC INFO: input %3d inactive\n",ii+1);

        if (rsIndices[ii] == 0 && designParams != NULL && 
            designParams[ii] == 1)
        {
          printOutTS(PL_ERROR,
               "MCMC ERROR: inactive input %d cannot be design parameter\n",
               ii+1);
          fclose(fp);
          return -1;
        }

        if (rsIndices[ii] < 0 || rsIndices[ii] > nInputs)
        {
          printOutTS(PL_ERROR,
               "MCMC readIndexFile ERROR: input %3d = %d invalid\n",ii+1,
               rsIndices[ii]);
          fclose(fp);
          return -1;
        }
        rsIndices[ii]--;
        fscanf(fp, "%lg", &rsValues[ii]);
      }
      fclose(fp);
      printOutTS(PL_INFO, "Response surface index information: \n");
      for (ii = 0; ii < nInputs; ii++)
      {
        if (rsIndices[ii] == -1)
          printOutTS(PL_INFO, "Input %4d: fixed at default value  = %e\n",
                     ii+1, rsValues[ii]);
        else if (rsIndices[ii] >= 1000)
          printOutTS(PL_INFO, "Input %4d: uncertain, sample index = %4d\n",
                     ii+1, rsIndices[ii]-999);
        else if (designParams != NULL && designParams[ii] == 1)
          printOutTS(PL_INFO, "Input %4d: design parameter\n", ii+1);
        else
          printOutTS(PL_INFO, "Input %4d: calibration parameter\n",ii+1);
      }
    }
  }
  return 0;
}

// ************************************************************************
// read in response surface index file
// rsIndices[ii] == -1  ==> fixed parameters
// rsIndices[ii] == -2  ==> uncertain parameters 
// samMatrix ==> if there is any uncertain parameter, it has the sample
// ------------------------------------------------------------------------
int MCMCAnalyzer::readIndexFile2(PsuadeData *dataPtr, psIVector &rsIndices, 
                          psVector &rsValues, psMatrix &samMatrix)
{
  int    kk, ii, nInputs, nSamp, numUParams=0, status;
  double ddata;
  char   *rsFile;
  pData  pPtr;
  string   sampleFileName, str2;
  ifstream infile;

  dataPtr->getParameter("input_ninputs", pPtr);
  nInputs = pPtr.intData_;
  dataPtr->getParameter("ana_rsindexfile", pPtr);
  rsFile = pPtr.strArray_[0];
  if (strcmp(rsFile, "NONE"))
  {
    printOutTS(PL_INFO,
         "MCMC: A response surface index file has been specified.\n");

    infile.open(rsFile, ios::in);
    if (!infile.is_open())
    {
      printOutTS(PL_ERROR,
         "MCMC ERROR: rs_index_file %s not found.\n",rsFile);
      return -1;
    }
    else
    {
      printOutTS(PL_INFO,"INFO: rs_index_file %s found.\n",rsFile);
      infile >> kk;
      if (kk != nInputs)
      {
        printOutTS(PL_ERROR,
             "MCMC ERROR: invalid nInputs in rs_index_file (%d != %d).\n",
             kk, nInputs);
        printOutTS(PL_ERROR,"  Data format should be: \n");
        printOutTS(PL_ERROR,
             "  line 1: nInputs in rs data (driver) file\n");
        printOutTS(PL_ERROR,
             "  line 2: 1 1 <Input 1 is a calibration parameter>\n");
        printOutTS(PL_ERROR,
             "  line 3: 2 0 val <Input 2 is a parameter fixed at val>\n");
        printOutTS(PL_ERROR,
             "  line 4: 3 999 1 uSamp <Input 3 is an uncertain parameter>\n");
        printOutTS(PL_ERROR,"  (999   means uncertain parameter)\n");
        printOutTS(PL_ERROR,"  (1     means first column in uSamp)\n");
        printOutTS(PL_ERROR,"  (uSamp is a sample file)\n");
        printOutTS(PL_ERROR,"  ...\n");
        infile.close();
        return -1;
      }
      rsIndices.setLength(nInputs);
      rsValues.setLength(nInputs);
      sampleFileName = "none";
      for (ii = 0; ii < nInputs; ii++)
      {
        rsIndices[ii] = 0;
        infile >> kk;
        if (kk != ii+1)
        {
          printOutTS(PL_ERROR,
               "MCMC ERROR: 1st index in indexFile = %d (must be %d)).\n",
               kk, ii+1);
          printOutTS(PL_ERROR,"  Data format should be: \n");
          printOutTS(PL_ERROR,
               "  line 1: nInputs in rs data (driver) file\n");
          printOutTS(PL_ERROR,
               "  line 2: 1 1 <Input 1 is a calibration parameter>\n");
          printOutTS(PL_ERROR,
               "  line 3: 2 0 val <Input 2 is a parameter fixed at val>\n");
          printOutTS(PL_ERROR,
               "  line 4: 3 999 1 uSamp <Input 3 is an uncertain parameter>\n");
          printOutTS(PL_ERROR,"  (999   means uncertain parameter)\n");
          printOutTS(PL_ERROR,"  (1     means first column in uSamp)\n");
          printOutTS(PL_ERROR,"  (uSamp is a sample file)\n");
          printOutTS(PL_ERROR,"  ...\n");
          infile.close();
          return -1;
        }
        infile >> kk;
        rsIndices[ii] = kk;
        if (rsIndices[ii] == 0)
          printOutTS(PL_INFO,"MCMC INFO: input %3d inactive\n",ii+1);

        if (rsIndices[ii] == 999)
        {
          infile >> kk;
          infile >> str2;
          if (sampleFileName.compare("none") == 0) sampleFileName = str2;
          else if (sampleFileName.compare(str2) != 0)
          {
            printOutTS(PL_ERROR,
                 "MCMC readIndexFile2 ERROR: sample file for all\n");
            printOutTS(PL_ERROR,
                 "         uncertain parameters must be the same.\n");
            exit(1);
          }
          numUParams++;
          rsIndices[ii] = kk + 1000;
        }
        else if (rsIndices[ii] < 0 || rsIndices[ii] > nInputs)
        {
          printOutTS(PL_ERROR,
               "MCMC readIndexFile2 ERROR: input %3d = %d invalid\n",ii+1,
               rsIndices[ii]);
          infile.close();
          return -1;
        }

        rsIndices[ii] = rsIndices[ii] - 1;
        if (rsIndices[ii] >= -1 && rsIndices[ii] < nInputs)
        {
          infile >> ddata;
          rsValues[ii] = ddata;
        }
      }
      infile.close();
      printOutTS(PL_INFO, "Response surface index information: \n");
      for (ii = 0; ii < nInputs; ii++)
      {
        if (rsIndices[ii] == -1)
          printOutTS(PL_INFO, "Input %4d: fixed at default value  = %e\n",
                     ii+1, rsValues[ii]);
        else if (rsIndices[ii] >= 1000)
          printOutTS(PL_INFO, "Input %4d: uncertain, sample index = %4d\n",
                     ii+1, rsIndices[ii]-999);
        else
          printOutTS(PL_INFO, "Input %4d: calibration/design parameter\n",
                     ii+1);
      }
      if (numUParams > 0)
      {
        status = readSampleInputFile(sampleFileName.c_str(),rsIndices,
                                     samMatrix);
        if (status < 0) return -1;
      }
    }
  }
  return 0;
}

// ************************************************************************
// clean up
// ------------------------------------------------------------------------
void MCMCAnalyzer::cleanUp()
{
  char charString[200];
  if (means_) delete[] means_;
  if (sigmas_) delete[] sigmas_;
  if (mostLikelyInput_) delete[] mostLikelyInput_;
  if (mostLikelyOutput_) delete[] mostLikelyOutput_;
  means_ = NULL;
  sigmas_ = NULL;
  mostLikelyInput_ = NULL;
  mostLikelyOutput_ = NULL;
  FILE *fp = fopen("psuade_stop", "r");
  if (fp != NULL)
  {
     fclose(fp);
     fp = NULL;
     printOutTS(PL_INFO,
          "MCMC INFO: psuade_stop FILE FOUND. WILL BE REMOVED\n");
     strcpy(charString, "psuade_stop");
     unlink(charString);
  }
  fp = fopen("psuade_master", "r");
  if (fp != NULL)
  {
     fclose(fp);
     fp = NULL;
     printOutTS(PL_INFO,
          "MCMC INFO: psuade_master FILE FOUND. WILL BE REMOVED\n");
     strcpy(charString, "psuade_master");
     unlink(charString);
  }
  fp = fopen("psuade_nomaster", "r");
  if (fp != NULL)
  {
     fclose(fp);
     fp = NULL;
     printOutTS(PL_INFO,
          "MCMC INFO: psuade_nomaster FILE FOUND. WILL BE REMOVED\n");
     strcpy(charString, "psuade_nomaster");
     unlink(charString);
  }
  fp = fopen("psuade_gm", "r");
  if (fp != NULL)
  {
     fclose(fp);
     fp = NULL;
     printOutTS(PL_INFO,
          "MCMC INFO: psuade_gm FILE FOUND. WILL BE REMOVED\n");
     strcpy(charString, "psuade_gm");
     unlink(charString);
  }
  fp = fopen("psuade_nogm", "r");
  if (fp != NULL)
  {
     fclose(fp);
     fp = NULL;
     printOutTS(PL_INFO,
          "MCMC INFO: psuade_nogm FILE FOUND. WILL BE REMOVED\n");
     strcpy(charString, "psuade_nogm");
     unlink(charString);
  }
}

// ************************************************************************
// display banner for brute force version 2
// ************************************************************************
void MCMCAnalyzer::displayBanner_bf2(int printLevel)
{
  if (printLevel > 0)
  {
    printOutTS(PL_INFO,"TO GAIN ACCESS TO DIFFERENT OPTIONS: TURN ON\n\n");
    printOutTS(PL_INFO,
         " * ana_expert to finetune MCMC parameters, \n");
    printOutTS(PL_INFO,
         "   (e.g. sample size for burn-in can be adjusted).\n");
    printOutTS(PL_INFO,
         " * rs_expert to customize response surface for MCMC,\n");
    printDashes(PL_INFO,0);

    printOutTS(PL_INFO,"FEATURES AVAILABLE IN THE CURRENT VERSION OF MCMC:\n");
    printOutTS(PL_INFO,
         " * Support likelihood functions from multiple outputs\n");
    printOutTS(PL_INFO," * Option to include response surface errors \n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode.\n");
    printOutTS(PL_INFO,
         " * Option to include model form errors in the form of\n");
    printOutTS(PL_INFO,"   discrepancy models.\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode.\n");
    printOutTS(PL_INFO,
         " * Option to set some inputs as design parameters\n");
    printOutTS(PL_INFO,
         "   - to be specified in the observation data spec file\n");
    printOutTS(PL_INFO,
         " * Option to disable some parameters (set to default)\n");
    printOutTS(PL_INFO,
         "   - in case these parameters are not to be calibrated\n");
    printOutTS(PL_INFO,
         "   - use rs_index_file in PSUADE's ANALYSIS section\n");
    printOutTS(PL_INFO,"   - not available with discrepancy modeling\n");
    printOutTS(PL_INFO,
         " * Option to set some parameters as uncertain but not calibrated\n");
    printOutTS(PL_INFO,
         "   - a sample file is to be provided to characterize uncertainty\n");
    printOutTS(PL_INFO,
         "   - use rs_index_file in PSUADE's ANALYSIS section\n");
    printOutTS(PL_INFO,
         " * Option to generate a posterior sample\n");
    printOutTS(PL_INFO,"   - can be selected in ana_expert mode\n");
    printOutTS(PL_INFO,
         " * MCMC can be terminated gracefully by creating a file\n");
    printOutTS(PL_INFO,
         "   named 'psuade_stop' in the same directory during the\n");
    printOutTS(PL_INFO,"   run (in case it takes too long).\n");
    printEquals(PL_INFO, 0);
  }
}

// ************************************************************************
// display banner for brute force version 2
// ************************************************************************
double MCMCAnalyzer::createDiscrepancyFunctions(int nInputs, int nOutputs,
                       double *lower, double *upper, psIVector &rsIndices, 
                       psVector &rsValues, psIVector &dParams,
                       int dnInputs, int dnSamples, psMatrix &dSamInputs,
                       psMatrix &dSamMeans,
                       PsuadeData *dataPtr, 
                       psVector &discFuncConstantMeans,
                       psVector &discFuncConstantStds,
                       FuncApprox **faPtrs,
                       FuncApprox ***faPtrs1_in, int printLevel)
{
  int    ii, kk, ii2, cnt, askFlag=0, *states=NULL, iOne=1;
  int    *ExpSamStates, ind, ExpNSamples, dfaType, dnPerDim=16;
  double *dOneSample, expdata, simdata, *tLowers=NULL, *tUppers=NULL;
  double *ExpSamOutputs, *ExpSamInputs, *tSamInputs, *settings;
  char   charString[1000], lineIn[10000];
  FuncApprox **faPtrs1;
  pData  qData;
  PsuadeData *dPtr;

  ExpNSamples   = dnSamples;
  ExpSamInputs  = new double[ExpNSamples*nInputs];
  ExpSamOutputs = new double[ExpNSamples];
  ExpSamStates  = new int[ExpNSamples];
  discFuncConstantMeans.setLength(nOutputs);
  discFuncConstantStds.setLength(nOutputs);
  for (ii2 = 0; ii2 < nOutputs; ii2++)
  {
    discFuncConstantMeans[ii2] = PSUADE_UNDEFINED;
    discFuncConstantStds[ii2] = PSUADE_UNDEFINED;
  }
  if (dataPtr != NULL) dataPtr->getParameter("input_names", qData);

  printOutTS(PL_INFO,
       "*** SELECT RESPONSE SURFACE TYPE FOR DISCREPANCY FUNCTION:\n");
  dfaType = -1;
  while (dfaType < 0 || dfaType >= PSUADE_NUM_RS)
  {
    writeFAInfo(-1);
    sprintf(charString, "===> Enter your choice : ");
    dfaType = getInt(0, PSUADE_NUM_RS-1, charString);
  }

  settings = new double[nInputs];
  for (ii2 = 0; ii2 < nInputs; ii2++)
  {
    settings[ii2] = 0.5*(lower[ii2] + upper[ii2]);
    if (rsIndices.length() > 0 && rsIndices[ii2] < 0)
      settings[ii2] = rsValues[ii2];
  }

  faPtrs1 = new FuncApprox*[nOutputs];
  (*faPtrs1_in) = faPtrs1;
  dOneSample = new double[nInputs];
  tSamInputs = new double[ExpNSamples*nInputs];
  tLowers = new double[nInputs];
  tUppers = new double[nInputs];
  dPtr    = new PsuadeData();
  char       **iNames;
  for (ii = 0; ii < nOutputs; ii++)
  {
    for (kk = 0; kk < dnSamples; kk++)
    {
      cnt = 0;
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
        if (dParams.length() > 0 && dParams[ii2] == 1)
        {
          dOneSample[ii2] = dSamInputs.getEntry(kk,cnt);
          cnt++;
        }
        else dOneSample[ii2] = settings[ii2];
      }

      simdata = 0.0;
      if (psAnaExpertMode_ == 1 && askFlag == 0)
      {
        printOutTS(PL_INFO,
             "To create discrepancy functions, the calibration\n");
        printOutTS(PL_INFO,
             "parameters need to be set to some nominal values.\n");
        printOutTS(PL_INFO,
             "You can choose the nominal values, or it will be\n");
        printOutTS(PL_INFO,
             "set to the input means or mid points of the ranges.\n");
        printf( "Set nomininal values yourself ? (y or n) ");
        scanf("%s", charString);
        fgets(lineIn,1000,stdin);
        if (charString[0] == 'y')
        {
          for (ii2 = 0; ii2 < nInputs; ii2++)
          {
            if ((rsIndices.length() == 0 || rsIndices[ii2] >= 0) &&
                (dParams.length() == 0 || dParams[ii2] == 0))
            {
              printOutTS(PL_INFO,
                  "Input %d has lower and upper bounds = %e %e\n",
                  ii2+1, lower[ii2], upper[ii2]);
              sprintf(charString, 
                      "Nominal value for input %d : ",ii2+1);
              dOneSample[ii2] = getDouble(charString);
              settings[ii2]   = dOneSample[ii2];
            }
          }
        }
        askFlag = 1;
      }
      simdata = faPtrs[ii]->evaluatePoint(dOneSample);
      expdata = dSamMeans.getEntry(kk, ii);

      if (printLevel >= 4)
      {
        printOutTS(PL_INFO, 
             "Experiment %4d (out of %d) : ",kk+1,dnSamples);
        for (ii2 = 0; ii2 < nInputs; ii2++)
           printOutTS(PL_INFO,
                "Input %7d = %12.4e ",ii2+1,dOneSample[ii2]);
        printOutTS(PL_INFO, 
             "simuation, experimental data = %12.4e %12.4e\n",
                   simdata, expdata);
      }
      ExpSamOutputs[kk] = expdata - simdata;
    }

    if (dnInputs > 0)
    {
      for (kk = 0; kk < ExpNSamples; kk++)
        for (ii2 = 0; ii2 < dnInputs; ii2++)
          tSamInputs[kk*dnInputs+ii2] = dSamInputs.getEntry(kk,ii2);
      iNames = new char*[dnInputs];
      cnt = 0;
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
        if (dParams[ii2] == 1)
        {
          iNames[cnt] = new char[100];
          tLowers[cnt] = lower[ii2];
          tUppers[cnt] = upper[ii2];
          if (qData.strArray_ == NULL)
               sprintf(iNames[cnt], "X%d", ii2+1);
          else strcpy(iNames[cnt], qData.strArray_[ii2]);
          cnt++;
        }
      }
      dPtr->updateInputSection(ExpNSamples,dnInputs,NULL,tLowers,
                     tUppers,tSamInputs, iNames, NULL,NULL,NULL,NULL);
      for (ii2 = 0; ii2 < dnInputs; ii2++) delete [] iNames[ii2];
      delete [] iNames;
    }
    else
    {
      iNames = new char*[1];
      iNames[0] = new char[100];
      sprintf(iNames[0], "X0");
      for (ii2 = 0; ii2 < ExpNSamples; ii2++) tSamInputs[ii2] = 0.5;
      tLowers[0] = 0.0;
      tUppers[0] = 1.0;
      dPtr->updateInputSection(ExpNSamples,iOne,NULL,tLowers,tUppers,
                          tSamInputs, iNames, NULL,NULL,NULL,NULL);
      delete [] iNames[0];
      delete [] iNames;
    }
    states = new int[ExpNSamples];
    for (kk = 0; kk < ExpNSamples; kk++) states[kk] = 1;
    iNames = new char*[1];
    iNames[0] = new char[100];
    sprintf(iNames[0], "Y%d", ii+1);
    dPtr->updateOutputSection(ExpNSamples,iOne,ExpSamOutputs,
                                 states,iNames);
    delete [] states;
    delete [] iNames[0];
    delete [] iNames;
    dPtr->updateMethodSection(PSUADE_SAMP_MC, ExpNSamples, 1, -1, -1);
    sprintf(charString, "psDiscrepancyModel%d", ii+1);
    dPtr->writePsuadeFile(charString, 0);

    printOutTS(PL_INFO,
         "Creating discrepancy response surface for output %d\n",ii+1);
    faPtrs1[ii] = NULL;
    if (dnInputs > 0 && dnSamples > 1)
    {
      faPtrs1[ii] = genFA(dfaType,dnInputs,iOne,ExpNSamples);
      if (faPtrs1[ii] == NULL)
      {
        printOutTS(PL_ERROR,
           "MCMC ERROR: cannot create discrepancy func for output %d.\n",
             ii+1);
        return -1.0;
      }
    }
    if (faPtrs1[ii] != NULL)
    {
      faPtrs1[ii]->setNPtsPerDim(dnPerDim);
      faPtrs1[ii]->setBounds(lower, upper);
      faPtrs1[ii]->setOutputLevel(0);
      faPtrs1[ii]->initialize(tSamInputs,ExpSamOutputs);
    }
    else
    {
      discFuncConstantMeans[ii] = 0.0;
      for (kk = 0; kk < ExpNSamples; kk++)
        discFuncConstantMeans[ii] += ExpSamOutputs[kk];
      discFuncConstantMeans[ii] /= (double) ExpNSamples;
      discFuncConstantStds[ii] = 0.0;
      for (kk = 0; kk < ExpNSamples; kk++)
        discFuncConstantStds[ii] +=
              pow(ExpSamOutputs[kk]-discFuncConstantMeans[ii],2.0);
      discFuncConstantStds[ii] = 
               sqrt(discFuncConstantStds[ii]/ExpNSamples);
    }
    printOutTS(PL_INFO,
         "Discrepancy response surface for output %d created.\n",ii+1);
  }
  delete [] ExpSamInputs;
  delete [] dOneSample;
  delete [] ExpSamOutputs;
  delete [] settings;
  delete [] tSamInputs;
  delete [] tLowers;
  delete [] tUppers;
  delete dPtr;
  return 0.0;
}

// ************************************************************************
// perform MCMC-like analysis (brute force): for direct call externally
// assume some input parameters are certain and some may be fixed uncertain
// parameters.
// ------------------------------------------------------------------------
double MCMCAnalyzer::analyzeDirect(McmcData &mdata)
{
  int    ii, ii2, jj, kk, kk2, status, iZero=0, cnt, combFlag=0;
  double ddata;

  //**/ ---------------------------------------------------------------
  //**/ nFUInputs - number of fixed uncertain parameters
  //**/ VecFUInputs - a list of fixed certain parameters
  //**/ MatFUInpSample - sample for the fixed uncertain parameters
  //**/ nCUInputs - number of calibration parameters
  //**/ VecCUInputs - a list of calibration parameters
  //**/ MatPriorSample - prior sample for the calibration parameters
  //**/ MatExpInputs - experimental sample inputs
  //**/ MatExpMeans - experimental sample output means
  //**/ MatExpStds - experimental sample output standard deviation
  //**/ MatPostSample - posterior sample
  //**/ ---------------------------------------------------------------

  //**/ ---------------------------------------------------------------
  //**/ error checking
  //**/ ---------------------------------------------------------------
  int nFUInputs = mdata.VecFUInputs_.length();
  if (nFUInputs > 0)
  {
    for (ii = 0; ii < nFUInputs; ii++)
    {
      if (mdata.VecFUInputs_[ii] < 0 || 
          mdata.VecFUInputs_[ii] >= mdata.nInputs_)
      {
        printf("MCMC ERROR: FixedUInput list problem.\n");
        printf("            Input = %d is invalid\n",
             mdata.VecFUInputs_[ii]+1);
        exit(1);
      }
    }
    if (mdata.MatFUInpSample_.ncols() != nFUInputs)
    {
      printf("MCMC ERROR: FixedUInput sample problem.\n");
      printf("            Sample matrix has = %d columns\n",
           mdata.MatFUInpSample_.ncols());
      printf("            Expected  = %d\n", nFUInputs);
      exit(1);
    }
  }

  int nCUInputs = mdata.VecCUInputs_.length();
  for (ii = 0; ii < nCUInputs; ii++)
  {
    if (mdata.VecCUInputs_[ii] < 0 || 
        mdata.VecCUInputs_[ii] >= mdata.nInputs_)
    {
      printf("MCMC ERROR: Calibration input list problem.\n");
      printf("            Input = %d is invalid\n",
             mdata.VecCUInputs_[ii]+1);
      exit(1);
    }
  }
  if (mdata.MatPriorSample_.ncols() != nCUInputs)
  {
    printf("MCMC ERROR: Prior sample problem.\n");
    printf("            Prior sample has = %d columns\n",
           mdata.MatPriorSample_.ncols());
    printf("            Expected  = %d\n", nCUInputs);
    exit(1);
  }

  int nDInputs = mdata.nInputs_ - nFUInputs - nCUInputs;
  int nExperiments = mdata.MatExpMeans_.nrows();
  if (mdata.MatExpInputs_.ncols() != nDInputs)
  {
    printf("MCMC ERROR: experiment sample problem.\n");
    printf("            Experiment matrix has = %d columns\n",
           mdata.MatExpInputs_.ncols());
    printf("            Expected  = %d\n", nDInputs);
    exit(1);
  }
  if (mdata.MatExpInputs_.ncols() > 0 &&
      mdata.MatExpInputs_.nrows() != mdata.MatExpMeans_.nrows())
  {
    printf("MCMC ERROR: experiment sample size and means mismatch\n");
    printf("            Experiment input matrix has = %d rows\n",
           mdata.MatExpInputs_.nrows());
    printf("            Experiment means matrix has = %d rows\n",
           mdata.MatExpMeans_.nrows());
    exit(1);
  }
  if (mdata.MatExpMeans_.nrows() != mdata.MatExpStds_.nrows())
  {
    printf("MCMC ERROR: MatExpMeans mismatch with MatExpStds\n");
    printf("            Experiment means matrix has = %d rows\n",
           mdata.MatExpMeans_.nrows());
    printf("            Experiment stdvs matrix has = %d rows\n",
           mdata.MatExpStds_.nrows());
    exit(1);
  }
  if (mdata.MatExpMeans_.ncols() != mdata.MatExpStds_.ncols())
  {
    printf("MCMC ERROR: MatExpMeans mismatch with MatExpStds\n");
    printf("            Experiment means matrix has = %d columns\n",
           mdata.MatExpMeans_.ncols());
    printf("            Experiment stdvs matrix has = %d columns\n",
           mdata.MatExpStds_.ncols());
    exit(1);
  }

  int nOutputs = mdata.nOutputs_;
  if (mdata.MatExpMeans_.ncols() != nOutputs)
  {
    printf("MCMC ERROR: experiment sample output problem.\n");
    printf("            Experiment means matrix has = %d columns\n",
           mdata.MatExpMeans_.ncols());
    printf("            nOutputs = %d\n", nOutputs_);
    exit(1);
  }
  if (mdata.VecLowerB_.length() != mdata.nInputs_)
  {
    printf("MCMC ERROR: input lower bound problem.\n");
    printf("            Lower bound vector has length %d\n",
           mdata.VecLowerB_.length());
    printf("            nInputs = %d\n", mdata.nInputs_);
    exit(1);
  }
  if (mdata.VecUpperB_.length() != mdata.nInputs_)
  {
    printf("MCMC ERROR: input upper bound problem.\n");
    printf("            Upper bound vector has length %d\n",
           mdata.VecUpperB_.length());
    printf("            nInputs = %d\n", mdata.nInputs_);
    exit(1);
  }

  //**/ ---------------------------------------------------------------
  //**/ create vecDesParams and vecInpTypes for later use
  //**/ vecInpTypes[ii]  == 1000 if ii is a fixed uncertain parameter
  //**/ vecInpTypes[ii]  == 2000 if ii is a design parameter
  //**/ vecInpTypes[ii]  == ii    if ii is a design parameter
  //**/ vecDesParams[ii] == 1 if ii is a design parameter
  //**/ ---------------------------------------------------------------
  psIVector vecInpTypes, vecDesParams;
  vecInpTypes.setLength(mdata.nInputs_);
  for (ii = 0; ii < mdata.nInputs_; ii++) vecInpTypes[ii] = ii;
  cnt = 0;
  for (ii = 0; ii < nFUInputs; ii++)
    vecInpTypes[mdata.VecFUInputs_[ii]] = 1000;
  cnt = 0;
  for (ii = 0; ii < nCUInputs; ii++)
  {
    kk = mdata.VecCUInputs_[ii];
    if (vecInpTypes[kk] == 1000)
    {
      printf("MCMC ERROR: input %d both calibration/uncertain input\n",
             kk+1);
      exit(1);
    }
    vecInpTypes[kk] = 2000;
  }
  vecDesParams.setLength(mdata.nInputs_);
  for (ii = 0; ii < mdata.nInputs_; ii++) 
    if (vecInpTypes[ii] < 1000) vecDesParams[ii] = 1;

  //**/ ---------------------------------------------------------------
  //**/ FUNumSampletoUse - if large FU sample, select a fixed subset
  //**/ ---------------------------------------------------------------
  int FUNumSampleToUse=1;
  FUNumSampleToUse = mdata.MatFUInpSample_.nrows();
  if (FUNumSampleToUse == 0) FUNumSampleToUse = 1;
  if (FUNumSampleToUse > 10000) FUNumSampleToUse = 10000;

  //**/ ---------------------------------------------------------------
  //**/ create response surface for use in computing likelihood
  //**/ ==> faPtrs
  //**/ ---------------------------------------------------------------
  FuncApprox **faPtrs = new FuncApprox*[mdata.nOutputs_];
  psVector vecYT;
  vecYT.setLength(mdata.nSamples_);
  for (ii = 0; ii < mdata.nOutputs_; ii++)
  {
    faPtrs[ii] = genFA(mdata.faType_,mdata.nInputs_,iZero,
                       mdata.nSamples_);
    faPtrs[ii]->setNPtsPerDim(16);
    faPtrs[ii]->setBounds(mdata.VecLowerB_.getDVector(),
                          mdata.VecUpperB_.getDVector());
    faPtrs[ii]->setOutputLevel(0);
    for (kk = 0; kk < mdata.nSamples_; kk++) 
      vecYT[kk] = mdata.VecSamOutputs_[kk*mdata.nOutputs_+ii];

    status = faPtrs[ii]->initialize(mdata.VecSamInputs_.getDVector(), 
                                    vecYT.getDVector());
    if (status != 0)
    {
      printOutTS(PL_ERROR,
           "MCMC ERROR: Unable to create response surface.\n");
      printOutTS(PL_ERROR,"            Consult PSUADE developers.\n");
      for (kk = 0; kk < ii; kk++) delete [] faPtrs[kk];
      delete [] faPtrs;
      return PSUADE_UNDEFINED;
    }
  }

  //**/ ---------------------------------------------------------------
  //**/ ==> max MCMC Samples
  //**/ ---------------------------------------------------------------
  int maxMCMCSamples = 1000000;
  if (mdata.nInputs_ >= 10) maxMCMCSamples = 2000000;

  //**/ ---------------------------------------------------------------
  //**/ option to add discrepancy function 
  //**/ if no design parameter, the best model form is a constant
  //**/ ---------------------------------------------------------------
  int modelFormConst=0, genPosteriors = 1;
  if (mdata.addDiscrepancy_ == 1 && nDInputs == 0) modelFormConst = 1;
  if (mdata.addDiscrepancy_ == 1 && nExperiments == 1) 
    modelFormConst = 1;

  //**/ ---------------------------------------------------------------
  //**/ create discrepancy function, if desired 
  //**/ ==> vecDiscSamOuts or
  //**/ ==> vecDiscFuncConstMeans, vecDistFuncConstStds
  //**/ ---------------------------------------------------------------
  psVector vecDiscFuncConstMeans, vecDiscFuncConstStds, vecDiscSamOuts;
  if (mdata.addDiscrepancy_ == 1)
  {
    double dstatus;
    psIVector vecINULL;
    psVector  vecDNULL;
    dstatus = createDiscrepancyFunctions2(mdata.nInputs_,mdata.nOutputs_,
              mdata.VecLowerB_.getDVector(), mdata.VecUpperB_.getDVector(),
              vecINULL,vecDNULL,vecDesParams,nDInputs,nExperiments,
              mdata.MatExpInputs_, mdata.MatExpMeans_,NULL,
              vecDiscFuncConstMeans,vecDiscFuncConstStds,vecDiscSamOuts,
              faPtrs, iZero, modelFormConst);
    if (dstatus < 0) return PSUADE_UNDEFINED;
  }

  //**/ ---------------------------------------------------------------
  //  set up for inference: vecXmax, vecRanges, vecMeans, vecSigmas
  //                        vecMostLikelyInputs, vecMostLikelyOutputs
  //**/ ---------------------------------------------------------------
  //**/ storage allocation
  psVector vecXmax, vecRanges;
  vecRanges.setLength(mdata.nInputs_);
  for (ii = 0; ii < mdata.nInputs_; ii++) 
    vecRanges[ii] = 1.0/(mdata.VecUpperB_[ii]-mdata.VecLowerB_[ii]);
  vecXmax.setLength(mdata.nInputs_);

  //**/ ---------------------------------------------------------------
  //**/ MCMC parameters (samInc, maxMCMCSamples)
  //**/ ---------------------------------------------------------------
  int nPriorSamples = mdata.MatPriorSample_.nrows();
  int nIncrements = 20;
  nPriorSamples = nPriorSamples / nIncrements * nIncrements;
  maxMCMCSamples = FUNumSampleToUse * nPriorSamples;
  
  //**/ ---------------------------------------------------------------
  //**/ generate a large sample ==> vecPriorSample
  //**/ ---------------------------------------------------------------
  psVector vecPriorSample; 
  mdata.MatPriorSample_.convert2Vector(vecPriorSample);

  //**/ ---------------------------------------------------------------
  //  duplicate the generated sample FUNumSampleToUse-1 times
  //**/ so inferenceSamIns has a sample of size nPriorSample
  //**/ replicated FUNumSampleToUse times
  //**/ ---------------------------------------------------------------
  psVector vecInfSamInp, vecInfSamOut, vecInfTmpOut;
  vecInfSamInp.setLength(maxMCMCSamples * mdata.nInputs_);
  double *inferenceSamIns = vecInfSamInp.getDVector();
  vecInfSamOut.setLength(maxMCMCSamples);
  double *inferenceSamOut = vecInfSamOut.getDVector();
  for (ii = 0; ii < nPriorSamples; ii++)
  {
    for (jj = 0; jj < FUNumSampleToUse; jj++)
    {
      cnt = (ii * FUNumSampleToUse + jj) * mdata.nInputs_;
      for (kk = 0; kk < mdata.VecCUInputs_.length(); kk++)
      {
        kk2 = mdata.VecCUInputs_[kk]; 
        inferenceSamIns[cnt+kk2] = 
          vecPriorSample[ii*mdata.VecCUInputs_.length()+kk]; 
      }
    }
  }
  vecPriorSample.clean();

  //**/ ---------------------------------------------------------------
  //**/ if there are uncertain parameters, create a random vector
  //**/ that selects a fixed sub-sample for the uncertain parameters,
  //**/ and fill in the slots in inferenceSamIns for uncertain parameters
  //**/ ---------------------------------------------------------------
  psIVector randIVec;
  if (mdata.MatFUInpSample_.ncols() > 0)
  {
    randIVec.setLength(FUNumSampleToUse);
    for (int mm = 0; mm < FUNumSampleToUse; mm++)
    {
      if (FUNumSampleToUse >= mdata.MatFUInpSample_.nrows())
           ii = mm;
      else ii = PSUADE_rand() % mdata.MatFUInpSample_.nrows();
      randIVec[mm] = ii;
    }
    for (ii = 0; ii < maxMCMCSamples; ii+=FUNumSampleToUse)
    {
      for (kk = 0; kk < nFUInputs; kk++)
      {
        cnt = mdata.VecFUInputs_[kk];
        for (jj = 0; jj < FUNumSampleToUse; jj++)
        {
          ii2 = randIVec[jj];     
          ddata = mdata.MatFUInpSample_.getEntry(ii2,kk);
          inferenceSamIns[(ii+jj)*mdata.nInputs_+cnt] = ddata;
        }
      }
    }
  }

  //**/ ---------------------------------------------------------------
  //  allocate storage for inference
  //**/ ---------------------------------------------------------------
  psVector vecXS, vecYS, vecSamStd, vecYDesign, vecYDesStd;
  cnt = maxMCMCSamples / nIncrements * nExperiments;
  vecXS.setLength(cnt * mdata.nInputs_);
  double *XSample = vecXS.getDVector();
  vecYS.setLength(cnt * mdata.nOutputs_);
  double *YSample = vecYS.getDVector();
  vecSamStd.setLength(cnt * mdata.nOutputs_);
  double *YSamStd = vecSamStd.getDVector();
  vecYDesign.setLength(cnt * mdata.nOutputs_);
  double *YDesign = vecYDesign.getDVector();
  vecYDesStd.setLength(cnt * mdata.nOutputs_);
  double *YDesStd = vecYDesStd.getDVector();

  //**/ ---------------------------------------------------------------
  //  perform inference
  //**/ ---------------------------------------------------------------
  int    ss, mm, runSize, index, ind1, ind2;
  double stdev, stdv2, YT, YT1, YT2;

  if (mdata.printLevel_ > 0)
    printOutTS(PL_INFO, "MCMC_BFS Inference begins ... \n");
  fflush(stdout);

  int samInc = maxMCMCSamples / (FUNumSampleToUse * nIncrements);
  int curMCMCSample = 0;
  while (curMCMCSample < maxMCMCSamples)
  {
    //**/ have one more loop for the uncertain parameters
    for (mm = 0; mm < samInc; mm++)
    {
      cnt = 0;
      index = mm * nExperiments * FUNumSampleToUse;
      //**/ first fill in all XSample slots
      for (ii2 = 0; ii2 < mdata.nInputs_; ii2++)
      {
        //**/ load prior samples from inferenceSamIns
        if (vecDesParams[ii2] == 0)
        {
          for (kk2 = 0; kk2 < nExperiments; kk2++)
          {
            for (jj = 0; jj < FUNumSampleToUse; jj++)
            {
              ind1 = (index+kk2*FUNumSampleToUse+jj)*mdata.nInputs_+ii2; 
              ind2 = (curMCMCSample+mm*FUNumSampleToUse+jj)*mdata.nInputs_+ii2; 
              XSample[ind1] = inferenceSamIns[ind2];
            }
          }
        }
        //**/ load design parameters
        if (vecDesParams[ii2] == 1)
        {
          for (kk2 = 0; kk2 < nExperiments; kk2++)
          {
            for (jj = 0; jj < FUNumSampleToUse; jj++)
            {
              ind1 = (index+kk2*FUNumSampleToUse+jj)*mdata.nInputs_+ii2; 
              XSample[ind1] = mdata.MatExpInputs_.getEntry(kk2, cnt);
            }
          }
          cnt++;
        }
      }
    }

    //**/ set design outputs and stds (discrepancy) to zero 
    //**/ since they will be used but may not be set later
    int totCnt = samInc * nExperiments * mdata.nOutputs_;
    for (jj = 0; jj < totCnt; jj++)
    {
      YDesign[jj] = YDesStd[jj] = 0.0;
      YSample[jj] = YSamStd[jj] = 0.0;
    }
    //**/ run XGuessS through response surfaces
    //**/ ==> YDesign, YDesStd, YSample, YSamStd
    runSize = samInc * nExperiments;
    for (ii2 = 0; ii2 < mdata.nOutputs_; ii2++)
    {
      //**/ case 1: if RS error is requested
      if (mdata.useRSUncertainties_ == 1)
      {
        faPtrs[ii2]->evaluatePointFuzzy(runSize, 
                             XSample,&YSample[ii2*runSize],
                             &YSamStd[ii2*runSize]);
        //**/ add discrepancy function, if available
        if (mdata.addDiscrepancy_ == 1 && modelFormConst == 0)
        {
          for (mm = 0; mm < samInc; mm++)
          {
            index = mm * nExperiments * FUNumSampleToUse;
            for (kk2 = 0; kk2 < nExperiments; kk2++)
            {
              for (jj = 0; jj < FUNumSampleToUse; jj++)
              {
                ind1 = ii2*runSize+index+kk2*FUNumSampleToUse+jj;
                YDesign[ind1] = vecDiscSamOuts[ii2*nExperiments+kk2];
                YDesStd[ind1] = 0.0;
              }
            }
          }
        }
        else if (mdata.addDiscrepancy_ == 1 && modelFormConst == 1 &&
                 vecDiscFuncConstMeans.length() > 0 &&
                 vecDiscFuncConstMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < runSize; kk2++)
          {
            YDesign[ii2*runSize+kk2] = vecDiscFuncConstMeans[ii2];
            YDesStd[ii2*runSize+kk2] = 0.0;
          }
        }
      }
      //**/ case 2: if RS error is to be turned off
      else
      {
        faPtrs[ii2]->evaluatePoint(runSize,XSample,
                                   &YSample[ii2*runSize]);
        //**/ add discrepancy function, if available
        if (mdata.addDiscrepancy_ == 1 && modelFormConst == 0)
        {
          for (mm = 0; mm < samInc; mm++)
          {
            index = mm * nExperiments * FUNumSampleToUse;
            for (kk2 = 0; kk2 < nExperiments; kk2++)
            {
              for (jj = 0; jj < FUNumSampleToUse; jj++)
              {
                ind1 = ii2*runSize+index+kk2*FUNumSampleToUse+jj;
                YDesign[ind1] = vecDiscSamOuts[ii2*nExperiments+kk2];
              }
            }
          }
        }
        else if (mdata.addDiscrepancy_ == 1 && modelFormConst == 1 &&
                 vecDiscFuncConstMeans.length() > 0 &&
                 vecDiscFuncConstMeans[ii2] != PSUADE_UNDEFINED)
        {
          for (kk2 = 0; kk2 < runSize; kk2++)
            YDesign[ii2*runSize+kk2] = vecDiscFuncConstMeans[ii2];
        }
        for (kk2 = 0; kk2 < runSize; kk2++)
          YSamStd[ii2*runSize+kk2] = YDesStd[ii2*runSize+kk2] = 0.0;
      }
    }

    //**/ compute vecXDist
    for (ii = 0; ii < samInc; ii++)
    {
      index = ii * FUNumSampleToUse;
      for (jj = 0; jj < FUNumSampleToUse; jj++)
      {
        inferenceSamOut[curMCMCSample+index+jj] = 0.0;
        for (ii2 = 0; ii2 < mdata.nOutputs_; ii2++)
        {
          for (kk2 = 0; kk2 < nExperiments; kk2++)
          {
            ind1 = ii2*runSize+index*nExperiments+kk2*FUNumSampleToUse+jj;
            YT1 = YSample[ind1] + YDesign[ind1];
            stdev = YSamStd[ind1] + YDesign[ind1];
            stdev = YSamStd[ind1];
            stdv2 = YDesStd[ind1];
            YT2 = pow((YT1-mdata.MatExpMeans_.getEntry(kk2,ii2)),2.0) /
                 (pow(mdata.MatExpStds_.getEntry(kk2,ii2),2.0) +
                  stdev*stdev + stdv2*stdv2);
            inferenceSamOut[curMCMCSample+index+jj] += YT2;
          }
        }
        //**/ more correct?
        if (combFlag == 1)
          inferenceSamOut[curMCMCSample+index+jj] /= (double) nExperiments;
      }
    }
    curMCMCSample += (samInc * FUNumSampleToUse);
  }

  //**/ ---------------------------------------------------------------
  //**/  generate the posterior sample 
  //**/ ---------------------------------------------------------------
  int imax=-1;
  if (genPosteriors == 1)
  {
    int    maxPostSam=nPriorSamples;
    double dmax = 0.0;
    for (ss = 0; ss < maxMCMCSamples; ss++) 
    {
      inferenceSamOut[ss] = exp(-0.5*inferenceSamOut[ss]);
      if (inferenceSamOut[ss] > dmax) 
      {
        imax = ss;
        dmax = inferenceSamOut[ss];
      }
    }
    if (imax != -1)
    {
      if (mdata.printLevel_ > 1)
      {
        for (ii = 0; ii < mdata.VecCUInputs_.length(); ii++) 
        {
          kk = mdata.VecCUInputs_[ii];
          printf("Input %4d: xmax = %e\n", kk+1, 
                 inferenceSamIns[imax*mdata.nInputs_+kk]); 
        }
      }
    }
    if (dmax == 0)
    {
      printOutTS(PL_ERROR,
         "MCMC_BFS: ERROR encountered in posterior sample generation.\n");
      printOutTS(PL_ERROR,
         "          Maybe due to empty posterior. Use larger data errors.\n");
      exit(1);
    }
    else 
    {
      for (ss = 0; ss < maxMCMCSamples; ss++) 
        inferenceSamOut[ss] = inferenceSamOut[ss] / dmax;
      //**/ this format is to be compatible with ProbMatrix
      mdata.MatPostSample_.setFormat(PS_MAT2D);
      mdata.MatPostSample_.setDim(maxPostSam,mdata.VecCUInputs_.length());
      mdata.VecPostLikelihoods_.setLength(maxPostSam);
      int nValid=0;
      psIVector vecValid;
      vecValid.setLength(maxMCMCSamples);
      for (ss = 0; ss < maxMCMCSamples; ss++) 
        if (inferenceSamOut[ss] > 1e-5) vecValid[nValid++] = ss;
      cnt = 0;
      while (cnt < maxPostSam)
      {
        kk2 = PSUADE_rand() % nValid;
        kk2 = vecValid[kk2];
        ddata = PSUADE_drand();
        if (inferenceSamOut[kk2] > ddata)
        {
          for (jj = 0; jj < mdata.VecCUInputs_.length(); jj++)
          {
            ind1 = mdata.VecCUInputs_[jj];
            ind2 = kk2 * mdata.nInputs_ + ind1;
            mdata.MatPostSample_.setEntry(cnt,jj,inferenceSamIns[ind2]);
            mdata.VecPostLikelihoods_[cnt] = (double) inferenceSamOut[kk2];
          }
          cnt++;
        }
      }
    }
  }

  //**/ ---------------------------------------------------------------
  // clean up
  //**/ ---------------------------------------------------------------
  if (faPtrs != NULL)
  {
    for (ii = 0; ii < mdata.nOutputs_; ii++)
      if (faPtrs[ii] != NULL) delete faPtrs[ii];
    delete [] faPtrs;
  }
  return 0.0;
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int MCMCAnalyzer::get_nInputs()
{
  return nInputs_;
}
int MCMCAnalyzer::get_nOutputs()
{
  return nOutputs_;
}
double *MCMCAnalyzer::get_means()
{
  int    ii;
  double *retVal = NULL;
  if (means_)
  {
    retVal = new double[nInputs_];
    for (ii = 0; ii < nInputs_; ii++) retVal[ii] = means_[ii];
  }
  return retVal;
}
double *MCMCAnalyzer::get_mostLikelyInput()
{
  int    ii;
  double *retVal = NULL;
  if (mostLikelyInput_)
  {
    retVal = new double[nInputs_];
    checkAllocate(retVal, "retVal in MCMC::get_mostLikelyInput");
    for (ii = 0; ii < nInputs_; ii++) retVal[ii] = mostLikelyInput_[ii];
  }
  return retVal;
}
double *MCMCAnalyzer::get_mostLikelyOutput()
{
  int    ii;
  double *retVal = NULL;
  if (mostLikelyOutput_)
  {
    retVal = new double[nOutputs_];
    checkAllocate(retVal, "retVal in MCMC::get_mostLikelyOutput");
    for (ii = 0; ii < nOutputs_; ii++) 
      retVal[ii] = mostLikelyOutput_[ii];
  }
  return retVal;
}

// ************************************************************************
// create discrepancy function (either constant or some RS)
// ************************************************************************
double MCMCAnalyzer::createDiscrepancyFunctions2(int nInputs, int nOutputs,
                       double *lower, double *upper, psIVector &rsIndices, 
                       psVector &rsValues, psIVector &dParams,
                       int dnInputs, int dnSamples, psMatrix &dSamInputs,
                       psMatrix &dSamMeans, PsuadeData *dataPtr, 
                       psVector &vecDiscFuncConstMeans,
                       psVector &vecDistFuncConstStds,
                       psVector &vecDiscSamOuts, FuncApprox **faPtrs,
                       int printLevel, int constFlag)
{
  int    ii, kk, ii2, cnt, askFlag=0, iOne=1, pOrder, rstype, status;
  double expdata, simdata, ddata;
  char   charString[1000], lineIn[10000];
  pData  qData;
  psVector   vecSettings, vecXT, vecYT;
  PsuadeData *dPtr;
  FuncApprox *localFaPtr;

  //**/ ---------------------------------------------------------------
  //**/ allocate space for discrepancy data ==> vecDiscSamOuts
  //**/ and initialize the constant discrepancy function
  //**/ qData is for fetching input names for storing discrepancy
  //**/ ---------------------------------------------------------------
  int ExpNSamples = dnSamples;
  vecDiscSamOuts.setLength(ExpNSamples*nOutputs);
  double *discSamOutputs = vecDiscSamOuts.getDVector();
  vecDiscFuncConstMeans.setLength(nOutputs);
  vecDistFuncConstStds.setLength(nOutputs);
  for (ii2 = 0; ii2 < nOutputs; ii2++)
  {
    vecDiscFuncConstMeans[ii2] = PSUADE_UNDEFINED;
    vecDistFuncConstStds[ii2]  = PSUADE_UNDEFINED;
  }
  if (dataPtr != NULL) dataPtr->getParameter("input_names", qData);

  //**/ ---------------------------------------------------------------
  //**/ set discrepancy function calibration parameter default 
  //**/ (to the middle of the range). Also set fixed parameters
  //**/ ==> vecSsettings
  //**/ ---------------------------------------------------------------
  vecSettings.setLength(nInputs);
  for (ii2 = 0; ii2 < nInputs; ii2++)
  {
    vecSettings[ii2] = 0.5*(lower[ii2] + upper[ii2]);
    if (rsIndices.length() > 0 && rsIndices[ii2] < 0)
      vecSettings[ii2] = rsValues[ii2];
  }

  //**/ ---------------------------------------------------------------
  //**/ create and store discrepancy sample
  //**/ (same as experimental sample but with adjustment)
  //**/ thus, dSamMeans may be changed (if RS is used for discrepancy)
  //**/ ---------------------------------------------------------------
  psVector vecOneSample, vecLowers, vecUppers, vecSamIns;
  vecOneSample.setLength(nInputs);
  vecSamIns.setLength(ExpNSamples*nInputs);
  vecLowers.setLength(nInputs);
  vecUppers.setLength(nInputs);
  dPtr = new PsuadeData();
  psStrings Xnames, Ynames;
  for (ii = 0; ii < nOutputs; ii++)
  {
    //**/ for each output, compute discrepancies between
    //**/ simulation (evaluated at the experimental points)
    //**/ and experiment for every experiment
    for (kk = 0; kk < ExpNSamples; kk++)
    {
      //**/ inject the design, calibration default and fixed 
      //**/ parameter values (dSamInputs has design parameters
      //**/ for each of the ExpNSamples experiments)
      cnt = 0;
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
        if (dParams.length() > 0 && dParams[ii2] == 1)
        {
          vecOneSample[ii2] = dSamInputs.getEntry(kk,cnt);
          cnt++;
        }
        else vecOneSample[ii2] = vecSettings[ii2];
      }

      //**/ users can choose whether to re-set the calibration
      //**/ parameters or use the ones in vecSettings
      simdata = 0.0;
      if (psAnaExpertMode_ == 1 && askFlag == 0)
      {
        printOutTS(PL_INFO,
             "To create discrepancy functions, the calibration\n");
        printOutTS(PL_INFO,
             "parameters need to be set to some nominal values.\n");
        printOutTS(PL_INFO,
             "You can choose the nominal values, or it will be\n");
        printOutTS(PL_INFO,
             "set to the input means or mid points of the ranges.\n");
        printf( "Set nomininal values yourself ? (y or n) ");
        scanf("%s", charString);
        fgets(lineIn,1000,stdin);
        if (charString[0] == 'y')
        {
          for (ii2 = 0; ii2 < nInputs; ii2++)
          {
            //**/ if it is not fixed, and it is not a design parameter
            if ((rsIndices.length() == 0 || rsIndices[ii2] >= 0) &&
                (dParams.length() == 0 || dParams[ii2] == 0))
            {
              printOutTS(PL_INFO,
                  "Input %d has lower and upper bounds = %e %e\n",
                  ii2+1, lower[ii2], upper[ii2]);
              sprintf(charString, 
                      "Nominal value for input %d : ",ii2+1);
              vecOneSample[ii2] = getDouble(charString);
              vecSettings[ii2]  = vecOneSample[ii2];
            }
          }
        }
        askFlag = 1;
      }
      //**/ at this point, vecOneSample has been set to calibration default
      //**/ and design default (at sample kk) and at fixed values
      //**/ now evaluate it (ii-th output and kk-th designSample)
      simdata = faPtrs[ii]->evaluatePoint(vecOneSample.getDVector());
      expdata = dSamMeans.getEntry(kk, ii);

      if (printLevel >= 4)
      {
        printOutTS(PL_INFO, 
             "Experiment %4d (out of %d) : \n",kk+1,ExpNSamples);
        printOutTS(PL_INFO, "Inputs = \n");
        for (ii2 = 0; ii2 < nInputs; ii2++)
           printOutTS(PL_INFO, "%12.4e ",vecOneSample[ii2]);
        printOutTS(PL_INFO, 
             "\nSimuation/experimental data = %12.4e %12.4e\n",
             simdata, expdata);
      }
      //**/ generate the difference between experiment and simulation
      discSamOutputs[kk+ii*ExpNSamples] = expdata - simdata;
    }
    for (ii2 = 0; ii2 < nInputs; ii2++)
    {
      if ((rsIndices.length() == 0 || rsIndices[ii2] >= 0) &&
          (dParams.length() == 0 || dParams[ii2] == 0))
        printf("Nominal value for input %d = %e\n",ii2+1,vecSettings[ii2]);
    }

    //**/ Next fit the discrepancy with a function, if requested
    //**/ and then modify dSamMeans accordingly
    //**/ (This is needed if dnInputs > 0)
    if (dnInputs > 0 && dnSamples > 1)
    {
      pOrder = 0;
      if (ExpNSamples >= dnInputs+1) pOrder++;
      if (ExpNSamples >= dnInputs*(dnInputs+1)/2+1) pOrder++;
      if (ExpNSamples >= dnInputs*(dnInputs+1)*(dnInputs+2)/6+1) pOrder++;
      printf("PSUADE provides the following options discrepancy modeling:\n"); 
      printf("0. Just use (simdata-expdata) as discrepancies at experimental\n");
      printf("   points. Users are to create their own RS from the stored\n");
      printf("   discrepancy data at each experimental points.\n");
      printf("1. Gaussian process/Kriging (If TPROS is available, it will\n");
      printf("   be considered first. Otherwise, use MARS. The reason the\n");
      printf("   internal GP is not an option is because it is the same as\n");
      printf("   option 0 - exact interpolation at experimental points.)\n");
      if (pOrder > 0)
      {
        printf("2. Use polynomial function of order %d.\n",pOrder);
        sprintf(charString, "Response surface for output %d? (0-2) ",ii+1);
        rstype = getInt(0,2,charString);
      }
      else
      {
        sprintf(charString, "Response surface for output %d? (0-1) ",ii+1);
        rstype = getInt(0,1,charString);
      }
      if (rstype > 0)
      {
#ifdef HAVE_TPROS
        if      (rstype  == 1) rstype = PSUADE_RS_GP1;
#else
        if      (rstype  == 1) rstype = PSUADE_RS_KR;
#endif
        else
        {
          rstype = PSUADE_RS_REGR1;
          if      (pOrder == 2) rstype = PSUADE_RS_REGR1 + 1;
          else if (pOrder == 3) rstype = PSUADE_RS_REGR1 + 2;
        }
        localFaPtr = genFA(rstype, dnInputs, iOne, ExpNSamples);
        if (localFaPtr == NULL)
        {
          printOutTS(PL_ERROR,
             "MCMCAnalyzer ERROR: failed to create RS for discrepancy.\n");
          exit(1);
        }
        localFaPtr->setNPtsPerDim(100);
        cnt = 0;
        for (ii2 = 0; ii2 < nInputs; ii2++)
        {
          if (dParams[ii2] == 1)
          {
            vecLowers[cnt] = lower[ii2];
            vecUppers[cnt] = upper[ii2];
            cnt++;
          }
        }
        localFaPtr->setBounds(lower, upper);
        localFaPtr->setOutputLevel(-1);
        vecXT.setLength(dnInputs*ExpNSamples);
        for (kk = 0; kk < ExpNSamples; kk++)
          for (ii2 = 0; ii2 < dnInputs; ii2++)
            vecXT[kk*dnInputs+ii2] = dSamInputs.getEntry(kk,ii2);
        vecYT.setLength(ExpNSamples);
        for (kk = 0; kk < ExpNSamples; kk++)
          vecYT[kk] = discSamOutputs[kk+ii*ExpNSamples];
        status = localFaPtr->initialize(vecXT.getDVector(), vecYT.getDVector());
        if (status != 0)
        {
          printOutTS(PL_ERROR,
             "MCMCAnalyzer ERROR: crash in RS initialize for discrepancy.\n");
          exit(1);
        }
        vecYT.setLength(ExpNSamples);
        localFaPtr->evaluatePoint(ExpNSamples,vecXT.getDVector(),
                                  vecYT.getDVector());
        double errNorm = 0;
        for (kk = 0; kk < ExpNSamples; kk++)
        {
          printf("Output %d experiment %d: data = %12.6e, RS eval = %12.6e\n",
                 ii+1,kk+1,discSamOutputs[kk+ii*ExpNSamples],vecYT[kk]);
          ddata = discSamOutputs[kk+ii*ExpNSamples] - vecYT[kk];
          if (discSamOutputs[kk+ii*ExpNSamples] != 0)
               errNorm += pow(ddata/discSamOutputs[kk+ii*ExpNSamples],2.0);
          else errNorm += pow(ddata,2.0);
          discSamOutputs[kk+ii*ExpNSamples] = vecYT[kk];
        }
        printf("MSE of normalized error (true - interpolated discrepancy) = %e\n",
               sqrt(errNorm / ExpNSamples));
        printf("Please make sure this MSE is acceptable.\n");
        delete localFaPtr;
      }
    }

    //**/ store discrepancy file
    if (dnInputs > 0)
    {
      for (kk = 0; kk < ExpNSamples; kk++)
        for (ii2 = 0; ii2 < dnInputs; ii2++)
          vecSamIns[kk*dnInputs+ii2] = dSamInputs.getEntry(kk,ii2);
      Xnames.setNumStrings(dnInputs);
      cnt = 0;
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
        if (dParams[ii2] == 1)
        {
          vecLowers[cnt] = lower[ii2];
          vecUppers[cnt] = upper[ii2];
          if (qData.strArray_ == NULL)
               sprintf(charString, "X%d", ii2+1);
          else strcpy(charString, qData.strArray_[ii2]);
          Xnames.loadOneString(cnt, charString);
          cnt++;
        }
      }
      dPtr->updateInputSection(ExpNSamples,dnInputs,NULL,
                  vecLowers.getDVector(),vecUppers.getDVector(),
                  vecSamIns.getDVector(), Xnames.getStrings(), 
                  NULL,NULL,NULL,NULL);
    }
    else
    {
      Xnames.setNumStrings(1);
      strcpy(charString, "X0");
      Xnames.loadOneString(0, charString);
      //**/ set input = 0.5 (and later output = constant)
      for (ii2 = 0; ii2 < ExpNSamples; ii2++) vecSamIns[ii2] = 0.5;
      vecLowers[0] = 0.0;
      vecUppers[0] = 1.0;
      dPtr->updateInputSection(ExpNSamples,iOne,NULL,
                  vecLowers.getDVector(),vecUppers.getDVector(),
                  vecSamIns.getDVector(), Xnames.getStrings(), 
                  NULL,NULL,NULL,NULL);
    }
    psIVector vecStates;
    vecStates.setLength(ExpNSamples);
    for (kk = 0; kk < ExpNSamples; kk++) vecStates[kk] = 1;
    Ynames.setNumStrings(1);
    sprintf(charString, "Y%d", ii+1);
    Ynames.loadOneString(0, charString);
    dPtr->updateOutputSection(ExpNSamples,iOne,
                 &discSamOutputs[ii*ExpNSamples],
                 vecStates.getIVector(),Ynames.getStrings());
    dPtr->updateMethodSection(PSUADE_SAMP_MC, ExpNSamples, 1, -1, -1);
    sprintf(charString, "psDiscrepancyModel%d", ii+1);
    dPtr->writePsuadeFile(charString, 0);

    //**/ constant discrepancy function is requested
    //**/ set the constant to be the mean
    if (constFlag == 1)
    {
      vecDiscFuncConstMeans[ii] = 0.0;
      for (kk = 0; kk < ExpNSamples; kk++)
        vecDiscFuncConstMeans[ii] += discSamOutputs[ii*ExpNSamples+kk];
      vecDiscFuncConstMeans[ii] /= (double) ExpNSamples;
      vecDistFuncConstStds[ii] = 0.0;
      for (kk = 0; kk < ExpNSamples; kk++)
        vecDistFuncConstStds[ii] +=
         pow(discSamOutputs[ii*ExpNSamples+kk]-vecDiscFuncConstMeans[ii],2.0);
      vecDistFuncConstStds[ii] = 
         sqrt(vecDistFuncConstStds[ii]/ExpNSamples);
    }
  }

  //**/ return the discrepancy sample 
  delete dPtr;
  return 0.0;
}
