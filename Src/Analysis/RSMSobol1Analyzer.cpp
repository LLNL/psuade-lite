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
// Functions for the class RSMSobol1Analyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Matrix.h"
#include "Vector.h"
#include "pData.h"
#include "RSMSobol1Analyzer.h"
#include "Sampling.h"
#include "PDFManager.h"
#include "PDFNormal.h"
#include "RSConstraints.h"
#include "Psuade.h"
#include "PsuadeData.h"
#include "PsuadeConfig.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobol1Analyzer::RSMSobol1Analyzer() : Analyzer(),nInputs_(0),
                     outputMean_(0), outputStd_(0)
{
  setName("RSMSOBOL1");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobol1Analyzer::~RSMSobol1Analyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
void RSMSobol1Analyzer::analyze(int nInps, int nSamp, double *lbs, 
                                double *ubs, double *X, double *Y)
{
  aData adata;
  adata.nInputs_ = nInps;
  adata.nOutputs_ = 1;
  adata.nSamples_ = nSamp;
  adata.iLowerB_ = lbs;
  adata.iUpperB_ = ubs;
  adata.sampleInputs_ = X;
  adata.sampleOutputs_ = Y;
  adata.outputID_ = 0;
  adata.printLevel_ = 0;
  analyze3(adata);
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobol1Analyzer::analyze(aData &adata)
{
  int    nInputs, nOutputs, nSamples, nSteps, outputID, noPDF=1;
  int    ii, jj, kk, iL, iR, sCnt, status, currNLevels, method=1, corFlag;
  int    nSubSamples=1000, nLevels=200, pdfFileFlag=0;
  int    printLevel, nSamp, scatterFileFlag=0, totalCnt, *pdfFlags;
  double *X, *Y2, *xLower, *xUpper, *inputMeans, *inputStdevs, *tempV;
  double *oneSamplePt, variance, ddata, dmean;
  char   winput1[500], winput2[500], pdfFile[500], scatterFile[500];
  char   *cString, pString[500];
  FILE   *fp;
  PsuadeData    *ioPtr;
  FuncApprox    *faPtr;
  RSConstraints *constrPtr;
  PDFManager    *pdfman, *pdfman1, *pdfman2;
  Sampling      *sampler;
  pData         pCorMat;
  psVector  vecIn, vecOut, vecUB, vecLB, vecY, vecXX, vecYY, vecZZ;
  psVector  vecLower2, vecUpper2, vecSamplePts, vecVars, vecEcvs; 
  psVector  vecInpMeans2, vecInpStdvs2, vecSamPts1D, vecMeans, vecVces;
  psVector  vecmSamPts;
  psIVector vecSS, vecBins, vecInpFlags2;
  psMatrix  *corMatp, corMat;

  if (method == 1)
  {
    analyze3(adata);
    return 0;
  }
  nInputs     = adata.nInputs_;
  nOutputs    = adata.nOutputs_;
  nSamples    = adata.nSamples_;
  xLower      = adata.iLowerB_;
  xUpper      = adata.iUpperB_;
  X           = adata.sampleInputs_;
  Y2          = adata.sampleOutputs_;
  outputID    = adata.outputID_;
  ioPtr       = adata.ioPtr_;
  printLevel  = adata.printLevel_;
  pdfFlags    = adata.inputPDFs_;
  inputMeans  = adata.inputMeans_;
  inputStdevs = adata.inputStdevs_;
  corFlag = 0;
  if (pdfFlags != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] != 0) noPDF = 0;
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] == PSUADE_PDF_SAMPLE) corFlag = 1;
  } 
  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO, "* RSMSobol1 constructor\n");
  printDashes(PL_INFO, 0);
  if (noPDF == 1) 
    printOutTS(PL_INFO,"* RSMSobol1 INFO: all uniform distributions.\n");
  else
  {
    printOutTS(PL_INFO,"* RSMSobol1 INFO: non-uniform distributions");
    printOutTS(PL_INFO," detected.\n");
  }

  if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
  {
    printOutTS(PL_ERROR, "RSMSobol1 ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (outputID >= nOutputs || outputID < 0)
  {
    printOutTS(PL_ERROR,"RSMSobol1 ERROR: invalid output ID (%d).\n", 
               outputID);
    return PSUADE_UNDEFINED;
  }
  if (nInputs <= 1)
  {
    printOutTS(PL_ERROR,"RSMSobol1 ERROR: nInputs=1 does not need");
    printOutTS(PL_ERROR," this analysis.\n");
    return PSUADE_UNDEFINED;
  }
  if (ioPtr == NULL)
  {
    printOutTS(PL_ERROR,
         "RSMSobol1 ERROR: no data object (PsuadeData) found.\n");
    return PSUADE_UNDEFINED;
  } 
  ioPtr->getParameter("input_cor_matrix", pCorMat);
  corMatp = (psMatrix *) pCorMat.psObject_;
  for (ii = 0; ii < nInputs; ii++)
  {
    for (jj = 0; jj < ii; jj++)
    {
      if (corMatp->getEntry(ii,jj) != 0.0)
      {
        printOutTS(PL_INFO,"RSMSobol1 INFO: Correlated inputs have\n");
        printOutTS(PL_INFO,"   been detected. PSUADE will use a\n");
        printOutTS(PL_INFO,"   variant of this method.\n");
        corFlag = 1;
      }
    }
  }
  status = 0;
  for (ii = 0; ii < nSamples; ii++)
    if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR, 
       "RSMSobol1 ERROR: Some outputs are undefined. Prune the\n");
    printOutTS(PL_ERROR,
       "                 undefined sample point first.\n");
    return PSUADE_UNDEFINED;
  }
  if (corFlag != 0) return analyze2(adata);

  vecY.setLength(nSamples);
  for (ii = 0; ii < nSamples; ii++) vecY[ii] = Y2[ii*nOutputs+outputID];

  constrPtr = new RSConstraints();
  constrPtr->genConstraints(ioPtr);

  faPtr = genFAInteractive(ioPtr, 0);
  status = faPtr->initialize(X, vecY.getDVector());

  if (psAnaExpertMode_ == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO, "* RSMSobol1 creates one sample of size M \n");
    printOutTS(PL_INFO, "* for each of the K input levels. Therefore,\n");
    printOutTS(PL_INFO, "* the total sample size is\n");
    printOutTS(PL_INFO, "* N = M * K * nInputs\n");
    nSubSamples = 1000;
    nLevels = 200;
    printOutTS(PL_INFO, "* nInputs m = %d, and\n", nInputs);
    printOutTS(PL_INFO, "* default M = %d\n", nSubSamples);
    printOutTS(PL_INFO, "* default K = %d\n", nLevels);
    printOutTS(PL_INFO, "* As a user, please decide on M and K.\n");
    printOutTS(PL_INFO, "* Note: large M and K may take a long time\n");
    printEquals(PL_INFO, 0);
    sprintf(pString,"Enter M (suggestion: 1000-10000) : ");
    nSubSamples = getInt(1000, 50000, pString);
    sprintf(pString,"Enter K (suggestion: 100 - 1000) : ");
    nLevels = getInt(100, 5000, pString);
    printAsterisks(PL_INFO, 0);
  }
  else
  {
    nSubSamples = 1000;
    nLevels = 200;
    if (psConfig_ != NULL)
    {
      cString = psConfig_->getParameter("RSMSobol1_nsubsamples");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
        if (nSubSamples < 1000)
        {
          printOutTS(PL_INFO,
             "RSMSobol1 INFO: nSubSamples should be >= 1000.\n");
          nSubSamples = 1000;
        }
        else
        {
          printOutTS(PL_INFO, 
             "RSMSobol1 INFO: nSubSamples = %d (config).\n",nSubSamples);
        }
      }
      cString = psConfig_->getParameter("RSMSobol1_nlevels");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
        if (nLevels < 200)
        {
          printOutTS(PL_INFO,
             "RSMSobol1 INFO: nLevels should be >= 200.\n");
          nLevels = 200;
        }
        else
        {
          printOutTS(PL_INFO,
             "RSMSobol1 INFO: nLevels = %d (config).\n",nLevels);
        }
      }
    }

    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"* RSMSobol1 creates one sample of size M \n");
    printOutTS(PL_INFO,"* for each of the K input levels. Therefore,\n");
    printOutTS(PL_INFO,"* the total sample size is\n");
    printOutTS(PL_INFO,"* N = M * K * nInputs.\n");
    printOutTS(PL_INFO,"* nInputs m = %d \n", nInputs);
    printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
    printOutTS(PL_INFO,"* default K = %d\n", nLevels);
    printOutTS(PL_INFO,
         "* To change these settings, turn on ana_expert mode and rerun.\n");
    printAsterisks(PL_INFO, 0);
  }

  if (psAnaExpertMode_ == 0 && psConfig_ != NULL)
  {
    cString = psConfig_->getParameter("RSMSobol1_pdffile");
    if (cString != NULL)
    {
      strcpy(pdfFile, cString);
      pdfFileFlag = 1; 
    }
  }
  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO,
         "RSMSobol1 will create a sample for basic statistics. You\n");
    printOutTS(PL_INFO,
         "have the option to plot the probability density function.\n");
    sprintf(pString, "Create a pdf (probability) bar graph? (y or n) ");
    getString(pString, winput1);
    if (winput1[0] == 'y')
    {
      sprintf(pString,"Enter the file name for the bar graph : ");
      getString(pString, pdfFile);
      pdfFile[strlen(pdfFile)-1] = '\0';
      fp = fopen(pdfFile, "w");
      if (fp != NULL)
      {
        fclose(fp);
        pdfFileFlag = 1; 
      }
      else 
      {
        printOutTS(PL_INFO, 
             "RSMSobol1 WARNING: cannot open file %s\n", pdfFile);
        pdfFileFlag = 0; 
      }
    }
    printEquals(PL_INFO, 0);
  }
  else pdfFileFlag = 0;

  if (psAnaExpertMode_ == 0 && psConfig_ != NULL)
  {
    cString = psConfig_->getParameter("RSMSobol1_scatterfile");
    if (cString != NULL)
    {
      strcpy(pdfFile, cString);
      scatterFileFlag = 1; 
    }
  }
  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO,
         "RSMSobol1 will create many samples for Sobol1 analysis.\n");
    printOutTS(PL_INFO,"You have the option to plot these sample data.\n");
    sprintf(pString, "Create a scatter plot for RSMSobol1? (y or n) ");
    getString(pString, winput1);
    if (winput1[0] == 'y')
    {
      sprintf(pString,"Enter the file name for the scatter plot : ");
      getString(pString, scatterFile);
      scatterFile[strlen(scatterFile)-1] = '\0';
      fp = fopen(scatterFile, "w");
      if (fp != NULL)
      {
        fclose(fp);
        scatterFileFlag = 1;
      }
      else 
      {
        printOutTS(PL_INFO, "RSMSobol1 ERROR: cannot open file %s\n",
               scatterFile); 
        scatterFileFlag = 0;
      }
    }
    printEquals(PL_INFO, 0);
  }
  else scatterFileFlag = 0;

  nSamp = 100000;
  if (printLevel > 1)
  {
    printOutTS(PL_INFO,
         "RSMSobol1 INFO: creating a sample for basic statistics.\n");
    printOutTS(PL_INFO,"                sample size = %d\n", nSamp);
  }

  vecXX.setLength(nSamp*nInputs);
  vecYY.setLength(nSamp);

  if (noPDF == 0)
  {
    pdfman = new PDFManager();
    pdfman->initialize(nInputs,pdfFlags,inputMeans,
                       inputStdevs,*corMatp,NULL,NULL);
    vecLB.load(nInputs, xLower);
    vecUB.load(nInputs, xUpper);
    vecOut.setLength(nSamp*nInputs);
    pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
    for (ii = 0; ii < nSamp*nInputs; ii++) vecXX[ii] = vecOut[ii];
    delete pdfman;
  }
  else
  {
    if (nInputs > 51) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
    else              sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    sampler->setInputBounds(nInputs, xLower, xUpper);
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nSamp, 1, 1);
    sampler->initialize(0);
    vecSS.setLength(nSamp);
    sampler->getSamples(nSamp,nInputs,1,vecXX.getDVector(),
                        vecYY.getDVector(),vecSS.getIVector());
    delete sampler;
  }

  printOutTS(PL_INFO,
          "RSMSobol1: running the sample with response surface...\n");
  faPtr->evaluatePoint(nSamp, vecXX.getDVector(), vecYY.getDVector());
  printOutTS(PL_INFO,
          "RSMSobol1: done running the sample with response surface.\n");

  for (ii = 0; ii < nSamp; ii++)
  {
    oneSamplePt = &(vecXX[ii*nInputs]);
    ddata = constrPtr->evaluate(oneSamplePt,vecYY[ii],status);
    if (status == 0) vecYY[ii] = PSUADE_UNDEFINED;
  }

  sCnt  = 0;
  dmean = 0.0;
  for (ii = 0; ii < nSamp; ii++)
  {
    if (vecYY[ii] != PSUADE_UNDEFINED)
    {
      dmean += vecYY[ii];
      sCnt++;
    }
  }
  if (sCnt > 1) dmean /= (double) sCnt;
  else
  {
    printOutTS(PL_ERROR, 
         "RSMSobol1 ERROR: too few samples that satisify the ");
    printOutTS(PL_ERROR, "constraints (%d out of %d).\n", sCnt, nSamp);
    delete faPtr;
    return PSUADE_UNDEFINED;
  }
  variance = 0.0;
  for (ii = 0; ii < nSamp; ii++)
  {
    if (vecYY[ii] != PSUADE_UNDEFINED)
      variance += (vecYY[ii] - dmean) * (vecYY[ii] - dmean) ;
  }
  variance /= (double) sCnt;
  printOutTS(PL_INFO, 
       "RSMSobol1: sample mean    (based on N = %d) = %10.3e\n",
       sCnt, dmean);
  printOutTS(PL_INFO, 
       "RSMSobol1: sample std dev (based on N = %d) = %10.3e\n",
       sCnt, sqrt(variance));
  if (variance == 0.0) variance = 1.0;

  nSteps = 1;
  if (nSamp > 2000) nSteps = nSamp / 1000;
  if (pdfFileFlag == 1)
  {
    fp = fopen(pdfFile, "w");
    if (fp != NULL)
    {
      if (psPlotTool_ == 1)
           fprintf(fp, "// Sample PDF based on response surface.\n");
      else fprintf(fp, "%% Sample PDF based on response surface.\n");
      fprintf(fp, "A = [\n");
      for (ii = 0; ii < nSamp; ii++)
      {
        if (vecYY[ii] != PSUADE_UNDEFINED)
        {
          if ((pdfFileFlag == 1) && (ii % nSteps) == 0) 
            fprintf(fp,"   %e\n", vecYY[ii]);
        }
      }
      fprintf(fp, "];\n");
      fprintf(fp, "amax = max(A);\n");
      fprintf(fp, "amin = min(A);\n");
      fprintf(fp, "step = (amax - amin) / 10;\n");
      fprintf(fp, "numa = size(A,1);\n");
      fprintf(fp, "B = zeros(numa,1);\n");
      fprintf(fp, "for ii = 1 : numa\n");
      fprintf(fp, "   B(ii) = ceil(((A(ii) - amin))/step);\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "X = zeros(10,1);\n");
      fprintf(fp, "Y = zeros(10,1);\n");
      fprintf(fp, "for ii = 1 : 10\n");
      fprintf(fp, "   X(ii) = amin + (ii - 1 + 0.5) * step;\n");
      fprintf(fp, "   [ia,ja,aa] = find(B==ii);\n");
      fprintf(fp, "   Y(ii) = length(ia);\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "Y = Y / numa;\n");
      fprintf(fp, "bar(X,Y,1.0)\n");
      fwritePlotXLabel(fp, "Output Values");
      fwritePlotYLabel(fp, "Probabilities");
      fwritePlotAxes(fp);
      fclose(fp);
    }
  }

  vecLower2.setLength(nInputs);
  vecUpper2.setLength(nInputs);
  nSamp  = nSubSamples;
  vecSamPts1D.setLength(nLevels);
  vecSamplePts.setLength(nInputs*nSubSamples);
  vecXX.setLength(nSubSamples*nInputs);
  vecYY.setLength(nSubSamples);
  vecMeans.setLength(nLevels);
  vecVars.setLength(nLevels);
  vecVces.setLength(nInputs);
  vecEcvs.setLength(nInputs);
  vecBins.setLength(nLevels);
  vecInpFlags2.setLength(nInputs-1);
  vecInpMeans2.setLength(nInputs-1);
  vecInpStdvs2.setLength(nInputs-1);

  fp = NULL;
  if (scatterFileFlag == 1)
  {
    fp = fopen(scatterFile, "w");
    if (psPlotTool_ == 1)
         fprintf(fp, "// scatter plot for RSMSobol1 sample\n");
    else fprintf(fp, "%% scatter plot for RSMSobol1 sample\n");
    fprintf(fp, "A = [ \n");
  }

  for (ii = 0; ii < nInputs; ii++)
  {
    printOutTS(PL_INFO, "RSMSobol1: processing input %d\n", ii+1);

    currNLevels = nLevels / 8;
    for (iR = 0; iR < 4; iR++)
    {
      if (noPDF == 0)
      {
        corMat.setDim(1,1);
        corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
        pdfman1 = new PDFManager();
        pdfman1->initialize(1,&pdfFlags[ii],&inputMeans[ii],
                            &inputStdevs[ii],corMat,NULL,NULL);
        vecLB.load(1, &xLower[ii]);
        vecUB.load(1, &xUpper[ii]);
        vecOut.setLength(currNLevels);
        pdfman1->genSample(currNLevels, vecOut, vecLB, vecUB);
        for (jj = 0; jj < currNLevels; jj++)
           vecSamPts1D[jj] = vecOut[jj];
        delete pdfman1;
      }
      else
      {
        sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
        sampler->setInputBounds(1, &xLower[ii], &xUpper[ii]);
        sampler->setOutputParams(1);
        sampler->setSamplingParams(currNLevels, 1, 1);
        sampler->initialize(0);
        vecSS.setLength(currNLevels);
        vecZZ.setLength(currNLevels);
        sampler->getSamples(currNLevels, 1, 1, vecSamPts1D.getDVector(), 
                            vecZZ.getDVector(), vecSS.getIVector());
        delete sampler;
      }

      if (noPDF == 0)
      {
        corMat.setDim(nInputs-1, nInputs-1);
        for (jj = 0; jj < ii; jj++)
        {
          vecLower2[jj] = xLower[jj];
          vecUpper2[jj] = xUpper[jj];
          vecInpFlags2[jj] = pdfFlags[jj];
          vecInpMeans2[jj] = inputMeans[jj];
          vecInpStdvs2[jj] = inputStdevs[jj];
          for (kk = 0; kk < ii; kk++)
            corMat.setEntry(jj, kk, corMatp->getEntry(jj,kk));
          for (kk = ii+1; kk < nInputs; kk++)
            corMat.setEntry(jj, kk-1, corMatp->getEntry(jj,kk));
        }
        for (jj = ii+1; jj < nInputs; jj++)
        {
          vecLower2[jj-1] = xLower[jj];
          vecUpper2[jj-1] = xUpper[jj];
          vecInpFlags2[jj-1] = pdfFlags[jj];
          vecInpMeans2[jj-1] = inputMeans[jj];
          vecInpStdvs2[jj-1] = inputStdevs[jj];
          for (kk = 0; kk < ii; kk++)
            corMat.setEntry(jj-1, kk, corMatp->getEntry(jj,kk));
          for (kk = ii+1; kk < nInputs; kk++)
            corMat.setEntry(jj-1, kk-1, corMatp->getEntry(jj,kk));
        }
        pdfman2 = new PDFManager();
        pdfman2->initialize(nInputs-1,vecInpFlags2.getIVector(),
                            vecInpMeans2.getDVector(),
                            vecInpStdvs2.getDVector(),corMat,NULL,NULL);
        vecLB.load(nInputs-1, vecLower2.getDVector());
        vecUB.load(nInputs-1, vecUpper2.getDVector());
        vecOut.setLength(nSubSamples*(nInputs-1));
        pdfman2->genSample(nSubSamples, vecOut, vecLB, vecUB);
        for (jj = 0; jj < nSubSamples*(nInputs-1); jj++) 
          vecXX[jj] = vecOut[jj];
        delete pdfman2;
      }
      else
      {
        for (jj = 0; jj < ii; jj++)
        {
          vecLower2[jj] = xLower[jj];
          vecUpper2[jj] = xUpper[jj];
        }
        for (jj = ii+1; jj < nInputs; jj++)
        {
          vecLower2[jj-1] = xLower[jj];
          vecUpper2[jj-1] = xUpper[jj];
        }
        if (nInputs-1 > 51)
             sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
        else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
        sampler->setInputBounds(nInputs-1, vecLower2.getDVector(), 
                                vecUpper2.getDVector());
        sampler->setOutputParams(1);
        sampler->setSamplingParams(nSubSamples, 1, 1);
        sampler->initialize(0);
        vecSS.setLength(nSubSamples);
        vecZZ.setLength(nSubSamples);
        sampler->getSamples(nSubSamples,nInputs-1,1,vecXX.getDVector(),
                            vecZZ.getDVector(),vecSS.getIVector());
        delete sampler;
      }

      for (iL = 0; iL < currNLevels; iL++)
      {
        for (jj = 0; jj < nSubSamples; jj++)
        {
          oneSamplePt = vecXX.getDVector();
          oneSamplePt = &(oneSamplePt[jj*(nInputs-1)]);
          for (kk = 0; kk < ii; kk++)
            vecSamplePts[jj*nInputs+kk] = oneSamplePt[kk];
          for (kk = ii+1; kk < nInputs; kk++)
            vecSamplePts[jj*nInputs+kk] = oneSamplePt[kk-1];
          vecSamplePts[jj*nInputs+ii] = vecSamPts1D[iL];
        }

        faPtr->evaluatePoint(nSubSamples,vecSamplePts.getDVector(), 
                             vecYY.getDVector());

        for (jj = 0; jj < nSubSamples; jj++)
        {
          oneSamplePt = &(vecSamplePts[jj*nInputs]);
          ddata = constrPtr->evaluate(oneSamplePt, vecYY[jj], status);
          if (status == 0) vecYY[jj] = PSUADE_UNDEFINED;

          if (scatterFileFlag == 1 && iR == 0 && ddata != PSUADE_UNDEFINED)
          {
            for (kk = 0; kk < nInputs; kk++)
              fprintf(fp, "%e ", vecSamplePts[jj*nInputs+kk]);
            fprintf(fp, " %e\n", vecYY[jj]);
          }
        }

        vecMeans[iL] = 0.0;
        sCnt = 0;
        for (jj = 0; jj < nSubSamples; jj++)
        {
          if (vecYY[jj] != PSUADE_UNDEFINED)
          {
            vecMeans[iL] += vecYY[jj];
            sCnt++;
          }
        }
        vecBins[iL] = sCnt;
        if (sCnt < nSubSamples/10 && printLevel >= 5)
          printOutTS(PL_WARN, 
               "RSMSobol1 WARNING: subsample size = %d\n",sCnt);
        if (sCnt >= 1) vecMeans[iL] /= (double) sCnt;
        else           vecMeans[iL] = PSUADE_UNDEFINED;
        if (printLevel > 3)
        {
          printOutTS(PL_INFO,"RSMSobol1: input %d :\n", ii+1);
          printOutTS(PL_INFO,
               "  refinement %2d, level %3d, size %d), mean = %e\n",iR,
               iL+1, nSubSamples, vecMeans[iL]);
        }
        vecVars[iL] = 0.0;
        if (sCnt > 1)
        {
          for (jj = 0; jj < nSubSamples; jj++)
             if (vecYY[jj] != PSUADE_UNDEFINED)
               vecVars[iL] += (vecYY[jj]-vecMeans[iL])*
                              (vecYY[jj]-vecMeans[iL]);
          vecVars[iL] /= (double) sCnt;
        }
        else vecVars[iL] = PSUADE_UNDEFINED;
      }

      totalCnt = 0;
      for (iL = 0; iL < currNLevels; iL++) totalCnt += vecBins[iL];
      if (totalCnt == 0) 
      {
        printOutTS(PL_ERROR, "RSMSobol1 ERROR: no feasible region.\n");
        exit(1);
      }

      ddata = 0.0;
      for (iL = 0; iL < currNLevels; iL++) 
      {
        if (vecMeans[iL] != PSUADE_UNDEFINED)
          ddata += vecMeans[iL] * vecBins[iL] / totalCnt;
      }
      vecVces[ii] = 0.0;
      for (iL = 0; iL < currNLevels; iL++)
        if (vecMeans[iL] != PSUADE_UNDEFINED)
          vecVces[ii] += (vecMeans[iL] - ddata) * (vecMeans[iL] - ddata) * 
                      vecBins[iL] / totalCnt;
      if (printLevel > 2 || iR == 3)
      {
        printOutTS(PL_INFO,
             "VCE(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
             ii+1, iR, vecVces[ii], vecVces[ii]/variance);
      }
      if (iR == 3) vecVces[ii] = vecVces[ii] / variance;

      vecEcvs[ii] = 0.0;
      for (iL = 0; iL < currNLevels; iL++)
      {
        if (vecVars[iL] != PSUADE_UNDEFINED)
          vecEcvs[ii] += vecVars[iL] * vecBins[iL] / totalCnt;
      }
      if (printLevel > 2)
      {
        printOutTS(PL_INFO,
             "ECV(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
             ii+1, iR, vecEcvs[ii], vecEcvs[ii]/variance);
      }
      currNLevels *= 2;
    }
  }
  if (printLevel > 0)
  {
    printAsterisks(PL_INFO, 0);
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, "RSMSobol1: Normalized VCE for input %3d = %e\n",
                 ii+1,vecVces[ii]);
    for (ii = 0; ii < nInputs; ii++) vecMeans[ii] = (double) ii;
    printEquals(PL_INFO, 0);
    sortDbleList2(nInputs, vecVces.getDVector(), vecMeans.getDVector());
    for (ii = nInputs-1; ii >= 0; ii--)
      printOutTS(PL_INFO,
           "RSMSobol1: Normalized VCE (ordered) for input %3d = %e\n",
           (int) vecMeans[ii]+1,vecVces[ii]);
    printAsterisks(PL_INFO, 0);
  }

  if (scatterFileFlag == 1)
  {
    fprintf(fp, "];\n");
    fprintf(fp, "for ii = 1 : %d\n", nInputs);
    fprintf(fp, "   plot(A(:,ii),A(:,%d),'x')\n", nInputs+1);
    fwritePlotAxes(fp);
    fwritePlotXLabel(fp, "['input ' int2str(ii)]");
    fwritePlotYLabel(fp, "Y");
    fprintf(fp, "   disp('Press enter to continue')\n");
    fprintf(fp, "   pause\n");
    fprintf(fp, "end\n");
    fclose(fp);
  }
    
  pData *pPtr;
  if (ioPtr != NULL)
  {
    pPtr = ioPtr->getAuxData();
    pPtr->nDbles_ = nInputs;
    pPtr->dbleArray_ = new double[nInputs];
    for (ii = 0; ii < nInputs; ii++)
      pPtr->dbleArray_[ii] = vecVces[ii] * variance;
    pPtr->dbleData_ = variance;
  }

  delete faPtr;
  delete constrPtr;
  return 0.0;
}

// ************************************************************************
// perform analysis (for problems with joint PDFs)
// ------------------------------------------------------------------------
double RSMSobol1Analyzer::analyze2(aData &adata)
{
  int    ii, jj, iL, iR, status, currNLevels, nSteps, bin, totalCnt;
  double ddata, width;
  char   pdfFile[500], scatterFile[500], pString[500], winput1[500];
  FILE   *fp;
  FuncApprox    *faPtr;
  RSConstraints *constrPtr;
  psVector      vecIn, vecOut, vecUB, vecLB, vecXX, vecYY, vecY;
  PDFManager    *pdfman;

  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO,"RSMSobol1: since joint PDFs have been specified,\n");
  printOutTS(PL_INFO,"           a different interaction analysis will\n");
  printOutTS(PL_INFO,"           be performed.\n");
  printEquals(PL_INFO, 0);
  int nInputs    = adata.nInputs_;
  int nOutputs   = adata.nOutputs_;
  int nSamples   = adata.nSamples_;
  int outputID   = adata.outputID_;
  int printLevel = adata.printLevel_;
  double *xLower = adata.iLowerB_;
  double *xUpper = adata.iUpperB_;
  double *XIn    = adata.sampleInputs_;
  double *YIn    = adata.sampleOutputs_;
  PsuadeData *ioPtr = adata.ioPtr_;
  vecY.setLength(nSamples);
  for (ii = 0; ii < nSamples; ii++) vecY[ii] = YIn[ii*nOutputs+outputID];

  constrPtr = new RSConstraints();
  constrPtr->genConstraints(ioPtr);

  faPtr = genFAInteractive(ioPtr, 0);
  if(faPtr == NULL)
  {
    printOutTS(PL_INFO, "faPtr is NULL in file %s line %d aborting. \n", 
               __FILE__, __LINE__);
    abort();
  }
  else status = faPtr->initialize(XIn, vecY.getDVector());

  int nSubSamples=2000, nLevels=500;
  printAsterisks(PL_INFO, 0);
  if (psAnaExpertMode_ == 1)
  {
    printAsterisks(PL_INFO, 0);
    printf("* RSMSobol1 creates one sample of size M \n");
    printf("* for each of the K input levels. Therefore,\n");
    printf("* the total sample size is\n");
    printf("* N = M * K.\n");
    printf("* NOW, nInputs = %d\n", nInputs);
    printf("* As a user, please decide on M and K.\n\n");
    printEquals(PL_INFO, 0);

    sprintf(pString,"Enter M (suggestion: 100-10000) : ");
    nSubSamples = getInt(100, 100000, pString);
    if (nSubSamples > 100000)
      printOutTS(PL_INFO, "An M of %d may take very long time.\n", 
                 nSubSamples);

    sprintf(pString,"Enter K (suggestion: 1000 - 10000) : ");
    nLevels = getInt(1000, 10000, pString);
    if (nLevels > 5000)
      printOutTS(PL_INFO, "* A K of %d may take very long time.\n", 
                 nLevels);
    printAsterisks(PL_INFO, 0);
  }
  else
  {
    nSubSamples = 100;
    nLevels     = 1000;
    printOutTS(PL_INFO,"* RSMSobol1: default M = %d.\n", nSubSamples);
    printOutTS(PL_INFO,"* RSMSobol1: default K = %d.\n", nLevels);
    printOutTS(PL_INFO,"* To change settings, rerun with ana_expert on.\n");
  }
  printEquals(PL_INFO, 0);

  int pdfFileFlag=0;
  if (psMasterMode_ == 1)
  {
    printOutTS(PL_INFO,"* RSMSobol1 will create a sample for basic\n");
    printOutTS(PL_INFO,"* statistics. You have the option to plot\n");
    printOutTS(PL_INFO,"* the PDF.\n");
    sprintf(pString, "Create a pdf (probability) bar graph? (y or n) ");
    getString(pString, winput1);
    if (winput1[0] == 'y')
    {
      sprintf(pString,"Enter the file name for the bar graph : ");
      getString(pString, pdfFile);
      pdfFile[strlen(pdfFile)-1] = '\0';
      fp = fopen(pdfFile, "w");
      if (fp != NULL)
      {
        fclose(fp);
        if (psPlotTool_ == 0) pdfFileFlag = 1; 
      }
      else 
      {
        printOutTS(PL_INFO,"RSMSobol1 WARNING: cannot open file %s\n", 
                   pdfFile);
        pdfFileFlag = 0; 
      }
    }
    printEquals(PL_INFO, 0);
  }
  else pdfFileFlag = 0;

  int scatterFileFlag=0;;
  if (psMasterMode_ == 1)
  {
    printOutTS(PL_INFO,"RSMSobol1 will create many samples for Sobol1\n");
    printOutTS(PL_INFO,"analysis. You have the option to plot these\n");
    printOutTS(PL_INFO,"sample data.\n");
    sprintf(pString, "Create a scatter plot for RSMSobol1? (y or n) ");
    getString(pString, winput1);
    if (winput1[0] == 'y')
    {
      sprintf(pString,"Enter the file name for the scatter plot : ");
      getString(pString, scatterFile);
      scatterFile[strlen(scatterFile)-1] = '\0';
      fp = fopen(scatterFile, "w");
      if (fp != NULL)
      {
        fclose(fp);
        scatterFileFlag = 1;
      }
      else 
      {
        printOutTS(PL_ERROR, "RSMSobol1 ERROR: cannot open file %s\n",
               scatterFile); 
        scatterFileFlag = 0;
      }
    }
    printEquals(PL_INFO, 0);
  }
  else scatterFileFlag = 0;

  int nSamp = nSubSamples * nLevels;
  if (nSamp < 100000) nSamp = 100000;
  if (printLevel > 1)
  {
    printOutTS(PL_INFO,"* RSMSobol1 INFO: creating a sample for basic\n");
    printOutTS(PL_INFO,"*           statistics. Sample size = %d\n",nSamp);
  }
  vecXX.setLength(nSamp*nInputs);
  vecYY.setLength(nSamp);
  pdfman = new PDFManager();
  pdfman->initialize(ioPtr);
  vecLB.load(nInputs, xLower);
  vecUB.load(nInputs, xUpper);
  vecOut.setLength(nSamp*nInputs);
  pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
  for (ii = 0; ii < nSamp*nInputs; ii++) vecXX[ii] = vecOut[ii];
  printOutTS(PL_INFO,"RSMSobol1: response surface evaluation begins...\n");
  faPtr->evaluatePoint(nSamp, vecXX.getDVector(), vecYY.getDVector());
  printOutTS(PL_INFO,"RSMSobol1: response surface evaluation ends...\n");
  double *oneSamplePt;
  for (ii = 0; ii < nSamp; ii++)
  {
    oneSamplePt = &(vecXX[ii*nInputs]);
    status = 1;
    if (constrPtr != NULL)
      ddata = constrPtr->evaluate(oneSamplePt,vecYY[ii],status);
    if (status == 0) vecYY[ii] = PSUADE_UNDEFINED;
  }
  double dmean=0;
  int sCnt = 0;
  for (ii = 0; ii < nSamp; ii++)
  {
    if (vecYY[ii] != PSUADE_UNDEFINED)
    {
      dmean += vecYY[ii];
      sCnt++;
    }
  }
  if (sCnt > 1) dmean /= (double) sCnt;
  else
  {
    printOutTS(PL_ERROR,"RSMSobol1 ERROR: too few samples that satisify");
    printOutTS(PL_ERROR,"the constraints (%d out of %d).\n", sCnt, nSamp);
    delete faPtr;
    delete pdfman;
    if (constrPtr != NULL) delete constrPtr;
    return PSUADE_UNDEFINED;
  }
  double variance = 0.0;
  for (ii = 0; ii < nSamp; ii++)
  {
    if (vecYY[ii] != PSUADE_UNDEFINED)
      variance += (vecYY[ii] - dmean) * (vecYY[ii] - dmean) ;
  }
  variance /= (double) sCnt;
  printOutTS(PL_INFO,
       "* RSMSobol1: sample mean    (based on %d points) = %e\n",
       sCnt,dmean);
  printOutTS(PL_INFO,
       "* RSMSobol1: sample std dev (based on %d points) = %e\n",
       sCnt, sqrt(variance));
  if (variance == 0.0) variance = 1.0;
  delete pdfman;
  printAsterisks(PL_INFO, 0);
  if (pdfFileFlag == 1 && (fp = fopen(pdfFile, "w")) != NULL)
  {
    if (psPlotTool_ == 1)
         fprintf(fp,"// sample PDF based on the response surface.\n");
    else fprintf(fp,"%% sample PDF based on the response surface.\n");
    fprintf(fp, "A = [\n");
    nSteps = 1;
    if (nSamp > 2000) nSteps = nSamp / 1000;
    for (ii = 0; ii < nSamp; ii++)
    {
      if (vecYY[ii] != PSUADE_UNDEFINED)
      {
        if ((pdfFileFlag == 1) && (ii % nSteps) == 0) 
          fprintf(fp,"   %e\n", vecYY[ii]);
      }
    }
    fprintf(fp, "];\n");
    fprintf(fp, "amax = max(A);\n");
    fprintf(fp, "amin = min(A);\n");
    fprintf(fp, "step = (amax - amin) / 10;\n");
    fprintf(fp, "numa = size(A,1);\n");
    fprintf(fp, "B = zeros(numa,1);\n");
    fprintf(fp, "for ii = 1 : numa\n");
    fprintf(fp, "   B(ii) = ceil(((A(ii) - amin))/step);\n");
    fprintf(fp, "end;\n");
    fprintf(fp, "X = zeros(10,1);\n");
    fprintf(fp, "Y = zeros(10,1);\n");
    fprintf(fp, "for ii = 1 : 10\n");
    fprintf(fp, "   X(ii) = amin + (ii - 1 + 0.5) * step;\n");
    fprintf(fp, "   [ia,ja,aa] = find(B==ii);\n");
    fprintf(fp, "   Y(ii) = length(ia);\n");
    fprintf(fp, "end;\n");
    fprintf(fp, "Y = Y / numa;\n");
    fprintf(fp, "bar(X,Y,1.0)\n");
    fwritePlotXLabel(fp, "Output Values");
    fwritePlotYLabel(fp, "Probabilities");
    fwritePlotAxes(fp);
    fclose(fp);
  }
  fp = NULL;
  if (scatterFileFlag == 1 && (fp = fopen(scatterFile, "w")) != NULL)
  {
    //fp = fopen(scatterFile, "w");
    fprintf(fp, "%% scatter plot for RSMSobol1 sample\n");
    fprintf(fp, "A = [ \n");
    nSteps = 1;
    if (nSamp > 2000) nSteps = nSamp / 1000;
    for (ii = 0; ii < nSamp; ii++)
    {
      if (vecYY[ii] != PSUADE_UNDEFINED)
      {
        if ((pdfFileFlag == 1) && (ii % nSteps) == 0) 
          fprintf(fp,"   %e\n", vecYY[ii]);
      }
    }
    fprintf(fp, "];\n");
    fprintf(fp, "for ii = 1 : %d\n", nInputs);
    fprintf(fp, "   plot(A(:,ii),A(:,%d),'x')\n", nInputs+1);
    fwritePlotXLabel(fp, "['input ' int2str(ii)]");
    fwritePlotYLabel(fp, "Y");
    fwritePlotAxes(fp);
    fprintf(fp, "   disp('Press enter to continue')\n");
    fprintf(fp, "   pause\n");
    fprintf(fp, "end\n");
    fclose(fp);
  }

  psVector  vecVces, vecVars, vecEcvs, vecMeans;
  psIVector vecBins;
  vecMeans.setLength(nLevels);
  vecVars.setLength(nLevels);
  vecVces.setLength(nInputs);
  vecEcvs.setLength(nInputs);
  vecBins.setLength(nLevels);

  for (ii = 0; ii < nInputs; ii++)
  {
    printOutTS(PL_INFO, "RSMSobol1: processing input %d\n", ii+1);

    currNLevels = nLevels / 8;
    for (iR = 0; iR < 4; iR++)
    {
      width = (xUpper[ii] - xLower[ii]) / currNLevels;
      for (iL = 0; iL < currNLevels; iL++)
      {
        vecBins[iL] = 0;
        vecMeans[iL] = 0.0;
        vecVars[iL] = 0.0;
      }
      for (jj = 0; jj < nSamp; jj++)
      {
        if (vecYY[jj] != PSUADE_UNDEFINED)
        {
          ddata = vecXX[jj*nInputs+ii];
          bin = (int) ((ddata - xLower[ii]) / width); 
          if (bin == currNLevels) bin = currNLevels - 1;
          vecMeans[bin] += vecYY[jj];
          vecBins[bin]++;
        }
      }
      for (iL = 0; iL < currNLevels; iL++)
      {
        if (vecBins[iL] < nSubSamples/10 && printLevel >= 5)
          printOutTS(PL_INFO,
               "RSMSobol1 WARNING: subsample size = %d\n",vecBins[iL]);
        if (vecBins[iL] >= 1) vecMeans[iL] /= (double) vecBins[iL];
        else                  vecMeans[iL] = PSUADE_UNDEFINED;
        if (printLevel > 3)
        {
          printOutTS(PL_INFO,"RSMSobol1: input %d :\n", ii+1);
          printOutTS(PL_INFO,
               "  refinement %2d, level %3d, size %d), mean = %e\n",
               iR, iL, nSubSamples, vecMeans[iL]);
        }
      }
      for (jj = 0; jj < nSamp; jj++)
      {
        ddata = vecXX[jj*nInputs+ii];
        bin = (int) ((ddata - xLower[ii]) / width); 
        if (bin == currNLevels) bin = currNLevels - 1;
        if (vecYY[jj] != PSUADE_UNDEFINED)
          vecVars[bin] += pow(vecYY[jj]-vecMeans[bin], 2.0); 
      }
      for (iL = 0; iL < currNLevels; iL++)
      {
        if (vecBins[iL] >= 1) vecVars[iL] /= (double) vecBins[iL];
        else                  vecVars[iL] = PSUADE_UNDEFINED;
      }

      totalCnt = 0;
      for (iL = 0; iL < currNLevels; iL++) totalCnt += vecBins[iL];
      if (totalCnt == 0) 
      {
        printOutTS(PL_ERROR, "RSMSobol1 ERROR: no feasible region.\n");
        exit(1);
      }

      ddata = 0.0;
      for (iL = 0; iL < currNLevels; iL++) 
      {
        if (vecMeans[iL] != PSUADE_UNDEFINED)
          ddata += vecMeans[iL] * vecBins[iL] / totalCnt;
      }
      vecVces[ii] = 0.0;
      for (iL = 0; iL < currNLevels; iL++)
        if (vecMeans[iL] != PSUADE_UNDEFINED)
          vecVces[ii] += (vecMeans[iL]-ddata) * (vecMeans[iL]-ddata) * 
                        vecBins[iL] / totalCnt;
      if (printLevel > 2 || iR == 3)
      {
        printOutTS(PL_INFO,
             "VCE(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
             ii+1, iR, vecVces[ii], vecVces[ii]/variance);
      }
      if (iR == 3) vecVces[ii] = vecVces[ii] / variance;

      vecEcvs[ii] = 0.0;
      for (iL = 0; iL < currNLevels; iL++)
      {
        if (vecVars[iL] != PSUADE_UNDEFINED)
          vecEcvs[ii] += vecVars[iL] * vecBins[iL] / totalCnt;
      }

      printOutTS(PL_INFO, 
           "ECV(%3d) (refinement=%3d) = %e, (normalized) = %e\n",
           ii+1, iR, vecEcvs[ii], vecEcvs[ii]/variance);

      currNLevels *= 2;
    }
  }
  for (ii = 0; ii < nInputs; ii++)
    printOutTS(PL_INFO,"RSMSobol1: Normalized VCE for input %3d = %e\n",
               ii+1,vecVces[ii]);

  delete faPtr;
  if (constrPtr != NULL) delete constrPtr;
  return 0.0;
}

// ************************************************************************
// perform analysis using fuzzy evaluation (12/2012)
// ------------------------------------------------------------------------
double RSMSobol1Analyzer::analyze3(aData &adata)
{
  int    ii, jj, kk, nn, ss, count, status, nSamp, rstype;
  int    nSubSamples=1000, nLevels=200, ntimes=1, nConstr, pdfNull=0;
  double dmean, dstds, ddata, ddata2, frac=0.8, *tempV, *oneSamplePt;
  char   pString[500], *cString, winput1[500], winput2[500];
  FuncApprox    *faPtr;
  RSConstraints *constrPtr;
  psVector      vecIn, vecOut, vecUB, vecLB, vecY;
  pData         pCorMat, pPtr;
  psMatrix      *corMatp=NULL, corMat;
  PDFManager    *pdfman;
  Sampling      *sampler;

  if (isScreenDumpModeOn())
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"*          RS-based First Order Sobol' Analysis \n");
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET\n");
    printOutTS(PL_INFO,"* - ana_expert to finetune RSMSobol1 parameters\n");
    printOutTS(PL_INFO,"*   (e.g. to adjust integration sample size).\n");
    printOutTS(PL_INFO,"* - rs_expert mode to finetune response surface\n");
    printOutTS(PL_INFO,"* - printlevel to display more information\n");
    printOutTS(PL_INFO,"* Or, use configure file to finetune parameters\n");
    printEquals(PL_INFO, 0);
  }
 
  int    nInputs, nOutputs, nSamples, outputID, printLevel;
  int    noPDF=1, corFlag=0, *pdfFlags;
  double *xLower, *xUpper, *X, *Y2, *inputMeans, *inputStdevs;
  PsuadeData *ioPtr;

  nInputs     = adata.nInputs_;
  nInputs_    = nInputs;
  nOutputs    = adata.nOutputs_;
  nSamples    = adata.nSamples_;
  xLower      = adata.iLowerB_;
  xUpper      = adata.iUpperB_;
  X           = adata.sampleInputs_;
  Y2          = adata.sampleOutputs_;
  outputID    = adata.outputID_;
  ioPtr       = adata.ioPtr_;
  printLevel  = adata.printLevel_;
  pdfFlags    = adata.inputPDFs_;
  inputMeans  = adata.inputMeans_;
  inputStdevs = adata.inputStdevs_;
  if (inputMeans == NULL || pdfFlags == NULL || inputStdevs == NULL)
  {
    pdfNull = 1;
    pdfFlags    = new int[nInputs];
    inputMeans  = new double[nInputs];
    inputStdevs = new double[nInputs];
    checkAllocate(inputStdevs, "inputStdevs in RSMSobol1::analyze3");
    for (ii = 0; ii < nInputs; ii++)
    {
      pdfFlags[ii] = 0;
      inputMeans[ii]  = 0;
      inputStdevs[ii] = 0;
    }
  }
  if (pdfFlags != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] != 0) noPDF = 0;
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] == PSUADE_PDF_SAMPLE) corFlag = 1;
  } 
  if (isScreenDumpModeOn())
  {
    if (noPDF == 1) 
      printOutTS(PL_INFO,"* RSMSobol1 INFO: all uniform distributions.\n");
    else
    {
      printOutTS(PL_INFO,"RSMSobol1 INFO: non-uniform distributions");
      printOutTS(PL_INFO," detected.\n");
    }
  }

  if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
  {
    printOutTS(PL_ERROR, "RSMSobol1 ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (outputID >= nOutputs || outputID < 0)
  {
    printOutTS(PL_ERROR,"RSMSobol1 ERROR: invalid output ID (%d).\n", 
               outputID);
    return PSUADE_UNDEFINED;
  }
  if (nInputs <= 1)
  {
    printOutTS(PL_ERROR,
         "RSMSobol1 INFO: analysis not needed for nInputs=1\n");
    return PSUADE_UNDEFINED;
  }

  if (ioPtr == NULL)
  {
    if (isScreenDumpModeOn())
    {
      printOutTS(PL_INFO,
        "RSMSobol1 INFO: no data object (PsuadeData) found.\n");
      printOutTS(PL_INFO,"   Several features will be turned off.\n");
    }
    corMatp = new psMatrix();
    corMatp->setDim(nInputs, nInputs);
    for (ii = 0; ii < nInputs; ii++) corMatp->setEntry(ii,ii,1.0e0);
  } 
  else
  {
    ioPtr->getParameter("input_cor_matrix", pCorMat);
    corMatp = (psMatrix *) pCorMat.psObject_;
    for (ii = 0; ii < nInputs; ii++)
    {
      for (jj = 0; jj < ii; jj++)
      {
        if (corMatp->getEntry(ii,jj) != 0.0)
        {
          if (isScreenDumpModeOn())
          {
            printOutTS(PL_INFO, 
              "RSMSobol1 INFO: Correlated inputs detected.\n");
            printOutTS(PL_INFO, 
              "   Alternative analyis is to be performed.\n");
          }
          corFlag = 1;
        }
      }
    }
  }

  status = 0;
  for (ii = 0; ii < nSamples; ii++)
    if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR,"RSMSobol1 ERROR: Some outputs are undefined.\n");
    printOutTS(PL_ERROR,"     Prune the undefined sample point first.\n");
    return PSUADE_UNDEFINED;
  }

  if (corFlag != 0) 
  {
    if (ioPtr == NULL) delete corMatp;
    if (pdfNull == 1)
    {
       delete [] pdfFlags;
       delete [] inputMeans;
       delete [] inputStdevs;
    }
    return analyze2(adata);
  }

  if (ioPtr != NULL)
  {
    constrPtr = new RSConstraints();
    constrPtr->genConstraints(ioPtr);
    nConstr = constrPtr->getNumConstraints();
  }
  else
  {
    nConstr = 0;
    constrPtr = NULL;
    if (isScreenDumpModeOn())
      printf("RSMSobolTSI INFO: no PsuadeData ==> no constraints.\n");
  }

  if (ioPtr == NULL)
  {
    if (rstype_ < 0)
    {
      printf("Select response surface. Options are: \n");
      writeFAInfo(0);
      strcpy(pString, "Choose response surface: ");
      rstype = getInt(0, PSUADE_NUM_RS, pString);
    }
    else rstype = rstype_;
  }
  else
  {
    ioPtr->getParameter("ana_rstype", pPtr);
    rstype = pPtr.intData_;
  }

  if (isScreenDumpModeOn() && psAnaExpertMode_ == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"* RSMSobol1 generates K levels for each input\n");
    printOutTS(PL_INFO,"* and creates a sample of size M for each level.\n");
    printOutTS(PL_INFO,"* Therefore, the total sample size for this\n");
    printOutTS(PL_INFO,"* analysis is:\n");
    printOutTS(PL_INFO,"*      N = M * K * nInputs\n");
    nSubSamples = 100;
    nLevels = 1000;
    printOutTS(PL_INFO,"* nInputs   = %d\n", nInputs);
    printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
    printOutTS(PL_INFO,"* default K = %d\n", nLevels);
    printOutTS(PL_INFO,"* Recommendation: K >> M\n");
    printOutTS(PL_INFO,"* Please select different M and K.\n");
    printOutTS(PL_INFO,"* NOTE: large M and K may take very long time.\n");
    printEquals(PL_INFO, 0);
    sprintf(pString,"Enter M (suggestion: 100 - 10000) : ");
    nSubSamples = getInt(100, 50000, pString);
    sprintf(pString,"Enter K (suggestion: 1000 - 10000) : ");
    nLevels = getInt(1000, 50000, pString);
    printOutTS(PL_INFO,"* To include response surface uncertainties in\n");
    printOutTS(PL_INFO,"* this analysis, the sensitivity calculation\n");
    printOutTS(PL_INFO,"* is to be repeated a number of times using\n");
    printOutTS(PL_INFO,"* different bootstrapped samples. Please specify\n");
    printOutTS(PL_INFO,"* the number of bootstrapped samples below.\n");
    printOutTS(PL_INFO,"* If you do not need error bars, set it to 1.\n");
    sprintf(pString,"Enter the number of bootstrapped samples (1 - 500) : ");
    ntimes = getInt(1, 500, pString);
    printAsterisks(PL_INFO, 0);
  }
  else
  {
    nSubSamples = 100;
    nLevels = 1000;
    if (psConfig_ != NULL)
    {
      cString = psConfig_->getParameter("RSMSobol1_nsubsamples");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
        if (nSubSamples < 100)
        {
          printOutTS(PL_INFO,
               "RSMSobol1 INFO: nSubSamples should be >= 100.\n");
          nSubSamples = 100;
        }
        else
        {
          printOutTS(PL_INFO,
               "RSMSobol1 INFO: nSubSamples = %d (config).\n",
               nSubSamples);
        }
      }
      cString = psConfig_->getParameter("RSMSobol1_nlevels");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
        if (nLevels < 1000)
        {
          printOutTS(PL_INFO,
               "RSMSobol1 INFO: nLevels should be >= 1000.\n");
          nLevels = 1000;
        }
        else
        {
          printOutTS(PL_INFO,
               "RSMSobol1 INFO: nLevels = %d (config).\n",nLevels);
        }
      }
    }
    if (isScreenDumpModeOn() && printLevel > 1)
    {
      printOutTS(PL_INFO,"* RSMSobol1 creates one sample of size M \n");
      printOutTS(PL_INFO,"* for each of the K levels of each input.\n");
      printOutTS(PL_INFO,"* Therefore, the total sample size is\n");
      printOutTS(PL_INFO,"* N = M * K * nInputs\n");
      printOutTS(PL_INFO,"* nInputs m = %d\n", nInputs);
      printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
      printOutTS(PL_INFO,"* default K = %d\n", nLevels);
      printOutTS(PL_INFO,
           "* To make changes, re-run with ana_expert mode on.\n");
      printEquals(PL_INFO, 0);
    }
  }
  if (ntimes > 1)
  {
    if (isScreenDumpModeOn())
    {
      printOutTS(PL_INFO,"RSMSobol1 INFO: number of bootstrapped samples\n");
      printOutTS(PL_INFO,"          greater than 1 is not recommended.\n");
      printOutTS(PL_INFO,"          Use the rssobol1b command instead.\n");
    }
    ntimes = 1;
  }
  if (ntimes == 1) frac = 1.0;

  nSamp = nLevels * nSubSamples;
  if (isScreenDumpModeOn())
  {
    printOutTS(PL_INFO,
      "RSMSobol1 INFO: creating a sample for basic statistics.\n");
    printOutTS(PL_INFO,"                sample size = %d\n", nSamp);
  }

  psIVector vecSS;
  psVector  vecXX, vecYY;
  vecXX.setLength(nSamp*nInputs);
  vecYY.setLength(nSamp);
  if (noPDF == 0)
  {
    pdfman = new PDFManager();
    pdfman->initialize(nInputs,pdfFlags,inputMeans,inputStdevs,*corMatp,
                       NULL,NULL);
    vecLB.load(nInputs, xLower);
    vecUB.load(nInputs, xUpper);
    vecOut.setLength(nSamp*nInputs);
    pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
    for (ii = 0; ii < nSamp*nInputs; ii++) vecXX[ii] = vecOut[ii];
    delete pdfman;
  }
  else
  {
    if (nInputs > 51)
       sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
    else
       sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    sampler->setInputBounds(nInputs, xLower, xUpper);
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nSamp, 1, 1);
    sampler->initialize(0);
    vecSS.setLength(nSamp);
    sampler->getSamples(nSamp, nInputs, 1, vecXX.getDVector(), 
                        vecYY.getDVector(), vecSS.getIVector());
    delete sampler;
  }

  psVector  vecBsX, vecBsY, vecBsMeans, vecBsStds;
  psIVector vecBsFlags;

  vecBsFlags.setLength(nSamples);
  vecBsX.setLength(nSamples*nInputs);
  vecBsY.setLength(nSamples);
  vecBsMeans.setLength(ntimes);
  vecBsStds.setLength(ntimes);
  count = (int) (frac * nSamples);
  faPtr = genFA(rstype, nInputs, 0, count);
  faPtr->setBounds(xLower, xUpper);
  faPtr->setOutputLevel(0);
  for (nn = 0; nn < ntimes; nn++)
  {
    if (ntimes == 1)
    {
      for (ss = 0; ss < nSamples*nInputs; ss++) vecBsX[ss] = X[ss];
      for (ss = 0; ss < nSamples; ss++) 
        vecBsY[ss] = Y2[ss*nOutputs+outputID];
    }
    else
    {
      for (ss = 0; ss < nSamples; ss++) vecBsFlags[ss] = 0;
      count = 0;
      while (count < frac * nSamples)
      {
        jj = PSUADE_rand() % nSamples;
        if (vecBsFlags[jj] == 0)
        {
          for (ii = 0; ii < nInputs; ii++)
            vecBsX[count*nInputs+ii] = X[jj*nInputs+ii];
          vecBsY[count] = Y2[jj*nOutputs+outputID];
          vecBsFlags[jj] = 1;
          count++;
        }
      }
    }
    status = faPtr->initialize(vecBsX.getDVector(), vecBsY.getDVector());
    if (status != 0)
    {
      printf("RSMSobol1 ERROR: in initializing response surface.\n");
      if (pdfNull == 1)
      {
        delete [] pdfFlags;
        delete [] inputMeans;
        delete [] inputStdevs;
      }
      return -1;
    }
    faPtr->evaluatePoint(nSamp,vecXX.getDVector(),vecYY.getDVector());
    count = 0;
    if (nConstr > 0)
    {
      tempV = vecXX.getDVector();
      for (kk = 0; kk < nSamp; kk++)
      {
        ddata = constrPtr->evaluate(&tempV[kk*nInputs],vecYY[kk],
                                    status);
        if (status == 0) 
        {
           vecYY[kk] = PSUADE_UNDEFINED;
           count++;
        }
      }
    }
    count = nSamp - count;
    if (nConstr > 0)
    {
      printOutTS(PL_INFO,
           "RSMSobol1 INFO: %6.2f percent passes the contraints.\n",
           (double) count * 100.0 /((double) nSamp));
    }
    if (count <= 1)
    {
      printf("RSMSobol1 ERROR: too few samples left after filtering\n");
      if (pdfNull == 1)
      {
        delete [] pdfFlags;
        delete [] inputMeans;
        delete [] inputStdevs;
      }
      return -1;
    }
    dmean = 0.0;
    for (kk = 0; kk < nSamp; kk++)
      if (vecYY[kk] != PSUADE_UNDEFINED) dmean += vecYY[kk];
    dmean /= (double) count;
    dstds = 0.0;
    for (kk = 0; kk < nSamp; kk++)
      if (vecYY[kk] != PSUADE_UNDEFINED) dstds += pow(vecYY[kk]-dmean,2.0);
    dstds /= (double) (count - 1);
    vecBsMeans[nn] = dmean;
    vecBsStds[nn]  = sqrt(dstds);
  }   
  dmean = 0.0; 
  for (nn = 0; nn < ntimes; nn++) dmean += vecBsMeans[nn];
  dmean /= (double) ntimes;
  dstds = 0.0;
  if (ntimes > 1)  
  {
    for (nn = 0; nn < ntimes; nn++) dstds += pow(vecBsMeans[nn]-dmean,2.0);
    dstds = sqrt(dstds/(double) (ntimes - 1));
  } 
  double smean=0, sstd=0;
  smean = 0.0; 
  for (nn = 0; nn < ntimes; nn++) smean += vecBsStds[nn];
  smean /= (double) ntimes;
  sstd = 0.0;
  if (ntimes > 1)  
  {
    for (nn = 0; nn < ntimes; nn++) sstd += pow(vecBsStds[nn]-smean,2.0);
    sstd = sqrt(sstd/(double) (ntimes - 1));
  } 
  if (isScreenDumpModeOn())
  {
    printOutTS(PL_INFO,
      "RSMSobol1: sample mean (std dev of mean) = %10.3e (%10.3e)\n",
      dmean, dstds);
    printOutTS(PL_INFO,
      "RSMSobol1: std dev (std dev of std dev)  = %10.3e (%10.3e)\n",
      smean, sstd);
  }
  if (smean == 0.0) smean = 1.0;
  //save mean & std
  outputMean_ = dmean;
  outputStd_  = smean;
  //cout << outputMean_ << ", " << outputStd_ << endl;

  int       nSteps=1;
  psVector  vecLower2, vecUpper2, vecZZ, samplePtsND;
  psVector  vecInpMeans2, vecInpStdvs2, vecSamPts1D;
  psIVector vecInpFlags2; 
  PDFManager *pdfman1, *pdfman2;

  if (nSamp > 2000) nSteps = nSamp / 1000;
  vecLower2.setLength(nInputs);
  vecUpper2.setLength(nInputs);
  nSamp  = nSubSamples;
  vecSamPts1D.setLength(nLevels*nInputs);
  samplePtsND.setLength(nSubSamples*nInputs*nInputs);
  vecInpFlags2.setLength(nInputs-1);
  vecInpMeans2.setLength(nInputs-1);
  vecInpStdvs2.setLength(nInputs-1);

  if (nLevels > nSubSamples) 
  {
    vecSS.setLength(nLevels);
    vecZZ.setLength(nLevels);
  }
  else
  {
    vecSS.setLength(nSubSamples);
    vecZZ.setLength(nSubSamples);
  }
  for (ii = 0; ii < nInputs; ii++)
  {
    if (noPDF == 0)
    {
      corMat.setDim(1,1);
      corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
      pdfman1 = new PDFManager();
      pdfman1->initialize(1,&pdfFlags[ii],&inputMeans[ii],
                           &inputStdevs[ii],corMat,NULL,NULL);
      vecLB.load(1, &xLower[ii]);
      vecUB.load(1, &xUpper[ii]);
      vecOut.setLength(nLevels);
      pdfman1->genSample(nLevels, vecOut, vecLB, vecUB);
      for (jj = 0; jj < nLevels; jj++)
        vecSamPts1D[ii*nLevels+jj] = vecOut[jj];
      delete pdfman1;
    }
    else
    {
      sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      sampler->setInputBounds(1, &xLower[ii], &xUpper[ii]);
      sampler->setOutputParams(1);
      sampler->setSamplingParams(nLevels, 1, 0);
      sampler->initialize(0);
      tempV = vecSamPts1D.getDVector();
      sampler->getSamples(nLevels, 1, 1, &(tempV[ii*nLevels]), 
                   vecZZ.getDVector(), vecSS.getIVector());
      delete sampler;
    }

    if (noPDF == 0)
    {
      corMat.setDim(nInputs-1, nInputs-1);
      for (jj = 0; jj < ii; jj++)
      {
        vecLower2[jj] = xLower[jj];
        vecUpper2[jj] = xUpper[jj];
        vecInpFlags2[jj] = pdfFlags[jj];
        vecInpMeans2[jj] = inputMeans[jj];
        vecInpStdvs2[jj] = inputStdevs[jj];
        for (kk = 0; kk < ii; kk++)
          corMat.setEntry(jj, kk, corMatp->getEntry(jj,kk));
        for (kk = ii+1; kk < nInputs; kk++)
          corMat.setEntry(jj, kk-1, corMatp->getEntry(jj,kk));
      }
      for (jj = ii+1; jj < nInputs; jj++)
      {
        vecLower2[jj-1] = xLower[jj];
        vecUpper2[jj-1] = xUpper[jj];
        vecInpFlags2[jj-1] = pdfFlags[jj];
        vecInpMeans2[jj-1] = inputMeans[jj];
        vecInpStdvs2[jj-1] = inputStdevs[jj];
        for (kk = 0; kk < ii; kk++)
          corMat.setEntry(jj-1, kk, corMatp->getEntry(jj,kk));
        for (kk = ii+1; kk < nInputs; kk++)
          corMat.setEntry(jj-1, kk-1, corMatp->getEntry(jj,kk));
      }
      pdfman2 = new PDFManager();
      pdfman2->initialize(nInputs-1,vecInpFlags2.getIVector(),
                    vecInpMeans2.getDVector(),vecInpStdvs2.getDVector(),
                    corMat,NULL,NULL);
      vecLB.load(nInputs-1, vecLower2.getDVector());
      vecUB.load(nInputs-1, vecUpper2.getDVector());
      vecOut.setLength(nSubSamples*(nInputs-1));
      pdfman2->genSample(nSubSamples, vecOut, vecLB, vecUB);
      for (jj = 0; jj < nSubSamples*(nInputs-1); jj++) 
        samplePtsND[ii*nSubSamples*nInputs+jj] = vecOut[jj];
      delete pdfman2;
    }
    else
    {
      for (jj = 0; jj < ii; jj++)
      {
        vecLower2[jj] = xLower[jj];
        vecUpper2[jj] = xUpper[jj];
      }
      for (jj = ii+1; jj < nInputs; jj++)
      {
        vecLower2[jj-1] = xLower[jj];
        vecUpper2[jj-1] = xUpper[jj];
      }
      if (nInputs-1 > 51)
           sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      sampler->setInputBounds(nInputs-1, vecLower2.getDVector(), 
                              vecUpper2.getDVector());
      sampler->setOutputParams(1);
      sampler->setSamplingParams(nSubSamples, 1, 0);
      sampler->initialize(0);
      tempV = samplePtsND.getDVector();
      sampler->getSamples(nSubSamples, nInputs-1, 1, 
                &(tempV[ii*nSubSamples*nInputs]), vecZZ.getDVector(), 
                vecSS.getIVector());
      delete sampler;
    }
  }
  if (pdfNull == 1)
  {
    delete [] pdfFlags;
    delete [] inputMeans;
    delete [] inputStdevs;
  }
  if (ioPtr == NULL) delete corMatp;

  int iL, sCnt, totalCnt, offset;
  psVector  vecVceMaxs,vecVceMins,vecVceMeds,vecVces,vecEcvs,vecVars;
  psVector  vecMeans,vecmSamPts;
  psIVector vecBins;

  vecYY.setLength(nSubSamples*nLevels);
  vecBsX.setLength(nSamples*nInputs);
  vecBsY.setLength(nSamples);
  vecBsFlags.setLength(nSamples);
  vecVceMaxs.setLength(nInputs);
  vecVceMins.setLength(nInputs);
  vecVceMeds.setLength(nInputs);
  vecVars.setLength(nLevels*ntimes);
  vecVces.setLength(nInputs*ntimes);
  vecEcvs.setLength(nInputs*ntimes);
  vecBins.setLength(nLevels*ntimes);
  vecMeans.setLength(nLevels*ntimes);
  vecmSamPts.setLength(nInputs*nSubSamples*nLevels);

  for (ii = 0; ii < nInputs; ii++)
  {
    if (isScreenDumpModeOn())
      printOutTS(PL_INFO, "RSMSobol1: processing input %d\n",ii+1);

    for (iL = 0; iL < nLevels; iL++)
    {
      offset = iL * nSubSamples * nInputs;
      tempV = samplePtsND.getDVector();
      for (jj = 0; jj < nSubSamples; jj++)
      {
        oneSamplePt = &(tempV[ii*nSubSamples*nInputs+jj*(nInputs-1)]);
        for (kk = 0; kk < ii; kk++)
          vecmSamPts[offset+jj*nInputs+kk] = oneSamplePt[kk];
        for (kk = ii+1; kk < nInputs; kk++)
          vecmSamPts[offset+jj*nInputs+kk] = oneSamplePt[kk-1];
        vecmSamPts[offset+jj*nInputs+ii] = vecSamPts1D[ii*nLevels+iL];
      }
    }

    for (nn = 0; nn < ntimes; nn++)
    {
      if (ntimes == 1)
      {
        for (ss = 0; ss < nSamples*nInputs; ss++) vecBsX[ss] = X[ss];
        for (ss = 0; ss < nSamples; ss++) 
          vecBsY[ss] = Y2[ss*nOutputs+outputID];
      }
      else
      {
        for (ss = 0; ss < nSamples; ss++) vecBsFlags[ss] = 0;
        count = 0;
        while (count < frac * nSamples)
        {
          jj = PSUADE_rand() % nSamples;
          if (vecBsFlags[jj] == 0)
          {
            for (kk = 0; kk < nInputs; kk++)
              vecBsX[count*nInputs+kk] = X[jj*nInputs+kk];
            vecBsY[count] = Y2[jj*nOutputs+outputID];
            vecBsFlags[jj] = 1;
            count++;
          }
        }
      }
      if (ii == 0 || ntimes > 1) 
        status = faPtr->initialize(vecBsX.getDVector(), vecBsY.getDVector());
      if (isScreenDumpModeOn() && printLevel > 3)
        printOutTS(PL_INFO,"RSMSobol1: function evaluations\n");
      faPtr->evaluatePoint(nSubSamples*nLevels,
                   vecmSamPts.getDVector(),vecYY.getDVector());
      if (nConstr > 0)
      {
        tempV = vecmSamPts.getDVector();
        for (kk = 0; kk < nLevels*nSubSamples; kk++)
        {
          ddata = constrPtr->evaluate(&tempV[kk*nInputs],
                                      vecYY[kk],status);
          if (status == 0) vecYY[kk] = PSUADE_UNDEFINED;
        }
      }
   
      if (isScreenDumpModeOn() && printLevel > 3 &&
         (iL % (nLevels/10) == 0))
        printOutTS(PL_INFO, "RSMSobol1: compute mean and std dev\n");
      for (iL = 0; iL < nLevels; iL++)
      {
        vecMeans[iL*ntimes+nn] = 0.0;
        sCnt = 0;
        for (jj = 0; jj < nSubSamples; jj++)
        {
          if (vecYY[iL*nSubSamples+jj] != PSUADE_UNDEFINED)
          {
            vecMeans[iL*ntimes+nn] += vecYY[iL*nSubSamples+jj];
            sCnt++;
          }
        }
        vecBins[iL*ntimes+nn] = sCnt;
        if (sCnt < nSubSamples/10 && printLevel >= 5)
          printOutTS(PL_INFO,"RSMSobol1 WARNING: subsample size = %d\n",
                     sCnt);
        if (sCnt >= 1) vecMeans[iL*ntimes+nn] /= (double) sCnt;
        else           vecMeans[iL*ntimes+nn] = PSUADE_UNDEFINED;
        vecVars[iL*ntimes+nn] = 0.0;
        if (sCnt > 1)
        {
          for (jj = 0; jj < nSubSamples; jj++)
            if (vecYY[iL*nSubSamples+jj] != PSUADE_UNDEFINED)
              vecVars[iL*ntimes+nn] += pow(vecYY[iL*nSubSamples+jj]-
                                        vecMeans[iL*ntimes+nn],2.0);
          vecVars[iL*ntimes+nn] /= (double) sCnt;
        }
        else vecVars[iL*ntimes+nn] = PSUADE_UNDEFINED;
      }
    }

    for (nn = 0; nn < ntimes; nn++)
    {
      totalCnt = 0;
      for (iL = 0; iL < nLevels; iL++) totalCnt += vecBins[iL*ntimes+nn];
      if (totalCnt == 0) 
      {
        printOutTS(PL_ERROR, "RSMSobol1 ERROR: no feasible region.\n");
        exit(1);
      }

      ddata = 0.0;
      for (iL = 0; iL < nLevels; iL++) 
      {
        if (vecMeans[iL*ntimes+nn] != PSUADE_UNDEFINED)
          ddata += vecMeans[iL*ntimes+nn]*vecBins[iL*ntimes+nn]/totalCnt;
      }
      vecVces[ii*ntimes+nn] = 0.0;
      for (iL = 0; iL < nLevels; iL++)
        if (vecMeans[iL*ntimes+nn] != PSUADE_UNDEFINED)
          vecVces[ii*ntimes+nn] += pow(vecMeans[iL*ntimes+nn]-ddata,2.0) * 
                        vecBins[iL*ntimes+nn] / totalCnt;
      vecEcvs[ii*ntimes+nn] = 0.0;
      for (iL = 0; iL < nLevels; iL++)
      {
        if (vecVars[iL*ntimes+nn] != PSUADE_UNDEFINED)
          vecEcvs[ii*ntimes+nn] += vecVars[iL*ntimes+nn]*
                     vecBins[iL*ntimes+nn]/totalCnt;
      }
    }
  }
  delete faPtr;
  if (constrPtr != NULL) delete constrPtr;

  psVector vecVceMedu;
  vecVceMedu.setLength(nInputs);
  for (ii = 0; ii < nInputs; ii++)
  {
    vecVceMaxs[ii] = -PSUADE_UNDEFINED;
    vecVceMins[ii] =  PSUADE_UNDEFINED;
    vecVceMeds[ii] =  0.0;
    vecVceMedu[ii] =  0.0;
    for (nn = 0; nn < ntimes; nn++)
    {
      if (vecVces[ii*ntimes+nn] > vecVceMaxs[ii])
        vecVceMaxs[ii] = vecVces[ii*ntimes+nn];
      if (vecVces[ii*ntimes+nn] < vecVceMins[ii]) 
        vecVceMins[ii] = vecVces[ii*ntimes+nn];
      vecVceMeds[ii] += vecVces[ii*ntimes+nn];
      vecVceMedu[ii] += 
        (vecVces[ii*ntimes+nn]-vecEcvs[ii*ntimes+nn]/nSubSamples);
    }
    vecVceMeds[ii] /= (double) ntimes;
    vecVceMedu[ii] /= (double) ntimes;
    if (smean != 0)
    {
      vecVceMeds[ii] /= (smean * smean);
      vecVceMedu[ii] /= (smean * smean);
      vecVceMaxs[ii] /= (smean * smean);
      vecVceMins[ii] /= (smean * smean);
    }
    else vecVceMeds[ii]=vecVceMedu[ii]=vecVceMaxs[ii]=vecVceMins[ii]=0;
  }

  vecVces_.setLength(nInputs_);
  if (isScreenDumpModeOn()) printAsterisks(PL_INFO, 0);
  for (ii = 0; ii < nInputs; ii++)
  {
    if (isScreenDumpModeOn())
    {
      printOutTS(PL_INFO,
        "RSMSobol1: Normalized mean VCE for input %3d = %12.4e",
        ii+1, vecVceMeds[ii]);
      if (ntimes > 1)
        printOutTS(PL_INFO,",bounds = [%12.4e, %12.4e]\n",vecVceMins[ii],
                   vecVceMaxs[ii]);
      else printOutTS(PL_INFO, "\n");
    }
    //save vce
    vecVces_[ii] = vecVceMeds[ii];
  }
  if (isScreenDumpModeOn() && printLevel >= 2)
  {
    ddata = ddata2 = 0.0;
    for (ii = 0; ii < nInputs; ii++)
    {
      printOutTS(PL_INFO,"Unnormalized VCE for input %3d = %12.4e\n",ii+1,
                 vecVceMeds[ii]*smean*smean);
      printOutTS(PL_INFO,
           "Unnormalized VCE for input %3d = %12.4e (unbiased)\n",
           ii+1, vecVceMedu[ii]*smean*smean);
      ddata  += vecVceMeds[ii]*smean*smean;
      ddata2 += vecVceMedu[ii]*smean*smean;
    }
    printOutTS(PL_INFO,"Sum of   biased VCEs = %12.4e\n",ddata);
    printOutTS(PL_INFO,"Sum of unbiased VCEs = %12.4e\n",ddata2);
    printOutTS(PL_INFO,"Total variance       = %12.4e\n",smean * smean);
  }

  pData *pObj = NULL;
  if (ioPtr != NULL)
  {
    pObj = ioPtr->getAuxData();
    if (ntimes == 1)
    {
      pObj->nDbles_ = nInputs;
      pObj->dbleArray_ = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++)
        pObj->dbleArray_[ii] = vecVceMeds[ii] * smean * smean;
      pObj->dbleData_ = smean * smean;
    }
    else
    {
      pObj->nDbles_ = 3*nInputs;
      pObj->dbleArray_ = new double[nInputs*3];
      for (ii = 0; ii < nInputs; ii++)
        pObj->dbleArray_[ii] = vecVceMeds[ii] * smean * smean;
      for (ii = 0; ii < nInputs; ii++)
        pObj->dbleArray_[nInputs+ii] = vecVceMins[ii] * smean * smean;
      for (ii = 0; ii < nInputs; ii++)
        pObj->dbleArray_[2*nInputs+ii] = vecVceMaxs[ii] * smean * smean;
      pObj->dbleData_ = smean * smean;
    }
  }

  if (isScreenDumpModeOn() && printLevel >= 1)
  {
    for (ii = 0; ii < nInputs; ii++) vecMeans[ii] = (double) ii;
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,"RSMSobol1: ordered normalized VCE : \n");
    printDashes(PL_INFO, 0);
    sortDbleList2(nInputs,vecVceMeds.getDVector(),vecMeans.getDVector());
    for (ii = nInputs-1; ii >= 0; ii--)
       printOutTS(PL_INFO,"RSMSobol1: Normalized VCE for input %3d = %e\n",
              (int) vecMeans[ii]+1,vecVceMeds[ii]);
    printAsterisks(PL_INFO, 0);
  }
  return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobol1Analyzer& RSMSobol1Analyzer::operator=(const RSMSobol1Analyzer &)
{
  printOutTS(PL_ERROR, 
       "RSMSobol1 operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int RSMSobol1Analyzer::get_nInputs()
{
  return nInputs_;
}
double RSMSobol1Analyzer::get_outputMean()
{
  return outputMean_;
}
double RSMSobol1Analyzer::get_outputStd()
{
  return outputStd_;
}
double RSMSobol1Analyzer::get_vce(int ind)
{
  if (ind < 0 || ind >= nInputs_)
  {
    printf("RSMSobol1 ERROR: get_vce index error.\n");
    return 0.0;
  }
  if (vecVces_.length() <= ind)
  {
    printf("RSMSobol1 ERROR: get_vce has not value.\n");
    return 0.0;
  }
  return vecVces_[ind];
}

