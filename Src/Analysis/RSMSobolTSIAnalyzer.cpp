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
// Functions for the class RSMSobolTSIAnalyzer  
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
#include "Psuade.h"
#include "FuncApprox.h"
#include "Sampling.h"
#include "RSConstraints.h"
#include "PDFManager.h"
#include "PDFNormal.h"
#include "PsuadeData.h"
#include "PsuadeConfig.h"
#include "RSMSobolTSIAnalyzer.h"
#include "sysdef.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer::RSMSobolTSIAnalyzer() : Analyzer(), nInputs_(0), 
                  outputMean_(0), outputStd_(0)
{
  setName("RSMSOBOLTSI");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer::~RSMSobolTSIAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
void RSMSobolTSIAnalyzer::analyze(int nInps, int nSamp, double *lbs,
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
double RSMSobolTSIAnalyzer::analyze(aData &adata)
{
  int    ii, jj, kk, iL, status, nSubSamples=10000, sCnt, nLevels=50;
  int    corFlag, method=1, noPDF, totalCnt, nSamp;
  double variance, dmean, dvar, ddata;
  char   pString[500], *cString, winput1[500], winput2[500];
  Sampling      *sampler;
  FuncApprox    *faPtr;
  RSConstraints *constrPtr;
  PDFManager    *pdfman, *pdfman1, *pdfman2;
  pData         pCorMat;
  psVector  vecIn, vecOut, vecUB, vecLB, vecYT, vecYT2, vecTSIRes;
  psVector  vecSamInps, vecSamOuts,vecVars, vecMeans, vecTSI, vecLower;
  psVector  vecmSamPts, vecInpMeans1, vecInpStdvs1, vecSam1D, vecUpper;
  psIVector vecSamStas, vecBins, vecPDFFlags1;
  psMatrix  *corMatp, corMat, corMat1, corMat2;

  if (method == 1)
  {
    analyze3(adata);
    return 0;
  }

  printAsterisks(PL_INFO, 0);
  printOutTS(PL_INFO, "*          RS-based Total Order Sobol' Indices \n");
  printEquals(PL_INFO, 0); 
  printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET \n");
  printOutTS(PL_INFO,"*\n");
  printOutTS(PL_INFO,
       "* - ana_expert mode to finetune RSMSobolTSI parameters, \n");
  printOutTS(PL_INFO,
       "*   (e.g. sample size for integration can be adjusted).\n");
  printOutTS(PL_INFO,
       "* - rs_expert mode to finetune response surface for RSMSobolTSI,\n");
  printOutTS(PL_INFO,"* - printlevel to display more information.\n");
  printEquals(PL_INFO,0);

  int printLevel  = adata.printLevel_;
  int nInputs     = adata.nInputs_;
  int nOutputs    = adata.nOutputs_;
  int nSamples    = adata.nSamples_;
  int outputID    = adata.outputID_;
  int *pdfFlags   = adata.inputPDFs_;
  double *xLower = adata.iLowerB_;
  double *xUpper = adata.iUpperB_;
  double *XIn    = adata.sampleInputs_;
  double *YIn    = adata.sampleOutputs_;
  double *inputMeans  = adata.inputMeans_;
  double *inputStdevs = adata.inputStdevs_;
  PsuadeData *ioPtr  = adata.ioPtr_;
  noPDF = 1;
  if (pdfFlags != NULL)
  {
    for (ii = 0; ii < nInputs; ii++) if (pdfFlags[ii] != 0) noPDF = 0;
  }
  if (noPDF == 1) 
    printOutTS(PL_INFO, "* RSMSobolTSI INFO: all uniform distributions.\n");
  else
  {
    printOutTS(PL_INFO,"* RSMSobolTSI INFO: non-uniform distributions\n");
    printOutTS(PL_INFO,
         "     detected which will be used in this analysis.\n");
  }

  if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
  {
    printOutTS(PL_ERROR, "RSMSobolTSI ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (nInputs <= 1)
  {
    printOutTS(PL_ERROR,
               "RSMSobolTSI: nInputs<=1 does not need this analysis.\n");
    return PSUADE_UNDEFINED;
  }
  if (outputID >= nOutputs || outputID < 0)
  {
    printOutTS(PL_ERROR,"RSMSobolTSI ERROR: invalid output ID (%d).\n", 
               outputID);
     return PSUADE_UNDEFINED;
  }
  if (ioPtr == NULL)
  {
    printOutTS(PL_ERROR, 
         "RSMSobolTSI ERROR: no data object (PsuadeData) found.\n");
    printOutTS(PL_ERROR, 
         "                   Consult PSUADE developers.\n");
    return PSUADE_UNDEFINED;
  }
  ioPtr->getParameter("input_cor_matrix", pCorMat);
  corMatp = (psMatrix *) pCorMat.psObject_;
  //for (ii = 0; ii < nInputs; ii++)
  //{
  //  for (jj = 0; jj < ii; jj++)
  //  {
  //    if (corMatp->getEntry(ii,jj) != 0.0)
  //    {
  //      printOutTS(PL_ERROR,
  //         "* RSMSobolTSI INFO: this method cannot handle correlated\n");
  //      printOutTS(PL_ERROR,
  //         "*           inputs using joint PDFs yet. Use group\n");
  //      printOutTS(PL_ERROR, "*           variance-based method.\n");
  //      return PSUADE_UNDEFINED;
  //    }
  //  }
  //}
  status = 0;
  for (ii = 0; ii < nSamples; ii++)
    if (YIn[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR, 
         "RSMSobolTSI ERROR: Some outputs are undefined. Prune the\n");
    printOutTS(PL_ERROR,
         "                   undefined sample points first.\n");
    return PSUADE_UNDEFINED;
  }
  vecYT.setLength(nSamples);
  for (ii = 0; ii < nSamples; ii++) vecYT[ii] = YIn[ii*nOutputs+outputID];

  constrPtr = new RSConstraints();
  if (ioPtr != NULL) constrPtr->genConstraints(ioPtr);

  faPtr = genFAInteractive(ioPtr, 0);
  status = faPtr->initialize(XIn, vecYT.getDVector());

  if (psAnaExpertMode_ == 1)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO, "\n");
    printOutTS(PL_INFO,
         "* RSMSobolTSIAnalyzer computes the total sensitivities one\n");
    printOutTS(PL_INFO,
         "*    input at a time. For each input, it first generates a\n");
    printOutTS(PL_INFO,
         "*    sample of size K (that is, K levels). For each level\n");
    printOutTS(PL_INFO,
         "*    a sample of size M is created from all other inputs.\n");
    printOutTS(PL_INFO,
         "*    The total sample size is thus: M * K * nInputs.\n");
    nSubSamples = 10000;
    nLevels = 50;
    printOutTS(PL_INFO,"* nInputs   = %d\n", nInputs);
    printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
    printOutTS(PL_INFO,"* default K = %d\n", nLevels);
    printOutTS(PL_INFO,"* Please enter your desired M and K below.\n");
    printOutTS(PL_INFO,
         "* Recommendation: K should be moderately large since the\n");
    printOutTS(PL_INFO,
         "*     samples are used for computing variances.\n");
    printOutTS(PL_INFO,"* NOTE: large M and K may take a long time\n");
    printEquals(PL_INFO, 0);
    sprintf(pString,"Enter M (suggestion: >= 1000) : ");
    nSubSamples = getInt(1000,100000,pString);
    sprintf(pString,"Enter K (suggestion: >= 50) : ");
    nLevels = getInt(50,1000,pString);
    printAsterisks(PL_INFO, 0);
  }
  else
  {  
    nSubSamples = 10000;
    nLevels = 50;
    if (psConfig_ != NULL)
    {
      cString = psConfig_->getParameter("RSMSoboltsi_nsubsamples");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
        if (nSubSamples < 1000)
        {
          printOutTS(PL_INFO, 
               "RSMSobolTSI INFO: nSubSamples should be >= 1000.\n");
          nSubSamples = 10000;
        }
        else
        {
          printOutTS(PL_INFO,
               "RSMSobolTSI INFO: nSubSamples = %d (config).\n",
               nSubSamples);
        }
      }
      cString = psConfig_->getParameter("RSMSoboltsi_nlevels");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
        if (nLevels < 50)
        {
          printOutTS(PL_INFO,
               "RSMSobolTSI INFO: nLevels should be >= 50.\n");
          nLevels = 50;
        }
        else
        {
          printOutTS(PL_INFO, 
              "RSMSobolTSI INFO: nLevels = %d (config).\n",nLevels);
        }
      }
    }
    if (printLevel > 0)
    {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO,"\n");
      printOutTS(PL_INFO,
         "* RSMSobolTSIAnalyzer computes the total sensitivities one\n");
      printOutTS(PL_INFO,
         "*    input at a time. For each input, it first generates a\n");
      printOutTS(PL_INFO,
         "*    sample of size K (that is, K levels). For each level\n");
      printOutTS(PL_INFO,
         "*    a sample of size M is created from all other inputs.\n");
      printOutTS(PL_INFO,
         "*    The total sample size is thus: M * K * nInputs.\n");
      printOutTS(PL_INFO,"* nInputs   = %d\n", nInputs);
      printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
      printOutTS(PL_INFO,"* default K = %d\n", nLevels);
      printOutTS(PL_INFO,
         "To change these settings, turn on ana_expert mode and rerun.\n");
      printEquals(PL_INFO, 0);
    }
  }

  nSamp = 100000;
  printOutTS(PL_INFO,
       "* RSMSobolTSI INFO: creating a sample for basic statistics.\n");
  printOutTS(PL_INFO,"*                   sample size = %d\n", nSamp);

  vecSamInps.setLength(nSamp*nInputs);
  vecSamOuts.setLength(nSamp);

  pdfman = new PDFManager();
  if (pdfFlags == NULL)
  {
    printOutTS(PL_ERROR,"pdfFlags is NULL in file %s line %d aborting\n",
               __FILE__, __LINE__);
    abort();
  }
  pdfman->initialize(nInputs,pdfFlags,inputMeans,inputStdevs,*corMatp,
                     NULL,NULL);
  vecLB.load(nInputs, xLower);
  vecUB.load(nInputs, xUpper);
  vecOut.setLength(nSamp*nInputs);
  pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
  for (ii = 0; ii < nSamp*nInputs; ii++) vecSamInps[ii] = vecOut[ii];

  printOutTS(PL_INFO, 
       "* RSMSobolTSI: running the sample with response surface...\n");
  faPtr->evaluatePoint(nSamp,vecSamInps.getDVector(),
                       vecSamOuts.getDVector());
  printOutTS(PL_INFO,
       "* RSMSobolTSI: done running the sample with response surface.\n");

  double *oneSamplePt;
  double *dPtr = vecSamInps.getDVector();
  for (ii = 0; ii < nSamp; ii++)
  {
    oneSamplePt = &(dPtr[nInputs*ii]);
    ddata = constrPtr->evaluate(oneSamplePt,vecSamOuts[ii],status);
    if (status == 0) vecSamOuts[ii] = PSUADE_UNDEFINED;
  }

  dmean = 0.0;
  sCnt = 0;
  for (ii = 0; ii < nSamp; ii++)
  {
    if (vecSamOuts[ii] != PSUADE_UNDEFINED)
    {
      dmean += vecSamOuts[ii];
      sCnt++;
    }
  }
  if (sCnt > 1) dmean /= (double) sCnt;
  else
  {
    printOutTS(PL_ERROR,
         "RSMSobolTSI ERROR: too few samples that satisify the ");
    printOutTS(PL_ERROR, "constraints (%d out of %d).\n", 
               sCnt, nSubSamples);
    delete faPtr;
    delete pdfman;
    return PSUADE_UNDEFINED;
  }
  variance = 0.0;
  for (ii = 0; ii < nSamp; ii++)
  {
    if (vecSamOuts[ii] != PSUADE_UNDEFINED)
      variance += (vecSamOuts[ii] - dmean) * (vecSamOuts[ii] - dmean);
  }
  variance /= (double) sCnt;
  printOutTS(PL_INFO,
       "* RSMSobolTSI: mean     (based on N = %d) = %10.3e\n",sCnt,dmean);
  printOutTS(PL_INFO,"* RSMSobolTSI: std dev. (based on N = %d) = %10.3e\n",
             sCnt, sqrt(variance));
  if (variance == 0.0) variance = 1.0;

  vecLower.setLength(nInputs);
  vecUpper.setLength(nInputs);
  vecSamInps.setLength(nSubSamples*nInputs);
  vecSamOuts.setLength(nSubSamples);
  vecSamStas.setLength(nSubSamples);
  vecYT2.setLength(nLevels);
  vecMeans.setLength(nSubSamples);
  vecBins.setLength(nSubSamples);
  vecVars.setLength(nSubSamples);
  vecTSI.setLength(nInputs);
  vecTSIRes.setLength(nInputs);
  vecmSamPts.setLength(nLevels*nInputs);
  vecSam1D.setLength(nLevels);
  vecInpMeans1.setLength(nInputs);
  vecInpStdvs1.setLength(nInputs);
  vecPDFFlags1.setLength(nInputs);

  for (ii = 0; ii < nInputs; ii++)
  {
    printOutTS(PL_INFO, "RSMSobolTSI: processing input %d\n", ii+1);

    corFlag = 0;
    for (jj = 0; jj < nInputs; jj++)
      if (ii != jj && corMatp->getEntry(ii,jj) != 0.0) corFlag = 1;
    if (corFlag == 0)
    {
      for (jj = 0; jj < nInputs; jj++) corFlag += pdfFlags[jj];
      if (corFlag == 0) corFlag = 1;
      else              corFlag = 0;
    }
      
    if (corFlag == 0)
    {
      printOutTS(PL_INFO, "RSMSobolTSI: create samples (1)\n");
      corMat1.setDim(nInputs-1, nInputs-1);
      for (jj = 0; jj < ii; jj++)
      {
        vecLower[jj] = xLower[jj];
        vecUpper[jj] = xUpper[jj];
        vecPDFFlags1[jj] = pdfFlags[jj];
        vecInpMeans1[jj] = inputMeans[jj];
        vecInpStdvs1[jj] = inputStdevs[jj];
        for (kk = 0; kk < ii; kk++)
          corMat1.setEntry(jj, kk, corMatp->getEntry(jj,kk));
        for (kk = ii+1; kk < nInputs; kk++)
          corMat1.setEntry(jj, kk-1, corMatp->getEntry(jj,kk));
      }
      for (jj = ii+1; jj < nInputs; jj++)
      {
        vecLower[jj-1] = xLower[jj];
        vecUpper[jj-1] = xUpper[jj];
        vecPDFFlags1[jj-1] = pdfFlags[jj];
        vecInpMeans1[jj-1] = inputMeans[jj];
        vecInpStdvs1[jj-1] = inputStdevs[jj];
        for (kk = 0; kk < ii; kk++)
          corMat1.setEntry(jj-1, kk, corMatp->getEntry(jj,kk));
        for (kk = ii+1; kk < nInputs; kk++)
          corMat1.setEntry(jj-1, kk-1, corMatp->getEntry(jj,kk));
      }
      pdfman1 = new PDFManager();
      pdfman1->initialize(nInputs-1,vecPDFFlags1.getIVector(),
                          vecInpMeans1.getDVector(),
                          vecInpStdvs1.getDVector(),corMat1,NULL,NULL);
      vecLB.load(nInputs-1, vecLower.getDVector());
      vecUB.load(nInputs-1, vecUpper.getDVector());
      vecOut.setLength(nSubSamples*(nInputs-1));
      pdfman1->genSample(nSubSamples, vecOut, vecLB, vecUB);
      for (jj = 0; jj < nSubSamples*(nInputs-1); jj++)
        vecSamInps[jj] = vecOut[jj];
      delete pdfman1;
      corMat2.setDim(1, 1);
      corMat2.setEntry(0, 0, corMatp->getEntry(0,0));
      pdfman2 = new PDFManager();
      pdfman2->initialize(1, &pdfFlags[ii], &inputMeans[ii],
                          &inputStdevs[ii],corMat2,NULL,NULL);
      vecLB.load(1, &xLower[ii]);
      vecUB.load(1, &xUpper[ii]);
      vecOut.setLength(nLevels);
      pdfman2->genSample(nLevels, vecOut, vecLB, vecUB);
      for (iL = 0; iL < nLevels; iL++) vecSam1D[iL] = vecOut[iL];
      delete pdfman2;
    }
    else
    {
      printOutTS(PL_INFO, "RSMSobolTSI: create samples (2)\n");
      if (nInputs > 51)
           sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      for (jj = 0; jj < ii; jj++)
      {
        vecLower[jj] = xLower[jj];
        vecUpper[jj] = xUpper[jj];
      }
      for (jj = ii+1; jj < nInputs; jj++)
      {
        vecLower[jj-1] = xLower[jj];
        vecUpper[jj-1] = xUpper[jj];
      }
      sampler->setInputBounds(nInputs-1,vecLower.getDVector(),
                              vecUpper.getDVector());
      sampler->setOutputParams(1);
      sampler->setSamplingParams(nSubSamples, 1, 0);
      sampler->initialize(0);
      sampler->getSamples(nSubSamples,nInputs-1,1,vecSamInps.getDVector(),
                          vecSamOuts.getDVector(),vecSamStas.getIVector());
      delete sampler;
      for (iL = 0; iL < nLevels; iL++)
        vecSam1D[iL] = (xUpper[ii] - xLower[ii]) / (nLevels-1) * 
                           iL + xLower[ii];
    }

    for (jj = 0; jj < nSubSamples; jj++)
    {
      for (iL = 0; iL < nLevels; iL++)
      {
        for (kk = 0; kk < ii; kk++)
          vecmSamPts[iL*nInputs+kk] = vecSamInps[jj*(nInputs-1)+kk];
        for (kk = ii+1; kk < nInputs; kk++)
          vecmSamPts[iL*nInputs+kk] = vecSamInps[jj*(nInputs-1)+kk-1];
        vecmSamPts[iL*nInputs+ii] = vecSam1D[iL];
      } 

      if (corFlag == 1)
      {
        vecIn.load(nLevels*nInputs, vecmSamPts.getDVector());
        vecOut.setLength(nLevels*nInputs);
        pdfman->invCDF(nLevels, vecIn, vecOut, vecLB, vecUB);
        for (kk = 0; kk < nLevels*nInputs; kk++) 
          vecmSamPts[kk] = vecOut[kk];
      }

      faPtr->evaluatePoint(nLevels,vecmSamPts.getDVector(),
                           vecYT2.getDVector());

      double *dPtr = vecmSamPts.getDVector();
      for (iL = 0; iL < nLevels; iL++)
      {
        oneSamplePt = &(dPtr[iL*nInputs]);
        ddata = constrPtr->evaluate(oneSamplePt,vecYT2[iL],status);
        if (status == 0) vecYT2[iL] = PSUADE_UNDEFINED;

        if (ddata == PSUADE_UNDEFINED) vecYT2[iL] = PSUADE_UNDEFINED;
      } 

      vecMeans[jj] = 0.0;
      vecVars[jj] = 0.0;
      sCnt = 0;
      for (iL = 0; iL < nLevels; iL++)
      {
        if (vecYT2[iL] != PSUADE_UNDEFINED)
        {
          vecMeans[jj] += vecYT2[iL];
          sCnt++;
        }
      }
      vecBins[jj] = sCnt;
      if (sCnt >= 1) vecMeans[jj] /= (double) sCnt;
      else           vecMeans[jj] = PSUADE_UNDEFINED;
      if (vecMeans[jj] == PSUADE_UNDEFINED) vecVars[jj] = PSUADE_UNDEFINED;
      else
      {
        for (iL = 0; iL < nLevels; iL++)
        {
          if (vecYT2[iL] != PSUADE_UNDEFINED)
            vecVars[jj] += (vecYT2[iL]-vecMeans[jj])*(vecYT2[iL]-vecMeans[jj]);
        }
        if (sCnt == 1) vecVars[jj] = 0.0;
        else           vecVars[jj] = vecVars[jj] / (double) sCnt;
      }
      if (sCnt < nSubSamples/10 && printLevel >= 5)
        printOutTS(PL_DUMP,
             "RSMSobolTSI WARNING: subsample size = %d\n", sCnt);
    }

    totalCnt = 0;
    for (jj = 0; jj < nSubSamples; jj++) totalCnt += vecBins[jj];
    dvar = 0.0;
    for (jj = 0; jj < nSubSamples; jj++)
    {
      if (vecVars[jj] != PSUADE_UNDEFINED)
        dvar += vecVars[jj] * vecBins[jj] / totalCnt;
    }
    vecTSI[ii] = dvar;

    dmean = 0.0;
    for (jj = 0; jj < nSubSamples; jj++)
    {
      if (vecMeans[jj] != PSUADE_UNDEFINED)
        dmean += vecMeans[jj] * vecBins[jj] / totalCnt;
    }
    dvar = 0.0;
    for (jj = 0; jj < nSubSamples; jj++)
      if (vecMeans[jj] != PSUADE_UNDEFINED)
        dvar += (vecMeans[jj] - dmean) * (vecMeans[jj] - dmean) * 
                 vecBins[jj] / totalCnt;
    vecTSIRes[ii] = dvar;

    printOutTS(PL_INFO,
         "RSMSobolTSI (unnormalized) for input %3d = %12.4e\n",
         ii+1,vecTSI[ii]);
    printOutTS(PL_INFO,
         "RSMSobolTSI (  normalized) for input %3d = %12.4e\n",
         ii+1, vecTSI[ii]/variance);
  }
  printAsterisks(PL_INFO, 0);
  for (ii = 0; ii < nInputs; ii++) vecTSI[ii] /= variance;
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO,
           "RSMSobolTSI (normalized) for input %3d = %12.4e\n",
           ii+1,vecTSI[ii]);
  printEquals(PL_INFO, 0);
  if (printLevel > 1)
  {
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, 
           "RSMSobolTSI (unnormalized) for input %3d = %12.4e\n",
           ii+1, vecTSI[ii] * variance);
    for (ii = 0; ii < nInputs; ii++) vecMeans[ii] = (double) ii;
    sortDbleList2(nInputs, vecTSI.getDVector(), vecMeans.getDVector());
    for (ii = nInputs-1; ii >= 0; ii--)
      printOutTS(PL_INFO,
           "RSMSobolTSI (normalized,ordered) for input %3d = %12.4e\n",
           (int) vecMeans[ii]+1,vecTSI[ii]);
    printAsterisks(PL_INFO, 0);
  }
  if (printLevel > 2)
  {
    printAsterisks(PL_INFO, 0);
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, 
           "RSMSobolTSI residual (normalized) for input %3d = %12.4e\n",
           ii+1, vecTSIRes[ii]/variance);
  }
  printAsterisks(PL_INFO, 0);
    
  pData *pPtr = ioPtr->getAuxData();
  pPtr->nDbles_ = nInputs;
  pPtr->dbleArray_ = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) 
    pPtr->dbleArray_[ii] = vecTSI[ii]*variance;
  pPtr->dbleData_ = variance;

  delete constrPtr;
  delete faPtr;
  delete pdfman;
  return 0.0;
}

// ************************************************************************
// perform analysis (version with error bars 12/2012)
// ------------------------------------------------------------------------
double RSMSobolTSIAnalyzer::analyze3(aData &adata)
{
  int    ii, jj, kk, iL, sCnt, status, nSubSamples=10000, nLevels=50;
  int    corFlag, ntimes=1, noPDF, offset, rstype, pdfNull=0, hasSPDF;
  int    totalCnt, nSamp;
  double variance, dmean, dvar, ddata;
  char   pString[500], *cString, winput1[500], winput2[500];
  Sampling      *sampler;
  FuncApprox    *faPtr;
  RSConstraints *constrPtr;
  PDFManager    *pdfman, *pdfman1, *pdfman2;
  psVector  vecIn, vecOut, vecUB, vecLB, vecY, vecInpMeans, vecInpStdvs;
  psVector  vecVars, vecMeans, vecTSI, vecTSIRes, vecSamInps, vecSamOuts;
  psIVector vecSamStas, vecBins, vecInpPDFs;
  pData     pCorMat, pSamFiles, pSamIndices, *pPtr;
  psMatrix  *corMatp, corMat, corMat1, corMat2;

  if (isScreenDumpModeOn())
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO, "*          RS-based Total Order Sobol' Indices \n");
    printEquals(PL_INFO, 0); 
    printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET \n");
    printOutTS(PL_INFO,"*\n");
    printOutTS(PL_INFO,"* - ana_expert mode to finetune internal parameters\n");
    printOutTS(PL_INFO,"*   (e.g. adjust sample size for integration).\n");
    printOutTS(PL_INFO,"* - rs_expert to mode finetune response surface\n");
    printOutTS(PL_INFO,"* - printlevel to 1 to display more information.\n");
    printOutTS(PL_INFO,"* - ntimes=100 to compute error bars for the results\n");
    printEquals(PL_INFO, 0);
  }
  int printLevel = adata.printLevel_;
  int nInputs    = nInputs_ = adata.nInputs_;
  int nOutputs   = adata.nOutputs_;
  int nSamples   = adata.nSamples_;
  int outputID   = adata.outputID_;
  int *pdfFlags  = adata.inputPDFs_;
  double *xLower = adata.iLowerB_;
  double *xUpper = adata.iUpperB_;
  double *XIn    = adata.sampleInputs_;
  double *YIn    = adata.sampleOutputs_;
  PsuadeData *ioPtr   = adata.ioPtr_;

  if (adata.inputMeans_ == NULL || adata.inputPDFs_ == NULL || 
      adata.inputStdevs_ == NULL)
  {
    pdfNull = 1;
    vecInpPDFs.setLength(nInputs);
    vecInpMeans.setLength(nInputs);
    vecInpStdvs.setLength(nInputs);
    for (ii = 0; ii < nInputs; ii++)
    {
      vecInpPDFs[ii] = 0;
      vecInpMeans[ii] = 0;
      vecInpStdvs[ii] = 0;
    }
  }
  else
  {
    if (adata.inputPDFs_ != NULL) 
      vecInpPDFs.load(nInputs,adata.inputPDFs_);
    if (adata.inputMeans_ != NULL) 
      vecInpMeans.load(nInputs,adata.inputMeans_);
    if (adata.inputStdevs_ != NULL) 
      vecInpStdvs.load(nInputs,adata.inputStdevs_);
  }

  noPDF = 1;
  hasSPDF = 0;
  if (pdfFlags != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] != 0) noPDF = 0;
    for (ii = 0; ii < nInputs; ii++)
      if (pdfFlags[ii] == PSUADE_PDF_SAMPLE) hasSPDF = 1;
  }
  if (hasSPDF > 0)
  {
    printOutTS(PL_ERROR,
         "RSMSobolTSI ERROR: S type PDF currently not supported.\n");
    if (ioPtr != NULL && nInputs > 0)
    {
      pPtr = ioPtr->getAuxData();
      pPtr->nDbles_ = nInputs;
      pPtr->dbleArray_ = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++) pPtr->dbleArray_[ii] = 0;
      pPtr->dbleData_ = 0;
    }
    return PSUADE_UNDEFINED;
  }
  if (isScreenDumpModeOn())
  {
    if (noPDF == 1) 
      printOutTS(PL_INFO,"* RSMSobolTSI INFO: all uniform distributions.\n");
    else
    {
      printOutTS(PL_INFO,"* RSMSobolTSI INFO: non-uniform distributions\n");
      printOutTS(PL_INFO,"* detected.\n");
    }
  }

  if (nInputs <= 0 || nSamples <= 0 || nOutputs <= 0) 
  {
    printOutTS(PL_ERROR, "RSMSobolTSI ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
    return PSUADE_UNDEFINED;
  } 
  if (nInputs <= 1)
  {
    printOutTS(PL_ERROR,
         "RSMSobolTSI: analysis not needed for nInputs<=1.\n");
    return PSUADE_UNDEFINED;
  }
  if (outputID >= nOutputs || outputID < 0)
  {
    printOutTS(PL_ERROR, "RSMSobolTSI ERROR: invalid output ID (%d).\n", 
               outputID);
    return PSUADE_UNDEFINED;
  }
  if (ioPtr == NULL)
  {
    if (hasSPDF > 0)
    {
      printOutTS(PL_ERROR,"RSMSobolTSI ERROR: no PsuadeData object found.\n");
      exit(1);
    }
    if (isScreenDumpModeOn())
    {
      printOutTS(PL_INFO,
        "RSMSobolTSI INFO: no data object (PsuadeData) found.\n");
      printOutTS(PL_INFO,"      Several features will be turned off.\n");
      printOutTS(PL_INFO,"      E.g. correlation matrix becomes identity.\n");
    }
    corMatp = new psMatrix();
    corMatp->setDim(nInputs, nInputs);
    for (ii = 0; ii < nInputs; ii++) corMatp->setEntry(ii,ii,1.0e0);
  }
  else
  {
    ioPtr->getParameter("input_cor_matrix", pCorMat);
    corMatp = (psMatrix *) pCorMat.psObject_;
  }
  status = 0;
  for (ii = 0; ii < nSamples; ii++)
    if (YIn[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR,
         "RSMSobolTSI ERROR: Some outputs are undefined. Prune\n");
    printOutTS(PL_ERROR,
         "                   the undefined sample points first.\n");
    return PSUADE_UNDEFINED;
  }

  if (ioPtr != NULL) 
  {
    constrPtr = new RSConstraints();
    constrPtr->genConstraints(ioPtr);
  }
  else
  {
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
    faPtr = genFA(rstype, nInputs, 0, nSamples);
  }
  else faPtr = genFAInteractive(ioPtr, 0);

  faPtr->setBounds(xLower, xUpper);
  vecY.setLength(nSamples);
  for (ii = 0; ii < nSamples; ii++) vecY[ii] = YIn[ii*nOutputs+outputID];
  status = faPtr->initialize(XIn, vecY.getDVector());
  if (status != 0)
  {
    printf("RSMSobolTSI ERROR: failed to build response surface.\n");
    printf("   A suggestion: re-run this with rs_expert mode on\n");
    printf("                 to examine what went wrong.\n");
    if (ioPtr == NULL) delete corMatp;
    return -1;
  }

  if (psAnaExpertMode_ == 1 && isScreenDumpModeOn())
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"\n");
    printOutTS(PL_INFO,
         "* RSMSobolTSIAnalyzer computes the total sensitivities one\n");
    printOutTS(PL_INFO,
         "*    input at a time. For each input, it first generates a\n");
    printOutTS(PL_INFO,
         "*    sample of size K (that is, K levels). For each level\n");
    printOutTS(PL_INFO,
         "*    a sample of size M is created from all other inputs.\n");
    printOutTS(PL_INFO,
         "*    The total sample size is thus: M * K * nInputs.\n");
    nSubSamples = 10000;
    nLevels = 50;
    printOutTS(PL_INFO,"* nInputs   = %d\n", nInputs);
    printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
    printOutTS(PL_INFO,"* default K = %d\n", nLevels);
    printOutTS(PL_INFO,"* Please enter your desired M and K below.\n");
    printOutTS(PL_INFO,
         "* Recommendation: K should be moderately large since the\n");
    printOutTS(PL_INFO,
         "*     samples are used for computing variances.\n");
    printOutTS(PL_INFO,"* NOTE: large M and K may take a long time\n");
    printEquals(PL_INFO, 0);
    sprintf(pString,"Enter M (1000 - 50000) : ");
    nSubSamples = getInt(1000,50000,pString);
    sprintf(pString,"Enter K (50 - 5000) : ");
    nLevels = getInt(50,5000,pString);
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
    nSubSamples = 10000;
    nLevels = 50;
    if (psConfig_ != NULL)
    {
      cString = psConfig_->getParameter("RSMSoboltsi_nsubsamples");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamples);
        if (nSubSamples < 1000)
        {
          printOutTS(PL_INFO,
               "RSMSobolTSI INFO: nSubSamples should be >= 1000.\n");
          nSubSamples = 1000;
        }
        else
        {
          printOutTS(PL_INFO, 
               "RSMSobolTSI INFO: nSubSamples = %d (config).\n",
               nSubSamples);
        }
      }
      cString = psConfig_->getParameter("RSMSoboltsi_nlevels");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nLevels);
        if (nLevels < 50)
        {
          printOutTS(PL_INFO,
               "RSMSobolTSI INFO: nLevels should be >= 50.\n");
          nLevels = 50;
        }
        else
        {
          printOutTS(PL_INFO,
               "RSMSobolTSI INFO: nLevels = %d (config).\n",nLevels);
        }
      }
    }
    if (isScreenDumpModeOn() && printLevel > 0)
    {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO,"\n");
      printOutTS(PL_INFO,
           "* RSMSobolTSIAnalyzer computes the total sensitivities one\n");
      printOutTS(PL_INFO,
           "*    input at a time. For each input, it first generates a\n");
      printOutTS(PL_INFO,
           "*    sample of size K (that is, K levels). For each level\n");
      printOutTS(PL_INFO,
           "*    a sample of size M is created from all other inputs.\n");
      printOutTS(PL_INFO,
           "*    The total sample size is thus: M * K * nInputs.\n");
      printOutTS(PL_INFO,"* nInputs   = %d\n", nInputs);
      printOutTS(PL_INFO,"* default M = %d\n", nSubSamples);
      printOutTS(PL_INFO,"* default K = %d\n", nLevels);
      printOutTS(PL_INFO,
           "To change settings, re-run with ana_expert mode on.\n");
      printAsterisks(PL_INFO, 0);
    }
  }

  nSamp = 200000;
  if (isScreenDumpModeOn())
  {
    printOutTS(PL_INFO,
      "RSMSobolTSI INFO: creating a sample for basic statistics.\n");
    printOutTS(PL_INFO,"                  sample size = %d\n", nSamp);
  }

  vecSamInps.setLength(nSamp*nInputs);
  vecSamOuts.setLength(nSamp);

  pdfman = new PDFManager();
  if (adata.inputPDFs_ == NULL)
  {
    printOutTS(PL_INFO, "pdfFlags is NULL in file %s line %d aborting\n", 
               __FILE__, __LINE__);
    abort();
  }
  if (ioPtr != NULL && hasSPDF > 0) 
  {
    ioPtr->getParameter("input_sample_files", pSamFiles);
    ioPtr->getParameter("input_sample_indices", pSamIndices);
    pdfman->initialize(nInputs,vecInpPDFs.getIVector(),
                       vecInpMeans.getDVector(),vecInpStdvs.getDVector(),
                       *corMatp,pSamFiles.strArray_,pSamIndices.intArray_);
  }
  else
  {
    pdfman->initialize(nInputs,vecInpPDFs.getIVector(),
                       vecInpMeans.getDVector(),vecInpStdvs.getDVector(),
                       *corMatp,NULL,NULL);
  }
  vecLB.load(nInputs, xLower);
  vecUB.load(nInputs, xUpper);
  vecOut.setLength(nSamp*nInputs);
  pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
  for (ii = 0; ii < nSamp*nInputs; ii++) vecSamInps[ii] = vecOut[ii];

  psVector vecSamStdvs;
  vecSamStdvs.setLength(nSamp);

  if (isScreenDumpModeOn() && printLevel > 1)
    printOutTS(PL_INFO,
         "RSMSobolTSI: response surface evaluation begins...\n");

  if (ntimes <= 1)
  {
    faPtr->evaluatePoint(nSamp,vecSamInps.getDVector(),
                         vecSamOuts.getDVector());
      for (ii = 0; ii < nSamp; ii++) vecSamStdvs[ii] = 0.0;
  }
  else
  {
    faPtr->evaluatePointFuzzy(nSamp,vecSamInps.getDVector(),
                   vecSamOuts.getDVector(),vecSamStdvs.getDVector());
  }
  if (isScreenDumpModeOn() && printLevel > 1)
    printOutTS(PL_INFO,
         "RSMSobolTSI: response surface evaluation ends...\n");

  dmean = 0.0;
  for (ii = 0; ii < nSamp; ii++) dmean += vecSamOuts[ii];
  dmean /= (double) nSamp;
  dvar = 0.0;
  for (ii = 0; ii < nSamp; ii++)
    dvar += (vecSamOuts[ii] - dmean) * (vecSamOuts[ii] - dmean) ;
  dvar /= (double) nSamp;
  if (isScreenDumpModeOn() && printLevel > 1)
  {
    printOutTS(PL_INFO,"RSMSobolTSI: sample mean (N=%d) = %e.\n",
               nSamp, dmean);
    printOutTS(PL_INFO,"RSMSobolTSI: sample std  (N=%d) = %e.\n",
               nSamp, sqrt(dvar));
  }

  int    nn, count;
  double d1, d2, *oneSamplePt, *dPtr;
  psIVector vecYcnts;
  psVector  vecYmeans, vecYstds, vecNewY;
  vecYcnts.setLength(ntimes);
  vecYmeans.setLength(ntimes);
  vecYstds.setLength(ntimes);
  vecNewY.setLength(ntimes);
  PDFNormal **rsPDFs = new PDFNormal*[nSamp];
  for (ii = 0; ii < nSamp; ii++)
  {
    if (vecSamStdvs[ii] != 0) 
         rsPDFs[ii] = new PDFNormal(vecSamOuts[ii],vecSamStdvs[ii]); 
    else rsPDFs[ii] = NULL;
  }
  for (nn = 0; nn < ntimes; nn++)
  {
    vecYmeans[nn] = 0.0;
    vecYstds[nn] = 0.0;
    vecYcnts[nn] = 0;
  }
  for (ii = 0; ii < nSamp; ii++)
  {
    if (rsPDFs[ii] != NULL)
    {
      d1 = vecSamOuts[ii] - 4 * vecSamStdvs[ii];
      d2 = vecSamOuts[ii] + 4 * vecSamStdvs[ii];
      rsPDFs[ii]->genSample(ntimes,vecNewY.getDVector(),&d1,&d2);
    }
    else 
    {
      for (nn = 0; nn < ntimes; nn++) vecNewY[nn] = vecSamOuts[ii];
    }
    oneSamplePt = vecSamInps.getDVector();
    oneSamplePt = &(oneSamplePt[ii*nInputs]);
    for (nn = 0; nn < ntimes; nn++)
    {
      status = 1;
      if (constrPtr != NULL)
        ddata = constrPtr->evaluate(oneSamplePt,vecNewY[nn],status);
      if (status != 0)
      {
        vecYmeans[nn] += vecNewY[nn];
        vecYcnts[nn]++;
      }
    }
  }

  for (nn = 0; nn < ntimes; nn++)
    if (vecYcnts[nn] != 0) vecYmeans[nn] /= (double) vecYcnts[nn];
  for (ii = 0; ii < nSamp; ii++)
  {
    if (rsPDFs[ii] != NULL)
    {
      d1 = vecSamOuts[ii] - 4 * vecSamStdvs[ii];
      d2 = vecSamOuts[ii] + 4 * vecSamStdvs[ii];
      rsPDFs[ii]->genSample(ntimes,vecNewY.getDVector(),&d1,&d2);
    } 
    else 
    {
      for (nn = 0; nn < ntimes; nn++) vecNewY[nn] = vecSamOuts[ii];
    }
    oneSamplePt = vecSamInps.getDVector();
    oneSamplePt = &(oneSamplePt[ii*nInputs]);
    for (nn = 0; nn < ntimes; nn++)
    {
      status = 1;
      if (constrPtr != NULL)
        ddata = constrPtr->evaluate(oneSamplePt,vecNewY[nn],status);
      if (status != 0)
        vecYstds[nn] += pow(vecNewY[nn] - vecYmeans[nn], 2.0);
    }
  }
  for (ii = 0; ii < nSamp; ii++)
    if (rsPDFs[ii] != NULL) delete rsPDFs[ii];
  delete [] rsPDFs;

  count  = 0;
  for (nn = 0; nn < ntimes; nn++) count += vecYcnts[nn];
  if (isScreenDumpModeOn())
    printOutTS(PL_INFO,
         "* RSMSobolTSI INFO: %6.2f percent passes the contraints\n",
         (double) count * 100.0 /((double) ntimes*nSamp));
  if (100.0 * count / ((double) ntimes*nSamp) < 10.0)
  {
    printOutTS(PL_ERROR,
         "RSMSobolTSI ERROR: too few samples that satisify the ");
    printOutTS(PL_ERROR, 
         "constraints (%d out of %d).\n", count, nSubSamples);
    delete faPtr;
    delete pdfman;
    return PSUADE_UNDEFINED;
  }
  for (nn = 0; nn < ntimes; nn++)
  {
    if (vecYcnts[nn] > 0) vecYstds[nn] = sqrt(vecYstds[nn]/nSamp);
    else                  vecYstds[nn] = 0.0;
  }
  double smean = 0.0;
  for (nn = 0; nn < ntimes; nn++) smean += vecYstds[nn];
  smean /= (double) ntimes;
  double sstd = 0.0;
  for (nn = 0; nn < ntimes; nn++)
    sstd += (vecYstds[nn] - smean) * (vecYstds[nn] - smean) ;
  sstd = sqrt(sstd/ntimes);

  dmean = 0.0;
  for (nn = 0; nn < ntimes; nn++) dmean += vecYmeans[nn];
  dmean /= (double) ntimes;
  variance = 0.0;
  for (nn = 0; nn < ntimes; nn++)
    variance += (vecYmeans[nn] - dmean) * (vecYmeans[nn] - dmean) ;
  variance /= (double) ntimes;

  if (isScreenDumpModeOn())
  {
    printOutTS(PL_INFO,
       "* RSMSobolTSI: sample mean (std dev of mean) = %10.3e (%10.3e)\n",
       dmean, sqrt(variance));
    printOutTS(PL_INFO,
       "* RSMSobolTSI: std dev (std dev of std dev)  = %10.3e (%10.3e)\n",
       smean, sstd);
  }

  outputMean_ = dmean;
  outputStd_ = smean;
  if (smean == 0.0) smean = 1.0;

  psVector  vecLower, vecUpper, vecSamPt1D, vecSamPtsND;
  psVector  vecInpMeans1, vecInpStdvs1;
  psIVector vecInpPDFs1;
  vecLower.setLength(nInputs);
  vecUpper.setLength(nInputs);
  vecSamOuts.setLength(nSubSamples);
  vecSamStas.setLength(nSubSamples);
  vecSamPt1D.setLength(nLevels*nInputs);
  vecInpMeans1.setLength(nInputs);
  vecInpStdvs1.setLength(nInputs);
  vecInpPDFs1.setLength(nInputs);
  vecSamPtsND.setLength(nSubSamples*nInputs*nInputs);
  int allUPDFs = 0;
  for (ii = 0; ii < nInputs; ii++)
  {
    if (isScreenDumpModeOn())
      printOutTS(PL_INFO,"RSMSobolTSI: processing input %d (phase 1)\n", 
                 ii+1);

    corFlag = 0;
    for (jj = 0; jj < nInputs; jj++)
      if (ii != jj && corMatp->getEntry(ii,jj) != 0.0) corFlag = 1;
    allUPDFs = 0;
    if (corFlag == 0)
    {
      for (jj = 0; jj < nInputs; jj++) allUPDFs += pdfFlags[jj];
      if (allUPDFs == 0) allUPDFs = 1;
    }
      
    if (corFlag == 0 && allUPDFs == 0)
    {
      if (isScreenDumpModeOn())
      {
        printOutTS(PL_DETAIL,
           "RSMSobolTSI: create sample (PDFs not uniform but )\n");
        printOutTS(PL_DETAIL, "no cross correlation)\n");
      }
      corMat1.setDim(nInputs-1, nInputs-1);
      for (jj = 0; jj < ii; jj++)
      {
        vecLower[jj] = xLower[jj];
        vecUpper[jj] = xUpper[jj];
        vecInpPDFs1[jj] = vecInpPDFs[jj];
        vecInpMeans1[jj] = vecInpMeans[jj];
        vecInpStdvs1[jj] = vecInpStdvs[jj];
        for (kk = 0; kk < ii; kk++)
          corMat1.setEntry(jj, kk, corMatp->getEntry(jj,kk));
        for (kk = ii+1; kk < nInputs; kk++)
          corMat1.setEntry(jj, kk-1, corMatp->getEntry(jj,kk));
      }
      for (jj = ii+1; jj < nInputs; jj++)
      {
        vecLower[jj-1] = xLower[jj];
        vecUpper[jj-1] = xUpper[jj];
        vecInpPDFs1[jj-1] = vecInpPDFs[jj];
        vecInpMeans1[jj-1] = vecInpMeans[jj];
        vecInpStdvs1[jj-1] = vecInpStdvs[jj];
        for (kk = 0; kk < ii; kk++)
          corMat1.setEntry(jj-1, kk, corMatp->getEntry(jj,kk));
        for (kk = ii+1; kk < nInputs; kk++)
          corMat1.setEntry(jj-1, kk-1, corMatp->getEntry(jj,kk));
      }
      pdfman1 = new PDFManager();
      pdfman1->initialize(nInputs-1,vecInpPDFs1.getIVector(),
                   vecInpMeans1.getDVector(),vecInpStdvs1.getDVector(),
                   corMat1,NULL,NULL);
      vecLB.load(nInputs-1, vecLower.getDVector());
      vecUB.load(nInputs-1, vecUpper.getDVector());
      vecOut.setLength(nSubSamples*(nInputs-1));
      pdfman1->genSample(nSubSamples, vecOut, vecLB, vecUB);
      for (jj = 0; jj < nSubSamples*(nInputs-1); jj++)
        vecSamPtsND[nSubSamples*ii*nInputs+jj] = vecOut[jj];
      delete pdfman1;
      corMat2.setDim(1, 1);
      corMat2.setEntry(0, 0, corMatp->getEntry(0,0));
      pdfman2 = new PDFManager();
      pdfman2->initialize(1, &(adata.inputPDFs_[ii]), 
                  &(adata.inputMeans_[ii]),&(adata.inputStdevs_[ii]),
                  corMat2,NULL,NULL);
      vecLB.load(1, &xLower[ii]);
      vecUB.load(1, &xUpper[ii]);
      vecOut.setLength(nLevels);
      pdfman2->genSample(nLevels, vecOut, vecLB, vecUB);
      for (iL = 0; iL < nLevels; iL++)
        vecSamPt1D[ii*nLevels+iL] = vecOut[iL];
      delete pdfman2;
    }
    else
    {
      if (isScreenDumpModeOn() && corFlag == 1)
        printOutTS(PL_DETAIL,
          "RSMSobolTSI: create sample (correlation: %d and the rest)\n",
          ii+1);
      if (isScreenDumpModeOn() && allUPDFs == 1)
        printOutTS(PL_DETAIL,
           "RSMSobolTSI: create sample (All uniform PDFs)\n");
      if (nInputs > 51)
           sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      for (jj = 0; jj < ii; jj++)
      {
        vecLower[jj] = xLower[jj];
        vecUpper[jj] = xUpper[jj];
      }
      for (jj = ii+1; jj < nInputs; jj++)
      {
        vecLower[jj-1] = xLower[jj];
        vecUpper[jj-1] = xUpper[jj];
      }
      sampler->setInputBounds(nInputs-1, vecLower.getDVector(), 
                              vecUpper.getDVector());
      sampler->setOutputParams(1);
      sampler->setSamplingParams(nSubSamples, 1, 0);
      sampler->initialize(0);
      dPtr = vecSamPtsND.getDVector();
      sampler->getSamples(nSubSamples, nInputs-1, 1, 
                     &(dPtr[ii*nSubSamples*nInputs]),
                     vecSamOuts.getDVector(), vecSamStas.getIVector());
      delete sampler;
      for (iL = 0; iL < nLevels; iL++)
        vecSamPt1D[iL+ii*nLevels] = (xUpper[ii]-xLower[ii])/(nLevels-1) * 
                                    iL + xLower[ii];
    }
  }

  psVector vecmSamPts, vecYZ, vecTSIMeds, vecTSIMins, vecTSIMaxs;
  psVector vecYFuzzy, vecSFuzzy;
  PDFNormal *rsPDF;
  vecYZ.setLength(nLevels*ntimes);
  vecTSIMeds.setLength(nInputs);
  vecTSIMins.setLength(nInputs);
  vecTSIMaxs.setLength(nInputs);
  vecYFuzzy.setLength(nLevels*nSubSamples);
  vecSFuzzy.setLength(nLevels*nSubSamples);
  vecmSamPts.setLength(nSubSamples*nLevels*nInputs);
  vecMeans.setLength(nSubSamples*ntimes);
  vecBins.setLength(nSubSamples*ntimes);
  vecVars.setLength(nSubSamples*ntimes);
  vecTSI.setLength(nInputs*ntimes);
  vecTSIRes.setLength(nInputs*ntimes);
  for (ii = 0; ii < nInputs; ii++)
  {
    if (isScreenDumpModeOn())
      printOutTS(PL_INFO,"RSMSobolTSI: processing input %d (phase 2)\n", 
                 ii+1);
    if (isScreenDumpModeOn())
      printOutTS(PL_DETAIL,"             Preparing sample\n"); 

    corFlag = 0;
    for (jj = 0; jj < nInputs; jj++)
      if (ii != jj && corMatp->getEntry(ii,jj) != 0.0) corFlag = 1;

    for (jj = 0; jj < nSubSamples; jj++)
    {
      offset = jj * nLevels * nInputs;
      for (iL = 0; iL < nLevels; iL++)
      {
        for (kk = 0; kk < ii; kk++)
          vecmSamPts[offset+iL*nInputs+kk] = 
              vecSamPtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)+kk];
        for (kk = ii+1; kk < nInputs; kk++)
          vecmSamPts[offset+iL*nInputs+kk] = 
              vecSamPtsND[ii*nSubSamples*nInputs+jj*(nInputs-1)+kk-1];
        vecmSamPts[offset+iL*nInputs+ii] = vecSamPt1D[iL+ii*nLevels];
      } 

      if (corFlag == 1)
      {
        dPtr = vecmSamPts.getDVector();
        vecIn.load(nLevels*nInputs, &(dPtr[offset]));
        vecOut.setLength(nLevels*nInputs);
        pdfman->invCDF(nLevels, vecIn, vecOut, vecLB, vecUB);
        for (kk = 0; kk < nLevels*nInputs; kk++) 
          vecmSamPts[offset+kk] = vecOut[kk];
      }
    }

    if (isScreenDumpModeOn())
      printOutTS(PL_DETAIL,"             Evaluating sample\n"); 
    if (ntimes == 1)
    {
      faPtr->evaluatePoint(nLevels*nSubSamples,vecmSamPts.getDVector(),
                           vecYFuzzy.getDVector());
      for (iL = 0; iL < nLevels*nSubSamples; iL++) vecSFuzzy[iL] = 0.0;
    }
    else
    {
      faPtr->evaluatePointFuzzy(nLevels*nSubSamples,vecmSamPts.getDVector(),
                          vecYFuzzy.getDVector(),vecSFuzzy.getDVector());
    }

    if (isScreenDumpModeOn())
      printOutTS(PL_DETAIL,"             Analyzing sample\n"); 
    for (jj = 0; jj < nSubSamples; jj++)
    {
      offset = jj * nLevels;
      dPtr = vecYZ.getDVector();
      for (iL = 0; iL < nLevels; iL++)
      {
        if (vecSFuzzy[offset+iL] != 0)
        {
          rsPDF = new PDFNormal(vecYFuzzy[offset+iL],vecSFuzzy[offset+iL]);
          d1 = vecYFuzzy[offset+iL] - 4 * vecSFuzzy[offset+iL];
          d2 = vecYFuzzy[offset+iL] + 4 * vecSFuzzy[offset+iL];
          rsPDF->genSample(ntimes,&(dPtr[iL*ntimes]),&d1,&d2);
          delete rsPDF;
        }
        else
        {
          for (nn = 0; nn < ntimes; nn++) 
            vecYZ[iL*ntimes+nn] = vecYFuzzy[offset+iL];
        }
      }

      dPtr = vecmSamPts.getDVector();
      for (iL = 0; iL < nLevels; iL++)
      {
        oneSamplePt = &(dPtr[jj*nLevels*nInputs+iL*nInputs]);
        for (nn = 0; nn < ntimes; nn++)
        {
          status = 1;
          if (constrPtr != NULL)
            ddata = constrPtr->evaluate(oneSamplePt,vecYZ[iL*ntimes+nn],
                                        status);
          if (status == 0) vecYZ[iL*ntimes+nn] = PSUADE_UNDEFINED;
        }
      }

      for (nn = 0; nn < ntimes; nn++)
      {
        vecMeans[jj*ntimes+nn] = 0.0;
        vecVars[jj*ntimes+nn] = 0.0;
        sCnt = 0;
        for (iL = 0; iL < nLevels; iL++)
        {
          if (vecYZ[iL*ntimes+nn] != PSUADE_UNDEFINED)
          {
            vecMeans[jj*ntimes+nn] += vecYZ[iL*ntimes+nn];
            sCnt++;
          }
        }
        vecBins[jj*ntimes+nn] = sCnt;
        if (sCnt >= 1) vecMeans[jj*ntimes+nn] /= (double) sCnt;
        else           vecMeans[jj*ntimes+nn] = PSUADE_UNDEFINED;
        if (vecMeans[jj*ntimes+nn] == PSUADE_UNDEFINED)
          vecVars[jj*ntimes+nn] = PSUADE_UNDEFINED;
        else
        {
          for (iL = 0; iL < nLevels; iL++)
          {
            if (vecYZ[iL*ntimes+nn] != PSUADE_UNDEFINED)
              vecVars[jj*ntimes+nn] += 
                 pow(vecYZ[iL*ntimes+nn]-vecMeans[jj*ntimes+nn],2.0);
          }
          if (sCnt == 1) vecVars[jj*ntimes+nn] = 0.0;
          else           vecVars[jj*ntimes+nn] /= (double) sCnt;
        }
      }
    }

    for (nn = 0; nn < ntimes; nn++)
    {
      totalCnt = 0;
      for (jj = 0; jj < nSubSamples; jj++) totalCnt += vecBins[jj*ntimes+nn];
      if (totalCnt == 0)
      {
        printOutTS(PL_ERROR, "RSMSobolTSI ERROR: no feasible region.\n");
        exit(1);
      }
      dvar = 0.0;
      for (jj = 0; jj < nSubSamples; jj++)
      {
        if (vecVars[jj*ntimes+nn] != PSUADE_UNDEFINED)
          dvar += vecVars[jj*ntimes+nn] * vecBins[jj*ntimes+nn] / totalCnt;
      }
      vecTSI[ii*ntimes+nn] = dvar;

      dmean = 0.0;
      for (jj = 0; jj < nSubSamples; jj++)
      {
        if (vecMeans[jj*ntimes+nn] != PSUADE_UNDEFINED)
          dmean += vecMeans[jj*ntimes+nn]*vecBins[jj*ntimes+nn]/totalCnt;
      }
      dvar = 0.0;
      for (jj = 0; jj < nSubSamples; jj++)
        if (vecMeans[jj*ntimes+nn] != PSUADE_UNDEFINED)
          dvar += pow(vecMeans[jj*ntimes+nn]-dmean, 2.0)*vecBins[jj]/totalCnt;
      vecTSIRes[ii*ntimes+nn] = dvar;
    }
  }

  for (ii = 0; ii < nInputs; ii++)
  {
    vecTSIMaxs[ii] = -PSUADE_UNDEFINED;
    vecTSIMins[ii] =  PSUADE_UNDEFINED;
    vecTSIMeds[ii] =  0.0;
    for (nn = 0; nn < ntimes; nn++)
    {
      if (vecTSI[ii*ntimes+nn] > vecTSIMaxs[ii]) 
        vecTSIMaxs[ii] = vecTSI[ii*ntimes+nn];
      if (vecTSI[ii*ntimes+nn] < vecTSIMins[ii]) 
        vecTSIMins[ii] = vecTSI[ii*ntimes+nn];
      vecTSIMeds[ii] += vecTSI[ii*ntimes+nn];
    }
    vecTSIMeds[ii] /= (double) ntimes;
    if (smean != 0)
    {
      vecTSIMeds[ii] /= (smean * smean);
      vecTSIMaxs[ii] /= (smean * smean);
      vecTSIMins[ii] /= (smean * smean);
    }
    else vecTSIMeds[ii] = vecTSIMaxs[ii] = vecTSIMins[ii] = 0.0;

    if (isScreenDumpModeOn())
    {
      printOutTS(PL_INFO, "RSMSobolTSI (normalized) for input %3d = %12.4e",
                 ii+1,vecTSIMeds[ii]);
      if (ntimes > 1)
        printOutTS(PL_INFO, ", bounds = [%12.4e, %12.4e]\n",
                   vecTSIMins[ii],vecTSIMaxs[ii]);
      else printOutTS(PL_INFO, "\n");
    }
  }
  if (ntimes == 1 && ioPtr != NULL)
  {
    pPtr = ioPtr->getAuxData();
    pPtr->nDbles_ = nInputs;
    pPtr->dbleArray_ = new double[nInputs];
    for (ii = 0; ii < nInputs; ii++) 
      pPtr->dbleArray_[ii] = vecTSIMeds[ii] * smean * smean;
    pPtr->dbleData_ = smean*smean;
  }
  if (isScreenDumpModeOn()) printAsterisks(PL_INFO, 0);
    
  if (ntimes == 1)
  {
    for (ii = 0; ii < nInputs; ii++) vecTSIMeds[ii] *= smean * smean;
    if (isScreenDumpModeOn())
      printResults(nInputs,smean*smean,vecTSIMeds.getDVector(),NULL,
                   NULL,ioPtr,0);
  }
  else
  {
    for (ii = 0; ii < nInputs; ii++)
    {
      vecTSIMeds[ii] *= smean * smean;
      vecTSIMins[ii] *= smean * smean;
      vecTSIMaxs[ii] *= smean * smean;
    }
    if (isScreenDumpModeOn())
      printResults(nInputs,smean*smean,vecTSIMeds.getDVector(),
            vecTSIMins.getDVector(),vecTSIMaxs.getDVector(),ioPtr,1);
  }

  vecTSIs_.setLength(nInputs_);
  for (ii = 0; ii < nInputs; ii++)
    vecTSIs_[ii] = vecTSIMeds[ii]/smean/smean;

  if (isScreenDumpModeOn() && printLevel > 1)
  {
    printEquals(PL_INFO, 0);
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, 
           "RSMSobolTSI (unnormalized) for input %3d = %12.4e\n",
           ii+1, vecTSIMeds[ii]);
    for (ii = 0; ii < nInputs; ii++) vecMeans[ii] = (double) ii;
    sortDbleList2(nInputs,vecTSIMeds.getDVector(),vecMeans.getDVector());
    for (ii = nInputs-1; ii >= 0; ii--)
      printOutTS(PL_INFO, 
           "RSMSobolTSI (unnormalized,ordered) for input %3d = %12.4e\n",
           (int) vecMeans[ii]+1,vecTSIMeds[ii]);
    printAsterisks(PL_INFO, 0);
  }

  if (constrPtr != NULL) delete constrPtr;
  delete faPtr;
  delete pdfman;
  if (ioPtr == NULL) delete corMatp;
  return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobolTSIAnalyzer& RSMSobolTSIAnalyzer::operator=(const RSMSobolTSIAnalyzer&)
{
   printOutTS(PL_ERROR,
              "RSMSobolTSI operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// print result
// ------------------------------------------------------------------------
int RSMSobolTSIAnalyzer::printResults(int nInputs, double variance,
                                      double *tsi, double *tsiMins, 
                                      double *tsiMaxs, PsuadeData *ioPtr,
                                      int flag)
{
   int   ii;
   FILE  *fp;
   char  **iNames;
   pData qData;

   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);
   if (qData.strArray_ != NULL) iNames = qData.strArray_;
   else                         iNames = NULL;
   printEquals(PL_INFO, 0);
   if (variance == 0.0)
   {
      printOutTS(PL_INFO, 
           "Total variance = 0. Hence, no total effect plot.\n");
      return 0;
   }

   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO, "Total Effect Statistics: \n");
   for (ii = 0; ii < nInputs; ii++)
   {
      printOutTS(PL_INFO, "Input %4d: Sobol' total sensitivity = %12.4e",
                 ii+1,tsi[ii]/variance);
      if (flag == 1)
           printOutTS(PL_INFO, ", bounds = [%12.4e, %12.4e]\n",
                      tsiMins[ii]/variance, tsiMaxs[ii]/variance);
      else printOutTS(PL_INFO, "\n");
   }
   if (psPlotTool_ == 1) fp = fopen("scilabrssoboltsi.sci", "w");
   else                  fp = fopen("matlabrssoboltsi.m", "w");
   if (fp != NULL)
   {
      if (psPlotTool_ == 1)
      {
         fprintf(fp, "// This file contains Sobol' total indices\n");
         fprintf(fp, "// set sortFlag = 1 and set nn to be the number\n");
         fprintf(fp, "// of inputs to display.\n");
      }
      else
      {
         fprintf(fp, "%% This file contains Sobol' total indices\n");
         fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
         fprintf(fp, "%% of inputs to display.\n");
      }
      fprintf(fp, "sortFlag = 0;\n");
      fprintf(fp, "nn = %d;\n", nInputs);
      fprintf(fp, "var = %e;\n", variance);
      fprintf(fp, "Mids = [\n");
      for (ii = 0; ii < nInputs; ii++) 
         fprintf(fp,"%24.16e\n",tsi[ii]/variance);
      fprintf(fp, "];\n");
      if (flag == 1)
      {
         fprintf(fp, "Mins = [\n");
         for (ii = 0; ii < nInputs; ii++)
            fprintf(fp,"%24.16e\n",tsiMins[ii]/variance);
         fprintf(fp, "];\n");
         fprintf(fp, "Maxs = [\n");
         for (ii = 0; ii < nInputs; ii++) 
            fprintf(fp,"%24.16e\n",tsiMaxs[ii]/variance);
         fprintf(fp, "];\n");
      }
      if (iNames == NULL)
      {
         fprintf(fp, "Str = {");
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
         fprintf(fp,"'X%d'};\n",nInputs);
      }
      else
      {
         fprintf(fp, "Str = {");
         for (ii = 0; ii < nInputs-1; ii++)
         {
            if (iNames[ii] != NULL) fprintf(fp,"'%s',",iNames[ii]);
            else                    fprintf(fp,"'X%d',",ii+1);
         }
         if (iNames[nInputs-1] != NULL) 
              fprintf(fp,"'%s'};\n",iNames[nInputs-1]);
         else fprintf(fp,"'X%d'};\n",nInputs);
      }
      fwritePlotCLF(fp);
      fprintf(fp, "if (sortFlag == 1)\n");
      if (psPlotTool_ == 1)
           fprintf(fp, "  [Mids, I2] = gsort(Mids);\n");
      else fprintf(fp, "  [Mids, I2] = sort(Mids,'descend');\n");
      if (flag == 1)
      {
         fprintf(fp, "  Maxs = Maxs(I2);\n");
         fprintf(fp, "  Mins = Mins(I2);\n");
      }
      fprintf(fp, "  Str  = Str(I2);\n");
      fprintf(fp, "  I2 = I2(1:nn);\n");
      fprintf(fp, "  Mids = Mids(1:nn);\n");
      if (flag == 1)
      {
         fprintf(fp, "  Maxs = Maxs(1:nn);\n");
         fprintf(fp, "  Mins = Mins(1:nn);\n");
      }
      fprintf(fp, "  Str  = Str(1:nn);\n");
      fprintf(fp, "end\n");
      if (flag == 1)
      {
         fprintf(fp, "ymin = min(Mins);\n");
         fprintf(fp, "ymax = max(Maxs);\n");
      }
      else
      {
         fprintf(fp, "ymin = min(Mids);\n");
         fprintf(fp, "ymax = max(Mids);\n");
      }
      fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
      fprintf(fp, "bar(Mids,0.8);\n");
      if (flag == 1)
      {
         fprintf(fp, "for ii = 1:nn\n");
         fprintf(fp, "   if (ii == 1)\n");
         if (psPlotTool_ == 1)
              fprintf(fp, "      set(gca(),\"auto_clear\",\"off\")\n");
         else fprintf(fp, "      hold on\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "   XX = [ii ii];\n");
         fprintf(fp, "   YY = [Mins(ii) Maxs(ii)];\n");
         fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
         fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',12)\n");
         fprintf(fp, "end;\n");
      }
      fwritePlotAxes(fp);
      fprintf(fp, "ymin=0;\n");
      if (psPlotTool_ == 1)
      {
         fprintf(fp, "a=gca();\n");
         fprintf(fp, "a.data_bounds=[0, 0; nn+1, ymax];\n");
         fprintf(fp, "newtick = a.x_ticks;\n");
         fprintf(fp, "newtick(2) = [1:nn]';\n");
         fprintf(fp, "newtick(3) = Str';\n");
         fprintf(fp, "a.x_ticks = newtick;\n");
         fprintf(fp, "a.x_label.font_size = 3;\n");
         fprintf(fp, "a.x_label.font_style = 4;\n");
      }
      else
      {
         fprintf(fp, "axis([0 nn+1 0 ymax])\n");
         fprintf(fp, "set(gca,'XTickLabel',[]);\n");
         fprintf(fp, "th=text(1:nn, repmat(ymin-0.07*(ymax-ymin),nn,1),Str,");
         fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
         fprintf(fp, "set(th, 'fontsize', 12)\n");
         fprintf(fp, "set(th, 'fontweight', 'bold')\n");
      }
      fwritePlotTitle(fp,"Sobol Total Order Indices");
      fwritePlotYLabel(fp, "Sobol Indices");
      if (psPlotTool_ == 1)
           fprintf(fp, "      set(gca(),\"auto_clear\",\"on\")\n");
      else fprintf(fp, "      hold off\n");
      fclose(fp);
      if (psPlotTool_ == 1)
           printOutTS(PL_INFO, 
                "RSMSobolTSI plot file = scilabrssoboltsi.sci\n");
      else printOutTS(PL_INFO, 
                "RSMSobolTSI plot file = matlabrssoboltsi.m\n");
      return 0;
   }
   else
   {
      printOutTS(PL_ERROR,
           "RSMSobolTSI ERROR: cannot create rssoboltsi graphics.\n");
      return 0;
   }
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int RSMSobolTSIAnalyzer::get_nInputs()
{
   return nInputs_;
}
double RSMSobolTSIAnalyzer::get_outputMean()
{
   return outputMean_;
}
double RSMSobolTSIAnalyzer::get_outputStd()
{
   return outputStd_;
}
double RSMSobolTSIAnalyzer::get_tsi(int ind)
{
   if (ind < 0 || ind >= nInputs_)
   {
      printf("RSMSobolTSI ERROR: get_tsi index error.\n");
      return 0.0;
   }
   if (vecTSIs_.length() == 0)
   {
      printf("RSMSobolTSI ERROR: get_tsi has not value.\n");
      return 0.0;
   }
   return vecTSIs_[ind];
}

