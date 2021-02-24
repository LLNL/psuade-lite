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
// Functions for the class RSMSobol2Analyzer  
// (Sobol' second order sensitivity analysis - with response surface)
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Vector.h"
#include "Matrix.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "Sampling.h"
#include "RSConstraints.h"
#include "PDFManager.h"
#include "pData.h"
#include "PsuadeData.h"
#include "PsuadeConfig.h"
#include "RSMSobol2Analyzer.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobol2Analyzer::RSMSobol2Analyzer() : Analyzer(),nInputs_(0),outputMean_(0),
                   outputStd_(0), vces_(0), ecvs_(0)
{
   setName("RSMSOBOL2");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobol2Analyzer::~RSMSobol2Analyzer()
{
   if (vces_ != NULL) delete [] vces_;
   if (ecvs_ != NULL) delete [] ecvs_;
}

// ************************************************************************
// perform analysis (this is intended for library calls)
// ------------------------------------------------------------------------
void RSMSobol2Analyzer::analyze(int nInps, int nSamp, double *lbs,
                                double *ubs, double *X, double *Y)
{
   int    ii, *pdfFlags;
   double *inputMeans, *inputStdevs;

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
   pdfFlags    = new int[nInps];
   inputMeans  = new double[nInps];
   inputStdevs = new double[nInps];
   checkAllocate(inputStdevs, "inputStdevs in RSMSobol2::analyze (lib)");
   for (ii = 0; ii < nInps; ii++)
   {
      pdfFlags[ii] = 0;
      inputMeans[ii]  = 0;
      inputStdevs[ii] = 0;
   }
   adata.inputPDFs_ = pdfFlags;
   adata.inputMeans_ = inputMeans;
   adata.inputStdevs_ = inputStdevs;
   analyze2(adata);
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobol2Analyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ii, jj, kk, status, outputID;
   int    nSubSamples=100, iL, sCnt, *SS, nLevels=100, noPDF=1, pdfNull=0;
   int    ii2, iR, currNLevels, nSamp, printLevel, *pdfFlags, rstype;
   int    *bins, totalCnt, *pdfFlags1, *pdfFlags2, selectedInput=-1;
   double *xLower, *xUpper, *X, *Y, *Y2, *YY, *cLower, *cUpper, *XX;
   double *oneSamplePt, *means, *vars, vce, ddata, variance, *ZZ;
   double dmean, ecv, *inputMeans, *inputStdevs, *mSamplePts;
   double *samplePts2D, *inputMeans1, *inputMeans2, *inputStdevs1;
   double *inputStdevs2;
   char   pString[500], *cString, winput[500], winput2[500];
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   psVector      vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat, pPDF;
   psMatrix      *corMatp, corMat;
   PDFManager    *pdfman, *pdfman1, *pdfman2;
   Sampling      *sampler;

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"*          RS-based Second Order Sobol' Indices \n");
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET\n");
   printOutTS(PL_INFO,"*\n");
   printOutTS(PL_INFO,
        "* - ana_expert mode to finetune RSMSobol2 parameters \n");
   printOutTS(PL_INFO,
        "*   (e.g. sample size for integration can be adjusted).\n");
   printOutTS(PL_INFO,"* - rs_expert mode to finetune response surface\n");
   printOutTS(PL_INFO,
        "* - printlevel to 1 or higher to display more information\n");
   printEquals(PL_INFO, 0);
   
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
      checkAllocate(inputStdevs, "inputStdevs in RSMSobol2::analyze");
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
   }
   if (noPDF == 1) 
      printOutTS(PL_INFO,"* RSMSobol2 INFO: all uniform distributions.\n");
   else
   {
      printOutTS(PL_INFO,
           "RSMSobol2 INFO: non-uniform distributions detected -\n");
      printOutTS(PL_INFO,"                will be used in this analysis.\n");
   }

   if (nInputs <= 1 || nSamples <= 0 || nOutputs <= 0)
   {
      printOutTS(PL_ERROR,"RSMSobol2 ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR,"   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR,"   nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR,"   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   }
   if (nInputs <= 2)
   {
      printOutTS(PL_ERROR,
           "RSMSobol2 ERROR: no need for this analysis (nInputs<=2).\n");
      return PSUADE_UNDEFINED;
   }
   if (outputID < 0 || outputID >= nOutputs)
   {
      printOutTS(PL_ERROR,
           "RSMSobol2 ERROR: invalid outputID (%d).\n",outputID);
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printOutTS(PL_INFO,
           "RSMSobol2 INFO: no data object (PsuadeData) found.\n");
      printOutTS(PL_INFO,"          Several features will be turned off.\n");
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
               printOutTS(PL_INFO,
                    "RSMSobol2 INFO: this method cannot handle correlated\n");
               printOutTS(PL_INFO,
                    "          inputs using joint PDFs. PSUADE will try\n");
               printOutTS(PL_INFO,
                    "          a variant of this method, or you can re-run\n");
               printOutTS(PL_INFO,
                    "          using the group variance-based method.\n");
               if (pdfNull == 1)
               {
                  delete [] pdfFlags;
                  delete [] inputMeans;
                  delete [] inputStdevs;
               }
               return analyze2(adata);
            }
         }
         if (pdfFlags[ii] == PSUADE_PDF_SAMPLE)
         {
            printOutTS(PL_ERROR,
                 "RSMSobol2 INFO: this method cannot handle S PDF type.\n");
            printOutTS(PL_INFO,
                 "          PSUADE will try a variant of this method.\n");
            if (pdfNull == 1)
            {
               delete [] pdfFlags;
               delete [] inputMeans;
               delete [] inputStdevs;
            }
            return analyze2(adata);
         }
      }
   }
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printOutTS(PL_ERROR,
                 "RSMSobol2 ERROR: Some outputs are undefined. Prune\n");
      printOutTS(PL_ERROR,
                 "                 the undefined sample points first.\n");
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
      printf("RSMSobolTSI INFO: no PsuadeData ==> no constraints.\n");
   }

   if (ioPtr == NULL)
   {
      printf("Select response surface. Options are: \n");
      writeFAInfo(0);
      strcpy(pString, "Choose response surface: ");
      rstype = getInt(0, PSUADE_NUM_RS, pString);
      faPtr = genFA(rstype, nInputs, 0, nSamples);
   }
   else faPtr = genFAInteractive(ioPtr, 0);
   faPtr->setBounds(xLower, xUpper);
   Y = new double[nSamples];
   checkAllocate(Y, "Y in RSMSobol2::analyze");
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];
   status = faPtr->initialize(X, Y);

   printAsterisks(PL_INFO, 0);
   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO,
                 "* RSMSobol2 generates a mesh of size K x K for every\n");
      printOutTS(PL_INFO,
                 "*   pair of inputs and then creates a sample of size\n");
      printOutTS(PL_INFO,
                 "*   M for each mesh point. The total sample size is:\n");
      printOutTS(PL_INFO,
                 "*      N = M * K * K * nInputs * (nInputs - 1) / 2.\n");
      printOutTS(PL_INFO,"* NOW, nInputs = %d\n", nInputs);
      printOutTS(PL_INFO,"* Please select M and K below.\n");
      printOutTS(PL_INFO,"* Recommendation: K x K >> M.\n");
      printOutTS(PL_INFO,"* NOTE: large M and K can take a long time.\n");
      printOutTS(PL_INFO,"* Default M = %d\n", nSubSamples);
      printOutTS(PL_INFO,"* Default K = %d\n", nLevels);
      printEquals(PL_INFO, 0);
      sprintf(pString,"Enter M (suggestion: 100 - 1000) : ");
      nSubSamples = getInt(100, 1000, pString);
      sprintf(pString, "Enter K (suggestion: 50 - 500) : ");
      nLevels = getInt(50, 500, pString);

   }
   else
   {
      nSubSamples = 100;
      nLevels = 100;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSobol2_nsubsamples");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput, winput2, &nSubSamples);
            if (nSubSamples < 100)
            {
               printOutTS(PL_INFO, 
                    "RSMSobol2 INFO: nSubSamples should be >= 100.\n");
               nSubSamples = 100;
            }
            else
            {
               printOutTS(PL_INFO, 
                  "RSMSobol2 INFO: nSubSamples = %d (config).\n",
                  nSubSamples);
            }
         }
         cString = psConfig_->getParameter("RSMSobol2_nlevels");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput, winput2, &nLevels);
            if (nLevels < 100)
            {
               printOutTS(PL_INFO, 
                    "RSMSobol2 INFO: nLevels should be >= 100.\n");
               nLevels = 100;
            }
            else
            {
               printOutTS(PL_INFO, 
                    "RSMSobol2 INFO: nLevels = %d (config).\n",nLevels);
            }
         }
      }
      if (printLevel > 0)
      {
         printOutTS(PL_INFO,"RSMSobol2: default M = %d.\n", nSubSamples);
         printOutTS(PL_INFO,"RSMSobol2: default K = %d.\n", nLevels);
         printOutTS(PL_INFO,
              "To change these settings, re-run with ana_expert mode on.\n");
      }
   }
   printEquals(PL_INFO, 0);

   nSamp = 100000;
   printOutTS(PL_INFO,
        "RSMSobol2 INFO: creating a sample for basic statistics.\n");
   printOutTS(PL_INFO,"                sample size = %d\n", nSamp);

   XX = new double[nSamp*nInputs];
   YY = new double[nSamp];
   checkAllocate(YY, "YY in RSMSobol2::analyze");

   if (noPDF == 0)
   {
      pdfman = new PDFManager();
      pdfman->initialize(nInputs,pdfFlags,inputMeans,
                         inputStdevs,*corMatp,NULL,NULL);
      vecLB.load(nInputs, xLower);
      vecUB.load(nInputs, xUpper);
      vecOut.setLength(nSamp*nInputs);
      pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
      for (ii = 0; ii < nSamp*nInputs; ii++) XX[ii] = vecOut[ii];
      delete pdfman;
   }
   else
   {
      if (nInputs < 51) sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      else              sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      sampler->setInputBounds(nInputs, xLower, xUpper);
      sampler->setOutputParams(1);
      sampler->setSamplingParams(nSamp, 1, 1);
      sampler->initialize(0);
      SS = new int[nSamp];
      checkAllocate(SS, "SS in RSMSobol2::analyze");
      sampler->getSamples(nSamp, nInputs, 1, XX, YY, SS);
      delete [] SS;
      delete sampler;
   }

   printOutTS(PL_INFO,
        "RSMSobol2: running the sample with response surface...\n");
   faPtr->evaluatePoint(nSamp, XX, YY);
   printOutTS(PL_INFO,
        "RSMSobol2: done running the sample with response surface.\n");
   
   for (ii = 0; ii < nSamp; ii++)
   {
      oneSamplePt = &(XX[ii*nInputs]);
      status = 1;
      if (constrPtr != NULL)
         ddata = constrPtr->evaluate(oneSamplePt,YY[ii],status);
      if (status == 0) YY[ii] = PSUADE_UNDEFINED;
   }
   
   dmean = 0.0;
   sCnt = 0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
      {
         dmean += YY[ii];
         sCnt++;
      }
   }
   if (sCnt > 1) dmean /= (double) sCnt;
   else
   {
      printOutTS(PL_ERROR, 
           "RSMSobol2 ERROR: too few samples that satisify\n");
      printOutTS(PL_ERROR, "constraints (%d out of %d).\n",sCnt,nSamp);
      delete [] XX;
      delete [] YY;
      delete faPtr;
      if (ioPtr == NULL) delete corMatp;
      if (constrPtr != NULL) delete constrPtr;
      if (pdfNull == 1)
      {
         delete [] pdfFlags;
         delete [] inputMeans;
         delete [] inputStdevs;
      }
      return PSUADE_UNDEFINED;
   }
   variance = 0.0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
         variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
   }
   variance /= (double) sCnt;
   printOutTS(PL_INFO,
        "RSMSobol2: sample mean    (based on N = %d) = %10.3e\n",
        sCnt, dmean);
   printOutTS(PL_INFO,
        "RSMSobol2: sample std dev (based on N = %d) = %10.3e\n",
        sCnt, sqrt(variance));
   if (variance == 0.0) variance = 1.0;
   delete [] XX;
   delete [] YY;

   //save mean & std
   outputMean_ = dmean;
   outputStd_ = sqrt(variance);

   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   nSamp  = nSubSamples;
   XX     = new double[nSubSamples*nInputs];
   YY     = new double[nSubSamples];
   means  = new double[nLevels*nLevels];
   vars   = new double[nLevels*nLevels];
   bins   = new int[nLevels*nLevels];
   pdfFlags1    = new int[2];
   inputMeans1  = new double[2];
   inputStdevs1 = new double[2];
   pdfFlags2    = new int[nInputs-2];
   inputMeans2  = new double[nInputs-2];
   inputStdevs2 = new double[nInputs-2];
   samplePts2D  = new double[nLevels*nLevels*2];
   mSamplePts   = new double[nInputs*nSubSamples];
   checkAllocate(mSamplePts, "mSamplePts in RSMSobol2::analyze");

   pData *pPtr = NULL;
   if (ioPtr != NULL)
   {
      pPtr = ioPtr->getAuxData();
      pPtr->nDbles_ = nInputs * nInputs;
      pPtr->dbleArray_ = new double[nInputs * nInputs];
      for (ii = 0; ii < nInputs*nInputs; ii++) pPtr->dbleArray_[ii] = 0.0;
      pPtr->dbleData_ = variance;
   }

   printAsterisks(PL_INFO, 0);
   vces_ = new double[nInputs*nInputs];
   ecvs_ = new double[nInputs*nInputs];
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         vce = 0.0;
         if (ii2 <= ii) continue;
         if (selectedInput != -1 && ii != selectedInput && 
             ii2 != selectedInput)
            continue;
         printOutTS(PL_DETAIL, "RSMSobol2: processing input pair %d, %d\n",
                    ii+1, ii2+1);

         currNLevels = nLevels / 2;
         for (iR = 0; iR < 2; iR++)
         {
            printOutTS(PL_DETAIL,
                 "RSMSobol2: processing refinement %d (4)\n",iR+1);

            cLower[0] = xLower[ii];
            cUpper[0] = xUpper[ii];
            cLower[1] = xLower[ii2];
            cUpper[1] = xUpper[ii2];
            if (noPDF == 0)
            {
               corMat.setDim(2,2);
               corMat.setEntry(0, 0, corMatp->getEntry(ii,ii));
               corMat.setEntry(1, 1, corMatp->getEntry(ii2,ii2));
               pdfFlags1[0] = pdfFlags[ii];
               pdfFlags1[1] = pdfFlags[ii2];
               inputMeans1[0] = inputMeans[ii];
               inputMeans1[1] = inputMeans[ii2];
               inputStdevs1[0] = inputStdevs[ii];
               inputStdevs1[1] = inputStdevs[ii2];
               pdfman1 = new PDFManager();
               pdfman1->initialize(2,pdfFlags1,inputMeans1,
                                   inputStdevs1,corMat,NULL,NULL);
               vecLB.load(2, cLower);
               vecUB.load(2, cUpper);
               vecOut.setLength(currNLevels*2);
               pdfman1->genSample(currNLevels*currNLevels,vecOut,vecLB,vecUB);
               for (jj = 0; jj < currNLevels*currNLevels*2; jj++) 
                  samplePts2D[jj] = vecOut[jj];
               delete pdfman1;
            }
            else
            {
               sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
               sampler->setInputBounds(2, cLower, cUpper);
               sampler->setOutputParams(1);
               sampler->setSamplingParams(currNLevels*currNLevels, 1, 0);
               sampler->initialize(0);
               SS = new int[currNLevels*currNLevels];
               ZZ = new double[currNLevels*currNLevels];
               checkAllocate(ZZ, "ZZ in RSMSobol2::analyze");
               sampler->getSamples(currNLevels*currNLevels,2,1,
                                   samplePts2D,ZZ,SS);
               delete [] SS;
               delete [] ZZ;
               delete sampler;
            }

            if (noPDF == 0)
            {
               corMat.setDim(nInputs-2, nInputs-2);
               for (jj = 0; jj < ii; jj++)
               {
                  cLower[jj] = xLower[jj];
                  cUpper[jj] = xUpper[jj];
                  corMat.setEntry(jj, jj, corMatp->getEntry(jj,jj));
                  pdfFlags2[jj] = pdfFlags[jj];
                  inputMeans2[jj] = inputMeans[jj];
                  inputStdevs2[jj] = inputStdevs[jj];
               }
               for (jj = ii+1; jj < ii2; jj++)
               {
                  cLower[jj-1] = xLower[jj];
                  cUpper[jj-1] = xUpper[jj];
                  corMat.setEntry(jj-1, jj-1, corMatp->getEntry(jj,jj));
                  pdfFlags2[jj-1] = pdfFlags[jj];
                  inputMeans2[jj-1] = inputMeans[jj];
                  inputStdevs2[jj-1] = inputStdevs[jj];
               }
               for (jj = ii2+1; jj < nInputs; jj++)
               {
                  cLower[jj-2] = xLower[jj];
                  cUpper[jj-2] = xUpper[jj];
                  corMat.setEntry(jj-2, jj-2, corMatp->getEntry(jj,jj));
                  pdfFlags2[jj-2] = pdfFlags[jj];
                  inputMeans2[jj-2] = inputMeans[jj];
                  inputStdevs2[jj-2] = inputStdevs[jj];
               }
               pdfman2 = new PDFManager();
               pdfman2->initialize(nInputs-2,pdfFlags2,inputMeans2,
                                   inputStdevs2,corMat,NULL,NULL);
               vecLB.load(nInputs-2, cLower);
               vecUB.load(nInputs-2, cUpper);
               vecOut.setLength(nSubSamples*(nInputs-2));
               pdfman2->genSample(nSubSamples, vecOut, vecLB, vecUB);
               for (jj = 0; jj < nSubSamples*(nInputs-2); jj++)
                  XX[jj] = vecOut[jj];
               delete pdfman2;
            }
            else
            {
               for (jj = 0; jj < ii; jj++)
               {
                  cLower[jj] = xLower[jj];
                  cUpper[jj] = xUpper[jj];
               }
               for (jj = ii+1; jj < ii2; jj++)
               {
                  cLower[jj-1] = xLower[jj];
                  cUpper[jj-1] = xUpper[jj];
               }
               for (jj = ii2+1; jj < nInputs; jj++)
               {
                  cLower[jj-2] = xLower[jj];
                  cUpper[jj-2] = xUpper[jj];
               }
               sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
               sampler->setInputBounds(nInputs-2, cLower, cUpper);
               sampler->setOutputParams(1);
               sampler->setSamplingParams(nSubSamples, 1, 1);
               sampler->initialize(0);
               SS = new int[nSubSamples];
               ZZ = new double[nSubSamples];
               checkAllocate(ZZ, "ZZ(2) in RSMSobol2::analyze");
               sampler->getSamples(nSubSamples,nInputs-2,1,XX,ZZ,SS);
               delete [] SS;
               delete [] ZZ;
               delete sampler;
            }

            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               for (jj = 0; jj < nSubSamples; jj++)
               {
                  oneSamplePt = &(XX[jj*(nInputs-2)]);
                  for (kk = 0; kk < ii; kk++)
                     mSamplePts[jj*nInputs+kk] = oneSamplePt[kk];
                  for (kk = ii+1; kk < ii2; kk++)
                     mSamplePts[jj*nInputs+kk] = oneSamplePt[kk-1];
                  for (kk = ii2+1; kk < nInputs; kk++)
                     mSamplePts[jj*nInputs+kk] = oneSamplePt[kk-2];
                  mSamplePts[jj*nInputs+ii] = samplePts2D[iL*2];
                  mSamplePts[jj*nInputs+ii2] = samplePts2D[iL*2+1];
               }

               faPtr->evaluatePoint(nSubSamples,mSamplePts,YY);

               for (jj = 0; jj < nSubSamples; jj++)
               {
                  oneSamplePt = &(mSamplePts[jj*nInputs]);
                  status = 1;
                  if (constrPtr != NULL)
                     ddata = constrPtr->evaluate(oneSamplePt,YY[jj],status);
                  if (status == 0) YY[jj] = PSUADE_UNDEFINED;
               }

               means[iL] = 0.0;
               sCnt = 0;
               for (jj = 0; jj < nSubSamples; jj++)
               {
                  if (YY[jj] != PSUADE_UNDEFINED)
                  {
                     means[iL] += YY[jj];
                     sCnt++;
                  }
               }
               bins[iL] = sCnt;
               if (sCnt < 1 && printLevel >= 5)
                  printOutTS(PL_INFO, 
                       "RSMSobol2 WARNING: subsample size = 0.\n");
               if (sCnt < 1) means[iL] = PSUADE_UNDEFINED;
               else          means[iL] /= (double) sCnt;

               vars[iL] = 0.0;
               ddata = means[iL];
               for (jj = 0; jj < nSubSamples; jj++)
               {
                  if (YY[jj] != PSUADE_UNDEFINED)
                     vars[iL] += (YY[jj]-ddata)*(YY[jj]-ddata);
               }
               if (sCnt < 1) vars[iL] = PSUADE_UNDEFINED;
               else          vars[iL] /= (double) sCnt;

            }

            dmean = 0.0;
            totalCnt = 0;
            for (iL = 0; iL < currNLevels*currNLevels; iL++) 
               totalCnt += bins[iL];
            if (totalCnt == 0)
            {
               printOutTS(PL_ERROR,
                    "RSMSobol2 ERROR: empty constrained space.\n");
               printOutTS(PL_ERROR,
                    "          Either try larger sample size or\n");
               printOutTS(PL_ERROR,"          use looser constraints.\n");
               exit(1);
            }
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               if (means[iL] != PSUADE_UNDEFINED)
                  dmean += means[iL] * bins[iL] / totalCnt;
            }
            vce = 0.0;
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               if (means[iL] != PSUADE_UNDEFINED)
                  vce += (means[iL]-dmean) * (means[iL]-dmean) * 
                         bins[iL] / totalCnt;
            }

            ecv = 0.0;
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
               if (vars[iL] != PSUADE_UNDEFINED) 
                  ecv += vars[iL] * bins[iL] / totalCnt;

            if (printLevel > 1 || iR == 1)
               printOutTS(PL_INFO, 
                    "VCE(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                    ii+1, ii2+1, vce, vce/variance);
            if (printLevel > 2)
               printOutTS(PL_INFO, 
                    "ECV(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                    ii+1, ii2+1, ecv, ecv/variance);
            currNLevels *= 2;
         }

         //save vces & ecvs
         vces_[ii*nInputs+ii2] = vce;
         ecvs_[ii*nInputs+ii2] = ecv;
         vces_[ii2*nInputs+ii] = vce;
         ecvs_[ii2*nInputs+ii] = ecv;
         if (pPtr != NULL)
         {
            pPtr->dbleArray_[ii*nInputs+ii2] = vce;
            pPtr->dbleArray_[ii2*nInputs+ii] = vce;
         }
      }
   }
   printAsterisks(PL_INFO, 0);

   delete [] cLower;
   delete [] cUpper;
   delete [] XX;
   delete [] YY;
   delete [] Y;
   delete [] means;
   delete [] vars;
   delete [] bins;
   delete [] mSamplePts;
   delete [] samplePts2D;
   delete [] pdfFlags1;
   delete [] pdfFlags2;
   delete [] inputMeans1;
   delete [] inputMeans2;
   delete [] inputStdevs1;
   delete [] inputStdevs2;
   delete faPtr;
   if (ioPtr == NULL) delete corMatp;
   if (constrPtr != NULL) delete constrPtr;
   if (pdfNull == 1)
   {
      delete [] pdfFlags;
      delete [] inputMeans;
      delete [] inputStdevs;
   }
   return 0.0;
}

// ************************************************************************
// perform analysis (for problems with joint PDFs)
// ------------------------------------------------------------------------
double RSMSobol2Analyzer::analyze2(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ii, ii2, jj, iR, status;
   int    nSubSamples=100, iL, sCnt, nLevels=100, outputID;
   int    currNLevels, nSamp, printLevel, *bins, totalCnt, bin1, bin2;
   double *xLower, *xUpper, *X, *Y, *XX, *YY, *Y2, dmean, ecv, ddata;
   double *oneSamplePt, *means, *vars, vce, variance, width1, width2;
   char   pString[500];
   pData         pPDF;
   PsuadeData    *ioPtr;
   FuncApprox    *faPtr;
   RSConstraints *constrPtr;
   psVector      vecIn, vecOut, vecUB, vecLB;
   PDFManager    *pdfman;

   if (isScreenDumpModeOn())
   {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO,"RSMSobol2: since joint PDFs have been specified, a\n");
      printOutTS(PL_INFO,"   different interaction analysis will be performed.\n");
   }
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
   Y = new double[nSamples];
   checkAllocate(Y, "Y in RSMSobol2::analyze2");
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   if (ioPtr != NULL)
   {
      constrPtr = new RSConstraints();
      constrPtr->genConstraints(ioPtr);
   }
   else constrPtr = NULL;

   if (ioPtr == NULL)
   {
      if (rstype_ < 0) jj = 0;
      else             jj  = rstype_; 
      faPtr = genFA(jj, nInputs, 0, nSamples);
      faPtr->setBounds(xLower, xUpper);
      faPtr->setOutputLevel(0);
   }
   else faPtr = genFAInteractive(ioPtr, 0);
   if (faPtr == NULL)
   {
      printOutTS(PL_ERROR,
           "RSMSobol2 ERROR: cannot create response surface.\n");
      delete [] Y;
      delete constrPtr;
      return 1.0e12;
   }
   status = faPtr->initialize(X, Y);

   if (isScreenDumpModeOn() || psAnaExpertMode_ == 1)
   {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO,
                 "* RSMSobol2 generates a mesh of size K x K for every\n");
      printOutTS(PL_INFO,
                 "*   pair of inputs and then creates a sample of size\n");
      printOutTS(PL_INFO,
                 "*   M for each mesh point. The total sample size is:\n");
      printOutTS(PL_INFO,
                 "*       N = M * K * K * nInputs * (nInputs - 1) / 2.\n");
      printOutTS(PL_INFO,"* NOW, nInputs = %d\n", nInputs);
      printOutTS(PL_INFO,"* Please select your desired M and K.\n");
      printOutTS(PL_INFO,"* Recommendation: K x K >> M.\n");
      printOutTS(PL_INFO,"* NOTE: large M and K can take a long time.\n");
      printEquals(PL_INFO, 0);
      sprintf(pString,"Enter M (suggestion: 100 - 1000) : ");
      nSubSamples = getInt(100, 1000, pString);
      sprintf(pString, "Enter nLevels (suggestion: 50 - 500) : ");
      nLevels = getInt(50, 500, pString);
      printAsterisks(PL_INFO, 0);
   }
   else
   {
      nSubSamples = 100;
      nLevels = 100;
      if (isScreenDumpModeOn())
      {
         printOutTS(PL_INFO,"* RSMSobol2: default M = %d.\n", nSubSamples);
         printOutTS(PL_INFO,"* RSMSobol2: default K = %d.\n", nLevels);
         printOutTS(PL_INFO,
           "* To change these settings, re-run with ana_expert mode on.\n");
         printAsterisks(PL_INFO, 0);
      }
   }

   nSamp = nLevels * nLevels * nSubSamples;
   if (isScreenDumpModeOn() && printLevel > 1)
   {
      printOutTS(PL_INFO,
           "* RSMSobol2 INFO: creating a sample for basic statistics.\n");
      printOutTS(PL_INFO,"*                 sample size = %d\n", nSamp);
   }

   XX = new double[nSamp*nInputs];
   YY = new double[nSamp];
   checkAllocate(YY, "YY in RSMSobol2::analyze2");

   pdfman = NULL;
   if (ioPtr != NULL)
   {
      pdfman = new PDFManager();
      pdfman->initialize(ioPtr);
   }
   else
   {
      pdfman = new PDFManager();
      int *inputPDFs = adata.inputPDFs_;
      double *inputMeans = adata.inputMeans_;
      double *inputStdevs = adata.inputStdevs_;
      if (inputPDFs == NULL || inputMeans == NULL || inputStdevs == NULL)
      {
         printf("RSMSobol2 ERROR: PDF information not provided.\n");
         exit(1);
      }
      psMatrix cMat;
      cMat.setDim(nInputs, nInputs);
      for (ii = 0; ii < nInputs; ii++) cMat.setEntry(ii,ii,1);
      char **snames = new char*[nInputs];
      for (ii = 0; ii < nInputs; ii++) 
      {
         snames[ii] = new char[100];
         sprintf(snames[ii], "X%d", ii+1);
      }
      pdfman->initialize(nInputs, inputPDFs, inputMeans, inputStdevs,
                         cMat, snames, NULL);
      for (ii = 0; ii < nInputs; ii++) delete [] snames[ii];
      delete [] snames;
   }
   vecLB.load(nInputs, xLower);
   vecUB.load(nInputs, xUpper);
   vecOut.setLength(nSamp*nInputs);
   pdfman->genSample(nSamp, vecOut, vecLB, vecUB);
   for (ii = 0; ii < nSamp*nInputs; ii++) XX[ii] = vecOut[ii];

   if (isScreenDumpModeOn())
      printOutTS(PL_INFO, 
        "RSMSobol2: running the sample with response surface...\n");
   faPtr->evaluatePoint(nSamp, XX, YY);
   if (isScreenDumpModeOn())
      printOutTS(PL_INFO, 
        "RSMSobol2: done running the sample with response surface.\n");
   
   if (constrPtr != NULL)
   {
      for (ii = 0; ii < nSamp; ii++)
      {
         oneSamplePt = &(XX[ii*nInputs]);
         ddata = constrPtr->evaluate(oneSamplePt,YY[ii],status);
         if (status == 0) YY[ii] = PSUADE_UNDEFINED;
      }
   }
   
   dmean = 0.0;
   sCnt = 0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
      {
         dmean += YY[ii];
         sCnt++;
      }
   }
   if (sCnt > 1) dmean /= (double) sCnt;
   else
   {
      printOutTS(PL_ERROR,"RSMSobol2 ERROR: too few samples that satisify\n");
      printOutTS(PL_ERROR,"constraints (%d out of %d).\n",sCnt,nSamp);
      delete [] XX;
      delete [] YY;
      delete [] Y;
      delete faPtr;
      delete pdfman;
      delete constrPtr;
      return PSUADE_UNDEFINED;
   }
   variance = 0.0;
   for (ii = 0; ii < nSamp; ii++)
   {
      if (YY[ii] != PSUADE_UNDEFINED)
         variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
   }
   variance /= (double) sCnt;
   if (isScreenDumpModeOn())
   {
      printOutTS(PL_INFO,
        "* RSMSobol2: sample mean    (based on %d points) = %e\n",
        sCnt, dmean);
      printOutTS(PL_INFO,
        "* RSMSobol2: total std dev  (based on %d points) = %e\n",
        sCnt, sqrt(variance));
      printAsterisks(PL_INFO, 0);
   }
   if (variance == 0.0) variance = 1.0;
   delete pdfman;

   //save mean & std
   outputMean_ = dmean;
   outputStd_ = sqrt(variance);

   pData *pPtr = NULL;
   if (ioPtr != NULL)
   {
      pPtr = ioPtr->getAuxData();
      pPtr->nDbles_ = nInputs * nInputs;
      pPtr->dbleArray_ = new double[nInputs * nInputs];
      for (ii = 0; ii < nInputs*nInputs; ii++) pPtr->dbleArray_[ii] = 0.0;
      pPtr->dbleData_ = variance;
   }

   vces_ = new double[nInputs*nInputs];
   ecvs_ = new double[nInputs*nInputs];
   means = new double[nLevels*nLevels];
   vars  = new double[nLevels*nLevels];
   bins  = new int[nLevels*nLevels];
   checkAllocate(bins, "bins in RSMSobol2::analyze2");

   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         if (isScreenDumpModeOn())
            printOutTS(PL_DETAIL, "RSMSobol2: processing input pair %d, %d\n",
                  ii+1, ii2+1);

         currNLevels = nLevels / 2;
         for (iR = 0; iR < 2; iR++)
         {
            if (isScreenDumpModeOn())
               printOutTS(PL_DETAIL, "RSMSobol2: processing refinement %d\n",iR);

            width1 = (xUpper[ii] - xLower[ii]) / currNLevels;
            width2 = (xUpper[ii2] - xLower[ii2]) / currNLevels;
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               means[iL] = 0.0;
               vars[iL] = 0.0;
               bins[iL] = 0;
            }
            for (jj = 0; jj < nSamp; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
               {
                  ddata = XX[jj*nInputs+ii];
                  bin1 = (int) ((ddata - xLower[ii]) / width1);
                  ddata = XX[jj*nInputs+ii2];
                  bin2 = (int) ((ddata - xLower[ii2]) / width2);
                  if (bin1 == currNLevels) bin1 = currNLevels - 1;
                  if (bin2 == currNLevels) bin2 = currNLevels - 1;
                  means[bin1*currNLevels+bin2] += YY[jj];
                  bins[bin1*currNLevels+bin2]++;
               }
            }
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               sCnt = bins[iL];
               if (sCnt < 1 && printLevel >= 5)
                  printOutTS(PL_DUMP,
                       "RSMSobol2 WARNING: subsample size = 0.\n");
               if (sCnt < 1) means[iL] = PSUADE_UNDEFINED;
               else          means[iL] /= (double) sCnt;
            }
            for (jj = 0; jj < nSamp; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
               {
                  ddata = XX[jj*nInputs+ii];
                  bin1 = (int) ((ddata - xLower[ii]) / width1);
                  ddata = XX[jj*nInputs+ii2];
                  bin2 = (int) ((ddata - xLower[ii2]) / width2);
                  if (bin1 == currNLevels) bin1 = currNLevels - 1;
                  if (bin2 == currNLevels) bin2 = currNLevels - 1;
                  ddata = means[bin1*currNLevels+bin2];
                  vars[bin1*currNLevels+bin2] += 
                         (YY[jj]-ddata)*(YY[jj]-ddata);
               }
            }
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               sCnt = bins[iL];
               if (sCnt < 1) vars[iL] = PSUADE_UNDEFINED;
               else          vars[iL] /= (double) sCnt;
            }

            totalCnt = 0;
            for (iL = 0; iL < currNLevels*currNLevels; iL++) 
               totalCnt += bins[iL];
            if (totalCnt == 0)
            {
               printOutTS(PL_ERROR, 
                    "RSMSobol2 ERROR: empty constrained space.\n");
               printOutTS(PL_ERROR, 
                    "          Either try larger sample size or\n");
               printOutTS(PL_ERROR, "          use looser constraints.\n");
               exit(1);
            }

            dmean = 0.0;
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               if (means[iL] != PSUADE_UNDEFINED)
                  dmean += means[iL] * bins[iL] / totalCnt;
            }
            vce = 0.0;
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
               if (means[iL] != PSUADE_UNDEFINED)
                  vce += (means[iL]-dmean) * (means[iL]-dmean) * 
                         bins[iL] / totalCnt;

            ecv = 0.0;
            for (iL = 0; iL < currNLevels*currNLevels; iL++)
            {
               if (vars[iL] != PSUADE_UNDEFINED)
                  ecv += vars[iL] * bins[iL] / totalCnt;
            }
            if (isScreenDumpModeOn() && (printLevel > 2 || iR == 1))
            {
               printOutTS(PL_INFO, 
                      "VCE(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                      ii+1, ii2+1, vce, vce/variance);
            }
            if (isScreenDumpModeOn() && (printLevel > 3 || iR == 1))
            {
               printOutTS(PL_DETAIL, 
                      "ECV(%3d,%3d) = %12.4e, (normalized) = %10.3e\n",
                      ii+1, ii2+1, ecv, ecv/variance);
            }
            currNLevels *= 2;
         }
         //save vces & ecvs
         vces_[ii*nInputs+ii2] = vces_[ii2*nInputs+ii] = vce;
         ecvs_[ii*nInputs+ii2] = ecvs_[ii2*nInputs+ii] = ecv;
         if (pPtr != NULL)
         {
            pPtr->dbleArray_[ii*nInputs+ii2] = vce;
            pPtr->dbleArray_[ii2*nInputs+ii] = vce;
         }
      }
   }
   if (isScreenDumpModeOn()) printAsterisks(PL_INFO, 0);
    
   delete constrPtr;
   delete faPtr;
   delete [] XX;
   delete [] YY;
   delete [] Y;
   delete [] means;
   delete [] vars;
   delete [] bins;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobol2Analyzer& RSMSobol2Analyzer::operator=(const RSMSobol2Analyzer &)
{
   printOutTS(PL_ERROR,"RSMSobol2 operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int RSMSobol2Analyzer::get_nInputs()
{
   return nInputs_;
}

double RSMSobol2Analyzer::get_outputMean()
{
   return outputMean_;
}

double RSMSobol2Analyzer::get_outputStd()
{
   return outputStd_;
}

double RSMSobol2Analyzer::get_vce(int ind1, int ind2)
{
   if (ind1 < 0 || ind1 >= nInputs_)
   {
      printf("RSMSobol2 ERROR: get_vce index 1 error.\n");
      return 0.0;
   }
   if (ind2 < 0 || ind2 >= nInputs_)
   {
      printf("RSMSobol2 ERROR: get_vce index 2 error.\n");
      return 0.0;
   }
   if (vces_ == NULL)
   {
      printf("RSMSobol2 ERROR: get_vce has no value.\n");
      return 0;
   }
   return vces_[ind1*nInputs_+ind2]/outputStd_/outputStd_;
}

double RSMSobol2Analyzer::get_ecv(int ind1, int ind2)
{
   if (ind1 < 0 || ind1 >= nInputs_)
   {
      printf("RSMSobol2 ERROR: get_ecv index 1 error.\n");
      return 0.0;
   }
   if (ind2 < 0 || ind2 >= nInputs_)
   {
      printf("RSMSobol2 ERROR: get_ecv index 2 error.\n");
      return 0.0;
   }
   if (ecvs_ == NULL)
   {
      printf("RSMSobol2 ERROR: get_ecv has no value.\n");
      return 0;
   }
   return ecvs_[ind1*nInputs_+ind2]/outputStd_/outputStd_;
}


