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
// Functions for the class RSMSobolGAnalyzer  
// Sobol' group main effect analysis (with response surface)
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Vector.h"
#include "Matrix.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "RSConstraints.h"
#include "Sampling.h"
#include "PDFManager.h"
#include "PsuadeData.h"
#include "pData.h"
#include "RSMSobolGAnalyzer.h"
#include "sysdef.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSMSobolGAnalyzer::RSMSobolGAnalyzer() : Analyzer()
{
   setName("RSMSOBOLG");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSMSobolGAnalyzer::~RSMSobolGAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double RSMSobolGAnalyzer::analyze(aData &adata)
{
   int        nInputs, nOutputs, nSamples, ii, ii2, jj, kk;
   int        nGroups, groupID, **groupMembers, length, status, sCnt;
   int        ir, index, currNSamples, nSubSamplesG, nInputsG, nInputsN;
   int        nSubSamplesN, nSamp, *pdfFlags, outputID, noPDF=1, *SS;
   int        printLevel, *bins, totalCnt, *iArray, *pdfFlagsG, *pdfFlagsN;
   size_t     compFlag;
   double     *xLower, *xUpper, *X, *Y, *Y2, *XX, *YY, *XXG, *XXN, ddata;
   double     *oneSamplePt, variance, *cLower, *cUpper, *means, *vars;
   double     *inputMeans, *inputStdevs, vce, dmean, ecv, *mSamplePts;
   double     *inputMeansG, *inputMeansN, *inputStdevsG, *inputStdevsN;
   PsuadeData *ioPtr;
   FILE       *fp=NULL;
   char       cfname[1001], pString[1001], *cString, winput1[1001];
   char       winput2[1001], lineIn[1001];
   RSConstraints *constrPtr;
   FuncApprox    *faPtr;
   psVector      vecIn, vecOut, vecUB, vecLB;
   pData         pCorMat;
   psMatrix      *corMatp, corMat;
   PDFManager    *pdfman, *pdfmanN, *pdfmanG;
   Sampling      *sampler;

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,
        "*          RS-based Group First Order Sobol' Indices \n");
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO,"* TO GAIN ACCESS TO DIFFERENT OPTIONS: SET\n");
   printOutTS(PL_INFO,"*\n");
   printOutTS(PL_INFO,
        "* - ana_expert mode to finetune RSMSobolG parameters, \n");
   printOutTS(PL_INFO,
        "*   (e.g. sample size for integration can be adjusted).\n");
   printOutTS(PL_INFO,
        "* - rs_expert to finetune response surface for RSMSobolG,\n");
   printOutTS(PL_INFO,
        "* - printlevel to 1 or higher to display more information.\n");
   printEquals(PL_INFO, 0);

   printLevel = adata.printLevel_;
   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   X        = adata.sampleInputs_;
   Y2       = adata.sampleOutputs_;
   outputID = adata.outputID_;
   ioPtr    = adata.ioPtr_;
   pdfFlags    = adata.inputPDFs_;
   inputMeans  = adata.inputMeans_;
   inputStdevs = adata.inputStdevs_;
   if (pdfFlags != NULL)
   {
      for (ii = 0; ii < nInputs; ii++)
         if (pdfFlags[ii] != 0) noPDF = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         if (pdfFlags[ii] == PSUADE_PDF_SAMPLE)
         {
            printOutTS(PL_ERROR,
                 "* RSMSobolG ERROR: S PDF type currently not supported.\n");
            return PSUADE_UNDEFINED;
         }
      }
   }
   if (noPDF == 1) 
        printOutTS(PL_INFO,"RSMSobolG INFO: all uniform distributions.\n");
   else
   {
      printOutTS(PL_INFO,
           "RSMSobolG INFO: non-uniform distributions detected,\n");
      printOutTS(PL_INFO,
           "                which will be used in this analysis.\n");
   }

   if (nInputs <= 1 || nSamples <= 0 || nOutputs <= 0)
   {
      printOutTS(PL_ERROR, "RSMSobolG ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nInputs <= 2)
   {
      printOutTS(PL_ERROR, "RSMSobolG INFO: nInputs == 2.\n");
      printOutTS(PL_ERROR, 
           "   You do not need to perform this when nInputs = 2.\n");
      return PSUADE_UNDEFINED;
   }
   if (outputID >= nOutputs || outputID < 0)
   {
      printOutTS(PL_ERROR,
           "RSMSobolG ERROR: invalid output ID (%d).\n",outputID);
      return PSUADE_UNDEFINED;
   }
   if (ioPtr == NULL)
   {
      printOutTS(PL_ERROR, "RSMSobolG ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   } 
   status = 0;
   for (ii = 0; ii < nSamples; ii++)
      if (Y2[nOutputs*ii+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printOutTS(PL_ERROR,
           "RSMSobolG ERROR: Some outputs are undefined. Prune\n");
      printOutTS(PL_ERROR,
           "                 the undefined sample points first.\n");
      return PSUADE_UNDEFINED;
   }
   Y = new double[nSamples];
   checkAllocate(Y, "Y in RSMSobolG::analyze");
   for (ii = 0; ii < nSamples; ii++) Y[ii] = Y2[ii*nOutputs+outputID];

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"To use this function, you need to provide a file\n");
   printOutTS(PL_INFO,"specifying group information, in the form of : \n");
   printOutTS(PL_INFO,"line 1: PSUADE_BEGIN\n");
   printOutTS(PL_INFO,"line 2: <d> specifying the number of groups\n");
   printOutTS(PL_INFO,
           "line 3 to line <d>+2: group number, size, input numbers\n");
   printOutTS(PL_INFO,"last line: PSUADE_END\n");
   while (1)
   {
      printOutTS(PL_INFO,"Enter the group file : ");
      scanf("%s", cfname);
      fp = fopen(cfname, "r");
      if (fp != NULL)  
         break;
      else 
         printOutTS(PL_ERROR,
              "ERROR : file not found (or file name too long).\n");
   }
   fp = fopen(cfname, "r");
   if (fp != NULL)
   {
      fgets(lineIn, 1000, fp);
      sscanf(lineIn, "%s", pString);
      if (!strcmp(pString, "PSUADE_BEGIN"))
      {
         fscanf(fp, "%d", &nGroups);
         if (nGroups <= 0)
         {
            printOutTS(PL_ERROR, "RSMSobolG ERROR: nGroups <= 0.\n");
            fclose(fp);
            exit(1);
         }
         groupMembers = new int*[nGroups];
         checkAllocate(groupMembers,"groupMembers in RSMSobolG::analyze");
         for (ii = 0; ii < nGroups; ii++)
         {
            fscanf(fp, "%d", &groupID);
            if (groupID != ii+1)
            {
               printOutTS(PL_ERROR,
                    "RSMSobolG ERROR: invalid groupID %d",groupID);
               printOutTS(PL_ERROR," should be %d\n", ii+1);
               fclose(fp);
               exit(1);
            }
            fscanf(fp, "%d", &length);
            if (length <= 0 || length >= nInputs)
            {
               printOutTS(PL_ERROR, 
                    "RSMSobolG ERROR: invalid group length.\n");
               fclose(fp);
               exit(1);
            }
            groupMembers[ii] = new int[nInputs];
            checkAllocate(groupMembers[ii],"gMembers in RSMSobolG::analyze");
            for (jj = 0; jj < nInputs; jj++) groupMembers[ii][jj] = 0;
            sCnt = 1;
            for (jj = 0; jj < length; jj++)
            {
               fscanf(fp, "%d", &index);
               if (index <= 0 || index > nInputs)
               {
                  printOutTS(PL_ERROR, 
                    "RSMSobolG ERROR: invalid group member.\n");
                  fclose(fp);
                  exit(1);
               }
               groupMembers[ii][index-1] = sCnt++;
            }
         }
         fgets(lineIn, 1000, fp);
         fgets(lineIn, 1000, fp);
         sscanf(lineIn, "%s", pString);
         if (strcmp(pString, "PSUADE_END"))
         {
            printOutTS(PL_ERROR, "RSMSobolG ERROR: PSUADE_END not found.\n");
            fclose(fp);
            exit(1);
         }
      }
      else
      {
         printOutTS(PL_ERROR, "RSMSobolG ERROR: PSUADE_BEGIN not found.\n");
         fclose(fp);
         exit(1);
      }
      fclose(fp);
   }

   ioPtr->getParameter("input_cor_matrix", pCorMat);
   corMatp = (psMatrix *) pCorMat.psObject_;
   for (ii = 0; ii < nGroups; ii++)
   {
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         for (jj = 0; jj < nInputs; jj++)
         {
            if (((groupMembers[ii][ii2] != 0 && groupMembers[ii][jj] == 0) ||
                 (groupMembers[ii][ii2] == 0 && groupMembers[ii][jj] != 0)) &&
                corMatp->getEntry(ii2,jj) != 0.0)
            {
               printOutTS(PL_ERROR,"RSMSobolG INFO: currently cannot handle\n");
               printOutTS(PL_ERROR,"          correlated inputs (joint PDF)\n");
               printOutTS(PL_ERROR,"          across different groups.\n");
               for (ir = 0; ir < nGroups; ir++) delete [] groupMembers[ir];
               delete [] groupMembers;
               return PSUADE_UNDEFINED;
            }
         }
      }
   }

   constrPtr = new RSConstraints();
   constrPtr->genConstraints(ioPtr);

   faPtr = genFAInteractive(ioPtr, 0);
   status = faPtr->initialize(X, Y);

   printAsterisks(PL_INFO, 0);
   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO,
           "* RSMSobolG creates a sample of size M1 for each subgroup\n");
      printOutTS(PL_INFO,
           "* of inputs and M2 for the other inputs when computing\n");
      printOutTS(PL_INFO,
           "* group sensitivity indices. The total sample size is thus:\n");
      printOutTS(PL_INFO,"*     N = M1 * M2 * nGroups.\n");
      nSubSamplesG = 50000;
      nSubSamplesN = 100;
      printOutTS(PL_INFO,"  default M1 = %d.\n", nSubSamplesG);
      printOutTS(PL_INFO,"  default M2 = %d.\n", nSubSamplesN);
      printOutTS(PL_INFO,"* Please select your desired M1 and M2.\n");
      printOutTS(PL_INFO,"* Recommendation: M1 > M2.\n");
      printOutTS(PL_INFO,"* Note: large M1 and M2 can take a long time.\n");
      printEquals(PL_INFO, 0);
      sprintf(pString,
             "Enter M1 (suggestion: 10000 - 100000, default = 50000) : ");
      nSubSamplesG = getInt(10000, 200000, pString);
      sprintf(pString, "Enter M2 (suggestion: 100 - 500, default = 100) : ");
      nSubSamplesN = getInt(100, 1000, pString);
      printAsterisks(PL_INFO, 0);
   }
   else
   {
      nSubSamplesG = 50000;
      nSubSamplesN = 100;
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("RSMSobolG_nsubsamples_ingroup");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamplesG);
            if (nSubSamplesG < 10000)
            {
               printOutTS(PL_INFO,
                    "RSMSobolG INFO: nSubSamplesG should be >= 10000.\n");
               nSubSamplesG = 10000;
            }
            else
            {
               printOutTS(PL_INFO,
                    "RSMSobolG INFO: nSubSamplesG = %d (config).\n",
                    nSubSamplesG);
            }
         }
         cString = psConfig_->getParameter("RSMSobolG_nsubsamples_outgroup");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %d", winput1, winput2, &nSubSamplesN);
            if (nSubSamplesN < 100)
            {
               printOutTS(PL_INFO,
                    "RSMSobolG INFO: nSubSamplesN should be >= 100.\n");
               nSubSamplesN = 100;
            }
            else
            {
               printOutTS(PL_INFO,
                    "RSMSobolG INFO: nSubSamplesN = %d (config).\n",
                    nSubSamplesN);
            }
         }
      }
      printOutTS(PL_INFO,"RSMSobolG: default M1 = %d.\n", nSubSamplesG);
      printOutTS(PL_INFO,"RSMSobolG: default M2 = %d.\n", nSubSamplesN);
      printOutTS(PL_INFO,
           "To change these settings, re-run with ana_expert mode on.\n");
   }
   printEquals(PL_INFO, 0);

   nSamp = 25000;
   printOutTS(PL_INFO,
        "RSMSobolG INFO: creating a sample for basic statistics.\n");
   printOutTS(PL_INFO,"                sample size = %d\n", nSamp);

   XX = new double[nSamp*8*nInputs];
   YY = new double[nSamp*8];
   checkAllocate(YY,"YY in RSMSobolG::analyze");

   for (ir = 0; ir < 3; ir++)
   {
      printOutTS(PL_INFO, 
           "RSMSobolG: compute statistics,sample size = %d\n",nSamp);
       
      if (noPDF == 0)
      {
         printOutTS(PL_INFO,"RSMSobolG INFO: non-uniform PDFs detected. \n");
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
         if (nInputs > 51)
              sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampler->setInputBounds(nInputs, xLower, xUpper);
         sampler->setOutputParams(1);
         sampler->setSamplingParams(nSamp, 1, 1);
         sampler->initialize(0);
         SS = new int[nSamp];
         checkAllocate(SS,"SS in RSMSobolG::analyze");
         sampler->getSamples(nSamp, nInputs, 1, XX, YY, SS);
         delete [] SS;
         delete sampler;
      }

      if (ir < 3 && printLevel > 1)
         printOutTS(PL_INFO,
              "RSMSobolG: running the sample with response surface...\n");
      faPtr->evaluatePoint(nSamp, XX, YY);
      if (ir < 3 && printLevel > 1)
         printOutTS(PL_INFO,
              "RSMSobolG: done running the sample with response surface.\n");

      for (ii = 0; ii < nSamp; ii++)
      {
         oneSamplePt = &(XX[ii*nInputs]);
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
              "RSMSobolG ERROR: too few samples that satisify the\n");
         printOutTS(PL_ERROR,"constraints (%d out of %d)\n",sCnt,nSamp);
         delete [] XX;
         delete [] YY;
         delete faPtr;
         return PSUADE_UNDEFINED;
      }
      variance = 0.0;
      for (ii = 0; ii < nSamp; ii++)
      {
         if (YY[ii] != PSUADE_UNDEFINED)
            variance += (YY[ii] - dmean) * (YY[ii] - dmean) ;
      }
      variance /= (double) sCnt;
      if (printLevel > 3 || ir == 2)
      {
         printOutTS(PL_INFO,
              "RSMSobolG: sample mean    (based on N = %d) = %10.3e\n",
              sCnt, dmean);
         printOutTS(PL_INFO,
              "RSMSobolG: sample std dev (based on N = %d) = %10.3e\n",
              sCnt, sqrt(variance));
      }
      nSamp *= 2;
   }
   if (variance == 0.0) variance = 1.0;
   delete [] XX;
   delete [] YY;

   cLower = new double[nInputs];
   cUpper = new double[nInputs];
   XXG    = new double[nSubSamplesG*nInputs];
   XXN    = new double[nSubSamplesN*nInputs];
   YY     = new double[nSubSamplesN+nSubSamplesG];
   means  = new double[nSubSamplesG];
   vars   = new double[nSubSamplesG];
   bins   = new int[nSubSamplesG];
   iArray = new int[nInputs];
   mSamplePts   = new double[nInputs*nSubSamplesN];
   inputMeansG  = new double[nInputs];
   inputMeansN  = new double[nInputs];
   inputStdevsG = new double[nInputs];
   inputStdevsN = new double[nInputs];
   pdfFlagsG    = new int[nInputs];
   pdfFlagsN    = new int[nInputs];
   checkAllocate(pdfFlagsN,"pdfFlagsN in RSMSobolG::analyze");

   printAsterisks(PL_INFO, 0);
   for (ii = 0; ii < nGroups; ii++)
   {
      if (printLevel > 1)
      {
         printOutTS(PL_INFO, "RSMSobolG: processing group %d\n", ii+1);
         printOutTS(PL_INFO, "           group members: ");
         for (jj = 0; jj < nInputs; jj++)
         {
            if (groupMembers[ii][jj] != 0)
               printOutTS(PL_INFO, "%d ", jj+1);
         }
         printOutTS(PL_INFO, "\n");
      }
      nInputsN = 0;
      for (jj = 0; jj < nInputs; jj++)
      {
         if (groupMembers[ii][jj] == 0)
         {
            cLower[nInputsN] = xLower[jj];
            cUpper[nInputsN] = xUpper[jj];
            iArray[nInputsN] = jj;
            nInputsN++;
         }
      }
      if (noPDF == 0) 
      {
         nInputsN = 0;
         for (jj = 0; jj < nInputs; jj++)
         {
            if (groupMembers[ii][jj] == 0)
            {
               pdfFlagsN[nInputsN] = pdfFlags[jj];
               inputMeansN[nInputsN] = inputMeans[jj];
               inputStdevsN[nInputsN] = inputStdevs[jj];
               nInputsN++;
            }
         }
         corMat.setDim(nInputsN, nInputsN);
         for (jj = 0; jj < nInputsN; jj++)
            for (ii2 = 0; ii2 < nInputsN; ii2++)
               corMat.setEntry(jj,ii2,corMatp->getEntry(iArray[jj],
                                                        iArray[ii2])); 
         pdfmanN = new PDFManager();
         pdfmanN->initialize(nInputsN, pdfFlagsN, inputMeansN,
                             inputStdevsN, corMat,NULL,NULL);
         vecLB.load(nInputsN, cLower);
         vecUB.load(nInputsN, cUpper);
         vecOut.setLength(nSubSamplesN*nInputsN);
         pdfmanN->genSample(nSubSamplesN, vecOut, vecLB, vecUB);
         for (jj = 0; jj < nSubSamplesN*nInputsN; jj++)
            XXN[jj] = vecOut[jj];
         delete pdfmanN;
      }
      else
      {
         if (nInputsN > 51)
              sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
         else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampler->setInputBounds(nInputsN, cLower, cUpper);
         sampler->setOutputParams(1);
         sampler->setSamplingParams(nSubSamplesN, 1, 1);
         sampler->initialize(0);
         SS = new int[nSubSamplesN];
         checkAllocate(SS,"SS(2) in RSMSobolG::analyze");
         sampler->getSamples(nSubSamplesN,nInputsN,1,XXN,YY,SS);
         delete [] SS;
         delete sampler;
      }

      currNSamples = nSubSamplesG / 8;
      nInputsG = 0;
      for (jj = 0; jj < nInputs; jj++)
      {
         if (groupMembers[ii][jj] != 0)
         {
            cLower[nInputsG] = xLower[jj];
            cUpper[nInputsG] = xUpper[jj];
            iArray[nInputsG] = jj;
            nInputsG++;
         }
      }
      if (noPDF == 0) 
      {
         nInputsG = 0;
         for (jj = 0; jj < nInputs; jj++)
         {
            if (groupMembers[ii][jj] != 0)
            {
               pdfFlagsG[nInputsG] = pdfFlags[jj];
               inputMeansG[nInputsG] = inputMeans[jj];
               inputStdevsG[nInputsG] = inputStdevs[jj];
               nInputsG++;
            }
         }
         corMat.setDim(nInputsG, nInputsG);
         for (jj = 0; jj < nInputsG; jj++)
            for (ii2 = 0; ii2 < nInputsG; ii2++)
               corMat.setEntry(jj,ii2,
                     corMatp->getEntry(iArray[jj],iArray[ii2])); 
         vecLB.load(nInputsG, cLower);
         vecUB.load(nInputsG, cUpper);
         vecOut.setLength(nSubSamplesG*nInputsG);
      }

      for (ir = 0; ir < 4; ir++)
      {

         printOutTS(PL_DETAIL,"RSMSobolG: processing refinement %d\n",ir+1);

         printOutTS(PL_DETAIL,"nSamplesG = %d, nSamplesN = %d\n",currNSamples,
                    nSubSamplesN);
         if (noPDF == 0) 
         {
            pdfmanG = new PDFManager();
            pdfmanG->initialize(nInputsG, pdfFlagsG, inputMeansG,
                                inputStdevsG, corMat,NULL,NULL);
            pdfmanG->genSample(nSubSamplesG, vecOut, vecLB, vecUB);
            for (jj = 0; jj < nSubSamplesG*nInputsG; jj++)
               XXG[jj] = vecOut[jj];
            delete pdfmanG;
         }
         else
         {
            if (nInputsN > 51)
                 sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
            else sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
            sampler->setInputBounds(nInputsG, cLower, cUpper);
            sampler->setOutputParams(1);
            sampler->setSamplingParams(nSubSamplesG, 1, 1);
            sampler->initialize(0);
            SS = new int[nSubSamplesG];
            checkAllocate(SS,"SS(3) in RSMSobolG::analyze");
            sampler->getSamples(nSubSamplesG,nInputsG,1,XXG,YY,SS);
            delete [] SS;
            delete sampler;
         }

         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
         {
            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               sCnt = 0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (groupMembers[ii][kk] != 0)
                  {
                     mSamplePts[jj*nInputs+kk] = XXG[ii2*nInputsG+sCnt];
                     sCnt++;
                  }
               }
            }

            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               sCnt = 0;
               for (kk = 0; kk < nInputs; kk++)
               {
                  if (groupMembers[ii][kk] == 0)
                  {
                     mSamplePts[jj*nInputs+kk] = XXN[jj*nInputsN+sCnt];
                     sCnt++;
                  }
               }
            }

            faPtr->evaluatePoint(nSubSamplesN, mSamplePts, YY);

            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               oneSamplePt = &mSamplePts[jj*nInputs];
               ddata = constrPtr->evaluate(oneSamplePt,YY[jj],status);
               if (status == 0) YY[jj] = PSUADE_UNDEFINED;
            }

            means[ii2] = 0.0;
            sCnt = 0;
            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
               {
                  means[ii2] += YY[jj];
                  sCnt++;
               }
            }
            bins[ii2] = sCnt;
            if (sCnt < 1 && printLevel >= 5)
               printOutTS(PL_DUMP, "RSMSobolG WARNING: subsample size = 0.\n");
            if (sCnt < 1) means[ii2] = PSUADE_UNDEFINED;
            else          means[ii2] /= (double) sCnt;

            vars[ii2] = 0.0;
            ddata = means[ii2];
            for (jj = 0; jj < nSubSamplesN; jj++)
            {
               if (YY[jj] != PSUADE_UNDEFINED)
                  vars[ii2] += (YY[jj] - ddata) * (YY[jj] - ddata);
            }
            if (sCnt < 1) vars[ii2] = PSUADE_UNDEFINED;
            else          vars[ii2] /= (double) sCnt;

            printOutTS(PL_DUMP, "RSMSobolG: Group %d\n", ii+1);
            printOutTS(PL_DUMP, "  refinement = %d, size = %d (%d),", ir, sCnt,
                      nSubSamplesN);
            printOutTS(PL_DUMP, 
                 " mean = %12.4e, var = %12.4e\n", means[ii2], vars[ii2]);
         }

         totalCnt = 0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++) totalCnt += bins[ii2];
         if (totalCnt == 0)
         {
            printOutTS(PL_ERROR, "RSMSobolG ERROR: empty constrained space.\n");
            exit(1);
         }
 
         dmean = 0.0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
         {
            if (means[ii2] != PSUADE_UNDEFINED)
               dmean += means[ii2] * bins[ii2] / totalCnt;
         }

         vce = 0.0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
            if (means[ii2] != PSUADE_UNDEFINED)
               vce += (means[ii2] - dmean) * (means[ii2] - dmean) *
                      bins[ii2] / totalCnt;

         ecv = 0.0;
         for (ii2 = 0; ii2 < nSubSamplesG; ii2++)
         {
            if (vars[ii2] != PSUADE_UNDEFINED)
               ecv += vars[ii2] * bins[ii2] / totalCnt;
         }

         if (ir == 3)
         {
            printOutTS(PL_INFO, "Unnormalized VCE (refinement=%3d) ", ir);
            printOutTS(PL_INFO, "for input group %3d = %12.4e\n", ii+1, vce);
         }
         printOutTS(PL_DETAIL, "Unnormalized ECV (refinement=%3d) ", ir);
         printOutTS(PL_DETAIL, "for input group %3d = %12.4e\n", ii+1, ecv);
         if (ir == 3)
            printOutTS(PL_INFO,
                 "** Normalized VCE for input group %3d = %12.4e\n",
                 ii+1, vce/variance);
         currNSamples *= 2;
      }
   }
   printAsterisks(PL_INFO, 0);
    
   delete constrPtr;
   delete faPtr;
   delete [] cLower;
   delete [] cUpper;
   delete [] XXG;
   delete [] XXN;
   delete [] inputMeansG;
   delete [] inputMeansN;
   delete [] inputStdevsG;
   delete [] inputStdevsN;
   delete [] pdfFlagsG;
   delete [] pdfFlagsN;
   delete [] mSamplePts;
   delete [] YY;
   delete [] Y;
   delete [] means;
   delete [] vars;
   delete [] bins;
   delete [] iArray;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSMSobolGAnalyzer& RSMSobolGAnalyzer::operator=(const RSMSobolGAnalyzer &)
{
   printOutTS(PL_ERROR, "RSMSobolG operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

