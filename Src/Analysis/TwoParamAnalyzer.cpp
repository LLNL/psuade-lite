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
// Functions for the class TwoParamAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : updated in 2006
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "TwoParamAnalyzer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Matrix.h"
#include "PsuadeData.h"
#include "pData.h"
#include "RSConstraints.h"
#include "PrintingTS.h"
using namespace std;

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
TwoParamAnalyzer::TwoParamAnalyzer(): Analyzer(),nInputs_(0),outputMean_(0),
                  outputStd_(0), sensitivity_(0), compareSensitivity_(0)
{
   setName("TwoParam");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
TwoParamAnalyzer::~TwoParamAnalyzer()
{
   if (sensitivity_) delete[] sensitivity_;
}

// ************************************************************************
// perform VCE analysis for two parameters
// ------------------------------------------------------------------------
double TwoParamAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, nSubSamples;
   int    ii, jj, ss, repID, subID, nReps1, ncount, index, status;
   int    symIndex, whichOutput, ii2, ss2, nReps, ir, *bins;
   int    errflag, nSubs, printLevel, iZero=0, iOne=1;
   int    totalCnt, validnReps, *pdfFlags;
   double *vce, *meanVCEVar, *vceMean, *vceVariance, *varVCEMean;
   double fixedInput, aMean, aVariance, **bsVCEs, ddata;
   double *txArray, *tyArray, *twArray, fixedInput2;
   double *varVCEVar, *mainEffects, *STI, *xLower, *xUpper;
   double *XIn, *YIn, *XX, *YY, *X, *Y;
   char   pString[500], winput[501], **inputNames;
   FILE   *fp, *fp1=NULL;
   PsuadeData    *ioPtr=NULL;
   RSConstraints *constrPtr=NULL;
   pData         pCorMat, pdata;
   psMatrix      *corMatp, corMat;

   nInputs  = adata.nInputs_;
   nInputs_ = nInputs;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   XIn      = adata.sampleInputs_;
   YIn      = adata.sampleOutputs_;
   outputID    = adata.outputID_;
   printLevel  = adata.printLevel_;
   nSubSamples = adata.nSubSamples_;
   ioPtr       = adata.ioPtr_;
   pdfFlags    = adata.inputPDFs_;

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR,"TwoParamEffect ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR,"    nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR,"    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nInputs <= 2)
   {
      printOutTS(PL_ERROR,"TwoParamEffect ERROR: there is no point in using\n");
      printOutTS(PL_ERROR,"              this method for nInputs <= 2.\n");
      return PSUADE_UNDEFINED;
   } 
   if (nSamples/nSubSamples*nSubSamples != nSamples)
   {
      printOutTS(PL_ERROR,"TwoParamEffect ERROR: nSamples != k*nSubSamples.\n");
      printOutTS(PL_ERROR,"    nSamples    = %d\n", nSamples);
      printOutTS(PL_ERROR,"    nSubSamples = %d\n", nSubSamples);
      return PSUADE_UNDEFINED;
   } 
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0)
   {
      printOutTS(PL_ERROR,"TwoParamEffect ERROR: invalid outputID.\n");
      printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR,"    outputID = %d\n", whichOutput+1);
      printOutTS(PL_ERROR,"    outputID reset to 1\n");
      whichOutput = 0;
   }
   if (ioPtr == NULL)
   {
      printOutTS(PL_ERROR,"TwoParamEffect ERROR: no data (PsuadeData).\n");
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
            printOutTS(PL_ERROR,
                 "* TwoParamEffect INFO: this method should not\n");
            printOutTS(PL_ERROR,
                 "*      be used if inputs are correlated with\n");
            printOutTS(PL_ERROR,
                 "*      joint PDFs. Use group variance-based method.\n");
            return PSUADE_UNDEFINED;
         }     
      }
   }
   status = 0;
   for (ss = 0; ss < nSamples; ss++)
      if (YIn[nOutputs*ss+whichOutput] > 0.9*PSUADE_UNDEFINED) status = 1;
   if (status == 1)
   {
      printOutTS(PL_ERROR, 
                 "TwoParamEffect ERROR: Some outputs are undefined.\n");
      printOutTS(PL_ERROR, 
                 "        Prune the undefined sample points first.\n");
      return PSUADE_UNDEFINED;
   }

   // clean up
   if (sensitivity_) delete[] sensitivity_;
   sensitivity_ = NULL;

   for (ii = 0; ii < nInputs; ii++)
   {
      if (pdfFlags != NULL && pdfFlags[ii] != PSUADE_PDF_UNIFORM)
      {
         printOutTS(PL_INFO,
              "* TwoParamEffect INFO: some inputs have non-uniform\n");
         printOutTS(PL_INFO,
              "*      PDFs. However, they will not be relevant in\n");
         printOutTS(PL_INFO,
              "*      this analysis (since the sample should have been\n");
         printOutTS(PL_INFO,
              "*      generated within the desired distributions.)\n");
         break;
      }
   }

   if (ioPtr != NULL)
   {
      constrPtr = new RSConstraints();
      constrPtr->genConstraints(ioPtr);
   }
   X = XIn;
   Y = new double[nSamples];
   checkAllocate(Y, "Y in TwoParam::analyze");
   ncount = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      Y[ss] = YIn[nOutputs*ss+whichOutput];
      ddata = constrPtr->evaluate(&(X[ss*nInputs]), Y[ss], status);
      if (status == 0) Y[ss] = PSUADE_UNDEFINED;
      else             ncount++;
   }
   if (ncount == 0)
   {
      printOutTS(PL_ERROR,"TwoParamEffect ERROR: no valid sample point.\n");
      printOutTS(PL_ERROR,"    nSamples before filtering = %d\n",nSamples);
      printOutTS(PL_ERROR,"    nSamples after  filtering = %d\n",ncount);
      printOutTS(PL_ERROR,
           "    INFO: check your data file for undefined's (1e35)\n");
      return 1.0;
   }
   if (ncount != nSamples)
   {
      printOutTS(PL_INFO,
           "* TwoParamEffect: nSamples before filtering = %d\n", nSamples);
      printOutTS(PL_INFO,
           "* TwoParamEffect: nSamples after filtering  = %d\n", ncount);
   }

   computeMeanVariance(nInputs,1,nSamples,Y,&aMean,&aVariance,0);
   outputMean_ = aMean;
   outputStd_  = sqrt(aVariance);

   if (PABS(aVariance) < 1.0e-15)
   {
      printOutTS(PL_INFO, 
           "* =====> TwoParamEffect: variance = %12.4e (sd =%12.4e)\n",
           aVariance, sqrt(aVariance));
      printOutTS(PL_ERROR, 
           "TwoParamEffect INFO: std dev too small ==> terminate.\n");
      return PSUADE_UNDEFINED;
   }
   if (psAnaExpertMode_ == 1 && printLevel > 4)
   {
      printOutTS(PL_INFO, "Perform crude analysis? (y or n) ");
      fgets(pString,400,stdin);
      if (pString[0] == 'y')
      {
         vce = new double[nInputs*(nInputs+1)];
         checkAllocate(vce, "vce in TwoParam::analyze");
         status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                  xLower, xUpper, vce);
         if (status == 0)
         {
            pData *pPtr = ioPtr->getAuxData();
            pPtr->nDbles_ = nInputs * nInputs;
            pPtr->dbleArray_ = new double[nInputs * (nInputs+1)];
            for (ii = 0; ii < nInputs*(nInputs+1); ii++)
               pPtr->dbleArray_[ii] = vce[ii];
            pPtr->dbleData_ = aVariance;
         }
         delete [] Y;
         delete [] vce;
         return 1.0;
      }
   }

   vceMean       = new double[nSubSamples];
   vceVariance   = new double[nSubSamples];
   bins          = new int[nSubSamples];
   vce           = new double[nInputs*(nInputs+1)];
   varVCEMean    = new double[nInputs*nInputs];
   varVCEVar     = new double[nInputs*nInputs];
   meanVCEVar    = new double[nInputs*nInputs];
   mainEffects   = new double[nInputs];
   STI           = new double[nInputs];
   checkAllocate(STI, "STI in TwoParam::analyze");
   errflag       = 0;

   printOutTS(PL_INFO,"\n");
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"* Second order Sensitivities on Raw Sample Data\n");
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO,
        "* Note: This method needs large samples for accurate results\n");
   printOutTS(PL_INFO,
        "*       (e.g. thousands or more depending on the functions).\n");
   printOutTS(PL_INFO,
        "*       For small to moderate sample sizes, use rssobol2.\n");
   printOutTS(PL_INFO,
        "* Note: This method works on replicated orthogonal arrays.\n");
   printOutTS(PL_INFO,
        "*       For random samples, a crude analysis will be performed.\n");
   printOutTS(PL_INFO,
        "*       (use analysis expert mode to set number of bins).\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,"* Number of sample points = %10d\n",nSamples);
   printOutTS(PL_INFO,"* Number of nInputs       = %10d\n",nInputs);
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,"Output %d\n", whichOutput+1);
   printOutTS(PL_INFO,"* ====> TwoParamEffect: mean     = %12.4e\n", aMean);
   printOutTS(PL_INFO,
        "* ====> TwoParamEffect: variance = %12.4e (sd =%12.4e)\n",
        aVariance, sqrt(aVariance));

   txArray  = new double[nSamples];
   twArray  = new double[nSamples];
   tyArray  = new double[nSamples];
   checkAllocate(tyArray, "tyArray in TwoParam::analyze");
   if (tyArray == NULL)
   {
      printOutTS(PL_ERROR,
           "TwoParamEffect ERROR:: memory allocation problem.\n");
      printOutTS(PL_ERROR,
           "                       Consult PSUADE developers.\n");
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         if (printLevel > 4 || nSamples > 100000)
            printOutTS(PL_DUMP,
                "TwoParamEffect:: input pairs = %d %d\n", ii+1, ii2+1);
         for (ss = 0; ss < nSamples; ss++)
         {
            txArray[ss] = X[nInputs*ss+ii];
            twArray[ss] = X[nInputs*ss+ii2];
            tyArray[ss] = Y[ss];
         }

         sortDbleList3(nSamples, txArray, twArray, tyArray);

         for (ss = 1; ss < nSamples; ss++)
            if (PABS(txArray[ss]-txArray[ss-1]) > 1.0E-10)
               break;
         nReps1 = ss;
         if (nReps1 <= 1)
         {
            printOutTS(PL_INFO,
               "TwoParamEffect INFO: nReps < 1 for input %d.\n",ii+1);
            printOutTS(PL_INFO,
               "        ==> not replicated orthogonal array nor Factorial\n");
            printOutTS(PL_INFO,
               "        ==> crude 2-way interaction analysis.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * (nInputs+1)];
               for (ii = 0; ii < nInputs*(nInputs+1); ii++) 
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps1)
         {
            fixedInput = txArray[ss];
            ncount = 0;
            for (ss2 = 0; ss2 < nReps1; ss2++)
               if (PABS(txArray[ss+ss2]-fixedInput) < 1.0E-10)
                  ncount++;
            if (ncount != nReps1) errflag = 1;
         }
         if (errflag != 0)
         {
            printOutTS(PL_WARN,
               "TwoParamEffect WARNING(2): replications not satisified.\n");
            printOutTS(PL_WARN,
               "    Are you using replicated orthogonal array or Factorial?\n");
            printOutTS(PL_WARN,
               "    If so, you need to use > 1 replications. A crude\n");
            printOutTS(PL_WARN,
               "    2-way interaction analysis will be done instead.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * (nInputs+1)];
               for (ii = 0; ii < nInputs*(nInputs+1); ii++) 
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps1)
            sortDbleList2(nReps1,&twArray[ss],&tyArray[ss]);

         for (ss = 1; ss < nReps1; ss++)
            if (PABS(twArray[ss]-twArray[ss-1]) > 1.0E-10)
               break;
         nReps = ss;
         if (nReps <= 1 || (nReps1/nReps*nReps != nReps1))
         {
            printOutTS(PL_WARN,
               "TwoParamEffect WARNING(3): replications not satisified.\n");
            printOutTS(PL_WARN,
               "    Are you using replicated orthogonal array or Factorial?\n");
            printOutTS(PL_WARN,
               "    If so, you need to use > 1 replications.\n");
            printOutTS(PL_WARN,
               "A Crude 2-way interaction analysis will be done instead.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * (nInputs+1)];
               for (ii = 0; ii < nInputs*(nInputs+1); ii++)   
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps1)
         {
            for (ss2 = 0; ss2 < nReps1; ss2+=nReps)
            {
               ncount = 1;
               fixedInput = twArray[ss+ss2];
               for (repID = 1; repID < nReps; repID++)
               {
                  if (PABS(twArray[ss+ss2+repID]-fixedInput)<1.0E-10)
                  ncount++;
               }
               if (ncount != nReps)
               {
                  printOutTS(PL_WARN,
                     "TwoParamEffect WARNING(3): sample not rOA/FACT.\n");
                  errflag = 1;
               }
            }
            if (errflag != 0) break;
         }
         if (errflag != 0)
         {
            printOutTS(PL_WARN,
               "TwoParamEffect WARNING(4): replications not satisified.\n");
            printOutTS(PL_WARN,
               "    Are you using replicated orthogonal array or Factorial?\n");
            printOutTS(PL_WARN,
               "    If so, you need to use > 1 replications.\n");
            printOutTS(PL_WARN,
               "A Crude 2-way interaction analysis will be done instead.\n");
            status = computeVCECrude(nInputs, nSamples, X, Y, aVariance,
                                     xLower, xUpper, vce);
            if (status == 0)
            {
               pData *pPtr = ioPtr->getAuxData();
               pPtr->nDbles_ = nInputs * nInputs;
               pPtr->dbleArray_ = new double[nInputs * (nInputs+1)];
               for (ii = 0; ii < nInputs*(nInputs+1); ii++)
                  pPtr->dbleArray_[ii] = vce[ii];
               pPtr->dbleData_ = aVariance;
            }
            delete [] vce;
            delete [] vceMean;
            delete [] vceVariance;
            delete [] varVCEMean;
            delete [] varVCEVar;
            delete [] meanVCEVar;
            delete [] mainEffects;
            delete [] bins;
            delete [] STI;
            delete [] txArray;
            delete [] tyArray;
            delete [] twArray;
            delete [] Y;
            return 1.0;
         }

         for (ss = 0; ss < nSamples; ss+=nReps)
         {
            symIndex = ss / nReps;
            fixedInput  = txArray[ss];
            fixedInput2 = twArray[ss];
            vceMean[symIndex] = 0.0;
            validnReps = 0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10 &&
                   PABS(twArray[ss+repID] - fixedInput2) < 1.0E-10)
               {
                  if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  {
                     validnReps++;
                     vceMean[symIndex] += tyArray[ss+repID];
                  }
               }
               else
               {
                  printOutTS(PL_ERROR,
                       "TwoParamEffect ERROR: consult developers.\n");
                  exit(1);
               }
            }
            bins[symIndex] = validnReps;
            if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
            vceVariance[symIndex] = 0.0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10 &&
                   PABS(twArray[ss+repID] - fixedInput2) < 1.0E-10)
               {
                  if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  {
                     vceVariance[symIndex] += 
                       ((tyArray[ss+repID]-vceMean[symIndex])*
                        (tyArray[ss+repID]-vceMean[symIndex]));
                  }
               }
            }
            if (validnReps > 0)
                 vceVariance[symIndex] /= ((double) validnReps);
            else vceVariance[symIndex] = PSUADE_UNDEFINED;
         }
         nSubs = nSamples / nReps;

         totalCnt = 0;
         for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];

         varVCEMean[ii*nInputs+ii2] = 0.0;
         meanVCEVar[ii*nInputs+ii2] = 0.0;
         aMean = 0.0;
         for (subID = 0; subID < nSubs; subID++)
         {
            if (vceVariance[subID] != PSUADE_UNDEFINED)
            {
               aMean += (vceMean[subID] / totalCnt * bins[subID]);
            }
         }
         for (subID = 0; subID < nSubs; subID++)
         {
            if (vceVariance[subID] != PSUADE_UNDEFINED)
            {
               varVCEMean[ii*nInputs+ii2] += ((vceMean[subID] - aMean)*
                                              (vceMean[subID] - aMean) *
                                              bins[subID] / totalCnt);
               meanVCEVar[ii*nInputs+ii2] += vceVariance[subID] *
                                             bins[subID] / totalCnt;
            }
         }
         vce[ii*nInputs+ii2] = varVCEMean[ii*nInputs+ii2];

         varVCEVar[ii*nInputs+ii2] = 0.0; 
         for (subID = 0; subID < nSubs; subID++)
            if (vceVariance[subID] != PSUADE_UNDEFINED)
               varVCEVar[ii*nInputs+ii2] += 
                  ((vceVariance[subID]-meanVCEVar[ii*nInputs+ii2])*
                   (vceVariance[subID]-meanVCEVar[ii*nInputs+ii2]) * 
                   bins[subID] / totalCnt);
         if (printLevel > 4 || nSamples > 100000) 
            printOutTS(PL_DUMP, "Interaction %3d %d = %e\n",ii+1,ii2+1,
                   vce[ii*nInputs+ii2]/aVariance);
      }
      if (errflag != 0) break;
   }

   for (ii = 0; ii < nInputs; ii++)
   {
      mainEffects[ii] = analyze1D(nInputs,nOutputs,nSamples,X,YIn,
                         ii,whichOutput,nReps1,aMean,aVariance,&STI[ii]);
      printOutTS(PL_INFO,"Main effect (Inputs %3d) = %12.4e\n", ii+1,
             mainEffects[ii]/aVariance);
      vce[nInputs*nInputs+ii] = mainEffects[ii];
   }

   pData *pPtr = ioPtr->getAuxData();
   pPtr->nDbles_ = nInputs * nInputs;
   pPtr->dbleArray_ = new double[nInputs * (nInputs+1)];
   for (ii = 0; ii < nInputs*(nInputs+1); ii++) pPtr->dbleArray_[ii] = vce[ii];
   pPtr->dbleData_ = aVariance;

#if 0
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         //varVCEMean[ii*nInputs+ii2] -= 
         //               (mainEffects[ii] + mainEffects[ii2]);
         varVCEMean[ii*nInputs+ii2] -= (STI[ii] + STI[ii2]);
         //if (varVCEMean[ii*nInputs+ii2] < 0.0) 
         //   varVCEMean[ii*nInputs+ii2] = 0.0;
      }
   }
#endif

   sensitivity_ = new double[nInputs_];
   checkAllocate(sensitivity_, "sensitivity_ in TwoParam::analyze");

   if (errflag == 0)
   {
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO,"* First order sensitivity indices (normalized)\n");
      printDashes(PL_INFO, 0);
      for (ii = 0; ii < nInputs; ii++)
      {
         if (mainEffects[ii] < 1.0e-10 * aVariance)
            mainEffects[ii] = 0.0;
         if (STI[ii] < 0.0) STI[ii] = 0.0;
         printOutTS(PL_INFO, 
                "* Sensitivity index (Inputs %3d) = %12.4e (raw = %12.4e)\n",
                ii+1,mainEffects[ii]/aVariance,mainEffects[ii]);

         //save sensitivity
         sensitivity_[ii] = mainEffects[ii];
         printf("sensitivity for input %d = %e\n", ii+1, sensitivity_[ii]);
      }
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO, 
           "* (First + second order) sensitivity indices (normalized)\n");
      printDashes(PL_INFO, 0);
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printOutTS(PL_INFO, 
             "* Sensitivity index ( Inputs %2d %2d ) = %12.4e (raw = %12.4e)\n",
             ii+1,ii2+1,vce[ii*nInputs+ii2]/aVariance,vce[ii*nInputs+ii2]);

            //save sensitivity comparisons
            record rec = {ii+1, ii2+1, vce[ii*nInputs+ii2]};
            compareSensitivity_.push_back(rec);
            printf("rec = %d %d %e\n", rec.index1, rec.index2, rec.value);
         }
      }
      if (constrPtr == NULL)
      {
         printEquals(PL_INFO, 0);
         printOutTS(PL_INFO,"* Just the second order sensitivity index\n");
         printOutTS(PL_INFO,
            "* (by subtracting mainEffects i, and just from VCE(i,j))\n");
         printOutTS(PL_INFO,
            "* Valid only for orthogonal (uncorrelated) inputs.\n");
         printDashes(PL_INFO, 0);
         for (ii = 0; ii < nInputs; ii++)
         {
            for (ii2 = ii+1; ii2 < nInputs; ii2++)
            {
               ddata = vce[ii*nInputs+ii2];
               ddata -= (mainEffects[ii] + mainEffects[ii2]);
               ddata /= aVariance;
               printOutTS(PL_INFO, 
                  "*2nd order sensitivity index (Inputs %2d %2d ) = %12.4e\n",
                  ii+1,ii2+1, ddata);
            }
         }
      }
   }

#if 0
   if (errflag == 0)
   {
      printOutTS(PL_INFO, 
         "*=======================================================**\n");
      printOutTS(PL_INFO, 
         "* Total variance due to the complementary set           **\n");
      printOutTS(PL_INFO, 
         "*-------------------------------------------------------**\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printOutTS(PL_INFO, 
               "*   Mean(Var) ratio   (Inputs %2d,%2d) = %12.4e\n",
               ii+1,ii2+1, meanVCEVar[ii*nInputs+ii2]/aVariance);
         }
      }
      printOutTS(PL_INFO, 
         "*=======================================================**\n");
      printOutTS(PL_INFO, 
         "* Total second order sensitivity index                  **\n");
      printOutTS(PL_INFO, 
         "*-------------------------------------------------------**\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printOutTS(PL_INFO, 
               "*   Sensitivity index (Inputs %2d,%2d) = %12.4e\n",
               ii+1,ii2+1, 1.0-meanVCEVar[ii*nInputs+ii2]/aVariance);
         }
      }
      printOutTS(PL_INFO, 
         "*=======================================================**\n");
      printOutTS(PL_INFO, 
         "* Strength of higher order interactions                 **\n");
      printOutTS(PL_INFO, 
         "*-------------------------------------------------------**\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
            printOutTS(PL_INFO, 
               "*   var VCE variance   (Inputs %2d,%2d) = %12.4e\n",
               ii+1,ii2+1,varVCEVar[ii*nInputs+ii2]);
         }
      }
   }
#endif
   printAsterisks(PL_INFO, 0);
   printAsterisks(PL_INFO, 0);

   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO, 
         "Bootstrap analysis draws n samples from the orignal sample\n");
      printOutTS(PL_INFO, 
         "and assess whether the sensitivity indices have converged.\n");
      printOutTS(PL_INFO, 
         "If you are performed iterative analysis with refinements,\n");
      printOutTS(PL_INFO, 
         "you will need to enter 'no index reuse' below at the first\n");
      printOutTS(PL_INFO, 
         "iteration and 'yes' afterward until the final refinement.\n");
      sprintf(pString,"Perform bootstrap interaction analysis? (y or n) ");
      getString(pString, winput);

      ncount = 0;
      if (winput[0] == 'y')
      {
         sprintf(pString,"Enter number of bootstrap samples (>=100): ");
         ncount = getInt(100, 2000, pString);
         XX = new double[nInputs*nSamples];
         YY = new double[nSamples];
         bsVCEs = new double*[ncount];
         for (ii = 0; ii < ncount; ii++)
            bsVCEs[ii] = new double[nInputs*nInputs];
         checkAllocate(bsVCEs[ncount-1], "bsVCEs in TwoParam::analyze");
         nReps = nSamples / nSubSamples;

         fp1 = fopen(".IE_bootstrap_indset", "r");
         if (fp1 != NULL)
         {
            printOutTS(PL_INFO, ".IE_bootstrap_indset file found.\n");
            sprintf(pString,"Re-use file? (y or n) ");
            getString(pString, winput);
            if (winput[0] == 'y')
            {
               fscanf(fp1, "%d", &ii);
               if (ii != nReps*ncount)
               {
                  printOutTS(PL_ERROR, 
                       "ERROR: expect the first line to be %d.\n",
                       nReps*ncount);
                  printOutTS(PL_ERROR, 
                       "       Instead found the first line to be %d\n",ii);
                  fclose(fp1);
                  exit(1);
               }
            }
            else
            {
               fclose(fp1);
               fp1 = fopen(".IE_bootstrap_indset", "w");
               if (fp1 == NULL)
               {
                  printOutTS(PL_ERROR, 
                     "ERROR: cannot open bootstrap_indset file (write).\n");
                  exit(1);
               }
               fprintf(fp1, "%d\n", nReps*ncount);
            }
         }
         else
         {
            fp1 = fopen(".IE_bootstrap_indset", "w");
            if (fp1 == NULL)
            {
               printOutTS(PL_ERROR, 
                    "ERROR: cannot open bootstrap_indset file (write).\n");
               exit(1);
            }
            fprintf(fp1, "%d\n", nReps*ncount);
            winput[0] = 'n';
         }
         fp = NULL;

         for (ii = 0; ii < ncount; ii++)
         {
            if (ii % (ncount / 10) == 0)
               printOutTS(PL_INFO, "Processing %d(a) of %d\n", ii+1, ncount);
            for (ir = 0; ir < nReps; ir++)
            {
               if (fp1 != NULL && winput[0] == 'y')
                    fscanf(fp1, "%d\n", &index);
               else index = PSUADE_rand() % nReps;
               if (fp1 != NULL && winput[0] != 'y')
                  fprintf(fp1, "%d\n", index);

               for (ss = 0; ss < nSubSamples; ss++)
               {
                  for (ii2 = 0; ii2 < nInputs; ii2++)
                     XX[(ir*nSubSamples+ss)*nInputs+ii2] =
                        XIn[(index*nSubSamples+ss)*nInputs+ii2];
                  YY[ir*nSubSamples+ss] = 
                     YIn[(index*nSubSamples+ss)*nOutputs+whichOutput];
               }
            }

            if (ii % (ncount / 10) == 0)
               printOutTS(PL_INFO, "Processing %d(b) of %d\n", ii+1, ncount);
            computeMeanVariance(nInputs,iOne,nSamples,YY,&aMean,&aVariance,
                                iZero);
            if (PABS(aVariance) < 1.0e-15)
            {
               printOutTS(PL_INFO,"INFO: variance too small ==> terminate.\n");
               continue;
            }
            computeVCE2(nInputs, nSamples, nSubSamples, XX, YY, fp, aMean, 
                        varVCEMean, meanVCEVar, varVCEVar, vce);
            for (ii2 = 0; ii2 < nInputs; ii2++)
               mainEffects[ii2] = analyze1D(nInputs,iOne,nSamples,XX,YY,
                                   ii2,iZero,nReps1,aMean,aVariance,&STI[ii2]);
            for (ii2 = 0; ii2 < nInputs; ii2++)
               for (ss = ii2+1; ss < nInputs; ss++)
                  bsVCEs[ii][ii2*nInputs+ss] = (vce[ii2*nInputs+ss])/aVariance; 
         }
         if (fp1 != NULL) fclose(fp1);
         if (psPlotTool_ == 0)
         {
           printf("Enter name of matlab/scilab file to store bootstrap info: ");
           scanf("%500s", winput);
           fgets(pString,500,stdin);
           fp = fopen(winput, "w");
         }
         else fp = NULL;
         if (fp != NULL)
         {
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "//bootstrap sample of interaction effects\n");
               fprintf(fp, "//VCE((i-1)*n+j,k) = VCE(i,j), bs sample k\n");
            }
            else
            {
               fprintf(fp, "%%bootstrap sample of interaction effects\n");
               fprintf(fp, "%%VCE((i-1)*n+j,k) = VCE(i,j), bs sample k\n");
            }
            fprintf(fp, "VCE = zeros(%d,%d);\n",nInputs*nInputs,ncount);
            if (ioPtr != NULL) ioPtr->getParameter("input_names",pdata);
            if (pdata.strArray_ != NULL) inputNames = pdata.strArray_;
            else                         inputNames = NULL;
            if (inputNames == NULL)
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
                  if (inputNames[ii] != NULL)
                     fprintf(fp,"'%s',",inputNames[ii]);
                  else fprintf(fp,"'X%d',",ii+1);
               }
               if (inputNames[nInputs-1] != NULL)
                  fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
               else fprintf(fp,"'X%d'};\n",nInputs);
            }
            for (ss = 0; ss < ncount; ss++)
            {
               for (ii = 0; ii < nInputs; ii++)
               {
                  for (ii2 = ii+1; ii2 < nInputs; ii2++)
                  {
                     if (ss == 0)
                     {
                        if (psPlotTool_ == 1)
                             fprintf(fp,"//Interaction(%d,%d) \n",ii+1,ii2+1);
                        else fprintf(fp,"%%Interaction(%d,%d) \n",ii+1,ii2+1);
                     }
                     fprintf(fp, "VCE(%d,%d) = %e;\n",ii*nInputs+ii2+1,ss+1,
                             bsVCEs[ss][ii*nInputs+ii2]);
                  }
               }
            }
            fprintf(fp,"VMM = sum(VCE')/%d;\n", ncount);
            fprintf(fp,"VMA = max(VCE');\n");
            fprintf(fp,"VMI = min(VCE');\n");
            fprintf(fp,"ymin = 0;\n");
            fprintf(fp,"ymax = max(VMA);\n");
            fprintf(fp,"hh = 0.05 * (ymax - ymin);\n");
            fprintf(fp,"ymax = ymax + hh;\n");
            fprintf(fp,"nn   = %d;\n",nInputs);
            if (psPlotTool_ == 1)
            {
               fprintf(fp,"VMM = matrix(VMM, nn, nn);\n");
               fprintf(fp,"VMM = VMM';\n");
               fprintf(fp,"VMA = matrix(VMA, nn, nn);\n");
               fprintf(fp,"VMA = VMA';\n");
               fprintf(fp,"VMI = matrix(VMI, nn, nn);\n");
               fprintf(fp,"VMI = VMI';\n");
               fprintf(fp,"drawlater\n");
               fprintf(fp,"hist3d(VMM);\n");
               fwriteHold(fp, 1);
               fprintf(fp,"a=gca();\n");
               fprintf(fp,"a.data_bounds=[0, 0, 0; nn, nn+1, ymax];\n");
               fprintf(fp,"newtick = a.x_ticks;\n");
               fprintf(fp,"newtick(2) = [1:nn]';\n");
               fprintf(fp,"newtick(3) = Str';\n");
               fprintf(fp,"a.x_ticks = newtick;\n");
               fprintf(fp,"a.x_label.font_size = 3;\n");
               fprintf(fp,"a.x_label.font_style = 4;\n");
               fprintf(fp,"a.y_ticks = newtick;\n");
               fprintf(fp,"a.y_label.font_size = 3;\n");
               fprintf(fp,"a.y_label.font_style = 4;\n");
               fprintf(fp,"a.rotation_angles = [5 -70];\n");
               fprintf(fp,"drawnow\n");
            }
            else
            {
               fprintf(fp,"VMM = reshape(VMM, nn, nn);\n");
               fprintf(fp,"VMM = VMM';\n");
               fprintf(fp,"VMA = reshape(VMA, nn, nn);\n");
               fprintf(fp,"VMA = VMA';\n");
               fprintf(fp,"VMI = reshape(VMI, nn, nn);\n");
               fprintf(fp,"VMI = VMI';\n");
               fprintf(fp,"hh = bar3(VMM,0.8);\n");
               fprintf(fp,"alpha = 0.2;\n");
               fprintf(fp,"set(hh,'FaceColor','b','facea',alpha);\n");
               fprintf(fp,"[X,Y] = meshgrid(1:nn,1:nn);\n");
               fprintf(fp,"for k = 1:nn\n");
               fprintf(fp,"  for l = k:nn\n");
               fprintf(fp,"    mkl = VMM(k,l);\n");
               fprintf(fp,"    ukl = VMA(k,l);\n");
               fprintf(fp,"    lkl = VMI(k,l);\n");
               fprintf(fp,"    if (mkl > .02 & (ukl-lkl)/mkl > .02)\n");
               fprintf(fp,"      xkl = [X(k,l), X(k,l)];\n");
               fprintf(fp,"      ykl = [Y(k,l), Y(k,l)];\n");
               fprintf(fp,"      zkl = [lkl, ukl];\n");
               fprintf(fp,"      plot3(xkl,ykl,zkl,'-mo',...\n");
               fprintf(fp,"        'LineWidth',5,'MarkerEdgeColor','k',...\n");
               fprintf(fp,"        'MarkerFaceColor','k','MarkerSize',10);\n");
               fprintf(fp,"    end\n");
               fprintf(fp,"  end\n");
               fprintf(fp,"end\n");
               fwriteHold(fp, 0);
               fprintf(fp,"axis([0.5 nn+0.5 0.5 nn+0.5 0 ymax])\n");
               fprintf(fp,"set(gca,'XTickLabel',Str);\n");
               fprintf(fp,"set(gca,'YTickLabel',Str);\n");
               fprintf(fp,"set(gca, 'fontsize', 12)\n");
               fprintf(fp,"set(gca, 'fontweight', 'bold')\n");
               fprintf(fp,"set(gca, 'linewidth', 2)\n");
            }
            fwritePlotAxes(fp);
            fwritePlotTitle(fp,"Sobol 1st+2nd Order Indices (with bootstrap)");
            fwritePlotZLabel(fp, "Sobol Indices (Normalized)");
            fwritePlotXLabel(fp, "Inputs");
            fwritePlotYLabel(fp, "Inputs");
            printOutTS(PL_INFO, "Total variance = %e\n", aVariance);
            fclose(fp);
            printOutTS(PL_INFO, 
               "Bootstrapped interaction effect plot in now in %s\n",winput);
         }
         else
         {
            printOutTS(PL_ERROR, "ERROR: cannot open file %s\n", winput);
         }
         for (ii = 0; ii < ncount; ii++) delete [] bsVCEs[ii];
         delete [] bsVCEs;
         delete [] XX;
         delete [] YY;
      }
   }

   delete [] vce;
   delete [] varVCEMean;
   delete [] meanVCEVar;
   delete [] varVCEVar;
   delete [] vceMean;
   delete [] vceVariance;
   delete [] mainEffects;
   delete [] STI;
   delete [] txArray;
   delete [] twArray;
   delete [] tyArray;
   delete [] Y;
   delete [] bins;
   if (constrPtr != NULL) delete constrPtr;
   // return 1.0 to facilitate continuous refinement
   return 1.0;
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int TwoParamAnalyzer::computeMeanVariance(int nInputs, int nOutputs, 
                                     int nSamples, double *y, double *aMean, 
                                     double *aVariance, int outputID)
{
   int    ss, count;
   double mean, variance;

   mean = 0.0;
   count = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      if (y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         mean += y[nOutputs*ss+outputID];
         count++;
      }

   }
   if (count <= 0)
   {
      printOutTS(PL_ERROR, "TwoParamEffect ERROR: no valid data.\n");
      exit(1);
   }
   mean /= (double) count;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      if (y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         variance += ((y[nOutputs*ss+outputID] - mean) *
                      (y[nOutputs*ss+outputID] - mean));
      }
   }
   variance /= (double) count;
   (*aMean) = mean;
   (*aVariance) = variance;
   return 0;
}

// ************************************************************************
// perform VCE (McKay) analysis
// ------------------------------------------------------------------------
double TwoParamAnalyzer::analyze1D(int nInputs, int nOutputs, int nSamples,
                           double *x, double *y, int inputID, int outputID,
                           int nReps, double aMean, double aVariance,
                           double *retData)
{
   int    ss, repID, subID, ncount, symIndex, nSubs;
   int    whichOutput, validnReps, totalCnt, *bins;
   double meanVCEVar, vce, *vceMean, *vceVariance, varVCEMean;
   double fixedInput, *txArray, *tyArray, tmean;

   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;

   if (nReps <= 1) return 1.0;
   nSubs         = nSamples / nReps;
   vceMean       = new double[nSubs];
   vceVariance   = new double[nSubs];
   txArray       = new double[nSamples];
   tyArray       = new double[nSamples];
   bins          = new int[nSubs];
   checkAllocate(bins, "bins in TwoParam::analyze1D");

   for (ss = 0; ss < nSamples; ss++)
   {
      txArray[ss] = x[nInputs*ss+inputID];
      tyArray[ss] = y[nOutputs*ss+whichOutput];
   }
   sortDbleList2(nSamples, txArray, tyArray);
   ncount = 0;
   for (ss = 0; ss < nSamples; ss+=nReps)
   {
      symIndex = ss / nReps;
      fixedInput = txArray[ss];
      vceMean[symIndex] = 0.0;
      validnReps = 0;
      for (repID = 0; repID < nReps; repID++)
      {
         if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10)
         {
            ncount++;
            if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
            {
               vceMean[symIndex] += tyArray[ss+repID];
               validnReps++;
            }
         }
      }
      if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
      vceVariance[symIndex] = 0.0;
      tmean = vceMean[symIndex];
      for (repID = 0; repID < nReps; repID++)
      {
         if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10)
         {
            if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
               vceVariance[symIndex] += ((tyArray[ss+repID]-tmean)*
                                         (tyArray[ss+repID]-tmean));
         }
      }
      if (validnReps > 0) vceVariance[symIndex] /= validnReps;
      else                vceVariance[symIndex] = PSUADE_UNDEFINED;
      bins[symIndex] = validnReps;
   }
   if (ncount != nSamples)
   {
      printOutTS(PL_ERROR,"TwoParamEffect ERROR: sample not valid, input %d\n",
             inputID);
      printOutTS(PL_ERROR, "       error data = %d (%d)\n",ncount, nSamples);
      exit(1);
   }

   varVCEMean = 0.0;
   meanVCEVar = 0.0;
   totalCnt = 0;
   for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];
   aMean = 0.0;
   for (subID = 0; subID < nSubs; subID++)
   {
      if (vceVariance[subID] != PSUADE_UNDEFINED)
         aMean += (vceMean[subID] / totalCnt * bins[subID]);
   }
   for (subID = 0; subID < nSubs; subID++)
   {
      if (vceVariance[subID] != PSUADE_UNDEFINED)
      {
         varVCEMean += ((vceMean[subID]-aMean) * (vceMean[subID]-aMean) *
                        bins[subID] / totalCnt);
         meanVCEVar += vceVariance[subID] * bins[subID] / totalCnt;
      }
   }
   vce = varVCEMean;

   delete [] txArray;
   delete [] tyArray;
   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   (*retData) = 1.0 - meanVCEVar / aVariance;
   return vce;
}

// *************************************************************************
// Compute two-way VCE 
// -------------------------------------------------------------------------
int TwoParamAnalyzer::computeVCE2(int nInputs, int nSamples, int nSubSamples,
                                 double *X, double *Y, FILE *fp, 
                                 double aMean, double *varVCEMean,
                                 double *meanVCEVar, double *varVCEVar,
                                 double *vce)
{
   int    ii, ii2, ss, nReps1, nReps, nSubs, symIndex, repID, subID;
   int    totalCnt, *bins, validnReps;
   double *txArray, *twArray, *tyArray, fixedInput, fixedInput2, *vceMean;
   double *vceVariance;


   txArray = new double[nSamples];
   twArray = new double[nSamples];
   tyArray = new double[nSamples];
   vceMean     = new double[nSubSamples];
   vceVariance = new double[nSubSamples];
   bins        = new int[nSubSamples];
   checkAllocate(bins, "bins in TwoParam::computeVCE2");

   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         for (ss = 0; ss < nSamples; ss++)
         {
            txArray[ss] = X[nInputs*ss+ii];
            twArray[ss] = X[nInputs*ss+ii2];
            tyArray[ss] = Y[ss];
         }

         sortDbleList3(nSamples, txArray, twArray, tyArray);

         for (ss = 1; ss < nSamples; ss++)
            if (PABS(txArray[ss]-txArray[ss-1]) > 1.0E-10)
               break;
         nReps1 = ss;

         for (ss = 0; ss < nSamples; ss+=nReps1)
            sortDbleList2(nReps1,&twArray[ss],&tyArray[ss]);

         for (ss = 1; ss < nReps1; ss++)
            if (PABS(twArray[ss]-twArray[ss-1]) > 1.0E-10)
               break;
         nReps = ss;

         for (ss = 0; ss < nSamples; ss+=nReps)
         {
            symIndex = ss / nReps;
            fixedInput  = txArray[ss];
            fixedInput2 = twArray[ss];
            vceMean[symIndex] = 0.0;
            validnReps = 0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-10 &&
                   PABS(twArray[ss+repID] - fixedInput2) < 1.0E-10)
               {
                  if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  {
                     vceMean[symIndex] += tyArray[ss+repID];
                     validnReps++;
                  }
               }
            }
            if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
            vceVariance[symIndex] = 0.0;
            for (repID = 0; repID < nReps; repID++)
            {
               if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
                  vceVariance[symIndex] += 
                       ((tyArray[ss+repID]-vceMean[symIndex])*
                        (tyArray[ss+repID]-vceMean[symIndex]));
            }
            if (validnReps > 0) 
                 vceVariance[symIndex] /= ((double) validnReps);
            else vceVariance[symIndex] = PSUADE_UNDEFINED;
         }
         nSubs = nSamples / nReps;

         varVCEMean[ii*nInputs+ii2] = 0.0;
         meanVCEVar[ii*nInputs+ii2] = 0.0;
         totalCnt = 0;
         for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];
         for (subID = 0; subID < nSubs; subID++)
         {
            if (vceVariance[subID] != PSUADE_UNDEFINED)
            {
               varVCEMean[ii*nInputs+ii2] += ((vceMean[subID] - aMean)*
                                              (vceMean[subID] - aMean) *
                                              bins[subID] / totalCnt);
               meanVCEVar[ii*nInputs+ii2] += vceVariance[subID]*bins[subID]/
                                             totalCnt;
            }
         }
         vce[ii*nInputs+ii2] = varVCEMean[ii*nInputs+ii2];

         varVCEVar[ii*nInputs+ii2] = 0.0; 
         for (subID = 0; subID < nSubs; subID++)
            if (vceVariance[subID] != PSUADE_UNDEFINED)
               varVCEVar[ii*nInputs+ii2] += 
                  ((vceVariance[subID]-meanVCEVar[ii*nInputs+ii2])*
                   (vceVariance[subID]-meanVCEVar[ii*nInputs+ii2]) *
                   bins[subID] / totalCnt);
      }
   }
   delete [] txArray;
   delete [] twArray;
   delete [] tyArray;
   delete [] vceMean;
   delete [] bins;
   delete [] vceVariance;
   return 0;
}

// *************************************************************************
// Compute VCE 
// -------------------------------------------------------------------------
int TwoParamAnalyzer::computeVCECrude(int nInputs, int nSamples, 
                             double *X, double *Y,double aVariance,
                             double *iLowerB, double *iUpperB, double *vce)
{
   int    ii, ii2, ss, nSize, nIntervals, ind1, ind2, ind;
   int    totalCnt, *bins, *tags, nIntervals1, nSize1;
   double *vceMean, *vceVariance, ddata, hstep1, hstep2, aMean, *vce1;
   char   pString[500];

   nSize = 50;
   ddata = 1.0 * nSamples / nSize;
   ddata = pow(ddata, 0.5);
   nIntervals = (int) ddata;
   nIntervals1 = (int) sqrt(1.0 * nSamples);
   if (nIntervals < 4)
   {
      printOutTS(PL_ERROR,"ERROR: sample size too small.\n");
      printOutTS(PL_ERROR,
         "       Need larger sample to give meaningful results.\n");
      return -1;
   }
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"           Crude 2-way Interaction Effect\n");
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO,
        "TwoParamEffect: default number of levels in each input  = %d\n",
          nIntervals);
   printOutTS(PL_INFO,
        "TwoParamEffect: default number of level for main effect = %d\n",
          nIntervals1);
   printOutTS(PL_INFO,
        "TwoParamEffect: sample size for each 2-dimensional box  = %d\n",
          nSize);
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,
        "* For small to moderate sample sizes, this method gives rough\n");
   printOutTS(PL_INFO,
        "* estimates of interaction effect (second order sensitivity).\n");
   printOutTS(PL_INFO,
        "* These estimates can vary with different choices of internal\n");
   printOutTS(PL_INFO,
        "* settings. For example, try different number of levels to\n");
   printOutTS(PL_INFO,
        "* assess sensitivity of the computed measures with respect to\n");
   printOutTS(PL_INFO,
        "* it. Turn on analysis expert mode to change the settings.\n");
   printOutTS(PL_INFO,
        "* Recommend sample size per 2D box to be at least 10 to\n");
   if (psAnaExpertMode_ == 1)
   {
      printEquals(PL_INFO, 0);
      sprintf(pString,
           "Number of levels for 2-way analysis (>= 4, default = %d): ",
           nIntervals);
      nIntervals = getInt(4, nSamples, pString);
      nSize = nSamples / nIntervals / nIntervals;
      if (nSize < 10)
      {
         printOutTS(PL_INFO,
              "* Using number of levels = %d for 2-way effect will give\n",
              nIntervals);
         printOutTS(PL_INFO,
            "  each 2D box less than 10 sample points. That is too small.\n");
         nSize = 10;
         ddata = 1.0 * nSamples / nSize;
         ddata = pow(ddata, 0.5);
         nIntervals = (int) ddata;
         printOutTS(PL_INFO, 
            "* Number of levels for 2-way analysis defaulted to %d\n", 
            nIntervals);
      } 
      sprintf(pString,
              "Number of levels for main effect (>= %d): ",nIntervals);
      nIntervals1 = getInt(nIntervals, nSamples, pString);
      nSize1 = nSamples / nIntervals1;
      if (nSize1 < 10)
      {
         printOutTS(PL_INFO, 
            "Sample size per level for main effect d too small (%d < 10).\n",
            nSize1);
         nSize1 = 10;
         nIntervals1 = nSamples / nSize1;
         printOutTS(PL_INFO, 
            "Default number of levels for main effect to %d\n",nIntervals1);
      } 
   }
   int nsize = nIntervals * nIntervals;
   if (nIntervals1 > nsize) nsize = nIntervals1; 
   vceMean     = new double[nsize];
   vceVariance = new double[nsize];
   bins        = new int[nsize];
   tags        = new int[nSamples];
   vce1        = new double[nInputs];
   checkAllocate(vce1, "vce1 in TwoParam::computeVCECrude");

   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO, "* First order sensitivity indices (normalized)\n");
   printDashes(PL_INFO, 0);

   sensitivity_ = new double[nInputs_];
   checkAllocate(sensitivity_, "sensitivity_ in TwoParam::computeVCECrude");
   for (ii = 0; ii < nInputs; ii++)
   {
      hstep1 = (iUpperB[ii] - iLowerB[ii]) / nIntervals1;
      for (ss = 0; ss < nIntervals1; ss++)
      {
         vceMean[ss] = 0.0;
         vceVariance[ss] = 0.0;
         bins[ss] = 0;
      }
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = (X[nInputs*ss+ii] - iLowerB[ii]) / hstep1;
         ind1 = (int) ddata;
         if (ind1 <  0) ind1 = 0;
         if (ind1 >= nIntervals1) ind1 = nIntervals1 - 1;
         tags[ss] = -1;
         if (Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceMean[ind1] += Y[ss];
            bins[ind1]++;
            tags[ss] = ind1;
         }
      }
      for (ss = 0; ss < nIntervals1; ss++)
         if (bins[ss] > 0) vceMean[ss] /= (double) bins[ss];

      for (ss = 0; ss < nSamples; ss++)
      {
         ind = tags[ss];
         if (ind >= 0 && Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceVariance[ind] += ((Y[ss]-vceMean[ind])*
                                 (Y[ss]-vceMean[ind]));
         }
      }
      for (ss = 0; ss < nIntervals1; ss++)
      {
         if (bins[ss] > 0)
            vceVariance[ss] /= (double) bins[ss];
         else vceVariance[ss] = PSUADE_UNDEFINED;
      }

      totalCnt = 0;
      for (ss = 0; ss < nIntervals1; ss++) totalCnt += bins[ss];
      aMean = 0.0;
      for (ss = 0; ss < nIntervals1; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
            aMean += (vceMean[ss] / totalCnt * bins[ss]);
      }
      vce1[ii] = 0.0;
      for (ss = 0; ss < nIntervals1; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
         {
            vce1[ii] += ((vceMean[ss] - aMean) *
                         (vceMean[ss] - aMean) * bins[ss]/totalCnt);
         }
      }
      printOutTS(PL_INFO,
           "* Sensitivity index (Inputs %3d) = %12.4e (raw = %12.4e)\n",
           ii+1,vce1[ii]/aVariance,vce1[ii]);

      //save sensitivity
      sensitivity_[ii] = vce1[ii];
      //cout << "sensitivity for input " << ii+1 << " = " << 
      //        sensitivity_[ii] << endl;
   }
   printEquals(PL_INFO, 0);
   printOutTS(PL_INFO, 
        "* (First + second order) sensitivity indices (normalized)\n");
   printDashes(PL_INFO, 0);
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = 0; ii2 < ii+1; ii2++) vce[ii*nInputs+ii2] = 0;
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         hstep1 = (iUpperB[ii] - iLowerB[ii]) / nIntervals;
         hstep2 = (iUpperB[ii2] - iLowerB[ii2]) / nIntervals;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            vceMean[ss] = 0.0;
            vceVariance[ss] = 0.0;
            bins[ss] = 0;
         }
         for (ss = 0; ss < nSamples; ss++)
         {
            ddata = (X[nInputs*ss+ii] - iLowerB[ii]) / hstep1;
            ind1 = (int) ddata;
            if (ind1 < 0) ind1 = 0;
            if (ind1 >= nIntervals) ind1 = nIntervals - 1;
            ddata = (X[nInputs*ss+ii2] - iLowerB[ii2]) / hstep2;
            ind2 = (int) ddata;
            if (ind2 < 0) ind2 = 0;
            if (ind2 >= nIntervals) ind2 = nIntervals - 1;
            ind = ind1 * nIntervals + ind2;
            tags[ss] = -1;
            if (Y[ss] < 0.9*PSUADE_UNDEFINED)
            {
               vceMean[ind] += Y[ss];
               bins[ind]++;
               tags[ss] = ind;
            }
         }

         for (ss = 0; ss < nIntervals*nIntervals; ss++)
            if (bins[ss] > 0) vceMean[ss] /= (double) bins[ss];

         for (ss = 0; ss < nSamples; ss++)
         {
            ind = tags[ss];
            if (ind >= 0 && Y[ss] < 0.9*PSUADE_UNDEFINED)
            {
               vceVariance[ind] += ((Y[ss]-vceMean[ind])*
                                    (Y[ss]-vceMean[ind]));
            }
         }
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            if (bins[ss] > 0)
                 vceVariance[ss] /= (double) bins[ss];
            else vceVariance[ss] = PSUADE_UNDEFINED;
         }

         totalCnt = 0;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
            totalCnt += bins[ss];
         aMean = 0.0;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            if (vceVariance[ss] != PSUADE_UNDEFINED)
               aMean += (vceMean[ss] / totalCnt * bins[ss]);
         }
         vce[ii*nInputs+ii2] = 0.0;
         for (ss = 0; ss < nIntervals*nIntervals; ss++)
         {
            if (vceVariance[ss] != PSUADE_UNDEFINED)
            {
               vce[ii*nInputs+ii2] += ((vceMean[ss] - aMean) *
                             (vceMean[ss] - aMean) * bins[ss]/totalCnt);
            }
         }
         if ( vce[ii*nInputs+ii2] < 0.0) vce[ii*nInputs+ii2] = 0.0;
         printOutTS(PL_INFO, 
            "* Sensitivity index ( Inputs %2d %2d ) = %12.4e (raw = %12.4e)\n",
            ii+1,ii2+1,vce[ii*nInputs+ii2]/aVariance,vce[ii*nInputs+ii2]);
         //           vce[ii*nInputs+ii2]/aVariance);

         //save sensitivity comparisons
         record rec = {ii+1, ii2+1, vce[ii*nInputs+ii2]};
         compareSensitivity_.push_back(rec);
      }
      vce[nInputs*nInputs+ii] = vce1[ii];
   }
   printAsterisks(PL_INFO, 0);
   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   delete [] tags;
   delete [] vce1;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
TwoParamAnalyzer& TwoParamAnalyzer::operator=(const TwoParamAnalyzer &)
{
   printOutTS(PL_ERROR, 
              "TwoParamEffect operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int TwoParamAnalyzer::get_nInputs()
{
   return nInputs_;
}
double TwoParamAnalyzer::get_outputMean()
{
   return outputMean_;
}
double TwoParamAnalyzer::get_outputStd()
{
   return outputStd_;
}
double *TwoParamAnalyzer::get_sensitivity()
{
   double* retVal = NULL;
   if (sensitivity_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in TwoParam::get_sensitivity");
      for (int ii = 0; ii < nInputs_; ii++) retVal[ii] = sensitivity_[ii];
   }
   return retVal;
}
std::vector<TwoParamAnalyzer::record> TwoParamAnalyzer::get_compareSensitivity()
{
   return compareSensitivity_;
}

