// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Functions for the class Optimizer
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Globals.h"
#include "PsuadeUtil.h"
#include "FuncApprox.h"
#include "pData.h"
#include "Optimizer.h"
#include "APPSPACKOptimizer.h"
#include "MMOptimizer.h"
#include "CobylaOptimizer.h"
#include "BobyqaOptimizer.h"
#include "SMOptimizer.h"
#include "MinpackOptimizer.h"
#include "TxMathOptimizer.h"
#include "SCEOptimizer.h"
#include "MultiObjectiveOptimizer.h"
#include "OUU1Optimizer.h"
#include "OUU2Optimizer.h"
#include "OUUOptimizer.h"
#include "LincoaOptimizer.h"
#include "NewuoaOptimizer.h"
#include "LBFGSOptimizer.h"
#include "NomadOptimizer.h"
#include "OUUMinlpOptimizer.h"

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
Optimizer::Optimizer()
{
   optimalX_ = NULL;
   optimalY_ = 1e35;
   numEvals_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Optimizer::~Optimizer()
{
   if (optimalX_ != NULL) delete [] optimalX_;
}

// ************************************************************************
// set objective function
// ------------------------------------------------------------------------
void Optimizer::setObjectiveFunction(void (*func)(int,double*,int,double*))
{
   objFunction_ = func;
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void Optimizer::optimize(oData *odata)
{
   (void) odata;
   return;
}

// ************************************************************************
// set parameter
// ------------------------------------------------------------------------
void Optimizer::setParam(char *)
{
   return;
}

// ************************************************************************
// get optimal X 
// ------------------------------------------------------------------------
double *Optimizer::getOptimalX()
{
   return optimalX_;
}

// ************************************************************************
// get optimal Y 
// ------------------------------------------------------------------------
double Optimizer::getOptimalY()
{
   return optimalY_;
}

// ************************************************************************
// get number of function evaluations 
// ------------------------------------------------------------------------
int Optimizer::getNumFuncEvals()
{
   return numEvals_;
}

// ************************************************************************
// perform optimization 
// ------------------------------------------------------------------------
int OptimizerSearch(PsuadeData *psuadeIO, FunctionInterface *funcIO,
                    double **optData, int *optCount)
{
   int    optimizeFlag, optimizeNumPts, optimizeUseRS, optMethod;
   int    minIndex, nPtsPerDim, faLeng, *dimIndices, *dimFactors;
   int    minCheck, returnFlag=0, ss, numYmin, istart, optimalCount;
   int    ii, ind1, ind2, nSamples, nInputs, nOutputs, totalEval;
   int    optID, optPrintLevel, optNewCnt, optNumFmin, M1=-1;
   double *tempOutputs, *dbleInd, *faXOut, *faYOut, optCutOff, optTol;
   double Ymin, *tmpOptX, *tmpOptY, *iLowerB, *iUpperB;
   double *sampleInputs, *sampleOutputs, optFmin, *tmpInitX, *tmpInitY;
   double *optInitXData, *optInitYData, *optimalXData, *optimalYData;
   char   specsFile[200], sparam[100];
   pData  pPtr, pLowerB, pUpperB, pInpData, pOutData;
   FuncApprox *faPtr;
#ifdef HAVE_TXMATH
   TxMathOptimizer   *TxMathPtr;
#endif
#ifdef HAVE_APPSPACK
   double *tempInputs=NULL;
   APPSPACKOptimizer *APPSPACKPtr;
#endif
#ifdef HAVE_MINPACK
   MinpackOptimizer  *MinpackPtr;
#endif
   CobylaOptimizer   *CobylaPtr;
   SMOptimizer       *SMPtr;
   MMOptimizer       *MMPtr;
   BobyqaOptimizer   *BobyqaPtr;
   SCEOptimizer      *SCEPtr;
   MultiObjectiveOptimizer *MooPtr;
   OUU1Optimizer     *OUU1Ptr;
   OUU2Optimizer     *OUU2Ptr;
   OUUOptimizer      *OUUPtr;
   LincoaOptimizer   *LincoaPtr;
   NewuoaOptimizer   *NewuoaPtr;
   LBFGSOptimizer    *lbfgsPtr;
   NomadOptimizer    *nomadPtr;
   OUUMinlpOptimizer *minlpPtr;
   oData             odata;

   psuadeIO->getParameter("input_ninputs", pPtr);
   nInputs = pPtr.intData_;
   psuadeIO->getParameter("output_noutputs", pPtr);
   nOutputs = pPtr.intData_;
   psuadeIO->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;
   psuadeIO->getParameter("input_sample", pInpData);
   sampleInputs = pInpData.dbleArray_;
   psuadeIO->getParameter("output_sample", pOutData);
   sampleOutputs = pOutData.dbleArray_;

   psuadeIO->getParameter("ana_opt_switch", pPtr);
   optimizeFlag = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_nlocalmin", pPtr);
   optimizeNumPts = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_rstype", pPtr);
   optimizeUseRS = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_outputid", pPtr);
   optID = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_printlevel", pPtr);
   optPrintLevel = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_method", pPtr);
   optMethod = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_cutoff", pPtr);
   optCutOff = pPtr.dbleData_;
   psuadeIO->getParameter("ana_opt_tolerance", pPtr);
   optTol = pPtr.dbleData_;
   psuadeIO->getParameter("ana_opt_fmin", pPtr);
   optFmin = pPtr.dbleData_;
   psuadeIO->getParameter("ana_opt_numfmin", pPtr);
   optNumFmin = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_maxfeval", pPtr);
   odata.maxFEval_ = pPtr.intData_;
   psuadeIO->getParameter("ana_opt_deltax", pPtr);
   odata.deltaX_ = pPtr.dbleData_;
   psuadeIO->getParameter("ana_opt_smtargetfile", pPtr);
   strcpy(specsFile, pPtr.strArray_[0]);
   optInitXData = optData[0];
   optInitYData = optData[1];
   optimalXData = optData[2];
   optimalYData = optData[3];
   optimalCount = optCount[0];
   psuadeIO->getParameter("app_maxparalleljobs", pPtr);
   odata.maxParallelJobs_ = pPtr.intData_;
   if (optimizeNumPts <= 0)
   {
      printf("PSUADE Optimizer: number of desired minimum <= 0.\n");
      return 0;
   }

   if ((optimizeFlag != 0) && (optimizeUseRS == 0) && optMethod != 10)
   {
      printAsterisks(PL_INFO, 0);
      printf("PSUADE : Search for local minima in sample space.\n");
      printf("         number of samples = %d\n", nSamples);

      tmpInitX = optInitXData;
      tmpInitY = optInitYData;
      tmpOptX  = optimalXData;
      tmpOptY  = optimalYData;
      optimalCount += optimizeNumPts;
      if (optimizeNumPts > 0) 
      {
         optInitXData = new double[optimalCount*nInputs];
         optInitYData = new double[optimalCount];
         optimalXData = new double[optimalCount*nInputs];
         optimalYData = new double[optimalCount];
         for (ind1 = 0; ind1 < optimalCount; ind1++)
            optimalYData[ind1] = PSUADE_UNDEFINED;
      }
      for (ind1 = 0; ind1 < optimalCount-optimizeNumPts; ind1++)
      {
         for (ii = 0; ii < nInputs; ii++)
            optInitXData[ind1*nInputs+ii] = tmpInitX[ind1*nInputs+ii]; 
         optInitYData[ind1] = tmpInitY[ind1]; 
         for (ii = 0; ii < nInputs; ii++)
            optimalXData[ind1*nInputs+ii] = tmpOptX[ind1*nInputs+ii]; 
         optimalYData[ind1] = tmpOptY[ind1]; 
      }

      tempOutputs = new double[nSamples];
      dbleInd = new double[nSamples];
      for (ss = 0; ss < nSamples; ss++)
      {
         tempOutputs[ss] = sampleOutputs[ss*nOutputs+optID];
         dbleInd[ss] = (double) ss;
      }
      if (nSamples > optimizeNumPts)
         sortDbleList2(nSamples, tempOutputs, dbleInd);

      optNewCnt = optimalCount - optimizeNumPts;
      ind1 = 0;
      while (optNewCnt < optimalCount && ind1 < nSamples)
      {
         minIndex = (int) dbleInd[ind1];
         for (ind2 = 0; ind2 < optNewCnt; ind2++)
         {
            for (ii = 0; ii < nInputs; ii++)
               if (sampleInputs[minIndex*nInputs+ii] != 
                  optInitXData[ind2*nInputs+ii]) break;
            if (ii == nInputs) break;
         }
         if (ind2 == optNewCnt)
         {
            for (ii = 0; ii < nInputs; ii++)
               optInitXData[optNewCnt*nInputs+ii] =
                  sampleInputs[minIndex*nInputs+ii]; 
            optInitYData[optNewCnt] = tempOutputs[ind1]; 
            optNewCnt++;
         }
         ind1++;
      }
      optimalCount = optNewCnt;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         printf("\t optimization starting point %d : \n", ind1+1);
         for (ii = 0; ii < nInputs; ii++)
            printf("\t\t %16.8e\n", optInitXData[ind1*nInputs+ii]);
         printf("\t\t\t\t Y = %16.8e\n", optInitYData[ind1]);
      }
      printAsterisks(PL_INFO, 0);
      delete [] tempOutputs;
      delete [] dbleInd;
      if (tmpInitX != NULL) delete [] tmpInitX;
      if (tmpInitY != NULL) delete [] tmpInitY;
      if (tmpOptX  != NULL) delete [] tmpOptX;
      if (tmpOptY  != NULL) delete [] tmpOptY;
      optData[0] = optInitXData;
      optData[1] = optInitYData;
      optData[2] = optimalXData;
      optData[3] = optimalYData;
   }

   else if ((optimizeFlag != 0) && (optimizeUseRS == 1))
   {
      if      (nInputs <=  3) nPtsPerDim = 32;
      else if (nInputs ==  4) nPtsPerDim = 24;
      else if (nInputs ==  5) nPtsPerDim = 16;
      else if (nInputs ==  6) nPtsPerDim = 10;
      else if (nInputs ==  7) nPtsPerDim =  8;
      else if (nInputs ==  8) nPtsPerDim =  6;
      else if (nInputs ==  9) nPtsPerDim =  5;
      else if (nInputs == 10) nPtsPerDim =  4;
      else if (nInputs == 11) nPtsPerDim =  3;
      else if (nInputs == 12) nPtsPerDim =  3;
      else if (nInputs == 13) nPtsPerDim =  3;
      else if (nInputs >= 14) nPtsPerDim =  2;
      faPtr = genFAInteractive(psuadeIO, 0);
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(iLowerB, iUpperB);
      tempOutputs = new double[nSamples];
      for (ss = 0; ss < nSamples; ss++) 
         tempOutputs[ss] = sampleOutputs[ss*nOutputs+optID];
      faPtr->genNDGridData(sampleInputs,tempOutputs, &faLeng,
                           &faXOut,&faYOut);
      dimIndices = new int[nInputs];
      dimFactors = new int[nInputs+1];
      tmpInitX = optInitXData;
      tmpInitY = optInitYData;
      tmpOptX  = optimalXData;
      tmpOptY  = optimalYData;
      optimalCount += optimizeNumPts;
      if (optimalCount > 0)
      {
         optInitXData = new double[optimalCount*nInputs];
         optInitYData = new double[optimalCount];
         optimalXData = new double[optimalCount*nInputs];
         optimalYData = new double[optimalCount];
         for (ind1 = 0; ind1 < optimalCount; ind1++)
            optimalYData[ind1] = PSUADE_UNDEFINED;
      }
      for (ind1 = 0; ind1 < optimalCount-optimizeNumPts; ind1++)
      {
         for (ii = 0; ii < nInputs; ii++)
            optInitXData[ind1*nInputs+ii] = tmpInitX[ind1*nInputs+ii]; 
         optInitYData[ind1] = tmpInitY[ind1]; 
         for (ii = 0; ii < nInputs; ii++)
            optimalXData[ind1*nInputs+ii] = tmpOptX[ind1*nInputs+ii]; 
         optimalYData[ind1] = tmpOptY[ind1]; 
      }

      printAsterisks(PL_INFO, 0);
      printf("PSUADE : Search for local minima in RS space.\n");
      printf("         number of samples = %d\n", nSamples);
      optNewCnt = optimalCount - optimizeNumPts;
      for (ind1 = 0; ind1 < faLeng; ind1++)
      {
         ind2 = ind1;
         dimFactors[0] = 1;
         for (ii = 0; ii < nInputs; ii++)
         {
            dimIndices[ii] = ind2 % nPtsPerDim;
            dimFactors[ii+1] = dimFactors[ii] * nPtsPerDim;
            ind2 = ind2 / nPtsPerDim;
         }
         minCheck = 1;
         for (ii = 0; ii < nInputs; ii++)
         {
            if ((dimIndices[ii]>0) &&  
                (faYOut[ind1]>=faYOut[ind1-dimFactors[ii]]))
               minCheck = 0;
            if ((dimIndices[ii]<nPtsPerDim-1) && 
                (faYOut[ind1]>=faYOut[ind1+dimFactors[ii]]))
               minCheck = 0;
         }
         if (minCheck == 1)
         {
            for (ii = 0; ii < nInputs; ii++)
               optInitXData[optNewCnt*nInputs+ii] = faXOut[ind1*nInputs+ii];  
            optInitYData[optNewCnt] = faYOut[ind1];  
            printf("PSUADE: Local minimum at X = \n");
            for (ii = 0; ii < nInputs; ii++)
               printf("\t\t %16.8e\n", faXOut[ind1*nInputs+ii]);  
            printf("\tY = %16.8e\n", faYOut[ind1]);  
            optNewCnt++;
         }
         if (optNewCnt >= optimalCount) break;
      }
      printAsterisks(PL_INFO, 0);
      delete [] tmpInitX;
      delete [] tmpInitY;
      delete [] tmpOptX;
      delete [] tmpOptY;
      delete [] dimIndices;
      delete [] dimFactors;
      delete [] tempOutputs;
      delete [] faXOut;
      delete [] faYOut;
      delete faPtr;
   }


   strcpy(odata.targetFile_, specsFile);
   odata.psIO_ = psuadeIO;
                                                                                
   if ((optimizeFlag == 1) && (optMethod == 1))
   {
#ifdef HAVE_TXMATH
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      TxMathPtr = new TxMathOptimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if (optInitYData[ind1] > optCutOff && ind1 > 0)
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optimalYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);
            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.optimalX_ = new double[nInputs];
            odata.tolerance_ = optTol;
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            TxMathPtr->optimize(&odata);
            printf("\t TxMath number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            for (ii = 0; ii < nInputs; ii++)
               printf("\t optimum  X(%6d) = %16.8e\n",ii+1,
                      odata.optimalX_[ii]);
            printf("\t\t\t optimum Y = %16.8e\n",odata.optimalY_);
            for (ii = 0; ii < nInputs; ii++)
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
            optimalYData[ind1] = odata.optimalY_;
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete TxMathPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
#else
      printf("TXMATH not linked.\n");
      exit(1);
#endif
   }

   else if ((optimizeFlag == 1) && (optMethod == 2))
   {
#ifdef HAVE_APPSPACK
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      APPSPACKPtr = new APPSPACKOptimizer();
      odata.optimalY_ = 1e35;
      odata.numFuncEvals_ = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);
            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.nSamples_ = nSamples;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            APPSPACKPtr->optimize(&odata);
            printf("\t APPSPACK number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("\t optimum  X(%6d) = %16.8e\n",ii+1,
                      odata.optimalX_[ii]);
            printf("\t\t\t optimum Y = %16.8e\n",odata.optimalY_);
            for (ii = 0; ii < nInputs; ii++)
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
            optimalYData[ind1] = odata.optimalY_;
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete APPSPACKPtr;
      delete [] tempInputs;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
#else
      printf("APPSPACK not linked.\n");
      exit(1);
#endif
   }

   else if ((optimizeFlag == 1) && (optMethod == 3))
   {
#ifdef HAVE_MINPACK
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      MinpackPtr = new MinpackOptimizer();
      odata.optimalY_ = 1e35;
      odata.numFuncEvals_ = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);
            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            MinpackPtr->optimize(&odata);
            printf("\t Minpack number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            printf("\t\t\t optimum Y = %16.8e\n",odata.optimalY_);
            for (ii = 0; ii < nInputs; ii++)
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
            optimalYData[ind1] = odata.optimalY_;
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete MinpackPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
#else
      printf("Minpack not linked.\n");
      exit(1);
#endif
   }

   else if ((optimizeFlag == 1) && (optMethod == 4))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      CobylaPtr = new CobylaOptimizer();
      odata.optimalY_ = 1e35;
      odata.numFuncEvals_ = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            CobylaPtr->optimize(&odata);
            printf("\t Cobyla number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete CobylaPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
   }

   else if ((optimizeFlag == 1) && (optMethod == 5))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      SMPtr = new SMOptimizer();
      odata.optimalY_ = 1e35;
      odata.numFuncEvals_ = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);
                                                                                
         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n",
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);
                                                                                
            odata.initialX_ = &(optInitXData[ind1*nInputs]);
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            SMPtr->optimize(&odata);
            printf("\t SM number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete SMPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0)
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
   }

   else if ((optimizeFlag == 1) && (optMethod == 6 || optMethod == 7))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      MMPtr = new MMOptimizer();
      odata.optimalY_ = 1e35;
      odata.numFuncEvals_ = 0;
      totalEval = 0;
      int totalCEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);
                                                                                
         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n",
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);
                                                                                
            odata.initialX_ = &(optInitXData[ind1*nInputs]);
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            if (optMethod == 7)
            {
               strcpy(sparam,"setAdaptive");
               MMPtr->setParam(sparam);
            } 
            MMPtr->optimize(&odata);
            printf("\t MM number of fine   function evaluations = %d\n",
                   odata.intData_);
            printf("\t MM number of coarse function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalCEval += odata.numFuncEvals_;
            totalEval  += odata.intData_;
            for (ii = 0; ii < nInputs; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      printf("\t MM total number of fine   function evaluations = %d\n",
             totalEval);
      printf("\t MM total number of coarse function evaluations = %d\n",
             totalCEval);
      delete MMPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0)
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
   }

   else if ((optimizeFlag == 1) && (optMethod == 8))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      BobyqaPtr = new BobyqaOptimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            BobyqaPtr->optimize(&odata);
            printf("\t Bobyqa number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            for (ii = 0; ii < nInputs; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete BobyqaPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t Bobyqa total number of function evaluations = %d\n",totalEval);
   }

   else if ((optimizeFlag == 1) && (optMethod == 9))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      SCEPtr = new SCEOptimizer();
      odata.optimalY_ = 1e35;
      odata.numFuncEvals_ = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         if (psInteractive_ == 1)
         {
           printAsterisks(PL_INFO, 0);
           printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);
         }
         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            if (psInteractive_ == 1)
            {
              for (ii = 0; ii < nInputs; ii++)
                 printf("\t starting X(%6d) = %16.8e\n",ii+1,
                        optInitXData[ind1*nInputs+ii]);
              printf("\t starting Y = %16.8e\n",optInitYData[ind1]);
            }
            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            SCEPtr->optimize(&odata);
            if (psInteractive_ == 1)
              printf("\t SCE number of function evaluations = %d\n",
                     odata.numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
            {
              optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
              if (psInteractive_ == 1)
                 printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            if (psInteractive_ == 1)
              printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete SCEPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
   }

   else if ((optimizeFlag == 1) && (optMethod == 10))
   {
      printAsterisks(PL_INFO, 0);
      printf("PSUADE OPTIMIZATION: \n");
      for (ii = 0; ii < nInputs; ii++)
         printf("\t starting X(%6d) = %16.8e\n",ii+1,optInitXData[ii]);
      printf("\t starting Y = %16.8e\n",optInitYData[0]);

      odata.initialX_ = optInitXData; 
      funcIO->setSynchronousMode();
      odata.funcIO_ = funcIO;
      odata.nInputs_ = nInputs;
      odata.nOutputs_ = nOutputs;
      odata.outputID_ = optID;
      odata.lowerBounds_ = iLowerB;
      odata.upperBounds_ = iUpperB;
      odata.outputLevel_ = optPrintLevel;
      odata.tolerance_   = optTol;
      odata.optimalX_ = new double[nInputs];
      MooPtr = new MultiObjectiveOptimizer();
      MooPtr->optimize(&odata);
      delete [] odata.optimalX_;
      odata.optimalX_ = NULL;
      odata.initialX_ = NULL;
      odata.lowerBounds_ = NULL;
      odata.upperBounds_ = NULL;
      odata.funcIO_ = NULL;
      delete MooPtr;
      printAsterisks(PL_INFO, 0);
   }

   else if ((optimizeFlag == 1) && (optMethod == 12))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      OUU1Ptr = new OUU1Optimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            odata.intData_ = -1;
            OUU1Ptr->optimize(&odata);
            printf("\t OUU1Optimizer number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            M1 = odata.intData_;
            if (M1 <= 0) M1 = nInputs;
            for (ii = 0; ii < M1; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete OUU1Ptr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t OUU1 total number of function evaluations = %d\n",
             totalEval);
   }

   else if ((optimizeFlag == 1) && (optMethod == 13))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      OUU2Ptr = new OUU2Optimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            odata.intData_ = -1;
            OUU2Ptr->optimize(&odata);
            printf("\t OUU2Optimizer number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            M1 = odata.intData_;
            if (M1 <= 0) M1 = nInputs;
            for (ii = 0; ii < M1; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete OUU2Ptr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t OUU2 total number of function evaluations = %d\n",
             totalEval);
   }
 
   else if ((optimizeFlag == 1) && ((optMethod == 11) || 
            (optMethod == 16) || (optMethod == 17) || 
            (optMethod == 18) || (optMethod == 20)))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      OUUPtr = new OUUOptimizer();
      if (optMethod == 11) OUUPtr->setLocalOptimizer(0);
      if (optMethod == 16) OUUPtr->setLocalOptimizer(1);
      if (optMethod == 17) OUUPtr->setLocalOptimizer(0);
      if (optMethod == 18) OUUPtr->setLocalOptimizer(2);
      if (optMethod == 20) OUUPtr->setLocalOptimizer(3);
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            odata.intData_ = -1;
            OUUPtr->optimize(&odata);
            printf("\t OUUOptimizer number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            M1 = odata.intData_;
            if (M1 <= 0) M1 = nInputs;
            for (ii = 0; ii < M1; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete OUUPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t OUU total number of function evaluations = %d\n",
             totalEval);
   }

   else if ((optimizeFlag == 1) && (optMethod == 14))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      LincoaPtr = new LincoaOptimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            LincoaPtr->optimize(&odata);
            printf("\t Lincoa number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            for (ii = 0; ii < nInputs; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete LincoaPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t Lincoa total number of function evaluations = %d\n",
             totalEval);
   }

   else if ((optimizeFlag == 1) && (optMethod == 15))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      NewuoaPtr = new NewuoaOptimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            NewuoaPtr->optimize(&odata);
            printf("\t Newuoa number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            for (ii = 0; ii < nInputs; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete NewuoaPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t Newuoa total number of function evaluations = %d\n",
             totalEval);
   }

   else if ((optimizeFlag == 1) && (optMethod == 19))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      lbfgsPtr = new LBFGSOptimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               printf("\t starting X(%6d) = %16.8e\n",ii+1,
                      optInitXData[ind1*nInputs+ii]);
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            lbfgsPtr->optimize(&odata);
            printf("\t LBFGS number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            for (ii = 0; ii < nInputs; ii++)
            {
               optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
               printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                      odata.optimalX_[ii]);
            }
            optimalYData[ind1] = odata.optimalY_;
            printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete lbfgsPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t LBFGS total number of function evaluations = %d\n",
             totalEval);
   }

   else if ((optimizeFlag == 1) && (optMethod == 21))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      nomadPtr = new NomadOptimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
            {
               ss = (int) optInitXData[ind1*nInputs+ii];
               if (optInitXData[ind1*nInputs+ii] - ss == 0)
                  printf("\t starting X(%6d) = %d\n",ii+1,ss);
               else
                  printf("\t starting X(%6d) = %16.8e\n",ii+1,
                         optInitXData[ind1*nInputs+ii]);
            }
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            nomadPtr->optimize(&odata);
            printf("\t Nomad number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            if (odata.optimalY_ != 1.0e50)
            {
               for (ii = 0; ii < nInputs; ii++)
               {
                  optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
                  ss = (int) optimalXData[ind1*nInputs+ii];
                  if (optimalXData[ind1*nInputs+ii] - ss == 0)
                     printf("\t optimum  X(%6d) = %d\n", ii+1, ss);
                  else
                     printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                         odata.optimalX_[ii]);
               }
               optimalYData[ind1] = odata.optimalY_;
               printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            }
            else
            {
               printf("** PSUADE NOMAD INFO: no feasible solution found.\n");
            }
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete nomadPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t Nomad total number of function evaluations = %d\n",
             totalEval);
   }

   else if ((optimizeFlag == 1) && (optMethod == 22))
   {
      numYmin = 0;
      istart = optimalCount - optimizeNumPts;
      if (istart < 0) istart = 0;
      minlpPtr = new OUUMinlpOptimizer();
      odata.optimalY_ = 1e35;
      totalEval = 0;
      for (ind1 = istart; ind1 < optimalCount; ind1++)
      {
         printAsterisks(PL_INFO, 0);
         printf("PSUADE OPTIMIZATION %d (%d) : \n", ind1+1, optimalCount);

         if ((optInitYData[ind1] > optCutOff) && (ind1 > 0))
         {
            printf("skip optimization (%16.8e > %16.8e) : \n", 
                   optInitYData[ind1], optCutOff);
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
            {
               ss = (int) optInitXData[ind1*nInputs+ii];
               if (optInitXData[ind1*nInputs+ii] - ss == 0)
                  printf("\t starting X(%6d) = %d\n",ii+1,ss);
               else
                  printf("\t starting X(%6d) = %16.8e\n",ii+1,
                         optInitXData[ind1*nInputs+ii]);
            }
            printf("\t starting Y = %16.8e\n",optInitYData[ind1]);

            odata.initialX_ = &(optInitXData[ind1*nInputs]); 
            funcIO->setSynchronousMode();
            odata.funcIO_ = funcIO;
            odata.nInputs_ = nInputs;
            odata.nOutputs_ = nOutputs;
            odata.outputID_ = optID;
            odata.lowerBounds_ = iLowerB;
            odata.upperBounds_ = iUpperB;
            odata.outputLevel_ = optPrintLevel;
            odata.tolerance_   = optTol;
            odata.optimalX_ = new double[nInputs];
            odata.numFuncEvals_ = 0;
            if      (optimalCount-istart == 1) odata.setOptDriver_ = 3;
            else if (ind1 == istart)           odata.setOptDriver_ = 1;
            else if (ind1 == optimalCount-1)   odata.setOptDriver_ = 2;
            else                               odata.setOptDriver_ = 0;
            minlpPtr->optimize(&odata);
            printf("\t OUU/MINLP number of function evaluations = %d\n",
                   odata.numFuncEvals_);
            totalEval += odata.numFuncEvals_;
            if (odata.optimalY_ != 1.0e50)
            {
               for (ii = 0; ii < nInputs; ii++)
               {
                  optimalXData[ind1*nInputs+ii] = odata.optimalX_[ii];
                  ss = (int) optimalXData[ind1*nInputs+ii];
                  if (optimalXData[ind1*nInputs+ii] - ss == 0)
                     printf("\t optimum  X(%6d) = %d\n", ii+1, ss);
                  else
                     printf("\t optimum  X(%6d) = %16.8e\n", ii+1,
                         odata.optimalX_[ii]);
               }
               optimalYData[ind1] = odata.optimalY_;
               printf("\t\t\t optimum Y = %16.8e\n", odata.optimalY_);
            }
            else
            {
               printf("** PSUADE NOMAD INFO: no feasible solution found.\n");
            }
            delete [] odata.optimalX_;
            odata.optimalX_ = NULL;
            odata.initialX_ = NULL;
            odata.lowerBounds_ = NULL;
            odata.upperBounds_ = NULL;
            odata.funcIO_ = NULL;
            if (optimalYData[ind1] <= optFmin) numYmin++;
         }
         printAsterisks(PL_INFO, 0);
      }
      delete minlpPtr;
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optFmin != 0.0) 
         {
            if (optimalYData[ind1] <= optFmin) returnFlag++;
         }
         else
         {
            if (PABS(optimalYData[ind1]) <= 1.0e-6) returnFlag++;
         }
      }
      printf("\t OUU/MINLP total number of function evaluations = %d\n",
             totalEval);
   }

   optCount[0] = optimalCount;
   int nTrack=0, nn, count;
   int *trackIndices = new int[optNumFmin];
   double minY, range, *trackYmin = new double[optNumFmin];
   if ((optimalCount > 0) && (optMethod != 0) && (optMethod != 10))
   {
      ind2 = 0;
      Ymin = 1e36;
      if (psInteractive_ != 0)
      {
        printAsterisks(PL_INFO, 0);
        printf("PSUADE Optimization Results.\n");
      }
      for (ind1 = 0; ind1 < optimalCount; ind1++)
      {
         if (optimalYData[ind1] != PSUADE_UNDEFINED)
         {
            if (psInteractive_ != 0)
              printf("PSUADE Optimization : local optima %d (%d) - \n",ind1+1,
                     optimalCount);
            ind2 = M1;
            if (ind2 < 0) ind2 = nInputs;
            for (ii = 0; ii < ind2; ii++)
            {
               ss = (int) optimalXData[ind1*nInputs+ii];
               if (psInteractive_ != 0)
               {
                 if (optimalXData[ind1*nInputs+ii] - ss == 0)
                    printf("\t\tX %5d = %d\n", ii+1,ss);
                 else 
                    printf("\t\tX %5d = %16.8e\n", ii+1,
                           optimalXData[ind1*nInputs+ii]);
               }
            }
            if (psInteractive_ != 0)
              printf("\t\t\tYmin = %16.8e\n",optimalYData[ind1]);
            if (nTrack < optNumFmin) 
            {
               count = 0;
               for (nn = 0; nn < nTrack; nn++)
               {
                  ind2 = trackIndices[nn];
                  count = 0;
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     range = iUpperB[ii] - iLowerB[ii];
                     if (range == 0) range = 1;
                     if (PABS(optimalXData[ind2*nInputs+ii]-
                              optimalXData[ind1*nInputs+ii])/range<1e-8) count++;
                  }
                  if (count == nInputs) break;
               }
               if (count < nInputs)
               {
                  trackIndices[nTrack] = ind1;
                  trackYmin[nTrack++] = optimalYData[ind1];
               }
            }
            else
            {
               ind2 = 0;
               minY = trackYmin[0];
               for (nn = 1; nn < nTrack; nn++)
               {
                  if (trackYmin[nn] > minY) 
                  {
                     ind2 = nn;
                     minY = trackYmin[nn];
                  }
               }
               if (optimalYData[ind1] < minY)
               {
                  count = 0;
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     range = iUpperB[ii] - iLowerB[ii];
                     if (range == 0) range = 1;
                     if (PABS(optimalXData[ind2*nInputs+ii]-
                              optimalXData[ind1*nInputs+ii])/range<1e-8) count++;
                  }
                  if (count < nInputs)
                  {
                     trackIndices[ind2] = ind1;
                     trackYmin[ind2] = optimalYData[ind1];
                  }
               }
            }
         }
      } 
      if (optimalCount > 0)
      {
         if (psInteractive_ != 0)
         {
           printf("##################################################\n");
           printf("PSUADE OPTIMIZATION : CURRENT GLOBAL MINIMUM - \n");
           for (nn = 0; nn < nTrack; nn++)
           {
              if (M1 < 0) M1 = nInputs;
              ind2 = trackIndices[nn];
              for (ii = 0; ii < M1; ii++)
              {
                 ind1 = (int) optimalXData[ind2*nInputs+ii];
                 if (optimalXData[ind2*nInputs+ii] - ind1 == 0)
                    printf("\t\tX %5d = %d\n",ii+1,ind1);
                 else
                    printf("\t\tX %5d = %16.8e\n",ii+1,
                           optimalXData[ind2*nInputs+ii]);
              }
              printf("\t\t\tYmin = %16.8e\n",optimalYData[ind2]);
              if (nn < nTrack-1)
                 printf("--------------------------------------------------\n");
           }
           printf("##################################################\n");
           printAsterisks(PL_INFO, 0);
         }
      }
   }
   delete [] trackIndices;
   delete [] trackYmin;
   if (optimizeFlag == 1 && returnFlag >= optNumFmin && optMethod != 10) 
   {
      printf("PSUADE Optimization : desired minimum found.\n");
      return returnFlag;
   }
   else returnFlag = 0;
   return returnFlag;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
Optimizer& Optimizer::operator=(const Optimizer &)
{
   printf("Optimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

