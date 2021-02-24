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
// Functions for the class PsuadeBase
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
extern "C" {
   void dsyev_(char *, char *, int *, double *, int *, double *,
               double *, int *, int *);
}

// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ------------------------------------------------------------------------
// local includes : class definition and utilities
// ------------------------------------------------------------------------
#include "Psuade.h"
#include "PsuadeBase.h"
#include "dtype.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PDFManager.h"
#include "JobCntl.h"
#include "PrintingTS.h"

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// run the special adaptive mode with gradient analyzer
// ------------------------------------------------------------------------
int PsuadeBase::runAdaptiveGradBased()
{
   int    nInputs, nOutputs, nSamples, ii, iR, iS, kk, iRov, nRefineCurr;
   int    *sampleROVInds, nROVMax, *states, nRefines, *refineSeps;
   int    maxParallelJobs, lwork, info, Ysize, count, status, flag;
   int    rovID, rovCnt, initNSamples, sampleIndex;
   int    refineSize, nRefineLocal;
   double *iLowerB, *iUpperB, *sampleInputs, *ranges, *Ystore, Ynew;
   double *sampleROVVals, *sampleROVErrs, error, anaThreshold, threshold;
   double *ROVStoreVals, *ROVStoreGrads, *ROVStorePts, *ROVStoreEigens;
   double *currentPt, YVal, *gradients, *hessian, errMax, dtemp;
   double diff, eigMax, dsum2, *work, *eigs, stdev, mean;
   double tolerance, tol2, eig2, threshDec, minEigen=0.0;
   char   jobz, uplo, **outNames, lineIn[500], fileName[500];
   char   **inNames;
   FILE   *fileIn;
   PsuadeData        *psIO;
   FunctionInterface *funcIO;
   pData             pPtr, pLowerB, pUpperB;
   pData             pInputs, pOutputs, pInputs2, pOutputs2;

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, 
           "PSUADE adaptiveGradBased: adaptive sampling based on using\n");
   printOutTS(PL_INFO, "       gradients to form regions of validity.\n");
   printOutTS(PL_INFO, "       (METIS sampling recommendated)\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO, "Note: Turn on rs_expert mode to set RS parameters.\n");
   printOutTS(PL_INFO, "      Turn on outputLevel (>0) to get error histogram.\n");
   printEquals(PL_INFO, 0);

   psuadeIO_->getParameter("input_ninputs", pPtr);
   nInputs  = pPtr.intData_;
   psuadeIO_->getParameter("output_noutputs", pPtr);
   nOutputs  = pPtr.intData_;
   if (nInputs > 12)
   {
      printOutTS(PL_ERROR, 
           "PSUADE adaptiveGradBased ERROR: nInputs should be <= 12.\n");
      exit(1);
   }
   if (nOutputs != 1)
   {
      printOutTS(PL_ERROR, 
           "PSUADE adaptiveGradBased ERROR: nOutputs should be 1.\n");
      exit(1);
   }

   psuadeIO_->getParameter("method_nsamples", pPtr);
   nSamples = pPtr.intData_;
   psuadeIO_->getParameter("method_nrefinements", pPtr);
   nRefines = pPtr.intData_;
   psuadeIO_->getParameter("method_refine_size", pPtr);
   refineSize = pPtr.intData_;
   nROVMax = nSamples + refineSize * nRefines;
   printOutTS(PL_INFO, 
          "PSUADE adaptiveGradBased: max no. of sample points = %d\n",
          nROVMax);
   printOutTS(PL_INFO, 
          "PSUADE adaptiveGradBased: no. of refinements (set at 4)   = %d\n",
          nRefines);

   psuadeIO_->getParameter("app_maxparalleljobs", pPtr);
   maxParallelJobs = pPtr.intData_;
   psuadeIO_->getParameter("ana_threshold", pPtr);
   anaThreshold = pPtr.dbleData_;
   printOutTS(PL_INFO, "PSUADE adaptiveGradBased: analysis threshold              = %e\n",
          anaThreshold);
   if (anaThreshold > 0.5)
   {
      printOutTS(PL_INFO, 
           "PSUADE adaptiveGradBased INFO: threshold must be <= 0.5.\n");
      printOutTS(PL_INFO, 
           "                               threshold reset to 0.5.\n");
      anaThreshold = 0.5;
   }
   psuadeIO_->getParameter("input_lbounds", pLowerB);
   iLowerB = pLowerB.dbleArray_;
   psuadeIO_->getParameter("input_ubounds", pUpperB);
   iUpperB = pUpperB.dbleArray_;
   Ysize = (nInputs + 1) * nInputs / 2 + nInputs + 2;
   currentPt = new double[nInputs];

   nSamples = 0;
   threshold = anaThreshold;
   if (psAnaExpertMode_ == 1)
   {
      printOutTS(PL_INFO, "The current cutoff threshold (2nd order) is %e.\n",
             anaThreshold);
      //threshold = 0.0;
      //while (threshold <= 0.0 || threshold < anaThreshold)
      //{
      //   printOutTS(PL_INFO, "For progressive thresholding, enter initial ");
      //   printOutTS(PL_INFO, "threshold (%e to 0.5): ", anaThreshold);
      //   scanf("%lg", &threshold);
      //}
      threshold = anaThreshold;
      if (threshold != anaThreshold)
      {
         threshDec = 0.0;
         while (threshDec <= 0.0 || threshDec >= 1.0)
         {
            printf("Enter a threshold decrement factor (> 0, < 1): ");
            scanf("%lg", &threshDec);
         }
      }
      printOutTS(PL_INFO, 
           "Max eigenvalue of Hessian controls extent of ROV.\n");
      printOutTS(PL_INFO, 
           "You can set the minimum value to control the extent.\n");
      minEigen = -1;
      while (minEigen < 0.0)
      {
         printf("Enter the minimum eigenvalue (0 if no minimum) : ");
         scanf("%lg", &minEigen);
      }
      printOutTS(PL_INFO, 
           "Do you have a partial PsuadeGradRSM.rov file? (y or n) ");
      scanf("%s", lineIn);
      if (lineIn[0] == 'y')
      {
         printAsterisks(PL_INFO, 0);
         printf("Enter ROV file name : ");
         scanf("%s", fileName);
         fileIn = fopen(fileName, "r");
         if (fileIn == NULL)
         {
            printOutTS(PL_INFO, "File name %s found.\n", fileName);
            return -1;
         }
         fclose(fileIn);
         psIO = new PsuadeData();
         status = psIO->readPsuadeFile(fileName);
         if (status != 0)
         {
            printOutTS(PL_ERROR, "ERROR: Problem reading file %s.\n", 
                       fileName);
            return -1;
         }
         psIO->getParameter("input_ninputs", pPtr);
         if (pPtr.intData_ != nInputs)
         {
            printOutTS(PL_INFO, 
                "nInputs in file %s does not match with current file.\n",
                   fileName);
            return -1;
         }
         psIO->getParameter("output_noutputs", pPtr);
         if (pPtr.intData_ != Ysize)
         {
            printOutTS(PL_INFO, 
                "nOutputs in file %s (%d) is not valid (should be %d).\n",
                   fileName, pPtr.intData_, Ysize);
            return -1;
         }
         psIO->getParameter("method_nsamples", pPtr);
         iRov = pPtr.intData_;
         if (iRov > nROVMax || iRov <= 0)
         {
            printOutTS(PL_ERROR, 
                "ERROR: nROVs in file %s too large (should be <= %d)\n",
                   fileName, nROVMax);
            return -1;
         }
         psIO->getParameter("input_sample", pInputs);
         ROVStorePts = new double[nROVMax*nInputs];
         for (iR = 0; iR < iRov*nInputs; iR++)
            ROVStorePts[iR] = pInputs.dbleArray_[iR];
         psIO->getParameter("output_sample", pOutputs);
         ROVStoreVals = new double[nROVMax];
         ROVStoreGrads = new double[nROVMax*nInputs];
         for (iR = 0; iR < iRov; iR++)
         {
            ROVStoreVals[iR] = pOutputs.dbleArray_[iR*Ysize];
            for (ii = 0; ii < nInputs; ii++)
               ROVStoreGrads[iR*nInputs+ii] =
                                    pOutputs.dbleArray_[iR*Ysize+ii+1];
         }
         ROVStoreEigens = new double[nROVMax];
         for (iR = 0; iR < iRov; iR++)
            ROVStoreEigens[iR] = pOutputs.dbleArray_[(iR+1)*Ysize-1];
         Ystore = new double[nROVMax*Ysize];
         for (iR = 0; iR < iRov*Ysize; iR++)
            Ystore[iR] = pOutputs.dbleArray_[iR];
         delete psIO;

         printAsterisks(PL_INFO, 0);
         printf("Enter the corresponding PsuadeGradRSM.sample file name : ");
         scanf("%s", fileName);
         fileIn = fopen(fileName, "r");
         if (fileIn == NULL)
         {
            printOutTS(PL_INFO, "Sample file name %s found.\n", fileName);
            return -1;
         }
         fclose(fileIn);
         psIO = new PsuadeData();
         status = psIO->readPsuadeFile(fileName);
         if (status != 0)
         {
            printOutTS(PL_ERROR, 
                "ERROR: Problem reading sample file %s.\n", fileName);
            return -1;
         }
         psIO->getParameter("input_ninputs", pPtr);
         if (pPtr.intData_ != nInputs)
         {
            printOutTS(PL_INFO, 
                "nInputs in file %s does not match with current file.\n",
                   fileName);
            return -1;
         }
         psIO->getParameter("output_noutputs", pPtr);
         if (pPtr.intData_ != nOutputs)
         {
            printOutTS(PL_INFO, 
                "nOutputs in file %s does not match with current file.\n",
                   fileName);
            return -1;
         }
         psIO->getParameter("method_nsamples", pPtr);
         nSamples = pPtr.intData_;
         if (nSamples < iRov)
         {
            printOutTS(PL_ERROR, 
                "ERROR: sample in file %s too small (should be > %d)\n",
                   fileName, iRov);
            return -1;
         }
         psIO->getParameter("input_sample", pInputs2);
         sampleInputs = new double[nSamples*nInputs];
         for (iR = 0; iR < nSamples*nInputs; iR++)
            sampleInputs[iR] = pInputs2.dbleArray_[iR];
         psIO->getParameter("output_sample", pOutputs2);
         sampleROVVals = new double[nSamples];
         for (iR = 0; iR < nSamples; iR++)
            sampleROVVals[iR] = pOutputs2.dbleArray_[iR];
         sampleROVErrs = new double[nSamples];
         for (iR = 0; iR < nSamples; iR++)
            sampleROVErrs[iR] = PSUADE_UNDEFINED;
         sampleROVInds = new int[nSamples];
         for (iR = 0; iR < nSamples; iR++) sampleROVInds[iR] = -1;
         delete psIO;

         for (iR = 0; iR < iRov; iR++)
         {
            YVal = ROVStoreVals[iR];
            eigMax = ROVStoreEigens[iR];
            for (iS = 0; iS < nSamples; iS++)
            {
               Ynew = YVal;
               dsum2 = 0.0;
               for (ii = 0; ii < nInputs; ii++)
               {
                  diff = sampleInputs[iS*nInputs+ii] - 
                         ROVStorePts[iR*nInputs+ii]; 
                  Ynew += ROVStoreGrads[iR*nInputs+ii] * diff;
                  dsum2 += diff * diff;
               }
               error = 0.5 * PABS(eigMax) * dsum2;
               if (error < sampleROVErrs[iS])
               {
                  sampleROVInds[iS] = iR;
                  sampleROVVals[iS] = Ynew;
                  sampleROVErrs[iS] = error;
               }
               if (error > anaThreshold) sampleROVInds[iS] = -1;
            }
         }

         error = sampleROVErrs[0];
         kk = 0;
         for (iS = 1; iS < nSamples; iS++)
         {
            if (sampleROVErrs[iS] > error)
            {
               kk = iS;
               error = sampleROVErrs[iS];
            }
         }
         for (ii = 0; ii < nInputs; ii++)
            currentPt[ii] = sampleInputs[kk*nInputs+ii];
         sampleIndex = -1;
         if (sampleROVInds[kk] < 0) sampleIndex = kk;
         if (sampleIndex == -1)
         {
            printOutTS(PL_INFO, 
              "PSUADE adaptiveGradBased: the ROV sample is good to go.\n");
            delete [] sampleInputs;
            delete [] sampleROVInds;
            delete [] sampleROVVals;
            delete [] sampleROVErrs;
            delete [] ROVStoreVals;
            delete [] ROVStoreGrads;
            delete [] ROVStorePts;
            delete [] ROVStoreEigens;
            delete [] Ystore;
            delete [] currentPt;
            return 0;
         }

         refineSeps = new int[2];
         refineSeps[0] = 0;
         refineSeps[1] = nSamples;
         nRefineCurr = 0;
         nRefines = 0;
         printAsterisks(PL_INFO, 0);
         iRov--;
      }
   }
 
   if (nSamples == 0)
   {
      psuadeIO_->getParameter("method_nsamples", pPtr);
      initNSamples = nSamples = pPtr.intData_; 
      nRefineLocal = 0;
      while (nSamples < 128000)
      {
         nRefineLocal++;
         nSamples *= 2;
      }
      nSamples = initNSamples; 
      psuadeIO_->getParameter("method_nrefinements", pPtr);
      nRefines = pPtr.intData_;
      refineSeps = new int[nRefineLocal+2];
      refineSeps[0] = 0;
      refineSeps[1] = nSamples;
      nRefineCurr = 1;
      while (nRefineCurr <= nRefineLocal)
      {
         sampler_->refine(2, 0, 0.0, nSamples, NULL);
         nSamples = sampler_->getNumSamples();
         refineSeps[nRefineCurr+1] = nSamples;
         printOutTS(PL_INFO, "PSUADE adaptiveGradBased: refinement %d (%d)\n",
                nRefineCurr, nSamples);
         nRefineCurr++;
      }
      nRefineCurr = 0;
      //cout << "PsuadeGradRSM: refinement done = " << nRefines
      //     << " (" << nSamples << ")\n";
      sampleInputs = new double[nSamples*nInputs];
      sampleROVVals = new double[nSamples];
      states = new int[nSamples];
      sampler_->getSamples(nSamples,nInputs,nOutputs,sampleInputs,
                           sampleROVVals,states);
      for (ii = 0; ii < nROVMax; ii++) states[ii] = 0;
      outNames = new char*[nOutputs];
      for (ii = 0; ii < nOutputs; ii++)
      {
         outNames[ii] = new char[100];
         sprintf(outNames[ii], "Y%d", ii+1);
      }
      inNames = new char*[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         inNames[ii] = new char[100];
         sprintf(inNames[ii], "X%d", ii+1);
      }
      psIO = new PsuadeData();
      psIO->updateMethodSection(PSUADE_SAMP_MC,nSamples,1,0,0);
      psIO->updateInputSection(nSamples,nInputs,NULL,iLowerB,iUpperB,
                               sampleInputs,inNames,NULL,NULL,NULL,NULL);
      psIO->updateOutputSection(nSamples,nOutputs,sampleROVVals,states,
                                outNames);
      psIO->writePsuadeFile("PsuadeGradBased.sample",0);
      delete psIO;
      for (ii = 0; ii < nInputs; ii++) delete [] inNames[ii];
      delete [] inNames;
      for (ii = 0; ii < nOutputs; ii++) delete [] outNames[ii];
      delete [] outNames;
      delete [] states;

      sampleROVInds = new int[nSamples];
      for (iS = 0; iS < nSamples; iS++) sampleROVInds[iS] = -1;
      sampleROVErrs = new double[nSamples];
      for (iS = 0; iS < nSamples; iS++) sampleROVErrs[iS] = PSUADE_UNDEFINED;

      ROVStorePts = new double[nROVMax*nInputs];
      ROVStoreVals = new double[nROVMax];
      ROVStoreGrads = new double[nROVMax*nInputs];
      ROVStoreEigens = new double[nROVMax];
      Ystore = new double[nROVMax*Ysize];
      iRov = 0;
      for (ii = 0; ii < nInputs; ii++) currentPt[ii] = sampleInputs[ii];
      sampleIndex = 0;
   }

   states = new int[nROVMax];
   for (ii = 0; ii < nROVMax; ii++) states[ii] = 1;
   outNames = new char*[Ysize];
   for (ii = 0; ii < Ysize; ii++)
   {
      outNames[ii] = new char[100];
      sprintf(outNames[ii], "Y%d", ii+1);
   }
   gradients = new double[nInputs];
   hessian = new double[nInputs*nInputs];
   jobz = 'N';
   uplo = 'U';
   lwork = 3 * nInputs;
   work = new double[lwork];
   eigs = new double[nInputs];

   funcIO = createFunctionInterface(psuadeIO_);
   ranges = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) ranges[ii] = (iUpperB[ii] - iLowerB[ii]);
   
   error = PSUADE_UNDEFINED;
   if (nRefines > 0)
   {
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO, "adaptiveGradBased: the first %d points will be evaluated.\n",
             initNSamples);
      printAsterisks(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
   }
   while (iRov < nROVMax && threshold >= anaThreshold)
   {
      iRov++;
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO, "ROV %4d = Sample %d\n", iRov, sampleIndex+1);

      evaluateFull(nInputs, currentPt, funcIO, ranges, maxParallelJobs, 
                   &YVal, gradients, hessian);
      if (outputLevel_ > 2) 
      {
         for (ii = 0; ii < nInputs; ii++)
         printOutTS(PL_INFO, "ROV %4d: input %2d value    = %e\n",iRov,ii+1,currentPt[ii]);
         printOutTS(PL_INFO, "ROV %4d: function value    = %e\n",iRov,YVal);
      }
      if (outputLevel_ > 3) 
      {
         for (ii = 0; ii < nInputs; ii++)
            printOutTS(PL_INFO, "ROV %4d: Gradient %3d      = %e\n",iRov,ii+1,
                   gradients[ii]);
         for (ii = 0; ii < nInputs; ii++)
            for (kk = 0; kk < nInputs; kk++)
            printOutTS(PL_INFO, "ROV %4d: Hessian (%3d,%3d) = %e\n",iRov,ii+1,kk+1,
                   hessian[ii*nInputs+kk]);
      }

      sampleROVInds[sampleIndex] = iRov;
      sampleROVErrs[sampleIndex] = 0.0;
      sampleROVVals[sampleIndex] = YVal;

      dsyev_(&jobz,&uplo,&nInputs,hessian,&nInputs,eigs,work,&lwork,&info);
      if (info != 0)
      {
         printOutTS(PL_INFO, "adaptiveGradBased INFO: dsyev returns a nonzero (%d).\n",
                info);
         exit(1);
      }
      eigMax = PABS(eigs[0]);
      for (ii = 1; ii < nInputs; ii++)
         if (PABS(eigs[ii]) > eigMax) eigMax = PABS(eigs[ii]);
      printOutTS(PL_INFO, "ROV %4d: max eigenvalue    = %e\n", iRov, eigMax);
      if (eigMax < minEigen)
      {
         eigMax = minEigen;
         printOutTS(PL_INFO, "ROV %4d: max value set to  = %e\n", iRov, eigMax);
      }

      for (ii = 0; ii < nInputs; ii++)
      {
         ROVStorePts[(iRov-1)*nInputs+ii] = currentPt[ii];
         ROVStoreGrads[(iRov-1)*nInputs+ii] = gradients[ii];
      }
      ROVStoreEigens[iRov-1] = eigMax;
      ROVStoreVals[iRov-1] = YVal;

      rovCnt = count = 0;
      for (iS = 0; iS < nSamples; iS++)
      {
         if (iS >= initNSamples || nRefines == 0) 
         {
            Ynew = YVal;
            dsum2 = 0.0;
            for (ii = 0; ii < nInputs; ii++)
            {
               diff = sampleInputs[iS*nInputs+ii] - currentPt[ii]; 
               Ynew += gradients[ii] * diff;
               dsum2 += diff * diff;
            }
            tolerance = 0.5 * dsum2 * eigMax;
            flag = 0;
            if (tolerance < sampleROVErrs[iS])
            {
               sampleROVInds[iS] = iRov - 1;
               sampleROVVals[iS] = Ynew;
               sampleROVErrs[iS] = tolerance;
               count++;
               flag = 1;
            }
            if (sampleROVErrs[iS] >= threshold)
            {
               sampleROVInds[iS] = - 1;
               if (flag == 1) count--;
            }
         }
         if (sampleROVInds[iS] >= 0) rovCnt++;
      }
      printOutTS(PL_INFO, "PSUADE adaptiveGradBased: region of validity limit not imposed.\n");
      printOutTS(PL_INFO, "PSUADE adaptiveGradBased: current threshold      = %e\n",
             threshold);
      printOutTS(PL_INFO, "PSUADE adaptiveGradBased: total number covered   = %d (%d)\n",
             rovCnt, nSamples);
      printOutTS(PL_INFO, "PSUADE adaptiveGradBased: number this ROV covers = %d\n",
             count+1);

      count = (iRov - 1) * Ysize;
      Ystore[count++] = YVal;
      for (ii = 0; ii < nInputs; ii++) Ystore[count++] = gradients[ii];
      for (ii = 0; ii < nInputs; ii++)
         for (kk = 0; kk <= ii; kk++)
            Ystore[count++] = hessian[ii*nInputs+kk];
      Ystore[count++] = eigMax;
      psuadeIO_->updateMethodSection(-1,iRov,1,0,0);
      psuadeIO_->updateInputSection(iRov,nInputs,NULL,NULL,NULL,
                                    ROVStorePts,NULL,NULL,NULL,NULL,NULL);
      psuadeIO_->updateOutputSection(iRov,Ysize,Ystore,states,outNames);
      psuadeIO_->writePsuadeFile("PsuadeGradBased.rov",0);

      for (iR = nRefineCurr; iR <= nRefineLocal; iR++)
      {
         for (kk = refineSeps[iR]; kk < refineSeps[iR+1]; kk++)
         {
            rovID = sampleROVInds[kk];
            if (rovID < 0) break;
         }
         sampleIndex = -1;
         if (kk < refineSeps[iR+1])
         {
            errMax = 0.0;
            for (kk = refineSeps[iR]; kk < refineSeps[iR+1]; kk++)
            {
               dtemp = sampleROVErrs[kk];
               if (dtemp > errMax)
               {
                  errMax = dtemp;
                  sampleIndex = kk;
               }
            }
         }
         if (sampleIndex >= 0) break;
         //printOutTS(PL_INFO, "PsuadeGradRSM: advance to refinement level %d (%d)\n",
         //       nRefineCurr, nRefines);
         nRefineCurr++;
      }

      if (sampleIndex < 0) 
      {
         error = 0.0;
         for (iS = 1; iS < nSamples; iS++)
            error += pow(sampleROVErrs[iS], 2.0);
         error = sqrt(error / (double) nSamples);
         printOutTS(PL_INFO, "PSUADE adaptiveGradBased: current mean error = %e\n", error);

         mean = 0.0;
         for (iS = 0; iS < nSamples; iS++) mean += sampleROVVals[iS];
         mean /= (double) nSamples;
         stdev = 0.0;
         for (iS = 0; iS < nSamples; iS++)
            stdev += pow(sampleROVVals[iS]-mean, 2.0);
         stdev /= (double) nSamples;
         stdev = sqrt(stdev);
         printOutTS(PL_INFO, "PSUADE adaptiveGradBased: Current sample mean   = %16.8e\n",
                mean);
         printOutTS(PL_INFO, "PSUADE adaptiveGradBased: Current std deviation = %16.8e\n",
                stdev);
      
         while (threshold > anaThreshold)
         {
            threshold *= threshDec;
            printOutTS(PL_INFO, "PSUADE adaptiveGradBased: threshold down to  = %e\n",
                   threshold);
            if (threshold <= anaThreshold) threshold = anaThreshold;
            for (iS = 0; iS < nSamples; iS++)
            {
               kk = sampleROVInds[iS];
               if (kk >= 0)
               {
                  Ynew = ROVStoreVals[kk];
                  eigMax = ROVStoreEigens[kk]; 
                  dsum2 = 0.0;
                  for (ii = 0; ii < nInputs; ii++)
                  {
                     diff = sampleInputs[iS*nInputs+ii] - 
                            ROVStorePts[kk*nInputs+ii]; 
                     Ynew += ROVStoreGrads[kk*nInputs+ii] * diff;
                     dsum2 += diff * diff;
                  }
                  tolerance = 0.5 * dsum2 * eigMax;
                  if (tolerance > threshold)
                  {
                     sampleROVErrs[iS] = PSUADE_UNDEFINED;
                     for (iR = 0; iR < iRov; iR++)
                     {
                        eig2 = ROVStoreEigens[iR]; 
                        dsum2 = 0.0;
                        Ynew = ROVStoreVals[iR];
                        for (ii = 0; ii < nInputs; ii++)
                        {
                           diff = sampleInputs[iS*nInputs+ii] - 
                                  ROVStorePts[iR*nInputs+ii]; 
                           Ynew += ROVStoreGrads[iR*nInputs+ii] * diff;
                           dsum2 += diff * diff;
                        }
                        tol2 = 0.5 * dsum2 * eig2;
                        if (tol2 < sampleROVErrs[iS])
                        {
                           sampleROVErrs[iS] = tol2;
                           sampleROVInds[iS] = iR;
                           sampleROVVals[iS] = Ynew;
                        }
                     }
                     if (sampleROVErrs[iS] > threshold) sampleROVInds[iS] = -1;
                  }
               }
            }
            nRefineCurr = 0;
            for (iR = nRefineCurr; iR <= nRefineLocal; iR++)
            {
               for (kk = refineSeps[iR]; kk < refineSeps[iR+1]; kk++)
               {
                  rovID = sampleROVInds[kk];
                  if (rovID < 0) break;
               }
               sampleIndex = -1;
               if (kk < refineSeps[iR+1])
               {
                  errMax = 0.0;
                  for (kk = refineSeps[iR]; kk < refineSeps[iR+1]; kk++)
                  {
                     dtemp = sampleROVErrs[kk];
                     if (dtemp > errMax)
                     {
                        errMax = dtemp;
                        sampleIndex = kk;
                     }
                  }
               }
               if (sampleIndex >= 0) break;
               nRefineCurr++;
            }
            if (sampleIndex >= 0) break;
         }
      }
      if (sampleIndex < 0)
      {
         printAsterisks(PL_INFO, 0);
         printOutTS(PL_INFO, "Cannot find next ROV to evaluate (%d).\n",sampleIndex);
         break;
      }
      if (iRov < initNSamples && nRefines > 0) sampleIndex = iRov;
      printOutTS(PL_INFO, "PSUADE adaptiveGradBased: next candidate ROV = %d\n",
             sampleIndex+1);
      for (ii = 0; ii < nInputs; ii++)
         currentPt[ii] = sampleInputs[sampleIndex*nInputs+ii];
      printAsterisks(PL_INFO, 0);
   }

   for (iS = 0; iS < nSamples; iS++)
   {
      sampleROVErrs[iS] = PSUADE_UNDEFINED;
      sampleROVInds[iS] = -1;
      for (iR = 0; iR < iRov; iR++)
      {
         eig2 = ROVStoreEigens[iR]; 
         dsum2 = 0.0;
         Ynew = ROVStoreVals[iR];
         for (ii = 0; ii < nInputs; ii++)
         {
            diff = sampleInputs[iS*nInputs+ii] - 
                   ROVStorePts[iR*nInputs+ii]; 
            Ynew += ROVStoreGrads[iR*nInputs+ii] * diff;
            dsum2 += diff * diff;
         }
         tol2 = 0.5 * dsum2 * eig2;
         if (tol2 < sampleROVErrs[iS])
         {
            sampleROVErrs[iS] = tol2;
            sampleROVInds[iS] = iR;
            sampleROVVals[iS] = Ynew;
         }
         if (sampleROVErrs[iS] >= anaThreshold) sampleROVInds[iS] = -1;
      }
      if (outputLevel_ > 4)
         printOutTS(PL_INFO, "Sample %d's (%d)index = %d, %e (%e)\n", iS+1,
                sampleROVInds[iS],nSamples,sampleROVErrs[iS],anaThreshold);
      if (sampleROVInds[iS] == -1)
      {
         printOutTS(PL_ERROR, "PSUADE adaptiveGradBased: fatal ERROR. Consult developers.\n");
         exit(1);
      }
   }

   mean = 0.0;
   for (iS = 0; iS < nSamples; iS++) mean += sampleROVVals[iS];
   mean /= (double) nSamples;
   stdev = 0.0;
   for (iS = 0; iS < nSamples; iS++)
      stdev += pow(sampleROVVals[iS]-mean, 2.0);
   stdev /= (double) nSamples;
   stdev = sqrt(stdev);
   printOutTS(PL_INFO, "PSUADE adaptiveGradBased: Final sample mean     = %16.8e\n", mean);
   printOutTS(PL_INFO, "PSUADE adaptiveGradBased: Final std deviation   = %16.8e\n", stdev);
   printOutTS(PL_INFO, "PSUADE adaptiveGradBased: Total number of ROVs  = %d\n", iRov);
   printOutTS(PL_INFO, "PSUADE adaptiveGradBased: Total number of runs  = %d\n",
          (Ysize-1)*iRov);
   if (iRov > 0)
   {
      printOutTS(PL_INFO, "PSUADE adaptiveGradBased: ROV info is in PsuadeGradBased.rov.\n");
   }

   if (sampler_ != NULL)
   {
      delete sampler_;
      sampler_ = NULL;
   }
   delete funcIO;
   delete [] work;
   delete [] eigs;
   delete [] sampleInputs;
   delete [] sampleROVInds;
   delete [] sampleROVVals;
   delete [] sampleROVErrs;
   delete [] ROVStoreVals;
   delete [] ROVStoreGrads;
   delete [] ROVStorePts;
   delete [] ROVStoreEigens;
   delete [] currentPt;
   delete [] gradients;
   delete [] hessian;
   delete [] refineSeps;
   for (ii = 0; ii < Ysize; ii++) delete [] outNames[ii];
   delete [] outNames;
   delete [] states;
   delete [] ranges;
   delete [] Ystore;
   return 0;
}

// ************************************************************************
// evaluate the function, its gradients, and Hessian
// ------------------------------------------------------------------------
int PsuadeBase::evaluateFull(int nInputs, double *currentPt, 
                             FunctionInterface *funcIO,
                             double *ranges, int maxJobs, double *YVal, 
                             double *gradients, double *hessian)
{
   int     launchInterval=1, ii, nJobs, iS, iS2, nOutputs=1, *sampleStates;
   int     jobCnt;
   double  *sampleInputs, *sampleOutputs, hstep, hstep2;
   JobCntl jcntl;

   if (maxJobs > 1) funcIO->setAsynchronousMode();
   funcIO->setLaunchInterval(launchInterval);

   nJobs = nInputs+1+(nInputs+1)*nInputs/2;
   sampleStates = new int[nJobs];
   sampleInputs = new double[nJobs*nInputs];
   sampleOutputs = new double[nJobs];
   for (iS = 0; iS < nJobs; iS++) sampleStates[iS] = 0;
   for (ii = 0; ii < nInputs; ii++) sampleInputs[ii] = currentPt[ii];
   jobCnt = 1;
   if (gradients != NULL)
   {
      for (iS = 0; iS < nInputs; iS++)
      {
         hstep = 1.0e-5 * ranges[iS];
         for (ii = 0; ii < nInputs; ii++)
            sampleInputs[jobCnt*nInputs+ii] = currentPt[ii];
         sampleInputs[jobCnt*nInputs+iS] = currentPt[iS] + hstep;
         jobCnt++;
      }
   }
   if (gradients != NULL && hessian != NULL)
   {
      for (iS = 0; iS < nInputs; iS++)
      {
         hstep = 1.0e-5 * ranges[iS];
         for (iS2 = 0; iS2 <= iS; iS2++)
         {
            hstep2 = 1.0e-5 * ranges[iS2];
            for (ii = 0; ii < nInputs; ii++)
               sampleInputs[jobCnt*nInputs+ii] = 
                                   sampleInputs[(1+iS)*nInputs+ii];
            sampleInputs[jobCnt*nInputs+iS2] += hstep2; 
            jobCnt++;
         }
      }
   }

   jcntl.setMaxParallelJobs(maxJobs);
   //jcntl.setOutputLevel(outputLevel_);
   jcntl.loadFuncIO(funcIO);
   jcntl.loadSample(jobCnt,nInputs,nOutputs,sampleInputs,sampleOutputs,
                    sampleStates);
   jcntl.execute();
   jcntl.getSampleOutputs(jobCnt, nOutputs, sampleOutputs); 

   (*YVal) = sampleOutputs[0]; 
   if (gradients != NULL)
   {
      for (iS = 0; iS < nInputs; iS++)
      {
         hstep = 1.0e-5 * ranges[iS];
         gradients[iS] = (sampleOutputs[1+iS] - sampleOutputs[0]) / hstep;
      }
      jobCnt = nInputs + 1;
   }
   if (gradients != NULL && hessian != NULL)
   {
      for (iS = 0; iS < nInputs; iS++)
      {
         hstep = 1.0e-5 * ranges[iS];
         for (iS2 = 0; iS2 <= iS; iS2++)
         {
            hstep2 = 1.0e-5 * ranges[iS2];
            hessian[iS*nInputs+iS2] = 
                ((sampleOutputs[jobCnt]-sampleOutputs[iS2+1]) - 
                 (sampleOutputs[iS+1]-sampleOutputs[0]))/(hstep*hstep2);
            hessian[iS2*nInputs+iS] = hessian[iS*nInputs+iS2]; 
            jobCnt++;
         }
      }
   }
   delete [] sampleInputs;
   delete [] sampleOutputs;
   delete [] sampleStates;
   return 0;
}

