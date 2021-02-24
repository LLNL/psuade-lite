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
// Ref: S. Lophaven H. B. Nielsen, and J. Sondergaard, "DACE: A Matlab
//      Kriging toolbox," Technical Report IMM-TR-2002-12, Informatics
//      and Mathematical Modelling, Technical University of Denmark.
// ************************************************************************
// Functions for the class PKriging
//   This module is intended primarily for parallel generation of 
//      Kriging coefficients. Evaluations should be done by Kriging module.
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "PKriging.h"
#include "sysdef.h"
#include "dtype.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PsuadeConfig.h"
#include "Sampling.h"
#include "PrintingTS.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// ************************************************************************
// external functions
// ------------------------------------------------------------------------
extern "C" {
   void dtrsv_(char *, char *, char *, int *, double *, int *, double *, int *);
   void dpotrf_(char *, int *, double *, int *, int *);
   void dpotrs_(char *, int *, int *, double *, int *, double *,
                int *, int *);
   void kbobyqa_(int *,int *, double *, double *, double *, double *,
                double *, int *, int *, double*);
}

// ************************************************************************
// ************************************************************************
// internal global variables
// ------------------------------------------------------------------------
int    PKRI_outputLevel=0;
int    PKRI_iter=-1;
int    PKRI_nInputs=-1;
int    PKRI_nSamples=-1;
int    PKRI_pOrder=-1;
double *PKRI_XDists=NULL;
double *PKRI_SMatrix=NULL;
double *PKRI_FMatrix=NULL;
double *PKRI_FMatTmp=NULL;
double *PKRI_MMatrix=NULL;
double *PKRI_X=NULL;
double *PKRI_Y=NULL;
double PKRI_currY=0.0;
double *PKRI_Ytmp=NULL;
double PKRI_OptY=0.0;
double *PKRI_OptThetas=NULL;
double *PKRI_dataStdDevs=NULL;
double PKRI_YStd=1.0;
int    PKRI_terminate=0;
int    PKRI_noProgressCnt=0;

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif
   void *pkribobyqaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    nBasis, count, Cleng, ii, jj, kk, status;
      double dist, ddata, ddata2;
      char   uplo='L';
      FILE   *fp=NULL;

      if (PKRI_noProgressCnt >= 200)
      {
         if (PKRI_outputLevel >= 1)
         {
            //PKRI_outputLevel = 0;
            printf("\t*** no progress for more than 200 iterations.\n");
         }
         (*YValue) = PSUADE_UNDEFINED;
         return NULL;
      }
      fp = fopen("psuade_stop", "r");
      if (fp != NULL)
      {
         fclose(fp);
         printf("Kriging: psuade_stop file found - terminating ....\n");
         (*YValue) = PKRI_OptY;
         return NULL;
      }
      PKRI_iter++;
      nBasis = 1;
      if (PKRI_pOrder == 1) nBasis = PKRI_nInputs + 1;
      Cleng = PKRI_nSamples + nBasis;

      count = 0;
      for (jj = 0; jj < PKRI_nSamples; jj++)
      {
         PKRI_SMatrix[jj*PKRI_nSamples+jj] = 1.0 + 1e-15 * PKRI_nSamples;
         if (PKRI_nInputs == 2) PKRI_SMatrix[jj*PKRI_nSamples+jj] += 0.01;

         if (PKRI_dataStdDevs != NULL) 
         {
            ddata = pow(PKRI_dataStdDevs[jj]/PKRI_YStd,2.0);
            PKRI_SMatrix[jj*PKRI_nSamples+jj] += ddata;
         }
         for (kk = jj+1; kk < PKRI_nSamples; kk++)
         {
            dist = 0.0;
            for (ii = 0; ii < PKRI_nInputs; ii++)
               dist += pow(PKRI_XDists[count*PKRI_nInputs+ii]/XValues[ii], 2.0);
            dist = exp(-dist);
            if (dist < 1.0e-16) dist = 0;
            PKRI_SMatrix[jj*PKRI_nSamples+kk] = dist;
            PKRI_SMatrix[kk*PKRI_nSamples+jj] = dist;
            count++;
         }
      }

      for (jj = 0; jj < PKRI_nSamples; jj++)
      {
         PKRI_FMatrix[jj] = 1.0;
         PKRI_FMatTmp[jj] = 1.0;
         if (PKRI_pOrder == 1)
         {
            for (kk = 1; kk < PKRI_nInputs+1; kk++)
            {
               ddata = PKRI_X[jj*PKRI_nInputs+kk-1];
               PKRI_FMatTmp[kk*PKRI_nSamples+jj] = ddata;
               PKRI_FMatrix[kk*PKRI_nSamples+jj] = ddata;
            }
         }
      }
      dpotrf_(&uplo, &PKRI_nSamples, PKRI_SMatrix, &PKRI_nSamples, &status);
      kk = nBasis;
      dpotrs_(&uplo, &PKRI_nSamples, &kk, PKRI_SMatrix, &PKRI_nSamples,
              PKRI_FMatTmp, &PKRI_nSamples, &status);
      for (jj = 0; jj < nBasis; jj++)
      {
         for (kk = 0; kk < nBasis; kk++)
         {
            ddata = 0.0;
            for (ii = 0; ii < PKRI_nSamples; ii++)
               ddata += PKRI_FMatrix[ii+jj*PKRI_nSamples]*
                        PKRI_FMatTmp[ii+kk*PKRI_nSamples];
            PKRI_MMatrix[jj+kk*nBasis] = ddata;
         }
      }
      if (nBasis > 1) dpotrf_(&uplo, &nBasis, PKRI_MMatrix, &nBasis, &status);
      for (jj = 0; jj < PKRI_nSamples; jj++) PKRI_Ytmp[jj] = PKRI_Y[jj];
      kk = 1;
      dpotrs_(&uplo,&PKRI_nSamples,&kk,PKRI_SMatrix,&PKRI_nSamples,PKRI_Ytmp,
              &PKRI_nSamples, &status);
      for (jj = 0; jj < nBasis; jj++)
      {
         ddata = 0.0;
         for (ii = 0; ii < PKRI_nSamples; ii++)
            ddata += PKRI_FMatrix[ii+jj*PKRI_nSamples]*PKRI_Ytmp[ii];
         PKRI_Ytmp[PKRI_nSamples+jj] = ddata;
      }
      if (nBasis == 1)
      {
         if (PKRI_MMatrix[0] == 0)
         {
            printf("PKriging ERROR: divide by 0.\n");
            exit(1);
         }
         PKRI_Ytmp[PKRI_nSamples] /= PKRI_MMatrix[0];
      }
      else
      {
         kk = 1;
         dpotrs_(&uplo,&kk,&kk,PKRI_MMatrix,&nBasis,&PKRI_Ytmp[PKRI_nSamples],
                 &kk, &status);
      }
      for (jj = 0; jj < PKRI_nSamples; jj++)
      {
         ddata = 0.0;
         for (ii = 0; ii < nBasis; ii++)
            ddata += PKRI_FMatrix[jj+ii*PKRI_nSamples]*PKRI_Ytmp[ii+PKRI_nSamples];
         PKRI_Ytmp[jj] = PKRI_Y[jj] - ddata;
      }
      kk = 1;
      dpotrs_(&uplo,&PKRI_nSamples,&kk,PKRI_SMatrix,&PKRI_nSamples,PKRI_Ytmp,
              &PKRI_nSamples, &status);

      ddata = 0.0;
      for (jj = 0; jj < PKRI_nSamples; jj++)
         ddata += PKRI_Ytmp[jj] * PKRI_Y[jj];
      ddata /= (double) PKRI_nSamples;

      for (jj = 0; jj < PKRI_nSamples; jj++) 
      {
         ddata2 = PKRI_SMatrix[PKRI_nSamples*jj+jj];
         if (ddata2 < 0.0) ddata2 = - ddata2;
         ddata *= pow(ddata2, 2.0/(double) PKRI_nSamples);
      }

      (*YValue) = PKRI_currY = ddata;

      if (ddata < PKRI_OptY || PKRI_iter == 1)
      {
         PKRI_noProgressCnt = 0;
         PKRI_OptY = ddata;
         for (ii = 0; ii < PKRI_nInputs; ii++)
            PKRI_OptThetas[ii] = XValues[ii];
         if (PKRI_outputLevel > 1)
         {
            printf("\t Kriging : iteration %d\n",PKRI_iter);
            for (ii = 0; ii < PKRI_nInputs; ii++)
               printf("\t    Current best theta %d = %e\n",ii+1,XValues[ii]);
            printf("\t    Current best objective value = %e\n",ddata);
         }
      }
      else PKRI_noProgressCnt++;

      if (psRSExpertMode_ == 1 && PKRI_outputLevel > 3)
      {
         printf("\t PKriging : iteration %d\n", PKRI_iter);
         for (ii = 0; ii < PKRI_nInputs; ii++)
            printf("\t    Current theta %d = %e\n", ii+1, XValues[ii]);
         printf("\t    Current Kriging objective value = %e\n", ddata);
         printf("\t* For early termination, just create a file\n");
         printf("\t* called 'psuade_stop' in your local directory.\n");
         fp = fopen("psuade_stop", "r");
         if (fp != NULL)
         {
            printf("\t*** psuade_stop file found.\n");
            (*YValue) = 0.0;
            fclose(fp);
         }
      }
      return NULL;
   }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// Constructor for object class Kriging
// ------------------------------------------------------------------------
PKriging::PKriging(int nInputs,int nSamples, CommManager *comm) : 
                   FuncApprox(nInputs,nSamples)
{
   int    ii, jj, pstatus, iOne=1;
   char   pString[500], winput[500], winput2[500], fname[500], *strPtr;
   FILE   *fp;

   // initialize variables and parameters
   faID_ = PSUADE_RS_KR;
   commMgr_ = comm;
   if (comm == NULL)
   {
      printOutTS(PL_ERROR,"PKriging ERROR: no communicator.\n");
      exit(1);
   }
   mypid_  = psCommMgr_->getPID();
   nprocs_ = psCommMgr_->getNumProcs();

   XNormalized_ = NULL;
   YNormalized_ = NULL;
   pOrder_  = 0;
   dataStdDevs_ = NULL;
   optTolerance_ = 1.0e-4;
   Thetas_ = new double[nInputs_+1];
   checkAllocate(Thetas_, "Thetas_ in PKriging::constructor");
   for (ii = 0; ii <= nInputs_; ii++) Thetas_[ii] = 0.01;

   // display banner and additonal information
   if (mypid_ == 0)
   {
      printAsterisks(PL_INFO, 0);
      printf("*                Kriging Analysis\n");
      printf("* Set printlevel to 1-4 to see Kriging details.\n");
      printf("* Create 'psuade_stop' file to gracefully terminate.\n");
      printf("* Create 'ps_print' file to set print level on the fly.\n");
      printEquals(PL_INFO, 0);
   }

   // read configure file, if any 
   if (mypid_ == 0 && psConfig_ != NULL)
   {
      strPtr = psConfig_->getParameter("KRI_DATA_STDEV_FILE");
      if (strPtr != NULL)
      {
         sscanf(strPtr, "%s %s %s", winput, winput2, fname);
         fp = fopen(fname, "r");
         if (fp == NULL)
         {
            printf("Kriging INFO: data variance file %s not found.\n",fname);
         }
         else
         {
            fscanf(fp, "%d", &ii); 
            if (ii != nSamples_)
            {
               printf("Kriging ERROR: std. dev. file should have %d entries.\n",
                      nSamples_);
               fclose(fp);
            }
            else
            {
               dataStdDevs_ = new double[nSamples_];
               checkAllocate(dataStdDevs_,
                             "dataStdDevs_ in PKriging::constructor");
               for (ii = 0; ii < nSamples_; ii++)
               {
                  fscanf(fp, "%d %lg", &jj, &dataStdDevs_[ii]); 
                  if (jj != ii+1)
                  {
                     printf("Kriging ERROR: line %d in the std dev file.\n",
                            jj+1);
                     delete [] dataStdDevs_;
                     dataStdDevs_ = NULL;
                     break;
                  }
               } 
               fclose(fp);
               if (dataStdDevs_ != NULL)
               {
                  printf("Kriging INFO: std. dev. file has been read.\n");
                  PKRI_dataStdDevs = dataStdDevs_; 
               }
            }
         }
      }
   }
   pstatus = nSamples_;
   if (mypid_ == 0 && dataStdDevs_ == NULL) pstatus = 0; 
   psCommMgr_->bcast((void *) &pstatus, iOne, INT, 0);
   if (mypid_ != 0 && pstatus != 0) dataStdDevs_ = new double[nSamples_];
   psCommMgr_->bcast((void *) dataStdDevs_, pstatus, DOUBLE, 0);
   PKRI_dataStdDevs = dataStdDevs_; 
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PKriging::~PKriging()
{
   if (Thetas_ != NULL) delete [] Thetas_;
   if (XNormalized_ != NULL) delete [] XNormalized_;
   if (YNormalized_ != NULL) delete [] YNormalized_;
   if (dataStdDevs_ != NULL) delete [] dataStdDevs_;
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int PKriging::initialize(double *X, double *Y)
{
   train(X,Y);
   return 0;
}
  
// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int PKriging::genNDGridData(double *, double *, int *, double **, double **)
{
   printf("PKriging genNDGridData INFO: use Kriging\n");
   return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int PKriging::gen1DGridData(double *, double *, int, double *, int *, 
                            double **, double **)
{
   printf("PKriging gen1DGridData INFO: use Kriging\n");
   return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int PKriging::gen2DGridData(double *, double *, int, int, double *, int *, 
                            double **, double **)
{
   printf("PKriging gen2DGridData INFO: use Kriging\n");
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int PKriging::gen3DGridData(double *, double *, int, int, int, double *, 
                            int *, double **, double **)
{
   printf("PKriging gen3DGridData INFO: use Kriging\n");
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int PKriging::gen4DGridData(double *, double *, int, int, int, int, double *,
                            int *, double **, double **)
{
   printf("PKriging gen4DGridData INFO: use Kriging\n");
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double PKriging::evaluatePoint(double *)
{
   printf("PKriging evaluatePoint INFO: use Kriging\n");
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double PKriging::evaluatePoint(int, double *, double *)
{
   printf("PKriging evaluatePoint INFO: use Kriging\n");
   return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double PKriging::evaluatePointFuzzy(double *, double &)
{
   printf("PKriging evaluatePointFuzzy INFO: use Kriging\n");
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double PKriging::evaluatePointFuzzy(int, double *, double *, double *)
{
   printf("PKriging evaluatePointFuzzy INFO: use Kriging\n");
   return 0.0;
}

// ************************************************************************
// training 
// ------------------------------------------------------------------------
double PKriging::train(double *X, double *Y)
{
   int    ii, jj, kk, nBasis, count, status, nDists, Cleng;
   double *XDists, dist, ddata; 
   char   pString[500];

   // clean up 
   if (XNormalized_ != NULL) delete [] XNormalized_;
   if (YNormalized_ != NULL) delete [] YNormalized_;
   XNormalized_ = NULL;
   YNormalized_ = NULL;

   // normalize the input and outputs 
   XNormalized_ = new double[nSamples_*nInputs_];
   checkAllocate(XNormalized_,"XNormalized_ in PKriging::train");
   initInputScaling(X, XNormalized_, 1);
   YNormalized_ = new double[nSamples_];
   checkAllocate(YNormalized_,"YNormalized_ in PKriging::train");
   initOutputScaling(Y, YNormalized_);
   PKRI_YStd = YStd_;

   // compute distances between all pairs of inputs (nDists, XDists)
   computeDistances(&XDists, &nDists);

   // slow mode: optimize
   int    maxfun, pLevel, nPts, iOne=1, iZero=0, nSamOpt, *samStates;
   double *TValues, *TUppers, *TLowers, rhobeg=1.0, rhoend=1.0e-4;
   double *work, *samInputs, *samOutputs, *optThetas, *MMatrix=NULL;
   double optY=PSUADE_UNDEFINED, *SMatrix=NULL, *FMatrix=NULL;
   double *FMatTmp;
   FILE   *fp=NULL;
   Sampling *sampler;

   status = 0;
   fp = fopen("psuade_stop", "r");
   if (mypid_ == 0 && fp != NULL)
   {
      printf("PKriging ERROR: remove the 'psuade_stop' file\n");
      printf("                first and re-do.\n");
      fclose(fp);
      status = 1;
   }
   psCommMgr_->bcast((void *) &status, iOne, INT, 0);
   if (status == 1) exit(1);

   if (nSamples_ <= nInputs_ + 1) pOrder_ = 0;
   if (pOrder_ == 1) nBasis = nInputs_ + 1;
   else              nBasis = 1;
   Cleng = nSamples_ + nBasis;

   TUppers = new double[nInputs_+1];
   TLowers = new double[nInputs_+1];
   checkAllocate(TLowers,"TLowers in PKriging::train");
   for (ii = 0; ii < nInputs_; ii++)
   {
      if (XMeans_[ii] == 0 && XStds_[ii] == 1.0)
      {
         TUppers[ii] = 20.0 * (upperBounds_[ii]-lowerBounds_[ii]);;
         TLowers[ii] = 0.1 * (upperBounds_[ii]-lowerBounds_[ii]);;
      }
      else
      {
         TUppers[ii] = 20.0;
         TLowers[ii] = 0.1;
      }
   }
   status = 0;
   if (mypid_ == 0 && psMasterMode_ == 1)
   {
      printf("PKriging: current optimization lower bound for input %d = %e",
              ii+1,TLowers[ii]);
      sprintf(pString,
              "PKriging: Enter optimization lower bound for input %d : ",ii+1);
      TLowers[ii] = getDouble(pString);
      if (TLowers[ii] <= 0.0)
      {
         printf("PKriging ERROR: lower bound <= 0\n");
         exit(1);
      }
      printf("PKriging: current optimization upper bound for input %d = %e",
             ii+1,TUppers[ii]);
      sprintf(pString,
             "PKriging: Enter optimization upper bound for input %d : ",ii+1);
      TUppers[ii] = getDouble(pString);
      if (TLowers[ii] > TUppers[ii])
      {
         printf("PKriging ERROR: lower bound >= upper bound\n");
         status = 1;
      }
   }
   psCommMgr_->bcast((void *) &status, iOne, INT, 0);
   if (status == 1) exit(1);
   kk = nInputs_ + 1;
   psCommMgr_->bcast((void *) TUppers, kk, DOUBLE, 0);
   psCommMgr_->bcast((void *) TLowers, kk, DOUBLE, 0);

   rhobeg = TUppers[0] - TLowers[0];
   for (ii = 1; ii < nInputs_; ii++)
   {
      ddata = TUppers[ii] - TLowers[ii];
      if (ddata < rhobeg) rhobeg = ddata;
   }
   rhobeg *= 0.5;
   rhoend = rhobeg * optTolerance_;
   TValues = new double[nInputs_+1];
   nPts = (nInputs_ + 1) * (nInputs_ + 2) / 2;
   work = new double[(nPts+5)*(nPts+nInputs_)+3*nInputs_*(nInputs_+5)/2+1];
   checkAllocate(work,"work in PKriging::train");
   if (mypid_ == 0)
   {
      printEquals(PL_INFO, 0);
      printf("* PKriging optimization tolerance = %e\n", rhoend);
   }

   if (mypid_ == 0)
      printf("PKriging training begins.... (order = %d)\n",pOrder_);
   nSamOpt = 10;
   if (nprocs_ > nSamOpt) nSamOpt = nprocs_;
   if (psRSExpertMode_ == 1 && mypid_ == 0)
   {
      printf("To thoroughly explore the parameter space, multi-start\n");
      printf("optimization is to be employed. Please enter the number\n");
      printf("of multi-starts (more the better, but also more expensive.\n");
      printf("Default is max(num. of processors, 10). \n");
      sprintf(pString,"Enter number of multi-starts (up to 100): ");
      nSamOpt = getInt(nSamOpt, 100, pString);
   }
   if (nInputs_ > 50)
        sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
   else sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
 
   sampler->setPrintLevel(0);
   sampler->setInputBounds(nInputs_, TLowers, TUppers);
   sampler->setOutputParams(iOne);
   sampler->setSamplingParams(nSamOpt, iOne, iZero);
   sampler->initialize(0);
   nSamOpt = sampler->getNumSamples();
   samInputs  = new double[nSamOpt * nInputs_];
   samOutputs = new double[nSamOpt];
   samStates  = new int[nSamOpt];
   checkAllocate(samStates,"samStates in PKriging::train");
   sampler->getSamples(nSamOpt, nInputs_, iOne, samInputs,
                       samOutputs, samStates);
   delete [] samOutputs;
   delete [] samStates;
   delete sampler;

   PKRI_XDists = XDists;
   PKRI_nSamples = nSamples_;
   PKRI_nInputs = nInputs_;
   PKRI_pOrder = pOrder_;
   PKRI_X = XNormalized_;
   PKRI_Y = YNormalized_;
   PKRI_Ytmp = new double[Cleng];
   PKRI_OptThetas = Thetas_;
   maxfun = 5000;
   PKRI_OptY = PSUADE_UNDEFINED;
   PKRI_outputLevel = outputLevel_;
   optThetas = new double[nInputs_];
   SMatrix = new double[nSamples_*nSamples_];
   FMatrix = new double[nSamples_*nBasis];
   FMatTmp = new double[nSamples_*nBasis];
   MMatrix = new double[nBasis*nBasis];
   checkAllocate(MMatrix,"MMatrix in PKriging::train");
   PKRI_SMatrix = SMatrix;
   PKRI_FMatrix = FMatrix;
   PKRI_FMatTmp = FMatTmp;
   PKRI_MMatrix = MMatrix;
   int    proc, csize;
   double dmean, dstd;
   double *commBuffer = new double[nInputs_+1];
   checkAllocate(commBuffer,"commBuffer in PKriging::train");

   for (kk = 0; kk < nSamOpt; kk++)
   {
      if (kk % nprocs_ == mypid_)
      {
        if (outputLevel_ >= 1) 
           printf("%4d: PKriging multi-start optimization: start = %d (%d)\n",
                  mypid_, kk+1, nSamOpt);
         PKRI_iter = 0;
         pLevel = 4444;
         for (ii = 0; ii < nInputs_; ii++) 
            TValues[ii] = samInputs[kk*nInputs_+ii];
#ifdef HAVE_BOBYQA
         PKRI_noProgressCnt = 0;
         kbobyqa_(&nInputs_,&nPts,TValues,TLowers,TUppers,&rhobeg,&rhoend,
                  &pLevel, &maxfun, work);
#else
         printf("PKriging ERROR: Bobyqa optimizer not installed.\n");
         exit(1);
#endif
         if (outputLevel_ >= 3) 
         {
            printf("%4d: multi-start optimization: iteration = %d (%d)\n",
                   mypid_, kk+1, nSamOpt);
            for (ii = 0; ii < nInputs_; ii++) 
               printf("%4d: Input %4d final length scale = %e\n",
                      mypid_, ii+1,TValues[ii]);
            printf("%4d: final objective value = %e\n", mypid_, PKRI_currY);
         }
         if (PKRI_OptY < optY)
         {
            optY = PKRI_OptY;
            for (ii = 0; ii < nInputs_; ii++)
               optThetas[ii] = PKRI_OptThetas[ii];
         }
         fp = fopen("psuade_stop", "r");
         if (fp != NULL)
         {
            fclose(fp);
            fp = NULL;
            printf("PKriging: psuade_stop file found - terminating ....\n");
            break;
         }
      }
   }
   for (kk = 0; kk < nSamOpt; kk++)
   {
      if (mypid_ == 0 && (kk % nprocs_ != mypid_))
      {
         proc = kk % nprocs_;
         csize = nInputs_ + 1; 
         psCommMgr_->recv((void *) commBuffer,csize,DOUBLE,kk,proc);
         if (commBuffer[nInputs_] < optY)
         {
            for (ii = 0; ii < nInputs_; ii++) optThetas[ii] = commBuffer[ii];
            optY = commBuffer[nInputs_];
         }
      }
      else if (mypid_ != 0 && kk % nprocs_ == mypid_)
      {
         for (ii = 0; ii < nInputs_; ii++) commBuffer[ii] = optThetas[ii];
         commBuffer[nInputs_] = optY;
         csize = nInputs_ + 1; 
         psCommMgr_->send((void *) optThetas,csize,DOUBLE,kk,0);
      }
   }
   if (mypid_ == 0)
   {
      for (ii = 0; ii < nInputs_; ii++) commBuffer[ii] = optThetas[ii];
      commBuffer[nInputs_] = optY;
   }
   csize = nInputs_ + 1;
   psCommMgr_->bcast((void *) commBuffer, csize, DOUBLE, 0);
   for (ii = 0; ii < nInputs_; ii++) optThetas[ii] = commBuffer[ii]; 
   optY = commBuffer[nInputs_];

   for (ii = 0; ii < nInputs_; ii++) Thetas_[ii] = optThetas[ii];
   if (mypid_ == 0) 
   {
      printAsterisks(PL_INFO,0);
      for (ii = 0; ii < nInputs_; ii++) 
         printf("PKriging: Input %4d optimal length scale = %e\n",
                ii+1,Thetas_[ii]);
      printf("PKriging: optimal objective value = %e\n",PKRI_OptY);
      printAsterisks(PL_INFO,0);
   }
   if (mypid_ == 0) 
   {
      fp = fopen("psuade_kriging_optdata","w");
      for (ii = 0; ii < nInputs_; ii++) 
         fprintf(fp, "%16.8e\n", Thetas_[ii]);
      printf("PKriging: length scales are saved in 'psuade_kriging_optdata'.\n");
      fclose(fp);
   }
   delete [] work;
   delete [] TValues;
   delete [] TUppers;
   delete [] TLowers;
   delete [] PKRI_Ytmp;
   delete [] samInputs;
   delete [] optThetas;
   PKRI_OptThetas = NULL;
   PKRI_XDists = NULL;
   PKRI_X = NULL;
   PKRI_Y = NULL;
   PKRI_Ytmp = NULL;
   if (SMatrix != NULL) delete [] SMatrix;
   if (FMatrix != NULL) delete [] FMatrix;
   if (FMatTmp != NULL) delete [] FMatTmp;
   if (MMatrix != NULL) delete [] MMatrix;
   PKRI_SMatrix = NULL;
   PKRI_FMatrix = NULL;
   PKRI_FMatTmp = NULL;
   PKRI_MMatrix = NULL;
   delete [] XDists;
   return 0.0;
}

// ************************************************************************
// compute distances between all pairs of inputs
// ------------------------------------------------------------------------
int PKriging::computeDistances(double **XDists, int *length)
{
   int    ii, jj, kk, count;
   double *LDists, dist;

   LDists = new double[(nSamples_*(nSamples_-1)/2)*nInputs_];
   checkAllocate(LDists,"LDists in PKriging::computeDistances");
   count = 0;
   for (jj = 0; jj < nSamples_; jj++)
   {
      for (kk = jj+1; kk < nSamples_; kk++)
      {
         dist = 0.0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            LDists[count*nInputs_+ii] = XNormalized_[jj*nInputs_+ii] - 
                                XNormalized_[kk*nInputs_+ii];
            if (LDists[count*nInputs_+ii] < 0) 
               LDists[count*nInputs_+ii] = - LDists[count*nInputs_+ii]; 
            dist += pow(LDists[count*nInputs_+ii], 2.0);
         }
         if (dist == 0.0)
         {
            printf("PKriging ERROR: repeated sample points.\n");
            printf("               Prune repeated points and re-run.\n");
            printf("Sample %d : (with sample point %d)\n", kk+1, jj+1);
            for (ii = 0; ii < nInputs_; ii++)
               printf("   Input %d : %e\n",ii+1,
                      XNormalized_[kk*nInputs_+ii]*XStds_[ii] +XMeans_[ii]);
            exit(1);
         }
         count++;
      }
   }
   (*length) = count;
   (*XDists) = LDists;
   return 0;
}

