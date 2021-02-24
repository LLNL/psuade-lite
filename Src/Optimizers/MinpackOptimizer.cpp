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
// Functions for the class MinpackOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************

#ifdef HAVE_MINPACK

#include <stdio.h>
#include <stdlib.h>

#include "MinpackOptimizer.h"
#include "Psuade.h"

#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// local variables
// ------------------------------------------------------------------------

int     M_nOutputs=0;
int     M_nFuncEvals=0;
int     M_outputLevel=0;
double *M_optimalX=NULL;
double  M_optimalY=0.0;
double *M_lowerBounds=NULL;
double *M_upperBounds=NULL;
double  M_deltaX=1.0e-8;
FunctionInterface *M_funcIO=NULL;
#define psMinpackMaxSaved_ 1000000
int     psMinpackNSaved_=0;
double  psMinpackSaveX_[psMinpackMaxSaved_];
double  psMinpackSaveY_[psMinpackMaxSaved_];
int     psMCurrDriver_=-1;

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------

extern "C" 
{
void *M_evaluateFunction(int *nOutputs,int *nInputs,double *XValues, 
                       double *YValues,double *YJac,int *ldYJac,int *flag)
{
   int    ii, jj, kk, mm, numLocalY, evalFlag, ignoreFlag, *statuses;
   int    count, *evalTrack, nInps, *sampleIDs;
   double *localY, *localX2, *localY2, deltaX, YData, dtemp;
   FILE   *infile;

   nInps = (*nInputs);
   if (psOptExpertMode_ != 0 && psMinpackNSaved_ == 0)
   {
      infile = fopen("psuade_minpack_data","r");
      if (infile != NULL)
      {
         fscanf(infile, "%d %d %d", &psMinpackNSaved_, &kk, &mm);
         if ((psMinpackNSaved_ <= 0) ||
             (psMinpackNSaved_+1+nInps) > psMinpackMaxSaved_/nInps)
         {
            printf("PSUADE Minpack: history file has too much data.\n");
            printf("                Please consult PSUADE developers.\n");
            fclose(infile);
         }
         else if (kk != nInps)
         {
            printf("PSUADE Minpack: history file has invalid nInputs.\n");
            fclose(infile);
         }
         else if (mm != M_nOutputs)
         {
            printf("PSUADE Minpack: history file has invalid nOutputs.\n");
            fclose(infile);
         }
         else
         {
            for (ii = 0; ii < psMinpackNSaved_; ii++)
            {
               fscanf(infile, "%d", &kk);
               if (kk != ii+1)
               {
                  printf("PSUADE Minpack: data index mismatch.\n");
                  fclose(infile);
                  psMinpackNSaved_ = 0;
               }
               for (jj = 0; jj < nInps; jj++)
                  fscanf(infile, "%lg", &psMinpackSaveX_[ii*nInps+jj]);
               for (jj = 0; jj < M_nOutputs; jj++)
                  fscanf(infile, "%lg", &psMinpackSaveY_[ii*M_nOutputs+jj]);
            }
            fclose(infile);
         }
      }
   }

   numLocalY = M_nOutputs;
   localY    = new double[numLocalY];
   evalTrack = new int[nInps + 1];
   sampleIDs = new int[nInps + 1];
   for (ii = 0; ii < nInps+1; ii++) evalTrack[ii] = 0;

   evalFlag = 1;
   for (ii = 0; ii < psMinpackNSaved_; ii++)
   {
      for (jj = 0; jj < nInps; jj++)
         if (PABS(psMinpackSaveX_[ii*nInps+jj]-XValues[jj])>1e-14) break;
      if (jj == nInps)
      {
         evalFlag = 0;
         break;
      }
   }

   statuses = new int[nInps+1];
   for (jj = 0; jj <= nInps; jj++) statuses[jj] = 1;
   if (evalFlag == 1)
   {
      statuses[0] = M_funcIO->evaluate(M_nFuncEvals,nInps,XValues,
                                       numLocalY,localY,0);
      sampleIDs[0] = M_nFuncEvals;
      M_nFuncEvals++;
      evalTrack[0] = 1;
   }
   else
   {
      statuses[0] = 0;
      for (jj = 0; jj < M_nOutputs; jj++)
         localY[jj] = psMinpackSaveY_[ii*M_nOutputs+jj];
   }

   if ((*flag) == 2)
   {
      localX2 = new double[nInps*nInps];
      localY2 = new double[(*nOutputs)*nInps];
      for (ii = 0; ii < nInps; ii++)
      {
         for (jj = 0; jj < nInps; jj++)
            localX2[ii*nInps+jj] = XValues[jj];
         deltaX = (M_upperBounds[ii] - M_lowerBounds[ii]) * M_deltaX;
         localX2[ii*nInps+ii] += deltaX;
         if (localX2[ii*nInps+ii] > M_upperBounds[ii])
         {
            deltaX = - deltaX;
            localX2[ii*nInps+ii] += 2.0 * deltaX;
         }

         evalFlag = 1;
         for (kk = 0; kk < psMinpackNSaved_; kk++)
         {
            for (jj = 0; jj < nInps; jj++)
            {
               dtemp = psMinpackSaveX_[kk*nInps+jj]-localX2[ii*nInps+jj];
               if (PABS(dtemp) > 1e-14) break;
            }
            if (jj == nInps)
            {
               evalFlag = 0;
               break;
            }
         }

         if (evalFlag == 1)
         {
            statuses[ii+1] = M_funcIO->evaluate(M_nFuncEvals,
                                   nInps,&localX2[ii*nInps],
                                   numLocalY,&localY2[ii*M_nOutputs],0);
            sampleIDs[ii+1] = M_nFuncEvals;
            M_nFuncEvals++;
            evalTrack[ii+1] = 1;
         }
         else
         {
            statuses[ii+1] = 0;
            for (jj = 0; jj < M_nOutputs; jj++)
               localY2[ii*M_nOutputs+jj] = psMinpackSaveY_[kk*M_nOutputs+jj];
         }
      }
   }

   while (statuses[0] != 0)
   {
      statuses[0] = M_funcIO->evaluate(sampleIDs[0],nInps,XValues,
                                       numLocalY,localY,2);
   } 
   if ((*flag) == 2)
   {
      count = nInps;
      while (count > 0)
      {
         count = nInps;
         for (ii = 1; ii <= nInps; ii++) 
         {
            if (statuses[ii] != 0)
               statuses[ii] = M_funcIO->evaluate(sampleIDs[ii],
                                      nInps,&localX2[(ii-1)*nInps],
                                      numLocalY,&localY2[(ii-1)*M_nOutputs],2);
            if (statuses[ii] == 0) count--;
         }
      } 
   }
   delete [] statuses;
         
   for (ii = 0; ii < nInps; ii++)
   {
      if (XValues[ii] < M_lowerBounds[ii])
         for (ii = 0; ii < M_nOutputs; ii++) 
            localY[ii] += 1.0e3 * pow(M_lowerBounds[ii] - XValues[ii],2.0);
      if (XValues[ii] > M_upperBounds[ii])
         for (ii = 0; ii < M_nOutputs; ii++) 
            localY[ii] += 1.0e3 * pow(XValues[ii] - M_upperBounds[ii],2.0);
   }
   for (ii = 0; ii < numLocalY; ii++) YValues[ii] = localY[ii];

   YData = 0.0;
   for (ii = 0; ii < numLocalY; ii++) YData += (localY[ii] * localY[ii]);

   ignoreFlag = 0;
   for (ii = 0; ii < nInps; ii++)
   {
      if (XValues[ii] < M_lowerBounds[ii]) ignoreFlag = 1;
      if (XValues[ii] > M_upperBounds[ii]) ignoreFlag = 1;
      if (ignoreFlag == 1) break;
   }
   if (ignoreFlag == 0 && YData < M_optimalY)
   {
      M_optimalY = YData;
      for (ii = 0; ii < nInps; ii++) M_optimalX[ii] = XValues[ii];
   }

   if (M_outputLevel > 1)
   {
      printf("==================================================>\n");
      count = 0;
      for (ii = 0; ii <= nInps; ii++) count += evalTrack[ii];
      count = M_nFuncEvals - count;
      if (evalTrack[0] == 1) count++;
      for (ii = 0; ii < nInps; ii++)
         printf("MinpackOptimizer(1) %6d : X%2d = %16.8e\n", 
                count, ii+1, XValues[ii]);
      YData = 0.0;
      for (ii = 0; ii < numLocalY; ii++)
         YData += (localY[ii] * localY[ii]);
      printf("MinpackOptimizer(1) %6d :        Y    = %16.8e\n", 
             count, YData);
      if ((*flag) == 2)
      {
         for (kk = 0; kk < nInps; kk++)
         {
            if (evalTrack[kk+1] == 1) count++;
            for (ii = 0; ii < nInps; ii++)
               printf("MinpackOptimizer(2) %6d : X%2d = %16.8e\n", 
                      count, ii+1, localX2[kk*nInps+ii]);
            YData = 0.0;
            for (ii = 0; ii < numLocalY; ii++)
               YData += pow(localY2[kk*M_nOutputs+ii], 2.0e0);
            printf("MinpackOptimizer(2) %6d :        Y    = %16.8e\n", 
                    count, YData);
         }
      }
      printf("<==================================================\n");
   }

   if (psOptExpertMode_ != 0 && evalFlag == 1)
   {
      if (evalTrack[0] == 1)
      {
         for (jj = 0; jj < nInps; jj++)
            psMinpackSaveX_[psMinpackNSaved_*nInps+jj] = XValues[jj];
         for (jj = 0; jj < M_nOutputs; jj++) 
            psMinpackSaveY_[psMinpackNSaved_*M_nOutputs+jj] = localY[jj];
         psMinpackNSaved_++;
      }
      if ((*flag) == 2)
      {
         for (jj = 0; jj < nInps; jj++)
         {
            if (evalTrack[jj+1] == 1)
            {
               for (kk = 0; kk < nInps; kk++)
                  psMinpackSaveX_[psMinpackNSaved_*nInps+kk] = 
                                             localX2[jj*nInps+kk];
               for (kk = 0; kk < M_nOutputs; kk++) 
                  psMinpackSaveY_[psMinpackNSaved_*M_nOutputs+kk] = 
                                             localY2[jj*M_nOutputs+kk];
               psMinpackNSaved_++;
            }
         }
      }
   }

   if ((*flag) == 2)
   {
      for (ii = 0; ii < nInps; ii++)
      {
         for (jj = 0; jj < M_nOutputs; jj++) 
            YJac[(*ldYJac)*ii+jj] = 
                  (localY2[ii*M_nOutputs+jj] - localY[jj]) / deltaX;
      }
      delete [] localX2;
      delete [] localY2;
   }

   if (psOptExpertMode_ != 0 && evalFlag == 1)
   {
      infile = fopen("psuade_minpack_data","w");
      if (infile != NULL)
      {
         fprintf(infile, "%d %d %d\n", psMinpackNSaved_, nInps, M_nOutputs);
         for (ii = 0; ii < psMinpackNSaved_; ii++)
         {
            fprintf(infile, "%d ", ii+1);
            for (jj = 0; jj < nInps; jj++)
               fprintf(infile, "%24.16e ", psMinpackSaveX_[ii*nInps+jj]);
            for (jj = 0; jj < M_nOutputs; jj++)
               fprintf(infile, "%24.16e\n", psMinpackSaveY_[ii*M_nOutputs+jj]);
         }
         fclose(infile);
      }
   }

   M_funcIO->setSynchronousMode();
   delete [] localY;
   delete [] evalTrack;
   delete [] sampleIDs;
   return NULL;
}
} /* extern C */

// ************************************************************************
// external function
// ------------------------------------------------------------------------
extern "C"
{
void lmder1_(void *,int *, int *, double *, double *, double *,
             int *, double *, int *, int *, double *, int *);
}

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MinpackOptimizer::MinpackOptimizer()
{
   M_optimalX    = NULL;
   M_lowerBounds = NULL;
   M_upperBounds = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MinpackOptimizer::~MinpackOptimizer()
{
   if (M_optimalX    != NULL) delete [] M_optimalX;
   if (M_lowerBounds != NULL) delete [] M_lowerBounds;
   if (M_upperBounds != NULL) delete [] M_upperBounds;
   M_optimalX    = NULL;
   M_lowerBounds = NULL;
   M_upperBounds = NULL;
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
// In order to use this optimizer, the number of output of interest should
// the number of terms in the least squares. 
// ------------------------------------------------------------------------
void MinpackOptimizer::optimize(oData *odata)
{
   int    nInputs, nFuns, lws, *iws, info, i, flag=0;
   double tol, *dws, *XValues, *YValues, *YJac, *initialX, ddata1, ddata2;

   if (M_optimalX    != NULL) delete [] M_optimalX;
   if (M_lowerBounds != NULL) delete [] M_lowerBounds;
   if (M_upperBounds != NULL) delete [] M_upperBounds;
   M_deltaX      = odata->deltaX_;
   M_optimalX    = NULL;
   M_lowerBounds = NULL;
   M_upperBounds = NULL;
   M_outputLevel = odata->outputLevel_;
   nInputs       = odata->nInputs_;
   M_nOutputs    = odata->nOutputs_;
   if (nInputs > M_nOutputs)
   {
      printf("MinpackOptimizer ERROR: nInputs should be <= nOutputs.\n");
      printf("    meaning that the number of terms in the least squares\n");
      printf("    should be at least as many as the number of inputs.\n");
      exit(1);
   }
   tol           = odata->tolerance_;
   M_optimalX    = new double[nInputs];
   M_lowerBounds = new double[nInputs];
   M_upperBounds = new double[nInputs];
   M_optimalY    = 1.0e50;
   M_funcIO      = odata->funcIO_;
   M_nFuncEvals  = 0;
   initialX      = odata->initialX_;
   for (i = 0; i < nInputs; i++)
   {
      M_lowerBounds[i] = odata->lowerBounds_[i];
      M_upperBounds[i] = odata->upperBounds_[i];
   }

   if ((odata->setOptDriver_ & 1))
   {
      printf("Minpack: setting optimization simulation driver.\n");
      psMCurrDriver_ = M_funcIO->getDriver();
      M_funcIO->setDriver(1);
   }
   if (odata->maxParallelJobs_ > 1) M_funcIO->setAsynchronousMode();

   XValues = new double[nInputs];
   YValues = new double[M_nOutputs];
   YJac    = new double[M_nOutputs*nInputs];
   lws     = 20 * nInputs;
   iws     = new int[lws];
   dws     = new double[lws];

   for (i = 0; i < nInputs; i++) XValues[i] = initialX[i];

   nFuns = M_nOutputs;
   i = 0;
   M_evaluateFunction(&M_nOutputs, &nInputs, XValues, YValues,
                      YJac, &i, &flag);
   ddata1 = pow(YValues[0], 2.0);
   for (i = 1; i < M_nOutputs; i++) ddata1 += pow(YValues[i], 2.0e0);
   ddata1 = sqrt(ddata1);

   ddata2 = M_upperBounds[0] - M_lowerBounds[0];
   for (i = 1; i < nInputs; i++) 
      if ((M_upperBounds[i]-M_lowerBounds[i]) < ddata2)
         ddata2 = M_upperBounds[i] - M_lowerBounds[i];

   if (ddata1 < ddata2) tol *= ddata1;
   else                 tol *= ddata2;

   printf("MinpackOptimizer tolerance = %e (Please check).\n", tol);
   
   lmder1_((void *) M_evaluateFunction, &nFuns, &nInputs, XValues,
           YValues, YJac, &nFuns, &tol, &info, iws, dws, &lws);

   odata->optimalY_ = M_optimalY;
   odata->numFuncEvals_ = M_nFuncEvals;
   for (i = 0; i < nInputs; i++) odata->optimalX_[i] = M_optimalX[i];

   if ((odata->setOptDriver_ & 2) && psMCurrDriver_ >= 0)
   {
      printf("Minpack: resetting optimization simulation driver.\n");
      M_funcIO->setDriver(psMCurrDriver_);
   }
   if (odata->maxParallelJobs_ > 1) M_funcIO->setSynchronousMode();

   delete [] XValues;
   delete [] YValues;
   delete [] YJac;
   delete [] iws;
   delete [] dws;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
MinpackOptimizer& MinpackOptimizer::operator=(const MinpackOptimizer &)
{
   printf("MinpackOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}
#else
   int minpack_bogus=0;
#endif


