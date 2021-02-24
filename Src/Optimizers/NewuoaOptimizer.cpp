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
// Functions for the class NewuoaOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>

#include "NewuoaOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

// ************************************************************************
// External functions
// ------------------------------------------------------------------------
extern "C" 
{
  void newuoa_(int *,int *,double *,double *,double *,int *,int *,double*);
}

// ************************************************************************
// Internal global variables
// ------------------------------------------------------------------------
#define psNewuoaMaxSaved_ 10000
int     psNewuoaNSaved_=0;
double  psNewuoaSaveX_[psNewuoaMaxSaved_*10];
double  psNewuoaSaveY_[psNewuoaMaxSaved_*10];
void    *psNewuoaObj_=NULL;
int     psNCurrDriver_=-1;
int     psNewuoaSaveHistory_=0;
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" 
{
#endif
   void *newuoaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, jj, kk, funcID, nInputs, nOutputs, outputID, found;
      double *localY, ddata;
      char   pString[1000], lineIn[1000];
      oData  *odata;
      FILE   *infile;

      nInputs = (*nInps);
      odata    = (oData *) psNewuoaObj_;
      nOutputs = odata->nOutputs_;
      localY   = (double *) malloc(nOutputs * sizeof(double));
      outputID = odata->outputID_;

      found = 0;
      for (ii = 0; ii < psNewuoaNSaved_; ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            if (PABS(psNewuoaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) 
               break;
         if (jj == nInputs)
         {
            found = 1;
            printf("Newuoa: simulation results reuse.\n");
            break;
         }
      }

      funcID = odata->numFuncEvals_;
      if (found == 0)
      {
         if (odata->optFunction_ != NULL)
            odata->optFunction_(nInputs, XValues, nOutputs, localY);
         else if (odata->funcIO_ != NULL)
            odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
         else
         {
            printf("NewuoaOptimizer ERROR: no function evaluator.\n");
            exit(1);
         }
         if (odata->outputLevel_ > 4)
         {
            printf("NewuoaOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            printf("    Y = %16.8e\n", localY[outputID]);
         }
         funcID = odata->numFuncEvals_++;
         (*YValue) = localY[outputID];
         if (odata->outputLevel_ > 4)
            printf("    ==> Objective function value = %16.8e\n", (*YValue));
      }
      else
      {
         localY[outputID] = psNewuoaSaveY_[ii];
         (*YValue) = localY[outputID];
      }

      if ((psNewuoaSaveHistory_ == 1 && found == 0 &&
          (psNewuoaNSaved_+1)*nInputs < psNewuoaMaxSaved_*10) && 
          psNewuoaNSaved_ < psNewuoaMaxSaved_)
      {
         for (jj = 0; jj < nInputs; jj++)
            psNewuoaSaveX_[psNewuoaNSaved_*nInputs+jj] = XValues[jj];
         psNewuoaSaveY_[psNewuoaNSaved_] = (*YValue);
         psNewuoaNSaved_++;
         infile = fopen("psuade_newuoa_history","w");
         if (infile != NULL)
         {
            for (ii = 0; ii < psNewuoaNSaved_; ii++)
            {
               fprintf(infile, "999 %d ", nInputs);
               for (kk = 0; kk < nInputs; kk++)
                  fprintf(infile, "%24.16e ", psNewuoaSaveX_[ii*nInputs+kk]);
               fprintf(infile, "%24.16e\n", psNewuoaSaveY_[ii]);
            }
            fclose(infile);
         }
      }

      if ((*YValue) < odata->optimalY_)
      {
         odata->optimalY_ = (*YValue);
         for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
         if (psOptExpertMode_ == 1 && odata->outputLevel_ > 1)
         {
            printf("NewuoaOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            if (found == 0) printf("    Y = %16.8e\n", localY[outputID]);
            printf("    *** Current best objective function value = %16.8e\n", 
                   odata->optimalY_);
         }
      }
      free(localY);
      return NULL;
   }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
NewuoaOptimizer::NewuoaOptimizer()
{
  if (isScreenDumpModeOn())
  {
    printAsterisks(PL_INFO, 0);
    printf("*    NEWUOA Optimizer Usage Information\n");
    printEquals(PL_INFO, 0);
    printf("* NEWUOA is an optimization software developed by Professor\n");
    printf("* Michael Powell to solve unconstrained optimization problems\n");
    printf("* without derrivatives suitable for problems with hundreds\n");
    printf("* of variables. The formulation of the problem is as follow:\n");
    printf("\n*                min_X F(X)\n");
    printf("*\n");
    printf("* To run this optimizer, do the following: \n");
    printf("* (1) Prepare a PSUADE input file with all variable defined.\n");
    printf("* (2) In the PSUADE input file, make sure opt_driver has been\n");
    printf("*     initialized to point your optimization objective function\n");
    printf("*      evaluator\n");
    printf("* (3) Since it is unconstrained, you must provide an initial\n");
    printf("      guess in the PSUADE_IO section, or you can generate an\n");
    printf("      initial guess using one of the sampling methods.\n");
    printf("* (4) Set optimization tolerance in PSUADE input file\n");
    printf("* (5) Set maximum number of iterations in PSUADE input file\n");
    printf("* (6) Set optimization print_level to give additonal outputs\n");
    printf("* (7) In Opt EXPERT mode, the optimization history log will be\n");
    printf("*     turned on automatically. Previous psuade_newuoa_history\n");
    printf("*     file will also be reused to save evaluations.\n");
    printf("* (8) If your opt_driver is a response surface which has more\n");
    printf("*     input than the number of optimization inputs, you can fix\n");
    printf("*     some driver inputs by using a rs_index_file (ANALYSIS).\n");
    printEquals(PL_INFO, 0);
    printf("To reuse the results (e.g. restart from abrupt termination),\n");
    printf("turn on save_history and use_history optimization options in\n");
    printf("the ANALYSIS section (e.g. optimization save_history). You will\n");
    printf("see a file created called 'psuade_newuoa_history' afterward.\n");
    printAsterisks(PL_INFO, 0);
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
NewuoaOptimizer::~NewuoaOptimizer()
{
}

// ************************************************************************
// optimize (this function should be used in the library call mode)
// ------------------------------------------------------------------------
void NewuoaOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                               double *ubds, int maxfun, double tol)
{
   double *optimalX;
   if (nInputs <= 0)
   {
      printf("NewuoaOptimizer ERROR: nInputs <= 0.\n");
      exit(1);
   }
   oData *odata = new oData();
   odata->outputLevel_ = 0;
   odata->nInputs_ = nInputs;
   optimalX = new double[nInputs];
   odata->optimalX_ = optimalX;
   odata->initialX_ = XValues;
   odata->lowerBounds_ = lbds;
   odata->upperBounds_ = ubds;
   odata->tolerance_ = tol;
   if (odata->tolerance_ <= 0) odata->tolerance_ = 1e-6;
   odata->nOutputs_ = 1;
   odata->outputID_ = 0;
   odata->maxFEval_ = maxfun;
   odata->numFuncEvals_ = 0;
   odata->tolerance_ = tol;
   odata->setOptDriver_ = 0;
   odata->optFunction_ = objFunction_;
   odata->funcIO_ = NULL;
   optimize(odata);
   odata->initialX_ = NULL;
   odata->lowerBounds_ = NULL;
   odata->upperBounds_ = NULL;
   odata->optFunction_ = NULL;
   delete [] optimalX;
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void NewuoaOptimizer::optimize(oData *odata)
{
   int    nInputs, printLevel=0, ii, kk, maxfun, nPts=0;
   int    nOutputs, outputID, newuoaFlag=0;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp, *workArray;
   char   cinput[5000], *cString, winput[5001];;
   FILE   *infile=NULL;

   printLevel = odata->outputLevel_;
   nInputs = odata->nInputs_;
   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
   odata->optimalY_ = 1.0e50;
   XValues = new double[nInputs+1];
   for (ii = 0; ii < nInputs; ii++) XValues[ii] = odata->initialX_[ii];
   rhobeg = odata->upperBounds_[0] - odata->lowerBounds_[0];
   for (ii = 1; ii < nInputs; ii++) 
   {
      dtemp = odata->upperBounds_[ii] - odata->lowerBounds_[ii];
      if (dtemp < rhobeg) rhobeg = dtemp;
   }
   rhobeg *= 0.5;
   rhoend = rhobeg * odata->tolerance_;
   if (rhobeg < rhoend)
   {
      printf("NewuoaOptimizer WARNING: tolerance too large.\n");
      printf("                         tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }
   nOutputs = odata->nOutputs_;
   outputID = odata->outputID_;
   maxfun = odata->maxFEval_;
   if ((odata->setOptDriver_ & 1))
   {
      printf("Newuoa: setting optimization simulation driver.\n");
      psNCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   psNewuoaObj_= (void *) odata;
   if (isScreenDumpModeOn())
   {
      printAsterisks(PL_INFO, 0);
      printf("Newuoa optimizer: max fevals = %d\n", odata->maxFEval_);
      printf("Newuoa optimizer: tolerance  = %e\n", odata->tolerance_);
      if (printLevel > 1)
         printf("Newuoa optimizer: rho1, rho2 = %e %e\n", rhobeg, rhoend);
      printEquals(PL_INFO, 0);
   }

   nPts = (nInputs + 1) * (nInputs + 2) / 2;
   kk   = (nPts + 13) * (nPts + nInputs) + 3*nInputs*(nInputs+3)/2;
   workArray = new double[kk];

#ifdef HAVE_NEWUOA
   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("opt_save_history");
      if (cString != NULL) psNewuoaSaveHistory_ = 1;
      cString = psConfig_->getParameter("opt_use_history");
      if (cString != NULL)
      {
         printf("Newuoa: use history has been turned on.\n");
         infile = fopen("psuade_newuoa_history","r");
         if (infile != NULL)
         {
            psNewuoaNSaved_ = 0;
            while (feof(infile) == 0)
            {
               fscanf(infile, "%d %d", &ii, &kk);
               if (ii != 999 || kk != nInputs)
               {
                  break;
               }
               else
               {
                  for (ii = 0; ii < nInputs; ii++) 
                     fscanf(infile, "%lg",
                        &psNewuoaSaveX_[psNewuoaNSaved_*nInputs+ii]);
                  fscanf(infile, "%lg",&psNewuoaSaveY_[psNewuoaNSaved_]);
                  psNewuoaNSaved_++;
               }
               if (((psNewuoaNSaved_+1)*nInputs > psNewuoaMaxSaved_*10) ||
                   psNewuoaNSaved_ > psNewuoaMaxSaved_) break;
            } 
            fclose(infile);
         }
      }
   }
   if (isScreenDumpModeOn())
   {
      for (ii = 0; ii < nInputs; ii++) 
         printf("Newuoa initial X %3d = %e\n", ii+1, XValues[ii]);
   }
   newuoaFlag = 7777;
#ifdef HAVE_NEWUOA
   newuoa_(&nInputs, &nPts, XValues, &rhobeg, &rhoend, &newuoaFlag, 
           &maxfun, workArray);
#else
      printf("Newuoa optimizer ERROR: newuoa not installed.\n");
      exit(1);
#endif

   if (isScreenDumpModeOn())
   {
      printf("Newuoa optimizer: total number of evaluations = %d\n",
              odata->numFuncEvals_);
   }
   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
   if (optimalX_ != NULL) delete [] optimalX_;
   optimalX_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) optimalX_[ii] = odata->optimalX_[ii];
   optimalY_ = odata->optimalY_;
   numEvals_ = odata->numFuncEvals_;

   if (psNewuoaSaveHistory_ == 1 && psNewuoaNSaved_ > 0)
   {
      infile = fopen("psuade_newuoa_history","w");
      if (infile != NULL)
      {
         for (ii = 0; ii < psNewuoaNSaved_; ii++)
         {
            fprintf(infile, "999 %d ", nInputs);
            for (kk = 0; kk < nInputs; kk++)
               fprintf(infile, "%24.16e ", psNewuoaSaveX_[ii*nInputs+kk]);
            fprintf(infile, "%24.16e\n", psNewuoaSaveY_[ii]);
         }
         fclose(infile);
      }
      printf("Newuoa: history saved in psuade_newuoa_history\n");
   }
#else
   printf("ERROR : Newuoa optimizer not installed.\n");
   exit(1);
#endif

   if ((odata->setOptDriver_ & 2) && psNCurrDriver_ >= 0)
   {
      printf("Newuoa: setting back to original simulation driver.\n");
      odata->funcIO_->setDriver(psNCurrDriver_);
   }
   delete [] XValues;
   delete [] workArray;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
NewuoaOptimizer& NewuoaOptimizer::operator=(const NewuoaOptimizer &)
{
   printf("NewuoaOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

