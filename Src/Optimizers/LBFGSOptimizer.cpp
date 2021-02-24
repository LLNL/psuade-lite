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
// Functions for the class LBFGSOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "LBFGSOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
int psLBCurrDriver_=-1;

#ifdef HAVE_LBFGS
extern "C" {
#include "../../External/L-BFGS-B-C/src/lbfgsb.h"
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
LBFGSOptimizer::LBFGSOptimizer()
{
  if (isScreenDumpModeOn())
  {
    printAsterisks(PL_INFO, 0);
    printf("*   LBFGS Optimizer Usage Information\n");
    printEquals(PL_INFO, 0);
    printf("* - To run this optimizer, first make sure opt_driver has\n");
    printf("*   been initialized to point to your optimization objective\n");
    printf("*   function evaluator\n");
    printf("* - Set optimization tolerance in your PSUADE input file\n");
    printf("* - Set maximum number of iterations in PSUADE input file\n");
    printf("* - Set num_local_minima to perform multistart optimization\n");
    printf("* - Set optimization print_level to give additonal outputs\n");
    printf("* - In Opt EXPERT mode, the optimization history log will be\n");
    printf("*   turned on automatically. Previous psuade_lbfgs_history\n");
    printf("*   file will also be reused.\n");
    printf("* - Your opt_driver must have the number of outputs equal to\n");
    printf("*   the number of inputs plus 1 whereby the last nInputs\n");
    printf("*   outputs are the derivatives of the outputs with respect\n");
    printf("*   to each input.\n");
    printAsterisks(PL_INFO, 0);
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
LBFGSOptimizer::~LBFGSOptimizer()
{
}

// ************************************************************************
// optimize (especially for library call mode)
// ------------------------------------------------------------------------
void LBFGSOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                              double *ubds, int maxfun, double tol)
{
   double *optimalX;
   if (nInputs <= 0)
   {
      printf("LBFGSOptimizer ERROR: nInputs <= 0.\n");
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
   odata->nOutputs_ = nInputs + 1;
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
   for (int ii = 0; ii < nInputs; ii++)
      XValues[ii] = optimalX_[ii] = odata->optimalX_[ii];
   delete [] optimalX;
   odata->optimalX_ = NULL;
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void LBFGSOptimizer::optimize(oData *odata)
{
#if HAVE_LBFGS
   int     funcID, ii, kk, nOuts, its, printLevel, nOne=1;
   double  *XValues, *lbounds, *ubounds, *GValues, FValue, *SValues;
   integer nInps, iprint=0, itask, *task=&itask, lsave[4], isave[44];
   integer *iwork, nCorr=5, *nbds, csave[60];
   double  factr, pgtol, *work, dsave[29]; 

   printLevel = odata->outputLevel_;
   nInps = odata->nInputs_;
   nOuts = odata->nOutputs_;
   if (nInps+1 != nOuts)
   {
      printf("WARNING: LBFGS optimizer needs nOuts = nInps+1.\n");
      if (nOuts != 1) exit(1);
      else
         printf("   INFO: will use finite difference.\n");
   }
   lbounds = odata->lowerBounds_;
   ubounds = odata->upperBounds_;
   XValues = new double[nInps];
   for (ii = 0; ii < nInps; ii++) XValues[ii] = odata->initialX_[ii];
   funcID = odata->numFuncEvals_;

   odata->optimalY_ = 1.0e50;
   if ((odata->setOptDriver_ & 1))
   {
      printf("LBFGS: setting optimization simulation driver.\n");
      psLBCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   GValues = new double[nInps];
   SValues = new double[nInps+1];
   for (ii = 0; ii < nInps; ii++) GValues[ii] = 0.0;
   for (ii = 0; ii < nInps+1; ii++) SValues[ii] = 0.0;
   nbds    = new integer[nInps];
   for (ii = 0; ii < nInps; ii++) nbds[ii] = 2;
   factr = 1e7;
   pgtol = 1e-6;
   kk = (2 * nCorr + 5) * nInps + 11 * nCorr * nCorr + 8 * nCorr;
   work  = new double[kk];
   iwork = new integer[3*nInps];
   *task = (integer) START;
   its   = 0;
 
   while (1)
   {
      its++;
      setulb(&nInps, &nCorr, XValues, lbounds, ubounds, nbds, &FValue,
          GValues, &factr, &pgtol, work, iwork, task, &iprint, csave,
          lsave, isave, dsave);
      if (IS_FG(*task))
      {
         if (isScreenDumpModeOn())
         {
            if (printLevel > 1)
               for (ii = 0; ii < nInps; ii++)
                  printf("   Current Input X %4d = %24.16e\n",ii+1,
                         XValues[ii]);
         }
         if (nOuts == (nInps+1))
         {
            if (odata->optFunction_ != NULL)
               odata->optFunction_(nInps, XValues, nOuts, SValues);
            else if (odata->funcIO_ != NULL)
               odata->funcIO_->evaluate(funcID,nInps,XValues,nOuts,SValues,0);
            else
            {
               printf("LBFGSOptimizer ERROR: no function evaluator.\n");
               exit(1);
            }
            funcID = odata->numFuncEvals_++;
         }
         else
         {
            if (odata->optFunction_ != NULL)
               odata->optFunction_(nInps, XValues, nOne, SValues);
            else if (odata->funcIO_ != NULL)
               odata->funcIO_->evaluate(funcID,nInps,XValues,nOne,SValues,0);
            funcID = odata->numFuncEvals_++;
            for (ii = 0; ii < nInps; ii++) 
            {
               XValues[ii] += 1.0e-8;
               if (odata->optFunction_ != NULL)
                  odata->optFunction_(nInps, XValues, nOne, &SValues[ii+1]);
               else if (odata->funcIO_ != NULL)
                  odata->funcIO_->evaluate(funcID,nInps,XValues,nOne,
                                           &SValues[ii+1],0);
               SValues[ii+1] = 1e8 * (SValues[ii+1] - SValues[0]);
               XValues[ii] -= 1.0e-8;
               funcID = odata->numFuncEvals_++;
            }
         }
         FValue = SValues[0];
         for (ii = 0; ii < nInps; ii++) GValues[ii] = SValues[ii+1];
         if (odata->outputLevel_ > 4)
         {
            printf("LBFGSOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInps; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            printf("    Y = %16.8e\n", FValue);
            for (ii = 1; ii < nInps; ii++)
               printf("    dY/dX %6d = %16.8e\n", ii, SValues[ii]);
         }
      }
      else if (*task != NEW_X)
      {
         if (isScreenDumpModeOn())
         {
            if (printLevel > 0)
               for (ii = 0; ii < nInps; ii++)
                  printf("Final Input X %4d = %24.16e\n", ii+1, XValues[ii]);
            printf("Final objective function value       = %16.8e\n",FValue);
            printf("Total number of function evaluations = %d\n", its);
         }
         break;
      }
   }

   if ((odata->setOptDriver_ & 2) && psLBCurrDriver_ >= 0)
   {
      printf("LBFGS: setting back to original simulation driver.\n");
      odata->funcIO_->setDriver(psLBCurrDriver_);
   }
   for (ii = 0; ii < nInps; ii++) odata->optimalX_[ii] = XValues[ii];
   odata->optimalY_ = FValue;
   if (optimalX_ != NULL) delete [] optimalX_;
   optimalX_ = new double[nInps];
   for (ii = 0; ii < nInps; ii++) optimalX_[ii] = XValues[ii];
   optimalY_ = odata->optimalY_;
   numEvals_ = odata->numFuncEvals_;
   delete [] work;
   delete [] iwork;
   delete [] nbds;
   delete [] XValues;
   delete [] GValues;
#else
   printf("ERROR : LBFGS optimizer not installed.\n");
   exit(1);
#endif
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
LBFGSOptimizer& LBFGSOptimizer::operator=(const LBFGSOptimizer &)
{
   printf("LBFGSOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

