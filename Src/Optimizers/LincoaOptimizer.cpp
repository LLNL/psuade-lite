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
// Functions for the class LincoaOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>

#include "LincoaOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

// ************************************************************************
// External functions
// ------------------------------------------------------------------------
extern "C" void lincoa_(int *,int *, int *, double *, int *, double *, 
                        double *,double *,double *, int *, int *, double*);

// ************************************************************************
// Internal global variables
// ------------------------------------------------------------------------
#define psLincoaMaxSaved_ 10000
int     psLincoaNSaved_=0;
double  psLincoaSaveX_[psLincoaMaxSaved_*10];
double  psLincoaSaveY_[psLincoaMaxSaved_*10];
void    *psLincoaObj_=NULL;
int     psLincoaNConstr_=0;
double  *psLincoaAmat_=NULL;
double  *psLincoaBvec_=NULL;
int     psLCurrDriver_=-1;
int     psLincoaSaveHistory_=0;
char    psLincoaConstraintFile_[1000]="NONE";
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" 
{
#endif
   void *lincoaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, jj, kk, funcID, nInputs, nOutputs, outputID, found;
      double *localY, ddata;
      char   pString[1000], lineIn[1000];
      oData  *odata;
      FILE   *infile;

      nInputs = (*nInps);
      odata    = (oData *) psLincoaObj_;
      nOutputs = odata->nOutputs_;
      localY   = (double *) malloc(nOutputs * sizeof(double));
      outputID = odata->outputID_;

      found = 0;
      for (ii = 0; ii < psLincoaNSaved_; ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            if (PABS(psLincoaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) break;
         if (jj == nInputs)
         {
            found = 1;
            printf("Lincoa: simulation results reuse.\n");
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
            printf("LincoaOptimizer ERROR: no function evaluator.\n");
            exit(1);
         }
         if (odata->outputLevel_ > 4)
         {
            printf("LincoaOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            for (ii = 0; ii < nOutputs; ii++) 
               printf("    Y %6d = %16.8e\n", ii+1, localY[ii]);
         }
         funcID = odata->numFuncEvals_++;
         (*YValue) = localY[outputID];
         if (odata->outputLevel_ > 4)
            printf("    ==> Objective function value = %16.8e\n", (*YValue));
      }
      else
      {
         localY[outputID] = psLincoaSaveY_[ii];
         (*YValue) = localY[outputID];
      }

      if ((psLincoaSaveHistory_ == 1 && found == 0 &&
          (psLincoaNSaved_+1)*nInputs < psLincoaMaxSaved_*10) && 
          psLincoaNSaved_ < psLincoaMaxSaved_)
      {
         for (jj = 0; jj < nInputs; jj++)
            psLincoaSaveX_[psLincoaNSaved_*nInputs+jj] = XValues[jj];
         psLincoaSaveY_[psLincoaNSaved_] = (*YValue);
         psLincoaNSaved_++;
         infile = fopen("psuade_lincoa_history","w");
         if (infile != NULL)
         {
            for (ii = 0; ii < psLincoaNSaved_; ii++)
            {
               fprintf(infile, "999 %d ", nInputs);
               for (kk = 0; kk < nInputs; kk++)
                  fprintf(infile, "%24.16e ", psLincoaSaveX_[ii*nInputs+kk]);
               fprintf(infile, "%24.16e\n", psLincoaSaveY_[ii]);
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
            printf("LincoaOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            if (found == 0)
            {
               for (ii = 0; ii < nOutputs; ii++) 
                  printf("    Y %6d = %16.8e\n", ii+1, localY[ii]);
            }
            printf("    *** Current best objective function value = %16.8e\n", 
                   odata->optimalY_);
            printf("    NOTE: This point may not satisfy the constraints.\n");
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
LincoaOptimizer::LincoaOptimizer()
{
  if (isScreenDumpModeOn())
  {
    printAsterisks(PL_INFO, 0);
    printf("*    LINCOA Optimizer Usage Information\n");
    printEquals(PL_INFO, 0);
    printf("* LINCOA is an optimization software developed by Professor\n");
    printf("* Michael Powell to solve linearly constrained optimization\n");
    printf("* without derrivatives suitable for problems with hundreds\n");
    printf("* of variables. The formulation of the problem is as follow:\n");
    printf("\n*                min_X F(X)\n");
    printf("*\n");
    printf("*            subject to A(:,j) X < b(j) for j = 1, 2, ... m\n");
    printf("* where m is the number of constraints.\n"); 
    printf("*\n");
    printf("* To run this optimizer, do the following: \n");
    printf("* (1) Prepare a PSUADE input file with all variable defined.\n");
    printf("* (2) In the PSUADE input file, make sure opt_driver has been\n");
    printf("*     initialized to point your optimization objective function\n");
    printf("*      evaluator\n");
    printf("* (3) Set optimization tolerance in PSUADE input file\n");
    printf("* (4) Set maximum number of iterations in PSUADE input file\n");
    printf("* (5) Set optimization print_level to give additonal outputs\n");
    printf("* (6) In Opt EXPERT mode, the optimization history log will be\n");
    printf("*     turned on automatically. Previous psuade_lincoa_history\n");
    printf("*     file will also be reused to save evaluations.\n");
    printf("* (7) If your opt_driver is a response surface which has more\n");
    printf("*     input than the number of optimization inputs, you can fix\n");
    printf("*     some driver inputs by using a rs_index_file (ANALYSIS).\n");
    printf("* (8) You also need to specific the linear constraints by\n");
    printf("*     creating a file called psuade_lincoa_constaints in your\n");
    printf("*     current directory having the following format:\n");
    printf("*     (If no such file exists, PSUADE treats it as if there is\n");
    printf("*      no constraint).\n");
    printf("*\n");
    printf("\tPSUADE_BEGIN\n");
    printf("\t<c> <m>         /* c in number of constraints, m - nInputs */\n");
    printf("\t1  <m+1 values> /* constraint 1: m for A(:,1), 1 for b(1) */\n");
    printf("\t2  <m+1 values> /* constraint 2: m for A(:,2), 1 for b(2) */\n");
    printf("\t...\n");
    printf("\tPSUADE_END\n\n");
    printf("where each line corresponds to one constraint.\n");
    printEquals(PL_INFO, 0);
    printf("To reuse the results (e.g. restart from abrupt termination),\n");
    printf("turn on save_history and use_history optimization options in\n");
    printf("the ANALYSIS section (e.g. optimization save_history). You will\n");
    printf("see a file created called 'psuade_lincoa_history' afterward.\n");
    printAsterisks(PL_INFO, 0);
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
LincoaOptimizer::~LincoaOptimizer()
{
}

// ************************************************************************
// set objective function
// ------------------------------------------------------------------------
void LincoaOptimizer::setConstraintFile(char *fname)
{
   strcpy(psLincoaConstraintFile_, fname);
}

// ************************************************************************
// optimize (this function should be used in the library mode)
// ------------------------------------------------------------------------
void LincoaOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                       double *ubds, int nOutputs, int maxfun, double tol)
{
   double *optimalX;
   if (nInputs <= 0)
   {
      printf("LincoaOptimizer ERROR: nInputs <= 0.\n");
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
   odata->nOutputs_ = nOutputs;
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
void LincoaOptimizer::optimize(oData *odata)
{
   int    nInputs, printLevel=0, ii, kk, maxfun, nPts=0, nConstr;
   int    nOutputs, outputID, lincoaFlag=0;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp, *workArray;
   char   cinput[5000], *cString, winput[5001], cfile[2000];
   FILE   *infile=NULL;
   FILE   *fp=NULL;

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
      printf("LincoaOptimizer WARNING: tolerance too large.\n");
      printf("                         tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }
   nOutputs = odata->nOutputs_;
   outputID = odata->outputID_;

   if (strcmp(psLincoaConstraintFile_, "NONE"))
        strcpy(cfile, psLincoaConstraintFile_);
   else strcpy(cfile, "psuade_lincoa_constraints");
   fp = fopen(cfile, "r");
   if (fp == NULL)
   {
      printf("Lincoa WARNING: constraint file %s not found.\n",cfile);
      printf("                Assume no constraint.\n");
      nConstr = 0;
      psLincoaNConstr_ = 0;
      psLincoaAmat_ = NULL;
      psLincoaBvec_ = NULL;
   } 
   else
   {
      fgets(winput, 5000, fp);
      if (!strcmp(winput,"PSUADE_BEGIN")) fgets(winput, 5000, fp);
      fscanf(fp,"%d %d", &nConstr, &kk);
      if (nConstr <= 0)
      {
         printf("Lincoa ERROR: nConstr <= 0.\n");
         fclose(fp);
         return;
      }
      printf("Lincoa: number of constraints read = %d\n", nConstr);
      if (kk != nInputs)
      {
         printf("Lincoa ERROR: nInputs do not match (%d %d).\n",kk,nInputs);
         fclose(fp);
         return;
      }
      psLincoaNConstr_ = nConstr;
      psLincoaAmat_ = new double[nConstr*nInputs];
      psLincoaBvec_ = new double[nConstr];
      kk = 0;
      while (kk < nConstr && feof(fp) == 0)
      {
         fscanf(fp, "%d", &ii);
         if (ii != kk+1)
         {
            printf("Lincoa ERROR: in reading constraint %d .\n", kk+1);
            fclose(fp);
            return;
         }
         for (ii = 0; ii < nInputs; ii++)
            fscanf(fp, "%lg", &psLincoaAmat_[kk*nInputs+ii]);
         fscanf(fp, "%lg", &psLincoaBvec_[kk]);
         if (isScreenDumpModeOn())
         {
            printf("Constr %4d: ",kk+1);
            for (ii = 0; ii < nInputs; ii++)
               printf(" %12.4e ",psLincoaAmat_[kk*nInputs+ii]);
            printf("%12.4e\n", psLincoaBvec_[kk]);
         }
         kk++;
      } 
      fclose(fp); 
      if (kk != nConstr)
      {
         printf("Lincoa ERROR: no. of constraints do not match (%d vs %d)\n",
                kk, nConstr);
         return;
      }
   }
   maxfun = odata->maxFEval_;
   if ((odata->setOptDriver_ & 1))
   {
      printf("Lincoa: setting optimization simulation driver.\n");
      psLCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   psLincoaObj_= (void *) odata;
   if (isScreenDumpModeOn())
   {
      printAsterisks(PL_INFO, 0);
      printf("Lincoa optimizer: max fevals = %d\n", odata->maxFEval_);
      printf("Lincoa optimizer: tolerance  = %e\n", odata->tolerance_);
      if (printLevel > 1)
         printf("Lincoa optimizer: rho1, rho2 = %e %e\n", rhobeg, rhoend);
      printEquals(PL_INFO, 0);
   }

   nPts = (nInputs + 1) * (nInputs + 2) / 2;
   kk = nConstr*(2+nInputs) + nPts*(4+nInputs+nPts) + 
        nInputs*(9+3*nInputs) + 2*nPts;
   workArray = new double[kk];

#ifdef HAVE_LINCOA
   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("opt_save_history");
      if (cString != NULL) psLincoaSaveHistory_ = 1;
      cString = psConfig_->getParameter("opt_use_history");
      if (cString != NULL)
      {
         printf("Lincoa: use history has been turned on.\n");
         infile = fopen("psuade_lincoa_history","r");
         if (infile != NULL)
         {
            psLincoaNSaved_ = 0;
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
                        &psLincoaSaveX_[psLincoaNSaved_*nInputs+ii]);
                  fscanf(infile, "%lg",&psLincoaSaveY_[psLincoaNSaved_]);
                  psLincoaNSaved_++;
               }
               if (((psLincoaNSaved_+1)*nInputs > psLincoaMaxSaved_*10) ||
                   psLincoaNSaved_ > psLincoaMaxSaved_) break;
            } 
            fclose(infile);
         }
      }
   }
   if (isScreenDumpModeOn())
   {
      for (ii = 0; ii < nInputs; ii++) 
         printf("Lincoa initial X %3d = %e\n", ii+1, XValues[ii]);
   }
   lincoaFlag = 9999;
   lincoa_(&nInputs, &nPts, &nConstr, psLincoaAmat_, &psLincoaNConstr_,
           psLincoaBvec_, XValues, &rhobeg, &rhoend, &lincoaFlag, 
           &maxfun, workArray);
   if (isScreenDumpModeOn())
   {
      printf("Lincoa optimizer: total number of evaluations = %d\n",
              odata->numFuncEvals_);
   }
   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
   odata->optimalY_ = workArray[0];
   if (optimalX_ != NULL) delete [] optimalX_;
   optimalX_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) optimalX_[ii] = XValues[ii];
   optimalY_ = odata->optimalY_ = workArray[0];
   numEvals_ = odata->numFuncEvals_;

   if (psLincoaSaveHistory_ == 1 && psLincoaNSaved_ > 0)
   {
      infile = fopen("psuade_lincoa_history","w");
      if (infile != NULL)
      {
         for (ii = 0; ii < psLincoaNSaved_; ii++)
         {
            fprintf(infile, "999 %d ", nInputs);
            for (kk = 0; kk < nInputs; kk++)
               fprintf(infile, "%24.16e ", psLincoaSaveX_[ii*nInputs+kk]);
            fprintf(infile, "%24.16e\n", psLincoaSaveY_[ii]);
         }
         fclose(infile);
      }
      printf("Lincoa: history saved in psuade_lincoa_history\n");
   }
#else
   printf("ERROR : Lincoa optimizer not installed.\n");
   exit(1);
#endif

   if ((odata->setOptDriver_ & 2) && psLCurrDriver_ >= 0)
   {
      printf("Lincoa INFO: reverting to original simulation driver.\n");
      odata->funcIO_->setDriver(psLCurrDriver_);
   }
   delete [] XValues;
   delete [] workArray;
   if (psLincoaAmat_ != NULL)
   {
      delete [] psLincoaAmat_;
      delete [] psLincoaBvec_;
      psLincoaAmat_ = NULL;
      psLincoaBvec_ = NULL;
      psLincoaNConstr_ = 0;
   }
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
LincoaOptimizer& LincoaOptimizer::operator=(const LincoaOptimizer &)
{
   printf("LincoaOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

