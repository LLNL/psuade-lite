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
// Functions for the class BobyqaOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2009
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "BobyqaOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

// ************************************************************************
// External functions
// ------------------------------------------------------------------------
extern "C" void bobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);

// ************************************************************************
// Internal global variables
// ------------------------------------------------------------------------
#define psBobyqaMaxSaved_ 10000
int     psBobyqaNSaved_=0;
double  psBobyqaSaveX_[psBobyqaMaxSaved_*10];
double  psBobyqaSaveY_[psBobyqaMaxSaved_*10];
void    *psBobyqaObj_=NULL;
int     psNumBOVars_=0;
int     psNumBLVars_=0;
int     *psBOVars_=NULL;
int     *psBLVars_=NULL;
int     psBCurrDriver_=-1;
int     psBobyqaSaveHistory_=0;
double  *psBOWghts_=NULL;
double  *psBLWghts_=NULL;
double  *psBLVals_=NULL;
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" 
{
#endif
   void *bobyqaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, jj, kk, funcID, nInputs, nOutputs, outputID, found, iOne=1;
      double *localY, ddata;
      char   pString[1000], lineIn[1000];
      oData  *odata;
      FILE   *infile;

      nInputs  = (*nInps);
      odata    = (oData *) psBobyqaObj_;
      nOutputs = odata->nOutputs_;
      outputID = odata->outputID_;
      localY   = (double *) malloc(nOutputs * sizeof(double));

      found = 0;
      for (ii = 0; ii < psBobyqaNSaved_; ii++)
      {
         for (jj = 0; jj < nInputs; jj++)
            if (PABS(psBobyqaSaveX_[ii*nInputs+jj]-XValues[jj])>1.0e-14) break;
         if (jj == nInputs)
         {
            found = 1;
            printf("Bobyqa: simulation results reuse.\n");
            break;
         }
      }

      funcID = odata->numFuncEvals_;
      if (found == 0)
      {
         if (odata->optFunction_ != NULL)
            odata->optFunction_(nInputs, XValues, iOne, localY);
         else if (odata->funcIO_ != NULL)
            odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
         else
         {
            printf("BobyqaOptimizer ERROR: no function evaluator.\n");
            exit(1);
         }
         if (odata->outputLevel_ > 4)
         {
            printf("BobyqaOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            for (ii = 0; ii < nOutputs; ii++) 
               printf("    Y %6d = %16.8e\n", ii+1, localY[ii]);
         }
         funcID = odata->numFuncEvals_++;
         if (psNumBOVars_ + psNumBLVars_ > 0)
         {
            ddata = 0.0;
            for (ii = 0; ii < psNumBOVars_; ii++)
            {
               kk = psBOVars_[ii];
               ddata += psBOWghts_[ii] * localY[kk];
            }
            (*YValue) = ddata;
            ddata = 0.0;
            for (ii = 0; ii < psNumBLVars_; ii++)
            {
               kk = psBLVars_[ii];
               ddata += pow(psBLVals_[ii] - localY[kk], 2.0) * psBLWghts_[ii];
            }
            (*YValue) += ddata;
         }
         else (*YValue) = localY[outputID];
         if (odata->outputLevel_ > 4)
            printf("    ==> Objective function value = %16.8e\n", (*YValue));
      }
      else
      {
         localY[outputID] = psBobyqaSaveY_[ii];
         (*YValue) = localY[outputID];
      }

      if ((psBobyqaSaveHistory_ == 1 && found == 0 &&
          (psBobyqaNSaved_+1)*nInputs < psBobyqaMaxSaved_*10) && 
          psBobyqaNSaved_ < psBobyqaMaxSaved_)
      {
         for (jj = 0; jj < nInputs; jj++)
            psBobyqaSaveX_[psBobyqaNSaved_*nInputs+jj] = XValues[jj];
         psBobyqaSaveY_[psBobyqaNSaved_] = (*YValue);
         psBobyqaNSaved_++;
         infile = fopen("psuade_bobyqa_history","w");
         if (infile != NULL)
         {
            for (ii = 0; ii < psBobyqaNSaved_; ii++)
            {
               fprintf(infile, "999 %d ", nInputs);
               for (kk = 0; kk < nInputs; kk++)
                  fprintf(infile, "%24.16e ", psBobyqaSaveX_[ii*nInputs+kk]);
               fprintf(infile, "%24.16e\n", psBobyqaSaveY_[ii]);
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
            printf("BobyqaOptimizer %6d : \n", odata->numFuncEvals_);
            for (ii = 0; ii < nInputs; ii++)
               printf("    X %6d = %16.8e\n", ii+1, XValues[ii]);
            if (found == 0)
            {
               for (ii = 0; ii < nOutputs; ii++) 
                  printf("    Y %6d = %16.8e\n", ii+1, localY[ii]);
            }
            printf("  *** Current optimal Objective function value = %16.8e\n", 
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
BobyqaOptimizer::BobyqaOptimizer()
{
  if (isScreenDumpModeOn())
  {
    printAsterisks(PL_INFO, 0);
    printf("*   BOBYQA Optimizer Usage Information\n");
    printEquals(PL_INFO, 0);
    printf("* - To run this optimizer in batch mode, first make sure\n");
    printf("*   opt_driver (in your PSUADE input file) has been set\n");
    printf("*   to point to your objective function evaluator.\n");
    printf("* - Set optimization tolerance in psuade.in file\n");
    printf("* - Set maximum number of iterations in psuade.in file\n");
    printf("* - Set num_local_minima to perform multistart optimization\n");
    printf("* - Set optimization print_level to give more screen outputs\n");
    printf("* - If opt_expert mode is turned on, the optimization history\n");
    printf("*   log will be turned on automatically. Previous history file\n");
    printf("*   (psuade_bobyqa_history) will also be reused.\n");
    printf("* - If your opt_driver is a response surface which has more\n");
    printf("*   inputs than the number of optimization inputs, you can fix\n");
    printf("*   some driver inputs by creating an rs_index_file and point\n");
    printf("*   to it in the ANALYSIS section (see user manual).\n");
    printf("* - In opt_expert mode, you can specialize the objective \n");
    printf("*   function by creating a file called psuade_bobyqa_special\n");
    printf("*   in your work directory. This will allow you to create\n");
    printf("*   your own objective function is in the following form: \n");
    printf("\n");
    printf("*        sum_{i=1}^m w_i O_i + sum_{j=1}^n (O_j - C_j)^2\n");
    printf("*\n");
    printf("*   where\n");
    printf("*   m   - number of outputs to be used to form linear sum.\n");
    printf("*   n   - number of outputs to form the squared term.\n");
    printf("*   w_i - weight of output i.\n");
    printf("*   C_j - constraint for output j.\n\n");
    printf("* psuade_bobyqa_special should have the following format: \n");
    printf("*\n");
    printf("\tPSUADE_BEGIN\n");
    printf("\t<m>         /* m in the above formula */\n");
    printf("\t1  <value>  /* the value of w_1 */\n");
    printf("\t3  <value>  /* the value of w_3 */\n");
    printf("\t...\n");
    printf("\t<n>         /* n in the above formula */\n");
    printf("\t2  <value>  /* the value of C_2 */\n");
    printf("\t...\n");
    printf("\tPSUADE_END\n");
    printEquals(PL_INFO, 0);
    printf("To reuse the simulation results (e.g. restart due to abrupt\n");
    printf("termination), turn on save_history and use_history optimization\n");
    printf("options in the ANALYSIS section. You will see a file created\n");
    printf("called 'psuade_bobyqa_history' afterward.\n");
    printAsterisks(PL_INFO, 0);
  }
  objFunction_ = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BobyqaOptimizer::~BobyqaOptimizer()
{
}

// ************************************************************************
// optimize (this function should be used in the library mode)
// ------------------------------------------------------------------------
void BobyqaOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                               double *ubds, int maxfun, double tol)
{
   double *optimalX;
   if (nInputs <= 0)
   {
      printf("BobyqaOptimizer ERROR: nInputs <= 0.\n");
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
   for (int ii = 0; ii < nInputs; ii++) 
      XValues[ii] = optimalX_[ii] = odata->optimalX_[ii];
   delete [] optimalX;
   odata->optimalX_ = NULL;
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void BobyqaOptimizer::optimize(oData *odata)
{
   int    nInputs, printLevel=0, ii, kk, maxfun, nPts=0;
   int    nOutputs, outputID, cmpFlag, bobyqaFlag=7777;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp, *workArray;
   char   filename[500], cinput[500], *cString;
   FILE   *infile=NULL;
   string   iLine;
   ifstream iFile;

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
      printf("BobyqaOptimizer WARNING: tolerance too large.\n");
      printf("                         tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }
   nOutputs = odata->nOutputs_;
   outputID = odata->outputID_;

   if (psOptExpertMode_ == 1)
   {
      strcpy(filename, "psuade_bobyqa_special");
      iFile.open(filename);
      if (iFile.is_open())
      {
         printf("INFO: *** psuade_bobyqa_special file found ***.\n");
         getline (iFile, iLine);
         cmpFlag = iLine.compare("PSUADE_BEGIN");
         if (cmpFlag != 0)
         {
            printf("Bobyqa ERROR: PSUADE_BEGIN not found.\n");
            iFile.close();
            return;
         }
         iFile >> psNumBOVars_;
         if (psNumBOVars_ < 0)
         {
            printf("Bobyqa ERROR: numOVars <= 0\n");
            iFile.close();
            return;
         }
         if (psNumBOVars_ > nOutputs)
         {
            printf("Bobyqa ERROR: numOVars > %d\n", nOutputs);
            iFile.close();
            return;
         }
         psBOVars_  = new int[psNumBOVars_];
         psBOWghts_ = new double[psNumBOVars_];
         for (ii = 0; ii < psNumBOVars_; ii++)
         {
            iFile >> psBOVars_[ii];
            if (psBOVars_[ii] <= 0 || psBOVars_[ii] > nOutputs)
            {
               printf("Bobyqa ERROR: invalid variable index %d.\n",
                      psBOVars_[ii]);
               iFile.close();
               return;
            }
            psBOVars_[ii]--;
            iFile >> psBOWghts_[ii];
         }
         printf("%d linear-term outputs selected: \n", psNumBOVars_);
         for (ii = 0; ii < psNumBOVars_; ii++)
            printf("%4d   weight = %16.8e\n",psBOVars_[ii]+1,psBOWghts_[ii]);
         iFile >> psNumBLVars_;
         if (psNumBLVars_ < 0)
         {
            printf("Bobyqa ERROR: numLVars < 0\n");
            iFile.close();
            return;
         }
         if (psNumBLVars_ > nOutputs)
         {
            printf("Bobyqa ERROR: numLVars > %d\n", nOutputs);
            iFile.close();
            return;
         }
         if (psNumBLVars_ > 0)
         {
            psBLVars_  = new int[psNumBLVars_];
            psBLWghts_ = new double[psNumBLVars_];
            psBLVals_  = new double[psNumBLVars_];
         }
         for (ii = 0; ii < psNumBLVars_; ii++)
         {
            iFile >> psBLVars_[ii];
            if (psBLVars_[ii] <= 0 || psBLVars_[ii] > nOutputs)
            {
               printf("Bobyqa ERROR: invalid variable index %d.\n",
                      psBLVars_[ii]);
               iFile.close();
               return;
            }
            psBLVars_[ii]--;
            iFile >> psBLVals_[ii];
            iFile >> psBLWghts_[ii];
         }
         printf("%d squared-term outputs selected: \n", psNumBLVars_);
         for (ii = 0; ii < psNumBLVars_; ii++)
            printf("%4d   value = %16.8e, weight = %16.8e\n", 
                   psBLVars_[ii]+1, psBLVals_[ii], psBLWghts_[ii]);
         getline (iFile, iLine);
         cmpFlag = iLine.compare("PSUADE_END");
         if (cmpFlag != 0)
         {
            getline (iFile, iLine);
            cmpFlag = iLine.compare("PSUADE_END");
            if (cmpFlag != 0)
            {
               printf("Bobyqa ERROR: PSUADE_END not found.\n");
               iFile.close();
               return;
            }
         }
         iFile.close();
      }
   }
   if (isScreenDumpModeOn() && (psNumBOVars_ + psNumBLVars_ == 0))
      printf("Bobyqa optimizer: selected output for optimization = %d\n", 
             outputID+1);

   maxfun = odata->maxFEval_;
   if ((odata->setOptDriver_ & 1))
   {
      printf("Bobyla: setting optimization simulation driver.\n");
      psBCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   psBobyqaObj_= (void *) odata;
   if (isScreenDumpModeOn())
   {
      printAsterisks(PL_INFO, 0);
      printf("Bobyqa optimizer: max fevals = %d\n", odata->maxFEval_);
      printf("Bobyqa optimizer: tolerance  = %e\n", odata->tolerance_);
      if (printLevel > 1)
         printf("Bobyqa optimizer: rho1, rho2 = %e %e\n", rhobeg, rhoend);
      printEquals(PL_INFO, 0);
   }
#ifdef HAVE_BOBYQA
   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("opt_save_history");
      if (cString != NULL) psBobyqaSaveHistory_ = 1;
      cString = psConfig_->getParameter("opt_use_history");
      if (cString != NULL)
      {
         printf("Bobyqa: use history has been turned on.\n");
         infile = fopen("psuade_bobyqa_history","r");
         if (infile != NULL)
         {
            psBobyqaNSaved_ = 0;
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
                     fscanf(infile, 
                        "%lg",&psBobyqaSaveX_[psBobyqaNSaved_*nInputs+ii]);
                  fscanf(infile, "%lg",&psBobyqaSaveY_[psBobyqaNSaved_]);
                  psBobyqaNSaved_++;
               }
               if (((psBobyqaNSaved_+1)*nInputs > psBobyqaMaxSaved_*10) ||
                   psBobyqaNSaved_ > psBobyqaMaxSaved_) break;
            } 
            fclose(infile);
         }
      }
   }

   nPts = (nInputs + 1) * (nInputs + 2) / 2;
   workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];
   if (isScreenDumpModeOn())
      for (ii = 0; ii < nInputs; ii++) 
         printf("Bobyqa initial X %3d = %e\n", ii+1, XValues[ii]);
   bobyqa_(&nInputs, &nPts, XValues, odata->lowerBounds_,
           odata->upperBounds_, &rhobeg, &rhoend, &bobyqaFlag, &maxfun, 
           workArray);
   if (isScreenDumpModeOn())
      printf("Bobyqa optimizer: total number of evaluations = %d\n",
              odata->numFuncEvals_);
   if (optimalX_ != NULL) delete [] optimalX_;
   optimalX_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) optimalX_[ii] = odata->optimalX_[ii];
   optimalY_ = odata->optimalY_;
   numEvals_ = odata->numFuncEvals_;

   if (psBobyqaSaveHistory_ == 1 && psBobyqaNSaved_ > 0)
   {
      infile = fopen("psuade_bobyqa_history","w");
      if (infile != NULL)
      {
         for (ii = 0; ii < psBobyqaNSaved_; ii++)
         {
            fprintf(infile, "999 %d ", nInputs);
            for (kk = 0; kk < nInputs; kk++)
               fprintf(infile, "%24.16e ", psBobyqaSaveX_[ii*nInputs+kk]);
            fprintf(infile, "%24.16e\n", psBobyqaSaveY_[ii]);
         }
         fclose(infile);
      }
      printf("Bobyqa: history saved in psuade_bobyqa_history\n");
   }
#else
   printf("ERROR : Bobyqa optimizer not installed.\n");
   exit(1);
#endif

   if ((odata->setOptDriver_ & 2) && psBCurrDriver_ >= 0)
   {
      printf("Bobyla INFO: reverting to original simulation driver.\n");
      odata->funcIO_->setDriver(psBCurrDriver_);
   }
   delete [] XValues;
   delete [] workArray;
   if (psBOVars_  != NULL) delete [] psBOVars_;
   if (psBOWghts_ != NULL) delete [] psBOWghts_;
   if (psBLVars_  != NULL) delete [] psBLVars_;
   if (psBLWghts_ != NULL) delete [] psBLWghts_;
   if (psBLVals_  != NULL) delete [] psBLVals_;
   psBLVals_ = NULL;
   psBLWghts_ = NULL;
   psBOWghts_ = NULL;
   psBOVars_ = NULL;
   psBLVars_ = NULL;
   psNumBOVars_ = psNumBLVars_ = 0;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
BobyqaOptimizer& BobyqaOptimizer::operator=(const BobyqaOptimizer &)
{
   printf("BobyqaOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

