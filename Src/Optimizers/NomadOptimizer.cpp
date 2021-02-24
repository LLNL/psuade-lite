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
// Functions for the class NomadOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "NomadOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

#ifdef HAVE_NOMAD
#include "../../External/NOMAD/src/nomad.hpp"
#include "../../External/NOMAD/src/Display.hpp"
using namespace NOMAD;
#endif

oData *psNomadoData_=NULL;
int    psNomadnOutputs_=-1;
double psNomadOptY_=1e50;
double *psNomadOptX_=NULL;
double psNomadTolerance_=1.0e-5;
double *psNomadObjFcnStore_=NULL;
int    psNomadPrintLevel_=-1;
int    psNomadStop_=0;
int    psNomadCurrDriver_=-1;
#define psNomadMaxSaved_ 10000
int     psNomadSaveHistory_=0;
int     psNomadNSaved_=0;
double  psNomadSaveX_[psNomadMaxSaved_*10];
double  psNomadSaveY_[psNomadMaxSaved_*10];
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
NomadOptimizer::NomadOptimizer()
{
   if (isScreenDumpModeOn())
   {
      printAsterisks(PL_INFO, 0);
      printf("*   NOMAD Optimizer Usage Information\n");
      printEquals(PL_INFO, 0);
      printf("* - To run this optimizer, first make sure opt_driver has\n");
      printf("*   been initialized to point to your optimization objective\n");
      printf("*   function evaluator\n");
      printf("* - Set maximum number of iterations in PSUADE input file\n");
      printf("* - Set num_local_minima to perform multi-start optimization\n");
      printf("* - Set optimization print_level to give additonal outputs\n");
// printf("* - In Opt EXPERT mode, the optimization history log will be\n");
// printf("*   turned on automatically. Previous psuade_nomad_history\n");
// printf("*   file will also be reused.\n");
      printAsterisks(PL_INFO, 0);
   }
   nInputs_ = 0;
   inputTypes_ = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
NomadOptimizer::~NomadOptimizer()
{
   if (inputTypes_ != NULL) delete [] inputTypes_;
}

// ************************************************************************
// set discrete variables 
// ------------------------------------------------------------------------
void NomadOptimizer::setDiscreteVariable(int index)
{
   char winput[1000];
   if (index <= 0)
   {
      printf("NOMAD setDiscreteVariable ERROR: variable number <= 0.\n");
      exit(1);
   }
   if (psConfig_ == NULL) psConfig_ = new PsuadeConfig();
   sprintf(winput,"iDiscrete%d", index);
   psConfig_->putParameter(winput);
}

// ************************************************************************
// optimize (this function should be used in the library mode)
// ------------------------------------------------------------------------
void NomadOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                          double *ubds, int nOuts, int maxfun, double tol)
{
   double *optimalX;
   if (nInputs <= 0)
   {
      printf("NomadOptimizer ERROR: nInputs <= 0.\n");
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
   odata->nOutputs_ = nOuts;
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
void NomadOptimizer::optimize(oData *odata)
{
#ifdef HAVE_NOMAD
   int    ii, kk, nInps, nOuts, printLevel, maxfun, idata, option;
   double *XValues, *lbounds, *ubounds, ddata, tol;
   char   pString[1000], *cString;
   FILE   *fp;

   printLevel = odata->outputLevel_;
   nInps  = odata->nInputs_;
   nOuts  = odata->nOutputs_;
   maxfun = odata->maxFEval_;
   psNomadTolerance_ = odata->tolerance_;
   lbounds = odata->lowerBounds_;
   ubounds = odata->upperBounds_;
   XValues = new double[nInps];
   for (ii = 0; ii < nInps; ii++) XValues[ii] = odata->initialX_[ii];
   for (ii = 0; ii < nInps; ii++) odata->optimalX_[ii] = XValues[ii];
   odata->optimalY_ = 1.0e50;
   if ((odata->setOptDriver_ & 1))
   {
      printf("NOMAD: setting optimization simulation driver.\n");
      psNomadCurrDriver_ = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }

   if (psConfig_ != NULL)
   {
      cString = psConfig_->getParameter("opt_save_history");
      if (cString != NULL) psNomadSaveHistory_ = 1;
      cString = psConfig_->getParameter("opt_use_history");
      if (cString != NULL)
      {
         printf("NOMAD: use history has been turned on.\n");
         fp = fopen("psuade_nomad_history","r");
         if (fp != NULL)
         {
            psNomadNSaved_ = 0;
            while (feof(fp) == 0)
            {
               fscanf(fp, "%d %d", &ii, &kk);
               if (ii != 999 || kk != nInps) break;
               else
               {
                  for (ii = 0; ii < nInps; ii++)
                     fscanf(fp, "%lg",
                            &psNomadSaveX_[psNomadNSaved_*nInps+ii]);
                  for (ii = 0; ii < nOuts; ii++)
                     fscanf(fp, "%lg",
                            &psNomadSaveY_[psNomadNSaved_*nOuts+ii]);
                  psNomadNSaved_++;
               }
               if (((psNomadNSaved_+1)*nInps > psNomadMaxSaved_*10) ||
                   ((psNomadNSaved_+1)*nOuts > psNomadMaxSaved_*10))
                  break;
            }
            fclose(fp);
         }
      }
      kk = 0;
      for (ii = 0; ii < nInps; ii++)
      {
         sprintf(pString,"iDiscrete%d", ii+1);
         cString = psConfig_->getParameter(pString);
         if (cString != NULL) kk++;
      }
      if (kk > 0)
      {
         if (inputTypes_ != NULL) delete [] inputTypes_;
         inputTypes_ = new int[nInps];
         for (ii = 0; ii < nInps; ii++) 
         {
            sprintf(pString,"iDiscrete%d", ii+1);
            cString = psConfig_->getParameter(pString);
            if (cString != NULL)
            {
               inputTypes_[ii] = 2;
               if (printLevel > 0)
                  printf("NOMAD input %4d is discrete\n",ii+1);
            }
            else inputTypes_[ii] = 1;
         }
         option = 2;
      }
   }

   psNomadPrintLevel_ = printLevel;
   Display dout ( std::cout );
   dout.precision( DISPLAY_PRECISION_STD );
   Parameters nomadp(dout);
   nomadp.set_DISPLAY_DEGREE(0);

   nomadp.set_DIMENSION(nInps);
   if (isInteractiveModeOn() && nInputs_ == 0 && inputTypes_ == NULL)
   {
      printf("NOMAD can solve either \n");
      printf("1. continuous (which is slow compared to bobyqa) or \n");
      printf("2. mixed-integer optimization.\n");
      sprintf(pString, "Please select (1) or (2) : ");
      option = getInt(1, 2, pString);

      if (option == 2)
      {
         inputTypes_ = new int[nInps];
         printf("NOMAD can handle 3 input types: \n");
         printf("1. Real (or R) - real/continuous\n");
         printf("2. Int  (or I) - integer\n");
         printf("3. Bin  (or B) - binary\n");
         //printf("4. Cat  (or C) - categorical\n");
         for (ii = 0; ii < nInps; ii++)
         {
            sprintf(pString, "Please enter type for input %d : ",ii+1);
            inputTypes_[ii] = getInt(1, 3, pString);
            idata = (int) lbounds[ii];
            if (inputTypes_[ii] == 3)
            {
               if (lbounds[ii] != 0 || lbounds[ii] != 1)
               {
                  printf("ERROR: input type has been set to binary but ");
                  printf("the bounds are not [0,1].\n");
                  printf("INFO: bounds reset to [0,1]\n"); 
                  lbounds[ii] = 0;
                  ubounds[ii] = 1;
               }
            }
            if (inputTypes_[ii] == 2)
            {
               idata = (int) lbounds[ii];
               kk    = (int) ubounds[ii];
               if ((lbounds[ii]-1.0*idata) != 0 || (ubounds[ii]-1.0*kk) != 0)
               {
                  printf("ERROR: input type has been set to integer but ");
                  printf("the bounds are not integers..\n");
                  printf("       Please set the bounds correctly.\n");
                  exit(1);
               }
            }
         }
      }
      if (nOuts > 1)
      {
         printf("The number of outputs = %d\n", nOuts);
         printf("NOMAD will treat the first one as objective function\n");
         printf("and the rest as inequality constraints.\n");
         printf("NOTE: For NOMAD to accept a design point, all inequality\n");
         printf("      constraints have to be <= 0.\n");
      }
   }
   nInputs_ = nInps;

   Point lbnds(nInps), ubnds(nInps), X0(nInps);
   for (ii = 0; ii < nInps; ii++)
   {
      lbnds[ii] = lbounds[ii];
      ubnds[ii] = ubounds[ii];
      if (inputTypes_ != NULL)
      {
         switch(inputTypes_[ii])  
         {
            case 1: nomadp.set_BB_INPUT_TYPE (ii, CONTINUOUS); 
                    X0[ii] = XValues[ii];
                    break;
            case 2: nomadp.set_BB_INPUT_TYPE (ii, INTEGER); 
                    X0[ii] = (int) (XValues[ii] + 0.5);
                    break;
            case 3: nomadp.set_BB_INPUT_TYPE (ii, BINARY); 
                    idata = (int) XValues[ii];
                    if (idata == 0) X0[ii] = 0;
                    else            X0[ii] = 1;
                    break;
            case 4: nomadp.set_BB_INPUT_TYPE (ii, CATEGORICAL); 
                    X0[ii] = (int) XValues[ii];
                    break;
         }
      }
      else X0[ii] = XValues[ii];
   }
   nomadp.set_LOWER_BOUND (lbnds);
   nomadp.set_UPPER_BOUND (ubnds);
   nomadp.set_X0(X0);

   psNomadnOutputs_ = nOuts;
   vector<bb_output_type> bbot (nOuts);
   bbot[0] = OBJ;
   for (ii = 1; ii < nOuts; ii++) bbot[ii] = PB;
   nomadp.set_BB_OUTPUT_TYPE ( bbot );

   nomadp.set_SPECULATIVE_SEARCH ( true );
   if (isInteractiveModeOn() && psOptExpertMode_ == 1)
   {
      sprintf(pString, 
              "Enter value for mesh update basis (default = 8) : ");
      ddata = getDouble(pString);
   } 
   else ddata = 8.0;
   nomadp.set_MESH_UPDATE_BASIS ( ddata );

   if (isInteractiveModeOn() && psOptExpertMode_ == 1)
   {
      sprintf(pString, 
              "Enter value for mesh coarsening exponent (default = 1) : ");
      idata = getInt(0, 10, pString);
   } 
   else idata = 1;
   nomadp.set_MESH_COARSENING_EXPONENT ( idata );
   nomadp.set_MAX_BB_EVAL (maxfun);
   for (ii = 0; ii < nInps; ii++)
   {
      ddata = 0.5 * (ubounds[ii] - lbounds[ii]);
      nomadp.set_INITIAL_MESH_SIZE ( ii, ddata, false );
   }
   //nomadp.set_LH_SEARCH (3, 3);

   nomadp.check();

   PsuadeNomadEvaluator madEvaluator(nomadp);   
   Mads mads(nomadp, &madEvaluator);
   psNomadoData_ = odata;
   psNomadoData_->numFuncEvals_ = 0;
   if (psNomadObjFcnStore_ != NULL) delete [] psNomadObjFcnStore_;
   psNomadObjFcnStore_ = new double[maxfun];
   for (ii = 0; ii < maxfun; ii++) psNomadObjFcnStore_[ii] = -1e35;

   psNomadStop_ = 0;
   mads.run();

   if (optimalX_ != NULL) delete [] optimalX_;
   optimalX_ = new double[nInps];
   for (ii = 0; ii < nInps; ii++) 
      optimalX_[ii] = psNomadoData_->optimalX_[ii];
   optimalY_ = psNomadoData_->optimalY_;
   numEvals_ = psNomadoData_->numFuncEvals_;

   char name[1000];
   if (isScreenDumpModeOn() && psOptExpertMode_ == 1)
   {
      fp = fopen("matlabnomad.m", "w");
      if (fp != NULL)
      {
         fprintf(fp, "A = [\n");
         for (ii = 0; ii < maxfun; ii++) 
         {
            if (psNomadObjFcnStore_[ii] <= -1e35) break;
            else fprintf(fp, "%e\n", psNomadObjFcnStore_[ii]);
         }
         fprintf(fp, "];\n");
         fprintf(fp, "plot(A,'linewidth',2)\n");
         fwritePlotAxes(fp);
         strcpy(name, "Iteration count");
         fwritePlotXLabel(fp, name);
         strcpy(name, "Objective function value");
         fwritePlotYLabel(fp, name);
         strcpy(name, "Nomad Convergence Plot (Some values may be infeasible)");
         fwritePlotTitle(fp, name);
         fclose(fp);
         printf("matlabnomad.m is now available.\n");
      }
   }
    
   if ((odata->setOptDriver_ & 2) && psNomadCurrDriver_ >= 0)
   {
      printf("NOMAD INFO: reverting to original simulation driver.\n");
      odata->funcIO_->setDriver(psNomadCurrDriver_);
   }
   delete [] XValues;
   delete [] psNomadObjFcnStore_;
   psNomadObjFcnStore_ = NULL;
#else
   printf("ERROR : NOMAD optimizer not installed.\n");
   exit(1);
#endif
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
NomadOptimizer& NomadOptimizer::operator=(const NomadOptimizer &)
{
   printf("NomadOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

#ifdef HAVE_NOMAD
// ************************************************************************
// PsuadeNomadEvaluator constructor
// ------------------------------------------------------------------------
PsuadeNomadEvaluator::PsuadeNomadEvaluator(const Parameters &pset): 
                      Evaluator(pset)
{
}

// ************************************************************************
// PsuadeNomadEvaluator destructor
// ------------------------------------------------------------------------
PsuadeNomadEvaluator::~PsuadeNomadEvaluator()
{
}

// ************************************************************************
// PsuadeNomadEvaluator evaluation 
// ------------------------------------------------------------------------
bool PsuadeNomadEvaluator::eval_x(Eval_Point &X,const Double &h_max,
                                  bool &count_eval)
{
   int    ii, jj, funcID, nInps = X.get_n(), goodFlag;
   double *XVals = new double[nInps], *YVals = new double[psNomadnOutputs_];
   bool   retFlag=true;
   FILE   *fp=NULL;

   if (psNomadStop_ == 1)
   {
      YVals[0] = -1e50;
      X.set_bb_output(0, YVals[0]);
      YVals[0] = 0.0;
      for (ii = 1; ii < psNomadnOutputs_; ii++) 
         X.set_bb_output(ii, YVals[0]);
      retFlag = false;
      return retFlag;
   }

   for (ii = 0; ii < nInps; ii++) XVals[ii] = X[ii].value();
   funcID = psNomadoData_->numFuncEvals_++;

   int found = 0;
   for (ii = 0; ii < psNomadNSaved_; ii++)
   {
      for (jj = 0; jj < nInps; jj++)
         if (PABS(psNomadSaveX_[ii*nInps+jj]-XVals[jj])>1.0e-14) break;
      if (jj == nInps)
      {
         found = 1;
         printf("NOMAD: simulation results reuse.\n");
         break;
      }
   }
   if (found == 0)
   {
      if (psNomadoData_->optFunction_ != NULL)
         psNomadoData_->optFunction_(nInps, XVals, psNomadnOutputs_,YVals);
      else if (psNomadoData_->funcIO_ != NULL)
         psNomadoData_->funcIO_->evaluate(funcID,nInps,XVals,
                                          psNomadnOutputs_,YVals,0);
      else
      {
         printf("NomadOptimizer ERROR: no function evaluator.\n");
         exit(1);
      }
      if (psNomadSaveHistory_ == 1 && found == 0 &&
          (psNomadNSaved_+1)*nInps < psNomadMaxSaved_*10 &&
          (psNomadNSaved_+1)*psNomadnOutputs_ < psNomadMaxSaved_*10) 
      {
         for (jj = 0; jj < nInps; jj++)
            psNomadSaveX_[psNomadNSaved_*nInps+jj] = XVals[jj];
         for (jj = 0; jj < psNomadnOutputs_; jj++)
            psNomadSaveY_[psNomadNSaved_*psNomadnOutputs_+jj] = YVals[jj];
         psNomadNSaved_++;
         fp = fopen("psuade_nomad_history","w");
         if (fp != NULL)
         {
            for (ii = 0; ii < psNomadNSaved_; ii++)
            {
               fprintf(fp, "999 %d ", nInps);
               for (jj = 0; jj < nInps; jj++)
                  fprintf(fp, "%24.16e ", psNomadSaveX_[ii*nInps+jj]);
               fprintf(fp, "%24.16e\n", psNomadSaveY_[ii]);
            }
            fclose(fp);
         }
      }
   }
   else
   {
      for (jj = 0; jj < psNomadnOutputs_; jj++)
         YVals[jj] = psNomadSaveY_[ii*psNomadnOutputs_+jj];
   }

   for (ii = 0; ii < psNomadnOutputs_; ii++) X.set_bb_output(ii, YVals[ii]);

   if (YVals[0] < psNomadoData_->optimalY_) 
   {
      goodFlag = 1;
      for (ii = 1; ii < psNomadnOutputs_; ii++) 
         if (YVals[ii] > 0) goodFlag = 0;
      if (goodFlag == 1)
      {
         for (ii = 0; ii < nInps; ii++) 
            psNomadoData_->optimalX_[ii] = XVals[ii];
         psNomadoData_->optimalY_ = YVals[0];
         if (psNomadPrintLevel_ > 0)
            printf("   NOMAD current best Y = %e (nfevals=%d)\n", 
                   YVals[0], funcID+1);
      }
      else
      {
         if (psNomadPrintLevel_ > 1)
         {
            printf("   NOMAD: ignore infeasible solution (nfeval=%d)\n",
                   funcID+1);
            //for (ii = 0; ii < nInps; ii++) 
            //   printf("   X[%5d] = %e\n", ii+1, XVals[ii]);
            //for (ii = 0; ii < psNomadnOutputs_; ii++) 
            //   printf("   Y[%5d] = %e\n", ii+1, YVals[ii]);
         }
      }
   }

   psNomadObjFcnStore_[funcID] = YVals[0];
   int    converged;
   double diff;
   if (funcID >= 30)
   {
      converged = 1;
      for (ii = funcID; ii > funcID-30; ii--) 
      {
         diff = psNomadObjFcnStore_[ii] - psNomadObjFcnStore_[ii-1]; 
         if (diff < 0) diff = -diff;
         if (diff > psNomadTolerance_) 
         {
            converged = 0;
            break;
         }
      }
      if (converged == 1)
      {
         printf("INFO: convergence detected at iteration %d ==> stop.\n",
                funcID);
         psNomadStop_ = 1;
      }
   }

   delete [] XVals;
   delete [] YVals;
   fp = fopen("psuade_stop", "r");
   if (fp != NULL)
   {
      printf("**************** ==> \n");
      printf("INFO: PSUADE found psuade_stop (REQUEST TO TERMINATE).\n");
      printf("      Remove this file if this is not what you desire.\n");
      printf("<== **************** \n");
      fclose(fp);
      psNomadStop_ = 1;
   }
   return retFlag;
}
#endif

