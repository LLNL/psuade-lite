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
// Definition for the class MMOptimizer
// AUTHOR : David Echeverria Ciaurri - Charles Tong
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "MMOptimizer.h"
#include "pData.h"
#include "Psuade.h"
#include "PsuadeUtil.h"

#ifdef HAVE_BOBYQA
extern "C" void bobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);
#endif

#define PSmin(X, Y)  ((X) < (Y) ? (X) : (Y))
#define PSmax(X, Y)  ((X) < (Y) ? (Y) : (X))

void    *psMMBobyqaObj_=NULL;
#define psMMMaxSaved_ 1000
int     psMMNSaved_=0;
double  psMMSaveX_[psMMMaxSaved_*10];
double  psMMSaveY_[psMMMaxSaved_];

#define psMMcMaxSaved_ 10000
int     psMMcNSaved_=0;
double  psMMcSaveX_[psMMcMaxSaved_*10];
double  psMMcSaveY_[psMMcMaxSaved_];

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// local variables
// ------------------------------------------------------------------------
double *MM_Yspec_=NULL;
double MM_normYspec_=1.0;
double *MM_OptimalY_=NULL;
extern "C"
{
  void dgesvd_(char *, char *, int *, int *, double *, int *, double *, 
               double *, int *, double *, int *, double *, int *, int*);
}

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif
  int evaluatefunctionmm_(int *nInps, double *XValues, double *YValue)
  {
    int    ii, jj, kk, funcID, nOutputs, outputID, nInputs;
    int    mm, found, index;
    double *localY, *aux, *xaux; /* *localX2, *localY2, deltaX;*/
    oData  *odata = (oData *) psMMBobyqaObj_;
    FILE   *infile;

    // ------ get history, if any ------
    nInputs  = *nInps;
    nOutputs = odata->nOutputs_;
    if(nOutputs <= 0)
    {
      printf("(file %s line %d) Invalid nOutputs = %d (should be > 0)",
             __FILE__, __LINE__, nOutputs);
      exit(1);
    }
    if (psOptExpertMode_ != 0 && psMMcNSaved_ == 0)
    {
      infile = fopen("psuade_mmc_data","r");
      if (infile != NULL)
      {
        fscanf(infile, "%d %d %d", &psMMcNSaved_, &kk, &mm);
        if ((psMMcNSaved_*nInputs > psMMcMaxSaved_*10) ||
            (psMMcNSaved_*nOutputs > psMMcMaxSaved_))
        {
          printf("PSUADE MM: history file has too much data.\n");
          printf("           History file not used.\n");
          psMMcNSaved_ = 0;
          fclose(infile);
        }
        else if (kk != nInputs)
        {
          printf("PSUADE MM: history file has invalid nInputs.\n");
          fclose(infile);
        }
        else if (mm != nOutputs)
        {
          printf("PSUADE MM: history file has invalid nOutputs.\n");
          fclose(infile);
        }
        else
        {
          for (ii = 0; ii < psMMcNSaved_; ii++)
          {
            fscanf(infile, "%d", &kk);
            if (kk != ii+1)
            {
              printf("PSUADE MM: index mismatch (actual %d vs expected %d).\n",
                     kk, ii+1);
              printf("INFO: check your psuade_mmc_data file.\n");
              fclose(infile);
              psMMcNSaved_ = 0;
              break;
            }
            for (jj = 0; jj < nInputs; jj++)
              fscanf(infile, "%lg", &psMMcSaveX_[ii*nInputs+jj]);
            for (jj = 0; jj < nOutputs; jj++)
              fscanf(infile, "%lg", &psMMcSaveY_[ii*nOutputs+jj]);
          }
          fclose(infile);
        }
      }
    }

    // ------ fetch data and allocate memory ------
    localY   = new double[nOutputs];
    aux      = new double[nInputs];
    xaux     = new double[nInputs];
    outputID = odata->outputID_;
    for (ii = 0; ii < nInputs; ii++) xaux[ii] = XValues[ii];

    // ------ search information in history ------
    found = 0;
    for (ii = 0; ii < psMMcNSaved_; ii++)
    {
      for (jj = 0; jj < nInputs; jj++)
        if (PABS(psMMcSaveX_[ii*nInputs+jj]-XValues[jj])>1e-14) break;
      if (jj == nInputs)
      {
        found = 1;
        printf("MMOptimizer: simulation results found in archive (fine).\n");
        index = ii;
        break;
      }
    }

    // ------ run simulation (always coarse model)------
    funcID = odata->numFuncEvals_;
    if (found == 0)
    {
      odata->funcIO_->evaluate(funcID,nInputs,xaux,nOutputs,localY,0);
      funcID = odata->numFuncEvals_++;
    }
    else
    {
      for (jj = 0; jj < nOutputs; jj++)
        localY[jj] = psMMcSaveY_[index*nOutputs+jj];
    }

    // ------ compute objective function ------
    (*YValue) = 0.0;
    for (ii = 0; ii < nOutputs; ii++) 
      (*YValue) += (localY[ii]-MM_Yspec_[ii])*(localY[ii]-MM_Yspec_[ii]);
    (*YValue) /= MM_normYspec_ * 100.0;

    // ------ update current optimal settings ------
    // ------ output optimization information ------
    if ((*YValue) < odata->optimalY_)
    {
      odata->optimalY_ = (*YValue);
      for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
      for (ii = 0; ii < nOutputs; ii++) MM_OptimalY_[ii] = localY[ii];
      if (odata->outputLevel_ > 2)
      {
        printDashes(PL_INFO, 0);
        for (ii = 0; ii < nInputs; ii++)
          printf("MMOptimizer %6d : X %3d = %16.8e\n",odata->numFuncEvals_,
                 ii+1, XValues[ii]);
        printf("MMOptimizer %6d : Ymin = %16.8e\n",
               odata->numFuncEvals_, odata->optimalY_);
        printf("MMOptimizer %6d : output 1 = %16.8e\n",
               odata->numFuncEvals_, localY[0]);
        printDashes(PL_INFO, 0);
      }
    }

    if (psOptExpertMode_ != 0 && found == 0)
    {
      if ((psMMcNSaved_+1)*nInputs < psMMcMaxSaved_*10 &&
          (psMMcNSaved_+1)*nOutputs < psMMcMaxSaved_)
      {
        for (jj = 0; jj < nInputs; jj++)
          psMMcSaveX_[psMMcNSaved_*nInputs+jj] = XValues[jj];
        for (jj = 0; jj < nOutputs; jj++)
          psMMcSaveY_[psMMcNSaved_*nOutputs+jj] = localY[jj];
        psMMcNSaved_++;
      }
    }
      
    if (psOptExpertMode_ != 0 && found == 0)
    {
      infile = fopen("psuade_mmc_data","w");
      if (infile != NULL)
      {
         fprintf(infile, "%d %d %d\n", psMMcNSaved_, nInputs, nOutputs);
         for (ii = 0; ii < psMMcNSaved_; ii++)
         {
           fprintf(infile, "%d ", ii+1);
           for (jj = 0; jj < nInputs; jj++)
             fprintf(infile, "%24.16e ", psMMcSaveX_[ii*nInputs+jj]);
           for (jj = 0; jj < nOutputs; jj++)
             fprintf(infile, "%24.16e\n", psMMcSaveY_[ii*nOutputs+jj]);
        }
        fclose(infile);
      }
    }
    delete [] localY;
    delete [] aux;
    delete [] xaux;
    return 0;
  }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MMOptimizer::MMOptimizer()
{
  adaptive_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MMOptimizer::~MMOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void MMOptimizer::optimize(oData *odata)
{
  int    nInputs, nOutputs, ii, jj, kk, nPts, maxfun, printLevel;
  int    nIter=1, nIterMAX=200, ncevals=0, nfevals=0, status;
  int    ldu, ldvt, lwork, info, minkn, minmn, *tempS, nAux;
  int    auxNInputs, auxNOutputs, *auxSStates, auxNSamples=0, match;
  double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp, *auxSInputs;
  double *zStar, *X0, *X1, *Yf, *Yfs, *Yc, *Ycs, *pS0, *auxSOutputs;
  double Fval, h=1e50, tolX=1e-3, *targetY, normTargetY, lastFval;
  double *Af, *AfT, *Ac, *AcT, *Uf, *UfT, *Vf, *VfT, *tempI, *tempO;
  double *Sf, *Uc, *VcT, *Vc, *Sc, *WORK, TOLpinv = 1e-10, Fh, tolF;
  double initFval, *workArray, optFval, *optX, *optY;
  char   jobu, jobvt, *auxDriverName, lineIn[501], inString[100];
  FILE   *targetFile, *fp;
  PsuadeData *psData=NULL;
  FunctionInterface *cfuncIO=NULL, *tempFuncIO;
  pData  pPtr, pInputs, pOutputs, pStates;

  // ------ prepare object variables ------
  nInputs  = odata->nInputs_;
  nOutputs = odata->nOutputs_;
  if(nOutputs <= 0)
  {
    printf("nOutputs in file %s line %d has nOutputs <= 0 (%d)", 
           __FILE__, __LINE__, nOutputs);
    exit(1);
  }
  for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
  odata->numFuncEvals_ = 0;
  odata->intData_ = 0;
  odata->optimalY_     = 1.0e50;
  tolX = odata->tolerance_;
  tolF = odata->tolerance_;

  // ------ prepare for adaptive mode if chosen ------
  auxNSamples = 0;
  psData = NULL;
  if (adaptive_ == 1)
  {
    printf("MMOptimizer INFO: adaptive on.\n");
    auxDriverName = odata->funcIO_->getAuxOptimizationDriver();
    if (!strcmp(auxDriverName, "NONE"))
    {
      printf("MMOptimizer ERROR: this optimizer needs aux_opt_driver.\n");
      exit(1);
    }
    fp = fopen(auxDriverName, "r");
    if (fp != NULL)
    {
      fgets(lineIn, 500, fp);
      sscanf(lineIn, "%10c", inString);
      if (strncmp(inString, "PSUADE_IO",9))
      {
        printf("MMOptimizer: aux opt driver is not a data file.\n");
        printf("             default back to non-adaptive.\n");
        adaptive_ = 0;
        printf("MMOptimizer INFO: adaptive turned off since aux_driver");
        printf(" is not a data file.\n");
      }
      else
      {
        fscanf(fp, "%d %d\n", &auxNInputs, &auxNOutputs);
        if (nInputs != auxNInputs || nOutputs != auxNOutputs)
        {
          printf("MMOptimizer ERROR: nInputs/nOutputs mismatch with\n");
          printf("                   auxiliary opt data file.\n");
          fclose(fp);
          exit(1);
        }
      }
      fclose(fp);
    }
    auxNSamples = 0;

    if (adaptive_ == 1)
    {
      psData = new PsuadeData();
      status = psData->readPsuadeFile(auxDriverName);
      if (status != 0) 
      {
        printf("MMOptimizer ERROR: cannot read aux opt data file %s.\n",
               auxDriverName);
        exit(1);
      }
      psData->getParameter("input_ninputs", pPtr);
      auxNInputs = pPtr.intData_;
      psData->getParameter("output_noutputs", pPtr);
      auxNOutputs = pPtr.intData_;
      psData->getParameter("method_nsamples", pPtr);
      auxNSamples = pPtr.intData_;
      auxSInputs  = new double[auxNSamples*auxNInputs];
      auxSOutputs = new double[auxNSamples*auxNOutputs];
      auxSStates  = new int[auxNSamples];
      psData->getParameter("output_sample", pOutputs);
      psData->getParameter("input_sample", pInputs);
      psData->getParameter("output_states", pStates);
      for (ii = 0; ii < auxNSamples*auxNInputs; ii++)
        auxSInputs[ii] = pInputs.dbleArray_[ii];
      for (ii = 0; ii < auxNSamples*auxNOutputs; ii++)
      {
        auxSOutputs[ii] = pOutputs.dbleArray_[ii];
        if (auxSOutputs[ii] == PSUADE_UNDEFINED)
        {
          printf("MMOptimizer ERROR: some sample outputs in %s",
                 auxDriverName);
          printf(" are undefined. Prune and re-run.\n");
          exit(1);
        }
      }
      for (ii = 0; ii < auxNSamples; ii++)
        auxSStates[ii] = pStates.intArray_[ii];
      pInputs.clean();
      pOutputs.clean();
      pStates.clean();
      cfuncIO = createFunctionInterfaceGivenAppDriver(auxNInputs,
                                  auxNOutputs, auxDriverName);
      tempFuncIO = odata->funcIO_;
      printf("MMOptimizer INFO: response surface has been created");
      printf(" from aux_driver (nSamples = %d).\n", auxNSamples);
    }
  }

  // read target (output target values) file, if any
  if (strcmp(odata->targetFile_, "NONE"))
  {
    targetFile = fopen(odata->targetFile_, "r");
    if (targetFile == NULL)
    {
      printf("MMOptimizer ERROR: cannot find target file %s.\n",
             odata->targetFile_);
      exit(1);
    }
    targetY = new double[nOutputs];
    for (ii = 0; ii < nOutputs; ii++)
      fscanf(targetFile,"%lg",&(targetY[ii]));
    fclose(targetFile);
    normTargetY = 0.0;
    for (ii = 0; ii < nOutputs; ii++) 
      normTargetY += targetY[ii] * targetY[ii];
    normTargetY = sqrt(normTargetY);
  }
  else
  {
    printf("MMOptimizer WARNING: cannot find target file.\n");
    printf("                     will assume target Fvalues = 0.\n");
    printf("NOTE: if this assumption is wrong, use target file.\n");
    targetY = new double[nOutputs];
    for (ii = 0; ii < nOutputs; ii++) targetY[ii] = 0.0;
    normTargetY = 1.0;
  }
  MM_Yspec_ = new double[nOutputs];
  for (ii = 0; ii < nOutputs; ii++) MM_Yspec_[ii] = targetY[ii];
  MM_normYspec_ = normTargetY;
  MM_OptimalY_ = new double[nOutputs];
  for (ii = 0; ii < nOutputs; ii++) MM_OptimalY_[ii] = 0.0;

  // ------ output target data information
  printAsterisks(PL_INFO, 0);
  printAsterisks(PL_INFO, 0);
  printf(" *         MANIFOLD MAPPING OPTIMIZATION BEGINS\n");
  printAsterisks(PL_INFO, 0);
  if (odata->outputLevel_ > 1 && MM_normYspec_ != PSUADE_UNDEFINED)
  {
    printEquals(PL_INFO, 0);
    printf(" MANIFOLD MAPPING: Specifications: (target output values)\n");
    for (ii = 0; ii < nOutputs; ii++)
      printf("   Output %3d = %20.14f\n", ii+1, MM_Yspec_[ii]);
    printEquals(PL_INFO, 0);
  }

  // ------ allocate temporary memory
  XValues = new double[nInputs+1];
  Yf      = new double[nOutputs];
  Yfs     = new double[nOutputs * (nIterMAX+1)];
  Af      = new double[nOutputs * nInputs];
  Uf      = new double[nOutputs*PSmax(nOutputs,nInputs)];
  AfT     = new double[nOutputs * nInputs];
  Vf      = new double[nInputs * nInputs];
  Sf      = new double[nInputs*PSmax(nOutputs,nInputs)];

  Yc      = new double[nOutputs];
  Ycs     = new double[nOutputs * (nIterMAX+1)];
  Ac      = new double[nOutputs * nInputs];
  Uc      = new double[nOutputs*PSmax(nOutputs,nInputs)];
  AcT     = new double[nOutputs * nInputs];
  Vc      = new double[nInputs * nInputs];
  Sc      = new double[nInputs*PSmax(nOutputs,nInputs)];

  X0      = new double[nInputs+1];
  pS0     = new double[nOutputs*nOutputs];
  UfT     = new double[nOutputs*PSmax(nOutputs,nInputs)];
  VfT     = new double[nInputs * nInputs];
  VcT     = new double[nInputs * nInputs];

  WORK    = new double[5 * PSmax(nOutputs,nInputs)];
  X1      = new double[nInputs + 1];
  zStar   = new double[nInputs];

  // Bobyqa STOP Criteria 
  rhobeg = odata->upperBounds_[0] - odata->lowerBounds_[0];
  for (ii = 1; ii < nInputs; ii++)
  {
    dtemp = odata->upperBounds_[ii] - odata->lowerBounds_[ii];
    if (dtemp > rhobeg) rhobeg = dtemp;
  }
  rhobeg *= 0.05;
  rhoend = rhobeg * odata->tolerance_;
  if (rhobeg < rhoend)
  {
    printf("MMOptimizer WARNING: tolerance too large.\n");
    printf("                     tolerance reset to 1.0e-6.\n");
    rhoend = rhobeg * 1.0e-6;
  }

  // Initializing state variables
  for (ii = 0; ii < nInputs; ii++) X0[ii]  = 0.0;

  for (ii = 0; ii < nOutputs; ii++)
  {
    for (jj = 0; jj < nOutputs; jj++)
    {
      if (ii==jj) pS0[ii + jj*nOutputs] = 1.0;
      else        pS0[ii + jj*nOutputs] = 0.0;
    }
  }

  for (ii = 0; ii < nInputs; ii++) XValues[ii] = odata->initialX_[ii];


  // ------ call optimizer (coarse model) ------
  if (odata->outputLevel_ > 1)
  {
    printEquals(PL_INFO, 0);
    printf(" MANIFOLD MAPPING: Initial coarse model optimization: \n");
    for (ii = 0; ii < nInputs; ii++)
      printf("   Initial X %3d = %20.14f\n", ii+1, XValues[ii]);
    printEquals(PL_INFO, 0);
  }

#ifdef HAVE_BOBYQA
  maxfun = 1000;
  printLevel = 6666;
  if (adaptive_ == 1) odata->funcIO_ = cfuncIO;
  else                odata->funcIO_->setDriver(2);
  psMMBobyqaObj_ = odata;
  nPts = (nInputs + 1) * (nInputs + 2) / 2;
  workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];
  bobyqa_(&nInputs, &nPts, XValues, odata->lowerBounds_, odata->upperBounds_, 
          &rhobeg, &rhoend, &printLevel, &maxfun, workArray);
  delete [] workArray;
  if (adaptive_ == 1) odata->funcIO_ = tempFuncIO;
#else
  printf("ERROR : Bobyqa optimizer not installed.\n");
  exit(1);
#endif
  ncevals = odata->numFuncEvals_;

  // zStar = XValues
  for (ii = 0; ii < nInputs; ii++) X0[ii] = zStar[ii] = XValues[ii];

  if (odata->outputLevel_ > 1)
  {
    printEquals(PL_INFO, 0);
    printf(" MANIFOLD MAPPING: Initial coarse model optimum: \n");
    for (ii = 0; ii < nInputs; ii++)
      printf("  Input %3d = %16.8e\n", ii+1, X0[ii]);
    printf("    Output   1 = %16.8e\n", MM_OptimalY_[0]);
    printf(" Optimal function value = %16.8e\n", odata->optimalY_);
    printEquals(PL_INFO, 0);
  }

  // 1.A.1 EVALUATING THE FINE MODEL AT THE COARSE OPTIMUM 
  odata->funcIO_->setDriver(1);      /* we evaluate the fine model */
  nfevals++;
  kk = fineSolve(odata, nInputs, X0, nOutputs, Yf, nfevals);
  if (kk == 1) nfevals--;

  Fval = 0.0;
  for (ii = 0; ii < nOutputs; ii++)
  {
    Fval += (Yf[ii] - targetY[ii]) * (Yf[ii] - targetY[ii]);
    Yfs[ii+(nIter-1)*nOutputs] = Yf[ii];
  }
  Fval = sqrt(Fval)/normTargetY * 100.00;
  optFval = Fval;
  optX = new double[nInputs];
  for (ii = 0; ii < nInputs; ii++) optX[ii] = X0[ii];
  optY = new double[nOutputs];
  for (ii = 0; ii < nOutputs; ii++) optY[ii] = Yf[ii];
  if (odata->outputLevel_ > 1)
  {
    printEquals(PL_INFO, 0);
    printf("Fine model evaluated initially at: \n");
    for (ii = 0; ii < nInputs; ii++)
      printf("   Input %d  = %16.8e\n", ii+1, X0[ii]);
    for (ii = 0; ii < nOutputs; ii++)
      printf("      Output %d = %16.8e\n", ii+1, Yf[ii]);
    printf("The corresponding cost function value = %16.8e\n", Fval);
    printEquals(PL_INFO, 0);
  }
  initFval = Fval;

  // 1.B EVALUATING THE COARSE MODEL AT THE COARSE OPTIMUM 
  if (adaptive_ == 1)
     cfuncIO->evaluate(odata->numFuncEvals_,nInputs,X0,nOutputs,Yc,0);
  else
  {
    odata->funcIO_->setDriver(2);      /* we evaluate the coarse model */
    odata->funcIO_->evaluate(odata->numFuncEvals_,nInputs,X0,nOutputs,Yc,0);
  }
  for (ii = 0; ii < nOutputs; ii++)
    Ycs[ii + (nIter-1)*nOutputs] = Yc[ii];
  ncevals = odata->numFuncEvals_;

  printEquals(PL_INFO, 0);
  printf("Iteration    # fEvals   # cEval   Fine Cost Func\n");
  printEquals(PL_INFO, 0);
  printf("%5d       %5d      %5d     %13.6e\n",
         nIter,nfevals,ncevals,Fval);
  nAux = 0;

  while (1)
  {
    // 2. COMPUTING THE NEXT SPECIFICATIONS 
    for (ii = 0; ii < nOutputs; ii++)
    {
      MM_Yspec_[ii] = 0.0;
      for (jj = 0; jj < nOutputs; jj++)
        MM_Yspec_[ii] += pS0[ii + jj*nOutputs]*(targetY[jj] - Yf[jj]);
      MM_Yspec_[ii] += Yc[ii];
    }

    MM_normYspec_ = 0.0;
    for (ii = 0; ii < nOutputs; ii++) 
      MM_normYspec_ += MM_Yspec_[ii]*MM_Yspec_[ii];
    MM_normYspec_ = sqrt(MM_normYspec_);
    lastFval = Fval;

    for (ii = 0; ii < nInputs; ii++) X1[ii] = X0[ii];

    if (odata->outputLevel_ > 2)
    {
      printDashes(PL_INFO, 0);
      printf("===> MM: coarse model initial X at iteration %d:\n",nfevals);
      for (ii = 0; ii < nInputs; ii++)
        printf("===>  Input %3d = %16.8e\n", ii+1, X1[ii]);
      printDashes(PL_INFO, 0);
    }

#ifdef HAVE_BOBYQA
    maxfun = 1000;
    if (adaptive_ == 1) odata->funcIO_ = cfuncIO;
    else                odata->funcIO_->setDriver(2);
    nPts = (nInputs + 1) * (nInputs + 2) / 2;
    workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];
    printLevel = 6666;
    bobyqa_(&nInputs, &nPts, X1, odata->lowerBounds_, odata->upperBounds_, 
            &rhobeg, &rhoend, &printLevel, &maxfun, workArray);
    delete [] workArray;
    if (adaptive_ == 1) odata->funcIO_ = tempFuncIO;
#else
    printf("ERROR : Bobyqa optimizer not installed.\n");
    exit(1);
#endif
    ncevals = odata->numFuncEvals_;

    if (odata->outputLevel_ > 2)
    {
      printDashes(PL_INFO, 0);
      printf("===> MM: coarse model optimum at iteration %d:\n", nfevals);
      for (ii = 0; ii < nInputs; ii++)
        printf("===>  Input %3d = %16.8e\n", ii+1, X0[ii]);
      printf("===>    Output   1 = %16.8e\n", MM_OptimalY_[0]);
      printf("===> Optimal function value = %16.8e\n", odata->optimalY_);
      printDashes(PL_INFO, 0);
    }

    // 3.A.1 EVALUATING THE FINE MODEL AT X1
    odata->funcIO_->setDriver(1);  /* we evaluate the fine model */
    nfevals++;
    kk = fineSolve(odata, nInputs, X1, nOutputs, Yf, nfevals);
    if (kk == 1) nfevals--;
    nAux++;
    Fval = 0.0;
    for (ii = 0; ii < nOutputs; ii++)
    {
      Fval += (Yf[ii] - targetY[ii])*(Yf[ii] - targetY[ii]);
      Yfs[ii + nIter*nOutputs] = Yf[ii];
    }
    Fval = sqrt(Fval) / normTargetY * 100.0;
    if (odata->outputLevel_ > 0)
    {
      printEquals(PL_INFO, 0);
      printf("===> Fine model evaluated at coarse optimum (iter = %d): \n",
             nfevals);
      for (ii = 0; ii < nInputs; ii++)
        printf("   Input %d  = %16.8e\n", ii+1, X1[ii]);
      for (ii = 0; ii < nOutputs; ii++)
        printf("      Output %d = %16.8e\n", ii+1, Yf[ii]);
      printf("The corresponding cost function value = %16.8e\n", Fval);
      printEquals(PL_INFO, 0);
    }
    if (Fval < optFval)
    {
      optFval = Fval;
      for (ii = 0; ii < nInputs; ii++) optX[ii] = X0[ii];
      for (ii = 0; ii < nOutputs; ii++) optY[ii] = Yf[ii];
    }

    // 3.A.2 store the new data for future use
    if (adaptive_ == 1 && auxNSamples > 0 && psData != NULL && nAux > 0)
    {
      match = 0;
      nAux = 0;
      //for (ii = 0; ii < auxNSamples; ii++)
      //{
      //   xDist = 0.0;
      //   for (jj = 0; jj < auxNInputs; jj++)
      //   {
      //      dtemp = odata->upperBounds_[jj] - odata->lowerBounds_[jj];
      //      dtemp = (auxSInputs[ii*auxNInputs+jj] - X1[jj]) / dtemp;
      //      xDist += dtemp * dtemp;
      //   }
      //   xDist = sqrt(xDist) / auxNInputs;
      //   if (xDist < 1.0e-4) {match = 1; break;}
      //}
      if (match == 0)
      {
        tempI = auxSInputs;
        tempO = auxSOutputs;
        tempS = auxSStates;
        auxSInputs  = new double[(auxNSamples+1)*auxNInputs];
        auxSOutputs = new double[(auxNSamples+1)*auxNOutputs];
        auxSStates  = new int[(auxNSamples+1)];
        for (ii = 0; ii < auxNSamples*auxNInputs; ii++)
          auxSInputs[ii] = tempI[ii];
        for (ii = 0; ii < auxNSamples*auxNOutputs; ii++)
          auxSOutputs[ii] = tempO[ii];
        for (ii = 0; ii < auxNSamples; ii++) auxSStates[ii] = tempS[ii];

        for (ii = 0; ii < auxNSamples; ii++)
        {
          for (jj = 0; jj < auxNInputs; jj++)
            if (PABS(tempI[ii*auxNInputs+jj]-X1[jj]) > 1.0e-15) break;
          if (jj == auxNInputs) break;
        }
        if (jj != auxNInputs)
        {
          for (ii = 0; ii < auxNInputs; ii++) 
            auxSInputs[auxNSamples*auxNInputs+ii] = X1[ii];
          for (ii = 0; ii < auxNOutputs; ii++) 
            auxSOutputs[auxNSamples*auxNOutputs+ii] = Yf[ii];
          auxSStates[auxNSamples] = 1;
          auxNSamples++;
        }
        psData->updateInputSection(auxNSamples, auxNInputs, NULL, NULL,
                         NULL,auxSInputs, NULL,NULL,NULL,NULL,NULL);
        psData->updateOutputSection(auxNSamples, auxNOutputs, auxSOutputs,
                                    auxSStates, NULL);
        psData->writePsuadeFile(auxDriverName,0);
        delete [] tempI;
        delete [] tempO;
        delete [] tempS;
        delete cfuncIO;
        cfuncIO = createFunctionInterfaceGivenAppDriver(auxNInputs,
                           auxNOutputs, auxDriverName);
      }
    }

    // 3.A.3 STOPPING CRITERIA
    h = 0.0;
    for (ii = 0; ii < nInputs; ii++)
      h += (X1[ii] - X0[ii]) * (X1[ii] - X0[ii]);
    h = sqrt(h);
    Fh = Fval - lastFval;
    if (Fh < 0) Fh = - Fh;
    printEquals(PL_INFO, 0);
    printf("Iteration  # fEvals  # cEval   Fine Cost Func  ||x_{k+1}-x_k||\n");
    printf("%5d     %5d     %5d     %13.6e       %4.2e\n", 
           nIter+1,nfevals,ncevals,Fval,h);
    printEquals(PL_INFO, 0);
    if (nIter == nIterMAX || Fh < tolF) break;

    // 3.B EVALUATING THE COARSE MODEL AT X1
    if (adaptive_ == 1)
      cfuncIO->evaluate(odata->numFuncEvals_,nInputs,X1,nOutputs,Yc,0);
    else
    {
      odata->funcIO_->setDriver(2);      /* we evaluate the coarse model */
      odata->funcIO_->evaluate(odata->numFuncEvals_,nInputs,X1,nOutputs,Yc,0);
    }
    for (ii = 0; ii < nOutputs; ii++) Ycs[ii+nIter*nOutputs] = Yc[ii];
    ncevals = odata->numFuncEvals_;

    if (odata->outputLevel_ > 1)
    {
      printDashes(PL_INFO, 0);
      printf("===> MM: coarse model evaluated at next guess (iter = %d)\n",
             nfevals-1);
      for (ii = 0; ii < nInputs; ii++)
        printf("===>    Input %3d = %16.8e\n", ii+1, X1[ii]);
      for (ii = 0; ii < nOutputs; ii++)
        printf("===>    Output %3d = %16.8e\n", ii+1, Yc[ii]);
      printDashes(PL_INFO, 0);
    }

    // 4. COMPUTING Af AND Ac
    // A. ESTABLISH DIMENSIONS (m x min(nIter,n)))
    minkn = PSmin(   nIter,nInputs);
    minmn = PSmin(nOutputs,  minkn);

    for (ii = 0; ii < minkn; ii++)
    {
      for (jj = 0; jj < nOutputs; jj++)
      {
        Af[jj + ii*nOutputs] = Yf[jj] - Yfs[jj + (nIter-ii-1)*nOutputs];
        Ac[jj + ii*nOutputs] = Yc[jj] - Ycs[jj + (nIter-ii-1)*nOutputs];
      }
    }

    // B. COPYING MATRICES FOR FORTRAN
    for (ii = 0; ii < nOutputs; ii++)
    {
      for(jj = 0; jj < minkn; jj++)
      {
        AfT[ii + nOutputs*jj] = Af[ii + jj*nOutputs];
        AcT[ii + nOutputs*jj] = Ac[ii + jj*nOutputs];
      }
    }

    // C. SVD DECOMPOSITION
    jobu  = 'S';
    jobvt = 'S';
    ldu   = nOutputs;
    ldvt  = minkn;
    lwork = 5 * PSmax(nOutputs,nInputs);
    dgesvd_(&jobu, &jobvt, &nOutputs, &minkn, AfT, &nOutputs, Sf, Uf, 
            &ldu, VfT, &ldvt, WORK, &lwork, &info);

    dgesvd_(&jobu, &jobvt, &nOutputs, &minkn, AcT, &nOutputs, Sc, Uc,
            &ldu, VcT, &ldvt, WORK, &lwork, &info);

    // D. COMPUTING pS0
    for (ii = 0; ii < minmn; ii++)
    {
      if (Sf[ii] < TOLpinv) Sf[ii] = 0.0;
      else                  Sf[ii] = 1/Sf[ii];;
    }

    // pinv(Sf)*U' U (in C):  nOutputs x minmn
    // U'
    for (ii = 0; ii < minmn; ii++)
      for(jj = 0; jj < nOutputs; jj++)
        UfT[ii + jj*minmn] = Uf[jj + ii*nOutputs];

    for (ii = 0; ii < minmn; ii++)
      for (jj = 0; jj < nOutputs; jj++)
        UfT[ii + jj*minmn] = Sf[ii]*UfT[ii + jj*minmn];

    // V
    for (ii = 0; ii < minkn; ii++)
      for(jj = 0; jj < minkn; jj++) Vf[jj + minkn*ii] = VfT[ii + jj*minkn];

    for (ii =0; ii <minkn; ii++)
    {
      for (jj=0; jj<nOutputs; jj++)
      {
        Af[ii + jj*minkn] = 0.0;
        for (kk = 0; kk < minmn; kk++)
          Af[ii + jj*minkn] += Vf[ii + kk*minkn]*UfT[kk + jj*minmn];
      }
    }

    // D2. Ac * pinv(Af) : Ac       (nOutputs x    minkn)
    //                     pinv(Af) (   minkn x nOutputs)
    for (ii = 0; ii < nOutputs; ii++)
    {
      for (jj =0; jj < nOutputs; jj++)
      {
        pS0[ii + jj*nOutputs] = 0.0;
        for (kk = 0; kk < minkn; kk++)
          pS0[ii + jj*nOutputs] += Ac[ii + kk*nOutputs]*Af[kk + jj*minkn];
      }
    }

    // D3. Ac * pinv(Af) + Im
    for (ii = 0; ii < nOutputs; ii++) pS0[ii + ii*nOutputs] += 1.0;

    // D4. Ac * pinv(Af) + Im - Uc*Uc'
    //                               Uc (in C) :  nOutputs x minmn
    for (ii = 0; ii < nOutputs; ii++)
    {
      for (jj = 0; jj < nOutputs; jj++)
      {
        for (kk = 0; kk < minmn; kk++)
          pS0[ii+jj*nOutputs] -= Uc[ii+kk*nOutputs]*Uc[jj+kk*nOutputs];
      }
    }
    nIter = nIter + 1;

    // updating variables 
    for (ii = 0; ii < nInputs; ii++) X0[ii] = X1[ii];
  }
  if (Fval == initFval)
  {
    printf(" ============================================================\n");
    printf(" ============================================================\n");
    printf("\n MANIFOLD MAPPING WARNING: \n");
    printf("     There seems to be no progress from the beginning.\n");
    printf("     This can be due to the initial guess being the optimal\n");
    printf("     or the coarse model is not good enough.\n\n");
    printf(" ============================================================\n");
    printf(" ============================================================\n");
  }

  // store data 
  if (psData != NULL)
  {
    psData->updateInputSection(auxNSamples, auxNInputs, NULL, 
                             NULL,NULL,auxSInputs,NULL,NULL,NULL,NULL,NULL);
    psData->updateOutputSection(auxNSamples, auxNOutputs, auxSOutputs,
                                auxSStates, NULL);
    psData->writePsuadeFile(auxDriverName,0);
  }

  // DIAGNOSTICS 
  if (nIter == nIterMAX)
    printf("\n Maximum number of %4d MM iterations reached.\n\n", nIterMAX);
  else if (Fval < initFval)
  {
    printAsterisks(PL_INFO, 0);
    printAsterisks(PL_INFO, 0);
    printf(" MANIFOLD MAPPING: Manifold-mapping solution:\n");
    for(ii = 0; ii < nInputs; ii++)
      printf("   Input %3d = %16.8e\n", ii+1, optX[ii]);
    for(ii = 0; ii < nOutputs; ii++)
      printf("  Output %3d = %16.8e\n", ii+1, optY[ii]);
    printf("  Function value = %16.8e\n", optFval);
    printf("  Number of fine   model evaluations = %d\n", nfevals);
    printf("  Number of coarse model evaluations = %d\n", ncevals);
    printAsterisks(PL_INFO, 0);
    printAsterisks(PL_INFO, 0);
  }
  printAsterisks(PL_INFO, 0);
  printf(" *       MANIFOLD MAPPING OPTIMIZATION ENDS\n");
  printAsterisks(PL_INFO, 0);
  printAsterisks(PL_INFO, 0);
  odata->optimalY_ = optFval;
  for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = optX[ii];

  if (adaptive_ == 1 && auxNSamples > 0 && psData != NULL)
  {
    delete [] auxSInputs;
    delete [] auxSOutputs;
    delete [] auxSStates;
    delete psData;
    delete cfuncIO;
  }
  delete [] X0;
  delete [] X1;
  delete [] XValues;
  delete [] zStar;
  delete [] Yf;
  delete [] Yc;
  delete [] Yfs;
  delete [] Ycs;
  delete [] pS0;
  delete [] targetY;
  delete [] MM_Yspec_;
  delete [] Af;
  delete [] AfT;
  delete [] Ac;
  delete [] AcT;
  delete [] Uf;
  delete [] UfT;
  delete [] Uc;
  delete [] Vc;
  delete [] Vf;
  delete [] VfT;
  delete [] VcT;
  delete [] Sf;
  delete [] Sc;
  delete [] WORK;
  delete [] optX;
  delete [] optY;
  delete [] MM_OptimalY_;
  MM_Yspec_ = NULL;
  MM_OptimalY_ = NULL;
}

// ************************************************************************
// set parameter
// ------------------------------------------------------------------------
void MMOptimizer::setParam(char *sparam)
{
  if (!strcmp(sparam, "setAdaptive")) adaptive_ = 1;
  return;
}

// ************************************************************************
// fine solve 
// ------------------------------------------------------------------------
int MMOptimizer::fineSolve(oData *odata, int nInputs, double *X0, 
                            int nOutputs, double *Yf, int runNumber)
{
  int  ii, jj, kk, mm, found;
  FILE *infile;

  if (psOptExpertMode_ != 0 && psMMNSaved_ == 0)
  {
    infile = fopen("psuade_mm_data","r");
    if (infile != NULL)
    {
      fscanf(infile, "%d %d %d", &psMMNSaved_, &kk, &mm);
      if ((psMMNSaved_*nInputs > psMMMaxSaved_*10) ||
          (psMMNSaved_*nOutputs > psMMMaxSaved_))
      {
        printf("PSUADE MM: history file has too much data.\n");
        printf("           History file is not used.\n");
        fclose(infile);
        psMMNSaved_ = 0;
      }
      else if (kk != nInputs)
      {
        printf("PSUADE MM: history file has invalid nInputs.\n");
        fclose(infile);
        psMMNSaved_ = 0;
      }
      else if (mm != nOutputs)
      {
        printf("PSUADE MM: history file has invalid nOutputs.\n");
        fclose(infile);
        psMMNSaved_ = 0;
      }
      else
      {
        for (ii = 0; ii < psMMNSaved_; ii++)
        {
          fscanf(infile, "%d", &kk);
          if (kk != ii+1)
          {
            printf("PSUADE MM: index mismatch (actual %d vs expected %d).\n",
                   kk, ii+1);
            printf("INFO: check your psuade_mm_data file.\n");
            psMMNSaved_ = 0;
            break;
          }
          for (jj = 0; jj < nInputs; jj++)
            fscanf(infile, "%lg", &psMMSaveX_[ii*nInputs+jj]);
          for (jj = 0; jj < nOutputs; jj++)
            fscanf(infile, "%lg", &psMMSaveY_[ii*nOutputs+jj]);
        }
        fclose(infile);
      }
    }
  }

  found = 0;
  for (ii = 0; ii < psMMNSaved_; ii++)
  {
    for (jj = 0; jj < nInputs; jj++)
      if (PABS(psMMSaveX_[ii*nInputs+jj]-X0[jj])>1.0e-14) break;
    if (jj == nInputs)
    {
      found = 1;
      printf("MMOptimizer: simulation results found in archive (fine).\n");
      break;
    }
  }

  if (found == 0)
  {
    odata->funcIO_->evaluate(runNumber,nInputs,X0,nOutputs,Yf,0);
    odata->intData_++;
  }
  else
  {
    for (jj = 0; jj < nOutputs; jj++)
      Yf[jj] = psMMSaveY_[ii*nOutputs+jj];
  }

  if (psOptExpertMode_ != 0 && found == 0)
  {
    if (((psMMNSaved_+1) * nInputs < psMMMaxSaved_*10) &&
        ((psMMNSaved_+1) * nOutputs < psMMMaxSaved_))
    {
      for (jj = 0; jj < nInputs; jj++)
        psMMSaveX_[psMMNSaved_*nInputs+jj] = X0[jj];
      for (jj = 0; jj < nOutputs; jj++)
        psMMSaveY_[psMMNSaved_*nOutputs+jj] = Yf[jj];
      psMMNSaved_++;
    }
  }

  if (psOptExpertMode_ != 0 && found == 0)
  {
    infile = fopen("psuade_mm_data","w");
    if (infile != NULL)
    {
      fprintf(infile, "%d %d %d\n", psMMNSaved_, nInputs, nOutputs);
      for (ii = 0; ii < psMMNSaved_; ii++)
      {
        fprintf(infile, "%d ", ii+1);
        for (jj = 0; jj < nInputs; jj++)
          fprintf(infile, "%24.16e ", psMMSaveX_[ii*nInputs+jj]);
        for (jj = 0; jj < nOutputs; jj++)
          fprintf(infile, "%24.16e ", psMMSaveY_[ii*nOutputs+jj]);
        fprintf(infile, "\n");
      }
      fclose(infile);
    }
  }
  return found;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
MMOptimizer& MMOptimizer::operator=(const MMOptimizer &)
{
  printf("MMOptimizer operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

