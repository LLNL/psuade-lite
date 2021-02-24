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
// Functions for the class SMOptimizer
// AUTHOR : David Echeverria Ciaurri - Charles Tong
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "PsuadeUtil.h"
#include "SMOptimizer.h"

#ifdef HAVE_BOBYQA
extern "C" void bobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);
#endif

void   *psSMBobyqaObj_=NULL;
double *px0, *B0, *x0, *fx, *y, normy;
int coarseoptimtype = 1;
//void savesmaux(int, int, double *, double *, double *);

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif
int evaluatefunctionsm_(int *nInps, double *XValues, double *YValue)
{
  int    i, j, funcID, nInputs, nOutputs, outputID, currDriver;
  double *localY, *aux, *xaux;
  oData  *odata;

  // ------ actual data from simulation should be 1 more ------
  odata    = (oData *) psSMBobyqaObj_;
  nInputs  = *nInps;
  nOutputs = odata->nOutputs_;
  if(nOutputs <= 0)
  {
    printf("(file %s line %d) Invalid nOutputs = %d (should be > 0)",
            __FILE__, __LINE__, nOutputs);
    exit(1);
  }

  // allocate memory
  localY   = new double[nOutputs];
  aux      = new double[nInputs];
  xaux     = new double[nInputs];
  funcID   = odata->numFuncEvals_++;
  outputID = odata->outputID_;
  for (i=0; i < nInputs; i++) xaux[i] = XValues[i];

  // modify values 
  if (coarseoptimtype == 1)
  {
    /* Computing pxk + Bk * (S - xk) */
    for (i = 0; i < nInputs; i++)
    {
      aux[i] = 0.0;
      for (j = 0; j < nInputs; j++)
        aux[i] = aux[i] + B0[i + j*nInputs]*(xaux[j]-x0[j]);
    };
    for (i = 0; i < nInputs; i++) xaux[i] = px0[i] + aux[i];
  }

  // ------ run simulation ------
  currDriver = odata->funcIO_->getDriver();
  odata->funcIO_->setDriver(2);      /* we always optimize the coarse model */
  odata->funcIO_->evaluate(funcID,nInputs,xaux,nOutputs,localY,0);
  odata->funcIO_->setDriver(currDriver);

  // ------ compute objective function ------
  (*YValue) = 0.0;
  if (coarseoptimtype == 1)
  {
     for (i=0; i<nOutputs; i++)
       (*YValue) += (localY[i] - y[i])*(localY[i] - y[i]);
     (*YValue) = sqrt((*YValue))/normy*100.00;
  }
  else
  {
    for (i=0; i<nOutputs; i++)
      (*YValue) += (localY[i] - fx[i])*(localY[i] - fx[i]);
    (*YValue) = sqrt((*YValue));
  }

  // ------ update current optimal settings ------
  // ------ output optimization information ------
  if ((*YValue) < odata->optimalY_)
  {
    odata->optimalY_ = (*YValue);
    for (i = 0; i < nInputs; i++) odata->optimalX_[i] = XValues[i];
    if (odata->outputLevel_ > 1)
    {
      printDashes(PL_INFO, 0);
      for (i = 0; i < nInputs; i++)
        printf("SMOptimizer %6d : X %3d = %16.8e\n",odata->numFuncEvals_,
               i+1, XValues[i]);
      printf("SMOptimizer %6d : Ymin = %16.8e\n",
              odata->numFuncEvals_, odata->optimalY_);
      printf("SMOptimizer %6d : output 1 = %16.8e\n",
                odata->numFuncEvals_, localY[0]);
      printDashes(PL_INFO, 0);
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
SMOptimizer::SMOptimizer()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SMOptimizer::~SMOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void SMOptimizer::optimize(oData *odata)
{
  int    nInputs, nOutputs, i, j, kk, maxfun=1000;
  int    nIter = 1, nIterMAX = 10, ncevals = 0, nfevals = 0;
  double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp;
  double *zstar, *px1, *B1, *x1, *B0h;
  double Fx, h = 1e50, tolX = 1e-4;
  FILE   *fspecs;
  extern double *px0, *B0, *x0, *fx, *y, normy;
  extern int coarseoptimtype;


  nInputs = odata->nInputs_;
  nOutputs = odata->nOutputs_;
  if(nOutputs <= 0)
  {
    printf("nOutputs in file %s line %d has nOutputs <= 0 (%d)",
           __FILE__, __LINE__, nOutputs);
    exit(1);
  }
  for (i = 0; i < nInputs; i++) odata->optimalX_[i] = 0.0;
  odata->numFuncEvals_ = 0;
  odata->optimalY_     = 1.0e50;
  tolX                 = odata->tolerance_;

  if (!strcmp(odata->targetFile_,"NONE"))
  {
    printf("SMOptimizer ERROR : no target file.\n");
    y = new double[nOutputs];
    for (i = 0; i < nOutputs; i++) y[i] = 0;
  }
  else
  {
    fspecs = fopen(odata->targetFile_, "r");
    y = new double[nOutputs];
    for (i = 0; i < nOutputs; i++) fscanf(fspecs, "%lg", &y[i]);
    fclose(fspecs);
  }
  normy = 0.0;
  for (i=0; i < nOutputs; i++) normy = normy+ y[i]*y[i];
  normy = sqrt(normy);
  if (normy == 0) normy = 1;

  printf("\n\n");
  printAsterisks(PL_INFO, 0);
  printf(" SPACE MAPPING: Specifications: \n");
  for(i=0; i<nOutputs; i++)
  printf("                                      %20.14f\n", y[i]);
  printAsterisks(PL_INFO, 0);

  // sm auxiliary variables: initialization + storing in file smaux 
  XValues = new double[nInputs+1];
  fx      = new double[nOutputs];
  x0      = new double[nInputs+1];
  px0     = new double[nInputs+1];
  B0      = new double[nInputs*nInputs];
  x1      = new double[nInputs+1];
  px1     = new double[nInputs+1];
  B1      = new double[nInputs*nInputs];
  B0h     = new double[nInputs];
  zstar   = new double[nInputs];

  /* Bobyqa STOP Criteria   */
  rhobeg = odata->upperBounds_[0] - odata->lowerBounds_[0];
  for (i = 0; i < nInputs; i++)
  {
    dtemp = odata->upperBounds_[i] - odata->lowerBounds_[i];
    if (dtemp > rhobeg) rhobeg = dtemp;
  }
  rhobeg *= 0.05;
  rhoend = rhobeg * odata->tolerance_;
  if (rhobeg < rhoend)
  {
    printf("SMOptimizer WARNING: tolerance too large.\n");
    printf("                      tolerance reset to 1.0e-6.\n");
    rhoend = rhobeg * 1.0e-6;
  }

  /* Initializing state variables */
  for (i = 0; i < nInputs; i++)
  {
    x0[i]  = 0.0;
    px0[i] = 0.0;
    for (j = 0; j < nInputs; j++)
    {
      if (i==j) B0[i + j*nInputs] = 1.0;
      else      B0[i + j*nInputs] = 0.0;
    }
  }

  /* next is optimization => 1 */
  //savesmaux(1, nInputs, px0, B0, x0);


  /* I. INITIAL GUESS */
  for (i = 0; i < nInputs; i++) XValues[i] = odata->initialX_[i];

  printAsterisks(PL_INFO, 0);
  printf(" SPACE MAPPING: Initial Guess: \n");
  for(i=0; i<nInputs; i++)
  printf("                                      %20.14f\n", XValues[i]);
  printAsterisks(PL_INFO, 0);

  /* 0. COARSE MODEL OPTIMIZATION */
  // ------ call optimizer ------

#ifdef HAVE_BOBYQA
  psSMBobyqaObj_ = (void *) odata;
  maxfun = 1000;
  int    printLevel=5555;
  int    nPts = (nInputs + 1) * (nInputs + 2) / 2;
  kk = (nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1;
  double *workArray = new double[kk];
  bobyqa_(&nInputs, &nPts, XValues, odata->lowerBounds_, odata->upperBounds_,
          &rhobeg, &rhoend, &printLevel, &maxfun, workArray);
  delete [] workArray;
#else
  printf("ERROR : Bobyqa optimizer not installed.\n");
  exit(1);
#endif
  ncevals += maxfun;
  coarseoptimtype = 0;

  for (i=0;i<nInputs;i++)
  {
    x0[i]    = XValues[i];
    zstar[i] = XValues[i];
  }

  printf("\n\n");
  printAsterisks(PL_INFO, 0);
  printf(" SPACE MAPPING: Coarse model optimum: \n");
  for(i=0; i<nInputs; i++)
  printf("                                      %20.14f\n", x0[i]);
  printAsterisks(PL_INFO, 0);
  printf("\n\n");

  /* 1. EVALUATING THE FINE MODEL AT THE COARSE OPTIMUM */
  odata->funcIO_->setDriver(1);      /* we evaluate the fine model */
  odata->funcIO_->evaluate(odata->numFuncEvals_++,nInputs,x0,nOutputs,fx,0);

  Fx = 0.0;
  for (i=0; i<nOutputs; i++)
    Fx += (fx[i] - y[i])*(fx[i] - y[i]);
  Fx = sqrt(Fx)/normy*100.00;

  nfevals += 1;
  printf("   Iteration      # f evals   # c eval   Fine Cost Function    ||x_{k+1} - x_k||\n");
  printf(" ---------------------------------------------------------------------------------\n");
  printf("   %5d           %5d       %5d         %8.3f\n",nIter,nfevals,
         ncevals,Fx);

  // Now f(x0) is stored in the file fileaux (format: m f(x0)) 
  /* 2. EVALUATING THE SPACE-MAPPING FUNCTION AT x0            */
  for (i=0;i<nInputs;i++)
    px0[i] = x0[i];         /* makes sense for parameter extraction */

#ifdef HAVE_BOBYQA
  maxfun=1000;
  printLevel=5555;
  nPts = (nInputs + 1) * (nInputs + 2) / 2;
  workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];
  bobyqa_(&nInputs, &nPts, px0, odata->lowerBounds_, odata->upperBounds_,
          &rhobeg, &rhoend, &printLevel, &maxfun, workArray);
  delete [] workArray;
#else
  printf("ERROR : Bobyqa optimizer not installed.\n");
  exit(1);
#endif
  ncevals += maxfun;
  coarseoptimtype = 1;

  while (1)
  {
    /* 2. OPTIMIZING THE MAPPED COARSE MODEL          */
    for (i=0;i<nInputs;i++) x1[i] = x0[i]; /* for parameter extraction */

#ifdef HAVE_BOBYQA
    maxfun=1000;
    printLevel=5555;
    nPts = (nInputs + 1) * (nInputs + 2) / 2;
    workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];
    bobyqa_(&nInputs, &nPts, x1, odata->lowerBounds_, odata->upperBounds_,
            &rhobeg, &rhoend, &printLevel, &maxfun, workArray);
    delete [] workArray;
#else
    printf("ERROR : Bobyqa optimizer not installed.\n");
    exit(1);
#endif
    ncevals += maxfun;
    coarseoptimtype = 0;

    /* 3. EVALUATING THE FINE MODEL AT x1 */
    odata->funcIO_->setDriver(1);      /* we evaluate the fine model */
    odata->funcIO_->evaluate(odata->numFuncEvals_++,nInputs,x1,nOutputs,fx,0);
    Fx   = 0.0;
    for (i=0; i<nOutputs; i++)
      Fx   += (fx[i] - y[i])*(fx[i] - y[i]);
    Fx   = sqrt(Fx)/normy*100.00;
    nfevals += 1;
    /*printf("                                      %20.14f\n", Fx);*/
    /* Now f(x1) is stored in the file fileaux (format: m f(x1)) */

    /* S. STOPPING CRITERIA              */
    h = 0.0;
    for (i=0;i<nInputs;i++) h += (x1[i]-x0[i])*(x1[i]-x0[i]);
    h = sqrt(h);
    /*printf(" h = %10.6f\n", h); */
    printf("   %5d           %5d       %5d         %8.3f              %4.2e \n", 
           nIter+1,nfevals,ncevals,Fx,h);
    if (nIter == nIterMAX || h < tolX) break;

    /* 4. EVALUATING THE SPACE-MAPPING FUNCTION AT x1            */
    for (i=0;i<nInputs;i++) px1[i] = px0[i];/* for parameter extraction */
#ifdef HAVE_BOBYQA
    maxfun=1000;
    printLevel=5555;
    nPts = (nInputs + 1) * (nInputs + 2) / 2;
    workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];
    bobyqa_(&nInputs, &nPts, px1, odata->lowerBounds_, odata->upperBounds_,
            &rhobeg, &rhoend, &printLevel, &maxfun, workArray);
    delete [] workArray;
#else
    printf("ERROR : Bobyqa optimizer not installed.\n");
    exit(1);
#endif
    ncevals += maxfun;
    coarseoptimtype = 1;

    nIter = nIter + 1;

    /* 5. UPDATING THE APPROX. OF p */
    /*    B1 = B0 + (px1-px0-B0*h)/(h'*h)*h' */
    for (i=0;i<nInputs;i++)
    {
      B0h[i] = 0.0;
      for (j=0;j<nInputs;j++)
        B0h[i] += B0[i+j*nInputs]*(x1[j]-x0[j]);
    }
    for (i=0;i<nInputs;i++)
      for (j=0;j<nInputs;j++)
        B1[i+j*nInputs] = B0[i+j*nInputs] + 
                          (px1[i]-px0[i]-B0h[i])/(h*h)*(x1[j]-x0[j]);

/* - CHECKINGS ------------------------------------------------------ */
/*
    for (i=0;i<nInputs;i++)
    {
      for (j=0;j<nInputs;j++)
        printf("%10.6f   ", B0[i+j*nInputs]);
      printf("\n");
    }
    printf("\n");
    for (i=0;i<nInputs;i++)
    {
      for (j=0;j<nInputs;j++) printf("%10.6f   ", B1[i+j*nInputs]);
      printf("\n");
    }
    printf("\n");
    for (j=0;j<nInputs;j++) printf("%10.6f   ", x0[j]);
    printf("\n\n");
    for (j=0;j<nInputs;j++) printf("%10.6f   ", x1[j]);
    printf("\n\n");
    for (j=0;j<nInputs;j++) printf("%10.6f   ", px0[j]);
    printf("\n\n");
    for (j=0;j<nInputs;j++) printf("%10.6f   ", px1[j]);
    printf("\n\n");
    for (j=0;j<nInputs;j++) printf("%10.6f   ", B0h[j]);
    printf("\n\n");
*/

    /* updating variables */
    for (i=0;i<nInputs;i++)
    {
      x0[i]  = x1[i];
      px0[i] = px1[i];
      for (j=0;j<nInputs;j++)
        B0[i+j*nInputs] = B1[i+j*nInputs];
    }
  }

  /* DIAGNOSTICS */
  if (nIter == nIterMAX)
    printf("\n Maximum number of %4d SM iterations reached.\n\n", nIterMAX);
  else
    printf("\n Normal termination.\n\n");

  printf("\n\n");
  printAsterisks(PL_INFO, 0);
  printf(" SPACE MAPPING: Space-mapping solution:\n");
  for(i=0; i<nInputs; i++)
  printf("                                      %20.14f\n", x1[i]);
  printAsterisks(PL_INFO, 0);
  printf("\n\n");
  odata->optimalY_ = Fx;
  for (i = 0; i < nInputs; i++) odata->optimalX_[i] = x1[i];

  // clean up
  delete [] XValues;
  delete [] fx;
  delete [] x0;
  delete [] px0;
  delete [] B0;
  delete [] x1;
  delete [] px1;
  delete [] B1;
  delete [] B0h;
  delete [] zstar;
  delete [] y;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
SMOptimizer& SMOptimizer::operator=(const SMOptimizer &)
{
  printf("SMOptimizer operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

