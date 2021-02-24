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
// Functions for the class Kriging
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
#ifdef WINDOWS
#include <windows.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "Kriging.h"
#include "sysdef.h"
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
  void dgeqrf_(int *, int *, double *, int *, double *, double *, int *, int *);
  void dpotrf_(char *, int *, double *, int *, int *);
  void dpotrs_(char *, int *, int *, double *, int *, double *,
               int *, int *);
  void kbobyqa_(int *,int *, double *, double *, double *, double *,
               double *, int *, int *, double*);
  void newuoa_(int *,int *,double *,double *,double *,int *,
               int *,double*);
  void dormqr_(char *, char *, int *, int *, int *, double *, int *, 
               double *, double *, int *, double *, int *, int *);
}

// ************************************************************************
// ************************************************************************
// internal global variables
// ------------------------------------------------------------------------
int    KRI_outputLevel=0;
int    KRI_iter=-1;
int    KRI_nInputs=-1;
int    KRI_nSamples=-1;
int    KRI_pOrder=-1;
double *KRI_XDists=NULL;
double *KRI_SMatrix=NULL;
double *KRI_FMatrix=NULL;
double *KRI_FMatTmp=NULL;
double *KRI_MMatrix=NULL;
double *KRI_X=NULL;
double *KRI_Y=NULL;
double KRI_currY=0.0;
double *KRI_Ytmp=NULL;
double KRI_OptY=0.0;
double *KRI_OptThetas=NULL;
double *KRI_dataStdDevs=NULL;
double KRI_Exponent=2.0;
double KRI_YStd=1.0;
double KRI_nugget=0.0;
int    KRI_terminate=0;
int    KRI_noProgressCnt=0;

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif
   void *kribobyqaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    nBasis, count, Cleng, ii, jj, kk, status;
      double dist, ddata, ddata2;
      char   uplo='L';
      FILE   *fp=NULL;

      for (ii = 0; ii < KRI_nInputs; ii++)
      {
         if (XValues[ii] < 0) 
         {
            (*YValue) = PSUADE_UNDEFINED;
            return NULL;
         }
      }
      if (KRI_noProgressCnt >= 200)
      {
         if (KRI_outputLevel >= 1)
         {
            //KRI_outputLevel = 0;
            printf("\t*** no progress for more than 200 iterations.\n");
         }
         (*YValue) = PSUADE_UNDEFINED;
         return NULL;
      }
      fp = fopen("psuade_stop", "r");
      if (fp != NULL)
      {
         fclose(fp);
         printf("Kriging: psuade_stop file found - terminating ...\n");
         (*YValue) = KRI_OptY;
         return NULL;
      }
      KRI_iter++;
      nBasis = 1;
      if (KRI_pOrder == 1) nBasis = KRI_nInputs + 1;
      Cleng = KRI_nSamples + nBasis;

      count = 0;
      for (jj = 0; jj < KRI_nSamples; jj++)
      {
         KRI_SMatrix[jj*KRI_nSamples+jj] = 1.0 + 1e-15 * KRI_nSamples;
         if (KRI_nugget != 0.0) 
            KRI_SMatrix[jj*KRI_nSamples+jj] += KRI_nugget;
         else if (KRI_nInputs == 2) 
            KRI_SMatrix[jj*KRI_nSamples+jj] += 0.01;

         if (KRI_dataStdDevs != NULL) 
         {
            ddata = pow(KRI_dataStdDevs[jj]/KRI_YStd,KRI_Exponent);
            KRI_SMatrix[jj*KRI_nSamples+jj] += ddata;
         }
         for (kk = jj+1; kk < KRI_nSamples; kk++)
         {
            dist = 0.0;
            for (ii = 0; ii < KRI_nInputs; ii++)
               dist += pow(KRI_XDists[count*KRI_nInputs+ii]/XValues[ii],
                           KRI_Exponent);
            dist = exp(-dist);
            if (dist < 1.0e-16) dist = 0;
            KRI_SMatrix[jj*KRI_nSamples+kk] = dist;
            KRI_SMatrix[kk*KRI_nSamples+jj] = dist;
            count++;
         }
      }

      for (jj = 0; jj < KRI_nSamples; jj++)
      {
         KRI_FMatrix[jj] = 1.0;
         KRI_FMatTmp[jj] = 1.0;
         if (KRI_pOrder == 1)
         {
            for (kk = 1; kk < KRI_nInputs+1; kk++)
            {
               ddata = KRI_X[jj*KRI_nInputs+kk-1];
               KRI_FMatTmp[kk*KRI_nSamples+jj] = ddata;
               KRI_FMatrix[kk*KRI_nSamples+jj] = ddata;
            }
         }
      }
      dpotrf_(&uplo, &KRI_nSamples, KRI_SMatrix, &KRI_nSamples, &status);
      kk = nBasis;
      dpotrs_(&uplo, &KRI_nSamples, &kk, KRI_SMatrix, &KRI_nSamples,
              KRI_FMatTmp, &KRI_nSamples, &status);
      for (jj = 0; jj < nBasis; jj++)
      {
         for (kk = 0; kk < nBasis; kk++)
         {
            ddata = 0.0;
            for (ii = 0; ii < KRI_nSamples; ii++)
               ddata += KRI_FMatrix[ii+jj*KRI_nSamples]*
                        KRI_FMatTmp[ii+kk*KRI_nSamples];
            KRI_MMatrix[jj+kk*nBasis] = ddata;
         }
      }
      if (nBasis > 1) 
         dpotrf_(&uplo,&nBasis,KRI_MMatrix,&nBasis,&status);
      for (jj = 0; jj < KRI_nSamples; jj++) KRI_Ytmp[jj] = KRI_Y[jj];
      kk = 1;
      dpotrs_(&uplo,&KRI_nSamples,&kk,KRI_SMatrix,&KRI_nSamples,KRI_Ytmp,
              &KRI_nSamples, &status);
      for (jj = 0; jj < nBasis; jj++)
      {
         ddata = 0.0;
         for (ii = 0; ii < KRI_nSamples; ii++)
            ddata += KRI_FMatrix[ii+jj*KRI_nSamples]*KRI_Ytmp[ii];
         KRI_Ytmp[KRI_nSamples+jj] = ddata;
      }
      if (nBasis == 1)
      {
         if (KRI_MMatrix[0] == 0)
         {
            printf("Kriging ERROR: divide by 0.\n");
            exit(1);
         }
         KRI_Ytmp[KRI_nSamples] /= KRI_MMatrix[0];
      }
      else
      {
         kk = 1;
         dpotrs_(&uplo,&kk,&kk,KRI_MMatrix,&nBasis,
                 &KRI_Ytmp[KRI_nSamples], &kk, &status);
      }
      for (jj = 0; jj < KRI_nSamples; jj++)
      {
         ddata = 0.0;
         for (ii = 0; ii < nBasis; ii++)
            ddata += KRI_FMatrix[jj+ii*KRI_nSamples]*
                     KRI_Ytmp[ii+KRI_nSamples];
         KRI_Ytmp[jj] = KRI_Y[jj] - ddata;
      }
      kk = 1;
      dpotrs_(&uplo,&KRI_nSamples,&kk,KRI_SMatrix,&KRI_nSamples,
              KRI_Ytmp,&KRI_nSamples, &status);

      ddata = 0.0;
      for (jj = 0; jj < KRI_nSamples; jj++)
         ddata += KRI_Ytmp[jj] * KRI_Y[jj];
      ddata /= (double) KRI_nSamples;

      for (jj = 0; jj < KRI_nSamples; jj++) 
      {
         ddata2 = KRI_SMatrix[KRI_nSamples*jj+jj];
         if (ddata2 < 0.0) ddata2 = - ddata2;
         ddata *= pow(ddata2, 2.0/(double) KRI_nSamples);
      }

      (*YValue) = KRI_currY = ddata;

      if (ddata < KRI_OptY || KRI_iter == 1)
      {
         KRI_noProgressCnt = 0;
         KRI_OptY = ddata;
         for (ii = 0; ii < KRI_nInputs; ii++)
            KRI_OptThetas[ii] = XValues[ii];
         if (KRI_outputLevel > 1)
         {
            printf("\t Kriging : iteration %d\n",KRI_iter);
            for (ii = 0; ii < KRI_nInputs; ii++)
               printf("\t    Current best theta %d = %e\n",ii+1,XValues[ii]);
            printf("\t    Current best objective value = %e\n",ddata);
         }
      }
      else KRI_noProgressCnt++;

      if (psRSExpertMode_ == 1 && KRI_outputLevel > 3)
      {
         printf("\t Kriging : iteration %d\n", KRI_iter);
         for (ii = 0; ii < KRI_nInputs; ii++)
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
Kriging::Kriging(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
  int    ii, jj;
  double ddata;
  char   pString[500], winput[500], winput2[500], fname[500], *strPtr;
  FILE   *fp;

  // initialize variables and parameters
  faID_ = PSUADE_RS_KR;
  XNormalized_ = NULL;
  YNormalized_ = NULL;
  Rmatrix_ = NULL;
  Mmatrix_ = NULL;
  V1_ = NULL;
  V2_ = NULL;
  pOrder_  = 0;
  initFlag_ = 0;
  workLength_ = 0;
  workArray_ = NULL;
  workX_ = NULL;
  dataStdDevs_ = NULL;
  optTolerance_ = 1.0e-4;
  fastMode_ = 3;
  Thetas_ = new double[nInputs_+1];
  checkAllocate(Thetas_, "Thetas_ in Kriging::constructor");

  for (ii = 0; ii <= nInputs_; ii++) Thetas_[ii] = 0.01;
  noReuse_ = 0;
  betas_ = NULL;
  gammas_ = NULL;
  betasOpt_ = NULL;
  gammasOpt_ = NULL;
  thetasOpt_ = NULL;

  // display banner and additonal information
  if (psInteractive_ == 1)
  {
    printAsterisks(PL_INFO, 0);
    printf("*                Kriging Analysis\n");
    printf("* Set printlevel to 1-4 to see Kriging details.\n");
    printf("* Turn on rs_expert mode to set slow or fast mode.\n");
    printf("*  + Fast mode: no optimization of hyperparameters.\n");
    printf("*      - turn on rs_expert to set hyperparameters.\n");
    printf("*      - default values = 1.0\n");
    printf("*  + Slow mode : hyperparameters are optimized.\n");
    printf("*      - to change optimization parameters, turn\n");
    printf("*        rs_expert mode.\n");
    printf("*  + Snail mode (DEFAULT): use multi-start optimization.\n");
    printf("*      - to change optimization parameters, turn\n");
    printf("*        rs_expert mode.\n");
    printf("* Create 'psuade_stop' file to gracefully terminate.\n");
    printf("* Create 'ps_print' file to set print level on the fly.\n");
    printEquals(PL_INFO, 0);
  }

  // read configure file, if any 
  if (psRSExpertMode_ == 0 && psConfig_ != NULL)
  {
    strPtr = psConfig_->getParameter("KRI_mode");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %d", winput, winput2, &ii);
      if (ii < 0 || ii > 3)
      {
        printf("Kriging INFO: mode from config not valid.\n");
        printf("              mode kept at %d.\n", fastMode_);
      }
      else
      {
        fastMode_ = ii;
        printf("Kriging INFO: mode from config = %d.\n",fastMode_);
      }
    }
    strPtr = psConfig_->getParameter("KRI_tol");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %lg", winput, winput2, &optTolerance_);
      if (optTolerance_ < 0.0 || optTolerance_ >= 1.0)
      {
        optTolerance_ = 1.0e-4;
        printf("Kriging INFO: tolerance from config not valid.\n");
        printf("              tolerance kept at %e.\n", optTolerance_);
      }
      else
      {
        printf("Kriging INFO: tolerance from config = %e.\n",
               optTolerance_);
      }
    }
    strPtr = psConfig_->getParameter("KRI_LENG_SCALE");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %d %s %lg", winput, &ii, winput2, &ddata);
      if (ii < 1 || ii > nInputs_)
      {
        printf("Kriging INFO: invalid input number for length scale.\n");
        printf("              Input number read = %d.\n", ii);
      }
      else
      {
        Thetas_[ii-1] = ddata;
        printf("Kriging INFO: length scale for input %d set to %e.\n",
               ii, ddata);
      }
    }
    strPtr = psConfig_->getParameter("KRI_DATA_STDEV_FILE");
    if (strPtr != NULL)
    {
      sscanf(strPtr, "%s %s %s", winput, winput2, fname);
      fp = fopen(fname, "r");
      if (fp == NULL)
      {
        printf("Kriging INFO: data variance file not found.\n");
      }
      else
      {
        fscanf(fp, "%d", &ii); 
        if (ii != nSamples_)
        {
          printf("Kriging ERROR: stdev file should have %d entries.\n",
                 nSamples_);
          fclose(fp);
        }
        else
        {
          dataStdDevs_ = new double[nSamples_];
          checkAllocate(dataStdDevs_,"dataStdDevs_ in Kriging::constructor");
          for (ii = 0; ii < nSamples_; ii++)
          {
            fscanf(fp, "%d %lg", &jj, &dataStdDevs_[ii]); 
            if (jj != ii+1)
            {
              printf("Kriging ERROR: line %d in the std dev file.\n",jj+1);
              delete [] dataStdDevs_;
              dataStdDevs_ = NULL;
              break;
            }
          } 
          fclose(fp);
          if (dataStdDevs_ != NULL)
          {
            printf("Kriging INFO: std. dev. file has been read.\n");
            KRI_dataStdDevs = dataStdDevs_; 
          }
        }
      }
    }
  }
  if (psInteractive_ == 1)
  {
    fp = fopen("psuade_kriging_optdata","r");
    if (fp != NULL)
    {
      printf("Kriging: psuade_kriging_optdata file found.\n");
      sprintf(pString,
              "Use Kriging length scales from the file? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
        for (ii = 0; ii < nInputs_; ii++) fscanf(fp,"%lg",&Thetas_[ii]); 
      fclose(fp);
      fastMode_ = 1;
    }
  }
   
  // if configure file is not used, ask for parameters if
  // response surface expert mode is on 
  if (psRSExpertMode_ == 1 && psInteractive_ == 1)
  {
    printf("There are three modes available: \n");
    printf("(1) fast mode with pre-specified thetas\n");
    printf("(2) slow mode with optimization on the thetas\n");
    printf("(3) very slow mode with multi-start optimization\n");
    //printf("(4) another fast mode using another optimization\n");
    sprintf(pString, "Please select mode (1 - 3) : ");
    fastMode_ = getInt(1,3,pString);
    if      (fastMode_ == 4) fastMode_ = 0;
    else if (fastMode_ == 1)
    {
      printf("Kriging: Length scales are correlation lengths\n");
      printf("         in the random parameter space.\n");
      printf("Current initial length scales are:\n");
      for (ii = 0; ii < nInputs_; ii++)
         printf("    Input %d: %e\n", ii+1, Thetas_[ii]);
      sprintf(pString, "Change length scales (thetas)? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
        for (ii = 0; ii < nInputs_; ii++)
        {
          sprintf(pString,"Enter theta for input %d (>0): ", ii+1);
          Thetas_[ii] = getDouble(pString);
          if (Thetas_[ii] <= 0.0)
          {
            printf("ERROR: theta <= 0 not valid.\n");
            exit(1);
          }
          if (ii == 0)
          {
            sprintf(pString,"Use %e for all other thetas? (y or n) ",
                    Thetas_[0]);
            getString(pString, winput);
            if (winput[0] == 'y')
            {
              for (jj = 1; jj < nInputs_; jj++) Thetas_[jj] = Thetas_[0];
              break;
            }
          }
        }
      }
    }
    else
    {
      sprintf(pString, "Enter optimization tolerance (default = 1e-4): ");
      optTolerance_ = getDouble(pString);
      if (optTolerance_ <= 0 || optTolerance_ > 0.5)
      {
        printf("Kriging INFO: optimization tol should be in (0,0.5]).\n");
        printf("              Tolerance set to default = 1.0e-4.\n");
        optTolerance_ = 1.0e-4;
      }
      if (fastMode_ == 2)
      {
        printf("Kriging: Current initial length scales (thetas) are:\n");
        for (ii = 0; ii < nInputs_; ii++)
          printf("     Input %d: %e\n", ii+1, Thetas_[ii]);
        printf("If some knowledge is available about the relative\n");
        printf("importance of some parameters, different initial\n");
        printf("thetas can be entered to reflect this knowledge \n");
        printf("(sensitive parameters have larger thetas.)\n");
        sprintf(pString, "Change initial thetas? (y or n) ");
        getString(pString, winput);
        if (winput[0] == 'y')
        {
          for (ii = 0; ii < nInputs_; ii++)
          {
            sprintf(pString,"Enter theta for input %d : ", ii+1);
            Thetas_[ii] = getDouble(pString);
            if (Thetas_[ii] <= 0.0)
              printf("Warning: theta < 0 not recommended.\n");
            if (ii == 0)
            {
              sprintf(pString,"Use %e for all other thetas? (y or n) ",
                      Thetas_[0]);
              getString(pString, winput);
              if (winput[0] == 'y')
              {
                for (jj = 1; jj < nInputs_; jj++) Thetas_[jj] = Thetas_[0];
                break;
              }
            }
          }
        }
      }
    }
    if (psMasterMode_ == 1)
    {
      sprintf(pString, "Add nugget? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y') 
      {
        KRI_nugget = 1.0;
        while (KRI_nugget >= 1.0 || KRI_nugget < 0.0)
        {
          sprintf(pString, "Enter nugget ([0,1)) : ");
          KRI_nugget = getDouble(pString);
        }
      }
    }
  }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Kriging::~Kriging()
{
  if (Thetas_ != NULL) delete [] Thetas_;
  if (XNormalized_ != NULL) delete [] XNormalized_;
  if (YNormalized_ != NULL) delete [] YNormalized_;
  if (Rmatrix_ != NULL) delete [] Rmatrix_;
  if (Mmatrix_ != NULL) delete [] Mmatrix_;
  if (V1_ != NULL) delete [] V1_;
  if (V2_ != NULL) delete [] V2_;
  if (workArray_ != NULL) delete [] workArray_;
  if (workX_ != NULL) delete [] workX_;
  if (dataStdDevs_ != NULL) delete [] dataStdDevs_;
  if (betas_ != NULL) delete [] betas_;
  if (betasOpt_ != NULL) delete [] betasOpt_;
  if (gammas_ != NULL) delete [] gammas_;
  if (gammasOpt_ != NULL) delete [] gammasOpt_;
  if (thetasOpt_ != NULL) delete [] thetasOpt_;
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int Kriging::initialize(double *X, double *Y)
{
  if (outputLevel_ >= 1) printf("Kriging training begins....\n");
  train(X,Y);
  if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  return 0;
}
  
// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Kriging::genNDGridData(double *XIn,double *YIn,int *NOut,double **XOut,
                           double **YOut)
{
  int totPts;

  if (outputLevel_ >= 1) printf("Kriging training begins....\n");
  train(XIn,YIn);
  if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  if ((*NOut) == -999) return 0;
 
  genNDGrid(NOut, XOut);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in Kriging::genNDGridData");

  if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
  predict(totPts, *XOut, *YOut, NULL);
  if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
  return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int Kriging::gen1DGridData(double *XIn,double *YIn,int ind1,
                           double *settings,int *NOut, double **XOut, 
                           double **YOut)
{
  int    ii, kk, totPts;
  double HX;
  psVector vecXT;

  if (outputLevel_ >= 1) printf("Kriging training begins....\n");
  train(XIn,YIn);
  if (outputLevel_ >= 1) printf("Kriging training completed.\n");
 
  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[totPts];
  checkAllocate(*XOut, "XOut in Kriging::gen1DGridData");
  vecXT.setLength(totPts*nInputs_);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (kk = 0; kk < nInputs_; kk++) 
      vecXT[ii*nInputs_+kk] = settings[kk]; 
    vecXT[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
    (*XOut)[ii] = HX * ii + lowerBounds_[ind1];
  }
    
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in Kriging::gen1DGridData");
  if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
  predict(totPts, vecXT.getDVector(), *YOut, NULL);
  if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int Kriging::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                           double *settings, int *NOut, double **XOut, 
                           double **YOut)
{
  int ii, jj, kk, totPts, index;
  psVector vecHX, vecXT;

  if (outputLevel_ >= 1) printf("Kriging training begins....\n");
  train(XIn, YIn);
  if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  
  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[2*totPts];
  checkAllocate(*XOut, "XOut in Kriging::gen2DGridData");
  vecXT.setLength(totPts*nInputs_);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      index = ii * nPtsPerDim_ + jj;
      for (kk = 0; kk < nInputs_; kk++) 
        vecXT[index*nInputs_+kk] = settings[kk]; 
      vecXT[index*nInputs_+ind1]  = vecHX[0] * ii + lowerBounds_[ind1];
      vecXT[index*nInputs_+ind2]  = vecHX[1] * jj + lowerBounds_[ind2];
      (*XOut)[index*2]   = vecHX[0] * ii + lowerBounds_[ind1];
      (*XOut)[index*2+1] = vecHX[1] * jj + lowerBounds_[ind2];
    }
  }
    
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in Kriging::gen2DGridData");
  if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
  predict(totPts, vecXT.getDVector(), *YOut, NULL);
  if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Kriging::gen3DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                           int ind3, double *settings, int *NOut, 
                           double **XOut, double **YOut)
{
  int ii, jj, ll, kk, totPts, index;
  psVector vecHX, vecXT;

  if (outputLevel_ >= 1) printf("Kriging training begins....\n");
  train(XIn, YIn);
  if (outputLevel_ >= 1) printf("Kriging training completed.\n");
 
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[3*totPts];
  checkAllocate(*XOut, "XOut in Kriging::gen3DGridData");
  vecXT.setLength(totPts*nInputs_);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      for (ll = 0; ll < nPtsPerDim_; ll++) 
      {
        index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
        for (kk = 0; kk < nInputs_; kk++) 
          vecXT[index*nInputs_+kk] = settings[kk]; 
        vecXT[index*nInputs_+ind1] = vecHX[0] * ii + lowerBounds_[ind1];
        vecXT[index*nInputs_+ind2] = vecHX[1] * jj + lowerBounds_[ind2];
        vecXT[index*nInputs_+ind3] = vecHX[2] * ll + lowerBounds_[ind3];
        (*XOut)[index*3]   = vecHX[0] * ii + lowerBounds_[ind1];
        (*XOut)[index*3+1] = vecHX[1] * jj + lowerBounds_[ind2];
        (*XOut)[index*3+2] = vecHX[2] * ll + lowerBounds_[ind3];
      }
    }
  }
   
  (*YOut) = new double[totPts];
  checkAllocate((*YOut), "YOut in Kriging::gen3DGridData");
  if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
  predict(totPts, vecXT.getDVector(), *YOut, NULL);
  if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Kriging::gen4DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                           int ind3, int ind4, double *settings, 
                           int *NOut, double **XOut, double **YOut)
{
  int ii, jj, ll, mm, kk, totPts, index;
  psVector vecHX, vecXT;

  if (outputLevel_ >= 1) printf("Kriging training begins....\n");
  train(XIn, YIn);
  if (outputLevel_ >= 1) printf("Kriging training completed.\n");
  
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  vecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[4*totPts];
  checkAllocate(*XOut, "XOut in Kriging::gen4DGridData");
  vecXT.setLength(totPts*nInputs_);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (jj = 0; jj < nPtsPerDim_; jj++) 
    {
      for (ll = 0; ll < nPtsPerDim_; ll++) 
      {
        for (mm = 0; mm < nPtsPerDim_; mm++) 
        {
          index = ii*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                  jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
          for (kk = 0; kk < nInputs_; kk++) 
            vecXT[index*nInputs_+kk] = settings[kk]; 
          vecXT[index*nInputs_+ind1] = vecHX[0] * ii + lowerBounds_[ind1];
          vecXT[index*nInputs_+ind2] = vecHX[1] * jj + lowerBounds_[ind2];
          vecXT[index*nInputs_+ind3] = vecHX[2] * ll + lowerBounds_[ind3];
          vecXT[index*nInputs_+ind4] = vecHX[3] * mm + lowerBounds_[ind4];
          (*XOut)[index*4]   = vecHX[0] * ii + lowerBounds_[ind1];
          (*XOut)[index*4+1] = vecHX[1] * jj + lowerBounds_[ind2];
          (*XOut)[index*4+2] = vecHX[2] * ll + lowerBounds_[ind3];
          (*XOut)[index*4+3] = vecHX[3] * mm + lowerBounds_[ind4];
        }
      }
    }
  }
    
  (*YOut) = new double[totPts];
  checkAllocate((*YOut), "YOut in Kriging::gen4DGridData");
  if (outputLevel_ >= 1) printf("Kriging interpolation begins....\n");
  predict(totPts, vecXT.getDVector(), *YOut, NULL);
  if (outputLevel_ >= 1) printf("Kriging interpolation completed.\n");
  (*NOut) = totPts;
  return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Kriging::evaluatePoint(double *X)
{
  int    iOne=1;
  double Y=0.0;
  predict(iOne, X, &Y, NULL);
  return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double Kriging::evaluatePoint(int npts, double *X, double *Y)
{
  predict(npts, X, Y, NULL);
  return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double Kriging::evaluatePointFuzzy(double *X, double &Ystd)
{
  int    iOne=1;
  double Y=0.0;
  predict(iOne, X, &Y, &Ystd);
  return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double Kriging::evaluatePointFuzzy(int npts, double *X, double *Y, 
                                   double *Ystds)
{
  predict(npts, X, Y, Ystds);
  return 0.0;
}

// ************************************************************************
// training 
// ------------------------------------------------------------------------
double Kriging::train(double *X, double *Y)
{
  int    ii, jj, kk, nBasis, count, status, nDists, Cleng;
  double *XDists, dist, ddata, *Cmatrix, *CFmatrix, *RFmatrix; 
  double *TUppers, *TLowers, *TValues;
  char   pString[500];

  // clean up 
  if (XNormalized_ != NULL) delete [] XNormalized_;
  if (YNormalized_ != NULL) delete [] YNormalized_;
  if (Rmatrix_ != NULL) delete [] Rmatrix_;
  if (Mmatrix_ != NULL) delete [] Mmatrix_;
  if (V1_ != NULL) delete [] V1_;
  if (V2_ != NULL) delete [] V2_;
  if (workArray_ != NULL) delete [] workArray_;
  if (workX_ != NULL) delete [] workX_;
  XNormalized_ = NULL;
  YNormalized_ = NULL;
  Mmatrix_ = NULL;
  Rmatrix_ = NULL;
  V1_ = NULL;
  V2_ = NULL;
  workArray_ = NULL;
  workX_ = NULL;

  // normalize the input and outputs 
  XNormalized_ = new double[nSamples_*nInputs_];
  checkAllocate(XNormalized_, "XNormalized_ in Kriging::train");
  initInputScaling(X, XNormalized_, 1);
  YNormalized_ = new double[nSamples_];
  checkAllocate(YNormalized_, "YNormalized_ in Kriging::train");
  initOutputScaling(Y, YNormalized_);
  KRI_YStd = YStd_;

  // compute distances between all pairs of inputs (nDists, XDists)
  if (fastMode_ != 0) 
  {
    status = computeDistances(&XDists, &nDists);
    if (status != 0)
    {
      printf("Kriging INFO: since there are repeated sample points.\n");
      printf("    Kriging will continue in fast mode taking all\n");
      printf("    length scales to be 1.\n");
      printf("    This may not give a good quality response surface.\n");
      printf("    So if you want to prune the sample and do it again,\n");
      printf("    terminate now. I am giving you 30 seconds to decide.\n");
      fastMode_ = 0;
      delete [] XDists;
      for (ii = 0; ii < nInputs_; ii++) Thetas_[ii] = 1.0;
#ifdef WINDOWS
      Sleep(20000);
#else
      sleep(20);
#endif
      printf("    10 more seconds.\n");
#ifdef WINDOWS
      Sleep(10000);
#else
      sleep(10);
#endif
    }
  }
  if (KRI_nugget != 0.0 && outputLevel_ > 0) 
    printf("Kriging INFO: nugget = %e\n",KRI_nugget);

  // fast mode = 0: nothing needs to be done
  if (fastMode_ == 0)
  {
    optimize();
    if (outputLevel_ >= 1) 
    {
      printf("Optimal length scales are:\n");
      for (ii = 0; ii < nInputs_; ii++)
        printf("     Input %d : %e\n", ii+1, Thetas_[ii]);
    }
    return 0.0;
  }
  // slower mode
  else if (fastMode_ == 1)
  {
    if (outputLevel_ > 0)
    {
      printEquals(PL_INFO,0);
      printf("Kriging training (1) begins.... (order = %d)\n",pOrder_);
      printf("Current length scales are:\n");
      for (ii = 0; ii < nInputs_; ii++)
        printf("     Input %d : %e\n", ii+1, Thetas_[ii]);
    }
    if (nSamples_ <= nInputs_ + 1) pOrder_ = 0;
    if (pOrder_ == 1) nBasis = nInputs_ + 1;
    else              nBasis = 1;
    if (outputLevel_ > 0) 
    {
      printf("Kriging training (1) ends.\n");
      printEquals(PL_INFO,0);
    }
  }
  else
  // slow mode: optimize
  {
    int    maxfun, pLevel, nPts, iOne=1, iZero=0, nSamOpt;
    int    *samStates, mode;
    double rhobeg=1.0, rhoend=1.0e-4;
    double *work, *samInputs, *samOutputs, *optThetas, *MMatrix=NULL;
    double optY=PSUADE_UNDEFINED, *SMatrix=NULL, *FMatrix=NULL;
    double *FMatTmp;
    FILE   *fp=NULL;
    Sampling *sampler;

    fp = fopen("psuade_stop", "r");
    if (fp != NULL)
    {
      printf("Kriging ERROR: remove the 'psuade_stop' file\n");
      printf("               first and re-do.\n");
      fclose(fp);
      exit(1);
    }

    if (nSamples_ <= nInputs_ + 1) pOrder_ = 0;
    if (pOrder_ == 1) nBasis = nInputs_ + 1;
    else              nBasis = 1;
    Cleng = nSamples_ + nBasis;

    TUppers = new double[nInputs_+1];
    TLowers = new double[nInputs_+1];
    checkAllocate(TLowers, "TLowers in Kriging::train");
    for (ii = 0; ii < nInputs_; ii++)
    {
      if (XMeans_[ii] == 0 && XStds_[ii] == 1.0)
      {
         TUppers[ii] = 30.0 * (upperBounds_[ii]-lowerBounds_[ii]);;
         TLowers[ii] = 0.1 * (upperBounds_[ii]-lowerBounds_[ii]);;
      }
      else
      {
         TUppers[ii] = 30.0;
         TLowers[ii] = 0.1;
      }
      if (psMasterMode_ == 1 && psInteractive_ == 1)
      {
        printf("Kriging: for input %d :\n", ii+1);
        printf("Kriging: current optimization lower bound = %e\n",
               TLowers[ii]);
        sprintf(pString,
           "Kriging: Enter new optimization lower bound : ");
        TLowers[ii] = getDouble(pString);
        if (TLowers[ii] <= 0.0)
        {
          printf("Kriging ERROR: lower bound <= 0\n");
          exit(1);
        }
        printf("Kriging: current optimization upper bound = %e\n",
               TUppers[ii]);
        sprintf(pString,
           "Kriging: Enter optimization upper bound : ");
        TUppers[ii] = getDouble(pString);
        if (TLowers[ii] > TUppers[ii])
        {
          printf("Kriging ERROR: lower bound >= upper bound\n");
          exit(1);
        }
      }
    }
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
    jj = (nPts+13) * (nPts+nInputs_) + 3*nInputs_*(nInputs_+3)/2;
    kk = (nPts+5)*(nPts+nInputs_)+3*nInputs_*(nInputs_+5)/2+1;
    if (jj > kk) work = new double[jj];
    else         work = new double[kk];
    checkAllocate(work, "work in Kriging::train");
    if (outputLevel_ > 0)
    {
      printEquals(PL_INFO, 0);
      printf("* Kriging optimization tolerance = %e\n", rhoend);
    }

    if (fastMode_ == 2)
    {
      if (outputLevel_ >= 1) 
        printf("Kriging training (2) begins.... (order = %d)\n",pOrder_);
      nSamOpt = 1;
      samInputs  = new double[nSamOpt * nInputs_];
      checkAllocate(samInputs, "samInputs_ in Kriging::train");
      for (ii = 0; ii < nInputs_; ii++) samInputs[ii] = Thetas_[ii];
    }
    else
    {
      if (outputLevel_ >= 1) 
        printf("Kriging training (3) begins.... (order = %d)\n",pOrder_);
      mode = 1;
      nSamOpt = 10;
      if (psGMMode_ == 1 && psInteractive_ == 1)
      {
        printf("Kriging: slow mode with multi-start optimization.\n");
        printf("Choose sampling method to generate multi-start.\n");
        sprintf(pString, "Sampling method (1-QMC, 2-LHS, 3-FF) : ");
        mode = getInt(1,3,pString);
        if (mode != 3)
        {
          sprintf(pString, "Enter sample size (1 - 100) : ");
          nSamOpt = getInt(1,100,pString);
        } 
      }
      if (mode == 1)
         sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      else if (mode == 2)
         sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
      else
         sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FF4);
 
      sampler->setPrintLevel(0);
      sampler->setInputBounds(nInputs_, TLowers, TUppers);
      sampler->setOutputParams(iOne);
      sampler->setSamplingParams(nSamOpt+2, iOne, iZero);
      sampler->initialize(0);
      nSamOpt = sampler->getNumSamples();
      samInputs  = new double[nSamOpt * nInputs_];
      samOutputs = new double[nSamOpt];
      samStates  = new int[nSamOpt];
      checkAllocate(samStates, "samStates in Kriging::train");
      sampler->getSamples(nSamOpt, nInputs_, iOne, samInputs,
                          samOutputs, samStates);
      nSamOpt -= 2;
      for (ii = 0; ii < nSamOpt*nInputs_; ii++)
         samInputs[ii] = samInputs[ii+2*nInputs_];
      delete [] samOutputs;
      delete [] samStates;
      delete sampler;
    }

    KRI_XDists = XDists;
    KRI_nSamples = nSamples_;
    KRI_nInputs = nInputs_;
    KRI_pOrder = pOrder_;
    KRI_X = XNormalized_;
    KRI_Y = YNormalized_;
    KRI_Ytmp = new double[Cleng];
    KRI_OptThetas = Thetas_;
    maxfun = 5000;
    KRI_OptY = PSUADE_UNDEFINED;
    KRI_outputLevel = outputLevel_;
    optThetas = new double[nInputs_];
    SMatrix = new double[nSamples_*nSamples_];
    FMatrix = new double[nSamples_*nBasis];
    FMatTmp = new double[nSamples_*nBasis];
    MMatrix = new double[nBasis*nBasis];
    checkAllocate(MMatrix, "MMatrix in Kriging::train");
    KRI_SMatrix = SMatrix;
    KRI_FMatrix = FMatrix;
    KRI_FMatTmp = FMatTmp;
    KRI_MMatrix = MMatrix;
    double *TValSave = new double[nSamOpt*nInputs_];
    checkAllocate(TValSave, "TValSave in Kriging::train");
    double dmean, dstd;
    int    stopFlag = 0;

    for (kk = 0; kk < nSamOpt; kk++)
    {
      fp = fopen("ps_print", "r");
      if (fp != NULL)
      {
        fscanf(fp, "%d", &outputLevel_);
        fclose(fp);
        fp = NULL;
        if (outputLevel_ > 0 && outputLevel_ <= 5)
          printf("Kriging: output level set to %d.\n", outputLevel_);
        else
        {
          outputLevel_ = 2;
          printf("Kriging: output level set to %d.\n", outputLevel_);
        }
        KRI_outputLevel = outputLevel_;
      }
      if (outputLevel_ >= 1) 
        printf("Kriging multi-start optimization: start = %d (%d)\n",
               kk+1, nSamOpt);
      KRI_iter = 0;
      for (ii = 0; ii < nInputs_; ii++) 
        TValues[ii] = samInputs[kk*nInputs_+ii];
      if (outputLevel_ >= 1) 
      {
        for (ii = 0; ii < nInputs_; ii++) 
          printf("Kriging: Input %4d initial length scale = %e\n",
                 ii+1,TValues[ii]);
      }
      KRI_noProgressCnt = 0;
#ifdef HAVE_BOBYQA
//       pLevel = 8888;
//       kbobyqa_(&nInputs_,&nPts,TValues,TLowers,TUppers,&rhobeg,&rhoend,
//                &pLevel, &maxfun, work);
#endif
#ifdef HAVE_NEWUOA
      pLevel = 6666;
      newuoa_(&nInputs_, &nPts, TValues, &rhobeg, &rhoend, &pLevel,
              &maxfun, work);
#endif
      if (outputLevel_ >= 1) 
      {
        printf("Kriging multi-start optimization: iteration = %d (%d)\n",
               kk+1, nSamOpt);
        for (ii = 0; ii < nInputs_; ii++) 
          printf("Kriging: Input %4d final length scale = %e\n",
                 ii+1,TValues[ii]);
        printf("Kriging final objective value = %e (ref = %e)\n", 
               KRI_currY, rhoend);
      }
      if (KRI_OptY < optY)
      {
        optY = KRI_OptY;
        for (ii = 0; ii < nInputs_; ii++)
          optThetas[ii] = KRI_OptThetas[ii];
        if (optY < rhoend)
        {
          if (outputLevel_ > 2) 
            printf("Kriging INFO: termination (sufficiently accurate)\n");
          break;
        }
      }
      fp = fopen("psuade_stop", "r");
      if (fp != NULL)
      {
        fclose(fp);
        fp = NULL;
        printf("Kriging: psuade_stop file found - terminating ....\n");
        break;
      }
      fp = fopen("ps_rs_expert", "r");
      if (fp != NULL)
      {
        fclose(fp);
        fp = NULL;
        printf("Kriging: turn on rs_expert mode.\n");
        psRSExpertMode_ = 1;
      }
      for (ii = 0; ii < nInputs_; ii++) 
        TValSave[kk*nInputs_+ii] = TValues[ii];
      if (kk >= 4)
      {
        for (ii = 0; ii < nInputs_; ii++)
        {
          dmean = 0.0;
          for (jj = 1; jj < kk+1; jj++) dmean += TValSave[jj*nInputs_+ii];
          dmean /= (double) kk; 
          dstd = 0.0;
          for (jj = 1; jj < kk+1; jj++)
            dstd += pow(TValSave[jj*nInputs_+ii]-dmean,2.0);
          dstd = sqrt(dstd/(double) (kk-1));
          if (dstd/dmean < 0.01) stopFlag++;
          if (outputLevel_ > 1) 
            printf("Opt convergence check: input %d: mean=%e, std=%e\n",
                   ii+1,dmean,dstd);
        }
        if (stopFlag == nInputs_) 
        {
          if (outputLevel_ >= 1) 
          {
            printf("Kriging INFO: same optimum after %d iterations.\n",
                   kk+1);
            printf("              Stop further processing.\n");
          }
          break;
        }  
      }  
    }
    for (ii = 0; ii < nInputs_; ii++) Thetas_[ii] = optThetas[ii];
    if (outputLevel_ >= 1) 
    {
      for (ii = 0; ii < nInputs_; ii++) 
        printf("Kriging: Input %4d optimal length scale = %e\n",
               ii+1,Thetas_[ii]);
      printf("Kriging: optimal objective value = %e\n",KRI_OptY);
    }
    delete [] work;
    delete [] TValues;
    delete [] TValSave;
    delete [] KRI_Ytmp;
    delete [] samInputs;
    delete [] optThetas;
    KRI_OptThetas = NULL;
    KRI_XDists = NULL;
    KRI_X = NULL;
    KRI_Y = NULL;
    KRI_Ytmp = NULL;
    if (SMatrix != NULL) delete [] SMatrix;
    if (FMatrix != NULL) delete [] FMatrix;
    if (FMatTmp != NULL) delete [] FMatTmp;
    if (MMatrix != NULL) delete [] MMatrix;
    KRI_SMatrix = NULL;
    KRI_FMatrix = NULL;
    KRI_FMatTmp = NULL;
    KRI_MMatrix = NULL;
    if (outputLevel_ > 0)
    {
      if (fastMode_ == 2)
           printf("Kriging training (2) ends.\n");
      else printf("Kriging training (3) ends.\n");
    }
    delete [] TUppers;
    delete [] TLowers;
  }

  char uplo='L';
  char trans='N';
  char diag='N';
  count = 0;
  Rmatrix_ = new double[nSamples_*nSamples_];
  checkAllocate(Rmatrix_, "Rmatrix_ in Kriging::train");
  for (jj = 0; jj < nSamples_; jj++)
  {
    Rmatrix_[jj*nSamples_+jj] = 1.0 + 1e-15 * nSamples_;
    if (KRI_nugget != 0.0)  Rmatrix_[jj*nSamples_+jj] += KRI_nugget;
    else if (nInputs_ == 2) Rmatrix_[jj*nSamples_+jj] += 0.01;

    if (dataStdDevs_ != NULL) 
      Rmatrix_[jj*nSamples_+jj] += pow(dataStdDevs_[jj]/YStd_,KRI_Exponent);
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
        dist += pow(XDists[count*nInputs_+ii]/Thetas_[ii],KRI_Exponent);
      dist = exp(-dist);
      if (dist < 1.0e-16) dist = 0.0;
      Rmatrix_[jj*nSamples_+kk] = dist;
      Rmatrix_[kk*nSamples_+jj] = dist;
      count++;
    }
  }
  if (outputLevel_ > 4)
  {
    printf("Kriging: covariance matrix for Cholesky decompositon.\n");
    for (jj = 0; jj < nSamples_; jj++)
    {
      for (kk = 0; kk < nSamples_; kk++)
        printf("%e ", Rmatrix_[jj*nSamples_+kk]);
      printf("\n");
    }
  }
  dpotrf_(&uplo, &nSamples_, Rmatrix_, &nSamples_, &status);
  if (status != 0) 
  {
    printf("Kriging ERROR: Cholesky decomposition not successful.\n");
    exit(1);
  }
  int    inc=1;
  double *CY = new double[nSamples_];
  checkAllocate(CY, "CY in Kriging::train");
  for (jj = 0; jj < nSamples_; jj++) CY[jj] = YNormalized_[jj];
  dtrsv_(&uplo,&trans,&diag,&nSamples_,Rmatrix_,&nSamples_,CY, &inc);
  CFmatrix = new double[nBasis*nSamples_];
  RFmatrix = new double[nBasis*nSamples_];
  checkAllocate(RFmatrix, "RFmatrix in Kriging::train");
  for (jj = 0; jj < nSamples_; jj++)
  {
    CFmatrix[jj] = RFmatrix[jj] = 1.0;
    if (pOrder_ == 1)
    {
      for (kk = 1; kk < nInputs_+1; kk++)
      {
        ddata = XNormalized_[jj*nInputs_+kk-1];
        CFmatrix[kk*nSamples_+jj] = RFmatrix[kk*nSamples_+jj] = ddata;
      }
    }
  }
  for (ii = 0; ii < nBasis; ii++)
  {
    dtrsv_(&uplo,&trans,&diag,&nSamples_,Rmatrix_,&nSamples_,
           &CFmatrix[ii*nSamples_], &inc);
  } 
  V1_ = new double[nBasis];
  checkAllocate(V1_, "V1_ in Kriging::train");
  for (jj = 0; jj < nBasis; jj++)
  {
    ddata = 0.0;
    for (ii = 0; ii < nSamples_; ii++)
      ddata += CFmatrix[ii+jj*nSamples_]*CY[ii];
    V1_[jj] = ddata;
  }
  Mmatrix_ = new double[nBasis*nBasis];
  checkAllocate(Mmatrix_, "Mmatrix_ in Kriging::train");
  for (jj = 0; jj < nBasis; jj++)
  {
    for (kk = jj; kk < nBasis; kk++)
    {
      ddata = 0.0;
      for (ii = 0; ii < nSamples_; ii++)
        ddata += CFmatrix[ii+jj*nSamples_]*
                 CFmatrix[ii+kk*nSamples_];
      Mmatrix_[jj+kk*nBasis] = Mmatrix_[kk+jj*nBasis] = ddata;
    }
  }
  dpotrf_(&uplo, &nBasis, Mmatrix_, &nBasis, &status);
  if (status != 0) 
  {
    printf("Kriging ERROR: Cholesky decomposition not successful.\n");
    exit(1);
  }
  dpotrs_(&uplo, &nBasis, &inc, Mmatrix_, &nBasis, V1_, &nBasis, &status);
  if (status != 0) 
  {
    printf("Kriging ERROR: LU solve not successful.\n");
    exit(1);
  }
  for (ii = 0; ii < nBasis; ii++)
  {
    for (jj = 0; jj < nSamples_; jj++) 
      CY[jj] -= CFmatrix[jj+ii*nSamples_] * V1_[ii];
  } 
  KrigingVariance_ = 0.0;
  for (jj = 0; jj < nSamples_; jj++) KrigingVariance_ += CY[jj] * CY[jj];
  KrigingVariance_ /= (double) nSamples_;
  KrigingVariance_ *= YStd_ * YStd_;
  if (outputLevel_ > 0) printf("Kriging variance = %e\n",KrigingVariance_);
  V2_ = new double[nSamples_];
  checkAllocate(V2_, "V2_ in Kriging::train");
  for (ii = 0; ii < nSamples_; ii++) V2_[ii] = CY[ii];
  trans = 'T';
  dtrsv_(&uplo,&trans,&diag,&nSamples_,Rmatrix_,&nSamples_,V2_, &inc);
  delete [] CFmatrix;
  delete [] RFmatrix;
  delete [] CY;
  delete [] XDists;
  //**  END

  //**  store results
  if (noReuse_ == 1 || psRSCodeGen_ == 0) return 0;
  FILE *fp = fopen("psuade_rs.info", "w");
  if (fp == NULL) return 0.0;
#if 1
  fprintf(fp,"This file contains information to re-construct Kriging\n");
  fprintf(fp,"response surface offline. Follow the steps below:\n");
  fprintf(fp,"1. Search for the keywords 'SPLIT HERE' in this file.\n");
  fprintf(fp,"2. Store the lines below keywords into main.c\n");
  fprintf(fp,"3. Add the to-be-evaluated points inside main() in main.c\n");
  fprintf(fp,"   (or, embed all except main() in your own program).\n");
  fprintf(fp,"4. Compile main.c (cc -o main main.c -lm) and run\n");
  fprintf(fp,"   (Do not change the psuade_rs.info file.\n");
  fprintf(fp,"Note: if std dev is desired, uncomment the corresponding\n");
  fprintf(fp,"      section and compile with: cc main.c -llapack -lm\n");
  fprintf(fp,"\n");
  fprintf(fp,"PSUADE_BEGIN\n");
  fprintf(fp, "%d\n", nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    fprintf(fp, "%24.16e %24.16e %24.16e\n", XMeans_[ii], XStds_[ii],
            Thetas_[ii]);
  fprintf(fp, "%24.16e %24.16e\n", YMean_, YStd_);
  fprintf(fp, "%d %d\n", nSamples_, nInputs_);
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, "%24.16e ", XNormalized_[jj*nInputs_+ii]);
    fprintf(fp, "\n");
  }
  fprintf(fp,"%d\n", pOrder_);
  fprintf(fp, "%24.16e \n", V1_[0]);
  if (pOrder_ == 1)
  {
    for (ii = 1; ii <= nInputs_; ii++) fprintf(fp, "%24.16e \n", V1_[ii]);
  }
  fprintf(fp, "%d\n", nSamples_);
  for (jj = 0; jj < nSamples_; jj++) fprintf(fp, "%24.16e \n", V2_[jj]);
  fprintf(fp,"%d\n", nSamples_);
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nSamples_; ii++)
      fprintf(fp, "%24.16e ", Rmatrix_[jj+ii*nSamples_]);
    fprintf(fp, "\n");
  }
  fprintf(fp,"%d\n", nBasis);
  for (jj = 0; jj < nBasis; jj++)
  {
    for (ii = 0; ii < nBasis; ii++)
      fprintf(fp, "%24.16e ", Mmatrix_[jj+ii*nBasis]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "%24.16e;\n", KrigingVariance_);
  fprintf(fp,"====================== SPLIT HERE =====================\n");
  fprintf(fp,"/* *******************************************/\n");
  fprintf(fp,"/* Kriging interpolator from PSUADE. */\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"#include <math.h>\n");
  fprintf(fp,"#include <stdlib.h>\n");
  fprintf(fp,"#include <stdio.h>\n");
  fprintf(fp,"int initialize();\n");
  fprintf(fp,"int finalize();\n");
  fprintf(fp,"int interpolate(int,double*,double*,double*);\n");
  fprintf(fp,"int CholSolve(int, double *, double *);\n");
  fprintf(fp,"main(int argc, char **argv) {\n");
  fprintf(fp,"  int    i, iOne=1, nInps;\n");
  fprintf(fp,"  double X[%d], Y, Std;\n",nInputs_);
  fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
  fprintf(fp,"  if (argc != 3) {\n");
  fprintf(fp,"     printf(\"ERROR: not enough argument.\\n\");\n");
  fprintf(fp,"     exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fIn = fopen(argv[1], \"r\");\n");
  fprintf(fp,"  if (fIn == NULL) {\n");
  fprintf(fp,"     printf(\"ERROR: cannot open input file.\\n\");\n");
  fprintf(fp,"     exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fscanf(fIn, \"%%d\", &nInps);\n");
  fprintf(fp,"  if (nInps != %d) {\n", nInputs_);
  fprintf(fp,"    printf(\"ERROR - wrong nInputs.\\n\");\n");
  fprintf(fp,"    exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  for (i=0; i<%d; i++) fscanf(fIn, \"%%lg\", &X[i]);\n",
          nInputs_);
  fprintf(fp,"  fclose(fIn);\n");
  fprintf(fp,"  initialize();\n");
  fprintf(fp,"  interpolate(iOne, X, &Y, &Std);\n");
  fprintf(fp,"  printf(\"Y = %%e (stdev = %%e)\\n\", Y, Std);\n");
  fprintf(fp,"  finalize();\n");
  fprintf(fp,"  fOut = fopen(argv[2], \"w\");\n");
  fprintf(fp,"  if (fOut == NULL) {\n");
  fprintf(fp,"     printf(\"ERROR: cannot open output file.\\n\");\n");
  fprintf(fp,"     exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fprintf(fOut,\" %%e\\n\", Y);\n");
  fprintf(fp,"  fclose(fOut);\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"/* Regression interpolation function         */\n");
  fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
  fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
  fprintf(fp," * ... */\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int    pOrder, nInps, nSamples, nBasis;\n");
  fprintf(fp,"double *XMeans=NULL,*XStds=NULL,YMean,YStd,*Thetas=NULL;\n");
  fprintf(fp,"double *V1=NULL,*V2=NULL,*Rmat=NULL,*Mmat=NULL,variance;\n");
  fprintf(fp,"double *XNorm=NULL;\n");
  fprintf(fp,"int interpolate(int npts,double *X,double *Y,\n");
  fprintf(fp,"                double *YStds) {\n");
  fprintf(fp,"  int    ss, ii, jj, iOne=1, status;\n");
  fprintf(fp,"  double dd, *x, *w1, *w2, *wx, sig;\n");
  fprintf(fp,"  char uplo='L';\n");
  fprintf(fp,"  w1 = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  w2 = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  wx = (double *) malloc(2*nBasis*sizeof(double));\n");
  fprintf(fp,"  x  = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  for (ss = 0; ss < npts; ss++) {\n");
  fprintf(fp,"    for (ii = 0; ii < nInps; ii++)\n");
  fprintf(fp,"      x[ii] = (X[ss*nInps+ii]-XMeans[ii])/XStds[ii];\n");
  fprintf(fp,"    for (jj = 0; jj < nSamples; jj++) {\n");
  fprintf(fp,"      sig = 0.0;\n");
  fprintf(fp,"      for (ii = 0; ii < nInps; ii++) {\n");
  fprintf(fp,"        dd = XNorm[jj*nInps+ii] - x[ii];\n");
  fprintf(fp,"        sig += dd * dd / Thetas[ii] / Thetas[ii];\n");
  fprintf(fp,"      }\n");
  fprintf(fp,"      dd = exp(-sig);\n");
  fprintf(fp,"      if (dd < 1e-16) dd = 0.0;\n");
  fprintf(fp,"      w1[jj] = w2[jj] = dd;\n");
  fprintf(fp,"    }\n");
  fprintf(fp,"    dd = V1[0];\n");
  fprintf(fp,"    if (pOrder == 1)\n");
  fprintf(fp,"      for (jj = 1; jj <= nInps; jj++) \n");
  fprintf(fp,"        dd += x[jj-1] * V1[jj];\n");
  fprintf(fp,"    for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"      dd += w2[jj] * V2[jj];\n");
  fprintf(fp,"    Y[ss] = dd * YStd + YMean;\n");
  fprintf(fp,"    /* ==== if need to compute std dev. =====*/\n");
  fprintf(fp,"    CholSolve(nSamples,Rmat,w2);\n");
  fprintf(fp,"    /*dpotrs_(&uplo,&nSamples,&iOne,Rmat,&nSamples,w2,\n");
  fprintf(fp,"            &nSamples,&status);*/\n");
  fprintf(fp,"    sig = 0.0;\n");
  fprintf(fp,"    for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"      sig += w1[jj] * w2[jj];\n");
  fprintf(fp,"    dd = 0.0;\n");
  fprintf(fp,"    for (jj = 0; jj < nSamples; jj++) dd += w2[jj];\n");
  fprintf(fp,"    wx[0] = wx[nBasis] = dd - 1;\n");
  if (nBasis > 1)
  {
    fprintf(fp,"    for (ii = 1; ii < nBasis; ii++) {\n");
    fprintf(fp,"      dd = 0.0;\n");
    fprintf(fp,"      for (jj = 0; jj < nSamples; jj++)\n");
    fprintf(fp,"        dd += x[ii-1] * w2[jj];\n");
    fprintf(fp,"      dd -= x[ii-1];\n");
    fprintf(fp,"      wx[ii] = wx[nBasis+ii] = dd;\n");
    fprintf(fp,"    }\n");
  }
  fprintf(fp,"    CholSolve(nBasis,Mmat,wx);\n");
  fprintf(fp,"    /* dpotrs_(&uplo,&nBasis,&iOne,Mmat,&nSamples,\n");
  fprintf(fp,"            wx, &nBasis, &status); */\n");
  fprintf(fp,"    dd = 0.0;\n");
  fprintf(fp,"    for (ii = 0; ii < nBasis; ii++) \n");
  fprintf(fp,"      dd += wx[ii] * wx[ii+nBasis];\n");
  fprintf(fp,"    dd = variance*(1+dd-sig);\n");
  fprintf(fp,"    if (dd < 0) dd = 0.0;\n");
  fprintf(fp,"    YStds[ss] = sqrt(dd) * YStd;\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  free(wx);\n");
  fprintf(fp,"  free(w1);\n");
  fprintf(fp,"  free(w2);\n");
  fprintf(fp,"  free(x);\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int CholSolve(int n, double *A, double *X) {\n");
  fprintf(fp,"  int    ii, jj, kk;\n");
  fprintf(fp,"  double *Y, ddata;\n");
  fprintf(fp,"  Y = (double *) malloc(n*sizeof(double));\n");
  fprintf(fp,"  Y[0] = X[0] / A[0];\n");
  fprintf(fp,"  for (ii = 1; ii < n; ii++) {\n");
  fprintf(fp,"     ddata = 0.0;\n");
  fprintf(fp,"     for (jj = 0; jj < ii; jj++) \n");
  fprintf(fp,"       ddata += A[jj*n+ii] * Y[jj];\n");
  fprintf(fp,"     Y[ii] = (X[ii] - ddata) / A[ii*n+ii];\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  X[n-1] = Y[n-1] / A[(n-1)*n+n-1];\n");
  fprintf(fp,"  for (ii = n-2; ii >= 0; ii--) {\n");
  fprintf(fp,"     ddata = 0.0;\n");
  fprintf(fp,"     for (jj = ii+1; jj < n; jj++) \n");
  fprintf(fp,"       ddata += A[ii*n+jj] * X[jj];\n");
  fprintf(fp,"     X[ii] = (Y[ii] - ddata) / A[ii*n+ii];\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  free(Y);\n");
  fprintf(fp,"  return 0;\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int initialize() {\n");
  fprintf(fp,"  int    ii, jj;\n");
  fprintf(fp,"  double ddata;\n");
  fprintf(fp,"  char   line[1001], word[1001];\n");
  fprintf(fp,"  FILE *fp = fopen(\"psuade_rs.info\", \"r\");\n");
  fprintf(fp,"  if (fp == NULL){\n");
  fprintf(fp,"     printf(\"Data file (psuade_rs.info) not found.\\n\");\n");
  fprintf(fp,"     exit(1);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  while (1) {\n");
  fprintf(fp,"    fgets(line, 500, fp);\n");
  fprintf(fp,"    sscanf(line, \"%%s\",word);\n");
  fprintf(fp,"    if (!strcmp(word, \"PSUADE_BEGIN\")) break;\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fscanf(fp, \"%%d\", &nInps);\n");
  fprintf(fp,"  XMeans = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  XStds  = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  Thetas = (double *) malloc(nInps*sizeof(double));\n");
  fprintf(fp,"  for (ii = 0; ii < nInps; ii++) {\n");
  fprintf(fp,"     fscanf(fp, \"%%lg\", &ddata);\n");
  fprintf(fp,"     XMeans[ii] = ddata;\n");
  fprintf(fp,"     fscanf(fp, \"%%lg\", &ddata);\n");
  fprintf(fp,"     XStds[ii] = ddata;\n");
  fprintf(fp,"     fscanf(fp, \"%%lg\", &ddata);\n");
  fprintf(fp,"     Thetas[ii] = ddata;\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fscanf(fp, \"%%lg %%lg\", &YMean, &YStd);\n");
  fprintf(fp,"  fscanf(fp, \"%%d %%d\", &nSamples, &nInps);\n");
  fprintf(fp,"  XNorm = (double *) malloc(nSamples*nInps*sizeof(double));\n");
  fprintf(fp,"  for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"    for (ii = 0; ii < nInps; ii++) \n");
  fprintf(fp,"      fscanf(fp, \"%%lg\", &XNorm[jj*nInps+ii]);\n");
  fprintf(fp,"  fscanf(fp, \"%%d\", &pOrder);\n");
  fprintf(fp,"  V1 = (double *) malloc((nInps+1)*sizeof(double));\n");
  fprintf(fp,"  fscanf(fp, \"%%lg\", &V1[0]);\n");
  fprintf(fp,"  if (pOrder == 1) { \n");
  fprintf(fp,"    for (ii = 1; ii <= nInps; ii++)\n");
  fprintf(fp,"        fscanf(fp, \"%%24.16e\", &V1[ii]);\n");
  fprintf(fp,"  }\n");
  fprintf(fp,"  fscanf(fp, \"%%d\", &nSamples);\n");
  fprintf(fp,"  V2 = (double *) malloc(nSamples*sizeof(double));\n");
  fprintf(fp,"  for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"    fscanf(fp, \"%%lg\", &V2[ii]);\n");
  fprintf(fp,"  fscanf(fp, \"%%d\", &nSamples);\n");
  fprintf(fp,"  Rmat=(double*) malloc(nSamples*nSamples*sizeof(double));\n");
  fprintf(fp,"  for (jj = 0; jj < nSamples; jj++) \n");
  fprintf(fp,"    for (ii = 0; ii < nSamples; ii++)\n");
  fprintf(fp,"      fscanf(fp, \"%%lg \", &Rmat[jj+ii*nSamples]);\n");
  fprintf(fp,"  fscanf(fp, \"%%d\", &nBasis);\n");
  fprintf(fp,"  Mmat=(double*) malloc(nBasis*nBasis*sizeof(double));\n");
  fprintf(fp,"  for (jj = 0; jj < nBasis; jj++) \n");
  fprintf(fp,"    for (ii = 0; ii < nBasis; ii++)\n");
  fprintf(fp,"      fscanf(fp, \"%%lg \", &Mmat[jj+ii*nBasis]);\n");
  fprintf(fp,"  fscanf(fp, \"%%lg\", &variance);\n");
  fprintf(fp,"  fclose(fp);\n");
  fprintf(fp,"  return 0;\n");
  fprintf(fp,"}\n");
  fprintf(fp,"/* ==========================================*/\n");
  fprintf(fp,"int finalize() {\n");
  fprintf(fp,"  if (V1 != NULL) free(V1);\n");
  fprintf(fp,"  if (V2 != NULL) free(V2);\n");
  fprintf(fp,"  if (Thetas != NULL) free(Thetas);\n");
  fprintf(fp,"  if (XMeans != NULL) free(XMeans);\n");
  fprintf(fp,"  if (XStds != NULL) free(XStds);\n");
  fprintf(fp,"  if (Mmat != NULL) free(Mmat);\n");
  fprintf(fp,"  if (Rmat != NULL) free(Rmat);\n");
  fprintf(fp,"  if (XNorm != NULL) free(XNorm);\n");
  fprintf(fp,"}\n");
#else
  fprintf(fp, "%% Kriging information after construction\n");
  fprintf(fp, "%% To interpolate with Kriging, do the following: \n");
  fprintf(fp, "%% 1. Take this file and rename it to be some .m file.\n");
  fprintf(fp, "%% 2. The, put the to-be-evaluated points in X\n");
  fprintf(fp, "%%    (X is a Nxm matrix of size N and m inputs)\n");
  fprintf(fp, "%% 3. Run this .m file. The results will be in Y and\n");
  fprintf(fp, "%%    the corresponding standard deviation is in Sd.\n");
  fprintf(fp, "%%\n");
  fprintf(fp, "%%In summary:\n");
  fprintf(fp, "%% Inputs : X\n");
  fprintf(fp, "%% Outputs: Y, Sd\n");
  fprintf(fp, "%%\n");
  fprintf(fp, "leng   = size(X,1);\n");
  fprintf(fp, "Y      = zeros(leng,1);\n");
  fprintf(fp, "Sd     = zeros(leng,1);\n\n");
  fprintf(fp, "nSamp  = %d;\n", nSamples_);
  fprintf(fp, "order  = %d;\n", pOrder_);
  fprintf(fp, "nbasis = %d;\n", nBasis);
  fprintf(fp, "NormalizedX = [\n");
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
      fprintf(fp, "%24.16e ", XNormalized_[jj*nInputs_+ii]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "];\n");
  fprintf(fp, "XMeans = [\n");
  for (ii = 0; ii < nInputs_; ii++) fprintf(fp,"%24.16e \n",XMeans_[ii]);
  fprintf(fp, "];\n");
  fprintf(fp, "XStds = [\n");
  for (ii = 0; ii < nInputs_; ii++) fprintf(fp,"%24.16e \n",XStds_[ii]);
  fprintf(fp, "];\n");
  fprintf(fp, "YMean = %24.16e ;\n", YMean_);
  fprintf(fp, "YStd  = %24.16e ;\n", YStd_);
  fprintf(fp, "Thetas = [\n");
  for (ii = 0; ii < nInputs_; ii++) fprintf(fp,"%24.16e \n",Thetas_[ii]);
  fprintf(fp, "];\n");
  fprintf(fp, "V1 = [\n");
  fprintf(fp, "%24.16e \n", V1_[0]);
  if (pOrder_ == 1)
  {
    for (ii = 1; ii <= nInputs_; ii++) fprintf(fp, "%24.16e \n", V1_[ii]);
  }
  fprintf(fp, "];\n");
  fprintf(fp, "V2 = [\n");
  for (jj = 0; jj < nSamples_; jj++)
    fprintf(fp, "%24.16e \n", V2_[jj]);
  fprintf(fp, "];\n");
  fprintf(fp, "Rmat = [\n");
  for (jj = 0; jj < nSamples_; jj++)
  {
    for (ii = 0; ii < nSamples_; ii++)
      fprintf(fp, "%24.16e ", Rmatrix_[jj+ii*nSamples_]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "];\n");
  fprintf(fp, "Mmat = [\n");
  for (jj = 0; jj < nBasis; jj++)
  {
    for (ii = 0; ii < nBasis; ii++)
      fprintf(fp, "%24.16e ", Mmatrix_[jj+ii*nBasis]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "];\n");
  fprintf(fp, "U = zeros(nbasis,1);\n");
  fprintf(fp, "KVar = %24.16e;\n", KrigingVariance_);
  fprintf(fp, "RR = zeros(size(NormalizedX,1),1);\n\n");
  fprintf(fp, "%% Processing \n");
  fprintf(fp, "for ii = 1 : leng\n");
  fprintf(fp, "   XN = (X(ii,:)' - XMeans) ./ XStds;\n");
  fprintf(fp, "   for jj = 1 : nSamp\n");
  fprintf(fp, "      XT = NormalizedX(jj,:)' - XN;\n");
  fprintf(fp, "      dd = sum((XT .* XT) ./ (Thetas .* Thetas));\n");
  fprintf(fp, "      RR(jj) = exp(-dd);\n");
  fprintf(fp, "   end;\n");
  fprintf(fp, "   YT = 1.0 * V1(1);\n");
  fprintf(fp, "   if (order == 1);\n");
  fprintf(fp, "      YT = YT + sum(V1(2:nbasis) .* XN);\n");
  fprintf(fp, "   end;\n");
  fprintf(fp, "   YT = YT + sum(V2 .* RR);\n");
  fprintf(fp, "   Y(ii) = YT * YStd + YMean;\n");
  fprintf(fp, "   SS = Rmat' \\ (Rmat \\ RR);\n");
  fprintf(fp, "   SR = SS' * RR;\n");
  fprintf(fp, "   U(1) = sum(RR) - 1;\n");
  fprintf(fp, "   for jj = 2 : nbasis\n");
  fprintf(fp, "      U(jj) = XN(jj-1) * (sum(RR) - 1); \n");
  fprintf(fp, "   end;\n");
  fprintf(fp, "   V = Mmat' \\ (Mmat \\ U);\n");
  fprintf(fp, "   Sd(ii) = KVar * (1 + U' * V - SR);\n");
  fprintf(fp, "   if (Sd(ii) < 0)\n");
  fprintf(fp, "      Sd(ii) = 0;\n");
  fprintf(fp, "   end\n");
  fprintf(fp, "end;\n");
  fprintf(fp, "display('The evaluated outputs and standard deviations are :\n");
  fprintf(fp, "[(1:leng)' Y Sd]\n");
  fprintf(fp, "%% an example for 2D\n");
  fprintf(fp, "%% X1 = 67:3:540; \n");
  fprintf(fp, "%% X2 = 0.4887:(1.0864-0.4887)/101:1.0864; \n");
  fprintf(fp, "%% [XX1,XX2] = meshgrid(X1,X2);\n");
  fprintf(fp, "%% MM = XX1;\n");
  fprintf(fp, "%% SD = XX1;\n");
  fprintf(fp, "%% for ii = 1 : size(XX1,2)\n");
  fprintf(fp, "%%    for kk = 1 : size(XX1,1)\n");
  fprintf(fp, "%%       XI = [XX1(kk,ii)\n");
  fprintf(fp, "%%             XX2(kk,ii)];\n");
  fprintf(fp, "%%       XN = (XI - XMeans) ./ XStds;\n");
  fprintf(fp, "%%       for jj = 1 : 31\n");
  fprintf(fp, "%%          XX = NormalizedX(jj,:)' - XN;\n");
  fprintf(fp, "%%          dd = sum((XX .* XX) ./ (Thetas .* Thetas));\n");
  fprintf(fp, "%%          RR(jj) = exp(-dd);\n");
  fprintf(fp, "%%       end;\n");
  fprintf(fp, "%%       YT = 1.0 * V1(1);\n");
  fprintf(fp, "%%       if (order == 1);\n");
  fprintf(fp, "%%          YT = YT + sum(V1(2:nbasis) .* XN);\n");
  fprintf(fp, "%%       end;\n");
  fprintf(fp, "%%       YT = YT + sum(V2 .* RR);\n");
  fprintf(fp, "%%       MM(kk,ii) = YT * YStd + YMean;\n");
  fprintf(fp, "%%       SS = Rmat' \\ (Rmat \\ RR);\n");
  fprintf(fp, "%%       SR = SS' * RR;\n");
  fprintf(fp, "%%       U(1) = sum(RR) - 1;\n");
  fprintf(fp, "%%       for jj = 2 : nbasis\n");
  fprintf(fp, "%%          U(jj) = XN(jj-1) * (sum(RR) - 1); \n");
  fprintf(fp, "%%       end;\n");
  fprintf(fp, "%%       V = Mmat' \\ (Mmat \\ U);\n");
  fprintf(fp, "%%       SD(kk,ii) = KVar * (1 + U' * V - SR);\n");
  fprintf(fp, "%%       if (SD(kk,ii) < 0)\n");
  fprintf(fp, "%%          SD(kk,ii) = 0;\n");
  fprintf(fp, "%%       end\n");
  fprintf(fp, "%%    end;\n");
  fprintf(fp, "%% end;\n");
  fprintf(fp, "%% mesh(MM)\n");
#endif
  fclose(fp);
  printf("Kriging response surface data file is in psuade_rs.info\n");
  fp = fopen("psuade_rs.py", "w");
  if (fp == NULL)
  {
    printf("ERROR: Cannot open file psuade_rs.py.\n");
    return 0;
  }
  fwriteRSPythonHeader(fp);
  fprintf(fp,"#==================================================\n");
  fprintf(fp,"# Kriging interpolation\n");
  fprintf(fp,"#==================================================\n");
  fwriteRSPythonCommon(fp);
  fprintf(fp, "nSamples = %d\n", nSamples_);
  fprintf(fp, "nInputs = %d\n", nInputs_);
  fprintf(fp, "nBasis = %d\n", nBasis);
  fprintf(fp, "KrigVariance = %24.16e\n", KrigingVariance_);
  fprintf(fp, "XParams = [\n");
  for (ii = 0; ii < nInputs_; ii++)
    fprintf(fp, "[%24.16e , %24.16e , %24.16e ],\n", XMeans_[ii], XStds_[ii],
            Thetas_[ii]);
  fprintf(fp, "]\n");
  fprintf(fp, "Ymean = %24.16e\n", YMean_);
  fprintf(fp, "Ystd  = %24.16e\n", YStd_);
  fprintf(fp, "SamInputs = [\n");
  for (jj = 0; jj < nSamples_; jj++)
  {
    fprintf(fp, "[ %24.16e ", XNormalized_[jj*nInputs_]);
    for (ii = 1; ii < nInputs_; ii++)
      fprintf(fp, ", %24.16e ", XNormalized_[jj*nInputs_+ii]);
    fprintf(fp, "], \n");
  }
  fprintf(fp, "]\n");
  fprintf(fp, "V1 = [\n");
  fprintf(fp, " %24.16e , \n", V1_[0]);
  if (pOrder_ == 1)
  {
    for (ii = 1; ii <= nInputs_; ii++)
      fprintf(fp, " %24.16e ,  \n", V1_[ii]);
  }
  fprintf(fp, "]\n");
  fprintf(fp, "V2 = [\n");
  for (jj = 0; jj < nSamples_; jj++)
    fprintf(fp, " %24.16e , \n", V2_[jj]);
  fprintf(fp, "]\n");
  fprintf(fp, "Rmat = [\n");
  for (jj = 0; jj < nSamples_; jj++)
  {
    fprintf(fp, " [ %24.16e ", Rmatrix_[jj]);
    for (ii = 1; ii < nSamples_; ii++)
      fprintf(fp, " , %24.16e ", Rmatrix_[jj+ii*nSamples_]);
    fprintf(fp, "], \n");
  }
  fprintf(fp, "]\n");
  fprintf(fp, "Mmat = [\n");
  for (jj = 0; jj < nBasis; jj++)
  {
    fprintf(fp, " [ %24.16e ", Mmatrix_[jj]);
    for (ii = 1; ii < nBasis; ii++)
      fprintf(fp, " , %24.16e ", Mmatrix_[jj+ii*nBasis]);
    fprintf(fp, "], \n");
  }
  fprintf(fp, "]\n");
  fprintf(fp,"###################################################\n");
  fprintf(fp,"# solve for A X = B\n");
  fprintf(fp,"#==================================================\n");
  fprintf(fp,"def MatSolve(n, Amat, B): \n");
  fprintf(fp,"  Y = n * [0.0]\n");
  fprintf(fp,"  Y[0] = B[0] / Amat[0][0]\n");
  fprintf(fp,"  for ii in range(n-1): \n");
  fprintf(fp,"    jj = ii + 1\n");
  fprintf(fp,"    ddata = 0.0\n");
  fprintf(fp,"    for kk in range(jj): \n");
  fprintf(fp,"      ddata += Amat[jj][kk] * Y[kk]\n");
  fprintf(fp,"    Y[jj] = (B[jj] - ddata)/Amat[jj][jj]\n");
  fprintf(fp,"  X = n * [0.0]\n");
  fprintf(fp,"  X[n-1] = Y[n-1] / Amat[n-1][n-1]\n");
  fprintf(fp,"  for ii in range(n-1): \n");
  fprintf(fp,"    jj = n - ii - 2\n");
  fprintf(fp,"    ddata = 0.0\n");
  fprintf(fp,"    for kk in range(ii+1): \n");
  fprintf(fp,"      kk2 = n - kk - 1\n");
  fprintf(fp,"      ddata += Amat[kk2][jj] * X[kk2]\n");
  fprintf(fp,"    X[jj] = (Y[jj] - ddata)/Amat[jj][jj]\n");
  fprintf(fp,"  return X\n");
  fprintf(fp,"###################################################\n");
  fprintf(fp,"# Interpolation function  \n");
  fprintf(fp,"# X[0], X[1],   .. X[m-1]   - first point\n");
  fprintf(fp,"# X[m], X[m+1], .. X[2*m-1] - second point\n");
  fprintf(fp,"# ... \n");
  fprintf(fp,"#==================================================\n");
  fprintf(fp,"def interpolate(XX): \n");
  fprintf(fp,"  npts = int(len(XX) / nInputs + 1e-8)\n");
  fprintf(fp,"  Xt = nInputs * [0.0]\n");
  fprintf(fp,"  Ys = 2 * npts * [0.0]\n");
  fprintf(fp,"  W1 = nSamples * [0.0]\n");
  fprintf(fp,"  W2 = nSamples * [0.0]\n");
  fprintf(fp,"  WX = (2 * nBasis) * [0.0]\n");
  fprintf(fp,"  for ss in range(npts) : \n");
  fprintf(fp,"    for ii in range(nInputs) : \n");
  fprintf(fp,
     "      Xt[ii] = (XX[ss*nInputs+ii]-XParams[ii][0])/XParams[ii][1]\n");
  fprintf(fp,"    for jj in range(nSamples) : \n");
  fprintf(fp,"      sig = 0.0\n");
  fprintf(fp,"      for ii in range(nInputs) : \n");
  fprintf(fp,"        dd = SamInputs[jj][ii] - Xt[ii];\n");
  fprintf(fp,"        sig = sig + dd*dd/(XParams[ii][2]*XParams[ii][2])\n");
  fprintf(fp,"      dd = math.exp(-sig)\n");
  fprintf(fp,"      if (dd < 1e-16) : \n");
  fprintf(fp,"        dd = 0.0\n");
  fprintf(fp,"      W1[jj] = dd\n");
  fprintf(fp,"      W2[jj] = dd\n");
  fprintf(fp,"    dd = V1[0]\n");
  if (pOrder_ == 1)
  {
    fprintf(fp,"    for ii in range(nInputs) : \n");
    fprintf(fp,"      dd = dd + Xt[ii] * V1[ii+1]\n");
  }
  fprintf(fp,"    for jj in range(nSamples) : \n");
  fprintf(fp,"      dd = dd + W2[jj] * V2[jj]\n");
  fprintf(fp,"    Ys[2*ss] = dd * Ystd + Ymean\n");
  fprintf(fp,"    W2a = MatSolve(nSamples,Rmat,W2)\n");
  fprintf(fp,"    sig = 0.0\n");
  fprintf(fp,"    for jj in range(nSamples) : \n");
  fprintf(fp,"      sig += W1[jj] * W2a[jj]\n");
  fprintf(fp,"    dd = 0.0;\n");
  fprintf(fp,"    for jj in range(nSamples) : \n");
  fprintf(fp,"      dd = dd + W2a[jj]\n");
  fprintf(fp,"    WX[0] = WX[nBasis] = dd - 1\n");
  if (nBasis > 1)
  {
    fprintf(fp,"    for ii in range(nBasis-1) : \n");
    fprintf(fp,"      dd = 0.0\n");
    fprintf(fp,"      for jj in range(nSamples) : \n");
    fprintf(fp,"        dd = dd + Xt[ii] * W2a[jj]\n");
    fprintf(fp,"      dd = dd - Xt[ii];\n");
    fprintf(fp,"      WX[ii+1] = dd\n");
    fprintf(fp,"      WX[nBasis+ii+1] = dd\n");
  }
  fprintf(fp,"    WXa = MatSolve(nBasis,Mmat,WX)\n");
  fprintf(fp,"    dd = 0.0\n");
  fprintf(fp,"    for ii in range(nBasis) : \n");
  fprintf(fp,"      dd = dd + WXa[ii] * WX[ii]\n");
  fprintf(fp,"    dd = KrigVariance*(1+dd-sig)\n");
  fprintf(fp,"    if (dd < 0) :\n"); 
  fprintf(fp,"      dd = 0.0\n");
  fprintf(fp,"    Ys[2*ss+1] = math.sqrt(dd) * Ystd\n");
  fprintf(fp,"  return Ys\n");
  fprintf(fp,"###################################################\n");
  fprintf(fp,"# main program\n");
  fprintf(fp,"#==================================================\n");
  fprintf(fp,"infileName  = sys.argv[1]\n");
  fprintf(fp,"outfileName = sys.argv[2]\n");
  fprintf(fp,"inputs = getInputData(infileName)\n");
  fprintf(fp,"outputs = interpolate(inputs)\n");
  fprintf(fp,"genOutputFile(outfileName, outputs)\n");
  fprintf(fp,"###################################################\n");
  fclose(fp);
  printf("FILE psuade_rs.py contains the Kriging interpolator.\n");
  return 0.0;
}

// ************************************************************************
// predict 
// ------------------------------------------------------------------------
double Kriging::predict(int length, double *X, double *Y, double *YStds)
{
  int    ii, jj, kk, ll, status, nRows, nBasis, leng,maxLeng=500;
  int    offset;
  double ddata, dist, mean, stdev;
  char   uplo='L';

  if (fastMode_ == 0)
  {
    predict0(length, X, Y, YStds);
    return 0.0;
  }
#if 1
  nRows = nSamples_ + 1;
  if (pOrder_ == 1) nRows = nRows + nInputs_;
  nBasis = nRows - nSamples_;
  if (workArray_ == NULL)
  {
    workArray_ = new double[2*nSamples_*maxLeng+maxLeng];
    workX_ = new double[maxLeng*nInputs_+2*maxLeng*nBasis];
    checkAllocate(workX_, "workX_ in Kriging::predict");
    workLength_ = maxLeng;
  }
  for (ll = 0; ll < length; ll+=maxLeng)
  {
    leng = maxLeng;
    if (ll+leng > length) leng = length - ll;
    for (ii = 0; ii < nInputs_; ii++)
    {
      mean  = XMeans_[ii];
      stdev = XStds_[ii];
      for (kk = 0; kk < leng; kk++)
        workX_[kk*nInputs_+ii] = (X[(ll+kk)*nInputs_+ii]-mean)/stdev;
    }
    for (kk = 0; kk < leng; kk++)
    {
      for (jj = 0; jj < nSamples_; jj++)
      {
        dist = 0.0;
        for (ii = 0; ii < nInputs_; ii++)
        {
          ddata = XNormalized_[jj*nInputs_+ii]-workX_[kk*nInputs_+ii];
          dist += ddata * ddata / (Thetas_[ii] * Thetas_[ii]);
        }
        ddata = exp(-dist);
        if (ddata < 1.0e-16) ddata = 0.0;
        workArray_[kk*nSamples_+jj] = ddata;
        workArray_[nSamples_*leng+kk*nSamples_+jj] = ddata;
      }
    }
    for (kk = 0; kk < leng; kk++)
    {
      ddata = 1.0 * V1_[0];
      if (pOrder_ == 1)
        for (ii = 1; ii <= nInputs_; ii++) 
          ddata += workX_[kk*nInputs_+ii-1]*V1_[ii];
      for (jj = 0; jj < nSamples_; jj++) 
        ddata += workArray_[nSamples_*leng+kk*nSamples_+jj] * V2_[jj];
      Y[ll+kk] = ddata * YStd_ + YMean_;
    }
    if (YStds != NULL)
    {
      dpotrs_(&uplo, &nSamples_, &leng, Rmatrix_, &nSamples_,
              &workArray_[leng*nSamples_], &nSamples_, &status);
      for (kk = 0; kk < leng; kk++)
      {
        ddata = 0.0;
        for (jj = 0; jj < nSamples_; jj++)
          ddata += workArray_[kk*nSamples_+jj] *
                   workArray_[leng*nSamples_+kk*nSamples_+jj];
        workArray_[2*leng*nSamples_+kk] = ddata; 
      }
      for (kk = 0; kk < leng; kk++)
      {
        ddata = 0.0;
        for (jj = 0; jj < nSamples_; jj++)
           ddata += workArray_[leng*nSamples_+kk*nSamples_+jj]; 
        workX_[leng*nInputs_+kk*nBasis] = ddata - 1.0; 
        workX_[leng*nInputs_+leng*nBasis+kk*nBasis] = ddata - 1.0; 
        for (ii = 1; ii < nBasis; ii++)
        {
          ddata = 0.0;
          for (jj = 0; jj < nSamples_; jj++)
            ddata += workX_[kk*nInputs_+ii-1] *
                     workArray_[leng*nSamples_+kk*nSamples_+jj]; 
          ddata -= workX_[kk*nInputs_+ii-1];
          workX_[leng*nInputs_+kk*nBasis+ii] = ddata; 
          workX_[leng*nInputs_+leng*nBasis+kk*nBasis+ii] = ddata; 
        }
      }
      dpotrs_(&uplo, &nBasis, &leng, Mmatrix_, &nSamples_,
              &workX_[leng*nInputs_+leng*nBasis], &nBasis, &status);

      for (kk = 0; kk < leng; kk++)
      {
        ddata = 0.0;
        for (ii = 0; ii < nBasis; ii++)
          ddata += workX_[leng*nInputs_+kk*nBasis+ii] *
                   workX_[leng*nInputs_+(leng+kk)*nBasis+ii];
        ddata = KrigingVariance_*
                (1.0 + ddata-workArray_[2*leng*nSamples_+kk]);
        if (ddata < 0.0)
        {
          printf("Kriging WARNING: prediction variance < 0\n");
          ddata = 0.0;
        }
        YStds[ll+kk] = sqrt(ddata);
      }
    }
  }
#else
  nRows = nSamples_ + 1;
  if (pOrder_ == 1) nRows = nRows + nInputs_;
  nBasis = nRows - nSamples_;
  if (workArray_ == NULL)
  {
    workArray_ = new double[2*nSamples_*length+length];
    workX_ = new double[length*nInputs_+2*length*nBasis];
    workLength_ = length;
    if (workX_ == NULL)
    {
      printOutTS(PL_ERROR,"ERROR: memory allocation in file %s line %d\n",
                 __FILE__, __LINE__);
      abort();
    }
  }
  else if (workArray_ != NULL && length > workLength_)
  {
    delete [] workArray_;
    delete [] workX_;
    workLength_ = length;
    workArray_ = new double[2*nSamples_*length+length];
    workX_ = new double[length*nInputs_+2*length*nBasis];
    if (workX_ == NULL)
    {
      printOutTS(PL_ERROR,"ERROR: memory allocation in file %s line %d\n",
                 __FILE__, __LINE__);
      abort();
    }
  }
  for (ii = 0; ii < nInputs_; ii++)
  {
    mean  = XMeans_[ii];
    stdev = XStds_[ii];
    for (jj = 0; jj < length; jj++)
      workX_[jj*nInputs_+ii] = (X[jj*nInputs_+ii] - mean) / stdev;
  }
  for (kk = 0; kk < length; kk++)
  {
    for (jj = 0; jj < nSamples_; jj++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = XNormalized_[jj*nInputs_+ii] - workX_[kk*nInputs_+ii];
        dist += ddata * ddata / (Thetas_[ii] * Thetas_[ii]);
      }
      ddata = exp(-dist);
      if (ddata < 1.0e-16) ddata = 0.0;
      workArray_[kk*nSamples_+jj] = ddata;
      workArray_[nSamples_*length+kk*nSamples_+jj] = ddata;
    }
  }
  for (kk = 0; kk < length; kk++)
  {
    ddata = 1.0 * V1_[0];
    if (pOrder_ == 1)
      for (jj = 1; jj <= nInputs_; jj++) 
        ddata += workX_[kk*nInputs_+jj-1]*V1_[jj];
    for (jj = 0; jj < nSamples_; jj++) 
      ddata += workArray_[nSamples_*length+kk*nSamples_+jj] * V2_[jj];
    Y[kk] = ddata * YStd_ + YMean_;
  }
  if (YStds != NULL)
  {
    dpotrs_(&uplo, &nSamples_, &length, Rmatrix_, &nSamples_,
            &workArray_[length*nSamples_], &nSamples_, &status);
    for (kk = 0; kk < length; kk++)
    {
      ddata = 0.0;
      for (ii = 0; ii < nSamples_; ii++)
        ddata += workArray_[kk*nSamples_+ii] *
                 workArray_[length*nSamples_+kk*nSamples_+ii];
      workArray_[2*length*nSamples_+kk] = ddata; 
    }
    for (kk = 0; kk < length; kk++)
    {
      ddata = 0.0;
      for (ii = 0; ii < nSamples_; ii++)
        ddata += workArray_[length*nSamples_+kk*nSamples_+ii]; 
      workX_[length*nInputs_+kk*nBasis] = ddata - 1.0; 
      workX_[length*nInputs_+length*nBasis+kk*nBasis] = ddata - 1.0; 
      for (ii = 1; ii < nBasis; ii++)
      {
        ddata = 0.0;
        for (jj = 0; jj < nSamples_; jj++)
          ddata += workX_[kk*nInputs_+ii-1] *
                   workArray_[length*nSamples_+kk*nSamples_+jj]; 
        ddata -= workX_[kk*nInputs_+ii-1];
        workX_[length*nInputs_+kk*nBasis+ii] = ddata; 
        workX_[length*nInputs_+length*nBasis+kk*nBasis+ii] = ddata; 
      }
    }
    dpotrs_(&uplo, &nBasis, &length, Mmatrix_, &nSamples_,
            &workX_[length*nInputs_+length*nBasis], &nBasis, &status);

    for (kk = 0; kk < length; kk++)
    {
      ddata = 0.0;
      for (ii = 0; ii < nBasis; ii++)
        ddata += workX_[length*nInputs_+kk*nBasis+ii] *
                 workX_[length*nInputs_+(length+kk)*nBasis+ii];
      ddata = KrigingVariance_*
              (1.0 + ddata-workArray_[2*length*nSamples_+kk]);
      if (ddata < 0.0)
      {
        printf("Kriging WARNING: prediction variance < 0\n");
        ddata = 0.0;
      }
      YStds[kk] = sqrt(ddata);
    }
  }
#endif
  return 0.0;
}

// ************************************************************************
// compute distances between all pairs of inputs
// ------------------------------------------------------------------------
int Kriging::computeDistances(double **XDists, int *length)
{
  int    ii, jj, kk, count, error=0;
  double *LDists, dist;

  LDists = new double[(nSamples_*(nSamples_-1)/2)*nInputs_];
  checkAllocate(LDists, "LDists in Kriging::computeDistances");
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
        printf("Kriging ERROR: repeated sample points.\n");
        printf("               Prune repeated points and re-run.\n");
        printf("Sample %d : (with sample point %d)\n", kk+1, jj+1);
        for (ii = 0; ii < nInputs_; ii++)
          printf("   Input %d : %e\n",ii+1,
                 XNormalized_[kk*nInputs_+ii]*XStds_[ii]+XMeans_[ii]);
        error = 1;
      }
      count++;
    }
  }
  (*length) = count;
  (*XDists) = LDists;
  return error;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double Kriging::setParams(int targc, char **targv)
{
  int    ii, *iArray;
  double *lengthScales, mmax, range;
  char   pString[500];
  FILE   *fp=NULL;

  if (targc > 0 && !strcmp(targv[0], "noReuse"))
     noReuse_ = 1;
  else if (targc > 0 && !strcmp(targv[0], "setMode0"))
     fastMode_ = 0;
  else if (targc > 0 && !strcmp(targv[0], "setMode1"))
     fastMode_ = 1;
  else if (targc > 0 && !strcmp(targv[0], "setMode2"))
     fastMode_ = 2;
  else if (targc > 0 && !strcmp(targv[0], "setMode3"))
     fastMode_ = 3;
  else if (targc > 0 && !strcmp(targv[0], "rank"))
  {
    lengthScales = new double[nInputs_];
    checkAllocate(lengthScales, "lengthScales in Kriging::setParams");
    mmax = 0.0;
    for (ii = 0; ii < nInputs_; ii++)
    {
      lengthScales[ii] = 1.0 / Thetas_[ii];
      if (XMeans_[ii] == 0 && XStds_[ii] == 1)
      {
        range = upperBounds_[ii] - lowerBounds_[ii];
        lengthScales[ii] *= range;
      }
      if (lengthScales[ii] > mmax) mmax = lengthScales[ii];
    }
    for (ii = 0; ii < nInputs_; ii++)
      lengthScales[ii] = lengthScales[ii] / mmax * 100.0;

    if (psPlotTool_ == 1)
         fp = fopen("scilabkrisa.sci", "w");
    else fp = fopen("matlabkrisa.m", "w");
    if (fp == NULL)
    {
      printf("Kriging ERROR: something wrong with opening a write file.\n");
    }
    else
    {
      fprintf(fp, "n = %d;\n", nInputs_);
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs_; ii++)
        fprintf(fp, "%24.16e \n", PABS(lengthScales[ii]));
      fprintf(fp, "]; \n");
      fprintf(fp, "ymax = max(Y);\n");
      fprintf(fp, "ymin = 0;\n");
      fprintf(fp, "if (ymax == ymin)\n");
      fprintf(fp, "   ymax = ymax * 0.1;\n");
      fprintf(fp, "end;\n");
      fwritePlotCLF(fp);
      fprintf(fp, "bar(Y,0.8);\n");
      fwritePlotAxes(fp);
      sprintf(pString, "Kriging Ranking");
      fwritePlotTitle(fp, pString);
      sprintf(pString, "Input Numbers");
      fwritePlotXLabel(fp, pString);
      sprintf(pString, "Kriging Measure");
      fwritePlotYLabel(fp, pString);
      if (psPlotTool_ == 1)
      {
        fprintf(fp,"a=gca();\n");
        fprintf(fp,"a.data_bounds=[0,ymin; n+1,ymax+0.01*(ymax-ymin)];\n");
      }
      else
      {
        fprintf(fp,"axis([0 n+1 ymin ymax+0.01*(ymax-ymin)])\n");
      }
      fclose(fp);
      if (psPlotTool_ == 1)
           printf("Kriging ranking in file scilabkrisa.sci\n");
      else printf("Kriging ranking in file matlabkrisa.m\n");
    }
    iArray = new int[nInputs_];
    checkAllocate(iArray, "iArray in Kriging::setParams");
    for (ii = 0; ii < nInputs_; ii++) iArray[ii] = ii;
    sortDbleList2a(nInputs_, lengthScales, iArray);
    printAsterisks(PL_INFO, 0);
    printf("* Kriging screening rankings\n");
    printAsterisks(PL_INFO, 0);
    for (ii = nInputs_-1; ii >= 0; ii--)
      printf("*  Rank %3d : Input = %4d (score = %5.1f) (ref = %e)\n",
             nInputs_-ii, iArray[ii]+1, lengthScales[ii], 
             0.01*lengthScales[ii]*mmax);
    printAsterisks(PL_INFO, 0);
    delete [] lengthScales;
  }
  return 0.0;
}

// ************************************************************************
// perform optimization evaluation 
// ------------------------------------------------------------------------
void Kriging::optimize()
{
  int    ii, jj, mm, nBasis, nDists, deg; 
  double *TUppers, *TLowers, *Thetas, *XDists, objfcn, objtmp, ddata;
  char   pString[1000];

  pOrder_ = 1;
  nBasis  = 1;
  if (pOrder_ == 1) nBasis = nInputs_ + 1;
  TUppers = new double[nInputs_];
  TLowers = new double[nInputs_];
  checkAllocate(TLowers, "TLowers in Kriging::optimize");
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
    if (psMasterMode_ == 1 && psInteractive_ == 1)
    {
      printf("Kriging: for input %d :\n", ii+1);
      printf("Kriging: current optimization lower bound = %e\n",
             TLowers[ii]);
      sprintf(pString,
         "Kriging: Enter new optimization lower bound : ");
      TLowers[ii] = getDouble(pString);
      if (TLowers[ii] <= 0.0)
      {
        printf("Kriging ERROR: lower bound <= 0\n");
        exit(1);
      }
      printf("Kriging: current optimization upper bound = %e\n",
             TUppers[ii]);
      sprintf(pString,
         "Kriging: Enter new optimization upper bound : ");
      TUppers[ii] = getDouble(pString);
      if (TUppers[ii] <= 0.0)
      {
        printf("Kriging ERROR: upper bound <= 0\n");
        exit(1);
      }
      if (TLowers[ii] >= TUppers[ii])
      {
        printf("Kriging ERROR: upper bound <= lower bound.\n");
        exit(1);
      }
    }
  }
  computeDistances(&XDists, &nDists);
  KRI_SMatrix = new double[nSamples_*nSamples_];
  KRI_FMatrix = new double[nSamples_*nBasis];
  KRI_FMatTmp = new double[nSamples_*nBasis];
  KRI_XDists  = XDists;
  KRI_iter    = 0;
  betas_      = new double[nBasis];
  betasOpt_   = new double[nBasis];
  gammas_     = new double[nSamples_];
  gammasOpt_  = new double[nSamples_];
  thetasOpt_  = new double[nInputs_];

  Thetas = new double[nInputs_+1];
  checkAllocate(Thetas, "Thetas in Kriging::optimize");
  for (ii = 0; ii < nInputs_; ii++) Thetas[ii] = 10.0;
  objfcn = evaluateFunction(Thetas);

  int maxIts = 4, flag;
  if (nInputs_ < maxIts) maxIts = nInputs_;
  if (nInputs_ == 1) maxIts = 2;
  double *scales = new double[nInputs_];
  checkAllocate(scales, "scales in Kriging::optimize");
  for (ii = 0; ii < nInputs_; ii++) 
    scales[ii] = pow(2.0, (1.0+ii)/(2.0+nInputs_));

  double *Thetas2 = new double[nInputs_];
  double *T1      = new double[nInputs_];
  double *S1      = new double[nInputs_];
  checkAllocate(S1, "S1 in Kriging::optimize");
  for (mm = 0; mm < maxIts; mm++)
  {
    for (ii = 0; ii < nInputs_; ii++) Thetas2[ii] = Thetas[ii];
     
    for (ii = 0; ii < nInputs_; ii++)
    {
      for (jj = 0; jj < nInputs_; jj++) T1[jj] = Thetas[jj];
      if (Thetas[ii] >= TUppers[ii])
      {
        flag = 1;
        T1[ii] = Thetas[ii] / sqrt(scales[ii]);
      }
      else if (Thetas[ii] == TLowers[ii])
      {
        flag = 1;
        T1[ii] = Thetas[ii] * sqrt(scales[ii]);
      }
      else
      {
        flag = 0;
        T1[ii] = Thetas[ii] * scales[ii];
        if (TUppers[ii] < T1[ii]) T1[ii] = TUppers[ii];
      }
      objtmp = evaluateFunction(T1);
      if (objtmp < objfcn)
      {
        for (jj = 0; jj < nInputs_; jj++) Thetas[jj] = T1[jj];
        objfcn = objtmp;
        for (jj = 0; jj < nBasis; jj++) betasOpt_[jj] = betas_[jj];
        for (jj = 0; jj < nSamples_; jj++) gammasOpt_[jj] = gammas_[jj];
        for (jj = 0; jj < nInputs_; jj++) thetasOpt_[jj] = Thetas[jj];
      }
      else
      {
        if (flag == 0)
        {
          T1[ii] = Thetas[ii] / scales[ii];
          if (T1[ii] < TLowers[ii]) T1[ii] = TLowers[ii];
          objtmp = evaluateFunction(T1);
          if (objtmp < objfcn)
          {
            for (jj = 0; jj < nInputs_; jj++) Thetas[jj] = T1[jj];
            objfcn = objtmp;
            for (jj = 0; jj < nBasis; jj++) betasOpt_[jj] = betas_[jj];
            for (jj = 0; jj < nSamples_; jj++) gammasOpt_[jj] = gammas_[jj];
            for (jj = 0; jj < nInputs_; jj++) thetasOpt_[jj] = Thetas[jj];
          }
        }
      }
    }

    flag = 1;
    deg  = 1;
    for (ii = 0; ii < nInputs_; ii++) S1[ii] = Thetas[ii]/Thetas2[ii];
    while (flag == 1)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = pow(S1[ii],1.0*deg) * Thetas[ii];
        if (ddata < TLowers[ii]) ddata = TLowers[ii]; 
        if (ddata > TUppers[ii]) ddata = TUppers[ii]; 
        T1[ii] = ddata;
      }
      objtmp = evaluateFunction(T1);
      if (objtmp < objfcn)
      {
        for (ii = 0; ii < nInputs_; ii++) Thetas[ii] = T1[ii];
        objfcn = objtmp;
        deg *= 2;
        for (jj = 0; jj < nBasis; jj++) betasOpt_[jj] = betas_[jj];
        for (jj = 0; jj < nSamples_; jj++) gammasOpt_[jj] = gammas_[jj];
        for (jj = 0; jj < nInputs_; jj++) thetasOpt_[jj] = Thetas[jj];
      }
      else flag = 0;
      for (ii = 0; ii < nInputs_; ii++) 
      {
        if (T1[ii] <= TLowers[ii]) flag = 0;
        if (T1[ii] >= TUppers[ii]) flag = 0;
      }
    }
    ddata = scales[0];
    for (ii = 0; ii < nInputs_-1; ii++) scales[ii] = scales[ii+1];
    scales[nInputs_-1] = ddata;
    for (ii = 0; ii < nInputs_; ii++) scales[ii] = pow(scales[ii],0.25);
  }

  delete [] T1;
  delete [] S1;
  delete [] scales;
  delete [] Thetas;
  delete [] Thetas2;
  delete [] TUppers;
  delete [] TLowers;
  delete [] KRI_SMatrix;
  delete [] KRI_FMatrix;
  delete [] KRI_FMatTmp;
  delete [] KRI_XDists;
  KRI_SMatrix = NULL;
  KRI_FMatrix = NULL;
  KRI_FMatTmp = NULL;
  KRI_XDists  = NULL;
}

// ************************************************************************
// function evaluation 
// ------------------------------------------------------------------------
double Kriging::evaluateFunction(double *thetas)
{
  int    nBasis, count, ii, jj, kk, status, inc=1;
  double dist, ddata, retVal;
  FILE   *fp=NULL;
  char   uplo='L', trans='N', diag='N', side='L';

  KRI_iter++;
  nBasis = 1;
  if (pOrder_ == 1) nBasis = nInputs_ + 1;

  count = 0;
  for (jj = 0; jj < nSamples_; jj++)
  {
    //KRI_SMatrix[jj*nSamples_+jj] = 1.0 + 1e-15 * nSamples_;
    KRI_SMatrix[jj*nSamples_+jj] = 1.0 + 2.22e-16 * (nSamples_ + 10);
    if (KRI_nugget != 0.0) 
      KRI_SMatrix[jj*nSamples_+jj] += KRI_nugget;
    else if (nInputs_ == 2) 
      KRI_SMatrix[jj*nSamples_+jj] += 0.01;
    if (KRI_dataStdDevs != NULL) 
    {
      ddata = pow(KRI_dataStdDevs[jj]/KRI_YStd,KRI_Exponent);
      KRI_SMatrix[jj*nSamples_+jj] += ddata;
    }
    for (kk = jj+1; kk < nSamples_; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
        dist += pow(KRI_XDists[count*nInputs_+ii],KRI_Exponent)*
                thetas[ii];
      dist = exp(-dist);
      //if (dist < 1.0e-16) dist = 0;
      KRI_SMatrix[jj*nSamples_+kk] = dist;
      KRI_SMatrix[kk*nSamples_+jj] = dist;
      count++;
    }
  }
  dpotrf_(&uplo, &nSamples_, KRI_SMatrix, &nSamples_, &status);
  if (status != 0) 
  {
    printf("Kriging ERROR: Cholesky decomposition not successful.\n");
    exit(1);
  }

  for (ii = 0; ii < nSamples_; ii++) 
  {
    KRI_FMatrix[ii] = 1.0;
    KRI_FMatTmp[ii] = 1.0;
    if (pOrder_ == 1)
    {
      for (kk = 1; kk < nBasis; kk++) 
      {
        ddata = XNormalized_[ii*nInputs_+kk-1];
        KRI_FMatTmp[kk*nSamples_+ii] = ddata;
        KRI_FMatrix[kk*nSamples_+ii] = ddata;
      }
    }
  }
  for (ii = 0; ii < nBasis; ii++) 
  { 
    dtrsv_(&uplo,&trans,&diag,&nSamples_,KRI_SMatrix,&nSamples_,
           &KRI_FMatTmp[ii*nSamples_], &inc);
  }
  for (ii = 0; ii < nSamples_*nBasis; ii++) 
    KRI_FMatrix[ii] = KRI_FMatTmp[ii];

  double *tau  = new double[nSamples_];
  double *work = new double[nSamples_];
  double *QMat = new double[nSamples_*nBasis];
  checkAllocate(QMat, "QMat in Kriging::evaluateFunction");
  dgeqrf_(&nSamples_, &nBasis, KRI_FMatTmp, &nSamples_, 
          tau, work, &nSamples_, &status);
  for (ii = 0; ii < nSamples_*nBasis; ii++) QMat[ii] = 0.0;
  for (ii = 0; ii < nBasis; ii++) QMat[ii*nSamples_+ii] = 1.0;
  dormqr_(&side, &trans, &nSamples_, &nBasis, &nBasis, KRI_FMatTmp, 
          &nSamples_, tau, QMat, &nSamples_, work, &nSamples_, 
          &status);
  delete [] work;
  delete [] tau;

  double *Ytmp = new double[nSamples_];
  checkAllocate(Ytmp, "Ytmp in Kriging::evaluateFunction");
  for (ii = 0; ii < nSamples_; ii++) Ytmp[ii] = YNormalized_[ii];
  dtrsv_(&uplo,&trans,&diag,&nSamples_,KRI_SMatrix,&nSamples_,
         Ytmp, &inc);
  work = new double[nBasis];
  checkAllocate(work, "work in Kriging::evaluateFunction");
  for (jj = 0; jj < nBasis; jj++)
  {
    ddata = 0.0;
    for (kk = 0; kk < nSamples_; kk++)
      ddata += QMat[kk+jj*nSamples_] * Ytmp[kk];
    work[jj] = ddata;
  }
  for (jj = nBasis-1; jj >= 0; jj--)
  {
    ddata = work[jj];
    for (kk = jj+1; kk < nBasis; kk++)
      ddata -= betas_[kk] * KRI_FMatTmp[kk*nSamples_+jj];
    betas_[jj] = ddata/ KRI_FMatTmp[jj*nSamples_+jj];
  }
  delete [] work;
  double *rho = new double[nSamples_];
  checkAllocate(rho, "rho in Kriging::evaluateFunction");
  for (ii = 0; ii < nSamples_; ii++) rho[ii] = Ytmp[ii];
  for (ii = 0; ii < nSamples_; ii++) 
  {
    ddata = 0.0;
    for (kk = 0; kk < nBasis; kk++) 
      ddata += KRI_FMatrix[kk*nSamples_+ii] * betas_[kk];
    rho[ii] -= ddata;
  }
  double sigma = 0.0;
  //for (ii = 0; ii < nSamples_; ii++) sigma += rho[ii] * rho[ii];
  for (ii = 0; ii < nSamples_; ii++) sigma += pow(rho[ii],2.0);
  sigma /= (double) nSamples_;
  double detR = 1.0;
  for (ii = 0; ii < nSamples_; ii++) 
  {
    ddata = KRI_SMatrix[ii*nSamples_+ii];
    detR *= pow(ddata, 2.0/nSamples_);
  }
  char transT = 'T';
  dtrsv_(&uplo,&transT,&diag,&nSamples_,KRI_SMatrix,&nSamples_,
         rho, &inc);
  for (ii = 0; ii < nSamples_; ii++) gammas_[ii] = rho[ii];
  delete [] QMat;
  delete [] Ytmp;
  delete [] rho;
  retVal = sigma * detR;
  return retVal;
}

// ************************************************************************
// predict 
// ------------------------------------------------------------------------
double Kriging::predict0(int length, double *X, double *Y, double *YStds)
{
  int    ii, jj, kk, nBasis;
  double *Xt = new double[nInputs_];
  double *FMat, ddata, Yt, *RMat;

  nBasis = 1;
  if (pOrder_ >= 1) nBasis += nInputs_;
  FMat = new double[nBasis];
  RMat = new double[nSamples_];
  checkAllocate(RMat, "RMat in Kriging::predict0");
  if (YStds != NULL)
    for (ii = 0; ii < length; ii++) YStds[ii] = 0;

  for (ii = 0; ii < length; ii++) 
  {
    for (jj = 0; jj < nInputs_; jj++) 
      Xt[jj] = (X[ii*nInputs_+jj] - XMeans_[jj])/XStds_[jj];

    //*/ construct regression matrix F
    FMat[0] = 1.0;
    if (pOrder_ == 1)
      for (kk = 1; kk < nBasis; kk++) FMat[kk] = Xt[kk-1];
    //*/ construct correlation R
    for (jj = 0; jj < nSamples_; jj++)
    {
      ddata = 0.0;
      for (kk = 0; kk < nInputs_; kk++)
        ddata += pow(Xt[kk]-XNormalized_[jj*nInputs_+kk],
                 KRI_Exponent)*thetasOpt_[kk];
      RMat[jj] = exp(-ddata);
    }
    Yt = 0.0;
    for (jj = 0; jj < nBasis; jj++) Yt += FMat[jj] * betasOpt_[jj];
    for (jj = 0; jj < nSamples_; jj++) Yt += RMat[jj] * gammasOpt_[jj];
    Y[ii] = Yt * YStd_ + YMean_;
  }
  delete [] FMat;
  delete [] RMat;
  delete [] Xt;
  return 0;
}

