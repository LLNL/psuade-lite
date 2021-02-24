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
// Functions for the class OUUOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <PsuadeCmakeConfig.h>

#ifdef WINDOWS
#include <windows.h>
#undef ERROR
#undef IS_ERROR
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "PDFManager.h"
#include "OUUOptimizer.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Sampling.h"
#include "FuncApprox.h"
#include "PrintingTS.h"
#ifdef HAVE_COBYLA
#include "cobyla.h"
#endif
#ifdef HAVE_LBFGS
extern "C" {
#include "../../External/L-BFGS-B-C/src/lbfgsb.h"
}
#endif
#ifdef HAVE_NOMAD
#include "../../External/NOMAD/src/nomad.hpp"
#include "../../External/NOMAD/src/Display.hpp"
using namespace NOMAD;
#endif

// ************************************************************************
// External functions
// ------------------------------------------------------------------------
extern "C" void bobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);
extern "C" void obobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);
extern "C"  void newuoa_(int *,int *,double *,double *,double *,int *,
                         int *,double*);
extern "C" void lincoa_(int *,int *, int *, double *, int *, double *,
                        double *,double *,double *, int *, int *, double*);

// ************************************************************************
// Internal 'global' variables (for passing parameters to Fortran and C
// functions - bobyqa, newuoa and cobyla)
// ------------------------------------------------------------------------
void  *psOUUObj_=NULL;
int psOUUM1_=-1;
int psOUUM2_=-1;
int psOUUM3_=-1;
int psOUUM4_=-1;
int psOUUM_=-1;
int psOUUnOutputs_=-1;

int psOUUUserOpt_ = 0;
int psOUUUseRS_=0;
int psOUUZ4RSType_=0;
int psOUUZ4RSAux_=0;
int psOUUValidateRS_=0;

double *psOUUM2Values_=NULL;
int psOUUZ3nSamples_=-1;
int psOUUZ4nSamples_=-1;
double *psOUUZ3SamInputs_=NULL;
double *psOUUZ4SamInputs_=NULL;
double *psOUUSamOutputs_=NULL;
double *psOUUZ4LBounds_=NULL;
double *psOUUZ4UBounds_=NULL;
double *psOUULBounds_=NULL;
double *psOUUUBounds_=NULL;
double *psOUUSamProbs_=NULL;

double *psOUUSamConstrNorms_=NULL;

double *psOUUXValues_=NULL;
double *psOUUWValues_=NULL;
double *psOUUOptimalX_=NULL;
double psOUUOptimalY_=0.0;

FuncApprox *psOUUfaPtr_=NULL;

int psOUULargeSampleSize_=0;
double *psOUULargeSamInputs_=NULL;
double *psOUULargeSamOutputs_=NULL;

int psOUUMode_=1;
double psOUUPercentile_=0.5;
double psOUUStdevMultiplier_=0;

int psOUUCMode_=1;
int psOUUOptCode_=-1;

int psOUUCounter_=0;
int psOUUParallel_=0;
int psOUUEnsembleEval_=0;
int psOUUStop_=0;
double psOUUTolerance_=1.0e-5;

#define psOUUMaxSaved_ 10000
int    psOUUSaveHistory_=0;
int    psOUUNSaved_=0;
double psOUUSaveX_[psOUUMaxSaved_*10];
double psOUUSaveY_[psOUUMaxSaved_*10];
double *psOUUNomadObjFcnStore_=NULL;

int *psOUUInputTypes_=NULL;
int psOUUMasterMode_=0;
int *psOUUDesignTypes_=NULL;
int psOUUNomadCnt_=0;

int psOUULincoaNConstr_ = 0;
double *psOUULincoaAmat_ = NULL;
double *psOUULincoaBvec_ = NULL;
char    psOUULincoaConstraintFile_[1000]="NONE";

#define PABS(x)  ((x) > 0 ? x : -(x))
#define psOUUType1 1
#define psOUUType2 2
#define psOUUType3 3
#define psOUUType4 4

// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" 
{
#endif
  void *ouuevalfunc2_(int *nInps, double *XValues, double *YValue)
  {
    int    ii, kk, M1, M2, M3, M4, M, index, iOne=1, found, funcID;
    double ddata;
    oData  *odata;
    FILE   *fp=NULL;

    odata = (oData *) psOUUObj_;
    M1    = psOUUM1_;
    M2    = psOUUM2_;
    M3    = psOUUM3_;
    M4    = psOUUM4_;
    M     = M1 + M2 + M3 + M4;
    funcID = odata->numFuncEvals_;

    index = 0;
    for (ii = 0; ii < M; ii++) 
    {
       if (psOUUInputTypes_[ii] == psOUUType2)
       {
          psOUUXValues_[ii] = XValues[index];
          index++;
       }
    }
    index = 0;
    for (ii = 0; ii < M; ii++) 
    {
      if (psOUUInputTypes_[ii] == psOUUType3)
      {
        psOUUXValues_[ii] = psOUUWValues_[index];
        index++;
      }
    }
    for (ii = 0; ii < M; ii++) 
    {
      if (psOUUInputTypes_[ii] == psOUUType4)
      {
        psOUUXValues_[ii] = psOUUWValues_[index];
        index++;
      }
    }

    found = 0;
    for (ii = 0; ii < psOUUNSaved_; ii++)
    {
      for (kk = 0; kk < M; kk++)
        if (PABS(psOUUSaveX_[ii*M+kk]-psOUUXValues_[kk])>1.0e-14) 
          break;
      if (kk == M)
      {
        found = 1;
        if (odata->outputLevel_ > 2)
          printf("OUUOptimizer: simulation results reuse.\n");
        ddata = psOUUSaveY_[ii];
        break;
      }
    }

    if (found == 0)
    {
      odata->funcIO_->evaluate(funcID,M,psOUUXValues_,iOne,&ddata,0);
      odata->numFuncEvals_++;
      if (psOUUSaveHistory_ == 1 && (psOUUNSaved_+1)*M < psOUUMaxSaved_*10)
      {
        for (kk = 0; kk < M; kk++)
           psOUUSaveX_[psOUUNSaved_*M+kk] = psOUUXValues_[kk];
        psOUUSaveY_[psOUUNSaved_] = ddata;
        psOUUNSaved_++;
        fp = fopen("psuade_ouu_history","w");
        if (fp != NULL)
        {
          for (ii = 0; ii < psOUUNSaved_; ii++)
          {
            fprintf(fp, "999 %d ", M);
            for (kk = 0; kk < M; kk++)
              fprintf(fp, "%24.16e ", psOUUSaveX_[ii*M+kk]);
            fprintf(fp, "%24.16e\n", psOUUSaveY_[ii]);
          }
          fclose(fp);
        }
      }
    }
    (*YValue) = ddata;

    if (ddata < psOUUOptimalY_)
    {
      psOUUOptimalY_ = ddata;
      for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = psOUUXValues_[ii];
      if (isScreenDumpModeOn() && odata->outputLevel_ > 2)
      {
        printf("    OUUOptimizer inner loop New Ymin = %16.8e (%d)\n",
               ddata, odata->numFuncEvals_);
        for (ii = 0; ii < M; ii++) 
          if (psOUUInputTypes_[ii] == psOUUType2)
            printf("     Input %3d = %e\n", ii+1, psOUUOptimalX_[ii]);
      }
    }
    return NULL;
  }

  //* -------------------------------------------------------------------
  void *ouuevalfunc_(int *nInps, double *XValues, double *YValue)
  {
    int    ii, kk, ss, funcID, M, M1, M2, M3, M4, bobyqaFlag=1112, nPts;
    int    maxfun, iOne=1, *readys, status, index, nSamp;
    double rhobeg, rhoend, ddata, *workArray, mean, stdev, *XLocal;
    double *lowers, *uppers;
    char   winput[1000];
    oData  *odata;
    FILE   *fp=NULL;

    odata = (oData *) psOUUObj_;
    M1    = psOUUM1_;
    M2    = psOUUM2_;
    M3    = psOUUM3_;
    M4    = psOUUM4_;
    M     = M1 + M2 + M3 + M4;

    fp = fopen("psuade_ouu_stop","r");
    if (fp != NULL)
    {
      unlink("psuade_ouu_stop");
      printf("OUUOptimizer: psuade_ouu_stop file found.\n");
      printf("              Abrupt termination.\n");
      if (psOUUNSaved_ > 0)
      {
        fp = fopen("psuade_ouu_history","w");
        if (fp != NULL)
        {
          for (ii = 0; ii < psOUUNSaved_; ii++)
          {
            fprintf(fp, "999 %d ", M);
            for (kk = 0; kk < M; kk++)
              fprintf(fp, "%24.16e ", psOUUSaveX_[ii*M+kk]);
            fprintf(fp, "%24.16e\n", psOUUSaveY_[ii]);
          }
          fclose(fp);
        }
        printf("OUUOptimizer: history saved in psuade_ouu_history\n");
      }
      exit(1);
    }

    nSamp = psOUUZ3nSamples_ * psOUUZ4nSamples_;
    XLocal = (double *) malloc(M*sizeof(double));
    readys = (int *) malloc(nSamp*sizeof(int));

    index = 0;
    for (ii = 0; ii < M; ii++)
    {
      if (psOUUInputTypes_[ii] == psOUUType1)
      {
        psOUUXValues_[ii] = XValues[index];
        index++;
      }
    }
    if (isScreenDumpModeOn())
    {
      printf("OUUOptimizer: Outer optimization iteration = %d\n",
             odata->numFuncEvals_/nSamp+1);
      if (odata->outputLevel_ > 1) 
      {
        printf("OUUOptimizer: Outer optimization loop FuncEval %d.\n",
               odata->numFuncEvals_+1);
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType1)
          {
            printf("    Current Level 1 input %3d = %e\n", ii+1, 
                   XValues[index]);
            index++;
          }
        }
      }
    }

    if (psOUUEnsembleEval_ == 0)
    {
      for (ss = 0; ss < nSamp; ss++)
      {
        if (isScreenDumpModeOn() && odata->outputLevel_ > 3) 
          printf("OUUOptimizer sample %d (of %d)\n",ss+1,nSamp); 
        readys[ss] = -1;
        psOUUOptimalY_ = PSUADE_UNDEFINED;
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType3)
          {
            kk = ss / psOUUZ4nSamples_;
            psOUUWValues_[index] = psOUUZ3SamInputs_[kk*M3+index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType4)
          {
            kk = ss % psOUUZ4nSamples_;
            psOUUWValues_[index] = psOUUZ4SamInputs_[kk*M4+index];
            index++;
          }
        }
        if (psOUUUserOpt_ == 0 && psOUUM2_ > 0)
        {
          lowers = (double *) malloc(M*sizeof(double));
          uppers = (double *) malloc(M*sizeof(double));
          rhobeg = 1.0e35;
          index = 0;
          for (ii = 0; ii < M; ii++)
          {
            if (psOUUInputTypes_[ii] == psOUUType2)
            {
              lowers[index] = odata->lowerBounds_[ii];
              uppers[index] = odata->upperBounds_[ii];
              //Nov 2017 - set to user-specified values
              //ddata = 0.5 * (lowers[index] + uppers[index]);
              XLocal[index] = psOUUM2Values_[index];
              ddata = uppers[index] - lowers[index];
              if (ddata < rhobeg) rhobeg = ddata;
              index++;
            }
          }
          rhobeg *= 0.5;
          rhoend = rhobeg * odata->tolerance_;
          if (rhobeg < rhoend) rhoend = rhobeg * 1.0e-6;
          maxfun = odata->maxFEval_;
          nPts = (M2 + 1) * (M2 + 2) / 2;
          kk = (nPts + 5) * (nPts + M2) + 3 * M2 * (M2 + 5) / 2 + 1;
          workArray = (double *) malloc(kk * sizeof(double));
          if (isScreenDumpModeOn() && odata->outputLevel_ > 3) 
            printf("OUU: inner optimization begins (sample %d of %d)\n",
                   ss+1,nSamp);
          bobyqaFlag = 1112;
          obobyqa_(&M2,&nPts,XLocal,lowers,uppers,&rhobeg,&rhoend, 
                   &bobyqaFlag, &maxfun, workArray);
          psOUUSamOutputs_[ss] = psOUUOptimalY_;
          if (isScreenDumpModeOn() && odata->outputLevel_ > 3) 
            printf("OUU: inner optimization %d ends, best Y = %e\n",
                   ss+1,psOUUSamOutputs_[ss]); 
          free(workArray);
          free(lowers);
          free(uppers);
        }
        else
        {
          index = 0;
          for (ii = 0; ii < M; ii++)
          {
            if (psOUUInputTypes_[ii] == psOUUType1)
            {
              XLocal[ii] = psOUUXValues_[index];
              index++;
            }
          }
          index = 0;
          for (ii = 0; ii < M; ii++)
          {
            if (psOUUInputTypes_[ii] == psOUUType2)
            {
              //Nov 2017 - set to user-specified values
              //XLocal[ii] = 0.5 * (odata->lowerBounds_[ii] + 
              //                    odata->upperBounds_[ii]);
              XLocal[ii] = psOUUM2Values_[index];
              index++;
            }
          }
          index = 0;
          for (ii = 0; ii < M; ii++)
          {
            if (psOUUInputTypes_[ii] == psOUUType3)
            {
              kk = ss / psOUUZ4nSamples_;
              XLocal[ii] = psOUUZ3SamInputs_[kk*M3+index];
              index++;
            }
          }
          index = 0;
          for (ii = 0; ii < M; ii++)
          {
            if (psOUUInputTypes_[ii] == psOUUType4)
            {
              kk = ss % psOUUZ4nSamples_;
              XLocal[ii] = psOUUZ4SamInputs_[kk*M4+index];
              index++;
            }
          }
          int found = 0;
          for (kk = 0; kk < psOUUNSaved_; kk++)
          {
            for (ii = 0; ii < M; ii++)
            {
              if (psOUUInputTypes_[ii] != psOUUType2 &&
                  (PABS(psOUUSaveX_[kk*M+ii]-XLocal[ii])>1.0e-14)) 
                break;
            }
            if (ii == M)
            {
              found = 1;
              if (isScreenDumpModeOn() && odata->outputLevel_ > 2) 
                printf("OUUOptimizer: simulation results reuse.\n");
              for (ii = 0; ii < M; ii++)
              {
                if (psOUUInputTypes_[ii] == psOUUType2)
                  XLocal[ii] = psOUUSaveX_[kk*M+ii];
              }   
              psOUUSamOutputs_[ss] = psOUUSaveY_[kk];
              for (ii = 0; ii < M; ii++) 
                psOUUOptimalX_[ii] = XLocal[ii];
              readys[ss] = 0;
              break;
            }
          }
          if (found == 0)
          {
            funcID = psOUUCounter_ * nSamp + ss;
            if (odata->optFunction_ != NULL)
            {
              odata->optFunction_(M, XLocal, iOne, &ddata);
              readys[ss] = 0;
            }
            else if (odata->funcIO_ != NULL)
            {
              readys[ss] = odata->funcIO_->evaluate(funcID,M,XLocal,
                                                    iOne,&ddata,0);
            }
            psOUUSamOutputs_[ss] = ddata;
            for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = XLocal[ii];
            if (isScreenDumpModeOn() && odata->outputLevel_ > 3 && 
                readys[ss] == 0) 
              printf("OUUOptimizer sample %d completed, best Y = %e\n",
                     ss+1,psOUUSamOutputs_[ss]); 
            if (psOUUSaveHistory_ == 1 && 
                (psOUUNSaved_+1)*M < psOUUMaxSaved_*10)
            {
              for (ii = 0; ii < M; ii++)
                psOUUSaveX_[psOUUNSaved_*M+ii] = XLocal[ii];
              psOUUSaveY_[psOUUNSaved_] = psOUUSamOutputs_[ss];
              psOUUNSaved_++;
            }
            if (psOUUParallel_ == 0) odata->numFuncEvals_++;
          }
        }
      }

      if (psOUUUserOpt_ == 1 && psOUUParallel_ == 1)
      {
        for (ss = 0; ss < nSamp; ss++)
        {
          funcID = psOUUCounter_ * nSamp + ss;
          if (readys[ss] != 0)
          {
            while (readys[ss] != 0)
            {
#ifdef WINDOWS
              Sleep(1000);
#else
              sleep(1);
#endif
              readys[ss] = odata->funcIO_->evaluate(funcID,M,XLocal,
                                                    iOne,&ddata,2);
            }
            psOUUSamOutputs_[ss] = ddata;
            if (isScreenDumpModeOn() && odata->outputLevel_ > 3) 
              printf("OUUOptimizer sample %d completed, best Y = %e\n",
                      ss+1,psOUUSamOutputs_[ss]); 
            odata->numFuncEvals_++;
          }
        }
      }
    }
    else
    {
      index = 0;
      for (ii = 0; ii < M; ii++)
      {
        if (psOUUInputTypes_[ii] == psOUUType1)
        {
          XLocal[index] = psOUUXValues_[ii];
          index++;
        }
      }
      for (ss = 0; ss < nSamp; ss++)
      {
        index = 0;
        for (ii = 0; ii < M; ii++) 
        {
          if (psOUUInputTypes_[ii] == psOUUType1)
          {
            psOUUXValues_[ss*M+ii] = XLocal[index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType2)
          {
            //Nov 2017 - set to user-specified values
            //psOUUXValues_[ss*M+ii] = 0.5*(odata->lowerBounds_[ii] + 
            //                              odata->upperBounds_[ii]);
            psOUUXValues_[ss*M+ii] = psOUUM2Values_[index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType3)
          {
            kk = ss / psOUUZ4nSamples_;
            psOUUXValues_[ss*M+ii] = psOUUZ3SamInputs_[kk*M3+index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType4)
          {
            kk = ss % psOUUZ4nSamples_;
            psOUUXValues_[ss*M+ii] = psOUUZ4SamInputs_[kk*M4+index];
            index++;
          }
        }
      }
      funcID = odata->numFuncEvals_;
      odata->funcIO_->ensembleEvaluate(nSamp,M,psOUUXValues_,
                                       iOne,psOUUSamOutputs_,funcID);
      odata->numFuncEvals_ += nSamp;
      for (ss = 0; ss < nSamp; ss++)
      {
        if (psOUUSaveHistory_ == 1 && 
           (psOUUNSaved_+1)*M < psOUUMaxSaved_*10)
        {
          for (kk = 0; kk < M; kk++)
            psOUUSaveX_[psOUUNSaved_*M+kk] = psOUUXValues_[ss*M+kk];
          psOUUSaveY_[psOUUNSaved_] = psOUUSamOutputs_[ss];
          psOUUNSaved_++;
        }
      }
      for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = psOUUXValues_[ii];
    }

    int failCnt=0;
    for (ss = 0; ss < nSamp; ss++) 
      if (psOUUSamOutputs_[ss] >= 0.98*PSUADE_UNDEFINED) failCnt++;
    if (failCnt != 0)
      printf("WARNING: there are %d failed runs out of %d\n",failCnt,nSamp);
    if (failCnt == nSamp || (failCnt > 0 && psOUUSamProbs_ != NULL))
    {   
      printf("ERROR: X3 sample cannot admit failures ==> terminate.\n");
      exit(1);
    }
    if (failCnt > 0 && psOUUMode_ > 2)
    {   
      printf("ERROR: OUU mode %d cannot admit failures ==> terminate.\n",
             psOUUMode_);
      exit(1);
    }
    if (psOUUUseRS_ == 0 || nSamp == 1)
    {
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1) 
        printf("OUUOpt: computing objective (no RS), nFuncEval = %d\n",
               odata->numFuncEvals_);
      if (psOUUMode_ == 1 || psOUUMode_ == 2)
      {
        mean = 0.0;
        if (psOUUZ3nSamples_ == 1)
        {
          for (ss = 0; ss < nSamp; ss++) 
            if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED) 
              mean += psOUUSamOutputs_[ss] / (double) (nSamp - failCnt);
        }
        else
        {
          for (ss = 0; ss < nSamp; ss++) 
          {
            index = ss / psOUUZ4nSamples_;
            mean += psOUUSamOutputs_[ss] / psOUUZ4nSamples_ * 
                    psOUUSamProbs_[index];
          }
        }
        (*YValue) = mean;
      }
      if (psOUUMode_ == 2 && psOUUStdevMultiplier_ != 0.0)
      {
        stdev = 0.0;
        for (ss = 0; ss < nSamp; ss++) 
        {
          index = ss / psOUUZ4nSamples_;
          if (psOUUZ3nSamples_ == 1)
          {
            if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED) 
              stdev += pow(psOUUSamOutputs_[ss]-mean,2.0)/
                          (double) (nSamp - failCnt);
          }
          else
          {
            stdev += pow(psOUUSamOutputs_[ss]-mean, 2.0)*
                     psOUUSamProbs_[index] / (double) psOUUZ4nSamples_;
          }
        }
        (*YValue) = mean + psOUUStdevMultiplier_ * sqrt(stdev);
      }
      if (psOUUMode_ == 3)
      {
        sortDbleList(nSamp, psOUUSamOutputs_);
        kk = (int) (psOUUPercentile_ * nSamp);
        if (kk >= nSamp) kk = nSamp - 1;
        (*YValue) = psOUUSamOutputs_[kk];
      }
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1) 
      {
        printf("OUUOptimizer: computed  objective (no RS) = %e.\n",
               (*YValue));
      }
    }
    else
    {
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1) 
        printf("OUUOpt: computing objective with RS, nFuncEval = %d\n",
               odata->numFuncEvals_);
        
      double totalMean=0.0, totalStdv=0.0, *resultStore, cverrors[3];
      resultStore=(double *) malloc(psOUUZ3nSamples_*sizeof(double));
      for (ii = 0; ii < psOUUZ3nSamples_; ii++)
      {
        if (psGMMode_ == 1)
        {
          fp = fopen("ouuZ4Sample", "w");
          if (fp != NULL)
          {
            fprintf(fp,"%d %d 1\n",psOUUZ4nSamples_,M4);
            for (ss = 0; ss < psOUUZ4nSamples_; ss++) 
            {
              for (kk = 0; kk < M4; kk++) 
                fprintf(fp,"%24.16e ",psOUUZ4SamInputs_[ss*M4+kk]);
              fprintf(fp,"%24.16e\n",
                      psOUUSamOutputs_[ii*psOUUZ4nSamples_+ss]);
            }
          }
          fclose(fp);
          printf("OUU INFO: a Z4 sample is ready for viewing in file");
          printf(" ouuZ4Sample.\n");
          printf("To continue, enter 0 (or 1 if no more interruption) : ");
          scanf("%d", &kk);
          if (kk == 1) psGMMode_ = 0;
        }
        if (isScreenDumpModeOn() && odata->outputLevel_ > 1) 
          printf("OUUOptimizer: constructing RS\n");
        status = psOUUfaPtr_->initialize(psOUUZ4SamInputs_,
                                &psOUUSamOutputs_[ii*psOUUZ4nSamples_]);
        if (psOUUValidateRS_ == 1)
        {
          validate(psOUUZ4nSamples_,psOUUZ4SamInputs_,
                   &psOUUSamOutputs_[ii*psOUUZ4nSamples_], cverrors);
          printf("OUUOptimizer: Z3 sample %d (of %d)\n",ii+1,
                 psOUUZ3nSamples_);
          printf("RS CV avg error = %e (scaled) \n", cverrors[0]);
          printf("RS CV rms error = %e (scaled) \n", cverrors[1]);
          printf("RS CV max error = %e (scaled) \n", cverrors[2]);
        }
        psOUUfaPtr_->evaluatePoint(psOUULargeSampleSize_,
                          psOUULargeSamInputs_,psOUULargeSamOutputs_);
        if (isScreenDumpModeOn() && odata->outputLevel_ > 1) 
          printf("OUUOptimizer: RS evaluation completed\n");
        if (psOUUMode_ == 1 || psOUUMode_ == 2)
        {
          mean = 0.0;
          for (ss = 0; ss < psOUULargeSampleSize_; ss++) 
            mean += psOUULargeSamOutputs_[ss];
          mean /= (double) psOUULargeSampleSize_;
          resultStore[ii] = mean;
        }
        if (psOUUMode_ == 2 && psOUUStdevMultiplier_ != 0.0)
        {
          stdev = 0.0;
          for (ss = 0; ss < psOUULargeSampleSize_; ss++) 
            stdev += pow(psOUULargeSamOutputs_[ss] - mean, 2.0);
          stdev /= (double) psOUULargeSampleSize_;
          resultStore[ii] = mean + psOUUStdevMultiplier_ * sqrt(stdev);
        }
        if (psOUUMode_ == 3)
        {
          sortDbleList(psOUULargeSampleSize_, psOUULargeSamOutputs_);
          kk = (int) (psOUUPercentile_ * psOUULargeSampleSize_);
          if (kk >= psOUULargeSampleSize_) 
            kk = psOUULargeSampleSize_ - 1;
          resultStore[ii] = psOUULargeSamOutputs_[kk];
        }
      }
      if (psOUUMode_ == 1 || psOUUMode_ == 2)
      {
        mean = 0.0;
        for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
        {
          if (psOUUSamProbs_ == NULL)
            mean += resultStore[ii] / (double) psOUUZ3nSamples_;
          else
            mean += resultStore[ii] * psOUUSamProbs_[ii];
        }
        (*YValue) = mean;
      }
      if (psOUUMode_ == 3)
      {
        mean = 0.0;
        for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
          mean += resultStore[ii] / (double) psOUUZ3nSamples_;
        (*YValue) = mean;
      }
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1) 
      {
        printf("OUUOptimizer: computed  objective (with RS) = %e.\n",
               (*YValue));
      }
      free(resultStore);
    }

    if ((*YValue) < odata->optimalY_)
    {
      odata->optimalY_ = (*YValue);
      for (ii = 0; ii < M; ii++) 
        odata->optimalX_[ii] = psOUUOptimalX_[ii];
      if (isScreenDumpModeOn() && odata->outputLevel_ > 0) 
      {
        printf("    OUUOptimizer outer loop new Ymin = %16.8e (***)\n",
               (*YValue));
        if (psOUUUserOpt_ == 0)
        {
          for (ii = 0; ii < M1+M2; ii++) 
            printf("        Input %3d at min = %e\n", ii+1, 
                   odata->optimalX_[ii]);
        }
      }
    }
    psOUUCounter_++;

    free(readys);
    free(XLocal);
    return NULL;
  }

  int ouuevalfunccobyla(int nInputs, int nConstraints, double *XValues,
                        double *YValue, double *constraints, void *odata2)
  {
    int    ii, kk, ss, funcID, M, M1, M2, M3, M4, NC;
    int    maxfun, *readys, status, index, nSamp, nOuts;
    double rhobeg, rhoend, ddata, *workArray, mean, stdev, *XLocal;
    double *lowers, *uppers, *outYs, *tempY;
    char   winput[1000];
    oData  *odata;
    FILE   *fp=NULL;

    odata = (oData *) psOUUObj_;
    nOuts = odata->nOutputs_;
    M1    = psOUUM1_;
    M2    = psOUUM2_;
    M3    = psOUUM3_;
    M4    = psOUUM4_;
    M     = M1 + M2 + M3 + M4;
    
    fp = fopen("psuade_ouu_stop","r");
    if (fp != NULL)
    {
      unlink("psuade_ouu_stop");
      printf("OUUOptimizer: psuade_ouu_stop file found.\n");
      printf("              Abrupt termination.\n");
      if (psOUUNSaved_ > 0)
      {
        fp = fopen("psuade_ouu_history","w");
        if (fp != NULL)
        {
          for (ii = 0; ii < psOUUNSaved_; ii++)
          {
            fprintf(fp, "999 %d ", M);
            for (kk = 0; kk < M; kk++)
              fprintf(fp, "%24.16e ", psOUUSaveX_[ii*M+kk]);
            for (kk = 0; kk < nOuts; kk++)
              fprintf(fp, "%24.16e\n",
                      psOUUSaveY_[ii*nOuts+kk]);
          }
          fclose(fp);
        }
        printf("OUUOptimizer: history saved in psuade_ouu_history\n");
      }
      exit(1);
    }

    nSamp = psOUUZ3nSamples_ * psOUUZ4nSamples_;
    XLocal = (double *) malloc(M*sizeof(double));
    readys = (int *) malloc(nSamp*sizeof(int));
    if (psOUUOptCode_ == 2) NC = nConstraints - 2 * M1;
    else                    NC = nConstraints;

    index = 0;
    for (ii = 0; ii < M; ii++)
    {
      if (psOUUInputTypes_[ii] == psOUUType1)
      {
        psOUUXValues_[ii] = XValues[index];
        index++;
      }
    }

    if (isScreenDumpModeOn())
    {
      printf("OUUOptimizer: Outer optimization iteration = %d\n",
             odata->numFuncEvals_/nSamp+1);
      if (odata->outputLevel_ > 1)
      {
        printf("OUUOptimizer: Outer optimization loop FuncEval %d.\n",
               odata->numFuncEvals_+1);
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType1)
          {
            printf("    Current Level 1 input %3d = %e\n", ii+1,
                   XValues[index]);
            index++;
          }
        }
      }
    }

    if (psOUUEnsembleEval_ == 0)
    {
      outYs = new double[nOuts];
      for (ss = 0; ss < nSamp; ss++)
      {
        if (isScreenDumpModeOn() && odata->outputLevel_ > 3)
           printf("OUUOptimizer sample %d (of %d)\n",ss+1,nSamp);
        readys[ss] = -1;
        psOUUOptimalY_ = PSUADE_UNDEFINED;
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType3)
          {
            kk = ss / psOUUZ4nSamples_;
            psOUUWValues_[index] = psOUUZ3SamInputs_[kk*M3+index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType4)
          {
            kk = ss % psOUUZ4nSamples_;
            psOUUWValues_[index] = psOUUZ4SamInputs_[kk*M4+index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType1)
          {
            XLocal[ii] = psOUUXValues_[index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType2)
          {
            //Nov 2017 - set to user-specified values
            //XLocal[ii] = 0.5 * (odata->lowerBounds_[ii] +
            //                    odata->upperBounds_[ii]);
            XLocal[ii] = psOUUM2Values_[index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType3)
          {
            kk = ss / psOUUZ4nSamples_;
            XLocal[ii] = psOUUZ3SamInputs_[kk*M3+index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType4)
          {
             kk = ss % psOUUZ4nSamples_;
             XLocal[ii] = psOUUZ4SamInputs_[kk*M4+index];
             index++;
          }
        }
        int found = 0;
        for (kk = 0; kk < psOUUNSaved_; kk++)
        {
          for (ii = 0; ii < M; ii++)
          {
            if (psOUUInputTypes_[ii] != psOUUType2 &&
                (PABS(psOUUSaveX_[kk*M+ii]-XLocal[ii])>1.0e-14))
              break;
          }
          if (ii == M)
          {
            found = 1;
            if (isScreenDumpModeOn() && odata->outputLevel_ > 2)
              printf("OUUOptimizer: simulation results reuse.\n");
            for (ii = 0; ii < M; ii++)
            {
              if (psOUUInputTypes_[ii] == psOUUType2)
                XLocal[ii] = psOUUSaveX_[kk*M+ii];
            }
            for (ii = 0; ii < nOuts; ii++)
              psOUUSamOutputs_[ss*nOuts+ii] = psOUUSaveY_[kk*nOuts+ii];
            for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = XLocal[ii];
            readys[ss] = 0;
            break;
          }
        }
        if (found == 0)
        {
          funcID = psOUUCounter_ * nSamp + ss;
          if (odata->optFunction_ != NULL)
          {
            odata->optFunction_(M, XLocal, nOuts, outYs);
            readys[ss] = 0;
          }
          else if (odata->funcIO_ != NULL)
          {
            readys[ss] = odata->funcIO_->evaluate(funcID,M,XLocal,
                                                  nOuts,outYs,0);
          }
          else
          {
            printf("OUUOptimizer ERROR: no function evaluator set.\n");
            exit(1);
          }
          for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = XLocal[ii];
          for (ii = 0; ii < nOuts; ii++)
            psOUUSamOutputs_[ss*nOuts+ii] = outYs[ii];
          if (psOUUSaveHistory_ == 1 && 
              (psOUUNSaved_+1)*M < psOUUMaxSaved_*10 &&
              (psOUUNSaved_+1)*nOuts < psOUUMaxSaved_*10)
          {
            for (ii = 0; ii < M; ii++)
              psOUUSaveX_[psOUUNSaved_*M+ii] = XLocal[ii];
            for (ii = 0; ii < nOuts; ii++)
              psOUUSaveY_[psOUUNSaved_*nOuts+ii] =
                        psOUUSamOutputs_[ss*nOuts+ii];
            psOUUNSaved_++;
          }
          if (psOUUParallel_ == 0) odata->numFuncEvals_++;
        }
      }

      if (psOUUUserOpt_ == 1 && psOUUParallel_ == 1)
      {
        for (ss = 0; ss < nSamp; ss++)
        {
          funcID = psOUUCounter_ * nSamp + ss;
          if (readys[ss] != 0)
          {
            while (readys[ss] != 0)
            {
#ifdef WINDOWS
              Sleep(1000);
#else
              sleep(1);
#endif
              readys[ss] = odata->funcIO_->evaluate(funcID,M,XLocal,
                                                    nOuts,outYs,2);
            }
            for (ii = 0; ii < nOuts; ii++)
              psOUUSamOutputs_[ss*nOuts+ii] = outYs[ii];
            odata->numFuncEvals_++;
          }
        }
      }
      delete [] outYs;
    }
    else
    {
      index = 0;
      for (ii = 0; ii < M; ii++)
      {
        if (psOUUInputTypes_[ii] == psOUUType1)
        {
          XLocal[index] = psOUUXValues_[ii];
          index++;
        }
      }
      for (ss = 0; ss < nSamp; ss++)
      {
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType1)
          {
            psOUUXValues_[ss*M+ii] = XLocal[index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType2)
          {
            //Nov 2017 - set to user-specified values
            //psOUUXValues_[ss*M+ii] = 0.5*(odata->lowerBounds_[ii] +
            //                              odata->upperBounds_[ii]);
            psOUUXValues_[ss*M+ii] = psOUUM2Values_[index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType3)
          {
            kk = ss / psOUUZ4nSamples_;
            psOUUXValues_[ss*M+ii] = psOUUZ3SamInputs_[kk*M3+index];
            index++;
          }
        }
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType4)
          {
            kk = ss % psOUUZ4nSamples_;
            psOUUXValues_[ss*M+ii] = psOUUZ4SamInputs_[kk*M4+index];
            index++;
          }
        }
      }
      funcID = odata->numFuncEvals_;
      odata->funcIO_->ensembleEvaluate(nSamp,M,psOUUXValues_,
                                       nOuts,psOUUSamOutputs_,funcID);
      odata->numFuncEvals_ += nSamp;
      for (ss = 0; ss < nSamp; ss++)
      {
        if (psOUUSaveHistory_ == 1 && 
            (psOUUNSaved_+1)*M < psOUUMaxSaved_*10 &&
            (psOUUNSaved_+1)*nOuts < psOUUMaxSaved_*10)
        {
          for (kk = 0; kk < M; kk++)
            psOUUSaveX_[psOUUNSaved_*M+kk] = psOUUXValues_[ss*M+kk];
          for (kk = 0; kk < nOuts; kk++)
            psOUUSaveY_[psOUUNSaved_*nOuts+kk] =
                psOUUSamOutputs_[ss*nOuts+kk];
          psOUUNSaved_++;
        }
      }
      for (ii = 0; ii < M; ii++) psOUUOptimalX_[ii] = psOUUXValues_[ii];
    }

    int failCnt=0;
    for (ss = 0; ss < nSamp; ss++)
       if (psOUUSamOutputs_[ss*nOuts] >= 0.98*PSUADE_UNDEFINED) failCnt++;
    if (failCnt != 0)
       printf("WARNING: there are %d failed runs out of %d\n",failCnt,nSamp);
    if (failCnt == nSamp || (failCnt > 0 && psOUUSamProbs_ != NULL))
    {
       printf("ERROR: X3 sample cannot admit any failures ==> terminate.\n");
       exit(1);
    }
    if (failCnt > 0 && psOUUMode_ > 2)
    {
      printf("ERROR: OUU mode %d cannot admit any failures ==> terminate.\n",
             psOUUMode_);
      exit(1);
    }
    for (ss = 0; ss < nConstraints; ss++) constraints[ss] = 0.0;
    if (psOUUUseRS_ == 0 || nSamp == 1)
    {
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1)
        printf("OUUOpt: computing objective (no RS), nFuncEval = %d\n",
               odata->numFuncEvals_);
      if (psOUUMode_ == 1 || psOUUMode_ == 2)
      {
        mean = 0.0;
        if (psOUUZ3nSamples_ == 1)
        {
          for (ss = 0; ss < nSamp; ss++)
            if (psOUUSamOutputs_[ss*nOuts] < 0.98*PSUADE_UNDEFINED)
              mean += psOUUSamOutputs_[ss*nOuts]/(1.0*nSamp-failCnt);
        }
        else
        {
          for (ss = 0; ss < nSamp; ss++)
          {
             index = ss / psOUUZ4nSamples_;
             mean += psOUUSamOutputs_[ss*nOuts] / psOUUZ4nSamples_ *
                     psOUUSamProbs_[index];
          }
        }
        (*YValue) = mean;
      }
      if (psOUUMode_ == 2 && psOUUStdevMultiplier_ != 0.0)
      {
        stdev = 0.0;
        for (ss = 0; ss < nSamp; ss++)
        {
           index = ss / psOUUZ4nSamples_;
          if (psOUUZ3nSamples_ == 1)
          {
            if (psOUUSamOutputs_[ss*nOuts] < 0.98*PSUADE_UNDEFINED)
              stdev += pow(psOUUSamOutputs_[ss*nOuts]-mean,2.0)/
                       (double) (nSamp - failCnt);
          }
          else
          {
            stdev += pow(psOUUSamOutputs_[ss*nOuts]-mean, 2.0)*
                     psOUUSamProbs_[index] / (double) psOUUZ4nSamples_;
          }
        }
        (*YValue) = mean + psOUUStdevMultiplier_ * sqrt(stdev);
      }
      if (psOUUMode_ == 3)
      {
        tempY = new double[nSamp];
        for (ss = 0; ss < nSamp; ss++)
          tempY[ss] = psOUUSamOutputs_[ss*nOuts];
        sortDbleList(nSamp, tempY);
        kk = (int) (psOUUPercentile_ * nSamp);
        if (kk >= nSamp) kk = nSamp - 1;
        (*YValue) = tempY[kk];
        delete [] tempY;
      }
      if (psOUUCMode_ == 1)
      {
        for (ii = 0; ii < NC; ii++)
        {
          constraints[ii] = 0;
          for (ss = 0; ss < nSamp; ss++)
            constraints[ii] += psOUUSamOutputs_[ss*nOuts+ii+1];
          constraints[ii] /= (double) nSamp;
        }
      }
      else if (psOUUCMode_ == 2)
      {
        for (ss = 0; ss < nSamp; ss++)
        {
          ddata = 0.0;
          for (ii = 0; ii < NC; ii++)
            if (psOUUSamOutputs_[ss*nOuts+ii+1] < 0)
              ddata += pow(psOUUSamOutputs_[ss*nOuts+ii+1],2.0);
          psOUUSamConstrNorms_[ss] = sqrt(ddata);
        }
        ddata = psOUUSamConstrNorms_[0];
        for (ss = 1; ss < nSamp; ss++)
        {
          if (psOUUSamConstrNorms_[ss] > ddata)
          {
            ddata = psOUUSamConstrNorms_[ss];
            kk = ss;
          }
        }
        for (ii = 0; ii < NC; ii++)
          constraints[ii] = psOUUSamOutputs_[kk*nOuts+ii+1];
      }
      else if (psOUUCMode_ == 3)
      {
        for (ss = 0; ss < nSamp; ss++)
        {
          ddata = 0.0;
          for (ii = 0; ii < NC; ii++)
            if (psOUUSamOutputs_[ss*nOuts+ii+1] < 0)
              ddata += psOUUSamOutputs_[ss*nOuts+ii+1];
          psOUUSamConstrNorms_[ss] = - ddata;
        }
        ddata = psOUUSamConstrNorms_[0];
        for (ss = 1; ss < nSamp; ss++)
        {
          if (psOUUSamConstrNorms_[ss] > ddata)
          {
            ddata = psOUUSamConstrNorms_[ss];
            kk = ss;
          }
        }
        for (ii = 0; ii < NC; ii++)
          constraints[ii] = psOUUSamOutputs_[kk*nOuts+ii+1];
      }
      else if (psOUUCMode_ == 4)
      {
        for (ss = 0; ss < nSamp; ss++)
        {
          ddata = 0;
          for (ii = 0; ii < NC; ii++)
            if (psOUUSamOutputs_[ss*nOuts+ii+1] < 0 &&
                psOUUSamOutputs_[ss*nOuts+ii+1] < ddata)
              ddata = psOUUSamOutputs_[ss*nOuts+ii+1];
          psOUUSamConstrNorms_[ss] = - ddata;
        }
        ddata = psOUUSamConstrNorms_[0];
        for (ss = 1; ss < nSamp; ss++)
        {
          if (psOUUSamConstrNorms_[ss] > ddata)
          {
            ddata = psOUUSamConstrNorms_[ss];
            kk = ss;
          }
        }
        for (ii = 0; ii < NC; ii++)
          constraints[ii] = psOUUSamOutputs_[kk*nOuts+ii+1];
      }
      if (psOUUOptCode_ == 2)
      {
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType1)
          {
            constraints[NC+2*index] = psOUULBounds_[index]-XValues[index];
            constraints[NC+2*index+1] = XValues[index]-
                                        psOUUUBounds_[index];
            index++;
          }
        }
      }
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1)
      {
        status = 1;
        for (ii = 0; ii < nConstraints; ii++)
          if (constraints[ii] > 0) status = 0;
        if (status == 1)
          printf("OUUOptimizer: computed  objective (no RS) = %e.\n",
                 (*YValue));
        else
          printf("OUUOptimizer: infeasible objective = %e.\n",(*YValue));
      }
      if (psOUUOptCode_ == 2)
      {
        for (ii = 0; ii < nConstraints; ii++)
          constraints[ii] = - constraints[ii];
      }
    }
    else
    {
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1)
         printf("OUUOpt: computing objective with RS, nFuncEval = %d\n",
                odata->numFuncEvals_);

      double totalMean=0.0, totalStdv=0.0, *resultStore, cverrors[3];
      resultStore=(double *) malloc(psOUUZ3nSamples_*sizeof(double));
      for (ii = 0; ii < psOUUZ3nSamples_; ii++)
      {
        status = psOUUfaPtr_->initialize(psOUUZ4SamInputs_,
                        &psOUUSamOutputs_[(ii*psOUUZ4nSamples_)*nOuts]);
        psOUUfaPtr_->evaluatePoint(psOUULargeSampleSize_,
                           psOUULargeSamInputs_,psOUULargeSamOutputs_);
        if (psOUUMode_ == 1 || psOUUMode_ == 2)
        {
          mean = 0.0;
          for (ss = 0; ss < psOUULargeSampleSize_; ss++)
            mean += psOUULargeSamOutputs_[ss*nOuts];
          mean /= (double) psOUULargeSampleSize_;
          resultStore[ii] = mean;
        }
        if (psOUUMode_ == 2 && psOUUStdevMultiplier_ != 0.0)
        {
          stdev = 0.0;
          for (ss = 0; ss < psOUULargeSampleSize_; ss++)
            stdev += pow(psOUULargeSamOutputs_[ss*nOuts]-mean, 2.0);
          stdev /= (double) psOUULargeSampleSize_;
          resultStore[ii] = mean + psOUUStdevMultiplier_ * stdev;
        }
        if (psOUUMode_ == 3)
        {
          tempY = new double[psOUULargeSampleSize_];
          for (ss = 0; ss < psOUULargeSampleSize_; ss++)
            tempY[ss] = psOUUSamOutputs_[ss*nOuts];
          sortDbleList(psOUULargeSampleSize_, tempY);
          kk = (int) (psOUUPercentile_ * psOUULargeSampleSize_);
          if (kk >= psOUULargeSampleSize_)
            kk = psOUULargeSampleSize_ - 1;
          resultStore[ii] = tempY[kk];
          delete [] tempY;
        }
      }
      if (psOUUMode_ == 1 || psOUUMode_ == 2)
      {
        mean = 0.0;
        for (ii = 0; ii < psOUUZ3nSamples_; ii++)
        {
          if (psOUUSamProbs_ == NULL)
            mean += resultStore[ii] / (double) psOUUZ3nSamples_;
          else
            mean += resultStore[ii] * psOUUSamProbs_[ii];
        }
        (*YValue) = mean;
      }
      if (psOUUMode_ == 3)
      {
        mean = 0.0;
        for (ii = 0; ii < psOUUZ3nSamples_; ii++)
          mean += resultStore[ii] / (double) psOUUZ3nSamples_;
        (*YValue) = mean;
      }
      if (isScreenDumpModeOn() && odata->outputLevel_ > 1)
      {
        status = 1;
        for (ii = 0; ii < nConstraints; ii++)
          if (constraints[ii] > 0) status = 0;
        if (status == 1)
          printf("OUUOptimizer: computed  objective (with RS) = %e.\n",
                  (*YValue));
        else
          printf("OUUOptimizer: infeasible objective = %e.\n",(*YValue));
      }
      if (psOUUOptCode_ == 2)
      {
        index = 0;
        for (ii = 0; ii < M; ii++)
        {
          if (psOUUInputTypes_[ii] == psOUUType1)
          {
            constraints[NC+2*index] = XValues[index]-psOUULBounds_[index];
            constraints[NC+2*index+1] = psOUUUBounds_[index] -
                                        XValues[index];
            index++;
          }
        }
        for (ii = 0; ii < NC; ii++)
          constraints[ii] = - constraints[ii];
      }
      free(resultStore);
    }

    if ((*YValue) < odata->optimalY_)
    {
      status = 1;
      if (psOUUOptCode_ == 2)
      {
        for (ii = 0; ii < nConstraints; ii++)
          if (constraints[ii] < 0) status = 0;
      }
      else
      {
        for (ii = 0; ii < nConstraints; ii++)
          if (constraints[ii] > 0) status = 0;
      }
      if (status == 1)
      {
        odata->optimalY_ = (*YValue);
        for (ii = 0; ii < M; ii++)
          odata->optimalX_[ii] = psOUUOptimalX_[ii];
        if (isScreenDumpModeOn() && odata->outputLevel_ > 0)
        {
          printf("    OUUOptimizer outer loop new Ymin = %16.8e (***)\n",
                 (*YValue));
          if (psOUUUserOpt_ == 0)
          {
            for (ii = 0; ii < M1+M2; ii++)
              printf("        Input %3d at min = %e\n", ii+1,
                     odata->optimalX_[ii]);
          }
        }
      }
    }
    psOUUCounter_++;

    free(readys);
    free(XLocal);
    return 0;
  }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
OUUOptimizer::OUUOptimizer()
{
  psOUUM1_ = -1;
  psOUUM2_ = -1;
  psOUUM3_ = -1;
  psOUUZ3nSamples_ = -1;
  psOUUZ4nSamples_ = -1;
  psOUUUserOpt_  = 0;
  psOUUObj_ = NULL;
  psOUUZ3SamInputs_ = NULL;
  psOUUZ4SamInputs_ = NULL;
  psOUUZ4LBounds_ = NULL;
  psOUUZ4UBounds_ = NULL;
  psOUUSamOutputs_ = NULL;
  psOUUSamConstrNorms_ = NULL;
  psOUUXValues_ = NULL;
  psOUUWValues_ = NULL;
  psOUUOptimalX_ = NULL;
  psOUUOptimalY_ = 0.0;
  psOUUUseRS_ = 0;
  optCode_ = -1;
  repeatFlag_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
OUUOptimizer::~OUUOptimizer()
{
  if (psOUUZ3SamInputs_ != NULL) delete [] psOUUZ3SamInputs_;
  if (psOUUZ4SamInputs_ != NULL) delete [] psOUUZ4SamInputs_;
  if (psOUUZ4LBounds_   != NULL) delete [] psOUUZ4LBounds_;
  if (psOUUZ4UBounds_   != NULL) delete [] psOUUZ4UBounds_;
  if (psOUUSamOutputs_  != NULL) delete [] psOUUSamOutputs_;
  if (psOUUSamProbs_    != NULL) delete [] psOUUSamProbs_;
  if (psOUUSamConstrNorms_  != NULL) delete [] psOUUSamConstrNorms_;
  if (psOUULargeSamInputs_  != NULL) delete [] psOUULargeSamInputs_;
  if (psOUULargeSamOutputs_ != NULL) delete [] psOUULargeSamOutputs_;
  if (psOUUXValues_  != NULL) delete [] psOUUXValues_;
  if (psOUUWValues_  != NULL) delete [] psOUUWValues_;
  if (psOUUOptimalX_ != NULL) delete [] psOUUOptimalX_;
  if (psOUUfaPtr_ != NULL) delete psOUUfaPtr_;
  if (psOUUInputTypes_ != NULL) delete psOUUInputTypes_;
  if (psOUUDesignTypes_ != NULL) delete psOUUDesignTypes_;
  if (psOUUM2Values_ != NULL) delete psOUUM2Values_;
  psOUUZ3SamInputs_ = NULL;
  psOUUZ4SamInputs_ = NULL;
  psOUUZ4LBounds_ = NULL;
  psOUUZ4UBounds_ = NULL;
  psOUUSamOutputs_ = NULL;
  psOUUSamConstrNorms_ = NULL;
  psOUUSamProbs_ = NULL;
  psOUULargeSamInputs_ = NULL;
  psOUULargeSamOutputs_ = NULL;
  psOUUXValues_ = NULL;
  psOUUWValues_ = NULL;
  psOUUOptimalX_ = NULL;
  psOUUfaPtr_ = NULL;
  psOUUInputTypes_ = NULL;
  psOUUDesignTypes_ = NULL;
  psOUUM2Values_ = NULL;
}

// ************************************************************************
// optimize (this function should be used in the library mode)
// ------------------------------------------------------------------------
void OUUOptimizer::optimize(int nInputs, double *XValues, double *lbds,
                      double *ubds, int nOutputs, int maxfun, double tol)
{
   double *optimalX;
   if (nInputs <= 0)
   {
      printf("OUUOptimizer ERROR: nInputs <= 0.\n");
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
   psOUUnOutputs_ = nOutputs;
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
   optimalX_ = new double[nInputs];
   for (int ii = 0; ii < nInputs; ii++)
      XValues[ii] = optimalX_[ii] = odata->optimalX_[ii];
   optimalY_ = odata->optimalY_;
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void OUUOptimizer::optimize(oData *odata)
{
   int    nInputs, printLevel=0, ii, kk, maxfun, nPts=0, bobyqaFlag=1111;
   int    M1, M2, M3, M4, index, iOne=1, iZero=0, hasColNums=0;
   int    count, *inputPDFs, nOutputs, itmp;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, ddata, *workArray, *tdata;
   char   pString[1000], *cString, lineIn[20001], filename[1000];
   FILE   *fp=NULL;
   pData  pdata;
   static int currDriver=-1;

   nInputs  = odata->nInputs_;
   nOutputs = odata->nOutputs_;
   printLevel = odata->outputLevel_;
   kk = psOUUM1_ + psOUUM2_ + psOUUM3_ + psOUUM4_;

   if (isLibraryModeOn())
   {
     if (nInputs != kk)
     {
       printf("OUUOptimizer ERROR: OUU set up incomplete.\n");
       printf("             Did you use setVariableType but you did \n");
       printf("                 set all variables?\n");
       printf("     current: M1 = %d (-1 means undefined) \n", psOUUM1_);
       printf("     current: M2 = %d\n", psOUUM2_);
       printf("     current: M3 = %d\n", psOUUM3_);
       printf("     current: M4 = %d\n", psOUUM4_);
       exit(1);
     }
     if (psOUUM1_ <= 0)
     {
       printf("OUUOptimizer ERROR: OUU set up invalid.\n");
       printf("             The number of design variables has to be >0\n");
       exit(1);
     } 
   } 
       
   if (repeatFlag_ == 1 && (nOutputs != psOUUnOutputs_ || kk != nInputs))
   {
     printf("OUUOptimizer WARNING: reuse but different information.\n");
     printf("             INFO: everything will be reset.\n");
     repeatFlag_ = 0;
     if (psOUUZ3SamInputs_ != NULL) delete [] psOUUZ3SamInputs_;
     if (psOUUZ4SamInputs_ != NULL) delete [] psOUUZ4SamInputs_;
     if (psOUUZ4LBounds_   != NULL) delete [] psOUUZ4LBounds_;
     if (psOUUZ4UBounds_   != NULL) delete [] psOUUZ4UBounds_;
     if (psOUUSamOutputs_  != NULL) delete [] psOUUSamOutputs_;
     if (psOUUSamProbs_    != NULL) delete [] psOUUSamProbs_;
     if (psOUUSamConstrNorms_  != NULL) delete [] psOUUSamConstrNorms_;
     if (psOUULargeSamInputs_  != NULL) delete [] psOUULargeSamInputs_;
     if (psOUULargeSamOutputs_ != NULL) delete [] psOUULargeSamOutputs_;
     if (psOUUXValues_  != NULL) delete [] psOUUXValues_;
     if (psOUUWValues_  != NULL) delete [] psOUUWValues_;
     if (psOUUOptimalX_ != NULL) delete [] psOUUOptimalX_;
     if (psOUUfaPtr_ != NULL) delete psOUUfaPtr_;
     if (psOUUDesignTypes_ != NULL) delete psOUUDesignTypes_;
     if (!isLibraryModeOn())
     {
       if (psOUUInputTypes_ != NULL) delete psOUUInputTypes_;
       psOUUInputTypes_ = NULL;
     }
     psOUUZ3SamInputs_ = NULL;
     psOUUZ4SamInputs_ = NULL;
     psOUUZ4LBounds_ = NULL;
     psOUUZ4UBounds_ = NULL;
     psOUUSamOutputs_ = NULL;
     psOUUSamConstrNorms_ = NULL;
     psOUUSamProbs_ = NULL;
     psOUULargeSamInputs_ = NULL;
     psOUULargeSamOutputs_ = NULL;
     psOUUXValues_ = NULL;
     psOUUWValues_ = NULL;
     psOUUOptimalX_ = NULL;
     psOUUfaPtr_ = NULL;
     psOUUDesignTypes_ = NULL;
     psOUUEnsembleEval_ = 0;
     psOUUParallel_ = 0;
     psOUUUseRS_ = 0;
   } 

   if (isScreenDumpModeOn() && repeatFlag_ == 0 && printLevel >= 0)
   {
     printAsterisks(PL_INFO, 0);
     printAsterisks(PL_INFO, 0);
     printf("*     1- OR 2-STAGE OPTIMIZATION UNDER UNCERTAINTY (OUU)\n");
     printEquals(PL_INFO, 0);
     printf("This optimization capability solves the following problem:\n");
     printf("\n   minimize_Z1 { Phi_{Z3,Z4} [ G(Z1,Z2,Z3,Z4) ] } \n\n");
     printf("   subject to either:\n");
     printf("    (a) bound constraints on Z1, Z2, and Z4 (use bobyqa); or \n");
     printf("    (b) no constraint on Z1, Z2, and Z4 (use newuoa); or\n");
     printf("    (c) inequality constraints on Z1,Z2,Z3,Z4 (use cobyla)\n");
     printf("    as the level-1 optimizer.)\n\n");
     printf("   Z3 is a set of discrete parameters for which a sample\n\n");
     printf("   is to be provided by the user. \n\n");
     printf("   (0) How to perform regular optimization? \n");
     printf("       In this case \n"); 
     printf("       - Z1 will be the optimization variables\n");
     printf("       - Z2 should be an empty set\n");
     printf("       - Z3 should be an empty set\n");
     printf("       - Z4 should be an empty set\n");
     printf("       - G(Z1,Z3,Z4)=G(Z1) is the simulator (opt_driver)\n");
     printf("   (1) How to perform 1-level OUU? \n");
     printf("       In this case \n"); 
     printf("       - Z1: the optimization variables\n");
     printf("       - Z2: should be an empty set\n");
     printf("       - Z3: parameters that you will provide a sample.\n");
     printf("       - Z4: parameters that you do not provide a sample.\n");
     printf("         (optionally, you can choose sampling scheme and\n");
     printf("          whether you want response surface for Z4).\n");
     printf("       - G(Z1,Z3,Z4) is the simulator (opt_driver)\n");
     printf("\n");
     printf("   (2) How to perform 2-level OUU? \n");
     printf("       In this case \n"); 
     printf("       - Z1: the optimization variables\n");
     printf("       - Z2: the set of level-2 optimization parameters.\n");
     printf("         If Z2 do not need to be published outside of the G\n");
     printf("         function, then it can be specified as an empty set.\n");
     printf("       - Z3: parameters that you will provide a sample for.\n");
     printf("       - Z4: parameters that you do not provide a sample.\n");
     printf("       - There are 2 options to how to set up G(Z1,Z2,Z3,Z4):\n");
     printf("         (a) A user-provided level 2 optimizer that does\n\n");
     printf("              G(Z1,Z2,Z3,Z4) = minimize_Z2{F(Z1,Z2,Z3,Z4)}\n\n");
     printf("             where F(Z1,Z2,Z3,Z4) is embedded in G(Z1,...)\n\n");
     printf("             In this case, the user will provide G(Z1,...)\n");
     printf("             via opt_driver (i.e. opt_driver = gfunction)\n");
     printf("         (b) A user-provided function F(Z1,...) such that\n\n");
     printf("              G(Z1,Z2,Z3,Z4) = minimize_Z2 {F(Z1,Z2,Z3,Z4)}\n\n");
     printf("             In this case, user is expected to provide the\n");
     printf("             F(Z1,...) function via 'opt_driver = ffunction'\n");
     printf("             and OUU provides the optimizer (BOBYQA).\n\n");
     printf("   In case 1 and 2, Phi_{Z3,Z4} is a functional on G(Z1,...)\n");
     printf("   with respect to Z3 and Z4, e.g. Phi_{Z3,Z4} may be:\n");
     printf("   1. mean of G(Z1,Z2,Z3,Z4) with respect to Z3,Z4 (default)\n"); 
     printf("   2. mean of G(Z1,Z2,Z3,Z4) + alpha * std dev of G(Z1,...)\n"); 
     printf("   3. G(Z1,Z2,Z3*,Z4*) such that \n");
     printf("         Prob(G(Z1,Z2,Z3,Z4)>G(Z1,Z2,Z3*,Z4*)) = epsilon\n");
     printf("   4. min_{Z3,Z4} G(Z1,Z2,Z3,Z4) given Z1,Z2 (robust opt.)\n");
     printEquals(PL_INFO, 0);
     printf("In the above formulation, the total no. of parameters M = %d\n", 
            nInputs);
     printf("These parameters are to be divided into four groups:\n");
     printf("(1) Stage 1 optimization parameters Z1 (M1 >= 1) \n");
     printf("(2) Stage 2 optimization (recourse) parameters Z2\n");
     printf("(3) uncertain parameters Z3 (with a user-provided sample) \n");
     printf("(4) uncertain parameters Z4 \n");
     printf("    - continuous parameters (PSUADE to generate sample), or\n");
     printf("    - a large sample (created by PSUADE from RS) \n");
     printf("Thus, the first M1+M2 parameters are considered to be\n");
     printf("optimization parameters, and M3+M4 are uncertain parameters\n");
     printf("so that M = M1 + M2 + M3 + M4.\n");
     printEquals(PL_INFO, 0);
     printf("To reuse simulation results (e.g. from before abrupt\n");
     printf("termination), turn on save_history and use_history\n");
     printf("optimization options in the ANALYSIS section (e.g.\n");
     printf("optimization save_history). You will see a text file\n");
     printf("called 'psuade_ouu_history' afterward.\n");
     printAsterisks(PL_INFO, 0);
     printf("IF YOU ARE READY TO MOVE ON, ENTER 'y' AND RETURN : ");
     lineIn[0] = '0';
     while (lineIn[0] != 'y' && lineIn[0] != 'n')
     {
       scanf("%s", lineIn);
       if (lineIn[0] == 'n') 
       {
         odata->optimalY_ = 0;
         for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
         printf("OUUOptimizer INFO: abort.\n");
         exit(1);
       }
     }
     printEquals(PL_INFO, 0);
   }

   if (psConfig_ != NULL)
   {
     kk = 0;
     for (ii = 0; ii < nInputs; ii++)
     {
       sprintf(pString,"iDiscrete%d", ii+1);
       cString = psConfig_->getParameter(pString);
       if (cString != NULL) kk++;
     }
     if (kk > 0)
     {
       if (psOUUDesignTypes_ != NULL) delete [] psOUUDesignTypes_;
       psOUUDesignTypes_ = new int[nInputs];
       for (ii = 0; ii < nInputs; ii++)
       {
         sprintf(pString,"iDiscrete%d", ii+1);
         cString = psConfig_->getParameter(pString);
         if (cString != NULL)
         {
           psOUUDesignTypes_[ii] = 1;
           if (printLevel > 0)
             printf("OUU input %4d is discrete\n",ii+1);
         }
         else psOUUDesignTypes_[ii] = 0;
       }
       optCode_ = 5;
     }
     sprintf(pString,"ouu_sample_has_nums");
     cString = psConfig_->getParameter(pString);
     if (cString != NULL) hasColNums = 1;
   }

   if (nOutputs > 1 && repeatFlag_ == 0 && optCode_ < 2)
   {
     if (isScreenDumpModeOn())
        printf("*** OUUOptimizer detects nOutputs > 1.\n");
     if (nOutputs != nInputs+1)
     {
       if (isScreenDumpModeOn())
       {
         printf("    This will be interpreted to mean using cobyla to\n");
         printf("*   solve inequality-constraint problems (G(x) < 0).\n");
       }
       optCode_ = 2;
     }
     else
     {
       if (psMasterMode_ == 1)
       {
         printf("   This will be interpreted to mean either:\n");
         printf("   (1) inequality-constraints are used (use cobyla), or\n");
         printf("   (2) derivative information has been provided\n");
         kk = 0;
         while (kk < 1 || kk > 2)
         {
           printf("Please select which is correct: (1 or 2) ");
           scanf("%d", &kk);
         }
         if (kk == 1) optCode_ = 2;
         if (kk == 2) optCode_ = 3;
       } 
       else
       {
         if (isScreenDumpModeOn())
         {
           printf("   This will be interpreted to mean using cobyla to\n");
           printf("   solve inequality-constraint problems (G(x) < 0).\n");
         }
         optCode_ = 2;
       } 
     } 
   } 
   if (nOutputs > 1 && optCode_ < 2)
   {
     printOutTS(PL_ERROR,"OUU ERROR: incompabible choice of optimizer.\n");
     printOutTS(PL_ERROR,"           Due to nOutputs>1 and solver choice.\n");
     exit(1);
   }

   if (nOutputs == 1 && repeatFlag_ == 0 && optCode_ == -1) optCode_ = 0;

   if (psMasterMode_ != 0)
   {
     printf("OUUOptimizer INFO: Master mode to be turned off.\n");
     psMasterMode_ = 0;
     psOUUMasterMode_ = 1;
   }

   if ((psOUUM1_+psOUUM2_+psOUUM3_+psOUUM4_) == nInputs && (psOUUM1_>0))
   {
     M1 = psOUUM1_;
     M2 = psOUUM2_;
     M3 = psOUUM3_;
     M4 = psOUUM4_;
   }
   else
   {
     M1 = M2 = M3 = M4 = 0;
     printf("M1 = number of design (level 1 optim.) parameters\n");
     while (M1 <= 0 || M1 > nInputs)
     {
       printf("Enter M1 (between 1 and %d) : ", nInputs);
       scanf("%d", &M1);
     }
     if (M1 < nInputs)
     {
       //if (nOutputs == 1)
       //{
         printf("M2 = no. of recourse (level 2 optim.) parameters\n");
         M2 = -1;
         while (M2 < 0 || M1+M2 > nInputs)
         {
            printf("Enter M2 (between 0 and %d) : ", nInputs-M1);
            scanf("%d", &M2);
         }
       //}
       if ((M1 + M2) < nInputs)
       {
         printf("M3 = number of discrete (scenario) parameters\n");
         M3 = -1;
         while (M3 < 0 || M1+M2+M3 > nInputs)
         {
           printf("Enter M3 (between 0 and %d) : ", nInputs-M1-M2);
           scanf("%d", &M3);
         }
         M4 = nInputs - M1 - M2 - M3;
       }
     }
     psOUUM1_ = M1;
     psOUUM2_ = M2;
     psOUUM3_ = M3;
     psOUUM4_ = M4;
   }

   if (isScreenDumpModeOn())
   {
     if (M1 == nInputs && repeatFlag_ == 0)
     {
       psOUUInputTypes_ = new int[nInputs];
       for (ii = 0; ii < nInputs; ii++) psOUUInputTypes_[ii] = 1;
     }
     else if (psOUUInputTypes_ == NULL)
     {
       printf("In the following, please select type for each variable:\n");
       printf("  1. design variable (level 1 optimization parameter)\n");
       printf("  2. operating variable (level 2 optimization parameter)\n");
       printf("  3. discrete uncertain variable (a sample will be given)\n");
       printf("  4. continuous uncertain variable\n");
       printf("NOTE: make sure your specification matches with above.\n");
       fgets(lineIn, 5000, stdin);
       psOUUInputTypes_ = new int[nInputs];
       for (ii = 0; ii < nInputs; ii++)
       {
         sprintf(pString, "Type for variable %d ? ", ii+1);
         psOUUInputTypes_[ii] = getInt(1,4,pString); 
       }
     }
     if (repeatFlag_ == 0)
     {
       printEquals(PL_INFO, 0);
       for (ii = 0; ii < nInputs; ii++)
       {
         if (psOUUInputTypes_[ii] == psOUUType1) 
         printf("Input %4d is a design parameter.\n",ii+1);
         else if (psOUUInputTypes_[ii] == psOUUType2) 
           printf("Input %4d is an operating parameter.\n",ii+1);
         else if (psOUUInputTypes_[ii] == psOUUType3) 
           printf("Input %4d is a discrete uncertain parameter.\n",ii+1);
         else if (psOUUInputTypes_[ii] == psOUUType4) 
           printf("Input %4d is a continuous uncertain parameter.\n",ii+1);
       }
     }
   }
   int c1=0, c2=0, c3=0, c4=0;
   for (ii = 0; ii < nInputs; ii++)
   {
     if (psOUUInputTypes_[ii] == psOUUType1) c1++;
     if (psOUUInputTypes_[ii] == psOUUType2) c2++;
     if (psOUUInputTypes_[ii] == psOUUType3) c3++;
     if (psOUUInputTypes_[ii] == psOUUType4) c4++;
   }
   M1 = psOUUM1_;
   M2 = psOUUM2_;
   M3 = psOUUM3_;
   M4 = psOUUM4_;
   if (c1 != M1 || c2 != M2 || c3 != M3 || c4 != M4)
   {
     printf("OUUOptimizer ERROR: input type counts do not match.\n");
     printf("    Number of type 1 = %d (expected = %d)\n",c1,M1);
     printf("    Number of type 2 = %d (expected = %d)\n",c2,M2);
     printf("    Number of type 3 = %d (expected = %d)\n",c3,M3);
     printf("    Number of type 4 = %d (expected = %d)\n",c4,M4);
     exit(1);
   }
   if (M4 > 0 && repeatFlag_ == 0)
   {
     odata->psIO_->getParameter("input_pdfs", pdata);
     inputPDFs = pdata.intArray_;
     pdata.intArray_ = NULL;
     for (ii = 0; ii < nInputs; ii++)
     {
       if (psOUUInputTypes_[ii] == psOUUType4)
       {
         if (inputPDFs[ii] == PSUADE_PDF_UNIFORM) 
           strcpy(pString, "Uniform");            
         if (inputPDFs[ii] == PSUADE_PDF_NORMAL) 
           strcpy(pString, "Normal");            
         if (inputPDFs[ii] == PSUADE_PDF_LOGNORMAL) 
           strcpy(pString, "Lognormal");            
         if (inputPDFs[ii] == PSUADE_PDF_TRIANGLE) 
           strcpy(pString, "Triangle");            
         if (inputPDFs[ii] == PSUADE_PDF_BETA) 
           strcpy(pString, "Beta");            
         if (inputPDFs[ii] == PSUADE_PDF_WEIBULL) 
           strcpy(pString, "Weibull");            
         if (inputPDFs[ii] == PSUADE_PDF_GAMMA) 
           strcpy(pString, "Gamma");            
         if (inputPDFs[ii] == PSUADE_PDF_EXPONENTIAL) 
           strcpy(pString, "Exponential");            
         if (inputPDFs[ii] == PSUADE_PDF_SAMPLE) 
           strcpy(pString, "User");            
         if (inputPDFs[ii] == PSUADE_PDF_F) 
           strcpy(pString, "F");            
         printf("PDF type for Input %5d = %s\n",ii+1,pString);
       }
     }
     delete [] inputPDFs;
   }
   if (isScreenDumpModeOn() && repeatFlag_ == 0) printEquals(PL_INFO,0);

   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
   odata->optimalY_ = 1.0e50;
   XValues = new double[M1+1];
   index = 0;
   for (ii = 0; ii < nInputs; ii++)
   {
      if (psOUUInputTypes_[ii] == psOUUType1)
      {
         XValues[index] = odata->initialX_[ii];
         index++;
      }
   }
   if (M2 > 0)
   {
      psOUUM2Values_ = new double[M2];
      index = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         if (psOUUInputTypes_[ii] == psOUUType2)
         {
            psOUUM2Values_[index] = odata->initialX_[ii];
            index++;
         }
      }
   }
   rhobeg = 1e35;
   index = 0;
   for (ii = 0; ii < nInputs; ii++) 
   {
      if (psOUUInputTypes_[ii] == psOUUType1)
      {
         ddata = odata->upperBounds_[ii] - odata->lowerBounds_[ii];
         if (ddata < rhobeg) rhobeg = ddata;
         index++;
      }
   }
   rhobeg *= 0.5;
   rhoend = rhobeg * odata->tolerance_;
   if (rhobeg < rhoend && repeatFlag_ == 0)
   {
      printf("OUUOptimizer WARNING: tolerance too large.\n");
      printf("                      tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }

   maxfun = odata->maxFEval_;
   if ((odata->setOptDriver_ & 1) && (odata->funcIO_ != NULL))
   {
      if (repeatFlag_ == 0)
      {
         printf("OUUOptimizer: setting optimization simulation driver.\n");
         printf("              driver = %s\n",
                odata->funcIO_->getOptimizationDriver());
      }
      currDriver = odata->funcIO_->getDriver();
      odata->funcIO_->setDriver(1);
   }
   psOUUObj_= (void *) odata;

   psOUUMode_ = 1;
   psOUUCMode_ = 1;
   if ((repeatFlag_ == 0) && isScreenDumpModeOn())
   {
     psOUUUserOpt_ = 0;
     psOUUUseRS_ = 0;
     if ((psOptExpertMode_ == 1) && (nInputs != (M1+M2)) && optCode_ != 3)
     {
       printEquals(PL_INFO, 0);
       printf("Select which functional Phi_{Z3,Z4} to use: \n");
       printf("  1. mean of G(Z1,...) with respect to Z3,Z4 (default)\n"); 
       printf("  2. mean of G(Z1,...) + beta * std dev of G(Z1,...)\n"); 
       printf("  3. G(Z1,Z2,Z3*,Z4*) such that \n");
       printf("         Prob(G(Z1,...)>G(Z1,Z2,Z3*,Z4*)) = 1 - alpha\n");
       printf("     This is also called value-at-risk with confidence\n");
       printf("     level alpha.\n");
       sprintf(pString,"Enter your choice of functional (1, 2 or 3) : ");
       psOUUMode_ = getInt(1, 3, pString);
       if (psOUUMode_ == 2)
       {
         sprintf(pString,"Enter beta (>= 0) : ");
         psOUUStdevMultiplier_ = -1;
         while (psOUUStdevMultiplier_ < 0)
           psOUUStdevMultiplier_ = getDouble(pString);
       } 
       if (psOUUMode_ == 3)
       {
         psOUUPercentile_ = 0.0;
         sprintf(pString,"Enter the confidence interval : [0.5 - 1.0] ");
         while ((psOUUPercentile_ < 0.5) || (psOUUPercentile_ > 1.0)) 
         { 
           psOUUPercentile_ = getDouble(pString);
         } 
       }
       if (nOutputs > 1)
       {
         printf("Please select below the scheme to compute the returned\n");
         printf("constraints from the ensemble constraints.\n");
         printf("The returned constraints are computed from:\n");
         printf(" 1. the means of constraints from all sample runs\n");
         printf(" 2. the simulation with largest 2-norm violations\n");
         printf(" 3. the simulation with largest 1-norm violations\n");
         printf(" 4. the simulation with largest infinity-norm violation\n");
         printf("NOTE: constraints are such that they should be >= 0.\n");
         sprintf(pString,"Enter your choice (1 - 4) : ");
         psOUUCMode_ = getInt(1, 4, pString);
       }
     }
     else if ((psOptExpertMode_==1) && (nInputs != (M1+M2)) && optCode_==3)
     {
       printEquals(PL_INFO, 0);
       printf("Select which functional Phi_{Z3,Z4} to use: \n");
       printf("  1. mean of G(Z1,...) with respect to Z3,Z4 (default)\n"); 
       printf("  2. mean of G(Z1,...) + beta * std dev of G(Z1,...)\n"); 
       sprintf(pString,"Enter your choice of functional (1 or 2) : ");
       psOUUMode_ = getInt(1, 2, pString);
       if (psOUUMode_ == 2)
       {
         sprintf(pString,"Enter beta (>= 0) : ");
         psOUUStdevMultiplier_ = -1;
         while (psOUUStdevMultiplier_ < 0)
           psOUUStdevMultiplier_ = getDouble(pString);
       } 
     }
   }

   if (repeatFlag_ == 0) psOUUZ3nSamples_ = 1;
   if (M3 > 0 && repeatFlag_ == 0)
   {
     printEquals(PL_INFO, 0);
     printf("A sample for Z3 is needed from you. Data format is:\n");
     printf("line 1: <nSamples> <nInputs> \n");
     printf("line 2: <sample 1 input 1> <input 2> ... <probability>\n");
     printf("line 3: <sample 2 input 1> <input 2> ... <probability>\n");
     printf("...\n");
     printf("Enter user sample file name : ");
     scanf("%s", filename);
     fgets(lineIn, 5000, stdin);
     fp = fopen(filename, "r");
     if (fp == NULL)
     {
       printf("OUUOptimizer ERROR: user sample file %s not found.\n",
              filename);
       exit(1);
     }
     fgets(lineIn, 10000, fp);
     sscanf(lineIn, "%d %d", &psOUUZ3nSamples_, &ii);
     if (ii != M3)
     {
       printf("OUUOptimizer ERROR: user sample nInputs %d != %d\n",
              ii, M3);
       fclose(fp);
       exit(1);
     } 
     if (psOUUZ3nSamples_ < M3+1)
     {
       printf("OUUOptimizer ERROR: user sample size must be >= %d\n",
              M3+1);
       fclose(fp);
       exit(1);
     } 
     psOUUZ3SamInputs_ = new double[psOUUZ3nSamples_ * M3];
     psOUUSamProbs_    = new double[psOUUZ3nSamples_];
     tdata             = new double[psOUUZ3nSamples_];
     ddata = 0.0;
     if (odata->outputLevel_ > 4)
       printf("Reading sample for X3.\n");
     for (ii = 0; ii < psOUUZ3nSamples_; ii++)
     {
       fgets(lineIn, 20000, fp);
       index = kk = 0;
       if (hasColNums == 1)
       {
         while (lineIn[index] == ' ') index++;
         sscanf(&lineIn[index],"%lg",&ddata);
         itmp = (int) (ddata + 1.0e-15);
         if (itmp != (ii+1))
         {
           printf("OUU ERROR: X3 sample file should have sample number ");
           printf("in the first column, but it is not valid (%d != %d)\n",
                  itmp, ii);
           exit(1);
         } 
         while (lineIn[index] != ' ') index++;
       }
       while (kk < M3)
       {
         while (lineIn[index] == ' ') index++;
         if (lineIn[index] == '\0' || lineIn[index] == '\n')
         {
           printf("ERROR: reading sample file %s line %d.\n",
                  filename,ii+2);
           fclose(fp);
           exit(1);
         }
         sscanf(&lineIn[index],"%lg",&psOUUZ3SamInputs_[ii*M3+kk]);
         while (lineIn[index] != ' ') index++;
         if (lineIn[index] == '\0' || lineIn[index] == '\n')
         {
           printf("ERROR: reading sample file %s line %d.\n",
                  filename,ii+2);
           fclose(fp);
           exit(1);
         }
         if (odata->outputLevel_ > 4)
           printf("%12.4e ",psOUUZ3SamInputs_[ii*M3+kk]);
         kk++;
       }
       while (lineIn[index] == ' ') index++;
       if (lineIn[index] == '\0' || lineIn[index] == '\n')
       {
         printf("WARNING: reading probability at line %d of %s.\n",
                ii+2, filename);
         printf("    Will set probabilities to be equiprobable.\n");
         psOUUSamProbs_[ii] = PSUADE_UNDEFINED;
         ddata += 1.0;
       }
       else sscanf(&lineIn[index],"%lg",&psOUUSamProbs_[ii]);
       if (odata->outputLevel_ > 4)
         printf("%12.4e\n",psOUUSamProbs_[ii]);
     }
     fclose(fp);
     if (ddata != 0)
     {
       for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
         psOUUSamProbs_[ii] = 1.0 / (double) psOUUZ3nSamples_;
     }
     else
     {
       ddata = 0.0;
       for (ii = 0; ii < psOUUZ3nSamples_; ii++) 
         ddata += psOUUSamProbs_[ii];
     }
     printf("User sample for Z3 has %d points\n", psOUUZ3nSamples_);
     printf("User sample for Z3 CDF = %e (should be ~1)\n", ddata);
   }

   int    methodZ4=1,*samStates,samplingOption,*iPdfs,index2;
   double *samOuts,*lowers,*uppers,*inputMeans=NULL,*inputStdevs=NULL;
   double *iMeans, *iStdvs;
   char   **targv;
   Sampling   *sampler = NULL;
   PDFManager *pdfman=NULL;
   psVector   vecLB, vecUB, vecOut;
   psMatrix   *corMat1, corMat2;
   if (repeatFlag_ == 0) psOUUZ4nSamples_ = 1;
   if (M4 > 0 && repeatFlag_ == 0)
   {
     printf("OUUOptimizer: generating a sample for Z4. Two options:\n");
     printf("(1) Users uploads a sample to PSUADE\n");
     printf("(2) PSUADE can internally create a sample\n");
     sprintf(pString, "Select option 1 or 2 : ");
     kk = getInt(1,2,pString);
     if (kk == 1)
     {
       printEquals(PL_INFO, 0);
       printf("A Z4 sample is needed from you. The file format should be:\n");
       printf("line 1: <nSamples> <nInputs> \n");
       printf("line 2: <sample 1 input 1> <input 2> \n");
       printf("line 3: <sample 2 input 1> <input 2> \n");
       printf("...\n");
       printf("Enter user sample file name : ");
       scanf("%s", filename);
       fgets(lineIn, 5000, stdin);
       fp = fopen(filename, "r");
       if (fp == NULL)
       {
         printf("OUUOptimizer ERROR: user sample file %s not found.\n",
                filename);
         exit(1);
       }
       fgets(lineIn, 5000, fp);
       sscanf(lineIn, "%d %d", &psOUUZ4nSamples_, &ii);
       if (ii != M4)
       {
         printf("OUUOptimizer ERROR: user sample nInputs %d != %d\n",
                ii, M4);
         fclose(fp);
         exit(1);
       }
       if (psOUUZ4nSamples_ < M4+1)
       {
         printf("OUUOptimizer ERROR: user sample size should be >= %d\n",
                M4+1);
         fclose(fp);
         exit(1);
       }
       psOUUZ4SamInputs_ = new double[psOUUZ4nSamples_ * M4];
       for (ii = 0; ii < psOUUZ4nSamples_; ii++)
       {
         if (hasColNums == 1)
         {
           fscanf(fp, "%lg", &ddata);
           itmp = (int) (ddata + 1.0e-15);
           if (itmp != (ii+1))
           {
             printf("OUU ERROR: X4 sample file should have sample number");
             printf(" in the first column, but it is not valid.\n");
           }
         }
         for (kk = 0; kk < M4; kk++)
           fscanf(fp,"%lg",&psOUUZ4SamInputs_[ii*M4+kk]);
       }
       fclose(fp);
       count = 0;
       for (ii = 0; ii < psOUUZ4nSamples_; ii++)
       {
         index = 0;
         for (kk = 0; kk < nInputs; kk++)
         {
           if (psOUUInputTypes_[kk] == psOUUType4)
           {
             ddata = psOUUZ4SamInputs_[ii*M4+index];
             if ((ddata-1e-10) > odata->upperBounds_[kk] ||
                 (ddata+1e-10) < odata->lowerBounds_[kk])
             {
               printf("WARNING: X4 sample not within input bounds.\n");
               printf("    Input %d Bounds are: [%e, %e]\n",kk+1,
                      odata->lowerBounds_[kk],odata->upperBounds_[kk]);
               count++;
             }
             index++;
           }
         }
       }
       if (count > 0)
       {
         count = 10;
         if (count > psOUUZ4nSamples_) count = psOUUZ4nSamples_;
         for (ii = 0; ii < count; ii++)
         {
           printf("Z4 %2d : ", ii+1);
           for (kk = 0; kk < M4; kk++)
             printf("%12.4e ", psOUUZ4SamInputs_[ii*M4+kk]);
           printf("\n");
         }
       }
       if (nOutputs == 1)
       {
         printf("The user sample for Z4 has %d points\n",psOUUZ4nSamples_);
         printf("You have the option to select a subset of Z4 to build a\n");
         printf("response surface and use the original larger Z4 sample\n");
         printf("to estimate the statistics from the response surface.\n");
         printf("Use response surface for Z4? (y or n) ");
         scanf("%s", lineIn);
         if (lineIn[0] == 'y') psOUUUseRS_ = 1;
       }
       else
       {
         psOUUUseRS_ = 0;
         printf("OUUOptimizer INFO: response surface not available \n");
         printf("             for problems with constraints.\n");
       }
       if (psOUUUseRS_ == 1)
       {
         fgets(lineIn, 5000, stdin);
         psOUULargeSampleSize_ = psOUUZ4nSamples_;
         psOUULargeSamInputs_= new double[psOUULargeSampleSize_*M4];
         psOUULargeSamOutputs_ = new double[psOUULargeSampleSize_*nOutputs];
         for (ii = 0; ii < psOUULargeSampleSize_*M4; ii++)
            psOUULargeSamInputs_[ii] = psOUUZ4SamInputs_[ii]; 
         for (ii = 0; ii < psOUULargeSampleSize_*nOutputs; ii++)
            psOUULargeSamOutputs_[ii] = 0.0;
         printf("Your Z4 sample size is %d.\n", psOUUZ4nSamples_);
         printf("This sample size may be too large for building a RS.\n");
         sprintf(pString,
                 "Number of points to use for building RS? (%d - %d) ",
                 M4+1, psOUUZ4nSamples_);
         psOUUZ4nSamples_ = getInt(M4+1, psOUUZ4nSamples_, pString);
         printf("You have 2 options on how to generate this RS set:\n");
         printf("(1) You upload another Z4 sample of size %d\n",
                psOUUZ4nSamples_);
         printf("(2) PSUADE randomly draws %d points from your sample\n",
                psOUUZ4nSamples_);
         sprintf(pString,"Select option 1 or 2 : ");
         kk = getInt(1, 2, pString);
         if (kk == 2)
         {
           printf("OUU will randomly select a subset of points from\n");
           printf("the original Z4 sample to build response surface.\n");
           count = 0;
           for (ii = 0; ii < psOUULargeSampleSize_; ii++)
           {
             if (psOUUZ4nSamples_ == psOUULargeSampleSize_) index = ii;
             else index = PSUADE_rand() % (psOUULargeSampleSize_-ii) + ii; 
             for (kk = 0; kk < M4; kk++)
               psOUUZ4SamInputs_[count*M4+kk] = 
                            psOUULargeSamInputs_[index*M4+kk];
             count++;
             if (count == psOUUZ4nSamples_) break;
           }
         }
         else if (kk == 1)
         {
           printEquals(PL_INFO, 0);
           printf("A Z4 RS sample is needed. The file format is:\n");
           printf("line 1: <nSamples> <nInputs> \n");
           printf("line 2: <sample 1 input 1> <input 2> \n");
           printf("line 3: <sample 2 input 1> <input 2> \n");
           printf("...\n");
           printf("Enter user sample file name : ");
           scanf("%s", filename);
           fgets(lineIn, 5000, stdin);
           fp = fopen(filename, "r");
           if (fp == NULL)
           {
             printf("OUUOptimizer ERROR: sample file %s not found.\n",
                    filename);
             exit(1);
           }
           fgets(lineIn, 5000, fp);
           sscanf(lineIn, "%d %d", &kk, &ii);
           if (ii != M4 || kk != psOUUZ4nSamples_)
           {
             printf("OUUOptimizer ERROR: parameter mismatch.\n");
             fclose(fp);
             exit(1);
           }
           for (ii = 0; ii < psOUUZ4nSamples_; ii++)
           {
             for (kk = 0; kk < M4; kk++)
               fscanf(fp,"%lg",&psOUUZ4SamInputs_[ii*M4+kk]);
           }
           fclose(fp);
         }
       }
     }
     else
     {
       psOUUZ4nSamples_ = 100;
       printEquals(PL_INFO, 0);
       printf("PSUADE will create a sample (size will be selected later)\n");
       printf("for Z4. Since the sample should be small for computational\n");
       printf("efficiency, and accuracy of statistics depends on sample\n");
       printf("size, you have the option to add one more step in OUU by\n");
       printf("estimating the statistics with a large sample evaluated on\n");
       printf("the response surfaces built from the small sample.\n");
       printf("Use response surface for Z4 to compute statistics? (y or n) ");
       scanf("%s", lineIn);
       if (lineIn[0] == 'y') psOUUUseRS_ = 1;
       if (psOUUMasterMode_ == 1 && psOUUUseRS_ == 1)
       {
         printf("For diagnostics purposes, you may choose to validate\n");
         printf("the response surface constructed at each OUU iteration.\n");
         printf("The more accurate the cross validation result is, the\n");
         printf("more accurate is the statistics.\n"); 
         printf("Validate the response surfaces for Z4? (y or n) ");
         scanf("%s", lineIn);
         if (lineIn[0] == 'y') psOUUValidateRS_ = 1;
       }
       fgets(lineIn, 5000, stdin);
       samplingOption = 0;
       if (psOUUUseRS_ == 0)
       {
         odata->psIO_->getParameter("input_pdfs", pdata);
         inputPDFs = pdata.intArray_;
         pdata.intArray_ = NULL;
         kk = index = 0;
         if (inputPDFs != NULL)
         {
           for (ii = 0; ii < nInputs; ii++) 
           {
             if (psOUUInputTypes_[ii] == psOUUType4)
             {
               kk += inputPDFs[ii];
               index++;
             }
           }
         }
         if (kk > 0) samplingOption = 1;
       }
       if (samplingOption == 0)
       {
         printEquals(PL_INFO, 0);
         printf("OUUOptimizer uses a Z4 sample to estimate the objective\n");
         printf("Default sampling method = Latin hypercube\n");
         printf("Default sample size     = %d\n",psOUUZ4nSamples_);
         printf("Available sampling method: \n");
         printf("   (1) LHS, \n");
         printf("   (2) factorial, or\n");
         printf("   (3) quasi-MC.\n");
         sprintf(pString,"Select sampling method (1, 2 or 3) : ");
         methodZ4 = getInt(1, 3, pString);
         if (methodZ4 == 1 || methodZ4 == 3)
         {
            kk = 2;
            if (psOUUUseRS_ == 1) kk = M4 + 1;
            sprintf(pString, "Enter sample size (>= %d, <= 1000) : ", kk);
            psOUUZ4nSamples_ = getInt(2, 10000, pString);
            printf("Latin hypercube/QMC has sample size = %d\n",
                   psOUUZ4nSamples_);  
         }
         else if (methodZ4 == 2)
         {
            sprintf(pString,
                    "Enter number of levels per variable (>=2, <=100) : ");
            psOUUZ4nSamples_ = getInt(2, 100, pString);
            kk = psOUUZ4nSamples_;
            for (ii = 1; ii < M4; ii++) psOUUZ4nSamples_ *= kk;
            printf("Factorial design has sample size = %d\n", 
                   psOUUZ4nSamples_);  
         }
         if (methodZ4 == 1)
              sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
         else if (methodZ4 == 2)
              sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FACT);
         else sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         sampler->setPrintLevel(0);
         lowers = new double[M4];
         uppers = new double[M4];
         index = 0;
         for (ii = 0; ii < nInputs; ii++)
         {
            if (psOUUInputTypes_[ii] == psOUUType4)
            {
               lowers[index] = odata->lowerBounds_[ii];
               uppers[index] = odata->upperBounds_[ii];
               index++;
            }
         }
         sampler->setInputBounds(M4, lowers, uppers);
         sampler->setOutputParams(iOne);
         sampler->setSamplingParams(psOUUZ4nSamples_, iOne, iZero);
         sampler->initialize(0);
         psOUUZ4nSamples_ = sampler->getNumSamples();
         psOUUZ4SamInputs_ = new double[psOUUZ4nSamples_ * M4];
         samStates = new int[psOUUZ4nSamples_];
         samOuts = new double[psOUUZ4nSamples_];
         sampler->getSamples(psOUUZ4nSamples_, M4, iOne, psOUUZ4SamInputs_, 
                             samOuts, samStates);
         delete sampler;
         delete [] lowers;
         delete [] uppers;
         delete [] samStates;
         delete [] samOuts;
         printEquals(PL_INFO, 0);
       }
       else
       {
         printEquals(PL_INFO, 0);
         printf("OUUOptimizer uses a Z4 sample to estimate the objective\n");
         printf("Default sample size     = %d\n",psOUUZ4nSamples_);
         sprintf(pString,"Enter your desired sample size (>=10, <=1000) : ");
         psOUUZ4nSamples_ = getInt(10, 10000, pString);
         odata->psIO_->getParameter("input_means", pdata);
         inputMeans = pdata.dbleArray_;
         if (inputMeans == NULL)
         {
           inputMeans = new double[nInputs];
           for (ii = 0; ii < nInputs; ii++) inputMeans[ii] = 0;
         }
         pdata.dbleArray_ = NULL;
         odata->psIO_->getParameter("input_stdevs", pdata);
         inputStdevs = pdata.dbleArray_;
         if (inputStdevs == NULL)
         {
           inputStdevs = new double[nInputs];
           for (ii = 0; ii < nInputs; ii++) inputStdevs[ii] = 1;
         }
         pdata.dbleArray_ = NULL;
         odata->psIO_->getParameter("input_cor_matrix", pdata);
         corMat1 = (psMatrix *) pdata.psObject_;
         pdata.psObject_ = NULL;

         corMat2.setDim(M4,M4);
         iPdfs  = new int[nInputs];
         iMeans = new double[nInputs];
         iStdvs = new double[nInputs];
         lowers = new double[nInputs];
         uppers = new double[nInputs];
         index = 0;
         for (ii = 0; ii < nInputs; ii++)
         {
           if (psOUUInputTypes_[ii] == psOUUType4)
           { 
             iMeans[index] = inputMeans[ii];
             iStdvs[index] = inputStdevs[ii];
             iPdfs[index]  = inputPDFs[ii];
             lowers[index] = odata->lowerBounds_[ii];
             uppers[index] = odata->upperBounds_[ii];
             index2 = 0;
             for (kk = 0; kk < nInputs; kk++)
             {
               if (psOUUInputTypes_[kk] == psOUUType4)
               {
                 ddata = corMat1->getEntry(ii,kk);
                 corMat2.setEntry(index,index2,ddata);
                 index2++;
               }
             }
             index++;
           }
         }
         pdfman = new PDFManager();
         pdfman->initialize(M4,iPdfs,iMeans,iStdvs,corMat2,NULL,NULL);
         vecLB.load(M4, lowers);
         vecUB.load(M4, uppers);
         vecOut.setLength(psOUUZ4nSamples_*M4);
         pdfman->genSample(psOUUZ4nSamples_, vecOut, vecLB, vecUB);
         psOUUZ4SamInputs_= new double[psOUUZ4nSamples_*M4];
         for (ii = 0; ii < psOUUZ4nSamples_*M4; ii++)
            psOUUZ4SamInputs_[ii] = vecOut[ii];
         //psOUUSamOutputs_ = new double[psOUUZ4nSamples_*nOutputs];
         delete [] inputPDFs;
         delete [] inputMeans;
         delete [] inputStdevs;
         delete [] iPdfs;
         delete [] iMeans;
         delete [] iStdvs;
         delete [] lowers;
         delete [] uppers;
         delete pdfman;
       }
       if (psOUUUseRS_ == 1)
       {
         odata->psIO_->getParameter("input_pdfs", pdata);
         inputPDFs = pdata.intArray_;
         pdata.intArray_ = NULL;
         if (inputPDFs == NULL)
         {
           inputPDFs = new int[nInputs];
           for (ii = 0; ii < nInputs; ii++) inputPDFs[ii] = 0;
         }
         kk = 0; index = 0;
         for (ii = 0; ii < nInputs; ii++)
         {
           if (psOUUInputTypes_[ii] == psOUUType4)
           {
             kk += inputPDFs[ii];
             index++;
           }
         }
         lowers = new double[M4];
         uppers = new double[M4];
         index = 0;
         for (ii = 0; ii < nInputs; ii++)
         {
           if (psOUUInputTypes_[ii] == psOUUType4)
           {
             lowers[index] = odata->lowerBounds_[ii];
             uppers[index] = odata->upperBounds_[ii];
             index++;
           }
         }
         if (kk == 0)
         {
           sampler = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
           sampler->setPrintLevel(0);
           sampler->setInputBounds(M4, lowers, uppers);
           sampler->setOutputParams(iOne);
           psOUULargeSampleSize_ = 100000;
           sampler->setSamplingParams(psOUULargeSampleSize_, iOne, iZero);
           sampler->initialize(0);
           psOUULargeSampleSize_ = sampler->getNumSamples();
           psOUULargeSamInputs_= new double[psOUULargeSampleSize_*M4];
           kk = psOUULargeSampleSize_*nOutputs;
           psOUULargeSamOutputs_ = new double[kk];
           samStates = new int[psOUULargeSampleSize_];
           sampler->getSamples(psOUULargeSampleSize_, M4, iOne,
                        psOUULargeSamInputs_,psOUULargeSamOutputs_,
                        samStates);
           delete [] samStates;
           delete sampler;
         }
         else
         {
           odata->psIO_->getParameter("input_means", pdata);
           inputMeans = pdata.dbleArray_;
           if (inputMeans == NULL)
           {
             inputMeans = new double[nInputs];
             for (ii = 0; ii < nInputs; ii++) inputMeans[ii] = 0;
           }
           pdata.dbleArray_ = NULL;
           odata->psIO_->getParameter("input_stdevs", pdata);
           inputStdevs = pdata.dbleArray_;
           if (inputStdevs == NULL)
           {
             inputStdevs = new double[nInputs];
             for (ii = 0; ii < nInputs; ii++) inputStdevs[ii] = 1;
           }
           pdata.dbleArray_ = NULL;
           odata->psIO_->getParameter("input_cor_matrix", pdata);
           corMat1 = (psMatrix *) pdata.psObject_;
           pdata.psObject_ = NULL;

           corMat2.setDim(M4,M4);
           iPdfs  = new int[nInputs];
           iMeans = new double[nInputs];
           iStdvs = new double[nInputs];
           index = 0;
           for (ii = 0; ii < nInputs; ii++)
           {
             if (psOUUInputTypes_[ii] == psOUUType4)
             {
               iMeans[index] = inputMeans[ii];
               iStdvs[index] = inputStdevs[ii];
               iPdfs[index]  = inputPDFs[ii];
               index2 = 0;
               for (kk = 0; kk < nInputs; kk++)
               {
                 if (psOUUInputTypes_[kk] == psOUUType4)
                 {
                   ddata = corMat1->getEntry(ii,kk);
                   corMat2.setEntry(index,index2,ddata);
                   index2++;
                 }
               }
               index++;
             }
           }
           pdfman = new PDFManager();
           pdfman->initialize(M4,iPdfs,iMeans,iStdvs,corMat2,NULL,NULL);
           vecLB.load(M4, lowers);
           vecUB.load(M4, uppers);
           psOUULargeSampleSize_ = 100000;
           vecOut.setLength(psOUULargeSampleSize_*M4);
           pdfman->genSample(psOUULargeSampleSize_, vecOut, vecLB, vecUB);
           psOUULargeSamInputs_= new double[psOUULargeSampleSize_*M4];
           kk = psOUULargeSampleSize_*nOutputs;
           psOUULargeSamOutputs_ = new double[kk];
           for (ii = 0; ii < psOUULargeSampleSize_*M4; ii++)
             psOUULargeSamInputs_[ii] = vecOut[ii];
           delete [] inputMeans;
           delete [] inputStdevs;
           delete [] iPdfs;
           delete [] iMeans;
           delete [] iStdvs;
           delete pdfman;
         }
         delete [] inputPDFs;
         delete [] lowers;
         delete [] uppers;
       }
     }
   }

   int nSamp=psOUUZ3nSamples_*psOUUZ4nSamples_;
   if (repeatFlag_ == 0)
   {
     psOUUSamOutputs_ = new double[nSamp*nOutputs];
     psOUUSamConstrNorms_ = new double[nSamp];
   }
   if (M2 > 0 && repeatFlag_ == 0)
   {
     printf("For 2-stage OUU, you have 2 options for inner optimization:\n");
     printf("(1) You can use BOBYQA available in PSUADE (in this case\n");
     printf("    opt_driver should point to your original function), or\n");
     printf("(2) You can provide your own optimizer (in opt_driver), or\n");
     printf("    it is one-level optimization.\n");
     printf("NOTE: Most likely you will choose (2) or answer y below.\n");
     printf("Use your own optimizer instead of BOBYQA? (y or n) ");
     scanf("%s", lineIn);
     if (lineIn[0] == 'y') 
     {
       psOUUUserOpt_ = 1;
       printf("NOTE: Make sure your optimizer executable has been\n");
       printf("      assigned to 'opt_driver' and it optimizes with\n");
       printf("      respect to the %d-th to %d-th parameters.\n", 
              M1+1, M1+M2);
     }
   }
   else if (repeatFlag_ == 0 && isScreenDumpModeOn())
   {
     printf("Since no recourse variable (for level 2 optimization) has\n");
     printf("been specified, PSUADE presumes that you are either doing\n");
     printf("1-level OUU, or you are providing the inner optimization\n");
     printf("solver in opt_driver.\n");
     psOUUUserOpt_ = 1;
   }
   if (repeatFlag_ == 0 && isScreenDumpModeOn()) printEquals(PL_INFO, 0);
 
   if (psOUUUserOpt_ == 1 && (M3+M4) > 0 && repeatFlag_ == 0)
   {
     printf("Each simulation calls the opt_driver with 1 sample point.\n");
     printf("For higher efficiency (less I/O), you have the option to\n");
     printf("provide in 'ensemble_opt_driver' an executable that can\n");
     printf("run multiple sample points. In this case, OUU calls the\n");
     printf("ensemble executable with the following sequence: \n");
     printf("    <ensemble_opt_driver> <sampleFile> <outputFile>\n\n");
     printf("where <sampleFile> is in the following format:\n");
     printf("   line 1: <nSamples>\n");
     printf("   line 2: Sample point 1 input values\n");
     printf("   line 3: Sample point 2 input values\n");
     printf("   line n: ...\n\n");
     printf("and <outputFile> should have all sample output values.\n\n");
     printf("Use ensemble opt driver for ensemble runs ? (y or n) ");
     lineIn[0] = '1';
     while (lineIn[0] != 'n' && lineIn[0] != 'y')
     {
       scanf("%s", lineIn);
       if (lineIn[0] == 'y') psOUUEnsembleEval_ = 1;
       printEquals(PL_INFO, 0);
     }
   }

   if ((psOUUEnsembleEval_ == 0) && ((M3+M4)>0) && repeatFlag_ == 0)
   {
     printf("You can configure OUU to run the ensemble simulations in\n");
     printf("parallel/asynchronous using the Linux fork/join. If 'n'\n");
     printf("is entered below, the opt_driver simulator will be\n");
     printf("evaluated sequentially (one sample point at a time).\n");
     printf("If 'y' is selected instead, be careful, because PSUADE\n");
     printf("PSUADE will launch %d jobs simultaneously, which may\n",
            psOUUZ3nSamples_*psOUUZ4nSamples_);
     printf("jam up the job queuing system.\n");
     printf("Turn on asynchronous mode ? (y or n) ");
     lineIn[0] = '1';
     while (lineIn[0] != 'n' && lineIn[0] != 'y')
     {
       scanf("%s", lineIn);
       if (lineIn[0] == 'y') psOUUParallel_ = 1;
       printEquals(PL_INFO, 0);
     }
   }
   if (((M3+M4) > 0 && (psOUUUserOpt_ == 1 || psOUUEnsembleEval_ == 0)) &&
       repeatFlag_ == 0)
     fgets(lineIn, 500, stdin);

   if (repeatFlag_ == 0)
   {
     psOUUWValues_ = new double[nInputs];
     psOUUXValues_ = new double[nInputs*nSamp];
     psOUUOptimalX_ = new double[nInputs];
   }

   int rstype=0, rsaux=0;
   if (repeatFlag_ == 0 && psOUUUseRS_ == 1 && M4 > 0)
   {
     psOUUfaPtr_ = NULL;
     if (psOUUMasterMode_ == 0)
     {
       if (printLevel > 2) 
          printf("OUUOptimizer: setting up response surface\n");
       if (psOUUZ4nSamples_ > 600)
       {
         rstype = PSUADE_RS_MARS;
         printf("OUUOptimizer: use MARS since nSamples > 400\n");
       }
       else if (psOUUZ4nSamples_ > 300)
       {
         rstype = PSUADE_RS_RBF;
         printf("OUUOptimizer: use RBF response surface\n");
       }
       else if (psOUUZ4nSamples_ > 200)
       {
         rstype = PSUADE_RS_KR;
         printf("OUUOptimizer: use Kriging (fast) response surface\n");
       }
       else
       {
         rstype = PSUADE_RS_KR;
         rsaux = 1;
         printf("OUUOptimizer: use Kriging (slow) response surface\n");
       }
     }
     else
     {
       printf("Response surface available: \n");
       printf("1. MARS\n");
       printf("2. Kriging (fast)\n");
       printf("3. Kriging (slow)\n");
       printf("4. Radial basis function\n");
       sprintf(pString, "Which response surface? ");
       rstype = getInt(1,4,pString); 
       if (rstype == 1) rstype = PSUADE_RS_MARS;
       if (rstype == 2)
       {
         rstype = PSUADE_RS_KR;
         rsaux = 0;
       }
       if (rstype == 3)
       {
          rstype = PSUADE_RS_KR;
          rsaux = 1;
       }
       if (rstype == 4) rstype = PSUADE_RS_RBF;
     }
     kk = psInteractive_;
     psInteractive_ = 0;
     psOUUfaPtr_ = genFA(rstype, M4, -1, psOUUZ4nSamples_);
     lowers = new double[M4];
     uppers = new double[M4];
     psOUUZ4RSType_ = rstype;
     psOUUZ4RSAux_ = rsaux;
     psOUUZ4LBounds_ = new double[M4];
     psOUUZ4UBounds_ = new double[M4];
     index = 0;
     for (ii = 0; ii < nInputs; ii++)
     {
       if (psOUUInputTypes_[ii] == psOUUType4)
       {
         lowers[index] = odata->lowerBounds_[ii];
         uppers[index] = odata->upperBounds_[ii];
         psOUUZ4LBounds_[index] = lowers[index];
         psOUUZ4UBounds_[index] = uppers[index];
         index++;
       }
     }
     psOUUfaPtr_->setBounds(lowers,uppers);
     psOUUfaPtr_->setOutputLevel(0);
     psInteractive_ = kk;
     if (rstype == PSUADE_RS_KR) 
     {
       targv = new char*[1];
       targv[0] = new char[100];
       if (rsaux == 0) strcpy(targv[0], "setMode2");
       else            strcpy(targv[0], "setMode3");
       psOUUfaPtr_->setParams(1, targv);
       delete [] targv[0];
       delete [] targv;
     }
     delete [] lowers;
     delete [] uppers;
   }

   nPts = (M1 + 1) * (M1 + 2) / 2;
   psOUUCounter_ = 0;
   if (psOUUParallel_ == 1 & odata->funcIO_ != NULL) 
     odata->funcIO_->setAsynchronousMode();

   if (psConfig_ != NULL)
   {
     cString = psConfig_->getParameter("opt_save_history");
     if (cString != NULL) psOUUSaveHistory_ = 1;
     cString = psConfig_->getParameter("opt_use_history");
     if (cString != NULL) 
     {
       printf("OUU: use history has been turned on.\n");
       fp = fopen("psuade_ouu_history", "r");
       if (fp != NULL)
       {
         printf("OUUOptimizer history file found.\n");
         psOUUNSaved_ = 0;
         while (feof(fp) == 0)
         {
           fscanf(fp, "%d %d", &ii, &kk);
           if (ii != 999 || kk != nInputs)
           {
             printf("OUU history file not used - nInputs mismatch.\n");
             break;
           }
           else
           {
             for (ii = 0; ii < nInputs; ii++)
               fscanf(fp, "%lg",&psOUUSaveX_[psOUUNSaved_*nInputs+ii]);
             for (ii = 0; ii < nOutputs; ii++)
               fscanf(fp, "%lg",&psOUUSaveY_[psOUUNSaved_*nOutputs+ii]);
             psOUUNSaved_++;
           }
           if (((psOUUNSaved_+1)*nInputs > psOUUMaxSaved_*10) ||
             psOUUNSaved_ > psOUUMaxSaved_) break;
         }
         fclose(fp);
       }
     }
   }

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"*** Optimization Under Uncertainty begins\n");
   if (isScreenDumpModeOn()) 
   {
     if (printLevel >= 0 && repeatFlag_ == 0)
     {
       printDashes(PL_INFO, 0);
       printf(" ** OUUOptimizer summary of optimization parameters: \n");
       printf("  * Number of first  stage optimization parameters = %d\n",M1);
       printf("  * Number of second stage optimization parameters = %d\n",M2);
       printf("  * Number of discrete   uncertain parameters      = %d\n",M3);
       printf("  * Number of continuous uncertain parameters      = %d\n",M4);
       printDashes(PL_INFO, 0);
       printf(" ** max fevals = %d\n", odata->maxFEval_);
       printf(" ** tolerance  = %e\n", odata->tolerance_);
       printDashes(PL_INFO, 0);
     }
   }
   switch (optCode_)
   {
     case 0: printOutTS(PL_INFO,
             " ** Optimizer = bobyqa (bound-constrained)\n");
             break;
     case 1: printOutTS(PL_INFO,
             " ** Optimizer = newuoa (unconstrained)\n");
             break;
     case 2: printOutTS(PL_INFO,
             " ** Optimizer = cobyla (inequality constrained)\n");
             break;
     case 3: printOutTS(PL_INFO,
             " ** Optimizer = LBFGS (derivative-based)\n");
             break;
     case 4: printOutTS(PL_INFO,
             " ** Optimizer = lincoa (linear inequality constrained)\n");
             break;
     case 5: printOutTS(PL_INFO,
             " ** Optimizer = nomad (mixed integer)\n");
             break;
   }
   psOUUOptCode_ = optCode_;
   if (isScreenDumpModeOn())
   {
     printDashes(PL_INFO, 0);
     for (ii = 0; ii < M1; ii++) 
        printf("OUUOptimizer initial X %3d = %e\n", ii+1, XValues[ii]);
     printEquals(PL_INFO, 0);
   }
   if (optCode_ == 0)
   {
#ifdef HAVE_BOBYQA
     kk = (nPts + 5) * (nPts + M1) + 3 * M1 * (M1 + 5) / 2 + 1;
     workArray = new double[kk];
     lowers = new double[M1];
     uppers = new double[M1];
     index  = 0;
     for (ii = 0; ii < nInputs; ii++)
     {
        if (psOUUInputTypes_[ii] == psOUUType1)
        {
           lowers[index] = odata->lowerBounds_[ii];
           uppers[index] = odata->upperBounds_[ii];
           index++;
        }
     }
     if (nOutputs > 1) 
        printf("OUUOptimizer WARNING: multiple outputs, use Output 1 only\n");
     if (isScreenDumpModeOn()) printf("OUUOptimizer: calling bobyqa\n");
     bobyqa_(&M1, &nPts, XValues, lowers, uppers,
             &rhobeg, &rhoend, &bobyqaFlag, &maxfun, workArray);
     delete [] workArray;
     delete [] lowers;
     delete [] uppers;
#else
     printf("OUUOptimizer ERROR: bobyqa not installed.\n");
     exit(1);
#endif
   }
   else if (optCode_ == 1)
   {
#ifdef HAVE_NEWUOA
     int newuoaFlag = 8888;
     kk = (nPts + 13) * (nPts + M1) + 3 * M1 * (M1 + 3) / 2;
     workArray = new double[kk];
     if (nOutputs > 1) 
        printf("OUUOptimizer WARNING: multiple outputs, use Output 1 only\n");
     if (isScreenDumpModeOn()) printf(" ** OUUOptimizer: calling newuoa\n");
     newuoa_(&M1, &nPts, XValues, &rhobeg, &rhoend, &newuoaFlag,
             &maxfun, workArray);
     delete [] workArray;
#else
     printf("OUUOptimizer ERROR: newuoa not installed.\n");
     exit(1);
#endif
   }
   else if (optCode_ == 2)
   {
#ifdef HAVE_COBYLA
     if (nOutputs == 1) 
        printf("OUUOptimizer INFO: no constraint declared.\n");
     if (isScreenDumpModeOn()) printf(" ** OUUOptimizer: calling cobyla\n");
     int cobylaNConstr = nOutputs - 1 + 2 * M1;
     psOUULBounds_ = odata->lowerBounds_;
     psOUUUBounds_ = odata->upperBounds_;
     cobyla(M1, cobylaNConstr, XValues, rhobeg, rhoend, 0,
            &maxfun, ouuevalfunccobyla, (void *) odata);
     psOUULBounds_ = NULL;
     psOUUUBounds_ = NULL;
#else
     printf("OUUOptimizer ERROR: cobyla not installed.\n");
     exit(1);
#endif
   }
   else if (optCode_ == 3)
   {
#ifdef HAVE_LBFGS
     if (nOutputs != M1+1) 
     {
       printf("OUUOptimizer ERROR: BFGS must have nOutputs=%d\n",M1+1);
       exit(1);
     }
     if (psOUUSaveHistory_ == 1)
     {
       psOUUSaveHistory_ = 0;
       printf("OUUOptimizer INFO: no optimization history will be saved.\n");
     }
     if (isScreenDumpModeOn()) printf("OUUOptimizer: calling lbfgs\n");

     integer nInps, iprint=0, itask, *task=&itask, lsave[4], isave[44];
     integer csave[60], *iwork, nCorr=5, *nbds, its;
     int     *readys, M, ss, funcID;
     double  factr, pgtol, *work, dsave[29], FValue, *YValues, *GValues;
     double  *XLocal, *means, mean, stdev;

     nInps = (integer) M1;
     nbds  = new integer[nInps];
     for (ii = 0; ii < nInps; ii++) nbds[ii] = 2;
     GValues = new double[nSamp*nInps];
     if (GValues == NULL)
     {
       printf("Fatal memory allocation error (OUU:GValues)\n");
       exit(1);
     }
     for (ii = 0; ii < nInps*nSamp; ii++) GValues[ii] = 0.0;
     factr = 1e7;
     pgtol = 1e-6;
     kk = (2 * nCorr + 5) * nInps + 11 * nCorr * nCorr + 8 * nCorr;
     work  = new double[kk];
     for (ii = 0; ii < kk; ii++) work[ii] = 0.0;
     iwork = new integer[3*nInps];
     for (ii = 0; ii < 3*nInps; ii++) iwork[ii] = 0;
     *task = (integer) START;
     its   = 0;
     M = M1 + M2 + M3 + M4;
     readys = new int[nSamp];
     means  = new double[M];
     lowers = new double[M1];
     uppers = new double[M1];
     index  = 0;
     for (ii = 0; ii < nInputs; ii++)
     {
       if (psOUUInputTypes_[ii] == psOUUType1)
       {
         lowers[index] = odata->lowerBounds_[ii];
         uppers[index] = odata->upperBounds_[ii];
         index++;
       }
     }
     YValues = new double[nOutputs];

     while (1)
     {
       its++;
       setulb(&nInps,&nCorr,XValues,lowers,uppers,nbds,&FValue,
              GValues, &factr, &pgtol, work, iwork, task, 
              &iprint, csave, lsave, isave, dsave);

       if (IS_FG(*task))
       {
         index = 0;
         for (ii = 0; ii < M; ii++)
         {
           if (psOUUInputTypes_[ii] == psOUUType1)
           {
             psOUUXValues_[ii] = XValues[index];
             index++;
           }
         }

         psOUUOptimalY_ = PSUADE_UNDEFINED;
         psOUUCounter_ = 0;
         if (psOUUEnsembleEval_ == 0)
         {
           for (ss = 0; ss < nSamp; ss++)
           {
             readys[ss] = -1;

             index = 0;
             for (ii = 0; ii < M; ii++)
             {
               if (psOUUInputTypes_[ii] == psOUUType3)
               {
                 kk = ss / psOUUZ4nSamples_;
                 psOUUXValues_[ii] = psOUUZ3SamInputs_[kk*M3+index];
                 index++;
               }
             }
             index = 0;
             for (ii = 0; ii < M; ii++)
             {
               if (psOUUInputTypes_[ii] == psOUUType4)
               {
                 kk = ss % psOUUZ4nSamples_;
                 psOUUXValues_[ii] = psOUUZ4SamInputs_[kk*M4+index];
                 index++;
               }
             }

             funcID = psOUUCounter_ * nSamp + ss;
             if (odata->optFunction_ != NULL)
             {
               odata->optFunction_(M,psOUUXValues_,nOutputs,YValues);
               readys[ss] = 0;
             }
             else if (odata->funcIO_ != NULL)
             {
               readys[ss] = odata->funcIO_->evaluate(funcID,M,
                                 psOUUXValues_,nOutputs,YValues,0);
             }
             psOUUSamOutputs_[ss] = YValues[0];
             for (ii = 0; ii < M1; ii++) 
               GValues[ss*M1+ii] = YValues[ii+1];
             if (psOUUParallel_ == 0) odata->numFuncEvals_++;
           }
           if (psOUUParallel_ == 1)
           {
             for (ss = 0; ss < nSamp; ss++)
             {
               funcID = psOUUCounter_ * nSamp + ss;
               if (readys[ss] != 0)
               {
                 while (readys[ss] != 0)
                 {
#ifdef WINDOWS
                   Sleep(1000);
#else
                   sleep(1);
#endif
                   readys[ss] = odata->funcIO_->evaluate(funcID,M,
                                psOUUXValues_,nOutputs,YValues,2);
                 }
                 psOUUSamOutputs_[ss] = YValues[0];
                 for (ii = 0; ii < M1; ii++)
                   GValues[ss*M1+ii] = YValues[ii+1];
                 odata->numFuncEvals_++;
               }
             }
           }
         }
         else
         {
           for (ss = 0; ss < nSamp; ss++)
           {
             index = 0;
             for (ii = 0; ii < M; ii++)
             {
               if (psOUUInputTypes_[ii] == psOUUType1)
               {
                 psOUUXValues_[ss*M+ii] = XValues[index];
                 index++;
               }
             }
             index = 0;
             for (ii = 0; ii < M; ii++)
             {
               if (psOUUInputTypes_[ii] == psOUUType3)
               {
                 kk = ss / psOUUZ4nSamples_;
                 psOUUXValues_[ss*M+ii] = 
                         psOUUZ3SamInputs_[kk*M3+index];
                 index++;
               }
             }
             index = 0;
             for (ii = 0; ii < M; ii++)
             {
               if (psOUUInputTypes_[ii] == psOUUType4)
               {
                 kk = ss % psOUUZ4nSamples_;
                 psOUUXValues_[ss*M+ii] = 
                          psOUUZ4SamInputs_[kk*M4+index];
                 index++;
               }
             }
           }
           funcID = odata->numFuncEvals_;
           odata->funcIO_->ensembleEvaluate(nSamp,M,psOUUXValues_,
                                  nOutputs,psOUUSamOutputs_,funcID);
           for (ss = 0; ss < nSamp; ss++)
             for (ii = 0; ii < M1; ii++)
               GValues[ss*M1+ii] = psOUUSamOutputs_[ss*M1+ii+1];
           odata->numFuncEvals_ += nSamp;
         }
         int failCnt=0;
         for (ss = 0; ss < nSamp; ss++)
           if (psOUUSamOutputs_[ss] >= 0.98*PSUADE_UNDEFINED) failCnt++;
         if (failCnt != 0)
           printf("WARNING: there are %d failed runs out of %d\n",
                  failCnt,nSamp);
         if (failCnt == nSamp || (failCnt > 0 && psOUUSamProbs_ != NULL))
         {
           printf("ERROR: X3 sample cannot admit failures => terminate.\n");
           exit(1);
         }
         if (psOUUMode_ == 1 || psOUUMode_ == 2)
         {
           if (psOUUZ3nSamples_ == 1)
           {
             mean = 0.0;
             for (ss = 0; ss < nSamp; ss++)
               if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED)
                 mean += psOUUSamOutputs_[ss] / 
                                  (double) (nSamp-failCnt);
             FValue = mean;
             for (ii = 0; ii < M1; ii++)
             {
               means[ii] = 0.0;
               for (ss = 0; ss < nSamp; ss++)
                 if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED)
                   means[ii] += GValues[ss*M1+ii] / 
                                    (double) (nSamp-failCnt);
               GValues[ii] = means[ii];
             }
           }
           else
           {
             mean = 0.0;
             for (ss = 0; ss < nSamp; ss++)
             {
               index = ss / psOUUZ4nSamples_;
               mean += psOUUSamOutputs_[ss] / psOUUZ4nSamples_ *
                         psOUUSamProbs_[index];
             }
             FValue = mean;
             for (ii = 0; ii < M1; ii++)
             {
               means[ii] = 0.0;
               for (ss = 0; ss < nSamp; ss++)
               {
                 index = ss / psOUUZ4nSamples_;
                 means[ii] += GValues[ss*M1+ii]/psOUUZ4nSamples_*
                                psOUUSamProbs_[index];
               }
               GValues[ii] = means[ii];
             }
           }
         }
         if (psOUUMode_ == 2 && psOUUStdevMultiplier_ != 0.0)
         {
           stdev = 0.0;
           for (ss = 0; ss < nSamp; ss++)
           {
             index = ss / psOUUZ4nSamples_;
             if (psOUUZ3nSamples_ == 1)
             {
               if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED)
                 stdev += pow(psOUUSamOutputs_[ss]-mean,2.0)/
                              (double) (nSamp - failCnt);
             }
             else
             {
               stdev += pow(psOUUSamOutputs_[ss]-mean, 2.0)*
                        psOUUSamProbs_[index] / (double) psOUUZ4nSamples_;
             }
           }
           FValue = mean + psOUUStdevMultiplier_ * sqrt(stdev);
           for (ii = 0; ii < M1; ii++)
           {
             stdev = 0.0;
             for (ss = 0; ss < nSamp; ss++)
             {
               index = ss / psOUUZ4nSamples_;
               if (psOUUZ3nSamples_ == 1)
               {
                 if (psOUUSamOutputs_[ss] < 0.98*PSUADE_UNDEFINED)
                   stdev += pow(GValues[ss*M1+ii]-means[ii],2.0)/
                                  (double) (nSamp - failCnt);
               }
               else
               {
                 stdev += pow(GValues[ss*M1+ii]-means[ii], 2.0)*
                             psOUUSamProbs_[index] / 
                            (double) psOUUZ4nSamples_;
               }
               GValues[ii] = means[ii] + psOUUStdevMultiplier_ * 
                                sqrt(stdev);
             }
           }
         }
         odata->numFuncEvals_ += nSamp;
         if (isScreenDumpModeOn() && odata->outputLevel_ > 1) 
         {
           printf("OUUOptimizer: Outer optimization iteration = %d\n",
                  odata->numFuncEvals_/nSamp);
           for (ii = 0; ii < nInps; ii++)
             printf("    Current Level 1 input %3d = %e\n", ii+1,
                    XValues[ii]);
           printf("OUUOpt: computing objective (no RS), nFuncEval = %d\n",
                  odata->numFuncEvals_);
           printf("OUUOptimizer: computed  objective (no RS) = %e.\n",
                  FValue);
         }
       }
       else if (*task != NEW_X)
       {
         if (isScreenDumpModeOn() && printLevel > 0)
         {
           for (ii = 0; ii < nInps; ii++)
             printf(" ** Final Input X %4d = %24.16e\n",ii+1,XValues[ii]);
           printf(" ** Final objective function value = %16.8e\n",FValue);
           printf(" ** Total number of iterations     = %ld\n", its);
         }
         break;
       }
     }
     odata->optimalY_ = FValue;
     for (ii = 0; ii < M1; ii++) odata->optimalX_[ii] = XValues[ii];
     delete [] YValues;
     delete [] GValues;
     delete [] nbds;
     delete [] iwork;
     delete [] work;
     delete [] readys;
     delete [] means;
     delete [] lowers;
     delete [] uppers;
#else
     printf("OUUOptimizer ERROR: lbfgs not installed.\n");
     exit(1);
#endif
   }
   else if (optCode_ == 4)
   {
#ifdef HAVE_LINCOA
     int lincoaFlag = 8888;
     kk = nOutputs*(2+M1) + nPts*(4+M1+nPts) + M1*(9+3*M1) + 2*nPts;
     workArray = new double[kk];
     if (psOUULincoaNConstr_ == 0) 
        printf(" ** OUUOptimizer INFO: no constraint declared.\n");
     if (isScreenDumpModeOn()) printf(" ** OUUOptimizer: calling lincoa\n");
     lincoa_(&M1,&nPts,&psOUULincoaNConstr_,psOUULincoaAmat_, 
             &psOUULincoaNConstr_,psOUULincoaBvec_,XValues, &rhobeg, 
             &rhoend, &lincoaFlag,&maxfun, workArray);
     for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
     odata->optimalY_ = psOUUOptimalY_ = workArray[0];
     delete [] workArray;
#else
     printf("OUUOptimizer ERROR: lincoa not installed.\n");
     exit(1);
#endif
   }
   else if (optCode_ == 5)
   {
#ifdef HAVE_NOMAD
     Display dout ( std::cout );
     dout.precision( DISPLAY_PRECISION_STD );
     Parameters nomadp(dout);
     nomadp.set_DISPLAY_DEGREE(0);
     nomadp.set_DIMENSION(M1);
     Point lbnds(M1), ubnds(M1), X0(M1);
     index  = 0;
     for (ii = 0; ii < nInputs; ii++)
     {
       if (psOUUInputTypes_[ii] == psOUUType1)
       {
         lbnds[index] = odata->lowerBounds_[ii];
         ubnds[index] = odata->upperBounds_[ii];
         if (psOUUDesignTypes_ != NULL)
         {
           switch(psOUUDesignTypes_[ii])
           {
             case 0: nomadp.set_BB_INPUT_TYPE (ii, CONTINUOUS);
                     X0[index] = XValues[index];
                     break;
             case 1: nomadp.set_BB_INPUT_TYPE (ii, INTEGER);
                     X0[index] = (int) (XValues[index] + 0.5);
                     break;
           }
         }
         else X0[index] = XValues[index];
         index++;
       }
     }
     nomadp.set_LOWER_BOUND (lbnds);
     nomadp.set_UPPER_BOUND (ubnds);
     nomadp.set_X0(X0);
     psOUUnOutputs_ = nOutputs;
     vector<bb_output_type> bbot (nOutputs);
     bbot[0] = OBJ;
     for (ii = 1; ii < nOutputs; ii++) bbot[ii] = PB;
     nomadp.set_BB_OUTPUT_TYPE ( bbot );
     nomadp.set_SPECULATIVE_SEARCH ( true );
     ddata = 8.0;
     nomadp.set_MESH_UPDATE_BASIS ( ddata );
     int idata = 1;
     nomadp.set_MESH_COARSENING_EXPONENT ( idata );
     nomadp.set_MAX_BB_EVAL (maxfun);
     for (ii = 0; ii < M1; ii++)
     {
        ddata = 0.5 * (ubnds[ii].value() - lbnds[ii].value());
        //nomadp.set_INITIAL_MESH_SIZE ( 0.4 );
        nomadp.set_INITIAL_MESH_SIZE ( ii, ddata, false );
     }
     nomadp.check();
     PsuadeOUUNomadEvaluator madEvaluator(nomadp);
     Mads mads(nomadp, &madEvaluator);
     psOUUObj_ = (void *) odata;
     odata->numFuncEvals_ = 0;
     if (psOUUNomadObjFcnStore_ != NULL) delete [] psOUUNomadObjFcnStore_;
     psOUUNomadObjFcnStore_ = new double[maxfun];
     for (ii = 0; ii < maxfun; ii++) psOUUNomadObjFcnStore_[ii] = -1e35;
     psOUUTolerance_ = odata->tolerance_;
     psOUUStop_ = 0;
     psOUUNomadCnt_ = 0;
     if (isScreenDumpModeOn()) printf(" ** OUUOptimizer: calling nomad\n");
     mads.run();
     delete [] psOUUNomadObjFcnStore_;
     psOUUNomadObjFcnStore_ = NULL;
#else
     printf("OUUOptimizer ERROR: nomad not installed.\n");
     exit(1);
#endif
   }
   else
   {
     printf("OUUOptimizer ERROR: optimizer not recognized.\n");
     exit(1);
   }
   printOutTS(PL_INFO,"*** Optimization Under Uncertainty ends\n");
   printAsterisks(PL_INFO, 0);
   if (isScreenDumpModeOn())
   {
     printf("OUUOptimizer: total number of evaluations = %d\n",
             odata->numFuncEvals_);
   }
   if (optimalX_ != NULL) delete [] optimalX_;
   optimalX_ = new double[nInputs];
   for (ii = 0; ii < nInputs; ii++) optimalX_[ii] = odata->optimalX_[ii];
   optimalY_ = odata->optimalY_;
   numEvals_ = odata->numFuncEvals_;
   
   if (psOUUSaveHistory_ == 1 && psOUUNSaved_ > 0)
   {
     fp = fopen("psuade_ouu_history","w");
     if (fp != NULL)
     {
       for (ii = 0; ii < psOUUNSaved_; ii++)
       {
         fprintf(fp, "999 %d ", nInputs);
         for (kk = 0; kk < nInputs; kk++)
           fprintf(fp, "%24.16e ", psOUUSaveX_[ii*nInputs+kk]);
         for (kk = 0; kk < nOutputs; kk++)
           fprintf(fp, "%24.16e\n", psOUUSaveY_[ii*nOutputs+kk]);
       }
       fclose(fp);
     }
     printf("OUUOptimizer: history saved in psuade_ouu_history\n");
   }
   if (psOUUParallel_ == 1 && odata->funcIO_ != NULL)
     odata->funcIO_->setSynchronousMode();

   if ((odata->setOptDriver_ & 2) && currDriver >= 0 &&
       odata->funcIO_ != NULL)
   {
      odata->funcIO_->setDriver(currDriver);
   }
   delete [] XValues;
   if (psOUUMasterMode_ != 0) psMasterMode_ = 1;

   psOUUnOutputs_ = nOutputs;
   repeatFlag_ = 1;
}

// ************************************************************************
// set local optimizer
// ------------------------------------------------------------------------
void OUUOptimizer::setLocalOptimizer(int optcode)
{
   int  nConstr, ii, kk, M;
   char winput[5000];
   FILE *fp;
   if      (optcode == 0) optCode_ = 0; 
   else if (optcode == 1) optCode_ = 1; 
   else if (optcode == 2) optCode_ = 2; 
   else if (optcode == 3) optCode_ = 3; 
   else if (optcode == 4)
   {
      if (psOUULincoaAmat_ != NULL) delete [] psOUULincoaAmat_;
      if (psOUULincoaBvec_ != NULL) delete [] psOUULincoaBvec_;
      psOUULincoaNConstr_ = 0;
      psOUULincoaAmat_ = NULL;
      psOUULincoaBvec_ = NULL;
      fp = fopen(psOUULincoaConstraintFile_, "r");
      if (fp == NULL)
      {
         printf("OUU WARNING: LINCOA constraints file not found.\n");
         printf("    Since there is no constraint, set default to bobyqa.\n");
      }
      else
      {
         fgets(winput, 5000, fp);
         if (!strcmp(winput,"PSUADE_BEGIN")) fgets(winput, 5000, fp);
         fscanf(fp,"%d %d", &nConstr, &kk);
         if (nConstr <= 0)
         {
            printf("OUU ERROR: nConstr <= 0 in psuade_lincoa_constraints.\n");
            fclose(fp);
            exit(1);
         }
         printf("OUU lincoa: number of constraints read = %d\n", nConstr);
         M = psOUUM1_ + psOUUM2_ + psOUUM3_ + psOUUM4_;
         if (kk != M)
         {
            printf("OUU Lincoa ERROR: nInputs do not match (%d %d).\n",kk,M);
            fclose(fp);
            exit(1);
         }
         psOUULincoaNConstr_ = nConstr;
         psOUULincoaAmat_ = new double[nConstr*M];
         psOUULincoaBvec_ = new double[nConstr];
         kk = 0;
         while (kk < nConstr && feof(fp) == 0)
         {
            fscanf(fp, "%d", &ii);
            if (ii != kk+1)
            {
               printOutTS(PL_ERROR,
                    "OUU Lincoa ERROR: in reading constraint %d .\n",kk+1);
               fclose(fp);
               return;
            }
            for (ii = 0; ii < M; ii++)
               fscanf(fp, "%lg", &psOUULincoaAmat_[kk*M+ii]);
            fscanf(fp, "%lg", &psOUULincoaBvec_[kk]);
            printf("Constr %4d: ",kk+1);
            for (ii = 0; ii < M; ii++)
               printf(" %12.4e ",psOUULincoaAmat_[kk*M+ii]);
            printf("%12.4e\n", psOUULincoaBvec_[kk]);
            kk++;
         }
         fclose(fp);
         if (kk != nConstr)
         {
            printOutTS(PL_ERROR,"OUU Lincoa ERROR: no. of constraints do ");
            printOutTS(PL_ERROR,"not match (%d vs %d)\n", kk, nConstr);
            exit(1);
         }
         optCode_ = 4; 
      }
   }
   else if (optcode == 5) optCode_ = 5; 
   else
   {
      printOutTS(PL_ERROR,"ERROR: invalid optimization code.\n");
      exit(1);
   }
   repeatFlag_ = 0;
   if (psOUUZ3SamInputs_ != NULL) delete [] psOUUZ3SamInputs_;
   if (psOUUZ4SamInputs_ != NULL) delete [] psOUUZ4SamInputs_;
   if (psOUUZ4LBounds_   != NULL) delete [] psOUUZ4LBounds_;
   if (psOUUZ4UBounds_   != NULL) delete [] psOUUZ4UBounds_;
   if (psOUUSamOutputs_  != NULL) delete [] psOUUSamOutputs_;
   if (psOUUSamProbs_    != NULL) delete [] psOUUSamProbs_;
   if (psOUUSamConstrNorms_  != NULL) delete [] psOUUSamConstrNorms_;
   if (psOUULargeSamInputs_  != NULL) delete [] psOUULargeSamInputs_;
   if (psOUULargeSamOutputs_ != NULL) delete [] psOUULargeSamOutputs_;
   if (psOUUXValues_  != NULL) delete [] psOUUXValues_;
   if (psOUUWValues_  != NULL) delete [] psOUUWValues_;
   if (psOUUOptimalX_ != NULL) delete [] psOUUOptimalX_;
   if (psOUUfaPtr_ != NULL) delete psOUUfaPtr_;
   psOUUZ3SamInputs_ = NULL;
   psOUUZ4SamInputs_ = NULL;
   psOUUZ4LBounds_ = NULL;
   psOUUZ4UBounds_ = NULL;
   psOUUSamOutputs_ = NULL;
   psOUUSamConstrNorms_ = NULL;
   psOUUSamProbs_ = NULL;
   psOUULargeSamInputs_ = NULL;
   psOUULargeSamOutputs_ = NULL;
   psOUUXValues_ = NULL;
   psOUUWValues_ = NULL;
   psOUUOptimalX_ = NULL;
   psOUUfaPtr_ = NULL;
   if (!isLibraryModeOn()) 
   {
      if (psOUUInputTypes_ != NULL) delete psOUUInputTypes_;
      psOUUInputTypes_ = NULL;
   }
}

// ************************************************************************
// set discrete inputs 
// ------------------------------------------------------------------------
int OUUOptimizer::setDiscreteInputs(int nInputs, int num, int *list)
{
   int  ii, kk, ignore=0;
   char *cString, pString[1000]; 
   if (psConfig_ != NULL)
   {
      kk = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         sprintf(pString,"iDiscrete%d", ii+1);
         cString = psConfig_->getParameter(pString);
         if (cString != NULL) kk++;
      }
      if (kk > 0)
      {
         printf("ERROR: psConfig already has discrete inputs defined.\n");
         printf("       Duplicates not allowed.\n");
         ignore = 1;
      }
   }
   if (ignore == 0)
   {
      if (psOUUDesignTypes_ == NULL) 
      {
         psOUUDesignTypes_ = new int[nInputs];
         for (ii = 0; ii < nInputs; ii++) psOUUDesignTypes_[ii] = 0;
      }
      for (ii = 0; ii < num; ii++)
      {
         if (list[ii] < 1 || list[ii] > nInputs)
         {
            printf("setDiscreteInputs ERROR: input not in range.\n");
            return -1;
         }
         psOUUDesignTypes_[list[ii]-1] = 1;
      }
   }
   return 0;
}

// ************************************************************************
// set variable type 
// ------------------------------------------------------------------------
int OUUOptimizer::setVariableType(int num, int vtype)
{
   int  ii, kk, ignore=0, *tmpInds=NULL;
   char *cString, pString[1000]; 

   if (num <= 0)
   {
      printf("setVariabletype ERROR: invalid argument.\n");
      printf("                Variable number %d must be > 0\n", num);
      exit(1);
   } 
   if (vtype <= 0 || vtype > 5)
   {
      printf("setVariabletype ERROR: invalid argument.\n");
      printf("                Variable type %d must be in [1,5]\n",vtype);
      exit(1);
   }
   if (psConfig_ != NULL && vtype == OUUDesignDisc)
   {
      sprintf(pString,"iDiscrete%d", num);
      psConfig_->putParameter(pString);
   }
   if (psOUUInputTypes_ == NULL) 
   {
      psOUUInputTypes_ = new int[2*num+1];
      psOUUM_ = 2 * num + 1;
      for (ii = 0; ii < psOUUM_; ii++) psOUUInputTypes_[ii] = 0;
   }
   if (psOUUInputTypes_ != NULL && num > psOUUM_)
   {
      tmpInds = psOUUInputTypes_;
      psOUUInputTypes_ = new int[2*num+1];
      for (ii = 0; ii < psOUUM_; ii++) psOUUInputTypes_[ii] = tmpInds[ii];
      for (ii = psOUUM_; ii < 2*num+1; ii++) psOUUInputTypes_[ii] = 0;
      psOUUM_ = 2 * num + 1;
   } 
   psOUUInputTypes_[num-1] = vtype;
   if (tmpInds != NULL) delete [] tmpInds;
   if (psOUUM1_ < 0) psOUUM1_ = psOUUM2_ = psOUUM3_ = psOUUM4_ = 0;
   switch (vtype)
   {
      case 1: psOUUM1_++; break;
      case 2: psOUUM2_++; break;
      case 3: psOUUM3_++; break;
      case 4: psOUUM4_++; break;
      case 5: psOUUM1_++; break;
   }
   return 0;
}

// ************************************************************************
// set objective function
// ------------------------------------------------------------------------
void OUUOptimizer::setConstraintFile(char *fname)
{
   int    ii, kk, M, nConstr;
   char   winput[5001];
   FILE   *fp;

   strcpy(psOUULincoaConstraintFile_, fname);
   if (psOUULincoaAmat_ != NULL) delete [] psOUULincoaAmat_;
   if (psOUULincoaBvec_ != NULL) delete [] psOUULincoaBvec_;
   psOUULincoaNConstr_ = 0;
   psOUULincoaAmat_ = NULL;
   psOUULincoaBvec_ = NULL;
   fp = fopen(psOUULincoaConstraintFile_, "r");
   if (fp == NULL)
   {
      printf("OUU WARNING: LINCOA constraints file not found.\n");
      printf("    Since there is no constraint, set default to bobyqa.\n");
   }
   else
   {
      fgets(winput, 5000, fp);
      if (!strcmp(winput,"PSUADE_BEGIN")) fgets(winput, 5000, fp);
      fscanf(fp,"%d %d", &nConstr, &kk);
      if (nConstr <= 0)
      {
         printf("OUU ERROR: nConstr <= 0 in psuade_lincoa_constraints.\n");
         fclose(fp);
         exit(1);
      }
      printf("OUU lincoa: number of constraints read = %d\n", nConstr);
      M = psOUUM1_ + psOUUM2_ + psOUUM3_ + psOUUM4_;
      if (kk != M)
      {
         printf("OUU Lincoa ERROR: nInputs do not match (%d %d).\n",kk,M);
         fclose(fp);
         exit(1);
      }
      psOUULincoaNConstr_ = nConstr;
      psOUULincoaAmat_ = new double[nConstr*M];
      psOUULincoaBvec_ = new double[nConstr];
      kk = 0;
      while (kk < nConstr && feof(fp) == 0)
      {
         fscanf(fp, "%d", &ii);
         if (ii != kk+1)
         {
            printOutTS(PL_ERROR,
               "OUU Lincoa ERROR: in reading constraint %d .\n",kk+1);
            fclose(fp);
            exit(1);
         }
         for (ii = 0; ii < M; ii++)
            fscanf(fp, "%lg", &psOUULincoaAmat_[kk*M+ii]);
         fscanf(fp, "%lg", &psOUULincoaBvec_[kk]);
         printf("Constr %4d: ",kk+1);
         for (ii = 0; ii < M; ii++)
            printf(" %12.4e ",psOUULincoaAmat_[kk*M+ii]);
         printf("%12.4e\n", psOUULincoaBvec_[kk]);
         kk++;
      }
      fclose(fp);
      if (kk != nConstr)
      {
         printOutTS(PL_ERROR,"OUU Lincoa ERROR: no. of constraints do ");
         printOutTS(PL_ERROR,"not match (%d vs %d)\n", kk, nConstr);
         exit(1);
      }
      optCode_ = 4; 
   }
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
void OUUOptimizer::setDiscreteVariable(int index)
{
   char winput[1000];
   if (index <= 0)
   {
      printf("OUU setDiscreteVariable ERROR: variable number <= 0.\n");
      exit(1);
   }
   if (psConfig_ == NULL) psConfig_ = new PsuadeConfig();
   sprintf(winput,"iDiscrete%d", index);
   psConfig_->putParameter(winput);
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
OUUOptimizer& OUUOptimizer::operator=(const OUUOptimizer &)
{
   printf("OUUOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************ 
// validate
// ------------------------------------------------------------------------
extern "C" 
{
void validate(int nSamps,double *samInputs,double *samOutputs,double *errors)
{
   int        nInps, nSubSamples, ii, ss, *iArray;
   double     *X, *Y;
   char       **targv;
   FuncApprox *fa;

   for (ii = 0; ii < 3; ii++) errors[ii] = 1e12;
   nInps = psOUUM4_;
   if (nSamps <= 20) nSubSamples = 1;
   else
   {
      nSubSamples = nSamps / 10;
      if (nSubSamples*10 != nSamps) nSubSamples++;
   }
   int rsState = psRSExpertMode_;
   int saveInteractive = psInteractive_;
   psRSExpertMode_ = 0;
   psInteractive_ = 0;
   fa = genFA(psOUUZ4RSType_, nInps, -1, nSamps-nSubSamples);
   if (fa == NULL)
   {
      printOutTS(PL_ERROR,"ERROR: cannot create function approximator.\n");
      return;
   }
   fa->setBounds(psOUUZ4LBounds_,psOUUZ4UBounds_);
   fa->setOutputLevel(0);
   if (psOUUZ4RSType_ == PSUADE_RS_KR) 
   {
      targv = new char*[1];
      targv[0] = new char[100];
      if (psOUUZ4RSAux_ == 0) strcpy(targv[0], "setMode2");
      else                    strcpy(targv[0], "setMode3");
      fa->setParams(1, targv);
      delete [] targv[0];
      delete [] targv;
   }
   iArray = new int[nSamps];
   if (nSubSamples > 1) generateRandomIvector(nSamps, iArray);
   else
   {
      for (ss = 0; ss < nSamps; ss++) iArray[ss] = ss;
   }
   X = new double[nSamps*nInps];
   for (ii = 0; ii < nInps; ii++)
   {
      for (ss = 0; ss < nSamps; ss++)
         X[iArray[ss]*nInps+ii] = samInputs[ss*nInps+ii];
   }
   Y = new double[nSamps];
   for (ss = 0; ss < nSamps; ss++) Y[iArray[ss]] = samOutputs[ss];
   oData *odata = (oData *) psOUUObj_;
   if (odata != NULL && odata->outputLevel_ > 2)
      printOutTS(PL_INFO, "Cross validation begins ...\n");
   int    status, ss2, count;
   double cvErr1s,cvErr2s,cvMaxs,CVErr1s=0,CVErr2s=0,CVMaxs=0,ddata;
   double *YT = new double[nSamps];
   double *X2 = new double[nSamps*nInps];
   double *Y2 = new double[nSamps];
   for (ss = 0; ss < nSamps; ss+=nSubSamples)
   {
      if (odata != NULL && odata->outputLevel_ > 2)
         printOutTS(PL_INFO, "Cross validation group %d (of 10)\n",
                    ss/nSubSamples+1);
      count = 0;
      for (ss2 = 0; ss2 < nSamps; ss2++)
      {
         if (ss2 < ss || ss2 >= (ss+nSubSamples))
         {
            for (ii = 0; ii < nInps; ii++)
               X2[count*nInps+ii] = X[ss2*nInps+ii];
            Y2[count++] = Y[ss2];
         }
      }
      status = fa->initialize(X2, Y2);
      if (status == -1) break;
      count = nSubSamples;
      if ((ss + nSubSamples) > nSamps) count = nSamps - ss;
      fa->evaluatePoint(count, &(X[ss*nInps]), YT);
      cvErr1s = cvErr2s = cvMaxs = 0.0;
      for (ss2 = 0; ss2 < count; ss2++)
      {
         ddata = YT[ss2] - Y[ss+ss2];
         if (Y[ss+ss2] != 0.0) ddata = ddata / PABS(Y[ss+ss2]);
         cvErr1s += ddata;
         cvErr2s += (ddata * ddata);
         if (PABS(ddata) > cvMaxs) cvMaxs = PABS(ddata);
      }
      CVErr1s += cvErr1s;
      CVErr2s += cvErr2s;
      if (cvMaxs > CVMaxs ) CVMaxs = cvMaxs; 
   }
   if (status >= 0)
   {
      errors[0] = CVErr1s / (double) nSamps;
      errors[1] = sqrt(CVErr2s / nSamps);
      errors[2] = CVMaxs;
   }
   psRSExpertMode_ = rsState;
   psInteractive_ = saveInteractive; 
   delete fa;
   delete [] X;
   delete [] Y;
   delete [] YT;
   delete [] X2;
   delete [] Y2;
   delete [] iArray;
}
}

#ifdef HAVE_NOMAD
// ************************************************************************
// PsuadeNomadEvaluator constructor
// ------------------------------------------------------------------------
PsuadeOUUNomadEvaluator::PsuadeOUUNomadEvaluator(const Parameters &pset):
                         Evaluator(pset)
{
}

// ************************************************************************
// PsuadeOUUNomadEvaluator destructor
// ------------------------------------------------------------------------
PsuadeOUUNomadEvaluator::~PsuadeOUUNomadEvaluator()
{
}

// ************************************************************************
// PsuadeNomadEvaluator evaluation 
// ------------------------------------------------------------------------
bool PsuadeOUUNomadEvaluator::eval_x(Eval_Point &X,const Double &h_max,
                                     bool &count_eval)
{
   int    ii, jj, M, index, funcID, nInps = X.get_n(), nOuts, goodFlag;
   double *XVals, *YVals = new double[psOUUnOutputs_], ddata;
   bool   retFlag=true;
   FILE   *fp=NULL;
   oData  *odata;

   if (psOUUStop_ == 1)
   {
      YVals[0] = -1e50;
      X.set_bb_output(0, YVals[0]);
      YVals[0] = 0.0;
      for (ii = 1; ii < psOUUnOutputs_; ii++)
         X.set_bb_output(ii, YVals[0]);
      retFlag = false;
      return retFlag;
   }
   XVals = new double[nInps];
   for (ii = 0; ii < nInps; ii++) XVals[ii] = X[ii].value();
   nOuts = psOUUnOutputs_;
   ouuevalfunccobyla(nInps, nOuts-1, XVals, YVals, &YVals[1], NULL);

   for (ii = 0; ii < psOUUnOutputs_; ii++) X.set_bb_output(ii, YVals[ii]);

   odata = (oData *) psOUUObj_;
   funcID = psOUUNomadCnt_;
   psOUUNomadCnt_++;
   if (YVals[0] < odata->optimalY_)
   {
      goodFlag = 1;
      for (ii = 1; ii < psOUUnOutputs_; ii++)
         if (YVals[ii] > 0) goodFlag = 0;
      if (goodFlag == 1)
      {
         for (ii = 0; ii < nInps; ii++)
            odata->optimalX_[ii] = XVals[ii];
         odata->optimalY_ = YVals[0];
         if (odata->outputLevel_ > 0)
            printf("   OUUNomad current best Y = %e (nfevals=%d)\n",
                   YVals[0], funcID+1);
      }
      //else
      //{
      //   if (odata->outputLevel_ > 1)
      //   {
      //      printf("OUUNomad INFO: ignore infeasible solution:\n");
      //      for (ii = 0; ii < nInps; ii++)
      //      printf("   X[%5d] = %e\n", ii+1, XVals[ii]);
      //      for (ii = 0; ii < psOUUnOutputs_; ii++)
      //         printf("   Y[%5d] = %e\n", ii+1, YVals[ii]);
      //   }
      //}
   }

   psOUUNomadObjFcnStore_[funcID] = YVals[0];
   int    converged, nn=20;
   double mean, stdev, dmax, track[30];
   if (funcID >= nn)
   {
      dmax = 0.0;
      for (ii = 0; ii <= funcID; ii++)
      {
         ddata = psOUUNomadObjFcnStore_[ii];
         if (ddata < 0) ddata = - ddata;
         if (ddata > dmax) dmax = ddata;
      }
      if (dmax == 0) dmax = 1;
      jj = 1;
      track[0] = psOUUNomadObjFcnStore_[funcID-nn+1];
      for (ii = funcID-nn+2; ii <= funcID; ii++)
      {
         if (psOUUNomadObjFcnStore_[ii] < track[jj-1])
            track[jj++] = psOUUNomadObjFcnStore_[ii];
      }
      if (jj >= 5)
      {   
         for (ii = 0; ii < 5; ii++) track[ii] = track[jj-5+ii];
         jj = 5;
         mean = 0.0;
         for (ii = 0; ii < jj; ii++) mean += track[ii];
         mean /= (double) jj;
         stdev = 0.0;
         for (ii = 0; ii < jj; ii++) stdev += pow(track[ii]-mean,2.0);
         stdev /= (double) (jj - 1);
         stdev = sqrt(stdev);
      }
      else stdev = 100 * dmax;
      converged = 0;
      if (odata->outputLevel_ > 1)
         printf("Convergence check: %e <? %e\n",stdev/dmax,psOUUTolerance_);
      if (stdev / dmax < psOUUTolerance_) 
      {
         converged = 1;
         printf("INFO: convergence detected at iteration %d ==> stop.\n",
                funcID);
         psOUUStop_ = 1;
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
      psOUUStop_ = 1;
   }
   return retFlag;
}
#endif

