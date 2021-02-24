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
// Functions for the class TGP
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_TGP
extern "C"
{
#include "matrix.h"
#include "rand_draws.h"
#include "rhelp.h"
#include "predict.h"
}
#include "model.h"
#include "params.h"
#include "mstructs.h"
#endif

#include "TBGP.h"
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class TGP
// ------------------------------------------------------------------------
TGP::TGP(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   faID_ = PSUADE_RS_TGP;
#ifdef HAVE_TGP
   tgp_ = NULL;
#endif
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
TGP::~TGP()
{
#ifdef HAVE_TGP
   if (tgp_ != NULL) delete tgp_;
#endif
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int TGP::initialize(double *X, double *Y)
{
   int status=-999;
   genNDGridData(X, Y, &status, NULL, NULL);
   if (psRSCodeGen_ == 1) 
      printf("TGP INFO: response surface stand-alone code not available.\n");
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int TGP::genNDGridData(double *X, double *Y, int *N, double **X2, 
                       double **Y2)
{
#ifdef HAVE_TGP
   int    BTE[3], R=1, linBurn=1, PS_FALSE=0, PS_TRUE=1;
   int    totPts, ii, ss, count, *stateIn, bte0, bte1;
   double *HX, *XX, *Xloc, *dparams, *itemps, *gpcs, *ZZMean;
   void   *tgp_state=NULL;
   char   lineOut[1000];

   bte0 = 2000;
   bte1 = 7000;
   if (psRSExpertMode_ == 1)
   {
      sprintf(lineOut, "TGP: enter burn-in sample size (e.g. 500-5000): ");
      bte0 = getInt(500, 5000, lineOut); 
      sprintf(lineOut,"TGP: enter MCMC sample size (e.g. %d-20000): ",2000+bte0);
      bte1 = getInt(2000+bte0, 20000, lineOut); 
   }

   if (outputLevel_ >= 1) printf("TGP training begins....\n");
   if (tgp_ != NULL) delete tgp_;

   stateIn = new int[3];
   stateIn[0] = 782;
   stateIn[1] = 267;
   stateIn[2] = 218;
   unsigned int lstate = three2lstate(stateIn);
   tgp_state = newRNGstate(lstate);
   delete [] stateIn;

   BTE[0] = bte0;
   BTE[1] = bte1;
   BTE[2] = 2;
   dparams = new double[(nInputs_+1)*(nInputs_+1)+nInputs_+45];
   checkAllocate(dparams, "dparams in TGP::genNDGridData");
   dparams[0] = 0.5;  // tree prior alpha
   dparams[1] = 2.0;  // tree prior beta
   dparams[2] = 10.0; // tree prior minpart
   if ((nInputs_ + 2) > 10) dparams[2] = 1.0 * (nInputs_ + 2.0);
   dparams[3] = 1.0;  // tree prior splitmin
   dparams[4] = 1.0 * nInputs_;  // tree prior basemax
   dparams[5] = 0.0;  // meanfn:  0 - linear, 1 : constant
   dparams[6] = 2.0;  // bflat in {b0, bmle, bflat, b0not,..}
   for (ii = 0; ii <= nInputs_; ii++) dparams[7+ii] = 0.0;
   count = nInputs_ + 8;
   dparams[count] = 1.0;  // Wi (1, (nInputs+1)^2 zeros)
   for (ii = 0; ii < (nInputs_+1)*(nInputs_+1)-1; ii++)
      dparams[count+1+ii] = 0.0;
   for (ii = 1; ii < nInputs_; ii++) dparams[count+ii*(nInputs_+1)] = 1.0;
   count += (nInputs_ + 1) * (nInputs_ + 1);
   dparams[count++] = 1.0;   // s2tau2 = c(1,1)
   dparams[count++] = 1.0;
   dparams[count++] = 5.0;   // s2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // s2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 5.0;   // tau2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // tau2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 0.1;
   dparams[count++] = 1;     // "expsep"
   dparams[count++] = 0.1;   // gd (1)
   dparams[count++] = 0.5;   // gd (2)
   dparams[count++] = 1;     // nug.p (1)
   dparams[count++] = 1;     // nug.p (2)
   dparams[count++] = 1;     // nug.p (3)
   dparams[count++] = 1;     // nug.p (4)
   dparams[count++] = -1;    // nug.lam (1)
   dparams[count++] = -1;    // nug.lam (2)
   dparams[count++] = -1;    // nug.lam (3)
   dparams[count++] = -1;    // nug.lam (4)
   dparams[count++] = 10.0;  // gamma(1)
   dparams[count++] = 0.2;   // gamma(2)
   dparams[count++] = 0.7;   // gamma(3)
   dparams[count++] = 1.0;   // d.p(1)
   dparams[count++] = 20.0;  // d.p(2)
   dparams[count++] = 10.0;  // d.p(3)
   dparams[count++] = 10.0;  // d.p(4)
   dparams[count++] = -1;    // d.lam (1)
   dparams[count++] = -1;    // d.lam (2)
   dparams[count++] = -1;    // d.lam (3)
   dparams[count++] = -1;    // d.lam (4)
   dparams[count++] = 0;     // nu?
   itemps = new double[7];
   itemps[0] = 1;
   itemps[1] = 0;
   itemps[2] = 0;
   itemps[3] = 1;
   itemps[4] = 1;
   itemps[5] = 0;
   itemps[6] = 1;

   if ((*N) != -999 && X2 != NULL && Y2 != NULL)
   {
      totPts = nPtsPerDim_;
      for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
      HX = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++) 
         HX[ii] = (upperBounds_[ii] - lowerBounds_[ii]) /
                  (double) (nPtsPerDim_ - 1); 
      XX = new double[totPts*nInputs_];
      Xloc = new double[nInputs_];
      checkAllocate(Xloc, "Xloc in TGP::genNDGridData");
      for (ii = 0; ii < nInputs_; ii++) Xloc[ii] = lowerBounds_[ii];
 
      for (ss = 0; ss < totPts; ss++)
      {
         for (ii = 0; ii < nInputs_; ii++ ) XX[ss*nInputs_+ii] = Xloc[ii];
         for (ii = 0; ii < nInputs_; ii++ ) 
         {
            Xloc[ii] += HX[ii];
            if (Xloc[ii] < upperBounds_[ii] || 
                PABS(Xloc[ii] - upperBounds_[ii]) < 1.0E-7) break;
            else Xloc[ii] = lowerBounds_[ii];
         }
      }
      delete [] HX;
      delete [] Xloc;
      tgp_ = new Tgp((void*) tgp_state, (int) nSamples_, (int) nInputs_, 
                     (int) totPts, (int) BTE[0], (int) BTE[1], (int) BTE[2], 
                     (int) R, (int) linBurn, (bool) PS_TRUE, (bool) PS_TRUE, 
                     (bool) PS_FALSE, (int) PS_FALSE, (bool) PS_FALSE, 
                     X, Y, XX, dparams, itemps, (bool) PS_FALSE, 
                     (int) PS_FALSE, (double *) NULL, (double *) NULL);
   }
   else
   {
      totPts = nSamples_;
      XX = new double[totPts*nInputs_];
      for (ss = 0; ss < totPts*nInputs_; ss++) XX[ss] = X[ss];
      tgp_ = new Tgp((void*) tgp_state, (int) nSamples_, (int) nInputs_, 
                     (int) totPts, (int) BTE[0], (int) BTE[1], (int) BTE[2], 
                     (int) R, (int) linBurn, (bool) PS_TRUE, (bool) PS_FALSE, 
                     (bool) PS_FALSE, (int) PS_FALSE, (bool) PS_FALSE, 
                     X, Y, XX, dparams, itemps, (bool) PS_FALSE, 
                     (int) PS_FALSE, (double *) NULL, (double *) NULL);
   }
   tgp_->Init();

   tgp_->Rounds();

#if 0
   ZZMean = NULL;
   if ((*N) != -999 && X2 != NULL && Y2 != NULL)
   {
      essOut = new double[3];  // itemps[1]*2+1
      ZpMean = new double[nSamples_];
      ZpQ = new double[nSamples_];
      ZpS2 = new double[nSamples_];
      ZpQ1 = new double[nSamples_];
      ZpMedian = new double[nSamples_];
      ZpQ2 = new double[nSamples_];
      ZZMean = new double[totPts];
      ZZQ = new double[totPts];
      ZZS2 = new double[totPts];
      ZZQ1 = new double[totPts];
      ZZMedian = new double[totPts];
      ZZQ2 = new double[totPts];
      tgp_->GetStats(iZero,ZpMean,ZZMean,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                     ZpS2,ZZS2,NULL,NULL,NULL,ZpQ1,ZpMedian,
                     ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);

      delete [] essOut;
      delete [] ZpMean;
      delete [] ZpQ;
      delete [] ZpS2;
      delete [] ZpQ1;
      delete [] ZpMedian;
      delete [] ZpQ2;
      delete [] ZZQ;
      delete [] ZZS2;
      delete [] ZZQ1;
      delete [] ZZMedian;
      delete [] ZZQ2;
   }
#else
   ZZMean = new double[totPts];
   evaluatePoint(totPts, XX, ZZMean);
#endif

   tgp_->GetPseudoPrior(itemps);
   delete [] dparams;
   delete [] itemps;

   gpcs = new double[4];
   tgp_->GetTreeStats(gpcs);
   delete [] gpcs;

   deleteRNGstate(tgp_state);
   tgp_state = NULL;

   if (outputLevel_ >= 1) printf("TGP training ends.\n");
   if ((*N) == -999 || X2 == NULL || Y2 == NULL) return 0;
   (*N) = totPts;
   (*X2) = XX;
   (*Y2) = ZZMean;
#else
   printf("PSUADE ERROR : TGP not installed.\n");
   (*N) = 0;
   (*X2) = NULL;
   (*Y2) = NULL;
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int TGP::gen1DGridData(double *X, double *Y, int ind1,
                      double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_TGP
   int    BTE[3], R=1, linBurn=1, PS_FALSE=0, PS_TRUE=1;
   int    totPts, ii, kk, count, *stateIn, bte0, bte1;
   double HX, *XX, *dparams, *itemps, *gpcs, *ZZMean;
   void   *tgp_state=NULL;
   char   lineOut[1000];

   bte0 = 2000;
   bte1 = 7000;
   if (psRSExpertMode_ == 1)
   {
      sprintf(lineOut, "TGP: enter burn-in sample size (e.g. 500-5000): ");
      bte0 = getInt(500, 5000, lineOut); 
      sprintf(lineOut,"TGP: enter MCMC sample size (e.g. %d-20000): ",2000+bte0);
      bte1 = getInt(2000+bte0, 20000, lineOut); 
   }

   if (outputLevel_ >= 1) printf("TGP training begins....\n");
   if (tgp_ != NULL) delete tgp_;

   stateIn = new int[3];
   stateIn[0] = 782;
   stateIn[1] = 267;
   stateIn[2] = 218;
   unsigned int lstate = three2lstate(stateIn);
   tgp_state = newRNGstate(lstate);
   delete [] stateIn;

   BTE[0] = bte0;
   BTE[1] = bte1;
   BTE[2] = 2;
   dparams = new double[(nInputs_+1)*(nInputs_+1)+nInputs_+45];
   checkAllocate(dparams, "dparams in TGP::genNDGridData");
   dparams[0] = 0.5;  // tree prior alpha
   dparams[1] = 2.0;  // tree prior beta
   dparams[2] = 10.0; // tree prior minpart
   if ((nInputs_ + 2) > 10) dparams[2] = 1.0 * (nInputs_ + 2.0);
   dparams[3] = 1.0;  // tree prior splitmin
   dparams[4] = 1.0 * nInputs_;  // tree prior basemax
   dparams[5] = 0.0;  // meanfn:  0 - linear, 1 : constant
   dparams[6] = 2.0;  // bflat in {b0, bmle, bflat, b0not,..}
   for (ii = 0; ii <= nInputs_; ii++) dparams[7+ii] = 0.0;
   count = nInputs_ + 8;
   dparams[count] = 1.0;  // Wi (1, (nInputs+1)^2 zeros)
   for (ii = 0; ii < (nInputs_+1)*(nInputs_+1)-1; ii++)
      dparams[count+1+ii] = 0.0;
   for (ii = 1; ii < nInputs_; ii++) dparams[count+ii*(nInputs_+1)] = 1.0;
   count += (nInputs_ + 1) * (nInputs_ + 1);
   dparams[count++] = 1.0;   // s2tau2 = c(1,1)
   dparams[count++] = 1.0;
   dparams[count++] = 5.0;   // s2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // s2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 5.0;   // tau2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // tau2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 0.1;
   dparams[count++] = 1;     // "expsep"
   dparams[count++] = 0.1;   // gd (1)
   dparams[count++] = 0.5;   // gd (2)
   dparams[count++] = 1;     // nug.p (1)
   dparams[count++] = 1;     // nug.p (2)
   dparams[count++] = 1;     // nug.p (3)
   dparams[count++] = 1;     // nug.p (4)
   dparams[count++] = -1;    // nug.lam (1)
   dparams[count++] = -1;    // nug.lam (2)
   dparams[count++] = -1;    // nug.lam (3)
   dparams[count++] = -1;    // nug.lam (4)
   dparams[count++] = 10.0;  // gamma(1)
   dparams[count++] = 0.2;   // gamma(2)
   dparams[count++] = 0.7;   // gamma(3)
   dparams[count++] = 1.0;   // d.p(1)
   dparams[count++] = 20.0;  // d.p(2)
   dparams[count++] = 10.0;  // d.p(3)
   dparams[count++] = 10.0;  // d.p(4)
   dparams[count++] = -1;    // d.lam (1)
   dparams[count++] = -1;    // d.lam (2)
   dparams[count++] = -1;    // d.lam (3)
   dparams[count++] = -1;    // d.lam (4)
   dparams[count++] = 0;     // nu?
   itemps = new double[7];
   itemps[0] = 1;
   itemps[1] = 0;
   itemps[2] = 0;
   itemps[3] = 1;
   itemps[4] = 1;
   itemps[5] = 0;
   itemps[6] = 1;


   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts];
   XX = new double[totPts*nInputs_];
   checkAllocate(XX, "XX in TGP::gen1DGridData");
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (kk = 0; kk < nInputs_; kk++) 
         XX[ii*nInputs_+kk] = settings[kk]; 
      XX[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
      (*X2)[ii] = HX * ii + lowerBounds_[ind1];
   }
    
   tgp_ = new Tgp((void*) tgp_state, (int) nSamples_, (int) nInputs_, 
                  (int) nSamples_, (int) BTE[0], (int) BTE[1], (int) BTE[2], 
                  (int) R, (int) linBurn, (bool) PS_TRUE, (bool) PS_TRUE, 
                  (bool) PS_FALSE, (int) PS_FALSE, (bool) PS_FALSE, 
                  X, Y, X, dparams, itemps, (bool) PS_FALSE, 
                  (int) PS_FALSE, (double *) NULL, (double *) NULL);
   
   tgp_->Init();

   tgp_->Rounds();

   if (outputLevel_ >= 1) printf("TGP interpolation begins....\n");
#if 0
   tgp_->LoadInterpolatedPts(totPts, XX);
   tgp_->Predict();
   essOut = new double[3];
   ZpMean = new double[nSamples_];
   ZpQ = new double[nSamples_];
   ZpS2 = new double[nSamples_];
   ZpQ1 = new double[nSamples_];
   ZpMedian = new double[nSamples_];
   ZpQ2 = new double[nSamples_];
   ZZMean = new double[totPts];
   ZZQ = new double[totPts];
   ZZS2 = new double[totPts];
   ZZQ1 = new double[totPts];
   ZZMedian = new double[totPts];
   ZZQ2 = new double[totPts];
   tgp_->GetStats(iZero,ZpMean,ZZMean,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                  ZpS2,ZZS2,NULL,NULL,NULL,ZpQ1,ZpMedian,
                  ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
   delete [] essOut;
   delete [] ZpMean;
   delete [] ZpQ;
   delete [] ZpS2;
   delete [] ZpQ1;
   delete [] ZpMedian;
   delete [] ZpQ2;
   delete [] ZZQ;
   delete [] ZZS2;
   delete [] ZZQ1;
   delete [] ZZMedian;
   delete [] ZZQ2;
#else
   ZZMean = new double[totPts];
   evaluatePoint(totPts, XX, ZZMean);
#endif

   tgp_->GetPseudoPrior(itemps);
   delete [] itemps;
   delete [] dparams;

   gpcs = new double[4];
   tgp_->GetTreeStats(gpcs);
   delete [] gpcs;

   deleteRNGstate(tgp_state);
   tgp_state = NULL;

   if (outputLevel_ >= 1) printf("TGP interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = ZZMean;
   delete [] XX;
#else
   printf("PSUADE ERROR : TGP not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int TGP::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                      double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_TGP
   int    BTE[3], R=1, linBurn=1, PS_FALSE=0, PS_TRUE=1, bte0, bte1;
   int    totPts, ii, jj, kk, count, *stateIn, index;
   double *HX, *XX, *dparams, *itemps, *gpcs, *ZZMean;
   void   *tgp_state=NULL;
   char   lineOut[1000];

   bte0 = 2000;
   bte1 = 7000;
   if (psRSExpertMode_ == 1)
   {
      sprintf(lineOut, "TGP: enter burn-in sample size (e.g. 500-5000): ");
      bte0 = getInt(500, 5000, lineOut); 
      sprintf(lineOut,"TGP: enter MCMC sample size (e.g. %d-20000): ",2000+bte0);
      bte1 = getInt(2000+bte0, 20000, lineOut); 
   }

   if (outputLevel_ >= 1) printf("TGP training begins....\n");
   if (tgp_ != NULL) delete tgp_;

   stateIn = new int[3];
   stateIn[0] = 782;
   stateIn[1] = 267;
   stateIn[2] = 218;
   unsigned int lstate = three2lstate(stateIn);
   tgp_state = newRNGstate(lstate);
   delete [] stateIn;

   BTE[0] = bte0;
   BTE[1] = bte1;
   BTE[2] = 2;
   dparams = new double[(nInputs_+1)*(nInputs_+1)+nInputs_+45];
   checkAllocate(dparams, "dparams in TGP::gen2DGridData");
   dparams[0] = 0.5;  // tree prior alpha
   dparams[1] = 2.0;  // tree prior beta
   dparams[2] = 10.0; // tree prior minpart
   if ((nInputs_ + 2) > 10) dparams[2] = 1.0 * (nInputs_ + 2.0);
   dparams[3] = 1.0;  // tree prior splitmin
   dparams[4] = 1.0 * nInputs_;  // tree prior basemax
   dparams[5] = 0.0;  // meanfn:  0 - linear, 1 : constant
   dparams[6] = 2.0;  // bflat in {b0, bmle, bflat, b0not,..}
   for (ii = 0; ii <= nInputs_; ii++) dparams[7+ii] = 0.0;
   count = nInputs_ + 8;
   dparams[count] = 1.0;  // Wi (1, (nInputs+1)^2 zeros)
   for (ii = 0; ii < (nInputs_+1)*(nInputs_+1)-1; ii++)
      dparams[count+1+ii] = 0.0;
   for (ii = 1; ii < nInputs_; ii++) dparams[count+ii*(nInputs_+1)] = 1.0;
   count += (nInputs_ + 1) * (nInputs_ + 1);
   dparams[count++] = 1.0;   // s2tau2 = c(1,1)
   dparams[count++] = 1.0;
   dparams[count++] = 5.0;   // s2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // s2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 5.0;   // tau2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // tau2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 0.1;
   dparams[count++] = 1;     // "expsep"
   dparams[count++] = 0.1;   // gd (1)
   dparams[count++] = 0.5;   // gd (2)
   dparams[count++] = 1;     // nug.p (1)
   dparams[count++] = 1;     // nug.p (2)
   dparams[count++] = 1;     // nug.p (3)
   dparams[count++] = 1;     // nug.p (4)
   dparams[count++] = -1;    // nug.lam (1)
   dparams[count++] = -1;    // nug.lam (2)
   dparams[count++] = -1;    // nug.lam (3)
   dparams[count++] = -1;    // nug.lam (4)
   dparams[count++] = 10.0;  // gamma(1)
   dparams[count++] = 0.2;   // gamma(2)
   dparams[count++] = 0.7;   // gamma(3)
   dparams[count++] = 1.0;   // d.p(1)
   dparams[count++] = 20.0;  // d.p(2)
   dparams[count++] = 10.0;  // d.p(3)
   dparams[count++] = 10.0;  // d.p(4)
   dparams[count++] = -1;    // d.lam (1)
   dparams[count++] = -1;    // d.lam (2)
   dparams[count++] = -1;    // d.lam (3)
   dparams[count++] = -1;    // d.lam (4)
   dparams[count++] = 0;     // nu?
   itemps = new double[7];
   itemps[0] = 1;
   itemps[1] = 0;
   itemps[2] = 0;
   itemps[3] = 1;
   itemps[4] = 1;
   itemps[5] = 0;
   itemps[6] = 1;


   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[2*totPts];
   checkAllocate(*X2, "X2 in TGP::gen2DGridData");
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         index = ii * nPtsPerDim_ + jj;
         for (kk = 0; kk < nInputs_; kk++) 
            XX[index*nInputs_+kk] = settings[kk]; 
         XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
         XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
         (*X2)[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         (*X2)[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }
   delete [] HX;
    
   tgp_ = new Tgp((void*) tgp_state, (int) nSamples_, (int) nInputs_, 
                  (int) nSamples_, (int) BTE[0], (int) BTE[1], (int) BTE[2], 
                  (int) R, (int) linBurn, (bool) PS_TRUE, (bool) PS_TRUE, 
                  (bool) PS_FALSE, (int) PS_FALSE, (bool) PS_FALSE, 
                  X, Y, X, dparams, itemps, (bool) PS_FALSE, 
                  (int) PS_FALSE, (double *) NULL, (double *) NULL);

   tgp_->Init();

   tgp_->Rounds();

   if (outputLevel_ >= 1) printf("TGP interpolation begins....\n");
#if 0
   tgp_->LoadInterpolatedPts(totPts, XX);
   tgp_->Predict();
   essOut = new double[3];
   ZpMean = new double[nSamples_];
   ZpQ = new double[nSamples_];
   ZpS2 = new double[nSamples_];
   ZpQ1 = new double[nSamples_];
   ZpMedian = new double[nSamples_];
   ZpQ2 = new double[nSamples_];
   ZZMean = new double[totPts];
   ZZQ = new double[totPts];
   ZZS2 = new double[totPts];
   ZZQ1 = new double[totPts];
   ZZMedian = new double[totPts];
   ZZQ2 = new double[totPts];
   tgp_->GetStats(iZero,ZpMean,ZZMean,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                  ZpS2,ZZS2,NULL,NULL,NULL,ZpQ1,ZpMedian,
                  ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
   delete [] essOut;
   delete [] ZpMean;
   delete [] ZpQ;
   delete [] ZpS2;
   delete [] ZpQ1;
   delete [] ZpMedian;
   delete [] ZpQ2;
   delete [] ZZQ;
   delete [] ZZS2;
   delete [] ZZQ1;
   delete [] ZZMedian;
   delete [] ZZQ2;
#else
   ZZMean = new double[totPts];
   evaluatePoint(totPts, XX, ZZMean);
#endif

   tgp_->GetPseudoPrior(itemps);
   delete [] itemps;
   delete [] dparams;

   gpcs = new double[4];
   tgp_->GetTreeStats(gpcs);
   delete [] gpcs;

   deleteRNGstate(tgp_state);
   tgp_state = NULL;

   if (outputLevel_ >= 1) printf("TGP interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = ZZMean;
   delete [] XX;
#else
   printf("PSUADE ERROR : TGP not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int TGP::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                      double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_TGP
   int    BTE[3], R=1, linBurn=1, PS_FALSE=0, PS_TRUE=1, bte0, bte1;
   int    totPts, ii, jj, ll, kk, count, *stateIn, index;
   double *HX, *XX, *dparams, *itemps, *gpcs, *ZZMean;
   void   *tgp_state=NULL;
   char   lineOut[1000];

   bte0 = 2000;
   bte1 = 7000;
   if (psRSExpertMode_ == 1)
   {
      sprintf(lineOut, "TGP: enter burn-in sample size (e.g. 500-5000): ");
      bte0 = getInt(500, 5000, lineOut); 
      sprintf(lineOut,"TGP: enter MCMC sample size (e.g. %d-20000): ",2000+bte0);
      bte1 = getInt(2000+bte0, 20000, lineOut); 
   }

   if (outputLevel_ >= 1) printf("TGP training begins....\n");
   if (tgp_ != NULL) delete tgp_;

   stateIn = new int[3];
   stateIn[0] = 782;
   stateIn[1] = 267;
   stateIn[2] = 218;
   unsigned int lstate = three2lstate(stateIn);
   tgp_state = newRNGstate(lstate);
   delete [] stateIn;

   BTE[0] = bte0;
   BTE[1] = bte1;
   BTE[2] = 2;
   dparams = new double[(nInputs_+1)*(nInputs_+1)+nInputs_+45];
   checkAllocate(dparams, "dparams in TGP::gen3DGridData");
   dparams[0] = 0.5;  // tree prior alpha
   dparams[1] = 2.0;  // tree prior beta
   dparams[2] = 10.0; // tree prior minpart
   if ((nInputs_ + 2) > 10) dparams[2] = 1.0 * (nInputs_ + 2.0);
   dparams[3] = 1.0;  // tree prior splitmin
   dparams[4] = 1.0 * nInputs_;  // tree prior basemax
   dparams[5] = 0.0;  // meanfn:  0 - linear, 1 : constant
   dparams[6] = 2.0;  // bflat in {b0, bmle, bflat, b0not,..}
   for (ii = 0; ii <= nInputs_; ii++) dparams[7+ii] = 0.0;
   count = nInputs_ + 8;
   dparams[count] = 1.0;  // Wi (1, (nInputs+1)^2 zeros)
   for (ii = 0; ii < (nInputs_+1)*(nInputs_+1)-1; ii++)
      dparams[count+1+ii] = 0.0;
   for (ii = 1; ii < nInputs_; ii++) dparams[count+ii*(nInputs_+1)] = 1.0;
   count += (nInputs_ + 1) * (nInputs_ + 1);
   dparams[count++] = 1.0;   // s2tau2 = c(1,1)
   dparams[count++] = 1.0;
   dparams[count++] = 5.0;   // s2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // s2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 5.0;   // tau2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // tau2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 0.1;
   dparams[count++] = 1;     // "expsep"
   dparams[count++] = 0.1;   // gd (1)
   dparams[count++] = 0.5;   // gd (2)
   dparams[count++] = 1;     // nug.p (1)
   dparams[count++] = 1;     // nug.p (2)
   dparams[count++] = 1;     // nug.p (3)
   dparams[count++] = 1;     // nug.p (4)
   dparams[count++] = -1;    // nug.lam (1)
   dparams[count++] = -1;    // nug.lam (2)
   dparams[count++] = -1;    // nug.lam (3)
   dparams[count++] = -1;    // nug.lam (4)
   dparams[count++] = 10.0;  // gamma(1)
   dparams[count++] = 0.2;   // gamma(2)
   dparams[count++] = 0.7;   // gamma(3)
   dparams[count++] = 1.0;   // d.p(1)
   dparams[count++] = 20.0;  // d.p(2)
   dparams[count++] = 10.0;  // d.p(3)
   dparams[count++] = 10.0;  // d.p(4)
   dparams[count++] = -1;    // d.lam (1)
   dparams[count++] = -1;    // d.lam (2)
   dparams[count++] = -1;    // d.lam (3)
   dparams[count++] = -1;    // d.lam (4)
   dparams[count++] = 0;     // nu?
   itemps = new double[7];
   itemps[0] = 1;
   itemps[1] = 0;
   itemps[2] = 0;
   itemps[3] = 1;
   itemps[4] = 1;
   itemps[5] = 0;
   itemps[6] = 1;


   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[3*totPts];
   checkAllocate(*X2, "X2 in TGP::gen3DGridData");
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         for (ll = 0; ll < nPtsPerDim_; ll++) 
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            for (kk = 0; kk < nInputs_; kk++) 
               XX[index*nInputs_+kk] = settings[kk]; 
            XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
            XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
            XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
            (*X2)[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            (*X2)[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            (*X2)[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }
   delete [] HX;
    
   tgp_ = new Tgp((void*) tgp_state, (int) nSamples_, (int) nInputs_, 
                  (int) nSamples_, (int) BTE[0], (int) BTE[1], (int) BTE[2], 
                  (int) R, (int) linBurn, (bool) PS_TRUE, (bool) PS_TRUE, 
                  (bool) PS_FALSE, (int) PS_FALSE, (bool) PS_FALSE, 
                  X, Y, X, dparams, itemps, (bool) PS_FALSE, 
                  (int) PS_FALSE, (double *) NULL, (double *) NULL);

   tgp_->Init();

   tgp_->Rounds();

   if (outputLevel_ >= 1) printf("TGP interpolation begins....\n");
#if 0
   tgp_->LoadInterpolatedPts(totPts, XX);
   tgp_->Predict();
   essOut = new double[3];
   ZpMean = new double[nSamples_];
   ZpQ = new double[nSamples_];
   ZpS2 = new double[nSamples_];
   ZpQ1 = new double[nSamples_];
   ZpMedian = new double[nSamples_];
   ZpQ2 = new double[nSamples_];
   ZZMean = new double[totPts];
   ZZQ = new double[totPts];
   ZZS2 = new double[totPts];
   ZZQ1 = new double[totPts];
   ZZMedian = new double[totPts];
   ZZQ2 = new double[totPts];
   tgp_->GetStats(iZero,ZpMean,ZZMean,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                  ZpS2,ZZS2,NULL,NULL,NULL,ZpQ1,ZpMedian,
                  ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
   delete [] essOut;
   delete [] ZpMean;
   delete [] ZpQ;
   delete [] ZpS2;
   delete [] ZpQ1;
   delete [] ZpMedian;
   delete [] ZpQ2;
   delete [] ZZQ;
   delete [] ZZS2;
   delete [] ZZQ1;
   delete [] ZZMedian;
   delete [] ZZQ2;
#else
   ZZMean = new double[totPts];
   evaluatePoint(totPts, XX, ZZMean);
#endif

   tgp_->GetPseudoPrior(itemps);

   gpcs = new double[4];
   tgp_->GetTreeStats(gpcs);
   delete [] gpcs;

   deleteRNGstate(tgp_state);
   tgp_state = NULL;
   delete [] itemps;
   delete [] dparams;

   if (outputLevel_ >= 1) printf("TGP interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = ZZMean;
   delete [] XX;
#else
   printf("PSUADE ERROR : TGP not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int TGP::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                       int ind4, double *settings, int *n, double **X2, 
                       double **Y2)
{
#ifdef HAVE_TGP
   int    BTE[3], R=1, linBurn=1, PS_FALSE=0, PS_TRUE=1, bte0, bte1;
   int    totPts, ii, jj, ll, mm, kk, count, *stateIn, index;
   double *HX, *XX, *dparams, *itemps, *gpcs, *ZZMean;
   void   *tgp_state=NULL;
   char   lineOut[1000];

   bte0 = 2000;
   bte1 = 7000;
   if (psRSExpertMode_ == 1)
   {
      sprintf(lineOut, "TGP: enter burn-in sample size (e.g. 500-5000): ");
      bte0 = getInt(500, 5000, lineOut); 
      sprintf(lineOut,"TGP: enter MCMC sample size (e.g. %d-20000): ",2000+bte0);
      bte1 = getInt(2000+bte0, 20000, lineOut); 
   }

   if (outputLevel_ >= 1) printf("TGP training begins....\n");
   if (tgp_ != NULL) delete tgp_;

   stateIn = new int[3];
   stateIn[0] = 782;
   stateIn[1] = 267;
   stateIn[2] = 218;
   unsigned int lstate = three2lstate(stateIn);
   tgp_state = newRNGstate(lstate);
   delete [] stateIn;

   BTE[0] = bte0;
   BTE[1] = bte1;
   BTE[2] = 2;
   dparams = new double[(nInputs_+1)*(nInputs_+1)+nInputs_+45];
   checkAllocate(dparams, "dparams in TGP::gen4DGridData");
   dparams[0] = 0.5;  // tree prior alpha
   dparams[1] = 2.0;  // tree prior beta
   dparams[2] = 10.0; // tree prior minpart
   if ((nInputs_ + 2) > 10) dparams[2] = 1.0 * (nInputs_ + 2.0);
   dparams[3] = 1.0;  // tree prior splitmin
   dparams[4] = 1.0 * nInputs_;  // tree prior basemax
   dparams[5] = 0.0;  // meanfn:  0 - linear, 1 : constant
   dparams[6] = 2.0;  // bflat in {b0, bmle, bflat, b0not,..}
   for (ii = 0; ii <= nInputs_; ii++) dparams[7+ii] = 0.0;
   count = nInputs_ + 8;
   dparams[count] = 1.0;  // Wi (1, (nInputs+1)^2 zeros)
   for (ii = 0; ii < (nInputs_+1)*(nInputs_+1)-1; ii++)
      dparams[count+1+ii] = 0.0;
   for (ii = 1; ii < nInputs_; ii++) dparams[count+ii*(nInputs_+1)] = 1.0;
   count += (nInputs_ + 1) * (nInputs_ + 1);
   dparams[count++] = 1.0;   // s2tau2 = c(1,1)
   dparams[count++] = 1.0;
   dparams[count++] = 5.0;   // s2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // s2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 5.0;   // tau2 prior (a0, g0) = c(5,10)
   dparams[count++] = 10.0;
   dparams[count++] = 0.2;   // tau2 hierar inv-gamma prior = c(5,10)
   dparams[count++] = 0.1;
   dparams[count++] = 1;     // "expsep"
   dparams[count++] = 0.1;   // gd (1)
   dparams[count++] = 0.5;   // gd (2)
   dparams[count++] = 1;     // nug.p (1)
   dparams[count++] = 1;     // nug.p (2)
   dparams[count++] = 1;     // nug.p (3)
   dparams[count++] = 1;     // nug.p (4)
   dparams[count++] = -1;    // nug.lam (1)
   dparams[count++] = -1;    // nug.lam (2)
   dparams[count++] = -1;    // nug.lam (3)
   dparams[count++] = -1;    // nug.lam (4)
   dparams[count++] = 10.0;  // gamma(1)
   dparams[count++] = 0.2;   // gamma(2)
   dparams[count++] = 0.7;   // gamma(3)
   dparams[count++] = 1.0;   // d.p(1)
   dparams[count++] = 20.0;  // d.p(2)
   dparams[count++] = 10.0;  // d.p(3)
   dparams[count++] = 10.0;  // d.p(4)
   dparams[count++] = -1;    // d.lam (1)
   dparams[count++] = -1;    // d.lam (2)
   dparams[count++] = -1;    // d.lam (3)
   dparams[count++] = -1;    // d.lam (4)
   dparams[count++] = 0;     // nu?
   itemps = new double[7];
   itemps[0] = 1;
   itemps[1] = 0;
   itemps[2] = 0;
   itemps[3] = 1;
   itemps[4] = 1;
   itemps[5] = 0;
   itemps[6] = 1;


   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[4*totPts];
   checkAllocate(*X2, "X2 in TGP::gen4DGridData");
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++) 
      {
         for (ll = 0; ll < nPtsPerDim_; ll++) 
         {
            for (mm = 0; mm < nPtsPerDim_; mm++) 
            {
               index = ii*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ +
                       jj* nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               for (kk = 0; kk < nInputs_; kk++) 
                  XX[index*nInputs_+kk] = settings[kk]; 
               XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
               XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
               XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
               XX[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
               (*X2)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               (*X2)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               (*X2)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               (*X2)[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }
   delete [] HX;
    
   tgp_ = new Tgp((void*) tgp_state, (int) nSamples_, (int) nInputs_, 
                  (int) nSamples_, (int) BTE[0], (int) BTE[1], (int) BTE[2], 
                  (int) R, (int) linBurn, (bool) PS_TRUE, (bool) PS_TRUE, 
                  (bool) PS_FALSE, (int) PS_FALSE, (bool) PS_FALSE, 
                  X, Y, X, dparams, itemps, (bool) PS_FALSE, 
                  (int) PS_FALSE, (double *) NULL, (double *) NULL);

   tgp_->Init();

   tgp_->Rounds();

   if (outputLevel_ >= 1) printf("TGP interpolation begins....\n");
#if 0
   tgp_->LoadInterpolatedPts(totPts, XX);
   tgp_->Predict();
   essOut = new double[3];
   ZpMean = new double[nSamples_];
   ZpQ = new double[nSamples_];
   ZpS2 = new double[nSamples_];
   ZpQ1 = new double[nSamples_];
   ZpMedian = new double[nSamples_];
   ZpQ2 = new double[nSamples_];
   ZZMean = new double[totPts];
   ZZQ = new double[totPts];
   ZZS2 = new double[totPts];
   ZZQ1 = new double[totPts];
   ZZMedian = new double[totPts];
   ZZQ2 = new double[totPts];
   tgp_->GetStats(iZero,ZpMean,ZZMean,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                  ZpS2,ZZS2,NULL,NULL,NULL,ZpQ1,ZpMedian,
                  ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
   delete [] essOut;
   delete [] ZpMean;
   delete [] ZpQ;
   delete [] ZpS2;
   delete [] ZpQ1;
   delete [] ZpMedian;
   delete [] ZpQ2;
   delete [] ZZQ;
   delete [] ZZS2;
   delete [] ZZQ1;
   delete [] ZZMedian;
   delete [] ZZQ2;
#else
   ZZMean = new double[totPts];
   evaluatePoint(totPts, XX, ZZMean);
#endif

   tgp_->GetPseudoPrior(itemps);
   delete [] itemps;
   delete [] dparams;

   gpcs = new double[4];
   tgp_->GetTreeStats(gpcs);
   delete [] gpcs;

   deleteRNGstate(tgp_state);
   tgp_state = NULL;

   if (outputLevel_ >= 1) printf("TGP interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = ZZMean;
   delete [] XX;
#else
   printf("PSUADE ERROR : TGP not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double TGP::evaluatePoint(double *X)
{
#ifdef HAVE_TGP
#if 1
   int    iOne=1;
   double Y;
   evaluatePoint(iOne, X, &Y);
   return Y;
#else
   int    iOne=1, iZero=0, PS_FALSE=0;
   double *essOut, ZpMean, ZpQ, ZpS2;
   double ZpQ1, ZpMedian, ZpQ2, ZZMean, ZZQ, ZZS2;
   double ZZQ1, ZZMedian, ZZQ2;
   tgp_->LoadInterpolatedPts(iOne, X);
   tgp_->Predict();
   essOut = new double[3];
   tgp_->GetStats(iZero,&ZpMean,&ZZMean,NULL,NULL,&ZpQ,&ZZQ,PS_FALSE,
                  &ZpS2,&ZZS2,NULL,NULL,NULL,&ZpQ1,&ZpMedian,
                  &ZpQ2,&ZZQ1,&ZZMedian,&ZZQ2,NULL,NULL,NULL,essOut);
   delete [] essOut;
   return ZZMean;
#endif
#else
   printf("PSUADE ERROR : TGP not installed.\n");
   return 0.0;
#endif
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double TGP::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_TGP
   int    iZero=0, PS_FALSE=0, chunksize;
   double *essOut, *ZpMean, *ZpQ, *ZpS2;
   double *ZpQ1, *ZpMedian, *ZpQ2, *ZZQ, *ZZS2;
   double *ZZQ1, *ZZMedian, *ZZQ2;
   essOut = new double[3];
   ZpMean = new double[nSamples_];
   ZpQ    = new double[nSamples_];
   ZpS2   = new double[nSamples_];
   ZpQ1   = new double[nSamples_];
   ZpMedian = new double[nSamples_];
   ZpQ2   = new double[nSamples_];
   ZZQ    = new double[nSamples_];
   ZZS2   = new double[nSamples_];
   ZZQ1   = new double[nSamples_];
   ZZMedian = new double[nSamples_];
   ZZQ2 = new double[nSamples_];

   double *XX, *YY;
   int    ii, jj, kk, nChunks = npts / nSamples_;
   if (nChunks == 0) nChunks = 1;
   XX = new double[nSamples_*nInputs_];
   YY = new double[nSamples_];
   checkAllocate(YY, "YY in TGP::evaluatePoint");
   
   for (ii = 0; ii < nChunks; ii++)
   {
      if (outputLevel_ > 1)
         printf("TGP: evaluate chunk %d (out of %d)\n",ii+1,nChunks);
      chunksize = nSamples_;
      if (ii*nSamples_+chunksize > npts) chunksize = npts - ii * nSamples_;
      for (jj = 0; jj < chunksize*nInputs_; jj++)
         XX[jj] = X[ii*nSamples_*nInputs_+jj]; 
      for (jj = chunksize; jj < nSamples_; jj++)
         for (kk = 0; kk < nInputs_; kk++)
            XX[jj*nInputs_+kk] = X[ii*nSamples_*nInputs_+(chunksize-1)*nInputs_+kk]; 
      tgp_->LoadInterpolatedPts(nSamples_, XX);
      tgp_->Predict();
      tgp_->GetStats(iZero,ZpMean,YY,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                     ZpS2,ZZS2,NULL,NULL,NULL,ZpQ1,ZpMedian,
                     ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
      for (jj = 0; jj < chunksize; jj++) Y[ii*nSamples_+jj] = YY[jj]; 
   }
   if (nSamples_*nChunks < npts)
   {
      for (ii = nSamples_*nChunks*nInputs_; ii < npts*nInputs_; ii++)
         XX[ii-nSamples_*nChunks*nInputs_] = X[ii];
      for (ii = npts; ii < (nChunks+1)*nSamples_; ii++)
         for (jj = 0; jj < nInputs_; jj++) 
            XX[ii*nInputs_+jj] = X[(npts-1)*nInputs_+jj];
      tgp_->LoadInterpolatedPts(nSamples_, XX);
      tgp_->Predict();
      tgp_->GetStats(iZero,ZpMean,YY,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                     ZpS2,ZZS2,NULL,NULL,NULL,ZpQ1,ZpMedian,
                     ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
      for (jj = 0; jj < npts-nChunks*nSamples_; jj++)
         Y[nChunks*nSamples_+jj] = YY[jj]; 
   }
   delete [] XX;
   delete [] YY;
   delete [] essOut;
   delete [] ZpMean;
   delete [] ZpQ;
   delete [] ZpS2;
   delete [] ZpQ1;
   delete [] ZpMedian;
   delete [] ZpQ2;
   delete [] ZZQ;
   delete [] ZZS2;
   delete [] ZZQ1;
   delete [] ZZMedian;
   delete [] ZZQ2;
#else
   printf("PSUADE ERROR : TGP not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// Evaluate a given point and return also the standard deviation 
// ------------------------------------------------------------------------
double TGP::evaluatePointFuzzy(double *X, double &std)
{
#ifdef HAVE_TGP
#if 1
   int    iOne=1;
   double Y;
   evaluatePointFuzzy(iOne, X, &Y, &std);
   return Y;
#else
   int    iOne=1, iZero=0, PS_FALSE=0;
   double *essOut, ZpMean, ZpQ, ZpS2;
   double ZpQ1, ZpMedian, ZpQ2, ZZMean, ZZQ, ZZS2;
   double ZZQ1, ZZMedian, ZZQ2;
   tgp_->LoadInterpolatedPts(iOne, X);
   tgp_->Predict();
   essOut = new double[3];
   tgp_->GetStats(iZero,&ZpMean,&ZZMean,NULL,NULL,&ZpQ,&ZZQ,PS_FALSE,
                  &ZpS2,&ZZS2,NULL,NULL,NULL,&ZpQ1,&ZpMedian,
                  &ZpQ2,&ZZQ1,&ZZMedian,&ZZQ2,NULL,NULL,NULL,essOut);
   delete [] essOut;
   std = ZZS2;
   return ZZMean;
#endif
#else
   printf("PSUADE ERROR : TGP not installed.\n");
   return 0.0;
#endif
}

// ************************************************************************
// Evaluate a number of points and return also the standard deviations 
// ------------------------------------------------------------------------
double TGP::evaluatePointFuzzy(int npts,double *X, double *Y, double *Ystd)
{
#ifdef HAVE_TGP
   int    iZero=0, PS_FALSE=0;
   double *essOut, *ZpMean, *ZpQ, *ZpS2;
   double *ZpQ1, *ZpMedian, *ZpQ2, *ZZQ;
   double *ZZQ1, *ZZMedian, *ZZQ2;
   essOut = new double[3];
   ZpMean = new double[nSamples_];
   ZpQ    = new double[nSamples_];
   ZpS2   = new double[nSamples_];
   ZpQ1   = new double[nSamples_];
   ZpMedian = new double[nSamples_];
   ZpQ2   = new double[nSamples_];
   ZZQ    = new double[nSamples_];
   ZZQ1   = new double[nSamples_];
   ZZMedian = new double[nSamples_];
   ZZQ2 = new double[nSamples_];

   double *XX, *YY, *YS;
   int    ii, jj, nChunks = npts / nSamples_;
   if (nChunks == 0) nChunks = 1;
   XX = new double[nSamples_*nInputs_];
   YY = new double[nSamples_];
   YS = new double[nSamples_];
   checkAllocate(YS, "YS in TGP::evaluatePointFuzzy");
   
   for (ii = 0; ii < nChunks; ii++)
   {
      for (jj = 0; jj < nSamples_*nInputs_; jj++)
         XX[jj] = X[ii*nSamples_*nInputs_+jj]; 
      tgp_->LoadInterpolatedPts(nSamples_, XX);
      tgp_->Predict();
      tgp_->GetStats(iZero,ZpMean,YY,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                     ZpS2,YS,NULL,NULL,NULL,ZpQ1,ZpMedian,
                     ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
      for (jj = 0; jj < nSamples_; jj++) Y[ii*nSamples_+jj] = YY[jj]; 
      for (jj = 0; jj < nSamples_; jj++) Ystd[ii*nSamples_+jj] = YS[jj]; 
   }
   if (nSamples_*nChunks < npts)
   {
      for (ii = nSamples_*nChunks*nInputs_; ii < npts*nInputs_; ii++)
         XX[ii-nSamples_*nChunks*nInputs_] = X[ii];
      for (ii = npts; ii < (nChunks+1)*nSamples_; ii++)
         for (jj = 0; jj < nInputs_; jj++) 
            XX[ii*nInputs_-nSamples_*nChunks*nInputs_+jj] = 
                                 X[(npts-1)*nInputs_+jj];
      tgp_->LoadInterpolatedPts(nSamples_, XX);
      tgp_->Predict();
      tgp_->GetStats(iZero,ZpMean,YY,NULL,NULL,ZpQ,ZZQ,PS_FALSE,
                     ZpS2,YS,NULL,NULL,NULL,ZpQ1,ZpMedian,
                     ZpQ2,ZZQ1,ZZMedian,ZZQ2,NULL,NULL,NULL,essOut);
      for (jj = 0; jj < npts-nChunks*nSamples_; jj++)
         Y[nChunks*nSamples_+jj] = YY[jj]; 
      for (jj = 0; jj < npts-nChunks*nSamples_; jj++)
         Ystd[nChunks*nSamples_+jj] = YS[jj]; 
   }
   delete [] XX;
   delete [] YY;
   delete [] YS;
   delete [] essOut;
   delete [] ZpMean;
   delete [] ZpQ;
   delete [] ZpS2;
   delete [] ZpQ1;
   delete [] ZpMedian;
   delete [] ZpQ2;
   delete [] ZZQ;
   delete [] ZZQ1;
   delete [] ZZMedian;
   delete [] ZZQ2;
#else
   printf("PSUADE ERROR : TGP not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double TGP::setParams(int targc, char **targv)
{
   if (targc > 0 && !strcmp(targv[0], "improv"))
   {
#ifdef HAVE_TGP
#else
      printf("PSUADE ERROR : TGP not installed.\n");
#endif
   }
   return 0.0;
}

