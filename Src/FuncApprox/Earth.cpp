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
// Functions for the class Earth
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Earth.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"

#ifdef HAVE_EARTH
#define STANDALONE 1
#include "earth.h"
#endif

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class Earth
// ------------------------------------------------------------------------
Earth::Earth(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
#ifdef HAVE_EARTH
   char pString[500];

   faID_ = PSUADE_RS_EARTH;

   if (nSamples_ > 200) maxTerms_ = 100;
   else                 maxTerms_ = nSamples_/2 + 1;

   if (nInputs >= 8) maxDegree_ = 4;
   else              maxDegree_ = nInputs;

   thresh_ = 0.00001;
   penalty_ = (maxDegree_ > 1) ? 3 : 2;
   penalty_ = 1.0;
   nMinSpan_ = 1;
   prune_ = true;
   nFastK_ = 20; 
   fastBeta_ = 0;
   newVarPenalty_ = 1.0;
   useBetaCache_ = true;
   nTrace_ =  0;
   numTerms_ = 0;

   if (psRSExpertMode_ == 1)
   {
      sprintf(pString,"Enter the number of basis functions (>10, <= %d): ",
              nSamples);
      maxTerms_ = getInt(10, nSamples, pString);
      sprintf(pString, "Enter the degree of interactions (<=%d) : ", nInputs);
      maxDegree_ = getInt(1, nInputs, pString);
   }

   bestSet_ = new bool[maxTerms_];
   dirs_    = new int[maxTerms_ * nInputs];
   cuts_    = new double[maxTerms_ * nInputs];
   betas_   = new double[maxTerms_];
   checkAllocate(betas_, "betas in Earth::constructor");
   wgts_    = NULL;
  
#else
   printf("PSUADE ERROR : Earth not installed.\n");
   exit(1);
#endif
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Earth::~Earth()
{
   if (bestSet_ != NULL) delete [] bestSet_;
   if (dirs_    != NULL) delete [] dirs_;
   if (cuts_    != NULL) delete [] cuts_;
   if (betas_   != NULL) delete [] betas_;
   if (wgts_    != NULL) delete [] wgts_;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Earth::initialize(double *X, double *Y)
{
#ifdef HAVE_EARTH
   int    *LinPreds, ii, jj, iOne=1;
   double *Residuals, *BX, BestGcv, *XX;

   Residuals = new double[nSamples_];
   BX = new double[nSamples_ * maxTerms_];
   LinPreds = new int[nInputs_];
   checkAllocate(LinPreds, "LinPreds in Earth::initialize");
   for (jj = 0; jj < nInputs_; jj++) LinPreds[jj] = 0;

   XX = new double[nSamples_ * nInputs_];
   for (ii = 0; ii < nSamples_; ii++)
      for (jj = 0; jj < nInputs_; jj++)
         XX[jj*nSamples_+ii] = X[ii*nInputs_+jj]; 

   if (outputLevel_ >= 2) printf("Entering Earth processing\n");
   if (outputLevel_ >= 1) 
      printf("If it crashes here, it is earth's problem (try larger N).\n");
   if (outputLevel_ >= 5) nTrace_ = 4;
   
   TrainEarth(&BestGcv, &numTerms_, bestSet_, BX, dirs_, cuts_, Residuals, 
         betas_, XX, Y, NULL, nSamples_, iOne, nInputs_, maxDegree_, maxTerms_,
         penalty_, thresh_, nMinSpan_, prune_, nFastK_, fastBeta_, 
         newVarPenalty_, LinPreds, useBetaCache_, nTrace_, NULL);
   delete [] XX;
   delete [] Residuals;
   delete [] BX;
   delete [] LinPreds;
   if (outputLevel_ >= 5) 
      FormatEarth(bestSet_,dirs_,cuts_,betas_,nInputs_,iOne,numTerms_,
                  maxTerms_,3,0);
   if (outputLevel_ >= 2) 
      printf("Returning from Earth, BestGcv = %e\n",BestGcv);
   if (psRSCodeGen_ == 1) 
      printf("Earth INFO: response surface stand-alone code not available.\n");
   return 0;
#else
   printf("PSUADE ERROR : Earth not used.\n");
   return -1;
#endif
}
  
// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Earth::genNDGridData(double *X, double *Y, int *N, double **X2, 
                        double **Y2)
{
#ifdef HAVE_EARTH
   int ss, iOne=1, totPts;
   initialize(X,Y);
   if ((*N) == -999) return 0;
  
   genNDGrid(N, X2);
   if ((*N) == 0) return 0;
   totPts = (*N);

   (*Y2) = new double[totPts];
   checkAllocate(*Y2, "Y2 in Earth::genNDGridData");
   (*N) = totPts;

   if (outputLevel_ >= 2) printf("Entering Earth prediction\n");
   if (outputLevel_ >= 1) 
      printf("If it crashes here, it is Earth prediction problem.\n");
   for (ss = 0; ss < totPts; ss++)
       PredictEarth(&((*Y2)[ss]), &((*X2)[ss*nInputs_]), bestSet_, 
             dirs_, cuts_, betas_, nInputs_, iOne, numTerms_, maxTerms_);

   if (outputLevel_ >= 2) printf("Returning from Earth\n");
#else
   printf("PSUADE ERROR : Earth not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Earth::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
                         int *N, double **X2, double **Y2)
{
#ifdef HAVE_EARTH
   int    ss, ii, iOne=1, totPts;
   double *XX, HX;
   initialize(X,Y);
   if ((*N) == -999) return 0;
  
   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts];
   (*Y2) = new double[totPts];
   XX = new double[totPts*nInputs_];
   checkAllocate(XX, "XX in Earth::gen1DGridData");
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XX[ss*nInputs_+ii] = settings[ii]; 
    
   for (ss = 0; ss < totPts; ss++) 
   {
      XX[ss*nInputs_+ind1] = HX * ss + lowerBounds_[ind1];
      (*X2)[ss] = HX * ss + lowerBounds_[ind1];
      (*Y2)[ii] = 0.0;
   }

   if (outputLevel_ >= 2) printf("Entering Earth prediction\n");
   if (outputLevel_ >= 1) 
      printf("If it crashes here, it is Earth prediction problem.\n");
   for (ss = 0; ss < totPts; ss++)
       PredictEarth(&((*Y2)[ss]), &(XX[ss*nInputs_]), bestSet_, 
             dirs_, cuts_, betas_, nInputs_, iOne, numTerms_, maxTerms_);
   if (outputLevel_ >= 2) printf("Returning from Earth\n");
   (*N) = totPts;
   delete [] XX;
#else
   printf("PSUADE ERROR : Earth not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Earth::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                         double *settings, int *N, double **X2, double **Y2)
{
#ifdef HAVE_EARTH
   int    ss, ii, jj, index, iOne=1, totPts;
   double *XX, *HX;

   initialize(X,Y);
   if ((*N) == -999) return 0;
  
   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts*2];
   (*Y2) = new double[totPts];
   XX = new double[totPts*nInputs_];
   checkAllocate(XX, "XX in Earth::gen2DGridData");
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XX[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
         XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
         (*X2)[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         (*X2)[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   if (outputLevel_ >= 2) printf("Entering Earth prediction\n");
   if (outputLevel_ >= 1) 
      printf("If it crashes here, it is Earth prediction problem.\n");
   for (ss = 0; ss < totPts; ss++)
       PredictEarth(&((*Y2)[ss]), &(XX[ss*nInputs_]), bestSet_, 
             dirs_, cuts_, betas_, nInputs_, iOne, numTerms_, maxTerms_);
   if (outputLevel_ >= 2) printf("Returning from Earth\n");
   (*N) = totPts;
   delete [] HX;
   delete [] XX;
#else
   printf("PSUADE ERROR : Earth not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Earth::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                         double *settings, int *N, double **X2, double **Y2)
{
#ifdef HAVE_EARTH
   int    ss, ii, jj, ll, index, iOne=1, totPts;
   double *XX, *HX;

   initialize(X,Y);
   if ((*N) == -999) return 0;
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts*3];
   (*Y2) = new double[totPts];
   XX = new double[totPts*nInputs_];
   checkAllocate(XX, "XX in Earth::gen3DGridData");
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XX[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            XX[index*nInputs_+ind1] = HX[0] * ii + lowerBounds_[ind1];
            XX[index*nInputs_+ind2] = HX[1] * jj + lowerBounds_[ind2];
            XX[index*nInputs_+ind3] = HX[2] * ll + lowerBounds_[ind3];
            (*X2)[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            (*X2)[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            (*X2)[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   if (outputLevel_ >= 2) printf("Entering Earth prediction\n");
   if (outputLevel_ >= 1) 
      printf("If it crashes here, it is Earth prediction problem.\n");
   for (ss = 0; ss < totPts; ss++)
       PredictEarth(&((*Y2)[ss]), &((*X2)[ss*nInputs_]), bestSet_, 
             dirs_, cuts_, betas_, nInputs_, iOne, numTerms_, maxTerms_);
   if (outputLevel_ >= 2) printf("Returning from Earth\n");
   (*N) = totPts;
   delete [] HX;
   delete [] XX;
#else
   printf("PSUADE ERROR : Earth not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Earth::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                        int ind4, double *settings, int *N, double **X2, 
                        double **Y2)
{
#ifdef HAVE_EARTH
   int    ss, ii, jj, ll, mm, index, iOne=1, totPts;
   double *XX, *HX;

   initialize(X,Y);
   if ((*N) == -999) return 0;
  
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts*4];
   (*Y2) = new double[totPts];
   XX = new double[totPts*nInputs_];
   checkAllocate(XX, "XX in Earth::gen4DGridData");
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XX[ss*nInputs_+ii] = settings[ii]; 

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            for (mm = 0; mm < nPtsPerDim_; mm++)
            {
               index = ii*nPtsPerDim_*nPtsPerDim_ * nPtsPerDim_ +
                       jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               XX[index*nInputs_+ind1] = HX[0] * ii + lowerBounds_[ind1];
               XX[index*nInputs_+ind2] = HX[1] * jj + lowerBounds_[ind2];
               XX[index*nInputs_+ind3] = HX[2] * ll + lowerBounds_[ind3];
               XX[index*nInputs_+ind4] = HX[3] * mm + lowerBounds_[ind4];
               (*X2)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               (*X2)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               (*X2)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               (*X2)[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }

   if (outputLevel_ >= 2) printf("Entering Earth prediction\n");
   if (outputLevel_ >= 1) 
      printf("If it crashes here, it is Earth prediction problem.\n");
   for (ss = 0; ss < totPts; ss++)
       PredictEarth(&((*Y2)[ss]), &((*X2)[ss*nInputs_]), bestSet_, 
             dirs_, cuts_, betas_, nInputs_, iOne, numTerms_, maxTerms_);
   if (outputLevel_ >= 2) printf("Returning from Earth\n");
   (*N) = totPts;
   delete [] XX;
   delete [] HX;
#else
   printf("PSUADE ERROR : Earth not used.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Earth::evaluatePoint(double *X)
{
   double Y=0.0;
#ifdef HAVE_EARTH
   int    iOne=1;
   PredictEarth(&Y, X, bestSet_, dirs_, cuts_, betas_, nInputs_, iOne, 
                numTerms_, maxTerms_);
#else
   printf("PSUADE ERROR : Earth not used.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double Earth::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_EARTH
   int ii, iOne=1;
   for (ii = 0; ii < npts; ii++)
      PredictEarth(&(Y[ii]), &(X[ii*nInputs_]), bestSet_, dirs_, cuts_, 
                   betas_, nInputs_, iOne, numTerms_, maxTerms_);
#else
   printf("PSUADE ERROR : Earth not used.\n");
#endif
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double Earth::evaluatePointFuzzy(double *X, double &std)
{
   double Y=0.0;
#ifdef HAVE_EARTH
   int    iOne=1;
   PredictEarth(&Y, X, bestSet_, dirs_, cuts_, betas_, nInputs_, iOne, 
                numTerms_, maxTerms_);
#else
   printf("PSUADE ERROR : Earth not used.\n");
#endif
   std = 0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double Earth::evaluatePointFuzzy(int npts,double *X, double *Y,double *Ystd)
{
#ifdef HAVE_EARTH
   int ii, iOne=1;
   for (ii = 0; ii < npts; ii++)
   {
      PredictEarth(&(Y[ii]), &(X[ii*nInputs_]), bestSet_, dirs_, cuts_, 
                   betas_, nInputs_, iOne, numTerms_, maxTerms_);
      Ystd[ii] = 0.0;
   }
#else
   printf("PSUADE ERROR : Earth not used.\n");
#endif
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
Earth& Earth::operator=(const Earth &)
{
   printf("Earth operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

