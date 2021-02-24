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
// Functions for the class SVM
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "sysdef.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "SVM.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#ifdef HAVE_SVM
extern "C" 
{
  void SVMSetGamma(double gamma, double tol);
  void SVMSetKernel(int kernel);
  void SVMTrain(int nInputs, int nTrains, double *trainInputs,
                double *trainOutput, int, double *, double *);

  void SVMInterp(int nTests, int nInputs, double *inputs, 
                 double *output, double *stds);
}
#endif

// ************************************************************************
// Constructor for object class SVM
// ------------------------------------------------------------------------
SVM::SVM(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
#ifdef HAVE_SVM
   int    idata;
   double ddata;
   char   *inStr, winput1[500], winput2[500];

   faID_ = PSUADE_RS_SVM;
   gamma_ = 1.0;
   tolerance_ = 1.0e-6;
   kernel_ = 2;
   if (psRSExpertMode_ != 1 && psConfig_ != NULL)
   {
      inStr = psConfig_->getParameter("SVM_tol");
      if (inStr != NULL)
      {
         sscanf(inStr, "%s %s %lg\n", winput1, winput2, &ddata);
         if (winput2[0] != '=')
            printf("SVM read config file syntax error : %s\n", inStr);
         else if (ddata < 0.0 || ddata < 1.0-6 || ddata >= 1e6)
         {
            printf("SVM read config file : read tol error : %e\n", ddata);
            printf("                       tol kept at %e\n", tolerance_);
         }
         else
         {
            tolerance_ = ddata;
            printf("SVM : tol   set to %e (config)\n", tolerance_);
         }
      }
      inStr = psConfig_->getParameter("SVM_gamma");
      if (inStr != NULL)
      {
         sscanf(inStr, "%s %s %lg\n", winput1, winput2, &ddata);
         if (winput2[0] != '=')
            printf("SVM read config file syntax error : %s\n", inStr);
         else if (ddata < 0.0 || ddata < 1.0-6 || ddata > 1.0e6)
         {
            printf("SVM read config file : gamma error : %e\n", ddata);
            printf("                       gamma kept at %e\n", gamma_);
         }
         else
         {
            gamma_ = ddata;
            printf("SVM : gamma set to %e (config)\n", gamma_);
         }
      }
      inStr = psConfig_->getParameter("SVM_kernel");
      if (inStr != NULL)
      {
         sscanf(inStr, "%s %s %d\n", winput1, winput2, &idata);
         if (winput2[0] != '=')
            printf("SVM read config file syntax error : %s\n", inStr);
         else if (idata < 1 || idata > 4)
         {
            printf("SVM read config file : kernel error : %d (must be 1-4)\n",
                   idata);
            printf("                       kernel kept at %d\n",kernel_+1);
         }
         else
         {
            kernel_ = idata - 1;
            printf("SVM : kernel set to %d (config)\n", kernel_);
         }
      }
   }
   if (psRSExpertMode_ == 1)
   {
      printf("SVM kernel options: \n");
      printf("1. linear\n");
      printf("2. third order polynomial\n");
      printf("3. radial basis function\n");
      printf("4. sigmoid function\n");
      printf("SVM: enter kernel type (1 - 4) : ");
      scanf("%d", &kernel_);
      if (kernel_ < 1 || kernel_ > 4)
      {
         printf("SVM ERROR : invalid kernel type, set to 3.\n");
         kernel_ = 3;
      }
      kernel_--;
      printf("SVM: enter tolerance (1.0e-6 to 1e6) : ");
      scanf("%lg", &tolerance_);
      if (tolerance_ < 1.0-6 || tolerance_ > 1e6)
      {
         printf("SVM ERROR : invalid tolerance, set to 1.0e-4\n");
         tolerance_ = 1.0e-4;
      }
      if (kernel_ == 2)
      {
         printf("SVM: enter RBF_gamma (1.0e-6 to 1.0e6) : ");
         scanf("%lg", &gamma_);
         if (gamma_ < 1.0-6 || gamma_ > 1.0e6)
         {
            printf("SVM ERROR : invalid RBF_gamma, set to 1.0\n");
            gamma_ = 1.0;
         }
      }
      fgets(winput1, 500, stdin);
   }
   if (kernel_ != -1) SVMSetKernel(kernel_);
   if (gamma_ != -1.0 || tolerance_ != -1.0) SVMSetGamma(gamma_, tolerance_);
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
SVM::~SVM()
{
#ifndef HAVE_SVM
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int SVM::initialize(double *X, double *Y)
{
#ifdef HAVE_SVM
   int    ss;
   double *stds;

   stds = new double[nSamples_];
   checkAllocate(stds, "stds in SVM::initialize");
   if (outputLevel_ >= 1)
   {
      printf("SVM training begins....\n");
      if (kernel_ == 0) printf("SVM kernel = linear\n");
      if (kernel_ == 1) printf("SVM kernel = third order polynomial\n");
      if (kernel_ == 2) printf("SVM kernel = radial basis function\n");
      if (kernel_ == 3) printf("SVM kernel = sigmoid function\n");
      printf("SVM epsilon   = %e\n",tolerance_);
      if (kernel_ == 2)
      printf("SVM RBF_gamma = %e\n",gamma_);
   }
   if (gamma_ != -1.0 || tolerance_ != -1.0) SVMSetGamma(gamma_, tolerance_);
   if (kernel_ != -1) SVMSetKernel(kernel_);
   SVMTrain(nInputs_, nSamples_, X, Y, 0, NULL, stds);
   for (ss = 0; ss < nSamples_; ss++) stds[ss] = 0.0;
   if (outputLevel_ >= 1) printf("SVM training completed.\n");
   if (psRSCodeGen_ == 1) 
      printf("SVM INFO: response surface stand-alone code not available.\n");
#else
   printf("PSUADE ERROR : SVM not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int SVM::genNDGridData(double *X, double *Y, int *N2, double **X2, 
                      double **Y2)
{
#ifdef HAVE_SVM
   int totPts;

   // initialization
   initialize(X,Y);
   if ((*N2) == -999 || X2 == NULL || Y2 == NULL) return 0;
  
   // generating regular grid data
   genNDGrid(N2, X2);
   if ((*N2) == 0) return 0;
   totPts = (*N2);

   // generate the data points 
   (*Y2) = new double[totPts];
   checkAllocate(*Y2, "Y2 in SVM::genNDGridData");
   if (outputLevel_ >= 1) printf("SVM interpolation begins....\n");
   SVMInterp(totPts, nInputs_, *X2, *Y2, NULL);
   if (outputLevel_ >= 1) printf("SVM interpolation completed.\n");
#else
   printf("PSUADE ERROR : SVM not installed.\n");
   return -1;
#endif
   return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int SVM::gen1DGridData(double *X, double *Y, int ind1,
                      double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_SVM
   int    ii, kk, totPts;
   double HX, *XX, *YY;

   // initialization
   initialize(X,Y);
  
   // generating regular grid data
   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[2*totPts];
   XX = new double[totPts*nInputs_];
   checkAllocate(XX, "XX in SVM::gen1DGridData");
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (kk = 0; kk < nInputs_; kk++) 
         XX[ii*nInputs_+kk] = settings[kk]; 
      XX[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
      (*X2)[ii] = HX * ii + lowerBounds_[ind1];
   }
    
   YY = new double[totPts];
   checkAllocate(YY, "YY in SVM::gen1DGridData");
   if (outputLevel_ >= 1) printf("SVM interpolation begins....\n");
   SVMInterp(totPts, nInputs_, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("SVM interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int SVM::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                       double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_SVM
   int    ii, jj, kk, totPts, index;
   double *HX, *XX, *YY;

   // initialization
   initialize(X,Y);
  
   // generating regular grid data
   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[2*totPts];
   checkAllocate(*X2, "X2 in SVM::gen2DGridData");
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
    
   YY = new double[totPts];
   checkAllocate(YY, "YY in SVM::gen2DGridData");
   if (outputLevel_ >= 1) printf("SVM interpolation begins....\n");
   SVMInterp(totPts, nInputs_, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("SVM interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int SVM::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                       double *settings, int *n, double **X2, double **Y2)
{
#ifdef HAVE_SVM
   int    ii, jj, ll, kk, totPts, index;
   double *HX, *XX, *YY;

   // initialization
   initialize(X,Y);
  
   // set up for generating regular grid data
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[3*totPts];
   checkAllocate(*X2, "X2 in SVM::gen3DGridData");
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
    
   YY = new double[totPts];
   checkAllocate(YY, "YY in SVM::gen3DGridData");
   if (outputLevel_ >= 1) printf("SVM interpolation begins....\n");
   SVMInterp(totPts, nInputs_, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("SVM interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int SVM::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                       int ind4, double *settings, int *n, double **X2, 
                       double **Y2)
{
#ifdef HAVE_SVM
   int    ii, jj, ll, mm, kk, totPts, index;
   double *HX, *XX, *YY;

   // initialization
   initialize(X,Y);
  
   // generating regular grid data
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   XX = new double[totPts*nInputs_];
   (*X2) = new double[4*totPts];
   checkAllocate(*X2, "X2 in SVM::gen4DGridData");
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
                  XX[index*nInputs_+kk] = settings[kk]; 
               XX[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
               XX[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
               XX[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
               XX[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind3];
               (*X2)[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               (*X2)[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               (*X2)[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               (*X2)[index*4+3] = HX[3] * mm + lowerBounds_[ind3];
            }
         }
      }
   }
    
   YY = new double[totPts];
   checkAllocate(YY, "YY in SVM::gen4DGridData");
   if (outputLevel_ >= 1) printf("SVM interpolation begins....\n");
   SVMInterp(totPts, nInputs_, XX, YY, NULL);
   if (outputLevel_ >= 1) printf("SVM interpolation completed.\n");
   (*n) = totPts;
   (*Y2) = YY;
   delete [] XX;
   delete [] HX;
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double SVM::evaluatePoint(double *X)
{
   double Y=0.0;
#ifdef HAVE_SVM
   int    iOne=1;
   SVMInterp(iOne, nInputs_, X, &Y, NULL);
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double SVM::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_SVM
   SVMInterp(npts, nInputs_, X, Y, NULL);
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// Evaluate a given point return also the standard deviation 
// ------------------------------------------------------------------------
double SVM::evaluatePointFuzzy(double *X, double &std)
{
   double Y=0.0;
#ifdef HAVE_SVM
   SVMInterp(1, nInputs_, X, &Y, &std);
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return Y;
}

// ************************************************************************
// Evaluate a number of points and return also the standard deviation 
// ------------------------------------------------------------------------
double SVM::evaluatePointFuzzy(int npts,double *X,double *Y,double *Ystd)
{
#ifdef HAVE_SVM
   SVMInterp(npts, nInputs_, X, Y, Ystd);
#else
   printf("PSUADE ERROR : SVM not installed.\n");
#endif
   return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double SVM::setParams(int targc, char **targv)
{
   if (targc > 0) gamma_ = *(double *) targv[0];
   if (targc > 1) tolerance_ = *(double *) targv[1];
   if (targc > 2) kernel_ = *(int *) targv[2];
   return 0.0;
}

