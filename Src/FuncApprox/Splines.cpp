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
// Functions for the class Splines
// AUTHOR : CHARLES TONG
// DATE   : 2013
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Splines.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "PsuadeConfig.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

extern "C" 
{
  void spline1d(int,double*,double*,int,double*,double*);
  void spline2d(int,int,double*,double*,double*,int,
                double*,double*,double*);
  void spline3d(int,int,int,double*,double*,double*,double*,int,
                double*,double*,double*,double*);
}

// ************************************************************************
// Constructor for object class Mars
// ------------------------------------------------------------------------
Splines::Splines(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   double ddata;

   faID_ = PSUADE_RS_SPLINES;

   if (nSamples_ <= 0)
   {
      printf("Splines ERROR: nSamples <= 0.\n");
      exit(1);
   }
   if (nInputs <= 0 || nInputs > 3)
   {
      printf("Splines ERROR: this method currently only supports 1, 2,\n");
      printf("               or 3 inputs on a regular grid.\n");
      exit(1);
   }
   samIns_ = NULL;
   samOut_ = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Splines::~Splines()
{
   if (samIns_ != NULL) delete [] samIns_;
   if (samOut_ != NULL) delete [] samOut_;
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int Splines::initialize(double *XIn, double *YIn)
{
   genNDGridData(XIn, YIn, NULL, NULL, NULL);
   if (psRSCodeGen_ == 1) 
      printf("Splines INFO: response surface stand-alone code not available.\n");
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Splines::genNDGridData(double *XIn, double *YIn, int *, double **, 
                           double **)
{
   if (nInputs_ == 1)
   {
      int ii;
      for (ii = 1; ii < nSamples_; ii++)
      {
         if (XIn[ii] <= XIn[ii-1])
         {
            printf("Splines initialize ERROR: sample not factorial (1).\n");
            printf("        Suggestion: re-order the inputs first.\n");
            exit(1);
         }
      }
      n1_ = nSamples_;
      n2_ = 1;
      n3_ = 1;
   }
   else if (nInputs_ == 2)
   {
      int ii, jj;
      for (ii = 1; ii < nSamples_; ii++)
         if (XIn[2*ii+1] != XIn[2*(ii-1)+1]) break;
      n1_ = ii;
      n2_ = nSamples_ / n1_;
      if (n1_ * n2_ != nSamples_)
      {
         printf("Splines initialize ERROR: sample not factorial (2a).\n");
         printf("   Note: factorial should be laid out in order ");
         printf("beginning with input 1.\n");
         exit(1);
      }
      n3_ = 1;
      for (ii = 1; ii < n1_; ii++)
      {
         if (XIn[2*ii] <= XIn[2*(ii-1)])
         {
            printf("Splines initialize ERROR: sample not factorial (2b).\n");
            exit(1);
         }
      }
      for (ii = 1; ii < n2_; ii++)
      {
         for (jj = 0; jj < n1_; jj++)
         {
            if (XIn[2*(n1_*ii+jj)] != XIn[2*jj])
            {
               printf("Splines initialize ERROR: sample not factorial (2c).\n");
               exit(1);
            }
         }
      }
      for (ii = 1; ii < n2_; ii++)
      {
         if (XIn[2*ii*n1_+1] <= XIn[2*(ii-1)*n1_+1])
         {
            printf("Splines initialize ERROR: sample not factorial (2d).\n");
            exit(1);
         }
      }
      for (ii = 0; ii < n2_; ii++)
      {
         for (jj = 1; jj < n1_; jj++)
         {
            if (XIn[2*(ii*n1_+jj)+1] != XIn[2*(ii*n1_+jj-1)+1])
            {
               printf("Splines initialize ERROR: sample not factorial (2e).\n");
               exit(1);
            }
         }
      }
      printf("Splines: number of distinct values for input 1 = %d\n", n1_);
      printf("Splines: number of distinct values for input 2 = %d\n", n2_);
   } 
   else if (nInputs_ == 3)
   {
      int    ii, jj, kk;
      double *dinput1, *dinput2, *dinput3;
      dinput1 = new double[nSamples_];
      dinput2 = new double[nSamples_];
      dinput3 = new double[nSamples_];
      checkAllocate(dinput3, "dinput3 in Splines::genNDGridData");
      for (ii = 0; ii < nSamples_; ii++) dinput1[ii] = XIn[3*ii];
      for (ii = 0; ii < nSamples_; ii++) dinput2[ii] = XIn[3*ii+1];
      for (ii = 0; ii < nSamples_; ii++) dinput3[ii] = XIn[3*ii+2];
      for (ii = 1; ii < nSamples_; ii++)
         if (dinput2[ii] != dinput2[ii-1]) break;
      n1_ = ii;
      printf("Splines: number of distinct values for input 1 = %d\n", n1_);
      for (ii = 1; ii < nSamples_; ii++)
         if (dinput3[ii] != dinput3[ii-1]) break;
      n2_ = ii / n1_;
      n3_ = nSamples_ / n1_ / n2_;
      printf("Splines: number of distinct values for input 2 = %d\n", n2_);
      printf("Splines: number of distinct values for input 3 = %d\n", n3_);
      if (n1_ * n2_ * n3_ != nSamples_)
      {
         printf("Splines initialize ERROR: sample not factorial (3a).\n");
         printf("   Note: factorial should be laid out in order ");
         printf("beginning with input 1.\n");
         exit(1);
      }
      for (ii = 1; ii < n1_; ii++)
      {
         if (dinput1[ii] <= dinput1[ii-1])
         {
            printf("Splines initialize ERROR: sample not factorial (3b).\n");
            exit(1);
         }
      }
      for (kk = 1; kk < n2_*n3_; kk++)
      {
         for (ii = 0; ii < n1_; ii++)
         {
            if (dinput1[kk*n1_+ii] != dinput1[ii])
            {
               printf("Splines initialize ERROR: sample not factorial (3c).\n");
               exit(1);
            }
         }
      }
      for (ii = 1; ii < n2_; ii++)
      {
         if (dinput2[ii*n1_] <= dinput2[(ii-1)*n1_])
         {
            printf("Splines initialize ERROR: sample not factorial (3d).\n");
            exit(1);
         }
      }
      for (kk = 0; kk < n2_*n3_; kk++)
      {
         for (ii = 1; ii < n1_; ii++)
         {
            if (dinput2[kk*n1_+ii] != dinput2[kk*n1_])
            {
               printf("Splines initialize ERROR: sample not factorial (3e).\n");
               exit(1);
            }
         }
      }
      for (kk = 1; kk < n3_; kk++)
      {
         for (ii = 0; ii < n1_*n2_; ii++)
         {
            if (dinput2[kk*n1_*n2_+ii] != dinput2[ii])
            {
               printf("Splines initialize ERROR: sample not factorial (3f).\n");
               exit(1);
            }
         }
      }
      for (ii = 1; ii < n3_; ii++)
      {
         if (dinput3[ii*n1_*n2_] <= dinput3[(ii-1)*n1_*n2_])
         {
            printf("Splines initialize ERROR: sample not factorial (3g).\n");
            exit(1);
         }
      }
      for (kk = 0; kk < n3_; kk++)
      {
         for (ii = 1; ii < n1_*n2_; ii++)
         {
            if (dinput3[kk*n1_*n2_+ii] != dinput3[kk*n1_*n2_])
            {
               printf("Splines initialize ERROR: sample not factorial (3h).\n");
               exit(1);
            }
         }
      }
      delete [] dinput1;
      delete [] dinput2;
      delete [] dinput3;
   }

   if (samIns_ != NULL) delete [] samIns_;
   if (samOut_ != NULL) delete [] samOut_;
   samIns_ = new double[n1_+n2_+n3_];
   checkAllocate(samIns_, "samIns in Splines::genNDGridData");
   int ii;
   for (ii = 0; ii < n1_; ii++) samIns_[ii] = XIn[ii*nInputs_];
   if (nInputs_ >= 2)
   {
      for (ii = 0; ii < n2_; ii++)
         samIns_[n1_+ii] = XIn[nInputs_*ii*n1_+1];
   }
   if (nInputs_ == 3)
   {
      for (ii = 0; ii < n3_; ii++)
         samIns_[n1_+n2_+ii] = XIn[nInputs_*ii*n1_*n2_+2];
   }
   samOut_ = new double[nSamples_];
   checkAllocate(samOut_, "samOut in Splines::genNDGridData");
   for (ii = 0; ii < nSamples_; ii++) samOut_[ii] = YIn[ii];
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Splines::gen1DGridData(double *XIn,double *YIn,int ind1,double *setting,
                           int *NOut,double **XOut,double **YOut)
{
   int    ii;
   double *X1=NULL, *X2=NULL, *X3=NULL, *YY, HX;
   if (ind1 < 0 || ind1 > 2)
   {
      printf("Splines gen1DGrid ERROR: invalid input.\n");
      printf("                         Input 1 = %d\n", ind1+1);
      (*NOut) = 0;
      (*XOut) = NULL;
      (*YOut) = NULL;
      return -1;
   }
   genNDGridData(XIn, YIn, NOut, XOut, YOut);
   if ((*NOut) == -999) return 0;

   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   X1 = new double[nPtsPerDim_];
   X2 = new double[nPtsPerDim_];
   X3 = new double[nPtsPerDim_];
   YY = new double[nPtsPerDim_];
   checkAllocate(YY, "YY in Splines::gen1DGridData");
   if (ind1 == 0)
   {
      for (ii = 0; ii < nPtsPerDim_; ii++)
      {
         X1[ii] = HX * ii + lowerBounds_[ind1];
         if (nInputs_ > 1) X2[ii] = setting[1];
         if (nInputs_ > 2) X3[ii] = setting[2];
      }
      (*XOut) = X1;
   }
   else if (ind1 == 1)
   {
      for (ii = 0; ii < nPtsPerDim_; ii++)
      {
         X2[ii] = HX * ii + lowerBounds_[ind1];
         X1[ii] = setting[0];
         if (nInputs_ > 2) X3[ii] = setting[2];
      }
      (*XOut) = X2;
   }
   else
   {
      for (ii = 0; ii < nPtsPerDim_; ii++)
      {
         X3[ii] = HX * ii + lowerBounds_[ind1];
         X1[ii] = setting[0];
         X2[ii] = setting[1];
      }
      (*XOut) = X3;
   }

   if (nInputs_ == 1) spline1d(nSamples_, XIn, YIn, nPtsPerDim_, X1, YY); 
   else if (nInputs_ == 2)
   {
      spline2d(n1_, n2_, samIns_, &samIns_[n1_],YIn,nPtsPerDim_,X1,X2,YY); 
   }
   else
   {
      spline3d(n1_,n2_,n3_,samIns_,&samIns_[n1_],&samIns_[n1_+n2_],YIn,
               nPtsPerDim_,X1,X2,X3,YY); 
   }

   if (ind1 == 0)
   {
      delete [] X2;
      delete [] X3;
   }
   else if (ind1 == 1)
   {
      delete [] X1;
      delete [] X3;
   }
   else if (ind1 == 2)
   {
      delete [] X1;
      delete [] X2;
   }
   (*YOut) = YY;
   (*NOut) = nPtsPerDim_;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Splines::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                           double *settings,int *NOut,double **XOut, 
                           double **YOut)
{
   int    totPts, ii, jj, index;
   double *HX, *X1, *X2, *X3, *XX, *YY;

   if (nInputs_ < 2)
   {
      printf("Splines gen2DGrid ERROR: nInputs < 2.\n");
      (*NOut) = 0;
      (*XOut) = NULL;
      (*YOut) = NULL;
      return -1;
   }
   if (ind1 < 0 || ind1 > 2 || ind2 < 0 || ind2 > 2 || ind1 >= nInputs_ ||
       ind2 >= nInputs_)
   {
      printf("Splines gen2DGrid ERROR: invalid input.\n");
      printf("                         Input 1 = %d\n", ind1+1);
      printf("                         Input 2 = %d\n", ind2+1);
      (*NOut) = 0;
      (*XOut) = NULL;
      (*YOut) = NULL;
      return -1;
   }
   if (ind1 == ind2)
   {
      printf("Splines gen2DGrid ERROR: the two inputs should be distinct.\n");
      printf("                         Input 1 = %d\n", ind1+1);
      printf("                         Input 2 = %d\n", ind2+1);
      (*NOut) = 0;
      (*XOut) = NULL;
      (*YOut) = NULL;
      return -1;
   }
   genNDGridData(XIn, YIn, NOut, XOut, YOut);
   if ((*NOut) == -999) return 0;

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1);
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1);

   X1 = new double[totPts];
   X2 = new double[totPts];
   X3 = new double[totPts];
   checkAllocate(X3, "X3 in Splines::gen2DGridData");
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         X1[index] = settings[0];
         X2[index] = settings[1];
         X3[index] = settings[2];
         if (ind1 == 0)
            X1[index] = HX[0] * ii + lowerBounds_[ind1];
         else if (ind1 == 1)
            X2[index] = HX[0] * ii + lowerBounds_[ind1];
         else if (ind1 == 2)
            X3[index] = HX[0] * ii + lowerBounds_[ind1];

         if (ind2 == 0)
            X1[index] = HX[1] * jj + lowerBounds_[ind2];
         else if (ind2 == 1)
            X2[index] = HX[1] * jj + lowerBounds_[ind2];
         else if (ind2 == 2)
            X3[index] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   YY = new double[totPts];
   checkAllocate(YY, "YY in Splines::gen2DGridData");
   if (nInputs_ == 2)
   {
      spline2d(n1_, n2_, samIns_, &samIns_[n1_],YIn,totPts,X1,X2,YY); 
   }
   else
   {
      spline3d(n1_,n2_,n3_,samIns_,&samIns_[n1_],&samIns_[n1_+n2_],YIn,
               totPts,X1,X2,X3,YY); 
   }

   XX = new double[2*totPts];
   checkAllocate(XX, "XX in Splines::gen2DGridData");
   (*XOut) = XX;
   (*YOut) = YY;
   (*NOut) = totPts;
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         XX[2*index]   = HX[0] * ii + lowerBounds_[ind1];
         XX[2*index+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   delete [] X1;
   delete [] X2;
   delete [] X3;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Splines::gen3DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                           int ind3, double *settings, int *NOut, 
                           double **XOut, double **YOut)
{
   int    totPts, ii, jj, ll, index;
   double *HX, *X1, *X2, *X3, *XX, *YY;

   if (nInputs_ < 3)
   {
      printf("Splines gen3DGrid ERROR: nInputs < 3.\n");
      (*NOut) = 0;
      (*XOut) = NULL;
      (*YOut) = NULL;
      return -1;
   }
   if (ind1 < 0 || ind1 > 2 || ind2 < 0 || ind2 > 2 || ind3 < 0 || ind3 > 2)
   {
      printf("Splines gen3DGrid ERROR: invalid input.\n");
      printf("                         Input 1 = %d\n", ind1+1);
      printf("                         Input 2 = %d\n", ind2+1);
      printf("                         Input 3 = %d\n", ind3+1);
      (*NOut) = 0;
      (*XOut) = NULL;
      (*YOut) = NULL;
      return -1;
   }
   if (ind1 == ind2 || ind2 == ind3 || ind1 == ind3)
   {
      printf("Splines gen3DGrid ERROR: the three inputs should be distinct.\n");
      printf("                         Input 1 = %d\n", ind1+1);
      printf("                         Input 2 = %d\n", ind2+1);
      printf("                         Input 3 = %d\n", ind3+1);
      (*NOut) = 0;
      (*XOut) = NULL;
      (*YOut) = NULL;
      return -1;
   }
   genNDGridData(XIn, YIn, NOut, XOut, YOut);
   if ((*NOut) == -999) return 0;

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1);
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1);
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1);

   X1 = new double[totPts];
   X2 = new double[totPts];
   X3 = new double[totPts];
   checkAllocate(X3, "X3 in Splines::gen3DGridData");
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            X1[index] = settings[0];
            X2[index] = settings[1];
            X3[index] = settings[2];
            if (ind1 == 0)
               X1[index] = HX[0] * ii + lowerBounds_[ind1];
            else if (ind1 == 1)
               X2[index] = HX[0] * ii + lowerBounds_[ind1];
            else if (ind1 == 2)
               X3[index] = HX[0] * ii + lowerBounds_[ind1];

            if (ind2 == 0)
               X1[index] = HX[1] * jj + lowerBounds_[ind2];
            else if (ind2 == 1)
               X2[index] = HX[1] * jj + lowerBounds_[ind2];
            else if (ind2 == 2)
               X3[index] = HX[1] * jj + lowerBounds_[ind2];

            if (ind3 == 0)
               X1[index] = HX[2] * ll + lowerBounds_[ind3];
            else if (ind3 == 1)
               X2[index] = HX[2] * ll + lowerBounds_[ind3];
            else if (ind3 == 2)
               X3[index] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   YY = new double[totPts];
   checkAllocate(YY, "YY in Splines::gen3DGridData");
   spline3d(n1_,n2_,n3_,samIns_,&samIns_[n1_],&samIns_[n1_+n2_],YIn,totPts,
            X1,X2,X3,YY); 

   XX = new double[3*totPts];
   checkAllocate(XX, "XX in Splines::gen3DGridData");
   (*XOut) = XX;
   (*YOut) = YY;
   (*NOut) = totPts;
   for (ii = 0; ii < nPtsPerDim_; ii++)
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            XX[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            XX[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            XX[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   delete [] X1;
   delete [] X2;
   delete [] X3;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Splines::gen4DGridData(double *, double *, int, int, int, int, double *, 
                           int *NOut, double **XOut, double **YOut)
{
   printf("Splines ERROR: not relevant since nInputs <= 2.\n");
   (*NOut) = 0;
   (*XOut) = NULL;
   (*YOut) = NULL;
   return -1;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Splines::evaluatePoint(double *X)
{
   int    iOne=1;
   double Y=0.0;
   if (nInputs_ == 1)
      spline1d(nSamples_,samIns_,samOut_,iOne, X, &Y); 
   else if (nInputs_ == 2)
   {
      spline2d(n1_,n2_,samIns_,&samIns_[n1_],samOut_,iOne,X,&X[1],&Y);
   }
   else 
   {
      spline3d(n1_,n2_,n3_,samIns_,&samIns_[n1_],&samIns_[n1_+n2_],samOut_,
               iOne,X,&X[1],&X[2],&Y);
   }
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double Splines::evaluatePoint(int npts, double *X, double *Y)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePoint(&X[nInputs_*kk]);
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double Splines::evaluatePointFuzzy(double *X, double &std)
{
   double Y = evaluatePoint(X);
   std = 0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double Splines::evaluatePointFuzzy(int npts,double *X,double *Y,double *Ystd)
{
   evaluatePoint(npts, X, Y);
   for (int kk = 0; kk < npts; kk++) Ystd[kk] = 0.0;
   return 0.0;
}

