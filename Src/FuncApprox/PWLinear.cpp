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
// Functions for the class PWLinear
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ------------------------------------------------------------------------
// local includes : class definition and utilities
// ------------------------------------------------------------------------
#include "Psuade.h"
#include "sysdef.h"
#include "pData.h"
#include "PsuadeUtil.h"
#include "PWLinear.h"

#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
PWLinear::PWLinear(int nInputs,int nSamples):FuncApprox(nInputs,nSamples)
{
   faID_ = PSUADE_RS_PWL;
   nROMs_ = nSamples;
   ROMStorePts_ = NULL;
   ROMStoreVals_ = NULL;
   ROMStoreGrads_ = NULL;
   ROMStoreEigens_ = NULL;
   threshold_ = 0;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PWLinear::~PWLinear()
{
   if (ROMStorePts_ != NULL) delete [] ROMStorePts_;
   if (ROMStoreVals_ != NULL) delete [] ROMStoreVals_;
   if (ROMStoreGrads_ != NULL) delete [] ROMStoreGrads_;
   if (ROMStoreEigens_ != NULL) delete [] ROMStoreEigens_;
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int PWLinear::initialize(double *X, double *Y)
{
   int     totPts, mm, iOne=1;
   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("PWLinear::genNDGridData ERROR - invalid argument.\n");
      exit(1);
   }
   if (nSamples_ <= nInputs_)
   {
      printf("PWLinear::genNDGridData INFO - not enough points.\n");
      return 0;
   }

   if (nInputs_ > 12)
   {
      printf("PWLinear::genNDGridData INFO - nInputs > 12.\n");
      printf("          No lattice points generated.\n");
      return 0;
   }

   if (ROMStorePts_ == NULL) setParams(iOne, NULL);
   if (psRSCodeGen_ == 1) 
      printf("PWLinear INFO: response surface stand-alone code not available.\n");
   return 0;
}

// ************************************************************************
// Generate lattice data based on the input set
// ------------------------------------------------------------------------
int PWLinear::genNDGridData(double *X, double *Y, int *N2, double **X2,
                            double **Y2)
{
   int     totPts, mm, iOne=1;

   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("PWLinear::genNDGridData ERROR - invalid argument.\n");
      exit(1);
   }
   if (nSamples_ <= nInputs_)
   {
      printf("PWLinear::genNDGridData INFO - not enough points.\n");
      return 0;
   }

   if (nInputs_ > 12)
   {
      printf("PWLinear::genNDGridData INFO - nInputs > 12.\n");
      printf("          No lattice points generated.\n");
      (*N2) = 0;
      (*X2) = 0;
      (*Y2) = 0;
      return 0;
   }

   if (ROMStorePts_ == NULL) setParams(iOne, NULL);
   if ((*N2) == -999) return 0;

   genNDGrid(N2, X2);
   if ((*N2) == 0) return 0;
   totPts = (*N2);

   (*Y2) = new double[totPts];
   checkAllocate(*Y2, "Y2 in PWLinear::genNDGridData");
   for (mm = 0; mm < totPts; mm++)
      (*Y2)[mm] = evaluatePoint(&((*X2)[mm*nInputs_]));

   return 0;
}

// ************************************************************************
// Generate 1D mesh results (setting others to some nominal values)
// ------------------------------------------------------------------------
int PWLinear::gen1DGridData(double *X, double *Y, int ind1,
                            double *settings, int *NN,
                            double **XX, double **YY)
{
   int    totPts, mm, nn;
   double HX, *Xloc;

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1);

   (*NN) = totPts;
   (*XX) = new double[totPts];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PWLinear::gen1DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn];

   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      Xloc[ind1] = HX * mm + lowerBounds_[ind1];
      (*XX)[mm] = Xloc[ind1];
      (*YY)[mm] = evaluatePoint(Xloc);
   }

   delete [] Xloc;
   return 0;
}

// ************************************************************************
// Generate 2D mesh results (setting others to some nominal values)
// ------------------------------------------------------------------------
int PWLinear::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                            double *settings, int *NN, double **XX, 
                            double **YY)
{
   int    totPts, mm, nn, index;
   double *HX, *Xloc;

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1);
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1);

   (*NN) = totPts;
   (*XX) = new double[totPts * 2];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PWLinear::gen2DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn];

   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         index = mm * nPtsPerDim_ + nn;
         Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
         Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
         (*XX)[index*2]   = Xloc[ind1];
         (*XX)[index*2+1] = Xloc[ind2];
         (*YY)[index] = evaluatePoint(Xloc);
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D mesh results (setting others to some nominal values)
// ------------------------------------------------------------------------
int PWLinear::gen3DGridData(double *X, double *Y, int ind1, int ind2,
                            int ind3, double *settings, int *NN,
                            double **XX, double **YY)
{
   int    totPts, mm, nn, pp, index;
   double *HX, *Xloc;

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1);
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1);
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1);

   (*NN) = totPts;
   (*XX) = new double[totPts * 3];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PWLinear::gen3DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn];

   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            index = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
            Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
            Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
            Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
            (*XX)[index*3]   = Xloc[ind1];
            (*XX)[index*3+1] = Xloc[ind2];
            (*XX)[index*3+2] = Xloc[ind3];
            (*YY)[index] = evaluatePoint(Xloc);
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D mesh results (setting others to some nominal values)
// ------------------------------------------------------------------------
int PWLinear::gen4DGridData(double *X, double *Y, int ind1, int ind2,
                            int ind3, int ind4, double *settings, int *NN,
                            double **XX, double **YY)
{
   int    totPts, mm, nn, pp, qq, index;
   double *HX, *Xloc;

   (*NN) = -999;
   genNDGridData(X, Y, NN, NULL, NULL);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1);
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1);
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1);
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1);

   (*NN) = totPts;
   (*XX) = new double[totPts * 4];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   checkAllocate(Xloc, "Xloc in PWLinear::gen4DGridData");
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn];

   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            for (qq = 0; qq < nPtsPerDim_; qq++)
            {
               index = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ + 
                       nn*nPtsPerDim_*nPtsPerDim_ + pp * nPtsPerDim_ + qq;
               Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
               Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
               Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
               Xloc[ind4] = HX[3] * qq + lowerBounds_[ind4];
               (*XX)[index*4]   = Xloc[ind1];
               (*XX)[index*4+1] = Xloc[ind2];
               (*XX)[index*4+2] = Xloc[ind3];
               (*XX)[index*4+3] = Xloc[ind4];
               (*YY)[index] = evaluatePoint(Xloc);
            }
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Evaluate a given point
// ------------------------------------------------------------------------
double PWLinear::evaluatePoint(double *X)
{
   int    iR, ii, count;
   double Y, Ytemp, dist2, diff, accum, dtemp, thresh, minDist2;

   thresh = threshold_ * 1.05;
   while (thresh < 2.0 * threshold_)
   {
      count = 0;
      accum = 0.0;
      Y = 0.0;
      minDist2 = PSUADE_UNDEFINED;
      for (iR = 0; iR < nROMs_; iR++)
      {
         Ytemp = ROMStoreVals_[iR];
         dist2 = 0.0;
         for (ii = 0; ii < nInputs_; ii++)
         {
            diff = X[ii] - ROMStorePts_[iR*nInputs_+ii]; 
            Ytemp += ROMStoreGrads_[iR*nInputs_+ii] * diff;
            dist2 += diff * diff;
         }
         dist2 = 0.5 * dist2 * ROMStoreEigens_[iR];
         if (dist2 <= 1.0e-12)
         { 
            Y = Ytemp;
            return Y;
         }
         if (dist2 < thresh)
         {
            dtemp = (thresh - dist2) / dist2;
            if (dtemp < 0.0) dtemp = 0.0;
            if (dtemp > 0.0)
            {
               Y += Ytemp * dtemp;
               accum += dtemp;
               count++;
            }
         }
         if (dist2 < minDist2) minDist2 = dist2;
      }
      if (count == 0)
      {
         printf("PWLinear ERROR: no ROV found, min dist2 = %e > %e.\n",
                minDist2, thresh);
         printf("PWLinear ERROR: no ROV found, reduce threshold to %e.\n",
                thresh);
         thresh = 1.1 * thresh;
      }
      else 
      {
         Y /= accum;
         return Y;
      }
   }
   printf("PWLinear ERROR: no ROV found for point, return 0.\n");
   printf("Data point is: \n");
   for (ii = 0; ii < nInputs_; ii++)
      printf("X %3d = %e\n", ii+1, X[ii]);
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double PWLinear::evaluatePoint(int npts, double *X, double *Y)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePoint(&X[kk*nInputs_]);
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviations
// ------------------------------------------------------------------------
double PWLinear::evaluatePointFuzzy(double *X, double &std)

{
   double Y;
   Y = evaluatePoint(X);
   std = 0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double PWLinear::evaluatePointFuzzy(int npts, double *X, double *Y,
                                    double *Ystd)
{
   for (int kk = 0; kk < npts; kk++)
      Y[kk] = evaluatePointFuzzy(&X[kk*nInputs_], Ystd[kk]);
   return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double PWLinear::setParams(int targc, char **targv)
{
   int        status, iR, ii, nOutputs, kk;
   double     *sampleInputs, *sampleOutputs;
   char       filename[500], lineIn[500];
   PsuadeData *psIO = NULL;
   pData      pPtr, pInputs, pOutputs; 

   if (targc == 1)
   {
      printf("PWLinear:: enter data file name : ");
      scanf("%s", filename);
      fgets(lineIn, 500, stdin);
      psIO = new PsuadeData();
      status = psIO->readPsuadeFile(filename);
      if (status != 0)
      {
         printf("ERROR: Problem reading file %s.\n", filename);
         exit(1);
      }
      psIO->getParameter("input_ninputs", pPtr);
      if (pPtr.intData_ != nInputs_)
      {
         printf("nInputs in file %s does not match with current file.\n",
                filename);
         exit(1);
      }
      psIO->getParameter("output_noutputs", pPtr);
      nOutputs = nInputs_ + nInputs_ * (nInputs_ + 1) / 2 + 2;
      if (pPtr.intData_ != nOutputs)
      {
         printf("nOutputs in file %s does not match with current file.\n",
                filename);
         exit(1);
      }
      psIO->getParameter("method_nsamples", pPtr);
      nROMs_ = pPtr.intData_;
      if (nROMs_ <= 0)
      {
         printf("ERROR: nROMs in file %s <= 0\n", filename);
         exit(1);
      }
      psIO->getParameter("input_sample", pInputs);
      ROMStorePts_ = new double[nROMs_*nInputs_];
      checkAllocate(ROMStorePts_,"ROMStorePts_ in PWLinear::setParams");
      for (iR = 0; iR < nROMs_*nInputs_; iR++)
         ROMStorePts_[iR] = pInputs.dbleArray_[iR];
      psIO->getParameter("output_sample", pOutputs);
      ROMStoreVals_ = new double[nROMs_];
      ROMStoreGrads_ = new double[nROMs_*nInputs_];
      checkAllocate(ROMStoreGrads_,"ROMStoreGrads_ in PWLinear::setParams");
      for (iR = 0; iR < nROMs_; iR++)
      {
         ROMStoreVals_[iR] = pOutputs.dbleArray_[iR*nOutputs];
         for (ii = 0; ii < nInputs_; ii++)
            ROMStoreGrads_[iR*nInputs_+ii] =
                                 pOutputs.dbleArray_[iR*nOutputs+ii+1];
      }
      ROMStoreEigens_ = new double[nROMs_];
      checkAllocate(ROMStoreEigens_,"ROMStoreEigens_ in PWLinear::setParams");
      for (iR = 0; iR < nROMs_; iR++)
         ROMStoreEigens_[iR] = pOutputs.dbleArray_[(iR+1)*nOutputs-1];
      psIO->getParameter("ana_threshold", pPtr);
      threshold_ = pPtr.dbleData_;
      if (threshold_ <= 0.0 || threshold_ >= 1.0)
      {
         printf("PWLinear ERROR: invalid analysis threshold %e.\n",threshold_);
      }
   }
   if (targc == 6)
   {
      printf("PWLinear: getting ROV information.\n");
      nROMs_ = *(int *) targv[0];
      nInputs_ = *(int *) targv[1];
      kk = *(int *) targv[2];
      sampleInputs = (double *) targv[3];
      sampleOutputs = (double *) targv[4];
      threshold_ = * (double *) targv[5];
      nOutputs = nInputs_ + nInputs_ * (nInputs_ + 1) / 2 + 2;
      if (kk != nOutputs)
      {
         printf("PWLinear ERROR: nOutputs mismatch (%d != %d).\n",
                kk, nOutputs);
         exit(1);
      }
      ROMStorePts_ = new double[nROMs_*nInputs_];
      checkAllocate(ROMStorePts_,"ROMStorePts_ in PWLinear::setParams");
      for (iR = 0; iR < nROMs_*nInputs_; iR++)
         ROMStorePts_[iR] = sampleInputs[iR];
      ROMStoreVals_ = new double[nROMs_];
      ROMStoreGrads_ = new double[nROMs_*nInputs_];
      checkAllocate(ROMStoreGrads_,"ROMStoreGrads_ in PWLinear::setParams");
      for (iR = 0; iR < nROMs_; iR++)
      {
         ROMStoreVals_[iR] = sampleOutputs[iR*nOutputs];
         for (ii = 0; ii < nInputs_; ii++)
            ROMStoreGrads_[iR*nInputs_+ii] =
                                 sampleOutputs[iR*nOutputs+ii+1];
      }
      ROMStoreEigens_ = new double[nROMs_];
      checkAllocate(ROMStoreEigens_,"ROMStoreEigens_ in PWLinear::setParams");
      for (iR = 0; iR < nROMs_; iR++)
         ROMStoreEigens_[iR] = sampleOutputs[(iR+1)*nOutputs-1];
      if (threshold_ <= 0.0 || threshold_ >= 1.0)
      {
         printf("PWLinear ERROR: invalid threshold %e.\n",threshold_);
         exit(1);
      }
   }
   if (psIO != NULL) delete psIO;
   return 0.0;
}

