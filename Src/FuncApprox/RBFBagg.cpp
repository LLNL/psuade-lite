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
// Functions for the class RBFBagg
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "RBFBagg.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class RBFBagg
// ------------------------------------------------------------------------
RBFBagg::RBFBagg(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   int  ii, itmp;
   char pString[500];

   faID_ = PSUADE_RS_RBFB;

   numRBFs_ = 25;

   usageIndex_ = 4;

   if (outputLevel_ > 1)
   {
      printAsterisks(PL_INFO, 0);
      printf("*                RBFBagg Analysis\n");
      printDashes(PL_INFO, 0);
      printf("Number of instantiations   = %d\n",numRBFs_);
      printEquals(PL_INFO, 0);
   }
   if (psRSExpertMode_ == 1 && psInteractive_ == 1)
   {
      sprintf(pString, 
              "How many instantiation of RBF (10-5000, default=100) ? ");
      numRBFs_ = getInt(10, 5000, pString);
   }

   itmp = psRSExpertMode_;
   psRSExpertMode_ = 0;
   psInteractive_ = 0;
   rbfObjs_ = new RBF*[numRBFs_];
   PsuadeConfig *tmpConfig = psConfig_;
   psConfig_ = NULL;
   for (ii = 0; ii < numRBFs_; ii++) 
      rbfObjs_[ii] = new RBF(nInputs_, nSamples_);
   psConfig_ = tmpConfig;
   psRSExpertMode_ = itmp;
   psInteractive_ = 1;

   dataSetX_ = NULL;
   dataSetY_ = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
RBFBagg::~RBFBagg()
{
   int ii;

   if (rbfObjs_ != NULL) 
   {
      for (ii = 0; ii < numRBFs_; ii++) delete rbfObjs_[ii];
      delete [] rbfObjs_;
   }
}

// ************************************************************************
// Set lower and upper bounds 
// ------------------------------------------------------------------------
int RBFBagg::setBounds( double *lower, double *upper )
{
   for (int ii=0 ; ii<numRBFs_; ii++) 
      rbfObjs_[ii]->setBounds(lower, upper);
   return 0;
}

// ************************************************************************
// set number of points to generate in each dimension
// ------------------------------------------------------------------------
void RBFBagg::setNPtsPerDim(int npoints)
{
   nPtsPerDim_ = npoints;
   for (int ii=0 ; ii<numRBFs_; ii++) rbfObjs_[ii]->setNPtsPerDim(npoints);
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int RBFBagg::initialize(double *XX, double *Y)
{
   int    ii, ss, jj, index, *iCnts, expertFlag;
   double *XB, *YB;
   FILE   *fp;

   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   iCnts = new int[nSamples_];
   checkAllocate(iCnts, "iCnts in RBFBagg::initialize");
   expertFlag = psRSExpertMode_;
   psRSExpertMode_ = 0;
   for (ii = 0; ii < numRBFs_; ii++)
   {
      if (outputLevel_ >= 2)
         printf("RBFBagg::initialize : creating RBF #%d (of %d)\n",
                ii+1, numRBFs_);
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++) iCnts[ss] = usageIndex_ * 2;
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = -1;
            while (index == -1)
            {
               index = PSUADE_rand() % nSamples_;
               if (iCnts[index] > 0 && (iCnts[index] % usageIndex_) == 0)
               {
                  iCnts[index]--;
               }
               else if (iCnts[index] > 0 && (iCnts[index] % usageIndex_) != 0)
               {
                  iCnts[index]--;
                  index = -1;
               }
               else if (iCnts[index] <= 0) index = -1;
            }
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = XX[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
         index = 0;
         for (ss = 0; ss < nSamples_; ss++) 
            if (iCnts[ss] < usageIndex_*2) index++;
         if (outputLevel_ >= 2)
            printf("     Number of sample points used = %d (out of %d)\n",
                   index,nSamples_);
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      rbfObjs_[ii]->setOutputLevel(0);
      jj = 0;
      for (ss = 0; ss < nSamples_; ss++)
         if (YB[ss] >= PSUADE_UNDEFINED) jj++;
      if (jj > 0)
      {
         printf("RBFBagg ERROR: some of the sample outputs are undefined.\n");
         exit(1);
      }
      rbfObjs_[ii]->initialize(XB, YB);
      fp = fopen("ps_print", "r");
      if (fp != NULL)
      {
         printf("RBFBagg: set print level to 2\n");
         outputLevel_ = 2;
         fclose(fp);
      }
   }
   delete [] iCnts;
   psRSExpertMode_ = expertFlag;
   delete [] XB;
   delete [] YB;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int RBFBagg::genNDGridData(double *XX, double *Y, int *N, double **XX2, 
                            double **Y2)
{
   int    totPts, ii, ss, jj, index, *iCnts, expertFlag;
   double *XXt, *Yt, *XB, *YB;
   FILE   *fp;

   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   checkAllocate(YB, "YB in RBFBagg::genNDGridData");
   if ((*N) == -999)
   {
      iCnts = new int[nSamples_];
      expertFlag = psRSExpertMode_;
      psRSExpertMode_ = 0;
      for (ii = 0; ii < numRBFs_; ii++)
      {
         if (outputLevel_ >= 2)
            printf("RBFBagg::genNDGridData : creating RBF_ #%d (of %d)\n",
                   ii+1, numRBFs_);
         if (dataSetX_ == NULL)
         {
            for (ss = 0; ss < nSamples_; ss++) iCnts[ss] = usageIndex_ * 2;
            for (ss = 0; ss < nSamples_; ss++)
            {
               index = -1;
               while (index == -1)
               {
                  index = PSUADE_rand() % nSamples_;
                  if (iCnts[index] > 0 && (iCnts[index] % usageIndex_) == 0)
                  {
                     iCnts[index]--;
                  }
                  else if (iCnts[index] > 0 && (iCnts[index] % usageIndex_) != 0)
                  {
                     iCnts[index]--;
                     index = -1;
                  }
                  else if (iCnts[index] <= 0) index = -1;
               }
               for (jj = 0; jj < nInputs_; jj++)
                  XB[ss*nInputs_+jj] = XX[index*nInputs_+jj]; 
               YB[ss] = Y[index]; 
            }
            index = 0;
            for (ss = 0; ss < nSamples_; ss++) 
               if (iCnts[ss] < usageIndex_*2) index++;
            if (outputLevel_ >= 2)
               printf("     Number of sample points used = %d (out of %d)\n",
                      index,nSamples_);
         }
         else
         {
            for (ss = 0; ss < nSamples_; ss++)
            {
               for (jj = 0; jj < nInputs_; jj++)
                  XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
               YB[ss] = dataSetY_[ii][ss]; 
            }
         }
         rbfObjs_[ii]->setOutputLevel(0);
         rbfObjs_[ii]->genNDGridData(XB, YB, N, NULL, NULL);
         fp = fopen("ps_print", "r");
         if (fp != NULL)
         {
            printf("RBFBagg: set print level to 2\n");
            outputLevel_ = 2;
            fclose(fp);
         }
      }
      delete [] iCnts;
      psRSExpertMode_ = expertFlag;
   }
   else
   {
      expertFlag = psRSExpertMode_;
      psRSExpertMode_ = 0;
      totPts = nPtsPerDim_;
      for (ii = 1; ii < nInputs_; ii++) totPts = totPts * nPtsPerDim_;
      (*XX2) = new double[nInputs_ * totPts];
      (*Y2)  = new double[totPts];
      checkAllocate(*Y2, "Y2 in RBFBagg::genNDGridData");
      for (ii = 0; ii < totPts; ii++) (*Y2)[ii] = 0.0;

      for (ii = 0; ii < numRBFs_; ii++)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = XX[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
         rbfObjs_[ii]->genNDGridData(XB, YB, N, &XXt, &Yt);

         for (ss = 0; ss < totPts; ss++) (*Y2)[ss] += Yt[ss];

         if (ii == numRBFs_-1)
            for (ss = 0; ss < totPts*nInputs_; ss++) (*XX2)[ss] = XXt[ss];
         delete [] XXt;
         delete [] Yt;
      }
      (*N) = totPts;

      for (ss = 0; ss < totPts; ss++) (*Y2)[ss] /= (double) numRBFs_;
      psRSExpertMode_ = expertFlag;
   }
   delete [] XB;
   delete [] YB;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int RBFBagg::gen1DGridData(double *X, double *Y, int ind1, 
                            double *settings, int *N, double **XX, 
                            double **YY)
{
   int    totPts, ii, ss, jj, index, expertFlag;
   double *XB, *YB, *XXt, *Yt;

   expertFlag = psRSExpertMode_;
   psRSExpertMode_ = 0;
   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[totPts];
   (*YY)  = new double[totPts];
   checkAllocate(*YY, "YY in RBFBagg::gen1DGridData");
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   for (ii = 0; ii < numRBFs_; ii++)
   {
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      rbfObjs_[ii]->gen1DGridData(XB, YB, ind1, settings, N, &XXt, &Yt);

      for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];
      if (ii == numRBFs_-1)
         for (ss = 0; ss < totPts; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }
   for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numRBFs_;
   delete [] XB;
   delete [] YB;
   psRSExpertMode_ = expertFlag;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int RBFBagg::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                            double *settings, int *N, double **XX, 
                            double **YY)
{
   int    totPts, ii, ss, jj, index, expertFlag;
   double *XB, *YB, *XXt, *Yt;

   expertFlag = psRSExpertMode_;
   psRSExpertMode_ = 0;
   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_ * nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[2 * totPts];
   (*YY)  = new double[totPts];
   checkAllocate(*YY, "YY in RBFBagg::gen2DGridData");
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   for (ii = 0; ii < numRBFs_; ii++)
   {
      if (outputLevel_ >= 1)
         printf("RBFBagg::gen2DGridData : creating RBF #%d (of %d)\n",
                ii+1, numRBFs_);
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      rbfObjs_[ii]->gen2DGridData(XB, YB, ind1, ind2, settings, N, 
                                   &XXt, &Yt);

      for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];

      if (ii == numRBFs_-1)
         for (ss = 0; ss < 2*totPts; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }
   for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numRBFs_;
   delete [] XB;
   delete [] YB;
   psRSExpertMode_ = expertFlag;
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int RBFBagg::gen3DGridData(double *X, double *Y, int ind1, int ind2, 
                            int ind3, double *settings, int *N, double **XX, 
                            double **YY)
{
   int    totPts, ii, ss, jj, index, expertFlag;
   double *XB, *YB, *XXt, *Yt;

   expertFlag = psRSExpertMode_;
   psRSExpertMode_ = 0;
   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[3 * totPts];
   (*YY)  = new double[totPts];
   checkAllocate(*YY, "YY in RBFBagg::gen3DGridData");
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   for (ii = 0; ii < numRBFs_; ii++)
   {
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      rbfObjs_[ii]->gen3DGridData(XB, YB, ind1, ind2, ind3, settings, 
                                   N, &XXt, &Yt);

      for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];

      if (ii == numRBFs_-1)
         for (ss = 0; ss < 3*totPts; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }

   for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numRBFs_;

   delete [] XB;
   delete [] YB;
   psRSExpertMode_ = expertFlag;
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int RBFBagg::gen4DGridData(double *X, double *Y, int ind1, int ind2, 
                            int ind3, int ind4, double *settings, int *N, 
                            double **XX, double **YY)
{
   int    totPts, ii, ss, jj, index, expertFlag;
   double *XB, *YB, *XXt, *Yt;

   expertFlag = psRSExpertMode_;
   psRSExpertMode_ = 0;
   XB = new double[nInputs_ * nSamples_];
   YB = new double[nSamples_];
   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   (*N)   = totPts;
   (*XX)  = new double[4 * totPts];
   (*YY)  = new double[totPts];
   checkAllocate(*YY, "YY in RBFBagg::gen4DGridData");
   for (ii = 0; ii < totPts; ii++) (*YY)[ii] = 0.0;

   for (ii = 0; ii < numRBFs_; ii++)
   {
      if (dataSetX_ == NULL)
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            index = PSUADE_rand() % nSamples_;
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = X[index*nInputs_+jj]; 
            YB[ss] = Y[index]; 
         }
      }
      else
      {
         for (ss = 0; ss < nSamples_; ss++)
         {
            for (jj = 0; jj < nInputs_; jj++)
               XB[ss*nInputs_+jj] = dataSetX_[ii][ss*nInputs_+jj]; 
            YB[ss] = dataSetY_[ii][ss]; 
         }
      }
      rbfObjs_[ii]->gen4DGridData(XB, YB, ind1, ind2, ind3, ind4,
                                   settings, N, &XXt, &Yt);

      for (ss = 0; ss < totPts; ss++) (*YY)[ss] += Yt[ss];

      if (ii == numRBFs_-1)
         for (ss = 0; ss < 4*totPts; ss++) (*XX)[ss] = XXt[ss];
      delete [] XXt;
      delete [] Yt;
   }

   for (ss = 0; ss < totPts; ss++) (*YY)[ss] /= (double) numRBFs_;

   delete [] XB;
   delete [] YB;
   psRSExpertMode_ = expertFlag;
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double RBFBagg::evaluatePoint(double *X)
{
   double Y=0.0;
   int    ii;
   double Yt;

   for (ii = 0; ii < numRBFs_; ii++) 
   {
      Yt = rbfObjs_[ii]->evaluatePoint(X);
      Y += Yt;
   }
   Y /= (double) numRBFs_;
   return Y;
}

// ************************************************************************
// Evaluate a number of points
// ------------------------------------------------------------------------
double RBFBagg::evaluatePoint(int npts, double *X, double *Y)
{
   int    in, ii;
   double YY, Yt;

   for (in = 0; in < npts; in++) 
   {
      YY = 0.0;
      for (ii = 0; ii < numRBFs_; ii++) 
      {
         Yt = rbfObjs_[ii]->evaluatePoint(&(X[in*nInputs_]));
         YY += Yt;
      }
      Y[in] = YY / (double) numRBFs_;
   }
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double RBFBagg::evaluatePointFuzzy(double *X, double &std)
{
   double Ymean=0.0;
   int    ii;
   double *Y;

   Y = new double[numRBFs_];
   checkAllocate(Y, "Y in RBFBagg::evaluatePointFuzzy");
   for (ii = 0; ii < numRBFs_; ii++) 
   {
      Y[ii] = rbfObjs_[ii]->evaluatePoint(X);
      Ymean += Y[ii];
   }
   Ymean /= (double) numRBFs_;
   std = 0.0;
   for (ii = 0; ii < numRBFs_; ii++) 
      std += (Y[ii] - Ymean) * (Y[ii] - Ymean);
   std = sqrt(std / (double) numRBFs_);
   delete [] Y;
   return Ymean;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double RBFBagg::evaluatePointFuzzy(int npts, double *X, double *Y,
                                    double *Ystd)
{
   int    in, ii;
   double *YY;

   YY = new double[numRBFs_];
   checkAllocate(YY, "YY in RBFBagg::evaluatePointFuzzy");
   for (in = 0; in < npts; in++) 
   {
      for (ii = 0; ii < numRBFs_; ii++) 
         YY[ii] = rbfObjs_[ii]->evaluatePoint(&(X[in*nInputs_]));
      Y[in] = 0.0;
      for (ii = 0; ii < numRBFs_; ii++) Y[in] += YY[ii];
      Y[in] /= (double) numRBFs_;
      Ystd[in] = 0.0;
      for (ii = 0; ii < numRBFs_; ii++) 
         Ystd[in] += (YY[ii] - Y[in]) * (YY[ii] - Y[in]);
      Ystd[in] = sqrt(Ystd[in] / (double) numRBFs_);
   }
   delete [] YY;
   return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double RBFBagg::setParams(int targc, char **targv)
{
   int    ii, itmp, leng;
   double *Xdata, *Ydata;
   char   cString[500], *argv[3];

   if (targc == 2 && !strcmp(targv[0], "num_rbf"))
   {
      if (rbfObjs_ != NULL) 
      {
         for (ii = 0; ii < numRBFs_; ii++) delete rbfObjs_[ii];
         delete [] rbfObjs_;
      }
      numRBFs_ = *(int *) targv[1];
      if (numRBFs_ < 2) numRBFs_ = 2;
      printf("RBF with bagging: no. of RBF set to = %d.\n", numRBFs_);
      itmp = psRSExpertMode_;
      psRSExpertMode_ = 0;
      rbfObjs_ = new RBF*[numRBFs_];
      for (ii = 0; ii < numRBFs_; ii++) 
         rbfObjs_[ii] = new RBF(nInputs_, nSamples_);
      psRSExpertMode_ = itmp;
   }
   else if (targc == 5 && !strcmp(targv[0], "rbf_sample"))
   {
      itmp = *(int *) targv[1];
      if (itmp < 0 || itmp >= numRBFs_)
      {
         printf("RBFBagg ERROR: in loading sample - invalid index.\n");
         exit(1);
      }
      leng = *(int *) targv[2];
      if (leng != nSamples_)
      {
         printf("RBFBagg ERROR: in loading sample - nSamples mismatch.\n");
         exit(1);
      }
      Xdata = (double *) targv[3];
      Ydata = (double *) targv[4];
      if (dataSetX_ == NULL)
      {
         dataSetX_ = new double*[numRBFs_];
         checkAllocate(dataSetX_, "dataSetX_ in RBFBagg::setParams");
         for (ii = 0; ii < numRBFs_; ii++) dataSetX_[ii] = NULL;
         dataSetY_ = new double*[numRBFs_];
         checkAllocate(dataSetY_, "dataSetY_ in RBFBagg::setParams");
         for (ii = 0; ii < numRBFs_; ii++) dataSetY_[ii] = NULL;
      }
      if (dataSetX_[itmp] == NULL)
      {
         dataSetX_[itmp] = new double[leng*nInputs_];
         dataSetY_[itmp] = new double[leng];
         checkAllocate(dataSetY_[itmp], "dataSetY_[itmp] in RBFBagg::setParams");
      }
      for (ii = 0; ii < leng*nInputs_; ii++) dataSetX_[itmp][ii] = Xdata[ii];
      for (ii = 0; ii < leng; ii++) dataSetY_[itmp][ii] = Ydata[ii];
   }
   else
   {
      printf("RBFBagg setParams ERROR: invalid command %s.\n", targv[0]);
   }
   return 0.0;
}

