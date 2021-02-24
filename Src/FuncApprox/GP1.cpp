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
// Functions for the class GP1
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "GP1.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "Vector.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#ifdef HAVE_TPROS
extern "C" 
{
  void TprosTrain(int nInputs, int nTrains, double *trainInputs,
                  double *trainOutput, int, double *, double *);

  void TprosInterp(int nTests, double *inputs, double *output, double *stds);

  void TprosGetLengthScales(int nInputs, double *lengthScales);
}
#endif

// ************************************************************************
// Constructor for object class GP1
// ------------------------------------------------------------------------
GP1::GP1(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   faID_ = PSUADE_RS_GP1;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GP1::~GP1()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int GP1::initialize(double *XIn, double *YIn)
{
#ifdef HAVE_TPROS
  int    ss, ii;
  char   pString[500], response[500], *cString;
  psVector vecStds, vecX, vecY;

  response[0] = 'n';
  if (psRSExpertMode_ != 1 && psConfig_ != NULL)
  {
    cString = psConfig_->getParameter("normalize_outputs");
    if (cString != NULL) response[0] = 'y';
  }
  if (psRSExpertMode_ == 1)
  {
    sprintf(pString, "GP1: normalize output? (y or n) ");
    getString(pString, response);
  }
   
  vecX.setLength(nSamples_*nInputs_);
  initInputScaling(XIn, vecX.getDVector(), 0);
  vecY.setLength(nSamples_);
  if (response[0] == 'y')
  {
    initOutputScaling(YIn, vecY.getDVector());
  }
  else
  {
    for (ii = 0; ii < nSamples_; ii++) vecY[ii] = YIn[ii];
    YMean_ = 0.0;
    YStd_ = 1.0;
  }
   
  psVector vecLScales;
  vecStds.setLength(nSamples_);
  if (outputLevel_ >= 1) printf("GP1 training begins....\n");
  TprosTrain(nInputs_,nSamples_,vecX.getDVector(),vecY.getDVector(),0,NULL, 
             vecStds.getDVector());
  for (ss = 0; ss < nSamples_; ss++) vecStds[ss] = 0.0;
  if (outputLevel_ >= 1) printf("GP1 training completed.\n");
  if (psRSExpertMode_)
  {
    vecLScales.setLength(nInputs_);
    TprosGetLengthScales(nInputs_, vecLScales.getDVector());
    printf("GP1 training information: \n");
    for (ii = 0; ii < nInputs_; ii++)
       printf("Input %d mean,std,length scale = %e %e %e\n",
              ii+1, XMeans_[ii], XStds_[ii], vecLScales[ii]); 
  }
  if (psRSCodeGen_ == 1) 
    printf("GP1 INFO: response surface stand-alone code not available.\n");
  return 0;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
  return -1;
#endif
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int GP1::genNDGridData(double *XIn, double *YIn, int *NOut, double **XOut, 
                      double **YOut)
{
#ifdef HAVE_TPROS
  int    totPts, ii, jj;
  double *XX;
  psVector vecXT;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  genNDGrid(NOut, &XX);
  if ((*NOut) == 0) return 0;
  totPts = (*NOut);

  if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in GP1::genNDGrid");
  vecXT.setLength(totPts*nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
  {
    for (jj = 0; jj < totPts; jj++)
      vecXT[jj*nInputs_+ii] = (XX[jj*nInputs_+ii] - XMeans_[ii]) /
                              XStds_[ii];
  } 
  TprosInterp(totPts, vecXT.getDVector(), *YOut, NULL);
  for (ii = 0; ii < totPts; ii++)
    (*YOut)[ii] = ((*YOut)[ii] * YStd_) + YMean_;
  if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
  (*XOut) = XX;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
  return -1;
#endif
  return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int GP1::gen1DGridData(double *XIn, double *YIn, int ind1,double *settings, 
                       int *NOut, double **XOut, double **YOut)
{
#ifdef HAVE_TPROS
  int    totPts, ii, jj, kk;
  double HX;
  psVector vecXT;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
 
  totPts = nPtsPerDim_;
  HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

  (*XOut) = new double[totPts];
  checkAllocate(*XOut, "XOut in GP1::gen1DGrid");
  vecXT.setLength(totPts*nInputs_);
  for (ii = 0; ii < nPtsPerDim_; ii++) 
  {
    for (kk = 0; kk < nInputs_; kk++) 
      vecXT[ii*nInputs_+kk] = settings[kk]; 
    vecXT[ii*nInputs_+ind1] = HX * ii + lowerBounds_[ind1];
    (*XOut)[ii] = HX * ii + lowerBounds_[ind1];
  }
    
  for (ii = 0; ii < nInputs_; ii++)
  {
    for (jj = 0; jj < totPts; jj++)
      vecXT[jj*nInputs_+ii] = (vecXT[jj*nInputs_+ii] - XMeans_[ii]) /
                               XStds_[ii];
  } 
  if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
  (*YOut) = new double[totPts];
  checkAllocate(*YOut, "YOut in GP1::gen1DGrid");
  TprosInterp(totPts, vecXT.getDVector(), *YOut, NULL);
  for (ii = 0; ii < totPts; ii++) 
    (*YOut)[ii] = ((*YOut)[ii] * YStd_) + YMean_;
  if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
  (*NOut) = totPts;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return 0;
}

// ************************************************************************
// Generate 2D results for display
// ------------------------------------------------------------------------
int GP1::gen2DGridData(double *XIn, double *YIn, int ind1, int ind2, 
                 double *settings, int *NOut, double **XOut, double **YOut)
{
#ifdef HAVE_TPROS
  int totPts, ii, jj, kk, index;
  psVector vecHX, vecXT;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(2);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1])/(nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2])/(nPtsPerDim_ - 1); 

  vecXT.setLength(totPts*nInputs_);
  (*XOut) = new double[2*totPts];
  checkAllocate(*XOut, "XOut in GP1::gen2DGrid");
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
  checkAllocate(*YOut, "YOut in GP1::gen2DGrid");
  for (ii = 0; ii < nInputs_; ii++)
  {
    for (jj = 0; jj < totPts; jj++)
      vecXT[jj*nInputs_+ii] = (vecXT[jj*nInputs_+ii] - XMeans_[ii]) /
                              XStds_[ii];
  } 
  if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
  TprosInterp(totPts, vecXT.getDVector(), *YOut, NULL);
  if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
  for (ii = 0; ii < totPts; ii++)
    (*YOut)[ii] = ((*YOut)[ii] * YStd_) + YMean_;
  (*NOut) = totPts;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int GP1::gen3DGridData(double *XIn,double *YIn,int ind1,int ind2,int ind3,
              double *settings, int *NOut, double **XOut, double **YOut)
{
#ifdef HAVE_TPROS
  int    totPts, ii, jj, kk, ll, index;
  psVector vecHX, vecXT;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(3);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

  vecXT.setLength(totPts*nInputs_);
  (*XOut) = new double[3*totPts];
  checkAllocate(*XOut, "XOut in GP1::gen3DGrid");
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
  checkAllocate(*YOut, "YOut in GP1::gen3DGrid");
  for (ii = 0; ii < nInputs_; ii++)
  {
    for (jj = 0; jj < totPts; jj++)
      vecXT[jj*nInputs_+ii] = (vecXT[jj*nInputs_+ii] - XMeans_[ii]) /
                              XStds_[ii];
  } 
  if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
  TprosInterp(totPts, vecXT.getDVector(), *YOut, NULL);
  if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
  for (ii = 0; ii < totPts; ii++)
    (*YOut)[ii] = ((*YOut)[ii] * YStd_) + YMean_;
  (*NOut) = totPts;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int GP1::gen4DGridData(double *XIn,double *YIn,int ind1,int ind2, int ind3,
                       int ind4,double *settings, int *NOut, double **XOut, 
                       double **YOut)
{
#ifdef HAVE_TPROS
  int    totPts, ii, jj, kk, ll, mm, index;
  psVector vecXT, vecHX;

  initialize(XIn, YIn);
  if ((*NOut) == -999 || XOut == NULL || YOut == NULL) return 0;
  
  totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
  vecHX.setLength(4);
  vecHX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
  vecHX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
  vecHX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
  vecHX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

  vecXT.setLength(totPts*nInputs_);
  (*XOut) = new double[4*totPts];
  checkAllocate(*XOut, "XOut in GP1::gen4DGrid");
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
  checkAllocate(*YOut, "YOut in GP1::gen4DGrid");
  for (ii = 0; ii < nInputs_; ii++)
  {
    for (jj = 0; jj < totPts; jj++)
      vecXT[jj*nInputs_+ii] = (vecXT[jj*nInputs_+ii] - XMeans_[ii]) /
                              XStds_[ii];
  } 
  if (outputLevel_ >= 1) printf("GP1 interpolation begins....\n");
  TprosInterp(totPts, vecXT.getDVector(), *YOut, NULL);
  if (outputLevel_ >= 1) printf("GP1 interpolation completed.\n");
  for (ii = 0; ii < totPts; ii++)
    (*YOut)[ii] = ((*YOut)[ii] * YStd_) + YMean_;
  (*NOut) = totPts;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double GP1::evaluatePoint(double *X)
{
  double Y=0.0;
#ifdef HAVE_TPROS
  int    ii, iOne=1;
  psVector vecXT;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  vecXT.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    vecXT[ii] = (X[ii] - XMeans_[ii]) / XStds_[ii];
  TprosInterp(iOne, vecXT.getDVector(), &Y, NULL);
  Y = Y * YStd_ + YMean_;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double GP1::evaluatePoint(int npts, double *X, double *Y)
{
#ifdef HAVE_TPROS
  int ii, jj;
  psVector vecXT;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  vecXT.setLength(npts*nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    for (jj = 0; jj < npts; jj++)
      vecXT[jj*nInputs_+ii] = (X[jj*nInputs_+ii]-XMeans_[ii])/XStds_[ii];
  TprosInterp(npts, vecXT.getDVector(), Y, NULL);
  for (jj = 0; jj < npts; jj++)
    Y[jj] = Y[jj] * YStd_ + YMean_;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double GP1::evaluatePointFuzzy(double *X, double &std)
{
  double Y=0.0;
#ifdef HAVE_TPROS
  int ii;
  psVector vecXT;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  vecXT.setLength(nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    vecXT[ii] = (X[ii] - XMeans_[ii]) / XStds_[ii];
  TprosInterp(1, vecXT.getDVector(), &Y, &std);
  Y = Y * YStd_ + YMean_;
  if (std < 0) printf("GP1 ERROR: variance < 0\n");
  else         std = sqrt(std) * YStd_;
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return Y;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double GP1::evaluatePointFuzzy(int npts,double *X, double *Y, double *Ystd)
{
#ifdef HAVE_TPROS
  int    ii, jj;
  psVector vecXT;
  if (XMeans_ == NULL)
  {
    printf("PSUADE ERROR : not initialized yet.\n");
    exit(1);
  }
  vecXT.setLength(npts*nInputs_);
  for (ii = 0; ii < nInputs_; ii++)
    for (jj = 0; jj < npts; jj++)
      vecXT[jj*nInputs_+ii] = (X[jj*nInputs_+ii]-XMeans_[ii])/XStds_[ii];
  TprosInterp(npts, vecXT.getDVector(), Y, Ystd);
  for (int ii = 0; ii < npts; ii++)
  {
    Y[ii] = Y[ii] * YStd_ + YMean_;
    if (Ystd[ii] < 0) printf("GP1 ERROR: variance < 0\n");
    else              Ystd[ii] = sqrt(Ystd[ii]) * YStd_;
  }
#else
  printf("PSUADE ERROR : GP1 not installed.\n");
#endif
  return 0.0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double GP1::setParams(int targc, char **targv)
{
  int    ii, ind;
  double mmax, ddata=0.0, range;
  char   pString[500];
  FILE   *fp;
  psVector  vecLScales;
  psIVector ivecT;
                                                                                
  if (targc > 0 && !strcmp(targv[0], "rank"))
  {
    vecLScales.setLength(nInputs_);
#ifdef HAVE_TPROS
    TprosGetLengthScales(nInputs_, vecLScales.getDVector());
#else
    printf("PSUADE ERROR : GP1 not installed.\n");
    return 0.0;
#endif
    mmax = 0.0;
    for (ii = 0; ii < nInputs_; ii++)
    {
      vecLScales[ii] = 1.0/vecLScales[ii];
      if (XMeans_[ii] == 0 && XStds_[ii] == 1)
      {
        range = upperBounds_[ii] - lowerBounds_[ii];
        vecLScales[ii] *= range;
      }
      if (vecLScales[ii] > mmax) mmax = vecLScales[ii];
    }
    for (ii = 0; ii < nInputs_; ii++)
      vecLScales[ii] = vecLScales[ii] / mmax * 100.0;
    if (psPlotTool_ == 1)
         fp = fopen("scilabgpsa.sci", "w");
    else fp = fopen("matlabgpsa.m", "w");
    if (fp == NULL)
    {
      printf("GP1 ERROR: something wrong with opening a write file.\n");
    }
    else
    {
      fprintf(fp, "n = %d;\n", nInputs_);
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs_; ii++)
        fprintf(fp, "%24.16e \n", PABS(vecLScales[ii]));
      fprintf(fp, "]; \n");
      fprintf(fp, "ymax = max(Y);\n");
      fprintf(fp, "ymin = 0;\n");
      fprintf(fp, "if (ymax == ymin)\n");
      fprintf(fp, "   ymax = ymax * 0.1;\n");
      fprintf(fp, "end;\n");
      fwritePlotCLF(fp);
      fprintf(fp, "bar(Y,0.8);\n");
      fwritePlotAxes(fp);
      sprintf(pString, "GP Ranking");
      fwritePlotTitle(fp, pString);
      sprintf(pString, "Input Numbers");
      fwritePlotXLabel(fp, pString);
      sprintf(pString, "GP Measure");
      fwritePlotYLabel(fp, pString);
      if (psPlotTool_ == 1)
      {
        fprintf(fp,"a=gca();\n");
        fprintf(fp,
           "a.data_bounds=[0, ymin; n+1, ymax+0.01*(ymax-ymin)];\n");
      }
      else
      {
        fprintf(fp,"axis([0 n+1 ymin ymax+0.01*(ymax-ymin)])\n");
      }
      fclose(fp);
      if (psPlotTool_ == 1)
           printf("GP ranking in file scilabgpsa.sci\n");
      else printf("GP ranking in file matlabgpsa.m\n");
    }
    ivecT.setLength(nInputs_);
    for (ii = 0; ii < nInputs_; ii++) ivecT[ii] = ii;
    sortDbleList2a(nInputs_, vecLScales.getDVector(), ivecT.getIVector());
    if (targc == 1)
    {
      printAsterisks(PL_INFO, 0);
      printf("* GP1 screening rankings \n");
      printAsterisks(PL_INFO, 0);
      for (ii = nInputs_-1; ii >= 0; ii--)
        printf("*  Rank %3d : Input = %3d (score = %5.1f) (ref = %e)\n", 
               nInputs_-ii, ivecT[ii]+1, vecLScales[ii], 
               vecLScales[ii]*mmax*0.01);
      printAsterisks(PL_INFO, 0);
    }
    if (targc > 1)
    {
      ind = *(int *) targv[1];
      if (ind >= 0 && ind < nInputs_) ddata = vecLScales[ind];
    }
  }
  return ddata;
}

