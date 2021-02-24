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
// Functions for the triangular distribution
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFTriangle.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFTriangle::PDFTriangle(double mean, double width)
{
  if (width <= 0.0)
  {
    printf("PDFTriangular ERROR: width = 0.0.\n");
    printf("     mode  = %e\n", mean);
    printf("     width = %e\n", width);
    exit(1);
  }
  mean_   = mean;
  width_  = width;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFTriangle::~PDFTriangle()
{
}

// ************************************************************************
// transformation 
// ------------------------------------------------------------------------
int PDFTriangle::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double height, ddata;

  if (psPDFDiagMode_ == 1)
    printf("PDFTriangle: getPDF begins (length = %d)\n",length);
  height = 1.0 / width_;
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata < (mean_-width_) || ddata > (mean_+width_)) 
      outData[ii] = 0.0;
    else if (ddata < mean_)
      outData[ii] = height - (mean_- ddata) * height / width_;
    else if (ddata >= mean_)
      outData[ii] = height - (ddata - mean_) * height / width_;
  }
  if (psPDFDiagMode_ == 1) printf("PDFTriangle: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFTriangle::getCDF(int length, double *inData, double *outData)
{
  int    ii, iOne=1;
  double ddata, ddata2;

  if (psPDFDiagMode_ == 1)
    printf("PDFTriangle: getCDF begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if      (ddata < (mean_ - width_)) outData[ii] = 0;
    else if (ddata > (mean_ + width_)) outData[ii] = 1;
    else if (ddata <= mean_)
    {
      getPDF(iOne, &ddata, &ddata2);
      outData[ii] = 0.5 * ddata2 * (ddata - mean_ + width_); 
    }
    else 
    {
      getPDF(iOne, &ddata, &ddata2);
      outData[ii] = 1.0 - 0.5 * ddata2 * (mean_ + width_ - ddata); 
    }
  }
  if (psPDFDiagMode_ == 1) printf("PDFTriangle: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation 
// ------------------------------------------------------------------------
int PDFTriangle::invCDF(int length, double *inData, double *outData,
                        double lower, double upper)
{
  int    ii, iOne=1;
  double range, low, ddata, scale, xlo, xhi, xmi, ylo, yhi, ymi;

  if (upper <= lower)
  {
    printf("PDFTriangle invCDF ERROR - lower bound >= upper bound.\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFTriangle: invCDF begins (length = %d)\n",length);
  getCDF(iOne, &lower, &low);
  getCDF(iOne, &upper, &range);
  range = range - low;
  scale = upper - lower;
  for (ii = 0; ii < length; ii++)
  {
    xlo = lower;
    getCDF(iOne, &xlo, &ylo);
    xhi = upper;
    getCDF(iOne, &xhi, &yhi);
    ddata = (inData[ii] - lower) / scale * range + low;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        getCDF(iOne, &xmi, &ymi);
        if (ddata > ymi)
        {
          xlo = xmi;
          ylo = ymi;
        }
        else
        {
          xhi = xmi;
          yhi = ymi;
        }
      }
      if (PABS(ddata-ylo) < PABS(ddata-yhi)) outData[ii] = xlo;
      else                                   outData[ii] = xhi;
    }
  }
   if (psPDFDiagMode_ == 1) printf("PDFTriangle: invCDF ends.\n");
   return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFTriangle::genSample(int length, double *outData, double *lowers,
                           double *uppers)
{
  int    ii, iOne=1;
  double range, low, xlo, xhi, xmi, ylo, yhi, ymi, UU, lower, upper;

  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFTriangle genSample ERROR - lower/upper bound unavailable.\n"); 
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (upper <= lower)
  {
    printf("PDFTriangle genSample ERROR - lower bound >= upper bound.\n");
    exit(1);
  }

  getCDF(iOne, &lower, &low);
  getCDF(iOne, &upper, &range);
  range = range - low;
  if (psPDFDiagMode_ == 1)
    printf("PDFTriangle: genSample begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    xlo = lower;
    getCDF(iOne, &xlo, &ylo);
    xhi = upper;
    getCDF(iOne, &xhi, &yhi);
    UU = PSUADE_drand() * range + low;
    if      (UU <= ylo) outData[ii] = xlo;
    else if (UU >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(UU-ylo) > 1.0e-12 || PABS(UU-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        getCDF(iOne, &xmi, &ymi);
        if (UU > ymi)
        {
          xlo = xmi;
          ylo = ymi;
        }
        else
        {
          xhi = xmi;
          yhi = ymi;
        }
      }
      if (PABS(UU-ylo) < PABS(UU-yhi)) outData[ii] = xlo;
      else                             outData[ii] = xhi;
    }
  }
  if (psPDFDiagMode_ == 1) printf("PDFTriangle: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFTriangle::getMean()
{
  return mean_;
}

