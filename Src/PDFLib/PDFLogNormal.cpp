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
// Functions for the lognormal distribution
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFLogNormal.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFLogNormal::PDFLogNormal(double lmean, double lstdev)
{
  if (lstdev < 0.0) 
  {
    printf("PDFLogNormal ERROR: stdev < 0.0.\n");
    exit(1);
  }
  stdev_ = lstdev;
  mean_  = lmean;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFLogNormal::~PDFLogNormal()
{
}

// ************************************************************************
// forward transformation to range
// ------------------------------------------------------------------------
int PDFLogNormal::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double denom, coef, xdata, expo;

  if (psPDFDiagMode_ == 1)
    printf("PDFLogNormal: getPDF begins (length = %d)\n", length);
  denom = 2.0 * stdev_ * stdev_;
  coef  = 1.0 / (stdev_ * sqrt(2.0*M_PI));
  for (ii = 0; ii < length; ii++)
  {
    xdata   = inData[ii];
    if (xdata <= 0.0)
    {
      printf("PDFLogNormal getPDF ERROR - data <= 0 (for log)\n");
      exit(1);
    }
    expo = - (log(xdata) - mean_) * (log(xdata) - mean_) / denom;
    outData[ii] = coef * exp(expo) / xdata;
  }
  if (psPDFDiagMode_ == 1) printf("PDFLogNormal: getPDF ends\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFLogNormal::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata, iroot2;

  if (psPDFDiagMode_ == 1)
    printf("PDFLogNormal: getCDF begins (length = %d)\n", length);
  iroot2 = sqrt(0.5)/stdev_;
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if      (ddata == 0.0) outData[ii] = 0;
    else if (ddata < 0.0)
    {
      printf("PDFLogNormal getCDF ERROR - data < 0 (for log)\n");
      exit(1);
    }
    else outData[ii] = 0.5 * (1.0 + erf((log(ddata)-mean_)*iroot2));
  }
  if (psPDFDiagMode_ == 1) printf("PDFLogNormal: getCDF ends\n");
  return 0;
}

// ************************************************************************
// look up cumulative distribution
// ------------------------------------------------------------------------
int PDFLogNormal::invCDF(int length, double *inData, double *outData,
                         double lower, double upper)
{
  int    ii;
  double ddata, scale, xlo, xhi, ylo, yhi, xmi, ymi, iroot2;

  if (lower < 0.0)
  {
    printf("PDFLogNormal invCDF ERROR - lower bound < 0.\n");
    printf("  Lower bound = %e\n", lower);
    exit(1);
  }
  if (upper <= lower)
  {
    printf("PDFLogNormal invCDF ERROR - lower bound >= upper bound.\n");
    printf("  Lower bound = %e\n", lower);
    printf("  Upper bound = %e\n", upper);
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFLogNormal: invCDF begins (length = %d)\n", length);
  scale = upper - lower;
  iroot2 = sqrt(0.5)/stdev_;
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata < 0.0)
    {
      printf("PDFLogNormal invCDF ERROR - data < 0.\n");
      exit(1);
    }
    xlo = lower;
    if (xlo == 0) ylo = 0.0;
    else          ylo = 0.5 * (1.0 + erf((log(xlo)-mean_)*iroot2));
    xhi = upper;
    yhi = 0.5 * (1.0 + erf((log(xhi)-mean_)*iroot2));
    ddata = (ddata - lower) / scale * (yhi - ylo) + ylo;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = 0.5 * (1.0 + erf((log(xmi)-mean_)*iroot2));
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
  if (psPDFDiagMode_ == 1) printf("PDFLogNormal: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFLogNormal::genSample(int length, double *outData, double *lowers,
                            double *uppers)
{
  int    ii, count, total;
  double U1, U2, R, theta, Z1, Z2, pi=3.14159, range, low, iroot2;
  double lower, upper, lower2, upper2;

  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFLogNormal genSample ERROR - lower/upper bound unavailable.\n");
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (length <= 0)
  {
    printf("PDFLogNormal genSample ERROR - length <= 0.\n");
    exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFLogNormal genSample ERROR - lower bound < 0.\n");
    printf("  Lower bound = %e\n", lower);
    exit(1);
  }
  if (upper <= lower)
  {
    printf("PDFLogNormal genSample ERROR - lower bound >= upper bound.\n");
    printf("  Lower bound = %e\n", lower);
    printf("  Upper bound = %e\n", upper);
    exit(1);
  }

  if (stdev_ == 0)
  {
    printf("PDFLogNormal: genSample WARNING - std. dev. = 0\n");
    for (ii = 0; ii < length; ii++) outData[ii] = exp(mean_);
    return 0;
  }

  iroot2 = sqrt(0.5)/stdev_;
  lower2 = exp(mean_-5*stdev_);
  upper2 = exp(mean_+5*stdev_);
  if (lower2 == 0) low = 0.0;
  else             low = 0.5 * (1.0 + erf((log(lower2)-mean_)*iroot2));
  range = 0.5 * (1.0 + erf((log(upper2)-mean_)*iroot2)) - low;
  count = total = 0;
  if (psPDFDiagMode_ == 1) 
    printf("PDFLogNormal: genSample begins (length = %d)\n",length);
  while (count < length)
  {
    U1 = PSUADE_drand() * range + low;
    U2 = PSUADE_drand() * range + low;
    R  = sqrt(-2.0 * log(U1));
    theta = 2 * pi * U2;
    Z1 = R * cos(theta);
    Z2 = R * sin(theta);
    outData[count] = exp(mean_ + stdev_ * Z1);
    if (outData[count] >= lower && outData[count] <= upper) count++;
    total++;
    if (count < length)
    {
      outData[count] = exp(mean_ + stdev_ * Z2);
      if (outData[count] >= lower && outData[count] <= upper) count++;
      total++;
    }
    if (total > length*1000)
    {
      printf("PDFLogNormal genSample ERROR - cannot generate enough.\n");
      printf("             sample points to be within range. Maybe\n");
      printf("             due to prescribed ranges too narro.\n");
      printf("      mean, stdev  = %e %e\n", mean_, stdev_);
      printf("      lower, upper = %e %e\n", lower, upper);
      exit(1);
    }
  }
  if (psPDFDiagMode_ == 1) 
  {
    printf("PDFLogNormal: genSample ends.\n");
    if (total > length)
      printf("PDFLogNormal Statistics: need %d to generate %d points.\n",
             total,length);
  }
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFLogNormal::getMean()
{
  return exp(mean_ + 0.5 * stdev_ * stdev_);
}

