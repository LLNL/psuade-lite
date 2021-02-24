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
// Functions for the normal distribution
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFNormal.h"
#define PABS(x) ((x >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFNormal::PDFNormal(double mean, double stdev)
{
  mean_   = mean;
  stdev_  = stdev;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PDFNormal::~PDFNormal()
{
}

// ************************************************************************
// forward transformation to range
// ------------------------------------------------------------------------
int PDFNormal::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double denom, coef, xdata, expo;

  if (psPDFDiagMode_ == 1) 
    printf("PDFNormal: getPDF begins (length = %d)\n",length);
  denom = 2.0 * stdev_ * stdev_;
  coef  = 1.0 / (stdev_ * sqrt(2.0*M_PI));
  for (ii = 0; ii < length; ii++)
  {
    xdata   = inData[ii];
    expo = - (xdata - mean_) * (xdata - mean_) / denom;
    outData[ii] = coef * exp(expo);
  }
  if (psPDFDiagMode_ == 1) printf("PDFNormal: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density 
// ------------------------------------------------------------------------
int PDFNormal::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata, iroot2;

  if (psPDFDiagMode_ == 1) 
    printf("PDFNormal: getCDF begins (length = %d)\n", length);
  iroot2 = sqrt(0.5)/stdev_;
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    outData[ii] = 0.5 * (1.0 + erf((ddata-mean_)*iroot2));
  }
  if (psPDFDiagMode_ == 1) printf("PDFNormal: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density in reverse
// ------------------------------------------------------------------------
int PDFNormal::invCDF(int length, double *inData, double *outData,
                      double lower, double upper)
{
  int    ii;
  double scale, ddata, iroot2, xlo, ylo, xhi, yhi, xmi, ymi;

  if (upper <= lower)
  {
    printf("PDFNormal invCDF ERROR - lower bound >= upper bound.\n");
    printf("                         lower bound = %e\n",lower);
    printf("                         upper bound = %e\n",upper);
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFNormal: invCDF begins (length = %d)\n",length);
  scale = upper - lower;
  iroot2 = sqrt(0.5)/stdev_;
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    xlo = lower;
    xhi = upper;
    ylo = 0.5 * (1.0 + erf((xlo-mean_)*iroot2));
    yhi = 0.5 * (1.0 + erf((xhi-mean_)*iroot2));
    ddata = (ddata - lower) / scale  * (yhi - ylo) + ylo;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = 0.5 * (1.0 + erf((xmi-mean_)*iroot2));
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
  if (psPDFDiagMode_ == 1) printf("PDFNormal: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFNormal::genSample(int length, double *outData, double *lowers,
                         double *uppers)
{
  int    ii, count, total;
  double U1, U2, R, theta, Z1, Z2, pi=3.14159, range, low, iroot2;
  double lower, upper, lower2, upper2;

  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFNormal genSample ERROR - lower/upper bound unavailable.\n"); 
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (length <= 0)
  {
    printf("PDFNormal genSample ERROR - length <= 0.\n");
    exit(1);
  }
  if (upper <= lower)
  {
    printf("PDFNormal genSample ERROR - lower bound >= upper bound.\n");
    printf("          lower = %24.16e\n", lower);
    printf("          upper = %24.16e\n", upper);
    exit(1);
  }

  if (stdev_ == 0)
  {
    printf("PDFNormal WARNING: genSample - std dev = 0.\n");
    for (ii = 0; ii < length; ii++) outData[ii] = mean_;
    return 0;
  }

  lower2 = mean_ - 5 * stdev_;
  upper2 = mean_ + 5 * stdev_;
  iroot2 = sqrt(0.5)/stdev_;
  low   = 0.5 * (1.0 + erf((lower2-mean_)*iroot2));
  range = 0.5 * (1.0 + erf((upper2-mean_)*iroot2)) - low;
  count = total = 0;
  if (psPDFDiagMode_ == 1)
    printf("PDFNormal: genSample begins (length = %d)\n",length);
  while (count < length)
  {
    U1 = PSUADE_drand() * range + low;
    U2 = PSUADE_drand() * range + low;
    R  = sqrt(-2.0 * log(U1));
    theta = 2 * pi * U2;
    Z1 = R * cos(theta);
    Z2 = R * sin(theta);
    outData[count] = mean_ + stdev_ * Z1;
    if (outData[count] >= lower && outData[count] <= upper) count++;
    if (count >= length) break;
    outData[count] = mean_ + stdev_ * Z2;
    if (outData[count] >= lower && outData[count] <= upper) count++;
    total += 2;
    if (total > length*100)
    {
      printf("PDFNormal genSample ERROR - Cannot generate enough\n");
      printf("          sample points to be within range. Maybe\n");
      printf("          due to prescribed ranges too narrow.\n");
      printf("     mean,  stdev = %e %e\n", mean_, stdev_);
      printf("     lower, upper = %e %e\n", lower, upper);
      printf("     ntrials, nsuccess = %d %d (%d)\n",total,count,length);
      exit(1);
    }
  }
  if (psPDFDiagMode_ == 1) printf("PDFNormal: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFNormal::getMean()
{
  return mean_;
}

