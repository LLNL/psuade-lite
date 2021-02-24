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
// Functions for the Weibull distribution
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFWeibull.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFWeibull::PDFWeibull(double lambda, double k)
{
  lambda_ = lambda;
  k_      = k;
  if (lambda_ < 0 || k_ < 0)
  {
    printf("PDFWeibull: lambda and k have to be > 0.\n");
    exit(1);
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFWeibull::~PDFWeibull()
{
}

// ************************************************************************
// forward
// ------------------------------------------------------------------------
int PDFWeibull::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double xdata, ddata, mult;

  if (psPDFDiagMode_ == 1)
    printf("PDFWeibull: getPDF begins (length = %d)\n",length);
  mult = k_ / lambda_;
  for (ii = 0; ii < length; ii++)
  {
    xdata = inData[ii];
    if (xdata < 0) outData[ii] = 0.0;
    else
    {
      ddata = xdata / lambda_;
      outData[ii] = mult * pow(ddata,k_-1) * exp(-pow(ddata,k_));
    }
  }
  if (psPDFDiagMode_ == 1) printf("PDFWeibull: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFWeibull::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata;

  if (psPDFDiagMode_ == 1)
    printf("PDFWeibull: getCDF begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if      (ddata < 0)
    {
      printf("PDFWeibull getCDF ERROR: incoming data < 0.\n");
      exit(1);
    }
    outData[ii] = 1.0 - exp(-pow(ddata/lambda_,k_));
  }
  if (psPDFDiagMode_ == 1) printf("PDFWeibull: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFWeibull::invCDF(int length, double *inData, double *outData,
                    double lower, double upper)
{
  int    ii;
  double scale, ddata, xlo, xhi, ylo, yhi, xmi, ymi;

  if (upper <= lower)
  {
    printf("PDFWeibull invCDF ERROR - lower bound >= upper bound.\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFWeibull: invCDF begins (length = %d)\n",length);
  scale = upper - lower;
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata < 0)
    {
      printf("PDFWeibull invCDF ERROR - data %e not in [0,infty).\n",
             ddata);
      exit(1);
    }
    xlo = lower;
    xhi = upper;
    ylo = 1.0 - exp(-pow(xlo/lambda_,k_));
    yhi = 1.0 - exp(-pow(xhi/lambda_,k_));
    ddata = (ddata - lower) / scale * (yhi - ylo) + ylo;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = 1.0 - exp(-pow(xmi/lambda_,k_));
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
  if (psPDFDiagMode_ == 1) printf("PDFWeibull: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFWeibull::genSample(int length, double *outData, double *lowers,
                          double *uppers)
{
  int    ii;
  double UU, xlo, xhi, xmi, ylo, yhi, ymi, lower, upper;

  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFWeibull genSample ERROR - lower/upper bound unavailable.\n");
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (upper <= lower)
  {
    printf("PDFWeibull genSample ERROR - lower bound >= upper bound.\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFWeibull: genSample begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    UU = PSUADE_drand();
    xlo = lower;
    xhi = upper;
    ylo = 1.0 - exp(-pow(xlo/lambda_,k_));
    yhi = 1.0 - exp(-pow(xhi/lambda_,k_));
    UU = UU * (yhi - ylo) + ylo;
    if      (UU <= ylo) outData[ii] = xlo;
    else if (UU >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(UU-ylo) > 1.0e-12 || PABS(UU-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = 1.0 - exp(-pow(xmi/lambda_,k_));
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
  if (psPDFDiagMode_ == 1) printf("PDFWeibull: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFWeibull::getMean()
{
  if (k_ > 0.0) return (lambda_ * Gamma_Function(1.0 + 1.0 / k_));
  else          return 0.0;
}

