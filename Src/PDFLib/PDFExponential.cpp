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
// Functions for the Exponential distribution
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFExponential.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFExponential::PDFExponential(double lambda, double)
{
  lambda_ = lambda;
  if (lambda_ <= 0)
  {
    printf("PDFExponential: lambda_ should be > 0.\n");
    exit(1);
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFExponential::~PDFExponential()
{
}

// ************************************************************************
// forward transformation 
// ------------------------------------------------------------------------
int PDFExponential::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double xdata;

  if (psPDFDiagMode_ == 1)
    printf("PDFExponential: getPDF begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    xdata   = inData[ii];
    if (xdata < 0)
    {
      printf("PDFExponential getPDF ERROR: data not in [0,infty).\n");
      exit(1);
    }
    outData[ii] = lambda_ * exp(-lambda_ * xdata);
  }
  if (psPDFDiagMode_ == 1) printf("PDFExponential: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFExponential::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata;

  if (psPDFDiagMode_ == 1)
    printf("PDFExponential: getCDF begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata < 0)
    {
      printf("PDFExponential getCDF ERROR - data >= 0.\n");
      exit(1);
    }
    outData[ii] = 1.0 - exp(-lambda_ * ddata);
  }
  if (psPDFDiagMode_ == 1) printf("PDFExponential: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFExponential::invCDF(int length, double *inData, double *outData,
                           double lower, double upper)
{
  int    ii;
  double scale, ddata, xlo, xhi, xmi, ylo, yhi, ymi;

  if (upper <= lower)
  {
    printf("PDFExponential invCDF ERROR - lower bound >= upper bound.\n");
    exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFExponential invCDF ERROR - lower bound should be >= 0.\n");
    exit(1);
  }
  if (upper <= 0.0)
  {
    printf("PDFExponential invCDF ERROR - upper bound should be > 0.\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFExponential: invCDF begins (length = %d)\n",length);
  scale = upper - lower;
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii]; 
    if (ddata < 0)
    {
      printf("PDFExponential invCDF ERROR - data %e not in [0,infty)).\n",
             ddata);
      exit(1);
    }
    xlo = lower;
    ylo = 1.0 - exp(-lambda_ * xlo);
    xhi = upper;
    yhi = 1.0 - exp(-lambda_ * xhi);
    ddata = (ddata - lower) / scale * (yhi - ylo) + ylo;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = 1.0 - exp(-lambda_ * xmi);
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
  if (psPDFDiagMode_ == 1) printf("PDFExponential: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFExponential::genSample(int length, double *outData, double *lowers,
                              double *uppers)
{
  int    ii;
  double UU, xlo, xhi, xmi, ylo, yhi, ymi, lower, upper;

  if (psPDFDiagMode_ == 1)
    printf("PDFExponential: genSample begins (length = %d)\n",length);
  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFExp genSample ERROR - lower/upper bound not available.\n"); 
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  for (ii = 0; ii < length; ii++)
  {
    UU = PSUADE_drand();
    xlo = lower;
    ylo = 1.0 - exp(-lambda_ * xlo);
    xhi = upper;
    yhi = 1.0 - exp(-lambda_ * xhi);
    UU = UU * (yhi - ylo) + ylo;
    if      (UU <= ylo) outData[ii] = xlo;
    else if (UU >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(UU-ylo) > 1.0e-12 || PABS(UU-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = 1.0 - exp(-lambda_ * xmi);
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
  if (psPDFDiagMode_ == 1) printf("PDFExponential: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFExponential::getMean()
{
  if (lambda_ > 0.0) return 1.0 / lambda_;
  else               return 0.0;
}

