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
// Functions for the Snedecor F distribution
// AUTHOR : CHARLES TONG
// DATE   : 2015
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFF.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFF::PDFF(double d1, double d2)
{
  d1_ = d1;
  d2_ = d2;
  if (d1 <= 0 || d2 <= 0)
  {
    printf("PDFF: d1 and d2 need to be > 0.\n");
    exit(1);
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFF::~PDFF()
{
}

// ************************************************************************
// forward transformation 
// ------------------------------------------------------------------------
int PDFF::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata, xdata, d2d2, d1d2, beta;

  if (psPDFDiagMode_ == 1)
    printf("PDFF: getPDF begins (length = %d)\n",length);
  d2d2 = pow(d2_, d2_);
  d1d2 = d1_ + d2_;
  beta = Beta_Function(0.5*d1_,0.5*d2_);
  for (ii = 0; ii < length; ii++)
  {
    xdata = inData[ii];
    if (xdata < 0.0)
    {
      printf("PDFF getPDF ERROR - data needs to be in [0,infty).\n");
      exit(1);
    }
    ddata = pow(d1_ * xdata, d1_) * d2d2;
    ddata = sqrt(ddata / pow(d1_ * xdata + d2_, d1d2));
    outData[ii] = ddata / (xdata * beta);
  }
  if (psPDFDiagMode_ == 1) printf("PDFF: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFF::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double xdata, mult;

  if (psPDFDiagMode_ == 1)
    printf("PDFF: getCDF begins (length = %d)\n",length);
  mult = 1.0 / Beta_Function(0.5*d1_,0.5*d2_);
  for (ii = 0; ii < length; ii++)
  {
    xdata = inData[ii];
    if      (xdata <= 0) outData[ii] = 0;
    else
    {
      xdata = d1_ * xdata / (d1_ * xdata + d2_);
      outData[ii] = mult*Incomplete_Beta_Function(xdata,0.5*d1_,0.5*d2_);
    }
  }
  if (psPDFDiagMode_ == 1) printf("PDFF: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFF::invCDF(int length, double *inData, double *outData,
                 double lower, double upper)
{
  int    ii;
  double scale, ddata, xlo, xhi, xmi, ylo, ymi, yhi, mult;

  if (upper <= lower)
  {
    printf("PDFF invCDF ERROR - lower bound >= upper bound.\n");
    exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFF invCDF ERROR - lower bound < 0.\n");
    printf("     INFO : support for F distribution is (0,1).\n");
    exit(1);
  }
  if (upper > 1.0)
  {
    printf("PDFF invCDF ERROR - upper bound > 1.\n");
    printf("     INFO : support for F distribution is (0,1).\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFF: invCDF begins (length = %d)\n",length);
  scale = upper - lower;
  mult = 1.0 / Beta_Function(0.5*d1_,0.5*d2_);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata < 0 || ddata > 1)
    {
      printf("PDFF invCDF ERROR - data %e not in (0,1).\n",ddata);
      printf("     INFO : support for F distribution is (0,1).\n");
      exit(1);
    }
    xlo = lower * d1_ / (d1_ * lower + d2_);
    ylo = mult*Incomplete_Beta_Function(xlo,0.5*d1_,0.5*d2_);
    xhi = upper * d1_ / (d1_ * upper + d2_);
    yhi = mult*Incomplete_Beta_Function(xhi,0.5*d1_,0.5*d2_);
    ddata = (ddata - lower) / scale * (yhi - ylo) + ylo;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        xmi = xmi * d1_ / (d1_ * xmi + d2_);
        ymi = mult*Incomplete_Beta_Function(xmi,0.5*d1_,0.5*d2_);
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
  if (psPDFDiagMode_ == 1) printf("PDFF: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFF::genSample(int length, double *outData, double *lowers, 
                    double *uppers)
{
  int    ii;
  double UU, xhi, xlo, xmi, yhi, ylo, ymi, mult, xlo2, xhi2, xmi2;
  double lower, upper;

  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFF genSample ERROR - lower/upper bound not available.\n"); 
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (upper <= lower)
  {
    printf("PDFF genSample ERROR - lower bound >= upper bound.\n");
    exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFF genSample  ERROR - lower bound < 0.\n");
    printf("     INFO : support for F distribution is (0,1).\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFF: genSample begins (length = %d)\n",length);
  mult = 1.0 / Beta_Function(0.5*d1_,0.5*d2_);
  for (ii = 0; ii < length; ii++)
  {
    UU   = PSUADE_drand();
    xlo  = lower;
    xlo2 = xlo * d1_ / (xlo * d1_ + d2_);
    ylo  = mult*Incomplete_Beta_Function(xlo2,0.5*d1_,0.5*d2_);
    xhi  = upper;
    xhi2 = xhi * d1_ / (xhi * d1_ + d2_);
    yhi  = mult*Incomplete_Beta_Function(xhi2,0.5*d1_,0.5*d2_);
    UU = UU * (yhi - ylo) + ylo;
    if      (UU <= ylo) outData[ii] = xlo;
    else if (UU >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(UU-ylo) > 1.0e-12 || PABS(UU-yhi) > 1.0e-12)
      {
        xmi  = 0.5 * (xhi + xlo);
        xmi2 = xmi * d1_ / (xmi * d1_ + d2_);
        ymi = mult*Incomplete_Beta_Function(xmi2,0.5*d1_,0.5*d2_);
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
  if (psPDFDiagMode_ == 1) printf("PDFF: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFF::getMean()
{
  if (d2_ > 2) return d2_ / (d2_ - 2);
  else
  {
    printf("PDFF INFO: getMean - d2 <= 2, return 0.\n");
    return 0.0;
  }
}

