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
// Functions for the inverse gamma distribution
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************

#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFInvGamma.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFInvGamma::PDFInvGamma(double alpha, double beta)
{
  alpha_ = alpha;
  beta_  = beta;
  if (alpha <= 0 || beta <= 0)
  {
    printf("PDFInvGamma: alpha and beta need to be > 0.\n");
    exit(1);
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFInvGamma::~PDFInvGamma()
{
}

// ************************************************************************
// forward transformation 
// ------------------------------------------------------------------------
int PDFInvGamma::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double mult, xdata;

  if (psPDFDiagMode_ == 1)
    printf("PDFInvGamma: getPDF begins (length = %d)\n",length);
  mult = pow(beta_, alpha_) / Gamma_Function(alpha_);
  for (ii = 0; ii < length; ii++)
  {
    xdata = inData[ii];
    if (xdata <= 0.0)
    {
      printf("PDFInvGamma getPDF ERROR - data needs to be in [0,infty).\n");
      exit(1);
    }
    outData[ii] = pow(xdata, -alpha_-1.0) * exp(-beta_ / xdata);
    outData[ii] *= mult;
  }
  if (psPDFDiagMode_ == 1) printf("PDFInvGamma: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFInvGamma::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata, mult;

  if (psPDFDiagMode_ == 1)
    printf("PDFInvGamma: getCDF begins (length = %d)\n",length);
  mult = 1.0 / Gamma_Function(alpha_);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if   (ddata <= 0) outData[ii] = 0;
    else outData[ii] = mult*Incomplete_Gamma_Function(beta_/ddata,alpha_);
  }
  if (psPDFDiagMode_ == 1) printf("PDFInvGamma: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFInvGamma::invCDF(int length, double *inData, double *outData,
                    double lower, double upper)
{
  int    ii;
  double scale, ddata, mult, xlo, xmi, xhi, ylo, ymi, yhi;

  if (upper <= lower)
  {
    printf("PDFInvGamma invCDF ERROR - lower bound >= upper bound.\n");
    exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFInvGamma invCDF ERROR - lower bound <= 0.\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFInvGamma: invCDF begins (length = %d)\n",length);
  scale = upper - lower;
  mult = 1.0 / Gamma_Function(alpha_);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata < 0.0)
    {
      printf("PDFInvGamma invCDF ERROR - data %e not in [0,infty).\n",
             ddata);
      exit(1);
    }
    xlo = lower;
    ylo = mult * Incomplete_Gamma_Function(beta_/xlo,alpha_);
    xhi = upper;
    yhi = mult * Incomplete_Gamma_Function(beta_/xhi,alpha_);
    ddata = (ddata - lower) / scale * (yhi - ylo) + ylo;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = mult * Incomplete_Gamma_Function(beta_/xmi,alpha_);
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
  if (psPDFDiagMode_ == 1) printf("PDFInvGamma: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFInvGamma::genSample(int length, double *outData, double *lowers,
                        double *uppers)
{
  int    ii;
  double UU, xhi, xlo, xmi, yhi, ylo, ymi, mult, lower, upper;

  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFInvGamma genSample ERROR - lower/upper bound unavailable.\n"); 
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (upper <= lower)
  {
     printf("PDFInvGamma genSample ERROR - lower bound >= upper bound.\n");
     exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFInvGamma genSample ERROR - lower bound < 0.\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFInvGamma: genSample begins (length = %d)\n",length);
  mult = 1.0 / Gamma_Function(alpha_);
  for (ii = 0; ii < length; ii++)
  {
    if ((length > 10000) && (ii % 1000 == 0))
      printf("InvGamma genSample: %d out of %d generated\n",ii,length);
    UU = PSUADE_drand();
    xlo = lower;
    ylo = mult * Incomplete_Gamma_Function(beta_/xlo,alpha_);
    xhi = upper;
    yhi = mult * Incomplete_Gamma_Function(beta_/xhi,alpha_);
    UU = UU * (yhi - ylo) + ylo;
    if      (UU <= ylo) outData[ii] = xlo;
    else if (UU >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(UU-ylo) > 1.0e-12 || PABS(UU-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = mult * Incomplete_Gamma_Function(beta_/xmi,alpha_);
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
  if (psPDFDiagMode_ == 1) printf("PDFInvGamma: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFInvGamma::getMean()
{
  if (alpha_ <= 1)
  {
    printf("PDFInvGamma ERROR: mean does not exist for alpha <= 1.\n");
    return 0.0;
  }
  return beta_ / (alpha_ - 1);
}

