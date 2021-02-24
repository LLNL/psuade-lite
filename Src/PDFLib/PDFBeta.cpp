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
// Functions for the beta distribution
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#include <stdio.h>
#include <math.h>
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFBeta.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PDFBeta::PDFBeta(double alpha, double beta)
{
  alpha_  = alpha;
  beta_   = beta;
  if (alpha < 0 || beta < 0)
  {
    printf("PDFBeta: currently only supporting alpha, beta > 0.\n");
    exit(1);
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFBeta::~PDFBeta()
{
}

// ************************************************************************
// forward transformation to range
// ------------------------------------------------------------------------
int PDFBeta::getPDF(int length, double *inData, double *outData)
{
  int    ii;
  double mult, xdata;

  if (psPDFDiagMode_ == 1)
    printf("PDFBeta: getPDF begins (length = %d)\n",length);
  mult = 1.0 / Beta_Function(alpha_,beta_);
  for (ii = 0; ii < length; ii++)
  {
    xdata = inData[ii];
    if (xdata <= 0 || xdata >= 1)
    {
      printf("PDFBeta getCDF ERROR: data not in (0, 1).\n");
      exit(1);
    }
    outData[ii] = pow(xdata, alpha_-1.0) * pow(1.0-xdata, beta_-1.0);
    outData[ii] *= mult;
  }
  if (psPDFDiagMode_ == 1) printf("PDFBeta: getPDF ends.\n");
  return 0;
}

// ************************************************************************
// look up cumulative density
// ------------------------------------------------------------------------
int PDFBeta::getCDF(int length, double *inData, double *outData)
{
  int    ii;
  double ddata, mult;

  if (psPDFDiagMode_ == 1)
    printf("PDFBeta: getCDF begins (length = %d)\n",length);
  mult = 1.0 / Beta_Function(alpha_,beta_);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if      (ddata <= 0) outData[ii] = 0;
    else if (ddata >= 1) outData[ii] = 1;
    else
      outData[ii] = mult*Incomplete_Beta_Function(ddata,alpha_,beta_);
  }
  if (psPDFDiagMode_ == 1) printf("PDFBeta: getCDF ends.\n");
  return 0;
}

// ************************************************************************
// transformation to range
// ------------------------------------------------------------------------
int PDFBeta::invCDF(int length, double *inData, double *outData,
                    double lower, double upper)
{
  int    ii;
  double scale, ddata, xlo, xhi, xmi, ylo, ymi, yhi, mult;

  if (upper <= lower)
  {
    printf("PDFBeta invCDF ERROR - lower bound >= upper bound.\n");
    exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFBeta invCDF ERROR - lower bound < 0.\n");
    printf("        INFO : support for beta distribution is (0,1).\n");
    exit(1);
  }
  if (upper > 1.0)
  {
    printf("PDFBeta invCDF ERROR - upper bound > 1.\n");
    printf("        INFO : support for beta distribution is (0,1).\n");
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
    printf("PDFBeta: invCDF begins (length = %d)\n",length);
  scale = upper - lower;
  mult = 1.0 / Beta_Function(alpha_,beta_);
  for (ii = 0; ii < length; ii++)
  {
    ddata = inData[ii];
    if (ddata < 0 || ddata > 1)
    {
      printf("PDFBeta invCDF ERROR - data %e not in (0,1).\n",ddata);
      printf("        INFO : support for beta distribution is (0,1).\n");
      exit(1);
    }
    xlo = lower;
    ylo = mult*Incomplete_Beta_Function(xlo,alpha_,beta_);
    xhi = upper;
    yhi = mult*Incomplete_Beta_Function(xhi,alpha_,beta_);
    ddata = (ddata - lower) / scale * (yhi - ylo) + ylo;
    if      (ddata <= ylo) outData[ii] = xlo;
    else if (ddata >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(ddata-ylo) > 1.0e-12 || PABS(ddata-yhi) > 1.0e-12)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = mult*Incomplete_Beta_Function(xmi,alpha_,beta_);
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
  if (psPDFDiagMode_ == 1) printf("PDFBeta: invCDF ends.\n");
  return 0;
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFBeta::genSample(int length,double *outData,double *lowers,
                       double *uppers)
{
  int    ii;
  double UU, xhi, xlo, xmi, yhi, ylo, ymi, mult, lower, upper;

  if (lowers == NULL || uppers == NULL)
  {
    printf("PDFBeta genSample ERROR - lower/upper bound not available.\n");
    exit(1);
  }
  lower = lowers[0];
  upper = uppers[0];
  if (upper <= lower)
  {
    printf("PDFBeta genSample ERROR - lower bound >= upper bound.\n");
    exit(1);
  }
  if (lower < 0.0)
  {
    printf("PDFBeta genSample  ERROR - lower bound < 0.\n");
    printf("        INFO : support for beta distribution is (0,1).\n");
    exit(1);
  }
  if (upper > 1.0)
  {
    printf("PDFBeta genSample  ERROR - upper bound > 1.\n");
    printf("        INFO : support for beta distribution is (0,1).\n");
    exit(1);
  }

  mult = 1.0 / Beta_Function(alpha_,beta_);
  if (psPDFDiagMode_ == 1)
    printf("PDFBeta: genSample begins (length = %d)\n",length);
  for (ii = 0; ii < length; ii++)
  {
    UU = PSUADE_drand();
    xlo = lower;
    ylo = mult*Incomplete_Beta_Function(xlo,alpha_,beta_);
    xhi = upper;
    yhi = mult*Incomplete_Beta_Function(xhi,alpha_,beta_);
    UU = UU * (yhi - ylo) + ylo;
    if      (UU <= ylo) outData[ii] = xlo;
    else if (UU >= yhi) outData[ii] = xhi;
    else
    {
      while (PABS(UU-ylo) > 1.0e-10 || PABS(UU-yhi) > 1.0e-10)
      {
        xmi = 0.5 * (xhi + xlo);
        ymi = mult*Incomplete_Beta_Function(xmi,alpha_,beta_);
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
  if (psPDFDiagMode_ == 1) printf("PDFBeta: genSample ends.\n");
  return 0;
}

// ************************************************************************
// get mean
// ------------------------------------------------------------------------
double PDFBeta::getMean()
{
  if (alpha_ + beta_ > 0.0) return alpha_ / (alpha_ + beta_);
  else                      return 0.0;
}

