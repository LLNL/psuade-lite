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
// Functions for the multivariate log normal distribution
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "sysdef.h"
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PDFNormal.h"
#include "PDFMVLogNormal.h"
#define PABS(x) (((x) >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
PDFMVLogNormal::PDFMVLogNormal(psVector &means, psMatrix &corMat)
{
  int ii, jj, length;
  double ddata, sigmaI, sigmaJ;

  if (means.length() <= 0)
  {
    printf("PDFMVLogNormal ERROR: number of inputs <= 0.\n");
    exit(1);
  }
  if (means.length() != corMat.ncols())
  {
    printf("PDFMVLogNormal ERROR: mismatch mean and cov dimensions.\n");
    printf("                      mean dimension =  %d\n",means.length());
    printf("                      cov  dimension =  %d\n",corMat.ncols());
    exit(1);
  }
  length = means.length();
  for (ii = 0; ii < length; ii++)
  {
    ddata = means[ii];
    if (ddata < 0)
    {
      printf("PDFMVLogNormal ERROR: mean should be >= 0.\n");
      exit(1);
    }
    ddata = corMat.getEntry(ii,ii);
    if (ddata <= 0)
    {
      printf("PDFMVLogNormal ERROR: stdev should be > 0.\n");
      exit(1);
    }
  }

  means_.setLength(length);
  covMat_.setDim(length,length);
  for (ii = 0; ii < length; ii++)
  {
    means_[ii] = means[ii];
    sigmaI = corMat.getEntry(ii,ii);
    for (jj = 0; jj < length; jj++)
    {
      sigmaJ = corMat.getEntry(jj,jj);
      if (ii == jj) ddata = 1.0;
      else          ddata = corMat.getEntry(ii,jj);
      ddata *= sigmaI*sigmaJ;
      covMat_.setEntry(ii,jj,ddata);
    }
  }
  covMat_.CholDecompose();
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PDFMVLogNormal::~PDFMVLogNormal()
{
}

// ************************************************************************
// generate a sample
// ------------------------------------------------------------------------
int PDFMVLogNormal::genSample(int length, psVector &vecOut, psVector &lower,
                              psVector &upper)
{
  int    ii, jj, nInputs, flag, count, total;
  double One=1.0, Zero=0.0, d1, d2;
  psVector vecData1, vecData2, vecTemp1, vecTemp2;
  PDFNormal *normalPtr;

  if (length <= 0)
  {
    printf("PDFMVLogNormal genSample ERROR - vec length <= 0.\n");
    exit(1);
  }
  nInputs = means_.length();
  if (vecOut.length()/length != nInputs)
  {
    printf("PDFMVLogNormal genSample ERROR - length mismatch.\n");
    printf("               lengths    = %d x %d.\n", length, nInputs);
    printf("               vec length = %d.\n", vecOut.length());
    exit(1);
  }

  if (psPDFDiagMode_ == 1)
     printf("PDFMVLogNormal: genSample begins (length = %d)\n",length);
  normalPtr  = new PDFNormal(Zero, One);
  vecTemp1.setLength(length);
  vecTemp2.setLength(length*nInputs);
  double *tmpData1 = vecTemp1.getDVector();
  double *tmpData2 = vecTemp2.getDVector();
  count = total = 0;
  d1 = -4.0 * One;
  d2 =  4.0 * One;
  while (count < length)
  {
    for (jj = 0; jj < nInputs; jj++)
    {
      normalPtr->genSample(length,tmpData1,&d1,&d2);
      for (ii = 0; ii < length; ii++)
        tmpData2[ii*nInputs+jj] = tmpData1[ii];
    }
    for (ii = 0; ii < length; ii++)
    {
      vecData1.load(nInputs, &(tmpData2[ii*nInputs]));
      covMat_.CholLMatvec(vecData1, vecData2);
      for (jj = 0; jj < nInputs; jj++)
        vecOut[count*nInputs+jj] = exp(means_[jj] + vecData2[jj]);
      flag = 0;
      for (jj = 0; jj < nInputs; jj++)
      {
        if (vecOut[count*nInputs+jj] < lower[jj] ||
           vecOut[count*nInputs+jj] > upper[jj]) {flag = 1; break;}
      }
      if (flag == 0) count++;
      total++;
      if (count >= length) break;
    }
    if (total > length*1000)
    {
      printf("PDFLogMVNormal ERROR : cannot generate enough sample points\n");
      printf("                       within the given ranges. Check your \n");
      printf("                       parameter lower and upper bounds.\n");
      for (jj = 0; jj < nInputs; jj++)
        printf("Input bounds = %e %e \n", lower[jj], upper[jj]);
      exit(1);
    }
  }
  if (psPDFDiagMode_ == 1)
  {
    printf("PDFMVLogNormal: genSample ends.\n");
    if (total > length)
      printf("PDFLogNormal Statistics: need %d to generate %d points.\n",
             total,length);
  }
  delete normalPtr;
  return 0;
}

