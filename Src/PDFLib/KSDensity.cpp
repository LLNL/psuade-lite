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
// Functions for the KSDensity
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "KSDensity.h"
#define PABS(x) ((x >= 0) ? x : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
KSDensity::KSDensity()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
KSDensity::~KSDensity()
{
}

// ************************************************************************
// get pdf
// ------------------------------------------------------------------------
void KSDensity::genDensity1D(psVector &dataSet, psVector &Xp, psVector &Pp)
{
  int    outLeng=101, ii, jj, kk;
  double lower, upper, width, median, accum, alpha, pi2Inv=0.5/3.1415928;
  double ddata, ddata2;
  psVector ds;

  ds = dataSet;
  ds.sort();
  kk = ds.length() / 2;
  median = ds[kk]; 
  for (ii = 0; ii < ds.length(); ii++) ds[ii] = PABS(ds[ii] - median);
  ds.sort();
  median = ds[kk]; 
  alpha = median / 0.6745;
  if (alpha == 0) alpha = dataSet.max() - dataSet.min();
  alpha *= pow(4.0 / (3.0 * dataSet.length()), 0.2);
  if (alpha == 0) alpha = 1.0;
  lower = dataSet.min() - 2.0 * alpha;
  upper = dataSet.max() + 2.0 * alpha;
  width = (upper - lower) / outLeng;

  Xp.setLength(outLeng);
  Pp.setLength(outLeng);
  for (ii = 0; ii < outLeng; ii++)
  {
    Xp[ii] = ddata = lower + width * (0.5 + ii);
    accum = 0.0;
    for (jj = 0; jj < dataSet.length(); jj++)
    {
      ddata2 = (ddata - dataSet[jj]) / alpha;
      accum += exp(-0.5 * ddata2 * ddata2) * pi2Inv;
    }
    Pp[ii] = accum / (double) dataSet.length();
  }
  ddata = Pp.sum();
  if (ddata != 0.0)
  {
    ddata = 1.0 / ddata;
    Pp.scale(ddata);
  }
}

