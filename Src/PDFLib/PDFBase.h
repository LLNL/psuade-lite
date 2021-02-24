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
// Definition for the abstract class PDFBase
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#ifndef __PDFBASEH__
#define __PDFBASEH__

#include "PsuadeData.h"

/**
 * @name Base class for probability density functions
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class PDFBase 
{

public:

  PDFBase();    
  virtual ~PDFBase();
  virtual int getPDF(int, double *, double *);
  virtual int getCDF(int, double *, double *);
  virtual int invCDF(int, double *, double *, double, double);
  virtual int genSample(int, double *, double *, double *);
  virtual double getMean();
  virtual int getMeans(double *);
  virtual int setParam(char *);
};

// ************************************************************************
// friend function
// ------------------------------------------------------------------------
                                                                                
extern "C"
{
  int PDFTransform(PsuadeData *psuadeIO, int nSamples, int nInputs,
                 double *sampleData);
  int PDFTransformSingle(int ptype, int nSamples, double *sampleData,
                 double lbound, double ubound, double mean, double stdev);
}

#endif // __PDFBASEH__

