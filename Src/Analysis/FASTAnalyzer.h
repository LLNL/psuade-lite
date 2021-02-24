// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team. 
// All rights reserved.
//
// Please see the COPYRIGHT and LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Definitions for the class FASTAnalyzer
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#ifndef __FASTANALYZERH__
#define __FASTANALYZERH__

#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class FASTAnalyzer : public Analyzer
{
private:

   int    nInputs_;
   int    M_;
   double *fourierCoefs_; //length nInputs_
   double FASTvariance_;

public:

   FASTAnalyzer();

   ~FASTAnalyzer();

   double analyze(aData &adata);

   FASTAnalyzer& operator=(const FASTAnalyzer &analyzer);

   int calculateOmegas(int nInputs, int nSamples, int *omegas);

   int computeCoefficents(int nSamples, int nInputs, double *Y, 
                          double **coefs, int outputLevel);

   /** Getters for analysis results */
   int    get_nInputs();
   int    get_M();
   double *get_fourierCoefs();
   double get_FASTvariance();
};

#endif // __FASTANALYZERH__

