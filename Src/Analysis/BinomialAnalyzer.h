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
// Definition for the class BinomialAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#ifndef __BINOMIALANALYZERH__
#define __BINOMIALANALYZERH__

#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
class BinomialAnalyzer : public Analyzer
{

public:

   BinomialAnalyzer();

   ~BinomialAnalyzer();

   double analyze(aData &adata);

   BinomialAnalyzer& operator=(const BinomialAnalyzer &analyzer);

   double *setupBinomialCDF(int n, double p0);

   double factorial(int n);

   /** Getters for analysis results */
   int    get_nSamples();
   double *get_Ybin();
   int    get_nBelow();
   double *get_BinomialCDF();
   double get_typeI();

private:
   int nSamples_;
   double *Ybin_; //length = nSamples_
   int    nBelow_;
   double *BinomialCDF_; //length = nSamples_+1
   double typeI_;
};

#endif // __BINOMIALANALYZERH__

