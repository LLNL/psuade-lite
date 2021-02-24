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
// Functions for the class BinomialAnalyzer (binomial test)
// Reference: "Hypothesis Testing Procedures for Eliminating Variables of 
//            Negligible Impact in Large Scale Computer Simulations"
//            by Jason Lenderman
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "BinomialAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
BinomialAnalyzer::BinomialAnalyzer() : Analyzer(), nSamples_(0), Ybin_(0), 
                     nBelow_(0), BinomialCDF_(0), typeI_(0)
{
   setName("BINOMIAL");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BinomialAnalyzer::~BinomialAnalyzer()
{
   if (Ybin_) delete [] Ybin_;
   if (BinomialCDF_) delete BinomialCDF_;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double BinomialAnalyzer::analyze(aData &adata)
{
   int     ss,  printLevel, info;
   double  *Y, p0=0.5, thresh;

   printLevel = adata.printLevel_;
   nSamples_   = adata.nSamples_;
   Y          = adata.sampleOutputs_;
   thresh     = adata.analysisThreshold_;

   if (nSamples_ <= 1)
   {
      printOutTS(PL_ERROR, "BinomialAnalyzer INFO: not meaningful to ");
      printOutTS(PL_ERROR, "do this when nSamples_ <= 1.\n");
      return PSUADE_UNDEFINED;
   } 
   if (Y == NULL)
   {
      printOutTS(PL_ERROR, "BinomialAnalyzer ERROR: no data.\n");
      return PSUADE_UNDEFINED;
   } 
   info = 0;
   for (ss = 0; ss < nSamples_; ss++)
      if (Y[ss] == PSUADE_UNDEFINED) info++;
   if (info > 0)
   {
      printOutTS(PL_ERROR,
           "BinomialAnalyzer ERROR: Some outputs are undefined.\n");
      printOutTS(PL_ERROR,
           "                        Prune them first before analyze.\n");
      return PSUADE_UNDEFINED;
   }

   if (Ybin_) delete [] Ybin_;
   if (BinomialCDF_) delete BinomialCDF_;
   Ybin_ = NULL;
   BinomialCDF_ = NULL;
   
   Ybin_ = new double[nSamples_];
   checkAllocate(Ybin_, "Ybin in BinomialAnalyzer::analyze");
   for (ss = 0; ss < nSamples_; ss++)
   {
      if (PABS(Y[ss]) <= thresh) Ybin_[ss] = 1.0;
      else                       Ybin_[ss] = 0.0;
   }

   nBelow_ = 0;
   for (ss = 0; ss < nSamples_; ss++) if (Ybin_[ss] == 1.0) nBelow_++;
   if (printLevel > 2)
      printOutTS(PL_INFO,
           "   ====> BinomialAnalyzer: number below threshold = %d (%d)\n",
           nBelow_, nSamples_);

   BinomialCDF_ = setupBinomialCDF(nSamples_, p0);
   for (ss = 0; ss <= nSamples_; ss++)
      BinomialCDF_[ss] = 1.0 - BinomialCDF_[ss];

   typeI_ = BinomialCDF_[nBelow_];

   //delete [] Ybin_;
   //delete [] BinomialCDF_;
   return typeI_;
}

// *************************************************************************
// create a binomial cumulative density function
// -------------------------------------------------------------------------
double *BinomialAnalyzer::setupBinomialCDF(int n, double p0)
{
   int    k;
   double *retValues, nFact, kFact, nkFact;

   nFact = factorial(n);
   retValues = new double[n+1];
   checkAllocate(retValues,"retValues in Binomial::setupBinomialCDF");
   for (k = 0; k <= n; k++)
   {
      kFact  = factorial(k);
      nkFact = factorial(n-k);
      retValues[k] = pow(p0,1.0*k) * pow(1-p0,1.0*(n-k));
      retValues[k] *= (nFact / (kFact * nkFact));
      if (k > 0) retValues[k] += retValues[k-1];
   }
   return retValues;
}

// *************************************************************************
// Compute factorial
// -------------------------------------------------------------------------
double BinomialAnalyzer::factorial(int n)
{
   int    ii;
   double fact;

   fact = 1.0;
   if (n == 0) return fact;
   for (ii = 2; ii <= n; ii++) fact *= (double) ii;
   return fact;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
BinomialAnalyzer& BinomialAnalyzer::operator=(const BinomialAnalyzer &)
{
   printOutTS(PL_ERROR, 
        "BinomialAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// Getters for analysis results
// ------------------------------------------------------------------------
int BinomialAnalyzer::get_nSamples()
{
   return nSamples_;
}
double *BinomialAnalyzer::get_Ybin()
{
   double* retVal = NULL;
   if (Ybin_)
   {
      retVal = new double[nSamples_];
      checkAllocate(retVal,"retVal in Binomial::get_nSamples");
      std::copy(Ybin_, Ybin_+nSamples_, retVal);
   }
   return retVal;
}
int BinomialAnalyzer::get_nBelow()
{
   return nBelow_;
}
double *BinomialAnalyzer::get_BinomialCDF()
{
   double* retVal = NULL;
   if (Ybin_)
   {
      retVal = new double[nSamples_+1];
      checkAllocate(retVal,"retVal in Binomial::get_BinomialCDF");
      std::copy(BinomialCDF_, BinomialCDF_+nSamples_+1, retVal);
   }
   return retVal;
}
double BinomialAnalyzer::get_typeI()
{
   return typeI_;
}

