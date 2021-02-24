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
// Definition for the class BootstrapAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#ifndef __BOOTSTRAPANALYZERH__
#define __BOOTSTRAPANALYZERH__

#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class BootstrapAnalyzer : public Analyzer
{
   int    nSteps_;
   double mean_;
   double stdev_;
   double *storedValues_; //length is nSteps_+1

public:

   BootstrapAnalyzer();

   BootstrapAnalyzer(const BootstrapAnalyzer &);

   ~BootstrapAnalyzer();

   double analyze(aData &adata);

   BootstrapAnalyzer& operator=(const BootstrapAnalyzer &analyzer);

   int computeMeanVariance(int nSamples, int nOutputs, 
                           double *sampleOutputs, double *mean,
                           double *variance,int outputID, int flag);

   int genBootstrap(int nSamples, double *Y, int ntimes, double *bmeans);

   double genJackknife(int nSamples, double *Y);

   int setupNormalCDF(double mean, double stdev);

   double normalCDFInv(double &value);

   /** Getters for analysis results */
   int get_nSteps();
   double get_mean();
   double get_stdev();
   double *get_storedValues();

};

#endif // __BOOTSTRAPANALYZERH__

