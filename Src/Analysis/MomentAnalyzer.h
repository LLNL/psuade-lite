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
// Definition for the class MomentAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************

#ifndef __MOMENTANALYZERH__
#define __MOMENTANALYZERH__

#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class MomentAnalyzer : public Analyzer
{
private:

	int nSamples_;
	int nGroups_;
	int nInputs_;
	int nOutputs_;
	double *moments_;

public:

   MomentAnalyzer();

   ~MomentAnalyzer();

   void analyze(int nInps, int nSamp, double *ilbs,
                double *iubs, double *sInps, double *sOuts);

   double analyze(aData &adata);

   MomentAnalyzer& operator=(const MomentAnalyzer &analyzer);

   int analyzeMore(int nInputs, int nOutputs, int nSamples, 
                   int nSubSamples, double *sampleOutputs, int outputID);

   int computeMeanVariance(int nInputs, int nOutputs, int nSamples,
                           double * sampleOutputs, double *mean,
                           double *variance,int outputID);

   int computeSkewnessKurtosis(int nInputs, int nOutputs, int nSamples,
                      double * sampleOutputs, double stdev,
                      double *skewness, double *kurtosis,int outputID);

   /** Getters for analysis results */
   int get_nSamples();
   int get_nGroups();
   int get_nInputs();
   int get_nOutputs();
   double get_mean();
   double get_stdev();
   double get_skewness();
   double get_kurtosis();
};

#endif // __MOMENTANALYZERH__

