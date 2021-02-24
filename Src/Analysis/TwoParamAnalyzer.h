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
// Definition for the class TwoParamAnalyzer (interaction analysis)
// AUTHOR : CHARLES TONG
// DATE   : updated in 2006
// ************************************************************************

#ifndef __TWOPARAMANALYZERH__
#define __TWOPARAMANALYZERH__

#include "Analyzer.h"
#undef max
#include <vector>

// ************************************************************************
// class definition
// ************************************************************************
class TwoParamAnalyzer : public Analyzer
{
private:

   int nInputs_;
   double outputMean_;
   double outputStd_;
   double *sensitivity_;	//length is nInputs_
   struct record
   {
      int index1;
      int index2;
      double value;
   };
   std::vector<record> compareSensitivity_;

   int computeMeanVariance(int,int,int,double*,double *,double *,int);
   int computeVCE2(int,int,int,double*,double*,FILE *,double,
                   double *, double *, double *, double *);
   int computeVCECrude(int, int, double *, double *, double, double *,
                       double *, double *);
public:

   TwoParamAnalyzer();

   ~TwoParamAnalyzer();

   double analyze(aData &adata);

   TwoParamAnalyzer& operator=(const TwoParamAnalyzer &analyzer);

   double analyze1D(int nInputs, int nOutputs, int nSamples,
                    double *sampleInputs, double *sampleOutputs,
                    int inputID, int outputID, int nSubSamples, 
                    double mean, double var, double *retdata);
   
   /** Getters for analysis results */
   int    get_nInputs();
   double get_outputMean();
   double get_outputStd();
   double *get_sensitivity();
   std::vector<record>  get_compareSensitivity();
};

#endif // __TWOPARAMANALYZERH__

