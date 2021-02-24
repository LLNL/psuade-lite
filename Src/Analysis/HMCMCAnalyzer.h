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
// Definition for the class HMCMCAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#ifndef __HMCMCANALYZERH__
#define __HMCMCANALYZERH__

#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
class HMCMCAnalyzer : public Analyzer
{
private:
   int    nInputs_;
   double *means_;           //length is nInputs_
   double *sigmas_;          //length is nInputs_
   double *mostLikelyInput_; //length is nInputs_
   double *mostLikelyOutput_;

public:

   HMCMCAnalyzer();

   ~HMCMCAnalyzer();

   double analyze();

   HMCMCAnalyzer& operator=(const HMCMCAnalyzer &analyzer);

   int readUserSpec(int *, int *, double **, int **, int ***, int **);

   int getHierarchicalPriors(int, int *, int **, double *, double *,
                             double *, int *, double *, double *);

   int checkConvergence(int leng, double *values, int step);

   int genMatlabFile(int samSize, int nParams, double *sample, 
                     double *inferenceStore);

   /** Getters for analysis results */
   int    get_nInputs();
   int    get_nOutputs();
   double *get_means();
   double *get_sigmas();
   double *get_mostLikelyInput();
   double *get_mostLikelyOutput();
};

#endif // __HMCMCANALYZERH__

