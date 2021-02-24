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
// Definition for the class MainEffectAnalyzer (McKay's VCE analysis)
// AUTHOR : CHARLES TONG
// DATE   : 2003 (updated 9/04)
// ************************************************************************

#ifndef __MAINEFFECTANALYZERH__
#define __MAINEFFECTANALYZERH__

#include "Analyzer.h"
#include "PsuadeData.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class MainEffectAnalyzer : public Analyzer
{
private:

   int    matlabPlotFlag_;
   int    nInputs_;
   int    outputID_;
   double *inputVCE_; 	//length is nInputs_
   double totalInputVCE_;
   double mainEffectMean_;
   double mainEffectStd_;

public:

   MainEffectAnalyzer();

   ~MainEffectAnalyzer();

   double analyze(aData &adata);

   int printResults(int, double, double *, PsuadeData *);

   MainEffectAnalyzer& operator=(const MainEffectAnalyzer &analyzer);

   int plotResponse(int nInputs, int nSamples, double *sampleInputs,
                    double *sampleOutputs, double *xLower, double *xUpper, 
                    int plotAxis1, int plotAxis2, double *settings);

   double analyzeOAT(aData &);

public:
   int computeMeanVariance(int,int,int,double*,double *,double *,int);
   int computeVCECrude(int,int,double *,double *,double *,double *,
                       double,double *);
   int computeVCE(int, int, int, double *, double *, int, FILE *,
                  double *, double *, double *, double *);

   /** Getters for analysis results */
   int    get_nInputs();
   int    get_outputID();
   double *get_inputVCE();
   double get_totalInputVCE();
   double get_mainEffectMean();
   double get_mainEffectStd();
};

#endif // __MAINEFFECTANALYZERH__

