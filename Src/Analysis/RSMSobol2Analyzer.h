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
// Definition for the class RSMSobol2Analyzer
// (Sobol variance decomposition for second order sensitivity)
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#ifndef __RSMSOBOL2ANALYZERH__
#define __RSMSOBOL2ANALYZERH__
#include <vector>
#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class RSMSobol2Analyzer : public Analyzer
{
private:

   int    nInputs_;
   double outputMean_;
   double outputStd_;
   double *vces_;
   double *ecvs_;

public:

   RSMSobol2Analyzer();

   ~RSMSobol2Analyzer();

   void analyze(int nInps, int nSamp, double *ilbs,
                double *iubs, double *sInps, double *sOuts);

   double analyze(aData &adata);

   double analyze2(aData &adata);

   RSMSobol2Analyzer& operator=(const RSMSobol2Analyzer &analyzer);

   /** Getters for analysis results */
   int    get_nInputs();
   double get_outputMean();
   double get_outputStd();
   double get_vce(int, int);
   double get_ecv(int, int);
};

#endif // __RSMSOBOL2ANALYZERH__

