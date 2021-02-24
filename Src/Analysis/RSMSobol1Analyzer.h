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
// Definition for the class RSSobol1Analyzer
// (Sobol variance decomposition for one parameter - with response surface)
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#ifndef __RSMSOBOL1ANALYZERH__
#define __RSMSOBOL1ANALYZERH__

#include "Analyzer.h"
#include "Vector.h"

// ************************************************************************
// class definition
// ************************************************************************
class RSMSobol1Analyzer : public Analyzer
{
private:

   int    nInputs_;
   double outputMean_;
   double outputStd_;
   psVector vecVces_;  //length is nInputs_

public:

   RSMSobol1Analyzer();

   ~RSMSobol1Analyzer();

   void analyze(int nInps, int nSamp, double *ilbs,
                double *iubs, double *sInps, double *sOuts);

   double analyze(aData &adata);

   double analyze2(aData &adata);

   double analyze3(aData &adata);

   RSMSobol1Analyzer& operator=(const RSMSobol1Analyzer &analyzer);

   /** Getters for analysis results */
   int    get_nInputs();
   double get_outputMean();
   double get_outputStd();
   double get_vce(int);
};

#endif // __RSMSOBOL1ANALYZERH__

