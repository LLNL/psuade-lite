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
// Definition for the class DeltaAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#ifndef __DELTAANALYZERH__
#define __DELTAANALYZERH__

#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************

class DeltaAnalyzer : public Analyzer
{
private:

   int    mode_;
   int    nBins_;
   int    nInputs_;
   int    nConfig_;
   double *minDeltas_; //length nBins_
   int    **deltaBins_;   //length nBins_, nInputs_
   double *dOrder_;    //length nInputs_
   int    *ranks_; 	   //length nInputs_

public:

   DeltaAnalyzer();

   ~DeltaAnalyzer();

   double analyze(aData &adata);

   DeltaAnalyzer& operator=(const DeltaAnalyzer &analyzer);

   int setParams(int nParams, char **params);

   /** Getters for analysis results */
   int    get_mode();
   int    get_nBins();
   int    get_nInputs();
   int    get_nConfig();
   double *get_minDeltas();
   int    **get_deltaBins();
   double *get_dOrder();
   int    *get_ranks();
};

#endif // __DELTAANALYZERH__

