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
// Definition for the class RSFuncApproxAnalyzer
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __RSFUNCAPPROXANALYZERH__
#define __RSFUNCAPPROXANALYZERH__

#include "Analyzer.h"

// ************************************************************************
// class definition
// ************************************************************************
class RSFuncApproxAnalyzer : public Analyzer
{
   int rsType_;
   int useCV_;
   int numCVGroups_;

public:

   RSFuncApproxAnalyzer();

   ~RSFuncApproxAnalyzer();

   double analyze(aData &adata);

   RSFuncApproxAnalyzer& operator=(const RSFuncApproxAnalyzer &analyzer);

   int setParams(int nParams, char **params);

private:

   double validate(aData &adata, char *datafile, char *errfile);
};

#endif // __RSFUNCAPPROXANALYZERH__

