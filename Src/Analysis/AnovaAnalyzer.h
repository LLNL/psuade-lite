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
// Definition for the class AnovaAnalyzer
// AUTHOR : CHARLES TONG
// DATE   : 2003 (updated 7/9/2004)
// ************************************************************************
#ifndef __ANOVAANALYZERH__
#define __ANOVAANALYZERH__

#include "Analyzer.h"

class AnovaAnalyzer : public Analyzer 
{

public:

   AnovaAnalyzer();

   ~AnovaAnalyzer();

   double analyze(aData &adata);

   AnovaAnalyzer& operator=(const AnovaAnalyzer &analyzer);

   /** Getters for analysis results */
   int    get_tableLeng();
   int    *get_dofs();
   double *get_sumSquare();
   double *get_meanSquares();
   double *get_fValues();
   int    **get_code();

private:

   double computeSumSquares1(int,int,int,int,int,int,double*,double*);
   double computeSumSquares2(int,int,int,int,int,int,int,int,double*,
                             double*,double,double);
   double computeSumSquares3(int,int,int,int,int,int,int,int,int,int,
                             double*,double*,double,double,double,double,
                             double,double);

   /*length of arrays is tableLeng_ */
   /*added - 4/30/14 - JMcE*/
   int    tableLeng_;
   int    *dofs_;
   double *sumSquares_;
   double *meanSquares_;
   double *fValues_;
   int    **code_;
};

#endif // __ANOVAANALYZERH__

