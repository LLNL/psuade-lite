// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Definition for the class SparseGridRegression
// AUTHOR : CHARLES TONG
// DATE   : 2011
// ************************************************************************

#ifndef __SGRIDREGRESSIONH__
#define __SGRIDREGRESSIONH__

#include "FuncApprox.h"
#include "Vector.h"

// ************************************************************************
// class definition
// ************************************************************************
class SparseGridRegression : public FuncApprox 
{
   int    pOrder_;
   int    numPerms_;
   int    **pcePerms_;
   psVector sampleInputs_;
   psVector sampleWeights_;
   psVector regCoefs_;
   psVector lBounds_;
   psVector uBounds_;
 
public:
   SparseGridRegression(int, int);
   ~SparseGridRegression();

   int    initialize(double*,double*);
   int    genNDGridData(double*,double*,int*,double**,double**);
   int    gen1DGridData(double*,double *,int,double*, 
                        int *, double **, double **);
   int    gen2DGridData(double*,double *,int,int,double*, 
                        int *, double **, double **);
   int    gen3DGridData(double*,double *,int,int,int,double*, 
                        int *, double **, double **);
   int    gen4DGridData(double*,double *,int,int,int,int,double*, 
                        int *, double **, double **);
   double evaluatePoint(double *);
   double evaluatePoint(int, double *, double *);
   double evaluatePointFuzzy(double *, double &);
   double evaluatePointFuzzy(int, double *, double *, double *);

private:
   int    analyze(double *, double *);
   int    GenPermutations();
   int    EvalLegendrePolynomials(double, double *);
};

#endif // __SGRIDREGRESSIONH__

