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
// Definition for the class Kriging
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
#ifndef __KRIGINGH__
#define __KRIGINGH__

#include <stdio.h>
#include "FuncApprox.h"

// ************************************************************************
// class definition
// ************************************************************************
class Kriging : public FuncApprox 
{
   double *Thetas_;
   double *XNormalized_;
   double *YNormalized_;
   double *Rmatrix_;
   double *Mmatrix_;
   double *V1_;
   double *V2_;
   int    workLength_;
   double *workArray_;
   double *workX_;
   double optTolerance_;
   int    pOrder_;
   int    initFlag_;
   int    fastMode_;
   double *dataStdDevs_;
   double KrigingVariance_;
   int    noReuse_;
   double *betas_;
   double *betasOpt_;
   double *gammas_;
   double *gammasOpt_;
   double *thetasOpt_;

public:

   Kriging(int, int);
   ~Kriging();

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

   double train(double *, double *);
   double predict(int, double *, double *, double *);
   int    computeDistances(double **, int *);
   double setParams(int, char **);
   double evaluateFunction(double *);
   void   optimize();
   double predict0(int, double *, double *, double *);
};

#endif // __KRIGINGH__

