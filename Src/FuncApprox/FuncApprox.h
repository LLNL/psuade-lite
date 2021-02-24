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
// Definition for the class FuncApprox
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __FUNCAPPROXH__
#define __FUNCAPPROXH__

#include "PsuadeData.h"

// ************************************************************************
// class definition
// ************************************************************************
class FuncApprox 
{
public:
   int    outputLevel_;
   int    faID_;
   int    nSamples_;
   int    nInputs_;
   int    nPtsPerDim_;
   double *lowerBounds_;
   double *upperBounds_;
   double *weights_;
   double *XMeans_;
   double *XStds_;
   double YMean_;
   double YStd_;

   FuncApprox(int, int);
   // Bill Oliver add prototype for copy constructor
   FuncApprox(const FuncApprox & fa);

   // Bill Oliver added overload operator=
   FuncApprox & operator=(const FuncApprox & fa);

   virtual ~FuncApprox();

   int            getID();
   virtual int    setOutputLevel(int);
   virtual int    setBounds(double *, double *);
   virtual int    loadWeights(int, double *);
   virtual void   setNPtsPerDim(int);
   virtual int    getNPtsPerDim();
   virtual int    initialize(double*,double*)=0;
   virtual int    genNDGridData(double*,double*,int*,double**,double**);
   virtual int    gen1DGridData(double*,double*,int,double*,int*,
                                double**,double**);
   virtual int    gen2DGridData(double*,double*,int,int,double*,int*,
                                double**,double**);
   virtual int    gen3DGridData(double*,double*,int,int,int,double*,int*,
                                double**,double**);
   virtual int    gen4DGridData(double*,double*,int,int,int,int,double*,
                                int*,double**,double**);
   virtual double evaluatePoint(double *)=0;
   virtual double evaluatePoint(int, double *, double *)=0;
   virtual double evaluatePointFuzzy(double *, double &);
   virtual double evaluatePointFuzzy(int, double *, double *, double *);
   virtual double setParams(int, char **);
   int            initInputScaling(double *, double *, int);
   int            initOutputScaling(double *, double *);

public:
   virtual int    genNDGrid(int*,double**);

};

// ************************************************************************
// friend function
// ------------------------------------------------------------------------
extern "C" 
{
   int        getFAType(char *pString);
   void       printThisFA(int faType);
   int        writeFAInfo(int);
   FuncApprox *genFA(int faType, int nInputs, int nOutputs, int nSamples);
   FuncApprox *genFAInteractive(PsuadeData *psuadeIO, int flag);
   FuncApprox *genFAFromFile(char *name, int outputID);
}

#endif // __FUNCAPPROXH__

