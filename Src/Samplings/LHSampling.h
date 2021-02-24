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
// Definition for the Latin hypercube class 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __LHSAMPLINGH_
#define __LHSAMPLINGH_

#include "Sampling.h"

/**
 * @name Latin Hypercube samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class LHSampling : public Sampling 
{
   int    trueRandom_;        // 0 - semi-randomize for main effect
   int    nSymbols_;

public:

   /** constructor */
   LHSampling();

   /** destructor */
   ~LHSampling();

   /** initialization 
       @param flag: flag to signal how far to initialize
    */
   int initialize(int flag);

   /** This function refines an incoming sample
       @param ratio: refinement ratio
       @param randomize: generate randomized sample
       @param thresh: threshold
       @param nSamples: sample size
       @param sampleErrs: errors for each sample point
    */
   int refine(int ratio,int randomize,double thresh,int nSamples,
              double *sampleErrs);

   /** This function overloads the assignment operator
       @param obj : Sampling object
    */
   LHSampling& operator=(const LHSampling &);

   /** This function sets internal input parameters
       @param nInputs : number of inputs
       @param counts : number of symbols per input
       @param settings : input settings
       @param symTable : symbol table
    */
   int setInputParams(int nInputs,int *counts,double **settings,int *symTable);
};

#endif // __LHSAMPLINGH_

