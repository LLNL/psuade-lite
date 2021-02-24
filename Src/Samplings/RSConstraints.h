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
// Definition for the abstract class RSConstraints
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#ifndef __RSCONSTRAINTS__
#define __RSCONSTRAINTS__

#include "FuncApproxFilter.h"
#include "PsuadeData.h"

/**
 * @name Defining output constraints
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class RSConstraints
{

   int              nConstraints_;
   int              *constraintNInputs_;
   int              **constraintInputIndices_;
   double           **constraintInputValues_;
   int              nInputs_;
   FuncApproxFilter **constraintFAs_;
   double           *YLBounds_;
   double           *YUBounds_;

public:

   /** constructor */
   RSConstraints();

   /** destructor */
   ~RSConstraints();

   /** initialization 
       @param flag: flag to signal how far to initialize
    */
   int initialize(int flag);

   /** get number of constraints 
    */
   int getNumConstraints();

   /** This function refines an incoming sample
       @param ratio: refinement ratio
       @param randomize: generate randomized sample
       @param thresh: threshold
       @param nSamples: sample size
       @param sampleErrs: errors for each sample point
    */
   int refine(int ratio,int randomize,double thresh,int nSamples,double *sampleErrs);

   /** This function overloads the assignment operator
       @param obj : Sampling object
    */
   RSConstraints& operator=(const RSConstraints &);

   /** This function generates the constraints
       @param psuadeIO : data container
    */
   int    genConstraints(PsuadeData *psuadeIO);

   /** This function evaluates the constraints
       @param samInputs : sample inputs
       @param samOutputs : sample output
       @param flag : return status flag
    */
   double evaluate(double *samInputs, double samOutput, int &flag);
};

#endif // __RSCONSTRAINTS__

