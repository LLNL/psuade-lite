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
// Definition for the SequentialSampling class
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************
#ifndef __SEQUENTIALSAMPLINGH__
#define __SEQUENTIALSAMPLINGH__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Sampling.h"
#include "Vector.h"
#include "Matrix.h"

/**
 * @name sequential samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class SequentialSampling: public Sampling 
{
   psVector XMeans_;
   psVector XStds_;

public:

   /** constructor */
   SequentialSampling();

   /** destructor */
   ~SequentialSampling();

   /** initialization 
       @param flag: flag to signal how far to initialize
    */
   int initialize(int flag);

   /** This function overloads the assignment operator
       @param obj : Sampling object
    */
   SequentialSampling& operator=(const SequentialSampling &);

private:
   void constructCMatrix(psMatrix &, psMatrix &, psVector);
   int  genDesigns(psVector&,psMatrix&,int,psMatrix &,psVector &);
   int  genDesigns2(psVector&,psMatrix&,int,psMatrix &,psVector &);
   int  genDesignsUltimate(psVector&,psMatrix&,int,psMatrix &);
   int  genInitialDesigns(psVector,psMatrix&,int,psMatrix &,psVector &);
   int  genInitialDesigns2(psVector,psMatrix&,int,psMatrix &,psVector &);
   int  evaluateUltimate(int, psIVector &, int, int, int, psVector &,
                         psMatrix &, psIVector &, double &, int &);
   int  printDesigns(psMatrix &, FILE *);
};

#endif // __SEQUENTIALSAMPLINGH__

