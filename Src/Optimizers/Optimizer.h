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
// Definitions for the abstract class Optimizer
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#ifndef __OPTIMIZER__
#define __OPTIMIZER__

#include <string.h>
//#include <iostream>
//using namespace std;

#include "Optimizer.h"
#include "oData.h"
#include "PsuadeData.h"
#include "FunctionInterface.h"

// ************************************************************************
// class definition
// ************************************************************************
class Optimizer
{
public:

   double *optimalX_;
   double optimalY_;
   int    numEvals_;
   void   (*objFunction_)(int, double *, int, double *);

   /** constructor */
   Optimizer();

   /** destructor */
   virtual ~Optimizer();

   /** send in the function to compute the objective */
   void setObjectiveFunction(void (*func)(int,double*,int,double*));

   /** run optimization 
     @param odata - an object that contains all needed data
     */
   virtual void optimize(oData *odata);

   /** set internal parameters
     @param sparam - string that contains the parameter information
     */
   void setParam(char *sparam);

   /** set optimal input values */
   double *getOptimalX();

   /** set optimal output */
   double getOptimalY();

   /** set number of function evaluations */
   int    getNumFuncEvals();

   /** assign operator override */
   Optimizer& operator=(const Optimizer &);
};

// ************************************************************************
// associated external function
// ************************************************************************
int OptimizerSearch(PsuadeData *, FunctionInterface *, double **, int *);

#endif // __OPTIMIZER__

