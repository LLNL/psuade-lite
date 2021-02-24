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
// Definitions for the class NomadOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#ifndef __NOMADOPTIMIZER__
#define __NOMADOPTIMIZER__

#include "Optimizer.h"
#include "oData.h"
#ifdef HAVE_NOMAD
#include "../../External/NOMAD/src/Evaluator.hpp"
#include "../../External/NOMAD/src/Eval_Point.hpp"
#include "../../External/NOMAD/src/Parameters.hpp"

// ************************************************************************
// class definition
// ************************************************************************
class PsuadeNomadEvaluator: public NOMAD::Evaluator
{
public:

   /** constructor */
   PsuadeNomadEvaluator(const NOMAD::Parameters &pset);

   /** destructor */
   ~PsuadeNomadEvaluator();

   /** evaluation of a point */
   bool eval_x(NOMAD::Eval_Point &,const NOMAD::Double &,bool &);
};
#endif

// ************************************************************************
// class definition
// ************************************************************************
class NomadOptimizer : public Optimizer
{
   int  nInputs_;
   int  *inputTypes_;

public:

   /** constructor */
   NomadOptimizer();

   /** destructor */
   ~NomadOptimizer();

   /** set discrete variable  
     @param index - rank of variable that is set to be discrete
     */
   void setDiscreteVariable(int index);

   /** run optimization in library mode (called by user directly)
     @param nInps - number of inputs
     @param XVals - initial guess and solution
     @param lbnds - lower bounds
     @param ubnds - upper bounds
     @param nOuts - number of outputs
     @param maxfe - maximum number of function evaluations
     @param tol   - tolerance
     */
   void optimize(int nInps, double *XVals, double *lbnds,
                 double *ubnds, int nOuts, int maxfe, double tol);

   /** run optimization 
     @param odata - an object that contains all needed data
     */
   void optimize(oData *odata);

   /** assign operator override */ 
   NomadOptimizer& operator=(const NomadOptimizer &);
};

#endif // __NOMADOPTIMIZER__

