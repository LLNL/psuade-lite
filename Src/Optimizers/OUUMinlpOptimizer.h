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
// Definitions for the class OUUMinlpOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#ifndef __OUUMINLPOPTIMIZER__
#define __OUUMINLPOPTIMIZER__

#include "Optimizer.h"
#include "oData.h"

#if HAVE_NOMAD
#include "../../External/NOMAD/src/Evaluator.hpp"
#include "../../External/NOMAD/src/Eval_Point.hpp"
#include "../../External/NOMAD/src/Parameters.hpp"

// ************************************************************************
// evaluator class definition
// ************************************************************************
class PsuadeMinlpEvaluator: public NOMAD::Evaluator
{
public:

   /** constructor */
   PsuadeMinlpEvaluator(const NOMAD::Parameters &pset);

   /** destructor */
   ~PsuadeMinlpEvaluator();

   /** evaluation of a point */
   bool eval_x(NOMAD::Eval_Point &,const NOMAD::Double &,bool &);
};
#endif

// ************************************************************************
// class definition
// ************************************************************************
class OUUMinlpOptimizer : public Optimizer
{
   int repeatFlag_;
   int optCode_;

public:

   /** constructor */
   OUUMinlpOptimizer();

   /** destructor */
   ~OUUMinlpOptimizer();

   /** run optimization 
     @param odata - an object that contains all needed data
     */
   void optimize(oData *odata);

   /** clean up */ 
   void cleanUp();

   /** assign operator override */ 
   OUUMinlpOptimizer& operator=(const OUUMinlpOptimizer &);

private:
   /** generate Z3 sample */ 
   void genZ3Sample();

   /** generate Z4 sample */ 
   void genZ4Sample();

   /** generate response surface */ 
   void genResponseSurface();
};

#endif // __OUUMINLPOPTIMIZER__

