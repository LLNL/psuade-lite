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
// Functions for the class TxMathOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
// ************************************************************************
#ifdef HAVE_TXMATH

#include "txc_streams.h"
#include "TxThroughStream.h"
#include "testOptimizer.h"
#include "TxPowellOpt.h"
#include "TxMathOptimizer.h"
#include "TxMathPSUADE.h"
#include "TxNoDerivPtrFunc.h"

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
TxMathOptimizer::TxMathOptimizer()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
TxMathOptimizer::~TxMathOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void TxMathOptimizer::optimize(oData *odata)
{
   int nInputs;
   TxOptimizer<double,TXSTD::vector<double> >* optimizer;
   TxThroughStream tts;

   nInputs = odata->nInputs_;
   TxMathPSUADEInitialize(odata);

   TxNoDerivPtrFunc <double, TXSTD::vector<double>,
         double (*)(const TXSTD::vector<double>&) >
         fpsuadeFunctor(TxMathPSUADEFunction, nInputs);

   optimizer = new
      TxPowellOpt<double,TXSTD::vector<double> >(&fpsuadeFunctor);
   optimizer->getFunctor().setClassName("PSUADE");
   testOptimizer(optimizer, TxMathPSUADEGetStartingPoint(), 
                 TxMathPSUADEGetDimension(), tts);
   odata->optimalY_ = TxMathPSUADEGetOptimalData(odata->optimalX_);
   odata->numFuncEvals_ = TxMathPSUADEGetNumFuncEval();
   delete optimizer;
   TxMathPSUADECleanUp();
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
TxMathOptimizer& TxMathOptimizer::operator=(const TxMathOptimizer &)
{
   printf("TxMathOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}
#else
   int txmath_bogus;
#endif


