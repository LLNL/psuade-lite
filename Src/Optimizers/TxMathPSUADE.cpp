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
// TxMathOptimizer routine : perform function evaluation 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifdef HAVE_TXMATH

#include <stdio.h>
#include "TxMathPSUADE.h"

//--------------------------------------------------------------------------
// local variables
//--------------------------------------------------------------------------
FunctionInterface *TxMath_funcIO_=NULL;
int TxMath_nInputs_=0;
int TxMath_nOutputs_=0;
int TxMath_outputID_=0;
int TxMath_outputLevel_=0;
int TxMath_sampleID_=0;
double *TxMath_startingPoint_=NULL;
double *TxMath_inputVector_=NULL;
double *TxMath_optimalX_=NULL;
double *TxMath_xLower_=NULL;
double *TxMath_xUpper_=NULL;
double TxMath_optimalY_=1.0e49;

//**************************************************************************
// function evaluation 
//--------------------------------------------------------------------------
void TxMathPSUADEInitialize(oData *odata)
{
   int    nInputs, nOutputs;
   double *xLower, *xUpper, *start;

   TxMath_funcIO_  = odata->funcIO_;
   nInputs  = odata->nInputs_;
   nOutputs = odata->nOutputs_;
   xLower   = odata->lowerBounds_;
   xUpper   = odata->upperBounds_;
   start    = odata->initialX_;
   TxMath_outputID_ = odata->outputID_;

   if (nInputs <= 0 || nOutputs != 1)
   {
      printf("TxMathPSUADE ERROR : invalid nInputs (%d), nOutputs (%d).\n",
             nInputs, nOutputs);
      exit(1);
   }
   TxMath_nInputs_ = nInputs;
   TxMath_nOutputs_ = nOutputs;
   if (TxMath_startingPoint_ != NULL) delete [] TxMath_startingPoint_;
   TxMath_startingPoint_ = new double[TxMath_nInputs_]; 
   for (int i = 0; i < TxMath_nInputs_; i++) 
      TxMath_startingPoint_[i] = start[i];
   if (TxMath_inputVector_ != NULL) delete [] TxMath_inputVector_;
   TxMath_inputVector_ = new double[TxMath_nInputs_]; 
   if (TxMath_optimalX_ != NULL) delete [] TxMath_optimalX_;
   TxMath_optimalX_ = new double[TxMath_nInputs_]; 
   for (int j = 0; j < TxMath_nInputs_; j++) TxMath_optimalX_[j] = start[j];
   if (TxMath_xLower_ != NULL) delete [] TxMath_xLower_;
   TxMath_xLower_ = new double[TxMath_nInputs_]; 
   if (TxMath_xUpper_ != NULL) delete [] TxMath_xUpper_;
   TxMath_xUpper_ = new double[TxMath_nInputs_]; 
   for (int k = 0; k < TxMath_nInputs_; k++) 
   {
      TxMath_xLower_[k] = xLower[k];
      TxMath_xUpper_[k] = xUpper[k];
   }
   TxMath_optimalY_ = 1.0e49;
   TxMath_sampleID_ = 0;
   TxMath_outputLevel_ = odata->outputLevel_;
}

//**************************************************************************
// function evaluation 
//--------------------------------------------------------------------------
double TxMathPSUADEFunction(const TXSTD::vector<double>& x) 
{
   int    i, evaluateFlag=1, currDriver;
   double *sampleOutputs, sampleOutput, addend;

   if (TxMath_funcIO_ == NULL)
   {
      printf("TxMathPSUADEFunction ERROR : function pointer missing.\n");
      exit(1);
   }
   for (i = 0; i < TxMath_nInputs_; i++) TxMath_inputVector_[i] = x[i];
   sampleOutput = addend = 0.0;
   sampleOutputs = new double[TxMath_nOutputs_];
   for (i = 0; i < TxMath_nInputs_; i++) 
   {
      if (TxMath_inputVector_[i] < TxMath_xLower_[i]) 
      {
         addend += (1.0e20 * (TxMath_xLower_[i] - TxMath_inputVector_[i]));
         evaluateFlag = 0;
      }
      if (TxMath_inputVector_[i] > TxMath_xUpper_[i]) 
      {
         addend += (1.0e20 * (TxMath_inputVector_[i] - TxMath_xUpper_[i]));
         evaluateFlag = 0;
      }
   }
   if (evaluateFlag == 1)
   {
      currDriver = TxMath_funcIO_->getDriver();
      TxMath_funcIO_->setDriver(1);
      TxMath_funcIO_->evaluate(TxMath_sampleID_,TxMath_nInputs_,
                        TxMath_inputVector_,TxMath_nOutputs_,
                        sampleOutputs,0);
      TxMath_funcIO_->setDriver(currDriver);
      TxMath_sampleID_++;
      sampleOutput = sampleOutputs[TxMath_outputID_];
   }
   sampleOutput += addend;
   if (sampleOutput < TxMath_optimalY_)
   {
      TxMath_optimalY_ = sampleOutput;
      for (i = 0; i < TxMath_nInputs_; i++) TxMath_optimalX_[i] = x[i];
   }
   if (TxMath_outputLevel_ > 0)
   {
      printf("TxMathOutput %6d = %16.8e\n", TxMath_sampleID_, sampleOutput);
      printf("TxMathOptimizer current Ymin = %16.8e\n", TxMath_optimalY_);
   }
   delete [] sampleOutputs;
   return sampleOutput;
}

//**************************************************************************
// define starting point
//--------------------------------------------------------------------------
TXSTD::vector<double> TxMathPSUADEGetStartingPoint() 
{
   TXSTD::vector<double> point(TxMath_nInputs_);
   for (int i = 0; i < TxMath_nInputs_; i++) 
      point[i] = TxMath_startingPoint_[i];
   return point;
}

//**************************************************************************
// define input dimension
//--------------------------------------------------------------------------
size_t TxMathPSUADEGetDimension()
{
   return TxMath_nInputs_;
}

//**************************************************************************
// fetch optimal data 
//--------------------------------------------------------------------------
double TxMathPSUADEGetOptimalData(double *optX)
{
   for (int i = 0; i < TxMath_nInputs_; i++ ) optX[i] = TxMath_optimalX_[i];
   return TxMath_optimalY_;
}

//**************************************************************************
// clean up 
//--------------------------------------------------------------------------
void TxMathPSUADECleanUp()
{
   if (TxMath_startingPoint_ != NULL) delete [] TxMath_startingPoint_;
   if (TxMath_inputVector_ != NULL) delete [] TxMath_inputVector_;
   if (TxMath_optimalX_ != NULL) delete [] TxMath_optimalX_;
   if (TxMath_xLower_ != NULL) delete [] TxMath_xLower_;
   if (TxMath_xUpper_ != NULL) delete [] TxMath_xUpper_;
   TxMath_startingPoint_ = NULL;
   TxMath_inputVector_ = NULL;
   TxMath_optimalX_ = NULL;
   TxMath_xLower_ = NULL;
   TxMath_xUpper_ = NULL;
}

//**************************************************************************
// get number of function evaluations
//--------------------------------------------------------------------------
int TxMathPSUADEGetNumFuncEval()
{
   return (TxMath_sampleID_+1);
}

#else
   int txmath_psuade_bogus=0;
#endif

