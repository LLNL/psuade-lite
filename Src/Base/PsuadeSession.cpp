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
// Functions for the class PsuadeSession  
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "PsuadeSession.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "PDFBase.h"
#include "pData.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
PsuadeSession::PsuadeSession()
{
   owned_ = 0;
   inputLBounds_ = NULL;
   inputUBounds_ = NULL;
   sampleInputs_ = NULL;
   sampleOutputs_ = NULL;
   sampleStates_  = NULL;
   inputPDFs_ = NULL;
   inputMeans_ = NULL;
   inputStdevs_ = NULL;
   inputNames_ = NULL;
   outputNames_ = NULL;
   outputLevel_ = 0;
   nSamples_ = 0;
   nInputs_ = 0;
   nOutputs_ = 0;
   psuadeIO_ = NULL;
}

// ************************************************************************
// Copy constructor by Bill Oliver 
// ------------------------------------------------------------------------
PsuadeSession::PsuadeSession(const PsuadeSession &ps)
{
   int ii;
   outputLevel_ = ps.outputLevel_;
   nSamples_ = ps.nSamples_;
   nInputs_ = ps.nInputs_;
   nOutputs_ = ps.nOutputs_;
   rsType_ = ps.rsType_;
   owned_ = 1;
   if (nInputs_ > 0 && ps.inputLBounds_ != NULL)
   {
      inputLBounds_ = new double[nInputs_];
      inputUBounds_ = new double[nInputs_];
      for (ii = 0; ii < nInputs_; ii++)
      {
         inputLBounds_[ii] = ps.inputLBounds_[ii];
         inputUBounds_[ii] = ps.inputUBounds_[ii];
      }
   }
   if (nSamples_ > 0 && nInputs_ > 0)
   {
      if (ps.sampleInputs_ != NULL)
      {
         sampleInputs_ = new double[nSamples_*nInputs_];
         for (ii = 0; ii < nSamples_*nInputs_; ii++)
            sampleInputs_[ii] = ps.sampleInputs_[ii];
      }
      if (ps.inputNames_ != NULL)
      {
         inputNames_ = new char*[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
         {
            inputNames_[ii] = new char[1001];
            strcpy(inputNames_[ii], ps.inputNames_[ii]);
         }
      }
      if (ps.inputPDFs_ != NULL)
      {
         inputPDFs_ = new int[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
            inputPDFs_[ii] = ps.inputPDFs_[ii];
      }
      if (ps.inputMeans_ != NULL)
      {
         inputMeans_ = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
            inputMeans_[ii] = ps.inputMeans_[ii];
      }
      if (ps.inputStdevs_ != NULL)
      {
         inputStdevs_ = new double[nInputs_];
         for (ii = 0; ii < nInputs_; ii++)
            inputStdevs_[ii] = ps.inputStdevs_[ii];
      }
      corMatrix_.load((psMatrix &) ps.corMatrix_); 
   }
   if (nSamples_ > 0 && nOutputs_ > 0 & ps.sampleOutputs_ != NULL)
   {
      sampleOutputs_ = new double[nSamples_*nOutputs_];
      for(ii = 0; ii < nSamples_*nOutputs_; ii++)
         sampleOutputs_[ii] = ps.sampleOutputs_[ii];
      if (ps.outputNames_ != NULL)
      {
         outputNames_ = new char*[nOutputs_];
         for(ii = 0; ii < nOutputs_; ii++)
         {
            outputNames_[ii] = new char[1001];
            strcpy(outputNames_[ii], ps.outputNames_[ii]);
         }
      }
   }
   if (nSamples_ > 0 && nOutputs_ > 0 & ps.sampleStates_ != NULL)
   {
      sampleStates_ = new int[nSamples_];
      for(ii = 0; ii < nSamples_; ii++)
         sampleStates_[ii] = ps.sampleStates_[ii];
   }
}

// ************************************************************************
// operator= 
// ------------------------------------------------------------------------
PsuadeSession & PsuadeSession::operator=(const PsuadeSession & ps)
{
   int ii;
   if(this == &ps)  return *this;
   owned_ = 1;
   inputLBounds_ = NULL;
   inputUBounds_ = NULL;
   sampleInputs_ = NULL;
   sampleOutputs_ = NULL;
   sampleStates_ = NULL;
   inputPDFs_ = NULL;
   inputMeans_ = NULL;
   inputStdevs_ = NULL;
   inputNames_ = NULL;
   outputNames_ = NULL;
   psuadeIO_ = NULL;
   outputLevel_ = ps.outputLevel_;
   nSamples_ = ps.nSamples_;
   nInputs_ = ps.nInputs_;
   nOutputs_ = ps.nOutputs_;
   rsType_ = ps.rsType_;
   if (nInputs_ > 0 && ps.inputLBounds_ != NULL)
   {
      inputLBounds_ = new double[nInputs_];
      inputUBounds_ = new double[nInputs_];
      for(int ii = 0; ii < nInputs_; ii++)
      {
         inputLBounds_[ii] = ps.inputLBounds_[ii];
         inputUBounds_[ii] = ps.inputUBounds_[ii];
      }
   }
   if (nSamples_ > 0 && nInputs_ > 0)
   {
      if (ps.sampleInputs_ != NULL)
      {
         sampleInputs_ = new double[nSamples_*nInputs_];
         for(ii = 0; ii < nSamples_*nInputs_; ii++)
            sampleInputs_[ii] = ps.sampleInputs_[ii];
      }
      if (ps.inputNames_ != NULL)
      {
         inputNames_ = new char*[nInputs_];
         for(ii = 0; ii < nInputs_; ii++)
         {
            inputNames_[ii] = new char[1001];
            strcpy(inputNames_[ii], ps.inputNames_[ii]);
         }
      }
   }
   if (nSamples_ > 0 && nOutputs_ > 0 & ps.sampleOutputs_ != NULL)
   {
      sampleOutputs_ = new double[nSamples_*nOutputs_];
      for(ii = 0; ii < nSamples_*nOutputs_; ii++)
         sampleOutputs_[ii] = ps.sampleOutputs_[ii];
      if (ps.outputNames_ != NULL)
      {
         outputNames_ = new char*[nOutputs_];
         for(ii = 0; ii < nOutputs_; ii++)
         {
            outputNames_[ii] = new char[1001];
            strcpy(outputNames_[ii], ps.outputNames_[ii]);
         }
      }
   }
   if (nSamples_ > 0 && nOutputs_ > 0 & ps.sampleStates_ != NULL)
   {
      sampleStates_ = new int[nSamples_];
      for(ii = 0; ii < nSamples_; ii++)
         sampleStates_[ii] = ps.sampleStates_[ii];
   }
   return *this;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
PsuadeSession::~PsuadeSession()
{
   int ii;
   if (owned_)
   {
      if (sampleInputs_ != NULL) delete [] sampleInputs_;
      if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
      if (sampleStates_ != NULL) delete [] sampleStates_;
      if (inputLBounds_ != NULL) delete [] inputLBounds_;
      if (inputUBounds_ != NULL) delete [] inputUBounds_;
      if (inputNames_ != NULL)
      {
         for(ii = 0; ii < nInputs_; ii++)
            if (inputNames_[ii] != NULL) delete [] inputNames_[ii];
         delete [] inputNames_[ii];
      }
      if (outputNames_ != NULL)
      {
         for(ii = 0; ii < nOutputs_; ii++)
            if (outputNames_[ii] != NULL) delete [] outputNames_[ii];
         delete [] outputNames_[ii];
      }
      if (inputPDFs_ != NULL) delete [] inputPDFs_;
      if (inputMeans_ != NULL) delete [] inputMeans_;
      if (inputStdevs_ != NULL) delete [] inputStdevs_;
   }
   psuadeIO_ = NULL;
}

