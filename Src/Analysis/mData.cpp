// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team. 
// All rights reserved.
//
// Please see the COPYRIGHT and LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free 
// Software Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU Lesser
// General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <assert.h>
#include "dtype.h"
#include "sysdef.h"
#include "mData.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------ 
McmcData::McmcData()
{
  clean();
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------ 
McmcData::~McmcData()
{ 
  clean();
}

// ************************************************************************
// clean
// ------------------------------------------------------------------------
void McmcData::clean()
{
  printLevel_ = 0;
  nSamples_ = 0;
  nInputs_ = 0;
  nOutputs_ = 0;
  faType_ = 0;
  useRSUncertainties_ = 0;
  addDiscrepancy_ = 0;
  VecFUInputs_.clean();
  MatFUInpSample_.cleanUp();
  VecCUInputs_.clean();
  VecSamInputs_.clean();
  VecSamOutputs_.clean();
  VecLowerB_.clean();
  VecUpperB_.clean();
  MatExpInputs_.cleanUp();
  MatExpMeans_.cleanUp();
  MatExpStds_.cleanUp();
  MatPriorSample_.cleanUp();
  MatPostSample_.cleanUp();
  VecRetValues_.clean();
}

