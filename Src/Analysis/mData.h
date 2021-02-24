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
// Definition for the class McmcData (for passing MCMC data)
// AUTHOR : CHARLES TONG
// DATE   : 2019
// ************************************************************************
#ifndef __MCMCDATAH__
#define __MCMCDATAH__

#include "PsuadeData.h"

// ************************************************************************
// class definition
// ************************************************************************
class McmcData 
{
public:

  int    printLevel_;
  //**/ sample size
  int    nSamples_;
  //**/ total number of inputs
  int    nInputs_;
  //**/ number of outputs
  int    nOutputs_;

  //**/ sample inputs lower and upper bounds
  //**/ sample inputs and outputs for response surface construction
  psVector VecLowerB_;
  psVector VecUpperB_;
  psVector VecSamInputs_;
  psVector VecSamOutputs_;

  //**/ number of fixed uncertain inputs
  psIVector VecFUInputs_;
  psMatrix MatFUInpSample_;

  //**/ number of calibration inputs
  psIVector VecCUInputs_;
  //**/ response surface type
  int    faType_;
  //**/ use response surface errors
  int    useRSUncertainties_;
  int    addDiscrepancy_;

  //**/ experimental data
  psMatrix MatExpInputs_;
  psMatrix MatExpMeans_;
  psMatrix MatExpStds_;

  //**/ input probability distribution functions
  psMatrix  MatPriorSample_;
  psMatrix  MatPostSample_;
  psVector  VecPostLikelihoods_;

  //**/ values returned from analyzer
  psVector VecRetValues_;

public:
  //**/ constructor
  McmcData();

  //**/ destructor
  ~McmcData();

  //**/ clean up memory
  void clean();
};

#endif // __MCMCDATAH__

