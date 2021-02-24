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
// Definition for the sampling class based on the Metis partitioning
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#ifndef __METISSAMPLINGH__
#define __METISSAMPLINGH__

#include <string.h>
#include "Sampling.h"

/**
 * @name Metis (domain decomposition) samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class MetisSampling : public Sampling 
{
  int refineType_;
  int refineSize_;
  int n1d_;
  int nAggrs_;
  psIVector vecAggrCnts_;
  psIVector *vecAggrLabels_;
  int graphN_;
  psIVector vecGraphI_;
  psIVector vecGraphJ_;
  psIVector vecCellsOccupied_;;
  int changeInfoName_;

public:

  /** constructor */
  MetisSampling();

  // Copy Constructor added by Bill Oliver
  MetisSampling(const MetisSampling & ms);

  /** destructor */
  ~MetisSampling();

  /** initialization 
      @param flag: flag to signal how far to initialize
   */
  int initialize(int flag);

  /** This function refines an incoming sample
      @param ratio: refinement ratio
      @param randomize: generate randomized sample
      @param thresh: threshold
      @param nSamples: sample size
      @param sampleErrs: errors for each sample point
   */
  int refine(int ratio,int randomize,double thresh,int nSamples,
             double *sampleErrs);

  /** This function overloads the assignment operator
      @param obj : Sampling object
   */
  MetisSampling& operator=(const MetisSampling &);

  /** This function sets internal parameters
      @param sparam : string containing parameter information
   */
  int setParam(char *sparam);
};

#endif // __METISSAMPLINGH__

