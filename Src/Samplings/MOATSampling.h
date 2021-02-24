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
// Definition for the Morris one-at-a-time sampling class
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#ifndef __MOATSAMPLINGH_
#define __MOATSAMPLINGH_

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Sampling.h"
#include "PsuadeData.h"

/**
 * @name Morris one-at-a-time samping method
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
class MOATSampling: public Sampling 
{
  int P_;
  psIVector vecInpSubset_;

public:

  /** constructor */
  MOATSampling();

  // Bill Oliver added a Copy Constructor
  MOATSampling(const MOATSampling &);

   /** destructor */
  ~MOATSampling();

  /** initialization 
      @param flag: flag to signal how far to initialize
   */
  int initialize(int flag);

  /** This function refines an incoming sample
      @param ratio: refinement ratio
      @param randomize: generate randomized sample
      @param thresh: threshold
      @param nSamples: sample size
      @param samErrs: errors for each sample point
   */
  int refine(int ratio,int randomize,double thresh,int nSamples,double *samErrs);

  /** This function overloads the assignment operator
      @param obj : Sampling object
   */
  MOATSampling& operator=(const MOATSampling &);

  /** This function modifies the current sample
      @param fname : name of file containing repair data
      @param start : the first sample location for repair
   */
  int repair(char *fname, int start);

  /** This function generates samples with constraints
      @param nInputs : number of inputs
      @param lbounds : input lower bounds
      @param ubounds : input upper bounds
   */
  int genRepair(int nInputs, double *lower, double *upper);

  /** This function generates samples with constraints
      @param psuadeIO : container for all available data
   */
  int genRepair(PsuadeData *);

private:

  int initializeHighDimension();
  int generate(double **);
  int merge();
  int checkSample(int, int, double *);
  int checkSample2(int, int, double *);
};

#endif // __MOATSAMPLINGH_

