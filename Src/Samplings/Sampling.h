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
// Definition for the abstract class Sampling
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __BASESAMPLINGH__
#define __BASESAMPLINGH__

#include <string.h>
#include "PsuadeData.h"
#include "Vector.h"
#include "Matrix.h"

/**
 * @name Sampling class (abstract)
 *
 * This is the base class on which all sampling methods are built.
 *
 **/
/*@{*/

// ************************************************************************
// class definition
// ************************************************************************
class Sampling 
{
public:

  int printLevel_;
  int samplingID_;
  int nSamples_;
  int nInputs_;
  int nOutputs_;
  int randomize_;
  int nReplications_;
  psVector  vecLBs_;
  psVector  vecUBs_;
  psVector  vecSamInps_;
  psVector  vecSamOuts_;
  psIVector vecSamStas_;
  psIVector vecSymTable_;
  psIVector vecInpNumSettings_;
  psVector  *vecInpSettings_;
  int maxNumSettings_;

  /** Constructor
   */
  Sampling();

  /** Copy Constructor created by Bill Oliver
   */
  Sampling(const Sampling &);

  /** Overload operator = function by Bill Oliver
   */
  Sampling & operator=(const Sampling &);

  /** Destructor
   */
  virtual ~Sampling();

  /** This function sets the input lower and upper bounds.
   *  @param printLevel : diagnostics level
   */
  int setPrintLevel(int printLevel);

  /** This function sets the input lower and upper bounds.
   *  @param nInputs : number of inputs
   *  @param iLowerB : input lower bounds
   *  @param iUpperB : input upper bounds
   */
  int setInputBounds(int nInputs, double *iLowerB, double *iUpperB);

  /** This function sets sampling parameters
   *  @param nSamples : sample size
   *  @param nReps : number of replications 
   *  @param randomize : randomize flags
   */
  int setSamplingParams(int nSamples, int nReps, int randomize);

  /** This function sets different input settings (user specified).
   *  @param nInputs : number of inputs
   *  @param nSettings : number of settings for each input
   *  @param settings : the settings for each input
   *  @param symTable : symbol table
   */
  virtual int setInputParams(int nInputs, int *nSettings, 
                             double **settings, int *symTable);

  /** This function sets output parameters
   *  @param nOutputs : number of outputs
   */
  virtual int setOutputParams(int nOutputs);

  /** This function loads externally created sample points.
   *  @param nInputs  : number of inputs
   *  @param nOutputs : number of outputs
   *  @param nSamples : number of samples
   *  @param sampleInputs : sample inputs
   *  @param sampleOutputs : sample outputs
   *  @param sampleStates : sample states
   */
  int loadSamples(int nInputs, int nOutputs, int nSamples,
                  double *sampleInputs, double *sampleOutputs, 
                  int *sampleStates);

  /** This function gets the number of sample points.
   */
  int getNumSamples();

  /** This function gets the number of input parameters.
   */
  int getNumInputs();

  /** This function gets the number of output parameters.
   */
  int getNumOutputs();

  /** This function gets the sample points.
   *  @param nInputs  : number of inputs
   *  @param nOutputs : number of outputs
   *  @param nSamples : number of samples
   *  @param sampleInputs : sample inputs
   *  @param sampleOutputs : sample outputs
   *  @param sampleStates : sample states
   */
  int getSamples(int nInputs, int nOutputs, int nSamples, 
                 double *sampleInputs, double *sampleOutputs, 
                 int *sampleStates);

  /** This function gets the sample points.
   *  @param nInputs  : number of inputs
   *  @param nOutputs : number of outputs
   *  @param nSamples : number of samples
   *  @param sampleInputs : sample inputs
   *  @param sampleOutputs : sample outputs
   */
  int getSamples(int nInputs, int nOutputs, int nSamples, 
                 double *sampleInputs, double *sampleOutputs); 

  /** This function stores the given sample outputs.
   *  @param nOutputs : number of outputs
   *  @param nSamples : number of samples
   *  @param sampleOutputs : sample outputs
   *  @param sampleStates : sample states
   */
  int storeResult(int nOuptuts, int nSamples, double *sampleOutputs, 
                  int *sampleStates);


  /** This function fully constructs a sampling from the psuadeIO_
   *  input.  All the common functionality it done in the Sampleing
   *  class, subclasses should overwrite this function.
   *  @param psuadeIO_ : A Psuade Data file, should have Input/Output 
   *  variables and Method defined.
   **/
  virtual int doSampling(PsuadeData* psuadeIO_);

  /** This function initializes the sample data.
   *  @param initLevel : 0 if we need to generate samples, 
   *                     1 if we will be loading external samples
   */
  virtual int initialize(int initLevel)=0;

  /** This function performs sample refinement.
   *  @param nLevels : degree of refinements
   *  @param randomize : randomize flag
   *  @param threshold : used to identify where to refine
   *  @param nSamples : sample size
   *  @param sampleErrors : errors for each one of the previous sample point
   */
  virtual int refine(int nLevels, int randomize, double threshold,
                     int nSamples, double *sampleErrors);

  /** This function performs sample repair (for MOAT mostly).
   *  @param fname : restart file name 
   *  @param start : start sample point 
   */
  virtual int repair(char *, int start);

  /** This function allocates storage for sample inputs and outputs.
   */
  int allocSampleData();

  /** This function sets internal sampling-specific parameters
      @param command : command string
      @param paramStrings : pointers to parameters
   */
  virtual int setParam(char *sparam);
};

/*@}*/

// ************************************************************************
// These functions are used to create and destroy sampling objects.
// ************************************************************************
extern Sampling *SamplingCreateFromID(int);
extern int SamplingDestroy(Sampling *);
extern int SamplingQuality(int,int,double *,double *,double *,double *);
extern int CheckSampleSmoothness(psVector,psVector,psVector,psVector,
                                 int, psVector &, psVector &, int);
extern int CheckSampleSmoothness2(psVector, psVector, psVector, psVector, 
               psVector, psVector, int, psVector &, psVector &, int);

extern int getSamplingMethod(char *pString);

#endif // __BASESAMPLINGH__

