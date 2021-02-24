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
// Definition for the class AnalysisManager
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#ifndef __ANALYSISMANAGERH__
#define __ANALYSISMANAGERH__

#ifdef HAVE_PYTHON
#include <Python.h>
#endif

#include "Sampling.h"
#include "Analyzer.h"
#include "PsuadeData.h"
#include "Vector.h"

#include "MOATAnalyzer.h"

// ************************************************************************
// class definition
// ************************************************************************

class AnalysisManager 
{
   int        numAnalyzers_;
   Analyzer   **analyzers_;
   Sampling   *sampler_;
   int        *logXsformFlags_;
   int        sizelogX;
   double     *analysisSampleErrors_;

public:

   AnalysisManager();
   ~AnalysisManager();

   int setSampler(Sampling *sampler);
   
   /**
    * Create the analyzers and do any necessary setup using the parameters 
    * from a PsuadeData object. Internally it just calls 
    * setup(int anaMethod, int faType)
    * @param psuadeIO: The PsuadeData object that holds all the information 
    * for setting up the analyzers
    */
   int setup(PsuadeData *psuadeIO); 

   /**
    * Create the analyzers and do any necessary setup.
    * @param anaMethod: anaMethod is a set of flags compressed into a int 
    *                   using &. An analyzer is turned on or off by a single 
    *                   bit in the int. see sysdef.h
    * @param faType:    Is only used if at least one of the used analyzers 
    *                   uses a response surface.  faType Defines the response 
    *                   surface type to use. 
    */
   int setup(int anaMethod, int faType);
   int clearAnalyzers();

   /**
    * analyze is run after setup.  analyze runs all the analyzers that were setup.
    * @param psuadeIO: The PsuadeData object that holds all the information 
    *                  for setting up all analyzers
    * @param nLevels: Only used when doing adaptive sampling/refinement.  
    *                 nLevels gives the ID of the the current refinment process.  
    *                 It varies from 0 to Total # of Refinements + 1!
    * @param levelSeps: An array of int, only used when doing adaptive 
    *                   sampling/refinement. levelSeps is the number of samples 
    *                   refined at each step.  It is indexed by nLevels.
    * @param analysisOutputID: Analysis is performed on a single output.  
    *                   analysisOutputID is the ID of that output.  It can take 
    *                   any value 1-numOutputs inclusive. 
    */
   int analyze(PsuadeData *psuadeIO, int nLevels, int *levelSeps, int analysisOutputID);


   /**
    * Set log transform flags (use logarithms of the inputs/outputs)
    * @param n: The size of the flags array
    * @param flags: The flags to set.  (Jim- I don't know anything about test flags.)
    */
   int loadLogXsformFlags(int n, int *flags);

   
   /**
    * This function performs analysis without a PsuadeData object, but 
    *               only with MomentAnalyzer.
    * 
    * @param anaMethod: Completely ignored.
    * @param nSamples:  The number of sample points 
    * @param vLower:    The lower bound for the inputs (the length is also take 
    *                      as the number of inputs)
    * @param vUpper:    The upper bound for the inputs 
    * @param vInputs:   The Inputs to the simulation (length = nSamples * nInputs) 
    * @param vOutputs:  The Outputs of the output to be anlyzed (Analysis only works 
    *                   on a single output, so length must = nSamples
    * @param outputLevel: The printLevel
    */
   int analyze(int anaMethod, int nSamples, psVector &vLower, 
               psVector &vUpper, psVector &vInputs, 
               psVector &vOutputs, int outputLevel);


   /**
    * analyze OneSampleTest and TwoSampleTest only
    * Only used interactively, the user must input the other parameters at the terminal
    * @param anaMethod: PSUADE_ANA_1SAMPLE or PSUADE_ANA_2SAMPLE (Can't 
    *                   do both in the same call either)
    */
   int analyze(int anaMethod);

   /**
    * get sample error bars for each individual sample point
    * Only two analyzers actually seem to produce this: OneSigmaAnalyzer and 
    *               RSFuncApproxAnalyzer -Jim
    * Return an array of doubles nSamples long.  Each double represents the 
    *                  error bars for that sample.
    */
   double *getSampleErrors();

   /**
    * special request (send specialized parameters to individual methods, 
    *                 analyzer specific)
    * @param anaMethod: The analyzer we're adding special paramters for 
    *                    (you must have already passed it to setup)
    * @param narg: The number of arguments being passed in (argc)
    * @param argv: An array of argument strings
    */
   int specialRequest(int anaMethod, int narg, char **argv);

   AnalysisManager& operator=(const AnalysisManager &);

   MOATAnalyzer *getMOATAnalyzer();
#ifdef HAVE_PYTHON
   PyObject *AnalysisDataList;
#endif
};

#endif // __ANALYSISMANAGERH__

