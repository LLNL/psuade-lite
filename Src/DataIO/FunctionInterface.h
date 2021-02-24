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
// Definition for the class FunctionInterface
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __FUNCTIONINTERFACEH__
#define __FUNCTIONINTERFACEH__

#include <math.h>
//#include <stream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sysdef.h"
#include "PsuadeData.h"
#include "FuncApprox.h"

// ************************************************************************
// Class Definition 
// ************************************************************************
class FunctionInterface 
{
   int    nInputs_;
   int    nOutputs_;
   char   **inputNames_;
   char   **outputNames_;
   char   appDriver_[2000];
   char   optDriver_[2000];
   char   auxOptDriver_[2000];
   char   ensembleDriver_[2000];
   char   ensembleOptDriver_[2000];
   char   appInputTemplate_[2000];
   char   appOutputTemplate_[2000];
   int    executionMode_;
   int    launchInterval_;
   int    appOptFlag_;
   int    useRSModel_;
   FuncApprox **rsPtrs_;
   int    printLevel_;
   int    rsRanFlag_;
   int    whichLocalFunction_;

public:
   int    *rsIndices_;
   double *rsValues_;
   int    rsnInps_;

public:

   /** constructor */
   FunctionInterface();

   /** Copy constructor by Bill Oliver */
   FunctionInterface(const FunctionInterface & fi);

   /** destructor */
   ~FunctionInterface();

   /** This function loads input data 
       @param nInputs - number of inputs
       @param iNames - names of each input
    */
   int  loadInputData(int nInputs, char **iNames);

   /** This function loads output data 
       @param nOutputs - number of outputs
       @param oNames - names of each output
    */
   int  loadOutputData(int nOutputs, char **oNames);

   /** This function loads output data 
       @param nFunc - number of function names
       @param fNames - names of each function
    */
   int  loadFunctionData(int nFunc, char **fNames);

   /** This function sets the run mode to be asynchronous */
   int  setAsynchronousMode();

   /** This function sets stochastic response surface evaluation */
   int  setStochasticRSMode(int);

   /** This function sets the run mode to be synchronous */
   int  setSynchronousMode();

   /** This function sets print level */
   int  setOutputLevel(int);

   /** This function sets interval between launching each run */
   int  setLaunchInterval(int);

   /** This function sets the simulator number
       @param num - simulator number
     */
   int  setDriver(int num);

   /** This function sets which local function to use
       @param option - local function identifier 
     */
   int  setLocalFunction(int num);

   /** This function calls the simulation runs */
   int  evaluate(int,int,double *,int,double *,int);

   /** This function calls ensemble simulation runs */
   int ensembleEvaluate(int,int,double *,int, double *, int);

   /** This function extracts the number of inputs */
   int  getNumInputs();

   /** This function extracts the number of outputs */
   int  getNumOutputs();

   /** This function gets the simulator number */
   int  getDriver();

   /** This function extracts the names of the inputs */
   char **getInputNames();

   /** This function extracts the names of the outputs */
   char **getOutputNames();

   /** This function extracts the name of the application driver */
   char *getApplicationDriver();

   /** This function extracts the name of the optimization driver */
   char *getOptimizationDriver();

   /** This function extracts the name of the auxiliary optimization driver */
   char *getAuxOptimizationDriver();

   /** This function overrides the assign operator */
   FunctionInterface& operator=(const FunctionInterface &);

   /** This function evaluates a test problem or a method */
   int psLocalFunction(int, double *, int, double *);

private :
   int psEnsembleLocalFunction(int, int, double *, int, double *);
   void *genRSModel(char *);
};

// ************************************************************************
// function to instantiate this class
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterface(PsuadeData *psuadeIO);

//** Blows away the appDriver for some reason. That's simpler?*/
FunctionInterface *createFunctionInterfaceSimplified(PsuadeData *psuadeIO);
FunctionInterface *createFunctionInterfaceGivenAppDriver(int, int, char *fname);

#endif //  __FUNCTIONINTERFACEH__

