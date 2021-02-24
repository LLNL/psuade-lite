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
// Definition for the PsuadeBase class
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __PSUADEBASEH__
#define __PSUADEBASEH__

#include "Matrix.h"
#include "PsuadeData.h"
#include "FunctionInterface.h"
#include "Sampling.h"
#include "AnalysisManager.h"
#include "Exceptions.h"
#include "PsuadeSession.h"

// ************************************************************************
// class definition 
// ************************************************************************
class PsuadeBase 
{
   int    outputLevel_;
   int    useRSModel_;
   int    plotFlag_;
   int    plotInput1_;
   int    plotInput2_;
   int    optimalCount_;
   double *optimalXData_;
   double *optimalYData_;
   double *optInitXData_;
   double *optInitYData_;
   Sampling    *sampler_;

public:
   double *sampleInputs_;
   double *sampleOutputs_;
   int    *sampleStates_;
   double *iLowerB_;
   double *iUpperB_;
   int    *inputPDFs_;
   double *inputMeans_;
   double *inputStds_;
   psMatrix *inputCMat_;
   char   **inputNames_;
   char   **outputNames_;
   int    *tagArray_;
   double *dataReg_;
   int    nInputs_;
   int    nSamples_;
   int    nOutputs_;

public:
   PsuadeData  *psuadeIO_;
   AnalysisManager analysisManager_;
   int jobsCompleted;
   int psuade_stop;

#ifdef HAVE_PYTHON
   PyObject* update_gui;
   PyObject* yesno_dialog;
#endif

public:

   // Constructor 
   PsuadeBase();

   // Destructor 
   ~PsuadeBase();

   // get input data from a given file
   // @param filename - name of file from which data are to be read
   // @param psuadeio_filename - optional - file to read psuade_io data from
   int  getInputFromFile(const char *filename,const char *psuadeio_filename=NULL);

   // run the experiments
   int  run() throw(Psuade_Stop_Exception);

   // run in interactive mode
   int  sessionInteractive();

   // run in parallel interactive mode
   int  sessionInteractiveParallel();
   //
   // assign operator 
   PsuadeBase& operator=(const PsuadeBase &);

private:

   int    runUniform();
   int    runEnsemble();
   int    runAdaptiveNN();
   int    runAdaptiveErrBased1();
   int    runAdaptiveErrBasedG();
   int    runAdaptivePRA();
   int    runAdaptiveOpt();
   int    runAdaptiveGradBased();
   void   pgPlotResponseSurface();
   void   cleanUp();
   int    interpretInteractive();
   int    interpretInteractiveParallel();
   int    RSBasedAnalysis(char *, PsuadeSession *);
   int    setupAnalysis(int, int);
   int    deleteAnalysis(int);
   void   *genRSModel(char *fname, int numInputs, int outputID);
   int    fillInSamples(int, int, double *, int *);
   int    setupGuide();
   int    genDriver(int);
   int    genBatchFile(int);
   int    genWorkFlow();
   int    genSetup(int, char *);
   int    evaluateFull(int, double *, FunctionInterface *, double *,
                       int, double *, double *, double *);
};

#endif // __PSUADEBASEH__

