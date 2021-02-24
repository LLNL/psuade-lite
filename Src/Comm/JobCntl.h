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
// Definition for the Job Control class 
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#ifndef __JOBCNTL__
#define __JOBCNTL__

#ifdef HAVE_PYTHON
#include <Python.h>
#endif
#include "FunctionInterface.h"

/**
 * @name Job control
 *
 **/

// ************************************************************************
// class definition
// ************************************************************************
                                                                                
class JobCntl
{
   double *sampleInputs_;
   double *sampleOutputs_;
   int    *sampleStates_;
   int    nSamples_;
   int    nInputs_;
   int    nOutputs_;
   int    maxPJobs_;
   int    jobWaitTime_;
   int    outputLevel_;
   FunctionInterface *funcIO_;

#ifdef HAVE_PYTHON
   PyObject* update_gui;
#endif

public:

   JobCntl();
   ~JobCntl();

   int setMaxParallelJobs(int);
   int setJobWaitTime(int);
   int setOutputLevel(int);
   int setGUIObject(void *);
   int loadFuncIO(FunctionInterface *);
   int loadSample(int,int,int,double*,double*,int*);
   int execute();
   int getSampleOutputs(int,int,double*);
};

#endif // __JOBCNTL__

