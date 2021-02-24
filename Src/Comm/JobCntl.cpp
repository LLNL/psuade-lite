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
// Functions for the Job Control class
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#ifdef WINDOWS
#define UNICODE
#include <windows.h>
//extern void Sleep(unsigned long milliseconds);
#endif //WINDOWS

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "JobCntl.h"

#ifdef HAVE_PYTHON
#include "Python.h"
#endif


// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
JobCntl::JobCntl()
{
   maxPJobs_ = 1;
   jobWaitTime_ = 0;
   funcIO_ = NULL;
   sampleInputs_ = NULL;
   sampleOutputs_ = NULL;
   sampleStates_ = NULL;
   nSamples_ = 0;
   nInputs_ = 0;
   nOutputs_ = 0;
   outputLevel_ = 0;
#ifdef HAVE_PYTHON
   update_gui = NULL;
#endif
}

// ************************************************************************
// Destructor
// ------------------------------------------------------------------------
JobCntl::~JobCntl()
{
   if (sampleInputs_  != NULL) delete [] sampleInputs_;
   if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
   if (sampleStates_  != NULL) delete [] sampleStates_;
}

// ************************************************************************
// load the function interface
// ------------------------------------------------------------------------
int JobCntl::loadFuncIO(FunctionInterface *funcIO)
{
   funcIO_ = funcIO;
   return 0;
}

// ************************************************************************
// set output level
// ------------------------------------------------------------------------
int JobCntl::setOutputLevel(int level)
{
   outputLevel_ = level;
   if (outputLevel_ < 0) outputLevel_ = 0;
   return 0;
}

// ************************************************************************
// set job wait time
// ------------------------------------------------------------------------
int JobCntl::setJobWaitTime(int wtime)
{
   jobWaitTime_ = wtime;
   if (jobWaitTime_ < 0) jobWaitTime_ = 0;
   if (jobWaitTime_ > 100000) jobWaitTime_ = 100000;
   return 0;
}

// ************************************************************************
// set GUI object
// ------------------------------------------------------------------------
int JobCntl::setGUIObject(void *obj)
{
#ifdef HAVE_PYTHON
   update_gui = (PyObject *) obj;
#endif
   return 0;
}

// ************************************************************************
// set maximum number of parallel jobs
// ------------------------------------------------------------------------
int JobCntl::setMaxParallelJobs(int pjobs)
{
   maxPJobs_ = pjobs;
   if (maxPJobs_ < 0) maxPJobs_ = 0;
   if (maxPJobs_ > 100000) maxPJobs_ = 100000;
   return 0;
}

// ************************************************************************
// load the sample data
// ------------------------------------------------------------------------
int JobCntl::loadSample(int nSamples, int nInputs, int nOutputs,
                        double *sampleInputs, double *sampleOutputs, 
                        int *sampleStates)
{
   int ii;

   if (sampleInputs_  != NULL) delete [] sampleInputs_;
   if (sampleOutputs_ != NULL) delete [] sampleOutputs_;
   if (sampleStates_  != NULL) delete [] sampleStates_;
   if (nSamples <= 0 || nInputs <= 0 || nOutputs <= 0)
   {
      printf("JobCntl ERROR: invalid parameters.\n");
      exit(1);
   }
   if (sampleInputs == NULL || sampleOutputs == NULL || sampleStates == NULL)
   {
      printf("JobCntl ERROR: invalid parameters.\n");
      exit(1);
   }
   nSamples_ = nSamples;
   nInputs_ = nInputs;
   nOutputs_ = nOutputs;
   sampleInputs_ = new double[nSamples*nInputs];
   sampleOutputs_ = new double[nSamples*nOutputs];
   sampleStates_ = new int[nSamples];
   for (ii = 0; ii < nSamples*nInputs; ii++)
      sampleInputs_[ii] = sampleInputs[ii];
   for (ii = 0; ii < nSamples*nOutputs; ii++)
      sampleOutputs_[ii] = sampleOutputs[ii];
   for (ii = 0; ii < nSamples; ii++) sampleStates_[ii] = sampleStates[ii];
   return 0;
}

// ************************************************************************
// run the sample
// ------------------------------------------------------------------------
int JobCntl::execute()
{
   int   ss, runJobs, pJobCnt, jobsCompleted;
   int   status, ii;

   if (nSamples_ <= 0 || nInputs_ <= 0 || nOutputs_ <= 0)
   {
      printf("JobCntl execute ERROR: scalar parameters not set.\n");
      exit(1);
   }
   if (sampleInputs_ == NULL || sampleOutputs_ == NULL || sampleStates_ == NULL)
   {
      printf("JobCntl execute ERROR: pointer parameters not set.\n");
      exit(1);
   }
   if (funcIO_ == NULL)
   {
      printf("JobCntl execute ERROR: funcIO not set.\n");
      exit(1);
   }

   runJobs = 0;
   for (ss = 0; ss < nSamples_; ss++) if (sampleStates_[ss] == 0) runJobs++;


   pJobCnt = 0;
   jobsCompleted = 0;
   while (jobsCompleted < runJobs)
   {
      for (ss = 0; ss < nSamples_; ss++)
      {
#ifdef HAVE_PYTHON
         PyObject* temp;
         if (update_gui != NULL) {
            temp = PyObject_CallObject(update_gui, NULL);
            if (temp != NULL) Py_DECREF(temp);
	 }
#endif
         if ((sampleStates_[ss] == 0) && (pJobCnt < maxPJobs_))
         {

            status = funcIO_->evaluate(ss,nInputs_,
                         &sampleInputs_[ss*nInputs_], nOutputs_, 
                         &sampleOutputs_[ss*nOutputs_],0);   

            if (outputLevel_ > 3)
            {
               printf("Completed: job = %6d\n", ss+1);
               for (ii = 0; ii < nInputs_; ii++)
                  printf("\t\tinput data %3d = %24.16e\n", ii+1,
                         sampleInputs_[ss*nInputs_+ii]);
               printf("\t\toutput data     = %24.16e\n",sampleOutputs_[ss]);
            }


            if (status == 0) 
            {
               jobsCompleted++;
               if (outputLevel_ > 0)
               {
                  if ((jobsCompleted % 11) == 1)
                  {
                     printf(".");
                     fflush(stdout);
                  }
               }
            }
            else
            {
               sampleStates_[ss] = status;
               pJobCnt++;
            }
         }
         else if (sampleStates_[ss] >= 2)
         {

            status = funcIO_->evaluate(ss,nInputs_,
                           &sampleInputs_[ss*nInputs_], nOutputs_, 
                           &sampleOutputs_[ss*nOutputs_],2);   
            if (status == 0) 
            {
               jobsCompleted++;
               if (outputLevel_ > 3)
               {
                  printf("Completed: job = %6d\n", ss+1);
                  for (ii = 0; ii < nInputs_; ii++)
                     printf("\t\tinput data %3d = %24.16e\n", ii+1,
                            sampleInputs_[ss*nInputs_+ii]);
                  printf("\t\toutput data     = %24.16e\n",sampleOutputs_[ss]);
               }
               pJobCnt--;
            }
            else sampleStates_[ss]++;
         }
      }


      if ((jobsCompleted < runJobs) && (jobWaitTime_ > 0))
      {
#ifdef WINDOWS
         Sleep(1000 * jobWaitTime_);
#else
         sleep(jobWaitTime_);
#endif

      }
      if (outputLevel_ > 0) 
         printf("\nPSUADE: jobs completed = %d(of %d)\n",
                jobsCompleted, runJobs);
   }
   return 0;
}

// ************************************************************************
// get sample outputs
// ------------------------------------------------------------------------
int JobCntl::getSampleOutputs(int nSamples, int nOutputs, double *outputs)
{
   int ii;

   if (nSamples_ <= 0 || nInputs_ <= 0 || nOutputs_ <= 0)
   {
      printf("JobCntl getSampleOutputs ERROR: scalar parameters not set.\n");
      exit(1);
   }
   if (sampleInputs_ == NULL || sampleOutputs_ == NULL || sampleStates_ == NULL)
   {
      printf("JobCntl getSampleOutputs ERROR: vector parameters not set.\n");
      exit(1);
   }
   if (nSamples_ != nSamples)
   {
      printf("JobCntl getSampleOutputs ERROR: nSamples mismatch.\n");
      exit(1);
   }
   if (nOutputs_ != nOutputs)
   {
      printf("JobCntl getSampleOutputs ERROR: nOutputs mismatch.\n");
      exit(1);
   }

   for (ii = 0; ii < nSamples*nOutputs; ii++) outputs[ii] = sampleOutputs_[ii];

   return 0;
}

