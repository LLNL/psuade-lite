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
// Functions for FunctionInterface 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <PsuadeCmakeConfig.h>

#ifdef WINDOWS
#include <windows.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "dtype.h"
#include "PsuadeUtil.h"
#include "FunctionInterface.h"
#include "FuncApprox.h"
#include "pData.h"
#include "PDFBase.h"
#include "PsuadeData.h"
#include "Psuade.h"
#include "ProbMatrix.h"
#include "MCMCAnalyzer.h"
#include "mData.h"
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
FunctionInterface::FunctionInterface()
{ 
  nInputs_ = 0;
  nOutputs_ = 0;
  inputNames_ = NULL;
  outputNames_ = NULL;
  strcpy(appDriver_, "true");
  strcpy(optDriver_, "true");
  strcpy(auxOptDriver_, "true");
  strcpy(ensembleDriver_, "true");
  strcpy(ensembleOptDriver_, "true");
  strcpy(appInputTemplate_, "psuadeApps_ct.in");
  strcpy(appOutputTemplate_, "psuadeApps_ct.out");
  executionMode_ = 0;
  launchInterval_ = 0;
  appOptFlag_ = 0;
  useRSModel_ = 0;
  rsPtrs_ = NULL;
  rsIndices_ = NULL;
  rsValues_ = NULL;
  rsRanFlag_ = 0;
  printLevel_ = 0;
  rsnInps_ = 0;
  whichLocalFunction_ = -1;
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
FunctionInterface::FunctionInterface(const FunctionInterface & fi)
{
  int ii;

  nInputs_ = fi.nInputs_;
  nOutputs_ = fi.nOutputs_;
  strcpy(appDriver_, fi.appDriver_);
  strcpy(optDriver_, fi.optDriver_);
  strcpy(auxOptDriver_, fi.auxOptDriver_);
  strcpy(ensembleDriver_, fi.ensembleDriver_);
  strcpy(ensembleOptDriver_, fi.ensembleOptDriver_);
  strcpy(appInputTemplate_, fi.appInputTemplate_);
  strcpy(appOutputTemplate_, fi.appOutputTemplate_);
  rsRanFlag_ = fi.rsRanFlag_;
  printLevel_ = fi.printLevel_;
  rsnInps_ = fi.rsnInps_;

  inputNames_ = new char*[nInputs_];
  for (ii = 0; ii < nInputs_; ii++)
  {
    inputNames_[ii] = new char[80];
    strcpy(inputNames_[ii], fi.inputNames_[ii]);
  }
 
  outputNames_ = new char *[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++)
  {
    outputNames_[ii] = new char[80];
    strcpy(outputNames_[ii], fi.outputNames_[ii]);
  }

  rsPtrs_ = new FuncApprox*[nOutputs_];
  for (ii = 0; ii < nOutputs_; ii++) rsPtrs_[ii] = fi.rsPtrs_[ii];

  rsIndices_ = new int[rsnInps_];
  rsValues_ = new double[rsnInps_];
  for (ii = 0; ii < rsnInps_; ii++)
  {
    rsIndices_[ii] = fi.rsIndices_[ii];
    rsValues_[ii] = fi.rsValues_[ii];
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FunctionInterface::~FunctionInterface()
{ 
  if (inputNames_ != NULL)
  {
    for (int i = 0; i < nInputs_; i++)
      if (inputNames_[i] != NULL) delete [] inputNames_[i];
    delete [] inputNames_;
  }
  if (outputNames_ != NULL)
  {
    for (int j = 0; j < nOutputs_; j++)
      if (outputNames_[j] != NULL) delete [] outputNames_[j];
    delete [] outputNames_;
  }
  if (rsPtrs_ != NULL)
  {
    for (int k = 0; k < nOutputs_; k++)
      if (rsPtrs_[k] != NULL) delete rsPtrs_[k];
    delete [] rsPtrs_;
  }
  if (rsIndices_ != NULL)
  {
    delete [] rsIndices_;
    rsIndices_ = NULL;
  }
  if (rsValues_ != NULL)
  {
    delete [] rsValues_;
    rsValues_ = NULL;
  }
}

// ************************************************************************
// load input data 
// ------------------------------------------------------------------------
int FunctionInterface::loadInputData(int nInputs, char **names)
{ 
  if (nInputs <= 0)
  {
    printf("FunctionInterface LoadInput ERROR: nInputs <= 0\n");
    return 1; 
  }
  nInputs_ = nInputs;
  inputNames_ = new char*[nInputs_];
  for (int i = 0; i < nInputs; i++)
  {
    inputNames_[i] = new char[80];
    if (names != NULL && names[i] != NULL)
         strcpy(inputNames_[i], names[i]);
    else sprintf(inputNames_[i], "X%d", i);
  }
  return 0;
}

// ************************************************************************
// load output data 
// ------------------------------------------------------------------------
int FunctionInterface::loadOutputData(int nOutputs, char **names)
{ 
  if (nOutputs <= 0)
  {
    printf("FunctionInterface loadOutput ERROR: nOutputs <= 0\n");
    return 1; 
  }
  nOutputs_ = nOutputs;
  outputNames_ = new char*[nOutputs_];
  for (int i = 0; i < nOutputs_; i++)
  {
    outputNames_[i] = new char[80];
    if (names != NULL && names[i] != NULL)
         strcpy(outputNames_[i], names[i]);
    else sprintf(outputNames_[i], "Y%d", i);
  }
  return 0;
}

// ************************************************************************
// load the Function information 
// ------------------------------------------------------------------------
int FunctionInterface::loadFunctionData(int length, char **names)
{ 
  int    ii, kk, nInps, status;
  char   inString[200], fname[200], *cString, pString[1001];
  char   winput1[1001], winput2[1001], winput3[1001];
  double ddata;
  FILE   *fp;
  PsuadeData *psIO;
  pData pPtr;

  if (length < 5)
  {
    printf("FunctionInterface loadFunction ERROR: length < 5\n");
    return 1; 
  }

  if (rsPtrs_ != NULL)
  {
    for (ii = 0; ii < nOutputs_; ii++)
      if (rsPtrs_[ii] != NULL) delete rsPtrs_[ii];
    delete [] rsPtrs_;
    rsPtrs_ = NULL;
  }
  if (rsIndices_ != NULL)
  {
    delete [] rsIndices_;
    rsIndices_ = NULL;
  }
  if (rsValues_ != NULL)
  {
    delete [] rsValues_;
    rsValues_ = NULL;
  }

  strcpy(appDriver_, names[0]);
  if (strcmp(names[1], "NONE")) strcpy(appInputTemplate_, names[1]);
  if (strcmp(names[2], "NONE")) strcpy(appOutputTemplate_, names[2]);
  strcpy(optDriver_, names[3]);
  strcpy(auxOptDriver_, names[4]);
  if (length >= 6) strcpy(ensembleDriver_, names[5]);
  if (length >= 7) strcpy(ensembleOptDriver_, names[6]);

  if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
       strcmp(optDriver_,"PSUADE_LOCAL"))
  {
    useRSModel_ = 2;
    return 0; 
  }
  if (strcmp(appDriver_,"PSUADE_LOCAL") &&
      !strcmp(optDriver_,"PSUADE_LOCAL"))
  {
    useRSModel_ = 3;
    return 0; 
  }
  if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
      !strcmp(optDriver_,"PSUADE_LOCAL"))
  {
    useRSModel_ = 4;
    return 0; 
  }

  useRSModel_ = 0;
  if      (appOptFlag_ == 0) strcpy(fname, appDriver_);
  else if (appOptFlag_ == 1) strcpy(fname, optDriver_);
  if      (appOptFlag_ == 2) strcpy(fname, auxOptDriver_);
  fp = fopen(fname, "r");
  if (fp != NULL)
  {
    fscanf(fp, "%10c", inString);
    if (!strncmp(inString, "PSUADE_IO",9)) useRSModel_ = 1;
    fclose(fp);
    if (useRSModel_ == 1)
    {
      psIO = new PsuadeData();
      psIO->setOutputLevel(0);
      status = psIO->readPsuadeFile(fname);
      if (status != 0)
      {
        printf("ERROR : cannot read file %s in PSUADE format\n",fname);
        exit(1);
      }
      psIO->getParameter("input_ninputs", pPtr);
      int rsModelnInps = pPtr.intData_;
      if (rsModelnInps < nInputs_ || rsModelnInps > nInputs_)
      {
        printAsterisks(PL_INFO,0);
        printf("FunctionInterface: setting up RS driver from %s.\n",
               fname);
        printEquals(PL_INFO,0);
        printf("WARNING: nInputs in RS driver does not match nInputs");
        printf(" in original psuade file.\n");
        printf("   nInputs in original psuade file = %d\n",nInputs_);
        printf("   nInputs in RS driver data file  = %d\n",rsModelnInps);
        psIO->getParameter("ana_rsindexfile", pPtr);
        if (!strcmp(pPtr.strArray_[0], "NONE"))
        {
          psIO->getParameter("input_names", pPtr);
          if (psConfig_ != NULL && pPtr.strArray_ != NULL)
          {
            cString = psConfig_->getParameter("num_fixed");
            if (cString != NULL)
            {
              sscanf(cString, "%s %s %d", winput1, winput2, &nInps);
              if (nInps > 0)
              {
                rsIndices_ = new int[rsModelnInps];
                rsValues_ = new double[rsModelnInps];
                for (ii = 0; ii < rsModelnInps; ii++) 
                  rsIndices_[ii] = ii;
                for (ii = 0; ii < nInps; ii++)
                {
                  sprintf(pString,"fixed-%d", ii+1);
                  cString = psConfig_->getParameter(pString);
                  if (cString != NULL)
                  {
                    sscanf(cString,"%s %s %s %lg",winput1,winput2,
                           winput3, &ddata);
                    for (kk = 0; kk < rsModelnInps; kk++)
                    {
                      if (!strcmp(winput2,pPtr.strArray_[kk]))
                      {
                        rsIndices_[kk] = -1;
                        rsValues_[kk] = ddata;
                        printf("Input %4d fixed at %e\n",kk+1,ddata);
                      }
                    }
                    if (kk == rsModelnInps) break;
                  }
                }
                if (ii == nInps)
                {
                  delete [] rsIndices_;
                  delete [] rsValues_;
                  rsIndices_ = NULL;
                  rsValues_ = NULL;
                }
                else rsnInps_ = rsModelnInps;
              }
            }
          }
          if (rsIndices_ == NULL)
          { 
            printf("ERROR: nInputs mismatch and missing rs_index_file.\n");
            printf("   nInputs in original psuade file = %d\n",nInputs_);
            printf("   nInputs in RS driver data file  = %d\n",
                   rsModelnInps);
            printf("ADVICE: Put rs_index_file in %s or in the",fname);
            printf(" original psuade file.\n");
            exit(1);
          }
        }
        else
        {
          printf("WARNING: rs_index_file found in the RS driver file.\n");
          fp = fopen(pPtr.strArray_[0], "r");
          if (fp == NULL)
          {
            printf("ERROR: missing rs_index_file %s in current folder.\n",
                   pPtr.strArray_[0]);
            exit(1);
          }
          else
          {
            fscanf(fp,"%d", &nInps);
            rsnInps_ = nInps;
            if (rsModelnInps != nInps)
            {
              printf("ERROR: invalid nInputs in rs index file.\n");
              printf("  It has to match nInputs in the RS data file.\n");
              printf("  nInputs read     = %d\n", nInps);
              printf("  nInputs expected = %d\n", rsModelnInps);
              printf("  Data format in rs index file should be: \n");
              printf("  line 1: nInputs in RS driver data file\n");
              printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
              printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
              printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
              printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
              printf("  ...\n");
              exit(1);
            }
            rsIndices_ = new int[nInps];
            rsValues_ = new double[nInps];
            for (ii = 0; ii < nInps; ii++) rsIndices_[ii] = 0;
            for (ii = 0; ii < nInps; ii++)
            {
              fscanf(fp, "%d", &kk);
              if (kk != ii+1)
              {
                printf("ERROR: first index in rs index file = %d.\n",
                       rsIndices_[ii]);
                printf("       Must be equal to %d.\n",ii+1);
                printf("  Data format in rs index file should be: \n");
                printf("  line 1: nInputs in RS driver data file\n");
                printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
                printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
                printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
                printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
                printf("  ...\n");
                exit(1);
              } 
              fscanf(fp, "%d", &rsIndices_[ii]);
              if (rsIndices_[ii] < 0 || rsIndices_[ii] > nInputs_)
              {
                printf("ERROR: input %3d = %d not valid\n",ii+1,
                       rsIndices_[ii]);
                printf("       Need to be between 1 and %d\n",nInputs_);
                exit(1);
              }
              rsIndices_[ii]--;
              fscanf(fp, "%lg", &rsValues_[ii]);
              if (rsIndices_[ii] == -1)
                printf("   RS Input %d inactive, fixed at %16.8e\n",
                       ii+1, rsValues_[ii]);
              else
              {
                printf("   RS Input %d   active, mapped to",ii+1);
                printf(" PSUADE input %d\n", rsIndices_[ii]+1);
              }
            }
            fclose(fp);
          }
        }
      }
      psIO->getParameter("input_names", pPtr);
      if (psConfig_ != NULL && rsIndices_ == NULL && pPtr.strArray_ != NULL)
      {
        cString = psConfig_->getParameter("num_fixed");
        if (cString != NULL)
        {
          sscanf(cString, "%s %s %d", winput1, winput2, &nInps);
          if (nInps > 0)
          {
            rsIndices_ = new int[rsModelnInps];
            rsValues_ = new double[rsModelnInps];
            for (ii = 0; ii < rsModelnInps; ii++) rsIndices_[ii] = ii;
            for (ii = 0; ii < nInps; ii++)
            {
              sprintf(pString,"fixed-%d", ii+1);
              cString = psConfig_->getParameter(pString);
              if (cString != NULL)
              {
                sscanf(cString,"%s %s %s %lg",winput1,winput2,
                       winput3, &ddata);
                for (kk = 0; kk < rsModelnInps; kk++)
                {
                  if (!strcmp(winput2,pPtr.strArray_[kk]))
                  {
                    rsIndices_[kk] = -1;
                    rsValues_[kk] = ddata;
                    break;
                  }
                }
                if (kk == rsModelnInps) break;
              }
            }
            if (ii != nInps)
            {
              printf("WARNING: Config info on fixed variables not used.\n");
              delete [] rsIndices_;
              delete [] rsValues_;
              rsIndices_ = NULL;
              rsValues_ = NULL;
            }
            else rsnInps_ = rsModelnInps;
          }
          for (kk = 0; kk < rsModelnInps; kk++)
          {
            if (rsIndices_[kk] == -1) 
              printf("Input %4d fixed at %e\n",kk+1,rsValues_[kk]);
          }
        }
      }
      delete psIO;
      rsPtrs_ = new FuncApprox*[nOutputs_];
      for (ii = 0; ii < nOutputs_; ii++)
      {
        printf("Creating response surface for output %d\n", ii+1);
        rsPtrs_[ii] = (FuncApprox *) genFAFromFile(fname,ii);
        if (rsPtrs_[ii] == NULL)
        {
          printf("FunctionInterface ERROR: no RS model given.\n");
          useRSModel_ = 0;
          exit(1);
        }
      }
    }
  }
  return 0;
}

// ************************************************************************
// set execution mode to be asynchronous
// ------------------------------------------------------------------------
int FunctionInterface::setAsynchronousMode()
{
  executionMode_ = 1;
  return 0;
}

// ************************************************************************
// set execution mode to be synchronous
// ------------------------------------------------------------------------
int FunctionInterface::setSynchronousMode()
{
  executionMode_ = 0;
  return 0;
}

// ************************************************************************
// set stochastic response surface mode (1:minus 3 sigma, 2: plus 3 sigma)
// ------------------------------------------------------------------------
int FunctionInterface::setStochasticRSMode(int indata)
{
  rsRanFlag_ = indata;
  return 0;
}

// ************************************************************************
// set time interval between asynchronous jobs
// ------------------------------------------------------------------------
int FunctionInterface::setLaunchInterval(int interval)
{
  launchInterval_ = interval;
  if (launchInterval_ < 0) launchInterval_ = 0;
  if (launchInterval_ > 100000) launchInterval_ = 0;
  return 0;
}

// ************************************************************************
// set print level
// ------------------------------------------------------------------------
int FunctionInterface::setOutputLevel(int level)
{
  printLevel_ = level;
  if (printLevel_ < 0)  printLevel_ = 0;
  if (printLevel_ > 10) printLevel_ = 0;
  return 0;
}

// ************************************************************************
// set which driver to use
// ------------------------------------------------------------------------
int FunctionInterface::setDriver(int which)
{
  int    ii, kk, nInps, status;
  char   inString[200], fname[200], winput1[1001], winput2[1001];
  char   winput3[1001], *cString, pString[1001];
  double ddata;
  FILE   *fp;
  PsuadeData *psIO;
  pData pPtr;

  if (rsPtrs_ != NULL)
  {
    for (ii = 0; ii < nOutputs_; ii++)
      if (rsPtrs_[ii] != NULL) delete rsPtrs_[ii];
    delete [] rsPtrs_;
    rsPtrs_ = NULL;
  }
  if (rsIndices_ != NULL)
  {
    delete [] rsIndices_;
    rsIndices_ = NULL;
  }
  if (rsValues_ != NULL)
  {
    delete [] rsValues_;
    rsValues_ = NULL;
  }

  if      (which >= 2) appOptFlag_ = 2;
  else if (which == 1) appOptFlag_ = 1;
  else                 appOptFlag_ = 0;
  useRSModel_ = 0;
  if      (appOptFlag_ == 0) strcpy(fname, appDriver_);
  else if (appOptFlag_ == 1) strcpy(fname, optDriver_);
  if      (appOptFlag_ == 2) strcpy(fname, auxOptDriver_);
  if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
       strcmp(optDriver_,"PSUADE_LOCAL"))
  {
    useRSModel_ = 2;
    return 0; 
  }
  if ( strcmp(appDriver_,"PSUADE_LOCAL") &&
      !strcmp(optDriver_,"PSUADE_LOCAL"))
  {
    useRSModel_ = 3;
    return 0; 
  }
  if (!strcmp(appDriver_,"PSUADE_LOCAL") &&
      !strcmp(optDriver_,"PSUADE_LOCAL"))
  {
    useRSModel_ = 4;
    return 0; 
  }

  fp = fopen(fname, "r");
  if (fp != NULL)
  {
    fscanf(fp, "%10c", inString);
    if (!strncmp(inString, "PSUADE_IO",9)) useRSModel_ = 1;
    fclose(fp);
    if (useRSModel_ == 1)
    {
      psIO = new PsuadeData();
      psIO->setOutputLevel(0);
      status = psIO->readPsuadeFile(fname);
      if (status != 0)
      {
        printf("ERROR : cannot read file %s in PSUADE format\n",fname);
        exit(1);
      }
      psIO->getParameter("input_ninputs", pPtr);
      int rsModelnInps = pPtr.intData_;
      if (rsModelnInps < nInputs_ || rsModelnInps > nInputs_)
      {
        printAsterisks(PL_INFO,0);
        printf("FunctionInterface: setting up RS driver from %s.\n",
               fname);
        printEquals(PL_INFO,0);
        printf("WARNING: nInputs in RS driver does not match nInputs");
        printf(" in original psuade file.\n");
        printf("   nInputs in original psuade file = %d\n",nInputs_);
        printf("   nInputs in RS driver data file  = %d\n",
               rsModelnInps);
        psIO->getParameter("ana_rsindexfile", pPtr);
        if (!strcmp(pPtr.strArray_[0], "NONE"))
        {
          printf("ERROR: nInputs mismatch and missing rs_index_file.\n");
          printf("   nInputs in original psuade file = %d\n",nInputs_);
          printf("   nInputs in RS driver data file  = %d\n",
                 rsModelnInps);
          printf("ADVICE: Put rs_index_file in %s or in the",fname);
          printf(" original psuade file.\n");
          exit(1);
        }
        printf("WARNING: rs_index_file found in the RS driver file.\n");
        fp = fopen(pPtr.strArray_[0], "r");
        if (fp == NULL)
        {
          printf("ERROR: missing rs_index_file %s in current folder.\n",
                 pPtr.strArray_[0]);
          exit(1);
        }
        else
        {
          fscanf(fp,"%d", &nInps);
          rsnInps_ = nInps;
          if (rsModelnInps != nInps)
          {
            printf("ERROR: invalid nInputs in rs index file.\n");
            printf("  nInputs read     = %d\n", nInps);
            printf("  nInputs expected = %d\n", rsModelnInps);
            printf("  Data format in rs index file should be: \n");
            printf("  line 1: nInputs in RS driver data file\n");
            printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
            printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
            printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
            printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
            printf("  ...\n");
            exit(1);
          }
          rsIndices_ = new int[nInps];
          rsValues_ = new double[nInps];
          for (ii = 0; ii < nInps; ii++) rsIndices_[ii] = 0;
          for (ii = 0; ii < nInps; ii++)
          {
            fscanf(fp, "%d", &kk);
            if (kk != ii+1)
            {
              printf("ERROR: first index in rs index file = %d.\n",
                     rsIndices_[ii]);
              printf("       Must be equal to %d.\n",ii+1);
              printf("  Data format in rs index file should be: \n");
              printf("  line 1: nInputs in RS driver data file\n");
              printf("  line 2: 1 <num> <num == 0 ==> fixed>\n");
              printf("  line 3: 2 <num> <0 if num != 0 (active)>\n");
              printf("  line 4: 3 <num> <num == 0 ==> fixed>\n");
              printf("  line 5: 4 <num> <0 if num != 0 (active)>\n");
              printf("  ...\n");
              exit(1);
            }
            fscanf(fp, "%d", &rsIndices_[ii]);
            if (rsIndices_[ii] < 0 || rsIndices_[ii] > nInputs_)
            {
              printf("INFO: input %3d = %d not valid\n",ii+1,
                     rsIndices_[ii]);
              printf("       Need to be between 1 and %d\n",nInputs_);
              exit(1);
            }
            rsIndices_[ii]--;
            fscanf(fp, "%lg", &rsValues_[ii]);
            if (rsIndices_[ii] == -1)
              printf("   RS Input %d inactive, fixed at %16.8e\n",
                     ii+1, rsValues_[ii]);
            else
            {
              printf("   RS Input %d   active, mapped to",ii+1);
              printf(" PSUADE input %d\n", rsIndices_[ii]+1);
            }
          }
          fclose(fp);
        }
      }
      psIO->getParameter("input_names", pPtr);
      if (psConfig_ != NULL && rsIndices_ == NULL && pPtr.strArray_ != NULL)
      {
        cString = psConfig_->getParameter("num_fixed");
        if (cString != NULL)
        {
          sscanf(cString, "%s %s %d", winput1, winput2, &nInps);
          if (nInps > 0)
          {
            rsIndices_ = new int[rsModelnInps];
            rsValues_ = new double[rsModelnInps];
            for (ii = 0; ii < rsModelnInps; ii++) rsIndices_[ii] = ii;
            for (ii = 0; ii < nInps; ii++)
            {
              sprintf(pString,"fixed-%d", ii+1);
              cString = psConfig_->getParameter(pString);
              if (cString != NULL)
              {
                sscanf(cString,"%s %s %s %lg",winput1,winput2,
                       winput3, &ddata);
                for (kk = 0; kk < rsModelnInps; kk++)
                {
                  if (!strcmp(winput2,pPtr.strArray_[kk]))
                  {
                    rsIndices_[kk] = -1;
                    rsValues_[kk] = ddata;
                    break;
                  }
                }
                if (kk == rsModelnInps) break;
              }
            }
            if (ii != nInps)
            {
              printf("WARNING: Config info on fixed variables not used.\n");
              delete [] rsIndices_;
              delete [] rsValues_;
              rsIndices_ = NULL;
              rsValues_ = NULL;
            }
            else rsnInps_ = nInps;
          }
          for (kk = 0; kk < rsModelnInps; kk++)
          {
            if (rsIndices_[kk] == -1) 
              printf("Input %4d fixed at %e\n",kk+1,rsValues_[kk]);
          }
        }
      }
      delete psIO;
      rsPtrs_ = new FuncApprox*[nOutputs_];
      for (ii = 0; ii < nOutputs_; ii++)
      {
        if (printLevel_ > 3)
          printf("Creating response surface for output %d.\n", ii+1);
        rsPtrs_[ii] = (FuncApprox *) genFAFromFile(fname,ii);
        if (rsPtrs_[ii] == NULL)
        {
          printf("FunctionInterface setDriver ERROR: no RS model.\n");
          useRSModel_ = 0;
          exit(1);
        }
      }
    }
  }
  return 0;
}

// ************************************************************************
// run application
// ------------------------------------------------------------------------
// flag = 0, this job has never been launched before 
// flag = 1, just create input file and return (for limited launch)
// flag = 2, check the status of this job
// ------------------------------------------------------------------------
int FunctionInterface::evaluate(int sampleID,int nInputs,double *inputs, 
                                int nOutputs, double *outputs, int flag)
{
  int    ii, outputCount, length, nfixed, status;
  double value, *myInputs, stdev;
  char   lineIn[500], command[500], winput1[500], winput2[500];
  char   outfile[500], infile[500], *cString, equal[100];
  FILE   *fp, *fIn, *fOut;

  if (nInputs_ != nInputs || nOutputs_ != nOutputs)
  {
    printf("FunctionInterface evaluate ERROR: nInputs/nOutputs mismatch.\n");
    printf("   nInputs  = %d versus %d (local)\n", nInputs, nInputs_);
    printf("   nOutputs = %d versus %d (local)\n", nOutputs, nOutputs_);
    exit(1);
  }

  if ((useRSModel_ == 0 || (useRSModel_ >= 2 && useRSModel_ <= 4)) && 
      rsPtrs_ == NULL)
  {
    if (appOptFlag_ == 0)
    {
      length = strlen(appInputTemplate_);
      for (ii = length-1; ii >= 0; ii--)
        if (appInputTemplate_[ii] == '/') break; 
      sprintf(infile, "%s.%d", &(appInputTemplate_[ii+1]), sampleID+1);
    }
    else sprintf(infile, "psuadeOpt.in.%d", sampleID+1);

    if (appOptFlag_ == 0)
    {
      length = strlen(appOutputTemplate_);
      for (ii = length-1; ii >= 0; ii--)
        if (appOutputTemplate_[ii] == '/') break; 
      sprintf(outfile, "%s.%d", &(appOutputTemplate_[ii+1]), sampleID+1);
    }
    else sprintf(outfile, "psuadeOpt.out.%d", sampleID+1);

    if ((fIn=fopen(outfile, "r")) == NULL) 
    {
      if (flag == 2) return 2;

      if (appOptFlag_ == 0 && (useRSModel_ == 2 || useRSModel_ == 4))
      {
        psLocalFunction(nInputs, inputs, nOutputs, outputs);
        return 0;
      }
      if (appOptFlag_ == 1 && (useRSModel_ == 3 || useRSModel_ == 4))
      {
        psLocalFunction(nInputs, inputs, nOutputs, outputs);
        return 0;
      }

      fOut = fopen(infile, "w");
      if (fOut == NULL)
      {
         printf("FunctionInterface ERROR: cannot open %s file\n",infile);
         exit(1);
      }
      fprintf(fOut, "%d\n", nInputs);
      for (ii = 0; ii < nInputs; ii++)
        fprintf(fOut, "%20.12e\n", inputs[ii]);
      if (psConfig_ != NULL)
      {
        cString = psConfig_->getParameter("num_fixed");
        if (cString != NULL)
        {
          sscanf(cString, "%s %s %d", winput1, winput2, &nfixed);
          fprintf(fOut,"num_fixed = %d\n", nfixed);
          ii = 0;
          while (ii < nfixed)
          {
            ii++;
            sprintf(winput1, "fixed-%d",ii);
            cString = psConfig_->getParameter(winput1);
            if (cString != NULL)
            {
              sscanf(cString, "%s %s %s %lg",winput1,winput2,equal,&value);
              fprintf(fOut,"fixed %d %s = %24.16e\n",ii,winput2,value);
            }
          }
        }
      }
      fclose(fOut);
      if (flag == 1) return 2;

      if (appOptFlag_ == 0)
      {
        if ((!strcmp(appDriver_, "NONE") || !strcmp(appDriver_, "true"))) 
        {
          printf("FunctionInterface ERROR: app driver not found.\n");
          exit(1);
        }
      }
      if (appOptFlag_ == 1)
      {
        if ((!strcmp(optDriver_, "NONE") || !strcmp(optDriver_, "true")))
        {
          printf("FunctionInterface ERROR: opt driver not found.\n");
          exit(1);
        }
      }
      if (appOptFlag_ == 2)
      {
        if ((!strcmp(auxOptDriver_, "NONE") || 
             !strcmp(auxOptDriver_, "true")))
        {
          printf("FunctionInterface ERROR: aux opt driver not found.\n");
          exit(1);
        }
      }

      if (appOptFlag_ == 0)
      {
        if (executionMode_ == 1)
          sprintf(command, "\"%s\" %s %s %d %d %d&",appDriver_,infile,
                  outfile, sampleID, flag, printLevel_);
        else
          sprintf(command, "\"%s\" %s %s %d %d %d", appDriver_,infile,
                   outfile, sampleID, flag, printLevel_);
        if (strstr((const char*) appDriver_, "rm ") != NULL ||
            strstr((const char*) appDriver_, "mv ") != NULL ||
            strstr((const char*) appDriver_, " -f ") != NULL ||
            strstr((const char*) appDriver_, "/bin/") != NULL) 
        {
          printf("FunctionInterface::evaluate ERROR: \n");
          printf("\t\t for security reason do not use rm in driver.\n");
          exit(1);
        }
      }
      else if (appOptFlag_ == 1)
      {
        if (executionMode_ == 1)
          sprintf(command, "%s %s %s %d %d %d&",optDriver_,infile,outfile,
                  sampleID+1, flag, printLevel_);
        else
          sprintf(command, "%s %s %s %d %d %d",optDriver_,infile,outfile,
                  sampleID+1, flag, printLevel_);
        if (strstr((const char*) optDriver_, "rm ") != NULL ||
            strstr((const char*) optDriver_, "mv ") != NULL ||
            strstr((const char*) optDriver_, " -f ") != NULL ||
            strstr((const char*) optDriver_, "/bin/") != NULL) 
        {
          printf("FunctionInterface::evaluate ERROR: \n");
          printf("\t\t for security reason do not use rm in driver.\n");
          exit(1);
        }
      }
      else if (appOptFlag_ == 2)
      {
        sprintf(command, "%s %s %s %d %d %d", auxOptDriver_, infile, 
                outfile, sampleID+1, flag, printLevel_);
        if (strstr((const char*) auxOptDriver_, "rm ") != NULL ||
            strstr((const char*) auxOptDriver_, "mv ") != NULL ||
            strstr((const char*) auxOptDriver_, " -f ") != NULL ||
            strstr((const char*) auxOptDriver_, "/bin/") != NULL) 
        {
          printf("FunctionInterface::evaluate ERROR: \n");
          printf("\t\t for security reason do not use rm in driver.\n");
          exit(1);
        }
      }
      status = system(command);   
      if (status != 0)
      {
        printf("FunctionInterface evaluate ERROR: system call returns %d.\n",
               status);
        printf("  INFO: system call return status should be 0.\n");
        printf("  INFO: check your simulation driver for correctiness.\n");
        exit(1);
      }
      if (executionMode_ == 1 && launchInterval_ > 0)
      {
#ifdef WINDOWS
        Sleep(1000 * launchInterval_);
#else
        sleep(launchInterval_);
#endif
      }
    }
    else if (flag == 0)
    {
      printf("WARNING: Output file %s exists before it is run.\n",outfile);
      printf("WARNING: PSUADE will use this output file.\n");
      printf("WARNING: If this is a mistake, stop PSUADE and clean up.\n");
      fclose(fIn);
    }
    else fclose(fIn);

    if (executionMode_ == 0)
    {
      //if (launchInterval_ == 0) launchInterval_ = 1;

      while ((fIn=fopen(outfile, "r")) == NULL)
      {
        if (printLevel_ > 2)
        {
          printf("Waiting for Job %d to complete.\n", sampleID+1);
          printf("If you run the simulator yourself, use the inputs\n");
          if (appOptFlag_ == 0)
          {
            printf("from psuadeApps_ct.in.%d for your simulation and\n", 
                   sampleID+1);
            printf("write the outputs to psuadeApps_ct.out.%d\n",
                   sampleID+1);
          }
          else if (appOptFlag_ == 1)
          {
            printf("from psuadeOpt.in.%d for your simulation and\n", 
                   sampleID+1);
            printf("write the outputs to psuadeOpt.out.%d\n",sampleID+1);
          }
        }
#ifdef WINDOWS
        Sleep(1000 * launchInterval_);
#else
        sleep(launchInterval_);
#endif
      }
      fclose(fIn);
    }
    if (executionMode_ == 1)
    {
      if ((fIn=fopen(outfile, "r")) == NULL) 
      {
        return 2;
      }
      fclose(fIn);
    }
    length = 0;
    for (ii = 0; ii < 100000; ii++) length = (length + ii) %  32768;

    fIn = fopen(outfile, "r");
    if (fIn == NULL) return 2;
    outputCount = 0;
    for (ii = 0; ii < nOutputs_; ii++)
    {
      fgets(lineIn, 100, fIn);
      sscanf(lineIn, "%lg", &(outputs[outputCount]));
      outputCount++;
      if (feof(fIn) == 1) break; 
    }
    fclose(fIn);

    if (outputCount != nOutputs) 
    {
      printf("\t\t output file %s found but nOutputs mismatch (%d,%d).\n", 
             outfile, outputCount, nOutputs);
      printf("\t\t (check your output format).\n");
      return 2;
    }

    if (strcmp(infile, "*"))
    {
      unlink(infile);
    }
    if (strcmp(outfile, "*"))
    {
      unlink(outfile);
    }
  }
  else if (useRSModel_ == 1 && rsPtrs_ != NULL)
  {
    if (rsIndices_ != NULL)
    {
      myInputs = new double[rsnInps_];
      for (ii = 0; ii < rsnInps_; ii++)
      {
        if (rsIndices_[ii] < 0) myInputs[ii] = rsValues_[ii];
        else                    myInputs[ii] = inputs[rsIndices_[ii]];
      }
      if (rsRanFlag_ == 0)
      {
        for (ii = 0; ii < nOutputs; ii++)
          outputs[ii] = rsPtrs_[ii]->evaluatePoint(myInputs);
      }
      else
      {
        for (ii = 0; ii < nOutputs; ii++)
        {
          outputs[ii] = rsPtrs_[ii]->evaluatePointFuzzy(myInputs,stdev);
          value = 3.0 * stdev;
          if (rsRanFlag_ == 1) outputs[ii] -= value;
          else                 outputs[ii] += value;
        } 
      }
      delete [] myInputs;
    }
    else
    {
      if (rsRanFlag_ == 0)
      {
        for (ii = 0; ii < nOutputs; ii++)
          outputs[ii] = rsPtrs_[ii]->evaluatePoint(inputs);
      }
      else
      {
        for (ii = 0; ii < nOutputs; ii++)
        {
          outputs[ii] = rsPtrs_[ii]->evaluatePointFuzzy(inputs,stdev);
          value = 3.0 * stdev;
          if (rsRanFlag_ == 1) outputs[ii] -= value;
          else                 outputs[ii] += value;
        }
      }
    }
  }
  else
  {
    printf("FunctionInterface ERROR: evaluate error.\n");
    if (useRSModel_ == 1)
      printf("       Did you forget to declare driver?\n");
    exit(1);
  } 
  return 0;
}

// ************************************************************************
// run ensemble simulation
// ------------------------------------------------------------------------
int FunctionInterface::ensembleEvaluate(int nSamp,int nInputs,double *inputs, 
                               int nOutputs, double *outputs, int ID)
{
  int    ii, ss, outputCount, nfixed, status;
  double value;
  char   outfile[500], infile[500], lineIn[5000], command[500];
  char   winput1[500], winput2[500], winput3[500], *cString;
  FILE   *fp, *fIn, *fOut;

  if (nInputs_ != nInputs || nOutputs_ != nOutputs)
  {
    printf("FunctionInterface ERROR: nInputs/nOutputs mismatch.\n");
    printf("   nInputs  = %d versus %d (local)\n", nInputs, nInputs_);
    printf("   nOutputs = %d versus %d (local)\n", nOutputs, nOutputs_);
    printf("NOTE: ensembleEvaluate expects driver is an actual simulator");
    printf(" and not RS data file.\n");
    printf("      Therefore, it does not take a rs_index_file.\n");
    exit(1);
  }

  sprintf(infile, "psuadeEval.in.%d", ID);
  sprintf(outfile, "psuadeEval.out.%d", ID);

  if ((fIn=fopen(outfile, "r")) == NULL) 
  {
    if (appOptFlag_ == 0)
    {
      if (!strcmp(ensembleDriver_,"PSUADE_LOCAL"))
      {
        psEnsembleLocalFunction(nSamp,nInputs,inputs,nOutputs,outputs);
        return 0;
      } 
    }
    else
    {
      if (!strcmp(ensembleOptDriver_,"PSUADE_LOCAL"))
      {
        psEnsembleLocalFunction(nSamp,nInputs,inputs,nOutputs,outputs);
        return 0;
      } 
    }

    fOut = fopen(infile, "w");
    if (fOut == NULL)
    {
      printf("FunctionInterface ERROR: cannot open %s file\n",infile);
      exit(1);
    }
    fprintf(fOut, "%d %d\n", nSamp, nInputs);
    for (ss = 0; ss < nSamp; ss++)
    {
      for (ii = 0; ii < nInputs; ii++)
        fprintf(fOut, "%20.12e ", inputs[ss*nInputs+ii]);
      fprintf(fOut, "\n");
    }
    if (psConfig_ != NULL)
    {
      cString = psConfig_->getParameter("num_fixed");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %d", winput1, winput2, &nfixed);
        fprintf(fOut,"num_fixed = %d\n", nfixed);
        ss = 0;
        for (ii = 0; ii < nfixed; ii++)
        {
          sprintf(winput1, "fixed-%d",ii+1);
          cString = psConfig_->getParameter(winput1);
          if (cString != NULL)
          {
            sscanf(cString, "%s %s %s %lg", winput1,winput2,winput3,&value);
            fprintf(fOut,"fixed %d %s = %24.16e\n", ii, winput2, value);
            ss++;
          }
          else printf("FuncInterface ERROR: %s not found.\n",winput1);
        }
        if (ss != nfixed)
        {
          printf("FuncInterface ERROR: fixed variables not found in config.\n");
          psConfig_->print();
          exit(1);
        }
      }
    }
    fclose(fOut);

    if (appOptFlag_ == 0)
    {
      if ((!strcmp(ensembleDriver_, "NONE") || 
           !strcmp(ensembleDriver_, "true")))
      {
        printf("FunctionInterface ERROR: ensemble driver not set.\n");
        exit(1);
      }
      fp = fopen(ensembleDriver_, "r");
      if (fp == NULL)
      {
        printf("FunctionInterface ERROR: ensemble driver %s not found.\n",
               ensembleDriver_);
        exit(1);
      }
      else fclose(fp);
      sprintf(command, "%s %s %s %d %d", ensembleDriver_,infile,outfile,ID, 
              printLevel_);
    }
    else
    {
      if ((!strcmp(ensembleOptDriver_, "NONE") || 
           !strcmp(ensembleOptDriver_, "true")))
      {
        printf("FunctionInterface ERROR: ensemble opt driver not set.\n");
        exit(1);
      }
      fp = fopen(ensembleOptDriver_, "r");
      if (fp == NULL)
      {
        printf("FunctionInterface ERROR: ensemble opt driver %s not found\n",
               ensembleOptDriver_);
        exit(1);
      }
      else fclose(fp);
      sprintf(command, "%s %s %s %d %d", ensembleOptDriver_,infile,outfile,
              ID, printLevel_);
    }

    status = system(command);   
    if (status != 0)
    {
      printf("FunctionInterface INFO: system call returns status = %d.\n",
             status);
      printf("  INFO: system call return status should be 0.\n");
      printf("  INFO: check your simulation driver for correctiness.\n");
      exit(1);
    }
  }
  else 
  {
    printf("WARNING: Output file %s exists before it is run.\n",outfile);
    printf("WARNING: PSUADE will use this output file.\n");
    printf("WARNING: If this is a mistake, stop PSUADE and clean up.\n");
    fclose(fIn);
  }

  while ((fIn=fopen(outfile, "r")) == NULL)
  {
    if (printLevel_ > 2)
    {
      printf("FunctionInterface: waiting for Job %d to complete.\n",ID);
#ifdef WINDOWS
      Sleep(1000 * launchInterval_);
#else
      sleep(launchInterval_);
#endif
    }
  }

  fIn = fopen(outfile, "r");
  if (fIn == NULL)
  {
    printf("FunctionInterface ERROR: output file %s not found.\n",outfile);
    exit(1);
  }
  outputCount = 0;
  for (ss = 0; ss < nSamp; ss++)
  {
    fgets(lineIn, 4000, fIn);
    for (ii = 0; ii < nOutputs_; ii++)
      sscanf(lineIn, "%lg", &(outputs[ss*nOutputs_+ii]));
    outputCount++;
    if (feof(fIn) == 1) break; 
  }
  fclose(fIn);

  if (outputCount != nSamp) 
  {
    printf("FunctionInterface ERROR: output file %s found but with\n",
           outfile);
    printf("                         insufficient data.\n");
    printf("Advice: Check the output format of your aux opt driver.\n");
    printf("        It should have %d lines each with %d output data.\n",
           nSamp, nOutputs_);
    exit(1);
  }

  if (strcmp(infile, "*")) unlink(infile);
  if (strcmp(outfile, "*")) unlink(outfile);
  return 0;
}

// ************************************************************************
// get number of input variables
// ------------------------------------------------------------------------
int FunctionInterface::getNumInputs()
{
  return nInputs_;
}

// ************************************************************************
// get number of output variables
// ------------------------------------------------------------------------
int FunctionInterface::getNumOutputs()
{
  return nOutputs_;
}

// ************************************************************************
// get current driver code
// ------------------------------------------------------------------------
int FunctionInterface::getDriver()
{
  return appOptFlag_;
}

// ************************************************************************
// get names of input variables
// ------------------------------------------------------------------------
char **FunctionInterface::getInputNames()
{
  return inputNames_;
}

// ************************************************************************
// get names of output variables
// ------------------------------------------------------------------------
char **FunctionInterface::getOutputNames()
{
  return outputNames_;
}

// ************************************************************************
// get application driver
// ------------------------------------------------------------------------
char *FunctionInterface::getApplicationDriver()
{
  return appDriver_;
}

// ************************************************************************
// get optimization driver
// ------------------------------------------------------------------------
char *FunctionInterface::getOptimizationDriver()
{
  return optDriver_;
}

// ************************************************************************
// get auxiliary optimization driver
// ------------------------------------------------------------------------
char *FunctionInterface::getAuxOptimizationDriver()
{
  return auxOptDriver_;
}

// ************************************************************************
// Local functions (use setLocalFunction to select)
// (some of these local functions are used for, for example, optimal
//  experimental design - coupled with some PSUADE's global optimizers). 
// ------------------------------------------------------------------------
int FunctionInterface::psLocalFunction(int nInputsIn, double *inputsIn,
                                       int nOutputs, double *outputs)
{
  int    ss, ii, jj, kk, ll, ind, lcnt, status, iOne=1, count;
  double dmean, ddata, dtmp, aggrVar, dvar, GMetric, IMetric, DMetric;
  double AMetric, EMetric;
  char   fname[2000], pString[2000], lineIn[2000];
  pData  pdata;
  psMatrix matCov;

  //**/ =============================================================
  //**/ make a local copy of inputs
  //**/ This is needed because user can call this function with
  //**/ nInputsIn=0 and inputsIn=NULL
  //**/ -------------------------------------------------------------
  int    nInputs = nInputsIn;
  double *inputs = NULL;
  psVector vecInps;
  if (nInputsIn > 0) 
  { 
    vecInps.setLength(nInputsIn);
    if (inputsIn != NULL)
      for (ii = 0; ii < nInputs; ii++) vecInps[ii] = inputsIn[ii];
    inputs = vecInps.getDVector();
  }
 
  //**/ =============================================================
  //**/ data structures to be used for optimal design
  //**/ =============================================================
  static int nHist, nInps=3, nOuts=1, maxHist=10000;
  static int ProblemInitialized=0;
  static psMatrix  matHistory;
  static psIVector vecUInputs, vecIT;
  static psMatrix matCandidates;
  static psMatrix matPriorSample;
  static psMatrix matEvalSet;
  static ProbMatrix **CandPostSamples=NULL;
  static FuncApprox **rsPtrs=NULL;
  static int nPreSelected=0;
  static double optimalVal;
  static McmcData mobj14;
  static psVector VecYPrior;
  //**/ -------------------------------------------------------------
  //**/ the following is used to select OUU test problem 
  //**/ -------------------------------------------------------------
  static int toyOption = -1;
 
  //**/ =============================================================
  //**/ option 999: clean up (after optimization is completed) 
  //**/ -------------------------------------------------------------
  //**/ The call sequence to use these optimiztions is:
  //**/ - instantiate FunctionInterface ==> funcIO
  //**/ - call funcIO->setLocalFunction(option) to select method
  //**/ - feed this funcIO to some optimization method, which
  //**/   calls this repeatedly for design evaluations
  //**/ - when completed, call funcIIO->setLocalFunction(999) 
  //**/ =============================================================
  if (whichLocalFunction_ == 999)
  {
    printf("FuncIO Local INFO: cleaning up.\n");
    if (CandPostSamples != NULL)
    {
      for (ii = 0; ii < matCandidates.nrows(); ii++)
        if (CandPostSamples[ii] != NULL)
          delete CandPostSamples[ii];
      delete [] CandPostSamples;
    } 
    CandPostSamples = NULL;
    if (rsPtrs != NULL)
    {
      for (ii = 0; ii < nOuts; ii++) 
        if (rsPtrs[ii] != NULL) delete rsPtrs[ii];
      delete [] rsPtrs;
    }
    rsPtrs = NULL;
    ProblemInitialized = 0;
    optimalVal = PSUADE_UNDEFINED;
    matPriorSample.cleanUp();
    matEvalSet.cleanUp();
    VecYPrior.clean();
  }

  //**/ =============================================================
  //**/ option: ODOE slow G/I/A/D/E-optimal with SCE optimization
  //**/ (this option, as opposed to LocalFunction 10-13, uses
  //**/ MCMC at every iteration and thus may be slower)
  //**/ -------------------------------------------------------------
  //**/ Bayes  20: G; 21: I; 22: A; 23: D; 24: E 
  //**/ =============================================================
  else if (whichLocalFunction_ >= 20 && whichLocalFunction_ <= 24)
  {
    //**/ -----------------------------------------------------------
    //**/ error checking (the objective function must have size 1
    //**/ even though the simulator has multiple outputs
    //**/ -----------------------------------------------------------
    if (nOutputs != 1 && nInputsIn > 0)
    {
      printf("FuncIO ERROR: nOutputs has to be = 1\n");
      exit(1);
    }
    psVector vecXT, vecYT;
    int nUInps;

    //**/ -----------------------------------------------------------
    //**/ initialization 
    //**/ -----------------------------------------------------------
    if (ProblemInitialized == 0)
    {
      printf("ODOE_(X)OPTIMAL: Initialization begins ..\n");

      //**/ --- read training sample ==> psuadeIO
      sprintf(pString,
          "Enter name of the training sample (for creating RS): ");
      getString(pString, fname);
      fname[strlen(fname)-1] = '\0';
      PsuadeData *psuadeIO = new PsuadeData();
      status = psuadeIO->readPsuadeFile(fname);
      if (status != 0)
      {
        printf("FuncIO ERROR when reading training sample\n");
        printf("       Maybe this file is in wrong format?\n");
        exit(1);
      }
      psuadeIO->getParameter("input_ninputs", pdata);
      nInps = pdata.intData_;
      psuadeIO->getParameter("output_noutputs", pdata);
      nOuts = pdata.intData_;

      //**/ --- user selects uncertain inputs ==> vecUInputs
      vecUInputs.setLength(nInps);
      sprintf(pString,
         "Enter uncertain input number (1 - %d, 0 to end) : ",nInps);
      ii = 0;
      while (1)
      {
        kk = getInt(0, nInps, pString);
        if (kk == 0 || kk > nInps) break;
        vecUInputs[ii] = kk - 1;
        ii++;
      }
      vecUInputs.subvector(0, ii-1);
      nUInps = vecUInputs.length();

      //**/ --- set uncertain parameter input to 1 in vecIT
      //**/ --- e.g. vecIT[kk]=1 if input kk is uncertain
      vecIT.setLength(nInps);
      kk = 1;
      for (ii = 0; ii < vecUInputs.length(); ii++)
      {
        vecIT[vecUInputs[ii]] = kk;
        kk++; 
      }

      //**/ --- read prior sample ==> matPriorSample
      printf("Uncertain parameters need a prior sample for inference.\n");
      sprintf(pString, "Enter the file name of your prior sample : ");
      getString(pString, fname);
      fname[strlen(fname)-1] = '\0';
      matPriorSample.setFormat(PS_MAT2D); // This is needed for later sort
      status = readIReadDataFile(fname, matPriorSample);
      if (status != 0)
      {
        printf("FuncIO ERROR when reading prior sample\n");
        printf("       Maybe this file is in wrong format?\n");
        exit(1);
      }
      if (matPriorSample.ncols() != vecUInputs.length())
      {
        printf("FuncIO ERROR: prior nInputs is not correct.\n");
        printf("       Should be equal to %d.\n",vecUInputs.length());
        exit(1);
      }

      //**/ --- read candidate set ==> matCandidates
      //**/ --- If no inputs is given, selected set is requested
      if (nInputs == 0)
        sprintf(pString,"Enter the file name of your selected set : ");
      else
        sprintf(pString,"Enter the file name of your candidate set : ");
      getString(pString, fname);
      fname[strlen(fname)-1] = '\0';
      status = readIReadDataFile(fname, matCandidates);
      if (status != 0)
      {
        printf("FuncIO ERROR when reading candidate set\n");
        printf("       Maybe this file is in wrong format?\n");
        exit(1);
      }
      int nCandidates = matCandidates.nrows();
      printf("Size of Candidate set = %d\n", nCandidates);
      if (nInputs != 0 && matCandidates.ncols() != 2*nOuts+nInps-nUInps)
      {
        printf("FuncIO ERROR: candidate file must have %d columns.\n",
               2*nOuts+nInps-nUInps);
        printf("Suggestion: use odoeu_rseval to append your ");
        printf("candidate set with output\n");
        printf("            means and weights.\n");
        exit(1);
      }
      if (nInputs == 0 && matCandidates.ncols() < nInps-nUInps)
      {
        printf("FuncIO ERROR: selected set file must have %d columns.\n",
               2*nOuts+nInps-nUInps);
        exit(1);
      }
#if 0
      if (matCandidates.ncols() == 2*nOuts+nInps-nUInps)
      {
        for (ii = 0; ii < nOuts; ii++)
        {
          for (jj = 0; jj < matCandidates.nrows(); jj++)
          {
            ddata = matCandidates.getEntry(jj,nInps-nUInps+2*ii);
            if (ddata <= 0)
            {
              printf("funcIO WARNING: candidate %d output %d <= 0\n",
                     jj+1, ii+1);
            }
            ddata = matCandidates.getEntry(jj,nInps-nUInps+2*ii+1);
            if (ddata < 0.1 || ddata > 1)
            {
              printf("funcIO WARNING: candidate %d weight not in [0.1,1]\n",
                     jj+1);
              printf("INFO: reset output %d weight to 1.\n",ii+1);
              matCandidates.setEntry(jj,nInps-nUInps+2*ii+1,1.0);
            }
          }
        }
      }
#endif

      //**/ --- read evaluation set (for G and I) ==> matEvalSet
      if (whichLocalFunction_ == 20 || whichLocalFunction_ == 21)
      {
        printf("An evaluation sample is needed to ");
        printf("compute the optimality\n");
        printf("metrics. This can be the same as the ");
        printf("candidate set.\n");
        sprintf(pString,
                "Enter the file name of your evaluation set : ");
        getString(pString, fname);
        fname[strlen(fname)-1] = '\0';
        status = readIReadDataFile(fname, matEvalSet);
        if (status != 0)
        {
          printf("FuncIO ERROR when reading evaluation set\n");
          printf("       Maybe this file is in wrong format?\n");
          exit(1);
        }
        if (matEvalSet.ncols() != 2*nOuts+nInps-nUInps &&
            matEvalSet.ncols() != nInps-nUInps)
        {
          printf("FuncIO ERROR: evaluation file must have %d or %d columns.\n",
                 2*nOuts+nInps-nUInps,nInps-nUInps);
          exit(1);
        }
      }

      //**/ --- construct response surface (psuadeIO ==> rsPtr)
      pData pInps, pOuts, pLBs, pUBs;
      psuadeIO->getParameter("method_nsamples", pdata);
      int nSamp = pdata.intData_;
      psuadeIO->getParameter("input_lbounds", pLBs);
      psuadeIO->getParameter("input_ubounds", pUBs);
      psuadeIO->getParameter("input_sample", pInps);
      psuadeIO->getParameter("output_sample", pOuts);
  
      int faFlag = 1, rsMethod=0;
      vecYT.setLength(nSamp);
      rsPtrs = new FuncApprox*[nOuts];
      int interactiveSave = psInteractive_;
      interactiveSave = 0;
      for (ii = 0; ii < nOuts; ii++)
      {
        if (ii == 0) 
        {
          rsPtrs[ii] = genFAInteractive(psuadeIO, faFlag);
          rsMethod = rsPtrs[ii]->getID();
        }
        else 
        {
          //**/ all outputs use the same RS method
          psuadeIO->updateAnalysisSection(-1,-1,rsMethod,-1,-1,-1);
          faFlag = 0;
          rsPtrs[ii] = genFAInteractive(psuadeIO, faFlag);
        }
        rsPtrs[ii]->setBounds(pLBs.dbleArray_,pUBs.dbleArray_);
        rsPtrs[ii]->setOutputLevel(0);
        for (jj = 0; jj < nSamp; jj++)
          vecYT[jj] = pOuts.dbleArray_[jj*nOuts+ii];
        status = rsPtrs[ii]->initialize(pInps.dbleArray_,
                                        vecYT.getDVector());
      }
      psInteractive_ = interactiveSave;

      //**/ --- the rest of the initialization
      delete psuadeIO;
      ProblemInitialized = 1;
      matHistory.setDim(maxHist,nInputs+1);
      nHist = 0;
      optimalVal = PSUADE_UNDEFINED;

      //**/ --- initialize MCMC object for subsequent processing
      mobj14.printLevel_ = 0;
      mobj14.nSamples_ = nSamp;
      mobj14.nInputs_ = nInps;
      mobj14.nOutputs_ = nOuts;
      mobj14.VecLowerB_.load(nInps, pLBs.dbleArray_);
      mobj14.VecUpperB_.load(nInps, pUBs.dbleArray_);
      mobj14.VecSamInputs_.load(nInps*nSamp, pInps.dbleArray_);
      mobj14.VecSamOutputs_.load(nOuts*nSamp, pOuts.dbleArray_);
      mobj14.VecCUInputs_ = vecUInputs;
      mobj14.faType_ = rsMethod;
      mobj14.MatPriorSample_ = matPriorSample;
      printf("ODOE_(X)OPTIMAL Initialization complete.\n");
    }

    //**/ -----------------------------------------------------------
    //**/ if inputs = NULL, create an inputs of all candidates
    //**/ (assume that user wants to evaluate with all candidates
    //**/ -----------------------------------------------------------
    if (nInputsIn == 0)
    {
      nInputs = matCandidates.nrows();
      vecInps.setLength(nInputs);
      inputs = vecInps.getDVector();
      for (ii = 0; ii < nInputs; ii++) inputs[ii] = ii + 1;
    }

    //**/ -----------------------------------------------------------
    //**/ check to see if there are duplicates (duplication selection
    //**/ is not allowed). If so, just return a large value
    //**/ -----------------------------------------------------------
    count = 0;
    for (ii = 0; ii < nInputs; ii++)
    {
      for (jj = ii+1; jj < nInputs; jj++)
      { 
        ind = inputs[ii];
        kk  = inputs[jj];
        if (ind == kk) count++;
      }
    }
    if (count > 0)
    {
      if (nInputsIn == 0)
      {
        for (jj = 0; jj < nOutputs; jj++) outputs[jj] = 0;
        printf(" ==> duplicate selection ");
      }
      if (optimalVal < 0) outputs[0] = 0.1 * optimalVal;
      else                outputs[0] = 10.0 * optimalVal;
      printf(" ==> duplicate selection ");
      for (ii = 0; ii < nInputs; ii++) printf("%d ", (int) inputs[ii]);
      printf(" ==> skip\n");
      return 0;
    }

    //**/ -----------------------------------------------------------
    //**/ error checking (whether inputs are valid) and display
    //**/ -----------------------------------------------------------
    if (nInputsIn > 0 && printLevel_ >= 0) 
    {
      if      (whichLocalFunction_ == 20) printf("GOPTIMAL inputs: ");
      else if (whichLocalFunction_ == 21) printf("IOPTIMAL inputs: ");
      else if (whichLocalFunction_ == 22) printf("DOPTIMAL inputs: ");
      else if (whichLocalFunction_ == 23) printf("AOPTIMAL inputs: ");
      else if (whichLocalFunction_ == 24) printf("EOPTIMAL inputs: ");
    }
    for (ii = 0; ii < nInputs; ii++)
    {
      ind = (int) inputs[ii];
      if (ind < 1 || ind > matCandidates.nrows())
      {
        printf("ODOE Function Evaluator ERROR: wrong input values.\n");
        printf("     Check candidate set size consistency.\n");
        printf("The erroneous inputs are:\n");
        for (jj = 0; jj < nInputs; jj++)
          printf("Candidate %3d for evaluation = %d\n", jj+1, 
                 (int) inputs[jj]);
        printf("But they should all be in the range of [1,%d]\n",
               matCandidates.nrows());
        exit(1);
      }
      if (nInputsIn > 0 && printLevel_ >= 0) printf("%5d ", ind);
    }

    //**/ -----------------------------------------------------------
    //**/ check history to see whether this has been evaluated before
    //**/ if so, just fetch it from history record
    //**/ -----------------------------------------------------------
    for (jj = 0; jj < nHist; jj++)
    {
      count = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
        for (kk = 0; kk < nInputs; kk++)
          if (matHistory.getEntry(jj,kk) == inputs[ii])
            break;
        if (kk != nInputs) count++;
      }
      if (count == nInputs)
      {
        outputs[0] = matHistory.getEntry(jj,nInputs);
        if (printLevel_ >= 0) 
          printf(" ===> output = %e (revisit)\n", outputs[0]);
        return 0; 
      }
    }

    //**/ -----------------------------------------------------------
    //**/ For the other Bayesian metrics (G, I, D, A, E)
    //**/ instantiate MCMC object 
    //**/ -----------------------------------------------------------
    int ind2, nDesInps = nInps - vecUInputs.length();
    psMatrix matExpInps, matExpMeans, matExpStds;
    matExpInps.setDim(nInputs, nInps-vecUInputs.length());
    matExpMeans.setDim(nInputs, nOuts);
    matExpStds.setDim(nInputs, nOuts);
    MCMCAnalyzer *mcmcAnalyzer = new MCMCAnalyzer();
    vecXT.setLength(nInps);

    //**/ -----------------------------------------------------------
    //**/ fill in the experimental data matrix for MCMC
    //**/ (matExpInps, matExpMeans, matExpStds)
    //**/ -----------------------------------------------------------
    for (ii = 0; ii < nInputs; ii++)
    {
      ind2 = (int) inputs[ii] - 1;
      lcnt = 0;
      for (jj = 0; jj < nInps; jj++)
      {
        //**/ if parameter is uncertain, compute vecXT = mean
        if (vecIT[jj] >= 1)
        {
          ind = vecIT[jj] - 1;
          dmean = 0;
          for (kk = 0; kk < matPriorSample.nrows(); kk++)
            dmean += matPriorSample.getEntry(kk,ind);
          vecXT[jj] = dmean / (double) matPriorSample.nrows();
        }
        //**/ else get the vecXT from candidate matrix
        else
        {
          vecXT[jj] = matCandidates.getEntry(ind2,lcnt);
          matExpInps.setEntry(ii, lcnt, vecXT[jj]);
          lcnt++;
        }
      }
      for (jj = 0; jj < nOuts; jj++)
      {
        //**/ get experimental mean/std from user 
        dmean = matCandidates.getEntry(ind2,nDesInps+jj*2);
        matExpMeans.setEntry(ii,jj,dmean);
        ddata = matCandidates.getEntry(ind2,nDesInps+jj*2+1);
        matExpStds.setEntry(ii,jj,ddata);
      }
    }

    //**/ -----------------------------------------------------------
    //**/ stuff the MCMC data object and call MCMC
    //**/ -----------------------------------------------------------
    mobj14.MatExpInputs_ = matExpInps;
    mobj14.MatExpMeans_ = matExpMeans;
    mobj14.MatExpStds_ = matExpStds;
    mcmcAnalyzer->analyzeDirect(mobj14);

    //**/ -----------------------------------------------------------
    //**/ extract posterior sample and evaluate
    //**/ -----------------------------------------------------------
    double **postSample = mobj14.MatPostSample_.getMatrix2D();
    GMetric = IMetric = AMetric = EMetric = 0;
    DMetric = 1;

    //**/ -----------------------------------------------------------
    //**/ only G- and I-optimal metrics need the following
    //**/ -----------------------------------------------------------
    if (whichLocalFunction_ == 20 || whichLocalFunction_ == 21)
    {
      //**/ #########################################################
      //**/ this portion calls rsPtr, and is done at each iteration.
      //**/ so if RS evaluation is expensive, this portion may be
      //**/ slow (2020)
      //**/ #########################################################
      vecXT.setLength(nInps);
      vecYT.setLength(mobj14.MatPostSample_.nrows());
      for (ii = 0; ii < matEvalSet.nrows(); ii++)
      {
        lcnt = 0;
        for (jj = 0; jj < nInps; jj++)
        {
          if (vecIT[jj] == 0)
          {
            vecXT[jj] = matEvalSet.getEntry(ii,lcnt);
            lcnt++;
          }
        }
        aggrVar = 0;
        for (jj = 0; jj < nOuts; jj++)
        {
          for (kk = 0; kk < mobj14.MatPostSample_.nrows(); kk++)
          {
            for (ll = 0; ll < vecUInputs.length(); ll++)
            {
              ind = vecUInputs[ll];
              vecXT[ind] = postSample[kk][ll];
            }
            ddata = rsPtrs[jj]->evaluatePoint(vecXT.getDVector());
            vecYT[kk] = ddata;
          }
          dmean = dvar = 0;
          for (kk = 0; kk < mobj14.MatPostSample_.nrows(); kk++)
            dmean += vecYT[kk];
          dmean /= (double) mobj14.MatPostSample_.nrows();
          for (kk = 0; kk < mobj14.MatPostSample_.nrows(); kk++)
            dvar += pow(vecYT[kk]-dmean,2);
          dvar = dvar / (double) mobj14.MatPostSample_.nrows();
          aggrVar += dvar;
        }
        aggrVar /= nOuts;
        if (aggrVar > GMetric) GMetric = aggrVar;
        IMetric += aggrVar;
      }
      IMetric /= (double) matEvalSet.nrows();
    }
    //**/ -----------------------------------------------------------
    //**/ D-, A-, and E-optimal metrics need the following
    //**/ It is done even G- or I-metrics are requested
    //**/ -----------------------------------------------------------
    psMatrix matCov, matEig;
    psVector vecEigs;
    matCov.setDim(vecUInputs.length(),vecUInputs.length());
    double dmean2, dcov;
    AMetric = 0.0;
    for (jj = 0; jj < vecUInputs.length(); jj++)
    {
      dmean = 0.0;
      for (kk = 0; kk < mobj14.MatPostSample_.nrows(); kk++)
        dmean += postSample[kk][jj];
      dmean /= (double) mobj14.MatPostSample_.nrows();
      for (kk = 0; kk < mobj14.MatPostSample_.nrows(); kk++)
        dvar += pow(postSample[kk][jj]-dmean,2.0);
      dvar /= (double) mobj14.MatPostSample_.nrows();
      matCov.setEntry(jj,jj,dvar);
      AMetric += dvar;
      for (ll = jj+1; ll < vecUInputs.length(); ll++)
      {
        dmean2 = 0;
        for (kk = 0; kk < mobj14.MatPostSample_.nrows(); kk++)
          dmean2 += postSample[kk][ll];
        dmean2 /= (double) mobj14.MatPostSample_.nrows();
        dcov = 0;
        for (kk = 0; kk < mobj14.MatPostSample_.nrows(); kk++)
          dcov += (postSample[kk][jj]-dmean)*
                  (postSample[kk][ll]-dmean2);
        dcov /= (double) mobj14.MatPostSample_.nrows();
        matCov.setEntry(jj,ll,dcov);
        matCov.setEntry(ll,jj,dcov);
      }
    }
    matCov.eigenSolve(matEig, vecEigs, 1);
    DMetric = 1.0;
    EMetric = 0;
    for (jj = 0; jj < vecUInputs.length(); jj++)
    {
      DMetric *= vecEigs[jj];
      if (vecEigs[jj] > EMetric) EMetric = vecEigs[jj];
    }

    //**/ -----------------------------------------------------------
    //**/ return all metric information
    //**/ -----------------------------------------------------------
    if (nInputsIn == 0)
    {
      if (nOutputs >= 1) outputs[0] = GMetric;
      if (nOutputs >= 2) outputs[1] = IMetric;
      if (nOutputs >= 3) outputs[2] = DMetric;
      if (nOutputs >= 4) outputs[3] = AMetric;
      if (nOutputs >= 5) outputs[4] = EMetric;
      return 0;
    }

    //**/ -----------------------------------------------------------
    //**/ return the relevant metric information
    //**/ -----------------------------------------------------------
    if      (whichLocalFunction_ == 20) outputs[0] = GMetric;
    else if (whichLocalFunction_ == 21) outputs[0] = IMetric;
    else if (whichLocalFunction_ == 22) outputs[0] = DMetric;
    else if (whichLocalFunction_ == 23) outputs[0] = AMetric;
    else if (whichLocalFunction_ == 24) outputs[0] = EMetric;
    if (outputs[0] < optimalVal) optimalVal = outputs[0];
    if (nInputsIn > 0 && printLevel_ >= 0) 
    {
      printf(" ===> output = %e (best so far = %e)\n",outputs[0],
             optimalVal);
    }

    //**/ -----------------------------------------------------------
    //**/ update history for faster evaluation for revisits
    //**/ -----------------------------------------------------------
    if (nHist >= maxHist)
    {
      for (jj = 0; jj < maxHist/2; jj++)
        for (ii = 0; ii < nInputs+1; ii++)
          matHistory.setEntry(jj,ii,
                    matHistory.getEntry(jj+maxHist/2,ii));
      nHist = maxHist/2;
    }
    for (ii = 0; ii < nInputs; ii++)
      matHistory.setEntry(nHist,ii,inputs[ii]);
    matHistory.setEntry(nHist,nInputs,outputs[0]);
    nHist++;
    return 0;
  }

  //**/ =============================================================
  //**/ option: G/I/A/D/E-optimal using the Fisher method and with 
  //**/ SCE optimization
  //**/ -------------------------------------------------------------
  //**/ Fisher 25: G; 26: I; 27: A; 28: D; 29: E 
  //**/ =============================================================
  else if (whichLocalFunction_ >= 25 && whichLocalFunction_ <= 29)
  {
    //**/ -----------------------------------------------------------
    //**/ error checking (the objective function must have size 1)
    //**/ But if nInputsIn == 0, it is asking for multiple metric
    //**/ evaluation (PSUADE interpreter is calling this from
    //**/ odoeu_feval) so it is acceptable to have nOutputs > 1 
    //**/ (in this case, nOutputs are metrics and not simulation
    //**/ outputs)
    //**/ -----------------------------------------------------------
    if (nInputsIn > 0 && nOutputs != 1)
    {
      printf("FuncIO ERROR: nOutputs has to be = 1\n");
      exit(1);
    }
    psVector vecXT, vecYT;
    int nUInps, nSamFisher;

    //**/ -----------------------------------------------------------
    //**/ initialization 
    //**/ -----------------------------------------------------------
    if (ProblemInitialized == 0)
    {
      if (printLevel_ > 0)
        printf("ODOE_(X)OPTIMAL (Fisher): Initialization begins ..\n");

      //**/ --- read training sample ==> psuadeIO
      printf("The key to this class of methods is to create ");
      printf("the Fisher information\n");
      printf("matrices and then derives metrics from them. In ");
      printf("order to create the\n");
      printf("Fisher information matrix, the following information ");
      printf("are needed:\n");
      printf("- A training sample to compute the Fisher matrix\n");
      printf("  (and identify design and uncertain parameters)\n");
      printf("- A prior distribution for the uncertain parameters\n");
      printf("- A candidate set (or a design set) of experimental designs\n");

      //**/ get the training sample
      sprintf(pString,
          "Enter name of the training sample (in PSUADE format): ");
      getString(pString, fname);
      fname[strlen(fname)-1] = '\0';
      PsuadeData *psuadeIO = new PsuadeData();
      status = psuadeIO->readPsuadeFile(fname);
      if (status != 0)
      {
        printf("FuncIO: ERROR encountered when reading the ");
        printf("training sample\n");
        printf("Possible reasons: \n");
        printf("- this file does not exist\n");
        printf("- this file is in wrong format (not PSUADE format)\n");
        exit(1);
      }
      psuadeIO->getParameter("input_ninputs", pdata);
      nInps = pdata.intData_;
      psuadeIO->getParameter("output_noutputs", pdata);
      nOuts = pdata.intData_;
      if (nOuts > 1)
      {
        printf("FuncIO ERROR: this method only works with nOutputs=1\n");
        printf("Suggestion: delete all but one output and re-run.\n");
        exit(1);
      }

      //**/ --- user selects uncertain inputs ==> vecUInputs
      printf("Out of the %d inputs, some should be design parameters ",
             nInputs_);
      printf("and some are\n");
      printf("uncertain parameters. In the following, please specify ");
      printf("which inputs\n");
      printf("are uncertain (and the rest are design parameters).\n");
      vecUInputs.setLength(nInps);
      sprintf(pString,
         "Enter uncertain input number (1 - %d, 0 to terminate) : ",nInps);
      ii = 0;
      while (1)
      {
        kk = getInt(0, nInps, pString);
        if (kk == 0 || kk > nInps) break;
        vecUInputs[ii] = kk - 1;
        ii++;
      }
      vecUInputs.subvector(0, ii-1);
      nUInps = vecUInputs.length();

      //**/ in order for the Fisher matrix to be well-conditioned, the
      //**/ following condition has to be met.
      if (nInputsIn != 0 && nInputsIn < nUInps)
      {
        printf("ERROR: design set must be larger than ");
        printf("uncertain parameter size.\n");
        printf("- design set has %d members\n",nInputsIn);
        printf("- number of uncertain parameters = %d\n",nUInps);
        exit(1);
      }

      //**/ --- set uncertain parameter input to 1 in vecIT
      //**/ --- e.g. vecIT[kk]=1 if input kk is uncertain
      vecIT.setLength(nInps);
      kk = 1;
      for (ii = 0; ii < vecUInputs.length(); ii++)
      {
        vecIT[vecUInputs[ii]] = kk;
        kk++; 
      }

      //**/ --- read prior sample ==> matPriorSample
      printf("Uncertain parameters need a prior sample. The sample ");
      printf("file should\n");
      printf("have the following format:\n");
      printf("Line 1: <nSamples> <nInputs>\n");
      printf("Line 2: 1 <sample 1 values>\n");
      printf("Line 3: 2 <sample 2 values>\n");
      printf("Line 4: ...\n");
      sprintf(pString, "Enter the file name of your prior sample : ");
      getString(pString, fname);
      fname[strlen(fname)-1] = '\0';
      matPriorSample.setFormat(PS_MAT2D); // This may not be needed
      status = readIReadDataFile(fname, matPriorSample);
      if (status != 0)
      {
        printf("FuncIO: ERROR encountered when reading the prior sample\n");
        printf("Possible reasons: \n");
        printf("- this file does not exist\n");
        printf("- this file is in wrong format\n");
        exit(1);
      }
      if (matPriorSample.ncols() != vecUInputs.length())
      {
        printf("FuncIO ERROR: prior sample nInputs=%d is not correct.\n",
               matPriorSample.ncols());
        printf("       Should be equal to %d.\n",vecUInputs.length());
        exit(1);
      }

      //**/ option to use smaller prior sample (this will not be 
      //**/ asked if this is called by odoeu_eval)
      if (nInputsIn > 0 && (matPriorSample.nrows() > 10))
      {
        printf("Fisher-based methods are computationally expensive, so ");
        printf("you may want\n");
        printf("to reduce the cost by collapsing the prior sample into ");
        printf("fewer sample\n");
        printf("points (prior sample size = %d).\n",matPriorSample.nrows());
        sprintf(pString,
           "Collapse prior sample into smaller sample ? (y or n) \n"); 
        getString(pString, lineIn);
        if (lineIn[0] == 'y')
        {
          printf("The size of the prior sample is %d.\n",
                 matPriorSample.nrows());
          sprintf(pString,
            "Use (1) the sample mean only or (2) a random sub-sample ? ");
          kk = getInt(1, 2, pString);
          if (kk == 1)
          {
            int nr = matPriorSample.nrows(), nc = matPriorSample.ncols();
            psMatrix tmpMat = matPriorSample;
            matPriorSample.setFormat(PS_MAT2D);
            matPriorSample.setDim(iOne, nc);
            for (ii = 0; ii < nc; ii++)
            {
              ddata = 0;
              for (jj = 0; jj < nr; jj++)
                ddata += tmpMat.getEntry(jj, ii);
              ddata /= (double) nr;
              matPriorSample.setEntry(0, ii, ddata); 
            }
          }
          else if (kk == 2)
          {
            sprintf(pString,"Enter sub-sample size (%d - %d) : ",
                    nInputs, matPriorSample.nrows()-1);
            nSamFisher = getInt(nInputs,matPriorSample.nrows()-1,pString);
            
            int nr = matPriorSample.nrows(), nc = matPriorSample.ncols();
            psMatrix tmpMat = matPriorSample;
            matPriorSample.setFormat(PS_MAT2D);
            matPriorSample.setDim(nSamFisher, nc);
            int nCurrent = 0, nTrials=0, checkFlag=1;
            while (nCurrent < nSamFisher)
            {
              nTrials++;
              if (nTrials > 20*nr)
              {
                printf("INFO: Unable to draw a sub-sample with ");
                printf("unique sample points.\n");
                printf("      Stop searching for a unique sub-sample.\n");
                checkFlag = 0;
              }
              kk = PSUADE_rand() % nr;
              //**/ check that there are no duplicates
              jj = 0;
              if (checkFlag == 1)
              {
                for (ii = 0; ii < nCurrent-1; ii++)
                {  
                  for (jj = 0; jj < nc; jj++)
                  {
                    ddata = matPriorSample.getEntry(ii,jj);
                    dtmp  = tmpMat.getEntry(kk,jj);
                    if (dtmp != ddata) break;
                  }
                  //**/ give up if the same sample point
                  if (jj == nc) break;
                }
              }
              if (jj != nc)
              {
                for (jj = 0; jj < nc; jj++)
                {
                  ddata = tmpMat.getEntry(kk,jj);
                  matPriorSample.setEntry(nCurrent,jj,ddata);
                }
                nCurrent++;
              }
            }
          }
          else nSamFisher = matPriorSample.nrows(); 
        }
      }

      //**/ --- read candidate set ==> matCandidates
      //**/ --- If no inputs is given (inputsIn = NULL), it means 
      //**/ --- the calling function is requesting using selected
      //**/ --- designs for analysis instead of parameters passed 
      //**/ --- to this function
      if (nInputs == 0 && nOutputs != 0)
      {
        printf("Next please provide provide a selected design set.\n");
        printf("The file should be in the following format:\n");
        printf("Line 1: <nSamples> <nInputs>\n");
        printf("Line 2: <selected design point 1 values>\n");
        printf("Line 3: <selected design point 2 values>\n");
        printf("Line 4: ...\n");
        sprintf(pString,
                "Enter the file name of the selected set of designs : ");
      }
      else
      {
        printf("Next please provide a candidate design list.\n");
        printf("The file should be in the following format:\n");
        printf("Line 1: <nSamples> <nInputs>\n");
        printf("Line 2: <candidate design 1 values>\n");
        printf("Line 3: <candidate design 2 values>\n");
        printf("Line 4: ...\n");
        sprintf(pString,
                "Enter the file name of your candidate design set : ");
      }
      getString(pString, fname);
      fname[strlen(fname)-1] = '\0';
      status = readIReadDataFile(fname, matCandidates);
      if (status != 0)
      {
        printf("FuncIO: ERROR encountered when reading candidate set\n");
        printf("Possible reasons: \n");
        printf("- this file does not exist\n");
        printf("- this file is in wrong format\n");
        exit(1);
      }
      int nCandidates = matCandidates.nrows();
      if (nInputsIn == 0 && nOutputs != 0)
      {
        printf("Size of the selected design set = %d\n", nCandidates);
        if (nCandidates < nUInps)
        {
          printf("ERROR: selected design set should have > %d members.\n",
                 nUInps);
          exit(1);
        }
      }
      else printf("Size of the candidate set = %d\n", nCandidates);

      //**/ --- read evaluation set (for G and I) ==> matEvalSet
      if (whichLocalFunction_ == 25 || whichLocalFunction_ == 26)
      {
        printf("An evaluation sample is needed to ");
        printf("compute the optimality metrics.\n");
        printf("This can be the same as the candidate set");
        printf(" (but not recommended).\n");
        printf("The file should be in the following format: \n");
        printf("Line 1: <numPoints> <nInputs>\n");
        printf("Line 2: 1 <sample 1 values> \n");
        printf("Line 3: 2 <sample 2 values> \n");
        printf("....\n");
        sprintf(pString,
                "Enter the file name of your evaluation set : ");
        getString(pString, fname);
        fname[strlen(fname)-1] = '\0';
        status = readIReadDataFile(fname, matEvalSet);
        if (status != 0)
        {
          printf("FuncIO ERROR when reading evaluation set\n");
          printf("Possible reasons: \n");
          printf("- this file does not exist\n");
          printf("- this file is in wrong format\n");
          exit(1);
        }
        if (matEvalSet.ncols() != 2*nOuts+nInps-nUInps &&
            matEvalSet.ncols() != nInps-nUInps)
        {
          printf("FuncIO ERROR: evaluation data must have %d or %d columns.\n",
                 2*nOuts+nInps-nUInps,nInps-nUInps);
          exit(1);
        }
      }

      //**/ --- construct response surface (psuadeIO ==> rsPtr)
      //**/ --- RS is needed to compute Fisher matrix
      pData pInps, pOuts, pLBs, pUBs;
      psuadeIO->getParameter("method_nsamples", pdata);
      int nSamp = pdata.intData_;
      psuadeIO->getParameter("input_lbounds", pLBs);
      psuadeIO->getParameter("input_ubounds", pUBs);
      psuadeIO->getParameter("input_sample", pInps);
      psuadeIO->getParameter("output_sample", pOuts);
  
      int faFlag = 1, rsMethod=0;
      vecYT.setLength(nSamp);
      rsPtrs = new FuncApprox*[1];
      rsPtrs[0] = genFAInteractive(psuadeIO, faFlag);
      rsPtrs[0]->setBounds(pLBs.dbleArray_,pUBs.dbleArray_);
      rsPtrs[0]->setOutputLevel(0);
      for (jj = 0; jj < nSamp; jj++)
        vecYT[jj] = pOuts.dbleArray_[jj];
      status = rsPtrs[0]->initialize(pInps.dbleArray_,
                                     vecYT.getDVector());

      //**/ --- the rest of the initialization
      delete psuadeIO;
      ProblemInitialized = 1;
      matHistory.setDim(maxHist,nInputs+1);
      nHist = 0;
      optimalVal = PSUADE_UNDEFINED;
      if (printLevel_ > 0)
        printf("ODOE_(X)OPTIMAL (Fisher) Initialization complete.\n");
    }

    //**/ -----------------------------------------------------------
    //**/ if nInputsIn = 0 and inputs = NULL, create an inputs of all 
    //**/ candidates because in this case, it is assumed that callers
    //**/ want to use the entire candidate set instead of 'inputs' 
    //**/ (e.g. called from odoeu_feval)
    //**/ Otherwise, copy inputsIn to vecInps
    //**/ -----------------------------------------------------------
    if (nInputsIn == 0 && inputsIn == NULL)
    {
      nInputs = matCandidates.nrows();
      vecInps.setLength(nInputs);
      inputs = vecInps.getDVector();
      for (ii = 0; ii < nInputs; ii++) inputs[ii] = ii + 1;
    }
    else
    {
      vecInps.setLength(nInputsIn);
      inputs = vecInps.getDVector();
      for (ii = 0; ii < nInputs; ii++) inputs[ii] = inputsIn[ii];
    }

    //**/ -----------------------------------------------------------
    //**/ check to see if there are duplicates (duplication selection
    //**/ is not allowed). If so, just return a large value
    //**/ -----------------------------------------------------------
    count = 0;
    for (ii = 0; ii < nInputs; ii++)
    {
      for (jj = ii+1; jj < nInputs; jj++)
      { 
        ind = (int) vecInps[ii];
        kk  = (int) vecInps[jj];
        if (ind == kk) count++;
      }
    }
    //**/ if the inputs is from the candidate set and it has 
    //**/ duplicates, compress it
    if (count > 0 && nInputsIn == 0)
    {
      if (nInputsIn == 0)
      {
        sortDbleList(nInputs, vecInps.getDVector());
        count = 1;
        for (ii = 1; ii < nInputs; ii++)
        {
          if (vecInps[count-1] != vecInps[ii]) 
          {
            vecInps[count] = vecInps[ii];
            count++;
          }
        }
        vecInps.subvector(0, nInputs-1);
      }
    }
    //**/ if inputs are from optimizer, return large values
    else if (count > 0) 
    {
      if (optimalVal < 0) outputs[0] = 0.1 * optimalVal;
      else                outputs[0] = 10.0 * optimalVal;
      printf(" ==> duplicate selection (not a user concern) : ");
      for (ii = 0; ii < nInputs; ii++) printf("%d ",(int) vecInps[ii]);
      printf(" ==> skip\n");
      return 0;
    }

    //**/ -----------------------------------------------------------
    //**/ error checking (whether inputs are valid) and display
    //**/ -----------------------------------------------------------
    if (nInputsIn > 0 && printLevel_ >= 0) 
    {
      if (whichLocalFunction_ == 25)
      {
        if (nInputsIn > 0)
          printf("GOPTIMAL (Fisher) inputs: ");
        else
          printf("<All>OPTIMAL (Fisher) inputs: ");
      }
      else if (whichLocalFunction_ == 26)
        printf("IOPTIMAL (Fisher) inputs: ");
      else if (whichLocalFunction_ == 27)
        printf("DOPTIMAL (Fisher) inputs: ");
      else if (whichLocalFunction_ == 28)
        printf("AOPTIMAL (Fisher) inputs: ");
      else if (whichLocalFunction_ == 29)
        printf("EOPTIMAL (Fisher) inputs: ");
    }
    for (ii = 0; ii < nInputs; ii++)
    {
      ind = (int) vecInps[ii];
      if (ind < 1 || ind > matCandidates.nrows())
      {
        printf("FuncIO ERROR: Wrong input values.\n");
        printf("              Check candidate set size consistency.\n");
        printf("The erroneous inputs are:\n");
        for (jj = 0; jj < nInputs; jj++)
          printf("Candidate %3d for evaluation = %d\n", jj+1, 
                 (int) vecInps[jj]);
        printf("But they should all be in the range of [1,%d]\n",
               matCandidates.nrows());
        exit(1);
      }
      if (nInputsIn > 0 && printLevel_ >= 0) printf("%5d ", ind);
    }

    //**/ -----------------------------------------------------------
    //**/ search history to see if this has been evaluated before
    //**/ -----------------------------------------------------------
    for (jj = 0; jj < nHist; jj++)
    {
      count = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
        for (kk = 0; kk < nInputs; kk++)
          if (matHistory.getEntry(jj,kk) == inputs[ii])
            break;
        if (kk != nInputs) count++;
      }
      if (count == nInputs)
      {
        outputs[0] = matHistory.getEntry(jj,nInputs);
        if (printLevel_ >= 0) 
          printf(" ===> output = %e (revisit)\n", outputs[0]);
        return 0; 
      }
    }

    //**/ -----------------------------------------------------------
    //**/ if nInputs = nOutputs = 0, it means odoeu_feval is calling
    //**/ this function to evaluate the GIDAE metrics for all the
    //**/ candidate points in the candidate set.
    //**/ vecInps (double) has the candidate indices (1-based)
    //**/ -----------------------------------------------------------
    vecXT.setLength(nInps); /* nInps = nInputs in training sample */
    vecYT.setLength(nOuts); /* nOuts = 1 */
    nUInps = vecUInputs.length(); /* number of uncertain inputs */

    //**/ -----------------------------------------------------------
    //**/ this 'else' segment is a different implementation than the
    //**/ 'if' segment, in the sense that this segment computes many
    //**/ M's = [dY(c_1)/d theta ... dY(c_n)/d theta] where M is the 
    //**/ fisher matrix, c_i is candidate i, and Y is simulation 
    //**/ output of interest. For each M, optimality metrics are
    //**/ computed and average over all prior sample points
    //**/ -----------------------------------------------------------
    psVector vecEig, vecCT, vecMCT;
    psMatrix matGrad, matGradT, matCovInv, matEig;
    int priorNR = matPriorSample.nrows(), cc, cand;
    int nCandidates = vecInps.length();
    double DMetric2=0,AMetric2=0,EMetric2=0;
    double GMetric2=0, IMetric2=0;
    vecMCT.setLength(nUInps);

    //**/ the process is repeated for each prior sample
    matGrad.setDim(nUInps, nCandidates);
    for (ss = 0; ss < priorNR; ss++)
    {
      //**/ --- first stuff the prior sample into vecXT
      for (jj = 0; jj < nInps; jj++)
      {
        //**/ if parameter is uncertain, use prior sample
        if (vecIT[jj] >= 1)
        {
          ind = vecIT[jj] - 1;
          vecXT[jj] = matPriorSample.getEntry(ss,ind);
        }
      }
      //**/ initialize matGrad

      //**/ fille matGrad with gradient wrt uncertain
      //**/ parameters for each candidate point
      for (cc = 0; cc < nCandidates; cc++)
      {
        cand = (int) (vecInps[cc]) - 1; 
        //**/ stuff candidate coordinate in vecXT
        //**/ (from candidate matrix)
        lcnt = 0;
        for (jj = 0; jj < nInps; jj++)
        {
          if (vecIT[jj] < 1) /* vecIT[jj] = 0 for design input jj */
          {
            vecXT[jj] = matCandidates.getEntry(cand,lcnt);
            lcnt++;
          } 
        }
        //**/ --- evaluate function and compute derivatives too
        rsPtrs[0]->evaluatePoint(iOne,vecXT.getDVector(),&ddata);
        vecYT[0] = ddata;

        //**/ compute partial y(theta) /partial theta_j
        for (jj = 0; jj < nInps; jj++)
        {
          //**/ if parameter is uncertain, perturb and evaluate
          if (vecIT[jj] >= 1)
          {
            ind = vecIT[jj] - 1;
            dtmp = vecXT[jj];
            vecXT[jj] *= (1.0 + 1e-6);
            rsPtrs[0]->evaluatePoint(iOne,vecXT.getDVector(),&ddata);
            //**/ finite difference delta Y wrt theta_jj for 
            //**/ candidate cc
            ddata = (ddata - vecYT[0]) / (vecXT[jj] - dtmp);
            vecXT[jj] = dtmp;
            ddata += matGrad.getEntry(ind, cc);
            matGrad.setEntry(ind, cc, ddata); 
          }
        }
      } /* cc for nCandidates */
    }
    //**/ Should be included
    for (cc = 0; cc < nCandidates; cc++)
    {
      for (jj = 0; jj < nUInps; jj++)
      {
        ddata = matGrad.getEntry(jj, cc);
        ddata /= (double) priorNR;
        matGrad.setEntry(jj, cc, ddata);
      }
    }

    //**/ --- now compute information matrix (Grad * Grad^T) 
    //**/ --- covariance matrix ~ inv(Grad * Grad^T)
    matGradT = matGrad;
    matGradT.transpose();
    matGrad.matmult(matGradT, matCovInv);
    ddata = matCovInv.computeDeterminant(); 
    if (PABS(ddata) < 1e-16) 
    {
      printf("INFO: Singular Fisher matrix ==> skip.\n");
      if (printLevel_ > 1)
      {
        printf("Current Fisher matrix:\n");
        matCovInv.print();
      }
      GMetric2 = PSUADE_UNDEFINED;
      IMetric2 = PSUADE_UNDEFINED;
      AMetric2 = PSUADE_UNDEFINED;
      DMetric2 = PSUADE_UNDEFINED;
      EMetric2 = PSUADE_UNDEFINED;
    }
    else
    {
      DMetric2 += (1.0 / matCovInv.computeDeterminant()); 

      //**/ --- now get back the covariance matrix by inversion
      //**/ --- and compute A metric
      matCovInv.computeInverse(matCov);
      for (kk = 0; kk < nUInps; kk++)
        AMetric2 += matCov.getEntry(kk, kk);

      //**/ --- compute EMetric (minimize the maximum
      //**/ --- eigenvalue
      matCov.eigenSolve(matEig, vecEig, 1);
      ddata = vecEig[0];
      for (kk = 1; kk < nUInps; kk++)
        if (vecEig[kk] > ddata) ddata = vecEig[kk];
      EMetric2 += ddata;

      //**/ --- compute G and I metrics, if needed
      if (matEvalSet.nrows() > 0)
      {
        //**/ --- should be added
        //**/ --- first stuff the mean prior sample into vecXT
        vecXT.setLength(nInps);
        for (jj = 0; jj < nInps; jj++)
        {
          for (ss = 0; ss < priorNR; ss++)
          {
            if (vecIT[jj] >= 1)
            {
              ind = vecIT[jj] - 1;
              vecXT[jj] += matPriorSample.getEntry(ss,ind);
            }
          }
          vecXT[jj] /= (double) priorNR;
        }
        GMetric2 = 0;
        int nGood=0;
        for (cc = 0; cc < matEvalSet.nrows(); cc++)
        {
          //**/ stuff evaluation coordinate in vecXT
          //**/ uncertain input values already stuffed
          lcnt = 0;
          for (jj = 0; jj < nInps; jj++)
          { 
            if (vecIT[jj] < 1)
            {
              vecXT[jj] = matEvalSet.getEntry(cc,lcnt);
              lcnt++;
            }
          }
          //**/ --- evaluate function and compute derivatives too
          rsPtrs[0]->evaluatePoint(iOne,vecXT.getDVector(),&ddata);
          vecYT[0] = ddata;
  
          //**/ compute partial y(theta) /partial theta_j for
          //**/ evaluation point cc
          vecCT.setLength(nUInps);
          for (jj = 0; jj < nInps; jj++)
          {
            //**/ if parameter is uncertain, perturb and evaluate
            if (vecIT[jj] >= 1)
            {
              ind = vecIT[jj] - 1;
              dtmp = vecXT[jj];
              vecXT[jj] *= (1.0 + 1e-6);
              rsPtrs[0]->evaluatePoint(iOne,vecXT.getDVector(),&ddata);
              ddata = (ddata - vecYT[0]) / (vecXT[jj] - dtmp);
              vecXT[jj] = dtmp;
              vecCT[ind] += ddata;
            }
          }
          //**/compute inner product c^t Cov c
          matCov.matvec(vecCT,vecMCT,0);
          ddata = 0;
          for (jj = 0; jj < nUInps; jj++) 
            ddata += (vecCT[jj] * vecMCT[jj]);
          if (ddata > 0)
          {
            IMetric2 += ddata;
            if (ddata > GMetric2) GMetric2 = ddata;
            nGood++;
          }
        }
        if (nGood == 0)
        {
          if (printLevel_ > 0)
            printf("\nINFO: evaluation set not informative ==> skip.\n");
          IMetric2 = PSUADE_UNDEFINED;
          GMetric2 = PSUADE_UNDEFINED;
        }   
        else
        { 
          IMetric2 /= (double) nGood;
        }
      }
    }

    //**/ if it is called by odoeu_feval, display information
    if (nInputsIn == 0)
    {
      printAsterisks(PL_INFO,0);
      printf("  *** Fisher metric summary: \n");
      printDashes(PL_INFO,0);
      printf("  G-metric = %e\n", GMetric2);
      printf("  I-metric = %e\n", IMetric2);
      printf("  D-metric = %e\n", DMetric2);
      printf("  A-metric = %e\n", AMetric2);
      printf("  E-metric = %e\n", EMetric2);
      printAsterisks(PL_INFO,0);
      if (nOutputs >= 1) outputs[0] = GMetric2;
      if (nOutputs >= 2) outputs[1] = IMetric2;
      if (nOutputs >= 3) outputs[2] = DMetric2;
      if (nOutputs >= 4) outputs[3] = AMetric2;
      if (nOutputs >= 5) outputs[4] = EMetric2;
      return 0;
    }

    //**/ --- return the proper metric
    if (whichLocalFunction_ == 25) outputs[0] = GMetric2;
    if (whichLocalFunction_ == 26) outputs[0] = IMetric2;
    if (whichLocalFunction_ == 27) outputs[0] = DMetric2;
    if (whichLocalFunction_ == 28) outputs[0] = AMetric2;
    if (whichLocalFunction_ == 29) outputs[0] = EMetric2;
    if (outputs[0] < optimalVal) optimalVal = outputs[0];
    if (printLevel_ >= 0) 
    {
      printf(" ===> output = %e (best so far = %e)\n",outputs[0],
             optimalVal);
    }

    //**/ --- update iteration history
    if (nHist >= maxHist)
    {
      for (jj = 0; jj < maxHist/2; jj++)
        for (ii = 0; ii < nInputs+1; ii++)
          matHistory.setEntry(jj,ii,
                    matHistory.getEntry(jj+maxHist/2,ii));
      nHist = maxHist/2;
    }
    for (ii = 0; ii < nInputs; ii++)
      matHistory.setEntry(nHist,ii,inputs[ii]);
    matHistory.setEntry(nHist,nInputs,outputs[0]);
    nHist++;
    return 0;
  }
  return 0;
}

// ************************************************************************
// Local function
// ------------------------------------------------------------------------
int FunctionInterface::psEnsembleLocalFunction(int nSamples, int nInputs, 
                              double *inputs, int nOutputs, double *outputs)
{
  int    ss, ii;
  double X1, X2, X3, X4, D, W1, W2, W3, T1, T2, Y, D1, D2, D3, D4, W4;
  double alpha=5, beta=1, delta=-10, gamma3, gamma4;
  static int initializeProblem=2;

  if (printLevel_ > 0)
     printf("FunctionInterface local function: Dowling's toy problem.\n");
  if (initializeProblem == -1)
  {
    printf("Available test problem :\n");
    printf(" 1. Dowling's toy problem - function in single-stage OUU.\n");
    printf(" 2. Dowling's toy problem - optimizer in two-stage OUU.\n");
    printf("Which problem to use ? (1 or 2) ");
    scanf("%d", &initializeProblem);
    if (initializeProblem != 2) initializeProblem = 1;
    if (initializeProblem == 1)
         printf("Dowling's toy function for single-stage OUU selected.\n");
    else printf("Dowling's toy optimizer for two-stage OUU selected.\n");
  }
  if (initializeProblem == 1)
  {
    for (ss = 0; ss < nSamples; ss++)
    {
      D  = inputs[ss*nInputs];
      X1 = inputs[ss*nInputs+1];
      X2 = inputs[ss*nInputs+2];
      W1 = inputs[ss*nInputs+3];
      W2 = inputs[ss*nInputs+4];
      W3 = inputs[ss*nInputs+5];
      Y = pow(X1 - D + W1, 2.0) + (1 + W1 * W1) * pow(X2 - D + W2, 2.0) +
          (W1 + W3) * X1 + (W2 + W3) * X2 + pow(W2*W2+(2+W3*W3)*D, 2.0);
      outputs[ss] = Y;
    }
  }
  else
  {
    if (nInputs != 12)
    {
      printf("psEnsembleLocalFunction ERROR: nInputs do not match.\n");
      exit(1);
    }
    for (ss = 0; ss < nSamples; ss++)
    {
      D1  = inputs[ss*nInputs];
      D2 = inputs[ss*nInputs+1];
      D3 = inputs[ss*nInputs+2];
      D4 = inputs[ss*nInputs+3];
      X1 = inputs[ss*nInputs+4];
      X2 = inputs[ss*nInputs+5];
      X3 = inputs[ss*nInputs+6];
      X4 = inputs[ss*nInputs+7];
      W1 = inputs[ss*nInputs+8];
      W2 = inputs[ss*nInputs+9];
      W3 = inputs[ss*nInputs+10];
      W4 = inputs[ss*nInputs+11];
      X1 = - W1;
      gamma3 = 1 + W3 * W3;
      gamma4 = 1 + W4 * W4;
      X2 = - (delta + beta * D2 * gamma4 + beta * gamma4 * W2)/(beta*gamma4);
      X3 = - (delta * W3 + D3 * gamma3 * gamma4 + gamma3 * gamma4 * W3) / 
           (gamma3 * gamma4);
      T1 = delta * gamma3 + beta*delta*W3*W3 + beta * gamma3 * gamma4 * W2;
      T2 = beta * gamma3 * gamma4 * W3 * W3 + beta * D2 * gamma3 * gamma4;
      T1 = T1 + T2;
      T2 = beta * D4 * gamma3 * gamma4;
      T1 = T1 - T2;
      T2 = beta * delta * gamma3 * gamma4 + beta*D3 * gamma3 * gamma4 * W3;
      X4 = (T1 + T2) / (beta * gamma3 * gamma4 * gamma4);

      Y = (X1 + W1) * (X1 + W1);
      Y = Y + beta * pow(X2 + W2 + D2, 2.0);
      Y = Y + (1 + W3 * W3) * pow(X3 + D3 + W3, 2.0);
      T1 = pow(D4 + X2 + W3 * X3 + X4 * (1 + W4 * W4), 2.0);
      Y = Y + 1.0/(1+W4*W4) * T1;
      Y = Y - 2.0 * delta * X4;
      Y = Y + alpha * pow(D1+W1, 2.0);
      Y = Y + (10 - alpha) * D1 * D1;
      Y = Y + (10 - beta) * D2 * D2;
      Y = Y + 3 * D3 * D3;
      Y = Y + D4 * D4 * sqrt(1+W3*W3+W4*W4);
      outputs[ss] = Y;
    }
  }
  return 0;
}

// ************************************************************************
// create an FunctionInterface instantiation
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterface(PsuadeData *psuadeIO)
{
  int   nInputs, nOutputs;
  char  **inputNames, **outputNames;
  pData pPtr, pINames, pONames, pAppFiles;
  FunctionInterface *funcIO;

  funcIO = new FunctionInterface();
  psuadeIO->getParameter("input_ninputs", pPtr);
  nInputs = pPtr.intData_;
  psuadeIO->getParameter("input_names", pINames);
  inputNames = pINames.strArray_;
  psuadeIO->getParameter("output_noutputs", pPtr);
  nOutputs = pPtr.intData_;
  psuadeIO->getParameter("output_names", pONames);
  outputNames = pONames.strArray_;
  psuadeIO->getParameter("app_files", pAppFiles);
  funcIO->loadInputData(nInputs, inputNames);
  funcIO->loadOutputData(nOutputs,outputNames);
  funcIO->loadFunctionData(pAppFiles.nStrings_, pAppFiles.strArray_);
  return funcIO;
}

// ************************************************************************
// create an FunctionInterface instantiation without setting up the
// driver (in case it is a response surface, this will save the time to
// set it up)
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterfaceSimplified(PsuadeData *psuadeIO)
{
  int   nInputs, nOutputs;
  char  **inputNames, **outputNames;
  pData pPtr, pINames, pONames, pAppFiles;
  FunctionInterface *funcIO;

  // ----------------------------------------------------------------
  // set up the FunctionInterface (for sample runs)
  // ----------------------------------------------------------------
  funcIO = new FunctionInterface();
  assert(psuadeIO->getParameter("input_ninputs", pPtr) == 0);
  nInputs = pPtr.intData_;
  assert(psuadeIO->getParameter("input_names", pINames) == 0);
  inputNames = pINames.strArray_;
  assert(psuadeIO->getParameter("output_noutputs", pPtr) == 0);
  nOutputs = pPtr.intData_;
  assert(psuadeIO->getParameter("output_names", pONames) == 0);
  outputNames = pONames.strArray_;
  assert(psuadeIO->getParameter("app_files", pAppFiles) == 0);
  funcIO->loadInputData(nInputs, inputNames);
  funcIO->loadOutputData(nOutputs,outputNames);
  strcpy(pAppFiles.strArray_[0], "NULL");
  funcIO->loadFunctionData(pAppFiles.nStrings_, pAppFiles.strArray_);
  return funcIO;
}

// ************************************************************************
// create an FunctionInterface instantiation using data file for
// response surface
// ------------------------------------------------------------------------
FunctionInterface *createFunctionInterfaceGivenAppDriver(int nInputs,
                                               int nOutputs, char *fname)
{
  int   ii, nAppFiles=5;
  char  **inputNames, **outputNames, **appFiles;
  FunctionInterface *funcIO;

  appFiles = new char*[nAppFiles];
  for (ii = 0; ii < nAppFiles; ii++)
  {
    appFiles[ii] = new char[500];
    strcpy(appFiles[ii], "NULL");
  }
  strcpy(appFiles[0], fname);

  inputNames = new char*[nInputs];
  for (ii = 0; ii < nInputs; ii++)
  {
    inputNames[ii] = new char[500];
    strcpy(inputNames[ii], "XX");
  }

  outputNames = new char*[nOutputs];
  for (ii = 0; ii < nOutputs; ii++)
  {
    outputNames[ii] = new char[500];
    strcpy(outputNames[ii], "YY");
  }

  funcIO = new FunctionInterface();
  funcIO->loadInputData(nInputs, inputNames);
  funcIO->loadOutputData(nOutputs,outputNames);
  funcIO->loadFunctionData(nAppFiles, appFiles);

  for (ii = 0; ii < nAppFiles; ii++) delete [] appFiles[ii];
  delete [] appFiles;
  for (ii = 0; ii < nInputs; ii++) delete [] inputNames[ii];
  delete [] inputNames;
  for (ii = 0; ii < nOutputs; ii++) delete [] outputNames[ii];
  delete [] outputNames;

  return funcIO;
}

// ************************************************************************
// choose between different local functions 
// ------------------------------------------------------------------------
int FunctionInterface::setLocalFunction(int problem)
{
  whichLocalFunction_ = problem;
  if (printLevel_ > 0)
  { 
    if (problem == 0)
      printf("FunctionInterface setLocalFunc: default OUU example\n");
    else if (problem == 10)
      printf("FunctionInterface setLocalFunc: ODOE_GOPTIMAL\n");
    else if (problem == 11)
      printf("FunctionInterface setLocalFunc: ODOE_IOPTIMAL\n");
    else if (problem == 12)
      printf("FunctionInterface setLocalFunc: ODOE_DOPTIMAL\n");
    else if (problem == 13)
      printf("FunctionInterface setLocalFunc: ODOE_AOPTIMAL\n");
    else if (problem == 14)
      printf("FunctionInterface setLocalFunc: ODOE_EOPTIMAL\n");
    else if (problem == 20)
      printf("FunctionInterface setLocalFunc: ODOE_GOPTIMAL (Bayes)\n");
    else if (problem == 21)
      printf("FunctionInterface setLocalFunc: ODOE_IOPTIMAL (Bayes)\n");
    else if (problem == 22)
      printf("FunctionInterface setLocalFunc: ODOE_DOPTIMAL (Bayes)\n");
    else if (problem == 23)
      printf("FunctionInterface setLocalFunc: ODOE_AOPTIMAL (Bayes)\n");
    else if (problem == 24)
      printf("FunctionInterface setLocalFunc: ODOE_EOPTIMAL (Bayes)\n");
    else if (problem == 25)
      printf("FunctionInterface setLocalFunc: ODOE_GOPTIMAL (Fisher)\n");
    else if (problem == 26)
      printf("FunctionInterface setLocalFunc: ODOE_IOPTIMAL (Fisher)\n");
    else if (problem == 27)
      printf("FunctionInterface setLocalFunc: ODOE_DOPTIMAL (Fisher)\n");
    else if (problem == 28)
      printf("FunctionInterface setLocalFunc: ODOE_AOPTIMAL (Fisher)\n");
    else if (problem == 29)
      printf("FunctionInterface setLocalFunc: ODOE_EOPTIMAL (Fisher)\n");
    else if (problem == 30)
      printf("FunctionInterface setLocalFunc: ODOE_MMD\n");
    else if (problem == 999)
      printf("FunctionInterface setLocalFunction: clean up\n");
  }
  if (problem == 999)
  {
    psLocalFunction(0, NULL, 0, NULL);
    whichLocalFunction_ = 0;
  }
  return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FunctionInterface& FunctionInterface::operator=(const FunctionInterface &)
{
  printf("FunctionInterface operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

