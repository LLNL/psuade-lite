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
// Functions for the class IntegrationAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "IntegrationAnalyzer.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0) ? x : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
IntegrationAnalyzer::IntegrationAnalyzer() : Analyzer(), integral_(0)
{
   setName("INTEGRATION");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
IntegrationAnalyzer::~IntegrationAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double IntegrationAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, ss, ii, whichOutput;
   int    outputID, nLevels, *levelSeps, ncount;
   double result, last, error, *X, *Y, *iLower, *iUpper;

   nInputs   = adata.nInputs_;
   nOutputs  = adata.nOutputs_;
   nSamples  = adata.nSamples_;
   iLower    = adata.iLowerB_;
   iUpper    = adata.iUpperB_;
   X         = adata.sampleInputs_;
   Y         = adata.sampleOutputs_;
   outputID  = adata.outputID_;
   nLevels   = adata.currRefineLevel_;
   levelSeps = adata.refineSeparators_;
   if (adata.inputPDFs_ != NULL)
   {
      ncount = 0;
      for (ii = 0; ii < nInputs; ii++) ncount += adata.inputPDFs_[ii];
      if (ncount > 0)
      {
         printOutTS(PL_INFO,
              "Integration INFO: some inputs have non-uniform PDFs,\n");
         printOutTS(PL_INFO,
              "            but they are not relevant in this analysis\n");
      }
   }

   if (nInputs <= 0)
   {
      printOutTS(PL_ERROR,  "Integration ERROR: invalid nInputs.\n");
      printOutTS(PL_ERROR,  "    nInputs  = %d\n", nInputs);
      return PSUADE_UNDEFINED;
   } 
   if (nOutputs <= 0)
   {
      printOutTS(PL_ERROR,  "Integration ERROR: invalid nOutputs.\n");
      printOutTS(PL_ERROR,  "    nOutputs = %d\n", nOutputs);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples <= 0)
   {
      printOutTS(PL_ERROR,  "Integration ERROR: invalid nSamples.\n");
      printOutTS(PL_ERROR,  "    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   whichOutput = outputID;
   if (whichOutput < 0 || whichOutput >= nOutputs)
   {
      printOutTS(PL_ERROR,  "Integration ERROR: invalid outputID (%d).\n",
             whichOutput+1);
      return PSUADE_UNDEFINED;
   }
   for (ss = 0; ss < nSamples; ss++)
   {
      if (Y[ss] == PSUADE_UNDEFINED)
      {
         printOutTS(PL_ERROR,  
              "Integration ERROR: some outputs are undefined.\n");
         printOutTS(PL_ERROR,  
              "                   Prune them before analyze.\n");
         return PSUADE_UNDEFINED;
      }
   }

   result = 0.0;
   for (ss = 0; ss < nSamples; ss++) result += Y[ss];
   result /= (double) nSamples;
   for (ii = 0; ii < nInputs; ii++) result *= (iUpper[ii] - iLower[ii]);
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"Integration: numerical integral = %14.4e\n",result);
   printAsterisks(PL_INFO, 0);

   integral_ = result;

   if (nLevels <= 0) return result;

   last = 0.0;
   for (ss = 0; ss < levelSeps[nLevels-1]; ss++) last += Y[ss];
   last /= (double) levelSeps[nLevels-1];
   for (ii = 0; ii < nInputs; ii++) last *= (iUpper[ii] - iLower[ii]);

   if (result == 0.0) error = result;
   else               error = PABS((last-result)/result); 
   if (adata.printLevel_ > 0)
      printOutTS(PL_INFO,"Integration: numerical error    = %14.4e\n",error);

   return error;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
IntegrationAnalyzer& IntegrationAnalyzer::operator=(const IntegrationAnalyzer &)
{
   printOutTS(PL_ERROR,"Integration operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
double IntegrationAnalyzer::get_integral()
{
   return integral_;
}

