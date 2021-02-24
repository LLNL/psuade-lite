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
// Functions for the class FORMAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FuncApprox.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "FORMAnalyzer.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FORMAnalyzer::FORMAnalyzer() : Analyzer()
{
   setName("FORM");
   rsType_ = PSUADE_RS_MARS; /* default is MARS */
   printOutTS(PL_ERROR, "FORMAnalyzer currently not supported.\n");
   exit(1);
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FORMAnalyzer::~FORMAnalyzer() 
{ 
} 

// ************************************************************************ 
// perform analysis 
// ------------------------------------------------------------------------
double FORMAnalyzer::analyze(aData &adata)
{
   int        nInputs, nOutputs, nSamples, outputID, sInd, ii, status;
   int        nPtsPerDim=64, wgtID, iter, printLevel, iOne=1;
   double     *X, *Y, *xLower, *xUpper, *YY, *wgts=NULL, *currZ, beta;
   double     ddata, gdata, Ldata, oldbeta, *alphas, thresh;
   FuncApprox *fa;

   printLevel = adata.printLevel_;
   nInputs    = adata.nInputs_;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   xLower     = adata.iLowerB_;
   xUpper     = adata.iUpperB_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   outputID   = adata.outputID_;
   wgtID      = adata.regWgtID_;
   thresh     = adata.analysisThreshold_;
   if (adata.inputPDFs_ != NULL)
   {
      printOutTS(PL_WARN, 
           "FORM INFO: some inputs have non-uniform PDFs, but\n");
      printOutTS(PL_WARN, 
           "           they will not be relevant in this analysis.\n");
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR, "FORMAnalyzer ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "    nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (printLevel > 0)
   {
      printAsterisks(PL_INFO, 0);
      printThisFA(rsType_);
      printEquals(PL_INFO, 0);
   }
   
   YY = new double[nSamples];
   checkAllocate(YY, "YY in FORM::analyze");
   for (sInd = 0; sInd < nSamples; sInd++)
      YY[sInd] = Y[sInd*nOutputs+outputID];
   if (wgtID >= 0)
   {
      wgts = new double[nSamples];
      checkAllocate(wgts, "wgts in FORM::analyze");
      for (sInd = 0; sInd < nSamples; sInd++)
         wgts[sInd] = Y[sInd*nOutputs+wgtID];
   }

   fa = genFA(rsType_, nInputs, iOne=1, nSamples);
   if (fa == NULL)
   {
      printOutTS(PL_INFO,"FORMAnalyzer ERROR: FuncApprox returns NULL.\n");
      delete [] YY;
      delete fa;
      if (wgts != NULL) delete [] wgts;
      return 1.0e12;
   }
   fa->setNPtsPerDim(nPtsPerDim);
   fa->setBounds(xLower, xUpper);
   fa->setOutputLevel(0);
   if (wgtID >= 0)
   {
      fa->loadWeights(nSamples, wgts);
      delete [] wgts;
   }
   status = fa->initialize(X, YY);
   if (status != 0)
   {
      printOutTS(PL_INFO, 
           "FORMAnalyzer ERROR: FuncApprox returns error.\n");
      delete [] YY;
      delete fa;
      return 1.0e12;
   }


   currZ = new double[nInputs];
   checkAllocate(currZ, "currZ in FORM::analyze");
   for (ii = 0; ii < nInputs; ii++)
      currZ[ii] = 0.5 * (xUpper[ii] + xLower[ii]);
   alphas = new double[nInputs];
   checkAllocate(alphas, "alphas in FORM::analyze");
   beta = 0.0;
   for (ii = 0; ii < nInputs; ii++) beta += currZ[ii] * currZ[ii];
   beta = sqrt(beta);

   iter = 1;
   while (iter < 1000)
   {
      gdata = fa->evaluatePoint(currZ) - thresh;
      for (ii = 0; ii < nInputs; ii++)
      {
         currZ[ii] += 1.0e-5;
         ddata = fa->evaluatePoint(currZ) - thresh;
         currZ[ii] -= 1.0e-5;
         alphas[ii] = (ddata - gdata) * 1.0e5;
      }
      Ldata = 0.0;
      for (ii = 0; ii < nInputs; ii++) Ldata += alphas[ii] * alphas[ii];
      Ldata = sqrt(Ldata);
      for (ii = 0; ii < nInputs; ii++) alphas[ii] /= Ldata;
      if (Ldata < 1.0e-12) break;
      for (ii = 0; ii < nInputs; ii++) 
         currZ[ii] = - alphas[ii] * (beta + gdata / Ldata);
      oldbeta = beta;
      beta = 0.0;
      for (ii = 0; ii < nInputs; ii++) beta += currZ[ii] * currZ[ii];
      beta = sqrt(beta);
      if (printLevel > 1)
         printOutTS(PL_INFO, 
              "FORM analyze: beta at iter %4d = %12.4e(%12.4e)\n",iter,
              beta, oldbeta);
      if (PABS((beta-oldbeta)/(beta+oldbeta)) < 1.0e-4) break;
      iter++;
   }
   if (printLevel > 1)
   {
      printOutTS(PL_INFO, "FORMAnalyzer: optimal point = \n");
      for (ii = 0; ii < nInputs; ii++) 
         printOutTS(PL_INFO, "\t X[%3d] = %12.4e\n",ii+1,currZ[ii]);
      printOutTS(PL_INFO, "optimal distance = %12.4e\n",beta);
      printOutTS(PL_INFO, "number of iterations = %d\n",iter);
   }

   delete fa;
   delete [] currZ;
   delete [] alphas;
   delete [] YY;
   return 0.0;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int FORMAnalyzer::setParams(int argc, char **argv)
{
   char  *request = (char *) argv[0];
   Analyzer::setParams(argc, argv);
   if (!strcmp(request, "rstype"))
   {
      if (argc != 2) printOutTS(PL_INFO,"FORMAnalyzer WARNING: setParams.\n");
      rsType_ = *(int *) argv[1];
      if (rsType_ < 0 || rsType_ >= PSUADE_NUM_RS) rsType_ = 0;
   }
   else
   {
      printOutTS(PL_ERROR, "FORMAnalyzer ERROR: setParams not valid.\n");
      exit(1);
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FORMAnalyzer& FORMAnalyzer::operator=(const FORMAnalyzer &)
{
   printOutTS(PL_ERROR, 
        "FORMAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

