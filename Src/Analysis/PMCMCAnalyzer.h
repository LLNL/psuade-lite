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
// Definition for the class PMCMCAnalyzer 
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#ifndef __PMCMCANALYZERH__
#define __PMCMCANALYZERH__

#include "pData.h"
#include "aData.h"
#include "FuncApprox.h"
#include "CommManager.h"

// ************************************************************************
// class definition
// ************************************************************************
class PMCMCAnalyzer 
{
private:
   char name_[100];
   int  mode_;
   int  nInputs_;
   int  nOutputs_;
   CommManager *commMgr_;
   int  mypid_;
   int  nprocs_;
   double *means_;
   double *sigmas_;

public:

   PMCMCAnalyzer(CommManager *);

   ~PMCMCAnalyzer();

   double analyze(aData &adata);

   PMCMCAnalyzer& operator=(const PMCMCAnalyzer &analyzer);

   double genMatlabFile(int nInputs, double *lower, double *upper,
                        double *ranges, int nPlots, int *plotIndices,
                        int nbins, int **bins, int ****bins2,
                        pData &pData, int nChains, int chainCnt,
                        double **XChains, int *chainStatus);

   int    genPostLikelihood(int nInputs, double *lower,
                 double *upper, double *XRange, int numChains,
                 int chainCnt, double **XChains, int *chainStatus,
                 int chainLimit, int *rsIndices, double *rsValues,
                 int *designParams, int dnInputs, int dnSample,
                 double *dSamInputs, FuncApprox **faPtrs,
                 FuncApprox **faPtrs1, int nOutputs, double *discOutputs,
                 double *discFuncConstantMeans, double *dSamMeans,
                 double *dSamStdevs);

   int setParams(int nParams, char **params);
};

#endif // __PMCMCANALYZERH__

