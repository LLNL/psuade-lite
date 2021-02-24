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
// Definition for the class MOATAnalyzer
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#ifndef __MOATANALYZERH__
#define __MOATANALYZERH__
#include <string>
#include <iostream>
#include <sstream>
#include "Analyzer.h"
#include "Vector.h"

// ************************************************************************
// class definition
// ************************************************************************
class MOATAnalyzer : public Analyzer
{

public:

   MOATAnalyzer();

   ~MOATAnalyzer();

   void analyze(int nInps, int nSamp, double *ilbs,
                double *iubs, double *sInps, double *sOuts);

   double analyze(aData &adata);

   MOATAnalyzer& operator=(const MOATAnalyzer &analyzer);

   /** Getters for analysis results */
   double get_mean(int ind);
   double get_stdev(int ind);
   double get_modifiedMean(int ind);
   double get_modifiedStdev(int ind);
 
private:
   int createScreenDiagramFile(int nSamples, int nInputs, double *Y, 
                               int *indices, double *modifiedMeans,
                               double *stds, int outputID, char **);
   int createScatterFile(int nSamples, int nInputs, double *Y, double *X,
                         int *indices, char **);
   int createBootstrapFile(int nSamples, int nInputs, double *Y, double *X,
                           int *indices, char **);

   /** All arrays are of length nInputs */
   int    nInputs_;
   int    nOutputs_;
   int    nSamples_;
   int    outputID_;
   psVector  means_;
   psVector  stds_;
   psVector  modifiedMeans_;
   psVector  modifiedStds_;
   psIVector indexesSortedByModifiedMeans_;

};

#endif // __MOATANALYZERH__

