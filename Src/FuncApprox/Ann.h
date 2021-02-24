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
// Definitions for the class ANN
// AUTHOR : Christopher K. Chan
// DATE   : 2005
// ************************************************************************

#ifndef __ANNH__
#define __ANNH__

#include <stdio.h>
#include "FuncApprox.h"

// ************************************************************************
// class definition
// ************************************************************************
class Ann : public FuncApprox
{
public:

   /** constructor */
   Ann(int nInputs, int nSamples);

   /** Copy constructor prototype created by Bill Oliver */
   Ann(const Ann &);	

   /** destructor */
   ~Ann();

   int initialize(double *, double *);
   int genNDGridData(double *, double *, int *, double **, double **);
   int gen1DGridData(double *, double *, int, double *, int *, 
                     double **, double **);
   int gen2DGridData(double *, double *, int, int, double *, int *, 
                     double **, double **);
   int gen3DGridData(double *, double *, int, int, int, double *, int *, 
                     double **, double **);
   int gen4DGridData(double *, double *, int, int, int, int, double *, 
                     int *, double **, double **);
   double evaluatePoint(double *);
   double evaluatePoint(int, double *, double *);
   double evaluatePointFuzzy(double *, double &);
   double evaluatePointFuzzy(int, double *, double *, double *);

private:

   // added by Bill Oliver required by the Copy Constructor
   /* Topology Section */
   int nNet_;	
   int nHidden_;

   /* Data Section */
   int samplingMethod_;
   bool normFlag_;	
   float normBounds_[2];

   /* Training Section */
   char *learnFuncName_;
   int nLearnFuncParam_;
   float aLearnFuncParam_[28];
   int nMaxRetrain_;
   int accuracyLevel_;
   bool trainPerformed_;

   /* Normalization */
   double maxMagnitudeOutput_;
   double *inputDataMax_;
   double *inputDataMin_;

   /* Output Options */
   bool verboseFlag_;
   bool saveNetFlag_;
   bool savePatFlag_;

   /* Internal data */

   double scalingFactor_;

   /* Data arrays */
   double *inputData_;
	
   /* Private functions */
   void checkClassVars(void);
   void readConfigFile(void);
   void readTopologySection(FILE *fp);
   void readDataSection(FILE *fp);
   void readTrainingSection(FILE *fp);
   void readOutputSection(FILE *fp);
   void generateNet(double *X, double *Y);	
   void cleanUp();
};

#endif
