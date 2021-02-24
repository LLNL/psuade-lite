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
// Functions for the class BootstrapAnalyzer (bootstrap confidence interval) 
// ------------------------------------------------------------------------
// Reference: "Hypothesis Testing Procedures for Eliminating Variables of 
//            Negligible Impact in Large Scale Computer Simulations"
//            by Jason Lenderman
// ************************************************************************
// AUTHOR : CHARLES TONG
// DATE   : 2007
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "BootstrapAnalyzer.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
BootstrapAnalyzer::BootstrapAnalyzer() : Analyzer()
{
   setName("BOOTSTRAP");
   nSteps_ = 0;
   mean_   = 0.0;
   stdev_  = 0.0;
   storedValues_ = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
BootstrapAnalyzer::~BootstrapAnalyzer()
{
   if (storedValues_ != NULL) delete [] storedValues_;
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double BootstrapAnalyzer::analyze(aData &adata)
{
   int     nOutputs, nSamples, outputID, ss, nB=5000, nPass,ind,printLevel;
   double  *Y, *YY, *bsmeans, dtemp, mean, var, z, z0, aVal;
   double  alpha=0.05, alpha_u, alpha_t;

   printLevel = adata.printLevel_;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   Y          = adata.sampleOutputs_;
   outputID   = adata.outputID_;

   if (nOutputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR, "BootstrapAnalyzer ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (outputID < 0 || outputID >= nOutputs)
   {
      printOutTS(PL_ERROR, "BootstrapAnalyzer ERROR: invalid outputID.\n");
      printOutTS(PL_ERROR, "    outputID = %d\n", outputID);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples <= 1)
   {
      printOutTS(PL_ERROR, "BootstrapAnalyzer INFO:not meaningful to ");
      printOutTS(PL_ERROR, "do this when nSamples <= 1.\n");
      return PSUADE_UNDEFINED;
   } 
   nPass = 0;
   for (ss = 0; ss < nSamples; ss++)
      if (Y[nOutputs*ss+outputID] == PSUADE_UNDEFINED) nPass++;
   if (nPass > 0)
   {
      printOutTS(PL_ERROR, 
           "BootstrapAnalyzer ERROR: Some outputs are undefined.\n");
      printOutTS(PL_ERROR, 
           "                         Prune the undefined's first.\n");
      return PSUADE_UNDEFINED;
   }
   if (storedValues_ != NULL) delete [] storedValues_;
   storedValues_ = NULL;

   if (printLevel > 1)
   {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO, "*                   Bootstrap Analysis\n");
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO, "*  Basic Statistics*\n");
      printDashes(PL_INFO, 0);
      printOutTS(PL_INFO, "* Output of interest = %d\n", outputID+1);
      printDashes(PL_INFO, 0);
   }
   computeMeanVariance(nSamples,nOutputs,Y,&mean,&var,outputID,printLevel);

   YY = new double[nSamples];
   checkAllocate(YY, "YY in Bootstrap::analyze");
   for (ss = 0; ss < nSamples; ss++) YY[ss] = Y[nOutputs*ss+outputID];

   bsmeans = new double[nB];
   checkAllocate(bsmeans, "bsmeans in Bootstrap::analyze");
   genBootstrap(nSamples, YY, nB, bsmeans);

   nPass = 0;
   for (ss = 0; ss < nB; ss++) if (bsmeans[ss] < mean) nPass++;

   setupNormalCDF(0.0e0, 1.0e0);
   dtemp = (double) nPass / (double) nB;
   z0    = normalCDFInv(dtemp);

   // ---------------------------------------------------------------
   aVal = genJackknife(nSamples, YY);

   alpha_t = 1.0 - alpha;
   z = normalCDFInv(alpha_t);
   dtemp = z0 + (z0 + z) / (1.0 - aVal * (z0 + z));
   if      (dtemp <= -10.0) alpha_u = 0.0;
   else if (dtemp >= 10.0)  alpha_u = 1.0;
   else
   {
       ind = (int) (0.05 * (dtemp + 10.0) * (double) nSteps_); 
       alpha_u = storedValues_[ind];
   } 
   if (aVal == 0.0) ind = 0;
   else             ind = (int) (alpha_u * nB);
   if (printLevel > 2)
      printOutTS(PL_DETAIL, 
         "Bootstrap confidence interval (%3.1f%% coverage) = {0, %e}\n",
         100.0*(1.0 - alpha), bsmeans[ind]);
   if (printLevel > 1)
      printAsterisks(PL_INFO, 0);

   delete [] YY;
   dtemp = bsmeans[ind];
   delete [] bsmeans;
   return dtemp;
}

// *************************************************************************
// Compute mean and variance
// -------------------------------------------------------------------------
int BootstrapAnalyzer::computeMeanVariance(int nSamples, int yDim, 
              double *Y, double *ymean, double *yvar, int yID, int flag)
{
   int    ss;
   double mean, variance;

   mean = 0.0;
   for (ss = 0; ss < nSamples; ss++) mean += Y[yDim*ss+yID];
   mean /= (double) nSamples;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
      variance += ((Y[yDim*ss+yID] - mean) * (Y[yDim*ss+yID] - mean));
   variance /= (double) (nSamples - 1);
   (*ymean) = mean;
   (*yvar)  = variance;

   printOutTS(PL_INFO, "BootstrapAnalyzer::mean     = %16.8e\n", mean);
   printOutTS(PL_INFO, "BootstrapAnalyzer::variance = %16.8e\n", variance);

   return 0;
}

// *************************************************************************
// get a bootstrap sample
// -------------------------------------------------------------------------
int BootstrapAnalyzer::genBootstrap(int nSamples, double *Y, int ntimes,
                                    double *bmeans)
{
   int    ss, ii, ind;
   double mean;
   
   for (ss = 0; ss < ntimes; ss++)
   {
      mean = 0.0;
      for (ii = 0; ii < nSamples; ii++)
      {
         ind   = PSUADE_rand() % nSamples;
         mean += Y[ind];
      }
      mean /= (double) nSamples;
      bmeans[ss] = mean;
   }
   return 0;
}

// *************************************************************************
// generate the measure based on the Jackknife technique
// -------------------------------------------------------------------------
double BootstrapAnalyzer::genJackknife(int nSamples, double *Y)
{
   int    ss, ii;
   double mean, dtemp, numer, denom;
   
   mean = 0.0;
   for (ss = 0; ss < nSamples; ss++) mean += Y[ss];
   mean /= (double) nSamples;
   denom = numer = 0.0; 
   for (ss = 0; ss < nSamples; ss++)
   {
      dtemp = 0.0;
      for (ii = 0; ii < nSamples; ii++)
         if (ii != ss) dtemp += Y[ii];
      dtemp /= (double) (nSamples - 1);
      numer += pow((mean - dtemp), 3.0);
      denom += pow((mean - dtemp), 2.0);
   }
   denom = 6.0 * pow(denom, 1.5);
   if (denom == 0.0) return 0.0;
   return (numer/denom);
}

// *************************************************************************
// create a CDF for a normal distribution
// -------------------------------------------------------------------------
int BootstrapAnalyzer::setupNormalCDF(double mean, double stdev)
{
   int    ii;
   double coef, denom, hstep, left, right, value, xdata, expo;

   mean_   = mean;
   stdev_  = stdev;
   nSteps_ = 100000;
   coef    = 1.0 / (stdev * sqrt(2.0*M_PI));
   denom   = 1.0 / (2.0 * stdev * stdev);
   left    = mean - 10.0 * stdev;
   right   = mean + 10.0 * stdev;
   hstep   = (right - left) / (double) nSteps_;


   value = 0.0;
   storedValues_ = new double[nSteps_+1];
   checkAllocate(storedValues_, "storedValues in Bootstrap::setupNormalCDF");
   for (ii = 0; ii < nSteps_; ii++)
   {
      xdata  = left + hstep * ii;
      expo   = - (xdata - mean) * (xdata - mean) * denom;
      value += coef * exp(expo) * hstep;
      storedValues_[ii] = value;
   }
   storedValues_[nSteps_] = 1.0;
   return 0;
}

// *************************************************************************
// perform CDF reverse for normal distribution
// -------------------------------------------------------------------------
double BootstrapAnalyzer::normalCDFInv(double &value)
{
   int    index;
   double hstep;

   hstep   = 20.0 * stdev_ / (double) nSteps_;
   if (value <= 0.0) return (-10.0*stdev_); 
   if (value >= 1.0) return (10.0*stdev_); 
   index = binarySearchDble(value, storedValues_, nSteps_+1);
   return (mean_-10.0*stdev_+hstep*index);
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
BootstrapAnalyzer& BootstrapAnalyzer::operator=(const BootstrapAnalyzer &)
{
   printOutTS(PL_ERROR, 
         "BootstrapAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}
// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int BootstrapAnalyzer::get_nSteps()
{
   return nSteps_;
}
double BootstrapAnalyzer::get_mean()
{
   return mean_;
}
double BootstrapAnalyzer::get_stdev()
{
   return stdev_;
}
double *BootstrapAnalyzer::get_storedValues()
{
   double* retVal = NULL;
   if (storedValues_)
   {
      retVal = new double[nSteps_];
      checkAllocate(retVal, "retVal in Bootstrap::get_storedValues");
      for (int ii = 0; ii < nSteps_; ii++) retVal[ii] = storedValues_[ii];
   }
   return retVal;
}

