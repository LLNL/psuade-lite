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
// Functions for the class GowerAnalyzer (Gower's distance analysis) 
// ************************************************************************
// AUTHOR : CHARLES TONG
// DATE   : 2009
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GowerAnalyzer.h"
#include "MainEffectAnalyzer.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "Sampling.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
GowerAnalyzer::GowerAnalyzer() : Analyzer()
{
   setName("GOWER");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GowerAnalyzer::~GowerAnalyzer()
{
}

// ************************************************************************
// perform mean/variance analysis
// ------------------------------------------------------------------------
double GowerAnalyzer::analyze(aData &adata)
{
   int     nInputs, nSamples, nInputs2, nSamples2, ss, ss2, printLevel;
   int     ii, status, nOutputs, outputID, iZero=0, nLHSample=100000;
   int     nLHSSub=500, *S2, iOne=1;
   double  *X, *X2, *ranges, gower, dmax, dmin, *dmeans, *dvars, ddata;
   double  *Y, *Y2, *vvm, *mvv, *vvv, *vce, *lowerB, *upperB, vsum;
   char    dataFile[500], lineIn[500];
   FILE    *fp;
   PsuadeData *pIO;
   pData   pPtr;
   MainEffectAnalyzer *mePtr;
   FuncApprox         *faPtr;
   Sampling           *sampPtr;

   printLevel = adata.printLevel_;
   nSamples   = adata.nSamples_;
   nInputs    = adata.nInputs_;
   nOutputs   = adata.nInputs_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   outputID   = adata.outputID_;
   lowerB     = adata.iLowerB_;
   upperB     = adata.iUpperB_;

   Y2 = new double[nSamples];
   checkAllocate(Y2, "Y2 in Gower::analyze");
   for (ii = 0; ii < nSamples; ii++) Y2[ii] = Y[ii*nOutputs+outputID];

   faPtr = genFA(PSUADE_RS_MARS, nInputs, iOne, nSamples);
   status = faPtr->initialize(X, Y2);
   delete [] Y2;

   sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
   sampPtr->setInputBounds(nInputs, lowerB, upperB);
   sampPtr->setOutputParams(1);
   sampPtr->setSamplingParams(nLHSample, 200, 0);
   sampPtr->initialize(0);
   X2 = new double[nLHSample*nInputs];
   Y2 = new double[nLHSample];
   S2 = new int[nLHSample];
   checkAllocate(S2, "S2 in Gower::analyze");
   sampPtr->getSamples(nLHSample, nInputs, 1, X2, Y2, S2);
   faPtr->evaluatePoint(nLHSample, X2, Y2);

   mePtr = new MainEffectAnalyzer();
   vce = new double[nInputs];
   vvm = new double[nInputs];
   vvv = new double[nInputs];
   mvv = new double[nInputs];
   checkAllocate(mvv, "mvv in Gower::analyze");
   mePtr->computeVCE((int) nInputs, (int) nLHSample, (int) nLHSSub, 
                     (double *) X2, (double *) Y2, (int) iZero, 
                     (FILE *) NULL, (double *) vvm, (double *) mvv, 
                     (double *) vvv, (double *) vce);
   vsum = 0.0;
   for (ii = 0; ii < nInputs; ii++) vsum += vvm[ii];
   if (vsum <= 0.0)
   {
      printOutTS(PL_INFO,"GowerAnalyzer INFO: vce sum <= 0.0.\n");
      printOutTS(PL_INFO,"                    vce not used for scaling.");
      for (ii = 0; ii < nInputs; ii++) vvm[ii] = 1.0;
   }
   else
   {
      printOutTS(PL_INFO,"GowerAnalyzer INFO: vce used for scaling.\n");
      for (ii = 0; ii < nInputs; ii++) vvm[ii] /= vsum;
      for (ii = 0; ii < nInputs; ii++)
         printOutTS(PL_INFO,  "VCE %4d = %e\n", ii+1, vvm[ii]);
   } 
   delete [] X2;
   delete [] Y2;
   delete [] S2;
   delete [] vce;
   delete [] vvv;
   delete [] mvv;
   delete faPtr;
   delete sampPtr;
   delete mePtr;

   printf("Enter file name for unobserved data (outside convex hull): ");
   scanf("%s", dataFile);
   fp = fopen(dataFile,"r");
   if (fp == NULL)
   {
      printOutTS(PL_ERROR,"GowerAnalyzer ERROR: cannot open file %s.\n", 
                 dataFile);
      delete [] vvm;
      return PSUADE_UNDEFINED;
   }
   fclose(fp);
   fgets(lineIn,500,stdin);

   pIO = new PsuadeData();
   pIO->setOutputLevel(0);
   status = pIO->readPsuadeFile(dataFile);
   if (status != 0)
   {
      printf("ERROR: cannot read file %s in PSUADE format.\n",dataFile);
      exit(1);
   }
   pIO->getParameter("input_ninputs", pPtr);
   nInputs2 = pPtr.intData_;
   pIO->getParameter("method_nsamples", pPtr);
   nSamples2 = pPtr.intData_;
   pIO->getParameter("input_sample", pPtr);
   X2 = pPtr.dbleArray_;
   if (nInputs != nInputs2)
   {
      printOutTS(PL_ERROR,
           "GowerAnalyzer ERROR: different input dimensions %d %d\n",
           nInputs, nInputs2);
      delete [] vvm;
      delete pIO;
      return PSUADE_UNDEFINED;
   }
   ranges = new double[nInputs];
   checkAllocate(ranges, "ranges in Gower::analyze");
   for (ii = 0; ii < nInputs; ii++)
   {
      dmax = - PSUADE_UNDEFINED;
      dmin =   PSUADE_UNDEFINED;
      for (ss = 0; ss < nSamples; ss++)
      {
         if (X[ss*nInputs+ii] < dmin) dmin = X[ss*nInputs+ii];
         if (X[ss*nInputs+ii] > dmax) dmax = X[ss*nInputs+ii];
      }
      ranges[ii] = dmax - dmin;
      if (ranges[ii] == 0.0)
      {
         printOutTS(PL_ERROR,
            "GowerAnalyzer ERROR: some input range = 0 (%d).\n",ii+1);
         delete [] ranges;
         delete [] vvm;
         delete pIO;
         return PSUADE_UNDEFINED;
      }
   }

   if (psPlotTool_ == 1)
   {
      fp = fopen("psuade_gower_data.sci", "w");
      if (fp != NULL) fprintf(fp, "// Gower distance statistics\n");
   }
   else
   {
      fp = fopen("psuade_gower_data.m", "w");
      if (fp != NULL) fprintf(fp, "%% Gower distance statistics\n");
   }
   if (fp == NULL)
   {
      printOutTS(PL_ERROR,  
         "GowerAnalyzer ERROR: cannot write to psuade_gower_data.m\n");
      delete [] vvm;
      delete [] ranges;
      delete pIO; 
      return 1.0;
   }

   fprintf(fp, "G = [\n");
   for (ss = 0; ss < nSamples; ss++)
   {
      for (ss2 = 0; ss2 < nSamples2; ss2++)
      {
         gower = 0.0;
         for (ii = 0; ii < nInputs; ii++)
            gower += PABS(X[ss*nInputs+ii]-X2[ss2*nInputs+ii])/ranges[ii];
         gower /= (double) nInputs;
         fprintf(fp, "%e ", gower);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp,"];\n");
   fwritePlotCLF(fp);
   fprintf(fp,"for ii = 1 : %d\n", nSamples2);
   fprintf(fp,"   X = G(:,ii);\n");
   if (psPlotTool_ == 1)
      fprintf(fp,"   X = gsort(X,'g','i');\n");
   else
      fprintf(fp,"   X = sort(X);\n");
   fprintf(fp,"   Y = [1:%d]' / %d;\n", nSamples, nSamples);
   fprintf(fp,"   plot(X,Y)\n");
   if (psPlotTool_ == 0)
   {
      fprintf(fp,"   if (ii == 1)\n");
      fprintf(fp,"      hold on\n");
      fprintf(fp,"   end;\n");
   }
   fprintf(fp,"end;\n");
   fwritePlotAxes(fp);
   fwritePlotXLabel(fp, "Gower distance");
   fwritePlotYLabel(fp, "Probabilities");

   fprintf(fp, "\n");
   if (psPlotTool_ == 0)
   {
      fprintf(fp, "figure(2);\n");
      fprintf(fp, "%% Scaled Gower distance statistics\n");
   }
   else
   {
      fprintf(fp, "scf(2);\n");
      fprintf(fp, "// Scaled Gower distance statistics\n");
   }
   fprintf(fp, "G2 = [\n");
   for (ss = 0; ss < nSamples; ss++)
   {
      for (ss2 = 0; ss2 < nSamples2; ss2++)
      {
         gower = 0.0;
         for (ii = 0; ii < nInputs; ii++)
            gower += PABS(X[ss*nInputs+ii]-X2[ss2*nInputs+ii]) / 
                     ranges[ii] * vvm[ii];
         gower /= (double) nInputs;
         fprintf(fp, "%e ", gower);
      }
      fprintf(fp, "\n");
   }
   fprintf(fp,"];\n");
   fwritePlotCLF(fp);
   fprintf(fp,"for ii = 1 : %d\n", nSamples2);
   fprintf(fp,"   X = G2(:,ii);\n");
   if (psPlotTool_ == 1)
      fprintf(fp,"   X = gsort(X,'g','i');\n");
   else
      fprintf(fp,"   X = sort(X);\n");
   fprintf(fp,"   Y = [1:%d]' / %d;\n", nSamples, nSamples);
   fprintf(fp,"   plot(X,Y)\n");
   if (psPlotTool_ == 0)
   {
      fprintf(fp,"   if (ii == 1)\n");
      fprintf(fp,"      hold on\n");
      fprintf(fp,"   end;\n");
   }
   fprintf(fp,"end;\n");
   fwritePlotAxes(fp);
   fwritePlotXLabel(fp, "Scaled Gower distance");
   fwritePlotYLabel(fp, "Probabilities");

   fprintf(fp, "\n");
   if (psPlotTool_ == 0)
   {
      fprintf(fp,"figure(3);\n");
      fprintf(fp,"%% Mahalanobis distance statistics\n");
      fprintf(fp,"%% showing distance from centers of clusters/std dev\n");
   }
   else
   {
      fprintf(fp,"scf(3);\n");
      fprintf(fp,"// Mahalanobis distance statistics\n");
      fprintf(fp,"// showing distance from centers of clusters/std dev\n");
   }
   fprintf(fp, "M = [\n");
   dmeans = new double[nInputs];
   dvars  = new double[nInputs];
   checkAllocate(dvars, "dvars in Gower::analyze");
   for (ii = 0; ii < nInputs; ii++)
   {
      dmeans[ii] = 0.0;
      for (ss = 0; ss < nSamples; ss++) dmeans[ii] += X[ss*nInputs+ii];
      dmeans[ii] /= nSamples;
      dvars[ii] = 0.0;
      for (ss = 0; ss < nSamples; ss++) 
         dvars[ii] += pow(X[ss*nInputs+ii] - dmeans[ii], 2.0);
      dvars[ii] /= nSamples;
      if (dvars[ii] == 0.0)
         printf("GowerAnalyzer WARNING: input %d has zero variance.\n",ii);
   }
   for (ss2 = 0; ss2 < nSamples2; ss2++)
   {
      ddata = 0.0;
      for (ii = 0; ii < nInputs; ii++)
      {
         if (dvars[ii] != 0.0)
            ddata += pow(X2[ss2*nInputs+ii]-dmeans[ii], 2.0) / dvars[ii];
         else
         {
            ddata = PSUADE_UNDEFINED;
            break;
         }
      }
      if (ddata != PSUADE_UNDEFINED) ddata = sqrt(ddata);
      fprintf(fp, "%e\n", ddata);
   }
   fprintf(fp,"];\n");
   fwritePlotCLF(fp);
   fprintf(fp,"X = [1:%d]';\n", ss2);
   fprintf(fp,"plot(X,M,'b*')\n");
   fwritePlotAxes(fp);
   fwritePlotXLabel(fp, "Sample Number");
   fwritePlotYLabel(fp, "Mahalanobis distance");
   fclose(fp);
   if (psPlotTool_ == 1)
      printOutTS(PL_INFO,
           "GowerAnalyzer: psuade_gower_data.sci file created.\n");
   else
      printOutTS(PL_INFO,
           "GowerAnalyzer: psuade_gower_data.m file created.\n");

   delete [] ranges; 
   delete [] vvm;
   delete [] dmeans;
   delete [] dvars;
   delete pIO;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
GowerAnalyzer& GowerAnalyzer::operator=(const GowerAnalyzer &)
{
   printOutTS(PL_ERROR,
              "GowerAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

