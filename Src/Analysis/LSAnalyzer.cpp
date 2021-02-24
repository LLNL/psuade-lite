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
// functions for the class LSAnalyzer (Local sensitivity analysis)  
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LSAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "pData.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
LSAnalyzer::LSAnalyzer() : Analyzer(), nInputs_(0), lsMeasures_(0)
{
   setName("LSA");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
LSAnalyzer::~LSAnalyzer()
{
   if (lsMeasures_ != NULL) delete [] lsMeasures_;
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double LSAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID;
   int    ii, ss, whichOutput, ncount, index;
   double *X, *Y, *mEffect, dtemp, *xLower, *xUpper;
   char   pString[500], **iNames=NULL;
   FILE   *fp=NULL;
   pData  qData;
   PsuadeData *ioPtr;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   outputID = adata.outputID_;
   X        = adata.sampleInputs_;
   Y        = adata.sampleOutputs_;
   if (adata.inputPDFs_ != NULL)
   {
      ncount = 0;
      for (ii = 0; ii < nInputs; ii++) ncount += adata.inputPDFs_[ii];
      if (ncount > 0)
      {
         printOutTS(PL_INFO, 
              "LocalSA INFO: some inputs have non-uniform PDFs, but\n");
         printOutTS(PL_INFO, 
              "              they are not relevant in this analysis.\n");
      }
   }
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;
   ioPtr = adata.ioPtr_;
   if (ioPtr != NULL)
   {
      ioPtr->getParameter("input_names", qData);
      iNames = qData.strArray_;
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR, "LocalSA ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples != nInputs + 1)
   {
      printOutTS(PL_ERROR, 
           "LSAnalyzer ERROR: nSamples should be equal to nInputs+1.\n");
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 

   printOutTS(PL_INFO, "\n");
   printAsterisks(PL_INFO, 0);
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"* Local Sensitivity Analysis\n");
   printOutTS(PL_INFO,"* This analysis works if the model output can be\n");
   printOutTS(PL_INFO,"* approximated by linear combination of the inputs.\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO, "* total number of samples = %d\n",nSamples);
   printOutTS(PL_INFO, "* number of Inputs        = %d\n",nInputs);
   printOutTS(PL_INFO, "* Output number           = %d\n", whichOutput+1);
   printDashes(PL_INFO, 0);
   if (psPlotTool_ == 1)
   {
      fp = fopen("scilablsa.sci", "w");
      if (fp != NULL)
      {
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// * Local Sensitivity Analysis                ** \n");
         fprintf(fp,"// *-------------------------------------------** \n");
         fprintf(fp,"// Output %d\n", whichOutput+1);
         fprintf(fp,"// *-------------------------------------------** \n");
      }
   }
   else
   {
      fp = fopen("matlablsa.m", "w");
      if (fp != NULL)
      {
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% * Local Sensitivity Analysis                ** \n");
         fprintf(fp,"%% *-------------------------------------------** \n");
         fprintf(fp,"%% Output %d\n", whichOutput+1);
         fprintf(fp,"%% *-------------------------------------------** \n");
      }
   }

   if (lsMeasures_) delete [] lsMeasures_;
   nInputs_ = nInputs;
   lsMeasures_ = new double[nInputs];
   mEffect = new double[nInputs];
   checkAllocate(mEffect, "mEffect in LSAnalyzer::analyze");
   for (ii = 0; ii < nInputs; ii++) 
      mEffect[ii] = lsMeasures_[ii] = PSUADE_UNDEFINED;

   for (ss = 1; ss < nSamples; ss++)
   {
      ncount = 0;
      for (ii = 0; ii < nInputs; ii++)
      {
         if (PABS(X[nInputs*ss+ii]-X[ii]) > 1.0e-15) 
         {
            ncount++;
            index = ii;
         }
      }
      if (ncount == 1)
      {
         mEffect[index] = Y[nOutputs*ss+whichOutput] - Y[whichOutput];
         dtemp = X[nInputs*ss+index] - X[index];
         mEffect[index] /= dtemp; 
         mEffect[index] *= (xUpper[index] - xLower[index]); 
         if (mEffect[index] < 0.0) mEffect[index] = - mEffect[index];
      }
   } 
   for (ii = 0; ii < nInputs; ii++)
   {
      if (mEffect[ii] == PSUADE_UNDEFINED)
      {
         printOutTS(PL_INFO,"LocalSA ERROR: sample not suitable for LSA.\n");
         delete [] mEffect;
         // Bill Oliver added code to close fp prior to returning
         if(fp != NULL) fclose(fp);
         return -1;
      }
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      lsMeasures_[ii] = mEffect[ii];
      printOutTS(PL_INFO,"* Input %3d (importance measure = %12.4e)\n", 
                 ii+1, mEffect[ii]);
   }

   if (fp != NULL)
   {
      fprintf(fp, "sortFlag = 0;\n");
      fprintf(fp, "nn = %d;\n", nInputs);
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%24.16e \n", PABS(mEffect[ii]));
      fprintf(fp, "]; \n");
      if (iNames == NULL)
      {
         fprintf(fp, "  Str = {");
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
         fprintf(fp,"'X%d'};\n",nInputs);
      }
      else
      {
         fprintf(fp, "  Str = {");
         for (ii = 0; ii < nInputs-1; ii++)
         {
            if (iNames[ii] != NULL) fprintf(fp,"'%s',",iNames[ii]);
            else                    fprintf(fp,"'X%d',",ii+1);
         }
         if (iNames[nInputs-1] != NULL) 
              fprintf(fp,"'%s'};\n",iNames[nInputs-1]);
         else fprintf(fp,"'X%d'};\n",nInputs);
      }
      fprintf(fp, "if (sortFlag == 1);\n");
      if (psPlotTool_ == 1) fprintf(fp,"  [AA,II] = gsort(Y);\n");
      else                  fprintf(fp,"  [AA,II] = sort(Y, 'descend');\n");
      fprintf(fp, "  II  = II(1:nn);\n");
      fprintf(fp, "  Str = Str(II);\n");
      fprintf(fp, "  Y   = Y(II);\n");
      fprintf(fp, "end\n");
      fprintf(fp, "ymax = max(Y);\n");
      fprintf(fp, "ymin = 0.0;\n");
      fprintf(fp, "if (ymax == ymin)\n");
      fprintf(fp, "   ymax = ymax * 1.01;\n");
      fprintf(fp, "   ymin = ymin * 0.99;\n");
      fprintf(fp, "end;\n");
      fprintf(fp, "if (ymax == ymin)\n");
      fprintf(fp, "   ymax = ymax + 0.01;\n");
      fprintf(fp, "   ymin = ymin - 0.01;\n");
      fprintf(fp,"end;\n");
      fwritePlotCLF(fp);
      fprintf(fp, "bar(Y,0.8);\n");
      fwritePlotAxes(fp);
      if (psPlotTool_ == 1)
      {
         fprintf(fp,"a=gca();\n");
         fprintf(fp,"a.data_bounds=[0, ymin-0.01*(ymax-ymin); nn+1, ");
         fprintf(fp,"ymax+0.01*(ymax-ymin)];\n");
         fprintf(fp,"newtick = a.x_ticks;\n");
         fprintf(fp,"newtick(2) = [1:nn]';\n");
         fprintf(fp,"newtick(3) = Str';\n");
         fprintf(fp,"a.x_ticks = newtick;\n");
         fprintf(fp,"a.x_label.font_size = 3;\n");
         fprintf(fp,"a.x_label.font_style = 4;\n");
      }
      else
      {
         fprintf(fp,
            "axis([0 nn+1 ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])\n");
         fprintf(fp,"set(gca,'XTickLabel',[]);\n");
         fprintf(fp,"th=text(1:nn, repmat(ymin-0.07*(ymax-ymin),nn,1),Str,");
         fprintf(fp,"'HorizontalAlignment','left','rotation',90);\n");
         fprintf(fp,"set(th, 'fontsize', 12)\n");
         fprintf(fp,"set(th, 'fontweight', 'bold')\n");
      }
      sprintf(pString, "Linear Sensitivity Measures");
      fwritePlotTitle(fp, pString);
      sprintf(pString, "Sensitivity Measure");
      fwritePlotYLabel(fp, pString);
      fclose(fp);
      if (psPlotTool_ == 1)
           printOutTS(PL_INFO, "LocalSA ranking result in scilablsa.sci\n");
      else printOutTS(PL_INFO, "LocalSA ranking result in matlablsa.m\n");
   }
   printAsterisks(PL_INFO, 0);

   delete [] mEffect;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
LSAnalyzer& LSAnalyzer::operator=(const LSAnalyzer &)
{
   printOutTS(PL_ERROR,"LocalSA operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int LSAnalyzer::get_nInputs()
{
   return nInputs_;
}
double *LSAnalyzer::get_lsMeasures()
{
   double* retVal = NULL;
   if (lsMeasures_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in LSAnalyzer::get_lsMeasures");
      for (int ii = 0; ii < nInputs_; ii++)
         retVal[ii] = lsMeasures_[ii];
   }
   return retVal;
}

