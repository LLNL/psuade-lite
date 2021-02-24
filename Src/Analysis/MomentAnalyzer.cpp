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
// Functions for the class MomentAnalyzer  
// Reference: from any introductory statistics book
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "MomentAnalyzer.h"
#include "Psuade.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "pData.h"
#include "RSConstraints.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MomentAnalyzer::MomentAnalyzer() : Analyzer(), nSamples_(0), nGroups_(0),
                                nInputs_(0), nOutputs_(0), moments_(0)
{
   setName("MOMENT");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MomentAnalyzer::~MomentAnalyzer()
{
   if (moments_) delete[] moments_;
}

// ************************************************************************
// perform analysis (library call mode)
// ------------------------------------------------------------------------
void MomentAnalyzer::analyze(int nInps, int nSamp, double *lbs,
                             double *ubs, double *X, double *Y)
{
   aData adata;
   adata.nInputs_ = nInps;
   adata.nOutputs_ = 1;
   adata.nSamples_ = nSamp;
   adata.iLowerB_ = lbs;
   adata.iUpperB_ = ubs;
   adata.sampleInputs_ = X;
   adata.sampleOutputs_ = Y;
   adata.outputID_ = 0;
   adata.printLevel_ = 0;
   analyze(adata);
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double MomentAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, nSubSamples, status, nValid;
   int    nGroups, groupID, whichOutput, ss, printLevel, nConstr, cnt, outID;
   int    errFlag, ii;
   double *gMeans, *gVariances, mean, variance, skewness, kurtosis, *X;
   double *Y2, *YY, *Y=NULL, variance2;
   char   winput[500], filename[500], pString[500];
   FILE   *fp=NULL;
   RSConstraints *constrPtr;
   PsuadeData    *ioPtr;
   pData         pPtr;

   // Providing initialization by Bill Oliver
   outID = 0;
   nValid = 1;  // nValid is set to 1 to avoid a divide by zero later

   printLevel = adata.printLevel_;
   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   X        = adata.sampleInputs_;
   YY       = adata.sampleOutputs_;
   outputID = adata.outputID_;
   nSubSamples = adata.nSubSamples_;
   ioPtr       = adata.ioPtr_;
   if (adata.inputPDFs_ != NULL)
   {
      cnt = 0;
      for (ii = 0; ii < nInputs; ii++) cnt += adata.inputPDFs_[ii];
      if (cnt > 0 && isScreenDumpModeOn())
      {
         printOutTS(PL_INFO, 
            "MomentAnalyzer INFO: non-uniform probability distributions\n");
         printOutTS(PL_INFO,
            "               have been defined in the data file, but\n");
         printOutTS(PL_INFO, 
            "               they will not be used in this analysis.\n");
      }
   }
   if (ioPtr != NULL) ioPtr->getParameter("output_names", pPtr);

   if (nInputs <= 0 || nOutputs <= 0)
   {
      printOutTS(PL_ERROR,
           "MomentAnalyzer ERROR: invalid nInputs or nOutputs.\n");
      printOutTS(PL_ERROR,"   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR,"   nOutputs = %d\n", nOutputs);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples <= 1)
   {
      printOutTS(PL_ERROR,"MomentAnalyzer ERROR: nSamples should be > 1.\n");
      printOutTS(PL_ERROR,"    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (psScreenOutput_ == 1 && nSubSamples <= 0)
   {
      printOutTS(PL_INFO,"MomentAnalyzer INFO: nSubSamples = 0.\n");
      printOutTS(PL_INFO,"    This affects additional group analysis only.\n");
      printOutTS(PL_INFO,"    Set nSubSamples to nSamples.\n");
      nSubSamples = nSamples;
   } 
   if (isScreenDumpModeOn() && (nSamples/nSubSamples*nSubSamples != nSamples))
   {
      printOutTS(PL_INFO,"MomentAnalyzer INFO: nSamples != k*nSubSamples.\n");
      printOutTS(PL_ERROR,"   nSamples    = %d\n", nSamples);
      printOutTS(PL_ERROR,"   nSubSamples = %d\n", nSubSamples);
      printOutTS(PL_INFO,"    Set nSubSamples to nSamples.\n");
      nSubSamples = nSamples;
   } 
   whichOutput = outputID;
   if (whichOutput < -1 || whichOutput >= nOutputs)
   {
      printOutTS(PL_ERROR,"MomentAnalyzer ERROR: invalid outputID (%d).\n",
                whichOutput+1);
      printOutTS(PL_ERROR, "   outputID = %d\n", outputID+1);
      return PSUADE_UNDEFINED;
   } 
   if (moments_) delete[] moments_;
   moments_ = NULL;
 
   if (whichOutput < 0) whichOutput = 0;
   errFlag = 0;
   for (ss = 0; ss < nSamples; ss++)
   {
      for (ii = 0; ii < nOutputs; ii++)
      {
         if (YY[ss*nOutputs+ii] == PSUADE_UNDEFINED)
         {
            errFlag++;
            break;
         }
      }
   }
   if (errFlag == 0) Y = YY;
   else
   {
      printOutTS(PL_INFO,"MomentAnalyzer INFO: Some outputs are undefined,\n");
      printOutTS(PL_INFO,
           "               which are purged before processing.\n");
      Y = new double[nSamples*nOutputs];
      checkAllocate(Y, "Y in Moment::analyze");

      cnt = 0;
      for (ss = 0; ss < nSamples; ss++)
      {
         for (ii = 0; ii < nOutputs; ii++)
            if (YY[ss*nOutputs+ii] == PSUADE_UNDEFINED) break;
         if (ii == nOutputs)
         {
            for (ii = 0; ii < nOutputs; ii++)
               Y[cnt*nOutputs+ii] = YY[ss*nOutputs+ii];
            cnt++;
         } 
      }
      if (nSubSamples == nSamples) nSubSamples = cnt;
      if (isScreenDumpModeOn())
         printOutTS(PL_INFO,"MomentAnalyzer INFO: original sample size = %d\n",
                    nSamples);
      nSamples = cnt;
      if (isScreenDumpModeOn())
         printOutTS(PL_INFO,"MomentAnalyzer INFO: purged   sample size = %d\n",
                    nSamples);
   }
   if (nSamples <= 0) 
   {
      printOutTS(PL_INFO,"INFO: Sample size = 0 ==> no analysis\n");
      return -1;
   }
   
   nConstr = 0;
   constrPtr = NULL;
   if (ioPtr != NULL)
   {
      constrPtr = new RSConstraints();
      constrPtr->genConstraints(ioPtr);
      nConstr = constrPtr->genConstraints(ioPtr);
   }
   if (nConstr > 0)
   {
      if (outputID < 0 || outputID >= nOutputs)
      {
         printOutTS(PL_INFO,
              "MomentAnalyzer: Generate output distribution plot\n");
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
         outID = getInt(1, nOutputs, pString);
         outID--;
      }
      else outID = outputID;
      if (psPlotTool_ == 1) strcpy(filename, "scilabua.sci");
      else                  strcpy(filename, "matlabua.m");
      fp = NULL;
      if (psScreenOutput_ == 1)
      {
         fp = fopen(filename, "w");
         if (fp == NULL)
            printOutTS(PL_ERROR, "ERROR: cannot open file %s.\n", filename);
      }
      Y2 = new double[nSamples];
      checkAllocate(Y2, "Y2 in Moment::analyze");
      nValid = 0;
      if (fp != NULL) fprintf(fp, "Y = [\n");
      for (ss = 0; ss < nSamples; ss++)
      {
         constrPtr->evaluate(&X[ss*nInputs],Y[ss*nOutputs+outputID],status);
         if (status == 1)
         {
            Y2[nValid++] = Y[ss*nOutputs+outputID];
            if (fp != NULL) fprintf(fp, "  %24.16e \n",Y2[nValid-1]);
         }
      }
      if (nValid == 0)
      {
         printOutTS(PL_INFO,"MomentAnalyzer INFO: no valid data points.\n");
         delete [] Y2;
         delete constrPtr;
         if (fp != NULL) fclose(fp);
         if (errFlag != 0 && Y != NULL) delete [] Y;
         return 0.0;
      }
      if (fp != NULL)
      {
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1)
         {
            fwritePlotCLF(fp);
            fprintf(fp, "ymin = min(Y);\n");
            fprintf(fp, "ymax = max(Y);\n");
            fprintf(fp, "ywid = 0.1 * (ymax - ymin);\n");
            fprintf(fp, "if (ywid < 1.0e-12)\n");
            fprintf(fp, "   disp('range too small.')\n");
            fprintf(fp, "   halt\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "histplot(10, Y/ywid, style=2);\n");
            fprintf(fp, "a = gce();\n");
            fprintf(fp, "a.children.fill_mode = \"on\";\n");
            fprintf(fp, "a.children.thickness = 2;\n");
            fprintf(fp, "a.children.foreground = 0;\n");
            fprintf(fp, "a.children.background = 2;\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Output Distribution");
            if (ioPtr != NULL)
                 fwritePlotXLabel(fp, pPtr.strArray_[outID]);
            else fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
         }
         else
         {
            fwritePlotCLF(fp);
            fprintf(fp, "[nk,xk]=hist(Y(:,1),10);\n");
            fprintf(fp, "NB = nk;\n");
            fprintf(fp, "N  = nk;\n");
            fprintf(fp, "X  = xk;\n");
            fprintf(fp, "bar(xk,nk/%d,1.0)\n",nValid);
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Output Distribution");
            if (ioPtr != NULL)
                 fwritePlotXLabel(fp, pPtr.strArray_[outID]);
            else fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
         }
         printOutTS(PL_INFO,"Output distribution plot is now in %s.\n",
                    filename);
         fclose(fp);
      }
      computeMeanVariance(nInputs,1,nValid,Y2,&mean,&variance,0);
      if (psScreenOutput_ == 1)
      {
         printAsterisks(PL_INFO, 0);
         printOutTS(PL_INFO, "*       Sample mean          = %12.4e\n",mean);
         printOutTS(PL_INFO, "*       Sample std dev       = %12.4e\n",
                    sqrt(variance));
         printOutTS(PL_INFO, "*       Based on nSamples    = %d\n",nValid);
         printOutTS(PL_INFO, "*       Output number        = %d\n", outID+1);
      }
      delete [] Y2;
      delete constrPtr;
      if (errFlag != 0 && Y != NULL) delete [] Y;
      return sqrt(variance/nValid);
   }

   if (nSubSamples <= 0) nSubSamples = nSamples;
   nGroups    = nSamples / nSubSamples;
   gMeans     = new double[nGroups];
   gVariances = new double[nGroups];
   checkAllocate(gVariances, "gVariances in Moment::analyze");

   if (isScreenDumpModeOn())
   {
      printOutTS(PL_INFO, "\n");
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO, "*             Basic Output Statistics\n");
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO, "* nSamples = %10d\n",nSamples);
      printOutTS(PL_INFO, "* nGroups  = %10d\n",nGroups);
      printOutTS(PL_INFO, "* nInputs  = %10d\n",nInputs);
      printOutTS(PL_INFO, "* nOutputs = %10d\n",nOutputs);
      printDashes(PL_INFO, 0);
   }
   if (outputID == -1)
   {
      for (ss = 0; ss < nOutputs; ss++)
      {
         if (isScreenDumpModeOn())
            printOutTS(PL_INFO, "* outputID = %10d\n", ss+1);
         computeMeanVariance(nInputs,nOutputs,nSubSamples, Y,
                             gMeans,gVariances,ss);
         if (isScreenDumpModeOn())
         {
            printOutTS(PL_INFO, "*       Sample mean          = %12.4e\n",
                       gMeans[0]);
            printOutTS(PL_INFO, "*       Sample std dev       = %12.4e\n",
                       sqrt(gVariances[0]));
         }
         if (gVariances[0] > 0)
         {
            computeSkewnessKurtosis(nInputs,nOutputs,nSubSamples, Y, 
                         sqrt(gVariances[0]), &skewness,&kurtosis,ss);
            if (isScreenDumpModeOn())
            {
               printOutTS(PL_INFO,"*       Sample skewness      = %12.4e\n", 
                          skewness);
               printOutTS(PL_INFO,"*       Sample kurtosis      = %12.4e\n", 
                       kurtosis);
            }
         }
         else
         {
            if (isScreenDumpModeOn())
               printOutTS(PL_INFO,
                    "*       Std dev=0, skeweness/kurtosis set to 0.\n");
            skewness = kurtosis = 0.0;
         }
         if (isScreenDumpModeOn() && ss < nOutputs-1 && printLevel >= 0) 
            printEquals(PL_INFO, 0);
      }
      variance = gVariances[0];
   }
   else
   {
      if (isScreenDumpModeOn())
         printOutTS(PL_INFO,"* outputID = %10d\n", outputID+1);
      computeMeanVariance(nInputs,nOutputs,nSubSamples,Y,gMeans,
                          gVariances,outputID);
      if (isScreenDumpModeOn())
      {
         printOutTS(PL_INFO,"*       Sample mean          = %12.4e\n",
                    gMeans[0]);
         printOutTS(PL_INFO,"*       Sample std dev       = %12.4e\n",
                    sqrt(gVariances[0]));
      }
      if (gVariances[0] > 0)
      {
         computeSkewnessKurtosis(nInputs,nOutputs,nSubSamples, Y, 
                      sqrt(gVariances[0]), &skewness,&kurtosis,outputID);
         if (isScreenDumpModeOn())
         {
            printOutTS(PL_INFO,"*       Sample skewness      = %12.4e\n", 
                       skewness);
            printOutTS(PL_INFO,"*       Sample kurtosis      = %12.4e\n", 
                       kurtosis);
         }
      }
      else
      {
         if (isScreenDumpModeOn())
            printOutTS(PL_INFO,
                 "*       Std dev=0, skeweness/kurtosis set to 0.\n");
         skewness = kurtosis = 0.0;
      }
      variance = gVariances[0];

      pData *pPtr=NULL;
      if (ioPtr != NULL)
      {
         pPtr = ioPtr->getAuxData();
         if (pPtr != NULL)
         {
            pPtr->nDbles_ = 4;
            pPtr->dbleArray_ = new double[4];
            pPtr->dbleArray_[0] = gMeans[0];
            pPtr->dbleArray_[1] = sqrt(gVariances[0]);
            pPtr->dbleArray_[2] = skewness;
            pPtr->dbleArray_[3] = kurtosis;
         }
      }

   }
   //save moments
   moments_ = new double[4];
   moments_[0] = gMeans[0];
   moments_[1] = sqrt(gVariances[0]);
   moments_[2] = skewness;
   moments_[3] = kurtosis;
   if (outputID == -1) 
      printf("INFO: calls to get moments return only the last output info.\n");
   if (isScreenDumpModeOn() && printLevel >= 0)
   {
      printDashes(PL_INFO, 0);
      printAsterisks(PL_INFO, 0);
   }

#ifdef HAVE_PYTHON
   PyObject *temp, *Ylist;

   PyDict_SetItemString( AnalysisDataDict, "nSamples",
			 temp=PyInt_FromLong(nSamples) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nGroups",
			 temp=PyInt_FromLong(nGroups) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nInputs",
			 temp=PyInt_FromLong(nInputs) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "nOutputs",
			 temp=PyInt_FromLong(nOutputs) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "outputID",
			 temp=PyInt_FromLong(whichOutput) ); Py_DECREF(temp);

   PyDict_SetItemString( AnalysisDataDict, "mean",
			 temp=PyFloat_FromDouble(gMeans[0]) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "std dev",
			 temp=PyFloat_FromDouble(sqrt(gVariances[0])) );
   Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "skewness",
			 temp=PyFloat_FromDouble(skewness) ); Py_DECREF(temp);
   PyDict_SetItemString( AnalysisDataDict, "kurtosis",
			 temp=PyFloat_FromDouble(kurtosis) ); Py_DECREF(temp);

   Ylist = PyList_New(0);
   for (ss = 0; ss < nSamples; ss++)
      //if (YG[ii] != PSUADE_UNDEFINED) {
         PyList_Append( Ylist,
			temp=PyFloat_FromDouble(Y[nOutputs*ss+outputID]) );
	 Py_DECREF(temp);
      //}
   PyDict_SetItemString( AnalysisDataDict, "Y",	Ylist ); Py_DECREF(Ylist);
#endif

   if (outputID < 0 || outputID >= nOutputs)
   {
      if (psInteractive_ == 1)
      {
         printf("MomentAnalyzer: Generating output distribution plot\n");
         sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
         outID = getInt(1, nOutputs, pString);
         outID--;
      }
      else
      {
         outID = 0;
         printOutTS(PL_INFO,"Invalid output ID. Set to default output 1.\n");
      }
   }
   else outID = outputID;
   if (psPlotTool_ == 1) strcpy(filename, "scilabua.sci");
   else                  strcpy(filename, "matlabua.m");
   fp = NULL;
   if (psScreenOutput_ == 1) fp = fopen(filename, "w");
   if (fp != NULL)
   {
      fprintf(fp, "Y = [\n");
      for (ss = 0; ss < nSamples; ss++)
         fprintf(fp, "  %24.16e\n", Y[nOutputs*ss+outID]);  
      fprintf(fp, "];\n");
      if (psPlotTool_ == 1)
      {
         fwritePlotCLF(fp);
         fprintf(fp, "ymin = min(Y);\n");
         fprintf(fp, "ymax = max(Y);\n");
         fprintf(fp, "ywid = 0.1 * (ymax - ymin);\n");
         fprintf(fp, "if (ywid < 1.0e-12)\n");
         fprintf(fp, "   disp('range too small.')\n");
         fprintf(fp, "   halt\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "histplot(10, Y/ywid, style=2);\n");
         fprintf(fp, "a = gce();\n");
         fprintf(fp, "a.children.fill_mode = \"on\";\n");
         fprintf(fp, "a.children.thickness = 2;\n");
         fprintf(fp, "a.children.foreground = 0;\n");
         fprintf(fp, "a.children.background = 2;\n");
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Output Distribution");
         if (ioPtr != NULL)
              fwritePlotXLabel(fp, pPtr.strArray_[outID]);
         else fwritePlotXLabel(fp, "Output Values");
         fwritePlotYLabel(fp, "Probabilities");
      }
      else
      {
         fwritePlotCLF(fp);
         fprintf(fp, "twoPlots = 0;\n");
         fprintf(fp, "if (twoPlots == 1)\n");
         fprintf(fp, "subplot(1,2,1)\n");
         fprintf(fp, "end;\n");
         if (nSamples > 500) fprintf(fp, "[nk,xk]=hist(Y(:,1),20);\n");
         else                fprintf(fp, "[nk,xk]=hist(Y(:,1),10);\n");
         fprintf(fp, "bar(xk,nk/%d,1.0)\n",nSamples);
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Probability Distribution");
         if (ioPtr != NULL)
              fwritePlotXLabel(fp, pPtr.strArray_[outID]);
         else fwritePlotXLabel(fp, "Output Value");
         fwritePlotYLabel(fp, "Probabilities");
         fprintf(fp, "if (twoPlots == 1)\n");
         fprintf(fp, "Yk = sort(Y(:,1));\n");
         fprintf(fp, "Xk = 1 : %d;\n", nSamples);
         fprintf(fp, "Xk = Xk / %d;\n", nSamples);
         fprintf(fp, "subplot(1,2,2)\n");
         fprintf(fp, "plot(Yk, Xk, 'lineWidth',3)\n");
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Cumulative Distribution");
         if (ioPtr != NULL)
              fwritePlotXLabel(fp, pPtr.strArray_[outID]);
         else fwritePlotXLabel(fp, "Output Values");
         fwritePlotYLabel(fp, "Probabilities");
         fprintf(fp, "end;\n");
      }
      printOutTS(PL_INFO,"Output distribution plot is now in %s.\n",filename);
      fclose(fp);
   }
   else 
   {
      if (psScreenOutput_ == 1)
         printOutTS(PL_ERROR,"ERROR: cannot open file %s.\n", filename);
   }

   if (outputID >= 0)
   {
      if (nGroups == 1)
      {
         if (isScreenDumpModeOn()) printAsterisks(PL_INFO, 0);
         variance2 = gVariances[0];
         delete [] gMeans;
         delete [] gVariances;
         if (isScreenDumpModeOn() && printLevel >= 0)
         {
            printOutTS(PL_INFO,"*       std error of mean       = %16.8e\n",
                   sqrt(variance2/(double) nSamples));
            printAsterisks(PL_INFO, 0);
         }
         if (errFlag != 0 && Y != NULL) delete [] Y;
         delete constrPtr;
         return sqrt(variance2/(double) nSamples);
      }
      for (groupID = 0; groupID < nGroups; groupID++)
      {
         computeMeanVariance(nInputs,nOutputs,nSubSamples,
                &Y[nSubSamples*groupID*nOutputs],
                &gMeans[groupID],&gVariances[groupID],whichOutput);
      }
      if (isScreenDumpModeOn())
         printOutTS(PL_INFO, "Output %d\n", whichOutput+1);

      mean = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         mean += gVariances[groupID];
      mean = mean / (double) nGroups;
      variance2 = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         variance2 += ((gVariances[groupID] - mean) *
                      (gVariances[groupID] - mean));
      variance2 = variance2 / (double) nGroups;

      if (isScreenDumpModeOn())
      {
         printOutTS(PL_INFO, "*       Mean of group variances = %16.8e\n",mean);
         printOutTS(PL_INFO, "*       Std error of variance   = %16.8e\n",
                    variance2);
      }
      mean = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         mean += gMeans[groupID];
      mean = mean / (double) nGroups;
      variance2 = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         variance2 += ((gMeans[groupID] - mean) * (gMeans[groupID] - mean));
      variance2 = variance2 / (double) nGroups;
      if (isScreenDumpModeOn())
      {
         printOutTS(PL_INFO, "*       Mean of group means     = %16.8e\n",mean);
         printOutTS(PL_INFO, "*       Std error of mean       = %16.8e\n", 
                    variance2);
         printAsterisks(PL_INFO, 0);
      }
   }

   if (isScreenDumpModeOn() && printLevel > 1 && outputID >= 0)
      analyzeMore(nInputs, nOutputs, nSamples, nSubSamples, Y, outputID);

   delete [] gMeans;
   delete [] gVariances;
   delete constrPtr;
   if (errFlag != 0 && Y != NULL) delete [] Y;
   return sqrt(variance/nValid);
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int MomentAnalyzer::computeMeanVariance(int nInputs, int nOutputs, 
                        int nSamples, double *Y, double *mean, 
                        double *variance, int outputID)
{
   int    ss;
   double meanTmp, varTmp;

   meanTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++) meanTmp += Y[nOutputs*ss+outputID];
   meanTmp /= (double) nSamples;
   varTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      varTmp += ((Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp));
   }
   varTmp /= (double) (nSamples - 1);
   (*mean)     = meanTmp;
   (*variance) = varTmp;
   return 0;
}

// *************************************************************************
// Compute skewness and kurtosis 
// kurtois of a standard normal distribution is 3.
// kurtois-3 measures excess kurtosis
// -------------------------------------------------------------------------
int MomentAnalyzer::computeSkewnessKurtosis(int nInputs, int nOutputs, 
                        int nSamples, double *Y, double stdev,
                        double *skewness, double *kurtosis, int outputID)
{
   int    ss;
   double meanTmp, skewTmp, kurTmp;

   meanTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++) meanTmp += Y[nOutputs*ss+outputID];
   meanTmp /= (double) nSamples;
   skewTmp = kurTmp = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      skewTmp += ((Y[nOutputs*ss+outputID] - meanTmp) *
                  (Y[nOutputs*ss+outputID] - meanTmp) *
                  (Y[nOutputs*ss+outputID] - meanTmp));
      kurTmp += ((Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp) *
                 (Y[nOutputs*ss+outputID] - meanTmp));
   }
   skewTmp /= (double) (nSamples - 1);
   kurTmp  /= (double) (nSamples - 1);
   if (stdev >= 0.0) (*skewness) = skewTmp / (stdev * stdev * stdev);
   else              (*skewness) = 0.0;
   if (stdev >= 0.0) (*kurtosis) = kurTmp / (stdev * stdev * stdev * stdev);
   else              (*kurtosis) = 0.0;
   return 0;
}

// ************************************************************************
// perform more analysis
// ------------------------------------------------------------------------
int MomentAnalyzer::analyzeMore(int nInputs, int nOutputs, int nSamples,
                              int nSubSamples, double *Y, int outputID)
{
   int    nGroups, groupID, myNSubSamples, nLevels=10, level, multiple;
   double *gMeans, *gVariances, mean, variance;

   for (level = 1; level <= nLevels; level++)
   {
      multiple = 1 << level;
      nGroups = nSamples / nSubSamples;
      if ((nGroups % multiple) != 0) break;
      myNSubSamples = nSubSamples * multiple;
      nGroups    = nSamples / myNSubSamples;
      gMeans     = new double[nGroups];
      gVariances = new double[nGroups];
      checkAllocate(gVariances, "gVariances in Moment::analyzeMore");


      for (groupID = 0; groupID < nGroups; groupID++)
      {
         computeMeanVariance(nInputs,nOutputs,myNSubSamples,
                &Y[myNSubSamples*groupID*nOutputs],
                &gMeans[groupID],&gVariances[groupID],outputID);
      }
      mean = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         mean += gMeans[groupID];
      mean = mean / (double) nGroups;
      variance = 0.0;
      for (groupID = 0; groupID < nGroups; groupID++)
         variance += ((gMeans[groupID] - mean) * (gMeans[groupID] - mean));
      variance = variance / (double) nGroups;
      printEquals(PL_INFO, 0);
      printOutTS(PL_INFO,"*       Agglomerate %d groups into 1 \n",multiple);
      printOutTS(PL_INFO,"*       Mean of group means     = %16.8e\n",mean);
      printOutTS(PL_INFO,"*       Std error of mean       = %16.8e\n",variance);
      printEquals(PL_INFO, 0);
      delete [] gMeans;
      delete [] gVariances;
   }
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MomentAnalyzer& MomentAnalyzer::operator=(const MomentAnalyzer &)
{
   printOutTS(PL_ERROR, 
        "MomentAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int MomentAnalyzer::get_nSamples()
{
   return nSamples_;
}
int MomentAnalyzer::get_nGroups()
{
   return nGroups_;
}
int MomentAnalyzer::get_nInputs()
{
   return nInputs_;
}
int MomentAnalyzer::get_nOutputs()
{
   return nOutputs_;
}
double MomentAnalyzer::get_mean()
{
   if (moments_) return moments_[0];
   else          return 0.0;
}
double MomentAnalyzer::get_stdev()
{
   if (moments_) return moments_[1];
   else          return 0.0;
}
double MomentAnalyzer::get_skewness()
{
   if (moments_) return moments_[2];
   else          return 0.0;
}
double MomentAnalyzer::get_kurtosis()
{
   if (moments_) return moments_[3];
   else          return 0.0;
}

