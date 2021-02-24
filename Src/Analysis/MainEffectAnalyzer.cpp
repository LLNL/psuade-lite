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
// functions for the class MainEffectAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include "MainEffectAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "pData.h"
#include "PsuadeData.h"
#include "RSConstraints.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MainEffectAnalyzer::MainEffectAnalyzer() : Analyzer(), nInputs_(0), 
                        outputID_(0), inputVCE_(0), totalInputVCE_(0), 
                        mainEffectMean_(0), mainEffectStd_(0)
{
   setName("MainEffect");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MainEffectAnalyzer::~MainEffectAnalyzer()
{
   if (inputVCE_) delete [] inputVCE_;
}

// ************************************************************************
// perform VCE (McKay) analysis (nSubSamples : number of LHS points)
// ------------------------------------------------------------------------
double MainEffectAnalyzer::analyze(aData &adata)
{
  int    nInputs, nOutputs, nSamples, outputID, nSubSamples, status;
  int    ii, jj, ir, ss, nReplications, ncount, index;
  int    whichOutput, nSubs, iZero=0, printLevel;
  double *txArray, *tyArray, *vce, *meanVCEVar, *varVCEMean;
  double aMean, aVariance, *varVCEVar, totalVCE;
  double *sampleInputs, *sampleOutputs, *X, *Y, ddata;
  double *XX, *YY, **bsVCEs, *iLowerB, *iUpperB;
  char   winput1[500], winput2[500], meFileName[500];
  char   pString[500], **inputNames;
  FILE   *fp=NULL, *fp1=NULL;
  pData  pdata;
  PsuadeData    *ioPtr;
  RSConstraints *constrPtr=NULL;

  nInputs       = adata.nInputs_;
  nInputs_	 = nInputs;
  nOutputs      = adata.nOutputs_;
  nSamples      = adata.nSamples_;
  outputID      = adata.outputID_;
  outputID_	 = outputID;
  sampleInputs  = adata.sampleInputs_;
  sampleOutputs = adata.sampleOutputs_;
  nSubSamples   = adata.nSubSamples_;
  iLowerB       = adata.iLowerB_;
  iUpperB       = adata.iUpperB_;
  printLevel    = adata.printLevel_;
  whichOutput   = outputID;
  ioPtr         = adata.ioPtr_;

  if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
  {
    printOutTS(PL_ERROR,"MainEffect ERROR: invalid arguments.\n");
    printOutTS(PL_ERROR,"    nInputs  = %d\n", nInputs);
    printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR,"    nSamples = %d\n", nSamples);
    exit(1);
  } 
  if (nSamples/nSubSamples*nSubSamples != nSamples)
  {
    printOutTS(PL_ERROR,"MainEffect ERROR: nSamples != k*nLevels.\n");
    printOutTS(PL_ERROR,"    nSamples = %d\n", nSamples);
    printOutTS(PL_ERROR,"    nLevels  = %d\n", nSubSamples);
    exit(1);
  } 
  if (whichOutput >= nOutputs || whichOutput < 0)
  {
    printOutTS(PL_ERROR,"MainEffect ERROR: invalid outputID.\n");
    printOutTS(PL_ERROR,"    nOutputs = %d\n", nOutputs);
    printOutTS(PL_ERROR,"    outputID = %d\n", whichOutput+1);
    exit(1);
  }
  if (ioPtr == NULL)
  {
    printOutTS(PL_ERROR,"MainEffect ERROR: no data (PsuadeData).\n");
    return PSUADE_UNDEFINED;
  }
  if (adata.inputPDFs_ != NULL)
  {
    ncount = 0;
    for (ii = 0; ii < nInputs; ii++) ncount += adata.inputPDFs_[ii];
    if (ncount > 0)
    {
      printOutTS(PL_INFO, 
           "MainEffect INFO: Some inputs have non-uniform PDFs.\n");
      printOutTS(PL_INFO, 
           "           However, they are relevant in this analysis\n");
      printOutTS(PL_INFO, 
           "           (the sample should have been generated\n");
      printOutTS(PL_INFO,"            with the desired distributions.)\n");
    }
  }
  if (inputVCE_) delete [] inputVCE_;
  inputVCE_ = NULL;

  if (ioPtr != NULL)
  {
    constrPtr = new RSConstraints();
    if(constrPtr != NULL) constrPtr->genConstraints(ioPtr);
    else 
    {
      printOutTS(PL_ERROR,"out of memory in file %s line %d, exiting.\n",
                 __FILE__, __LINE__);
      exit(1);
    }
  }
   
  X = sampleInputs;
  Y = new double[nSamples];
  checkAllocate(Y, "Y in MainEffect::analyze");

  ncount = 0;
  for (ss = 0; ss < nSamples; ss++)
  {
    Y[ss] = sampleOutputs[nOutputs*ss+whichOutput];
    ddata = constrPtr->evaluate(&(X[ss*nInputs]), Y[ss], status);
    if (status == 0) Y[ss] = PSUADE_UNDEFINED; 
    else             ncount++;
  }
  if (ncount == 0)
  {
    printOutTS(PL_ERROR,"MainEffect ERROR: no valid sample point.\n");
    printOutTS(PL_ERROR,"    nSamples before filtering = %d\n", nSamples);
    printOutTS(PL_ERROR,"    nSamples after  filtering = %d\n", ncount);
    printOutTS(PL_ERROR, 
         "    INFO: check your data file for undefined's (1e35)\n");
    delete [] Y;
    if (constrPtr != NULL) delete constrPtr;
    return 1.0;
  } 
  if (ncount != nSamples)
  {
    printOutTS(PL_INFO, 
         "MainEffect: nSamples before filtering = %d\n", nSamples);
    printOutTS(PL_INFO, 
         "MainEffect: nSamples after  filtering  = %d\n", ncount);
  }
  if (nSamples < 1000)
  {
    printOutTS(PL_INFO,"MainEffect INFO: nSamples may be too small to\n");
    printOutTS(PL_INFO, 
         "                 give results with acceptable accuracy.\n");
  }

  computeMeanVariance(nInputs,1,nSamples,Y,&aMean,&aVariance,0);
  if (printLevel >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"*              Main Effect Analysis\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO, 
         "* Turn on higher printlevel to display more information\n");
    printOutTS(PL_INFO,"* Turn on ana_expert mode for more plots\n");
    printDashes(PL_INFO, 0);
    printOutTS(PL_INFO,"* Number of sample points = %10d\n",nSamples);
    printOutTS(PL_INFO,"* Number of Inputs        = %10d\n",nInputs);
    printEquals(PL_INFO, 0);
    printOutTS(PL_INFO,"Output %d\n", whichOutput+1);
    printOutTS(PL_INFO,"====> MainEffect: mean               = %12.4e\n",
               aMean);
    printOutTS(PL_INFO,"====> MainEffect: standard deviation = %12.4e\n",
                 sqrt(aVariance));
  }
  mainEffectMean_ = aMean;
  mainEffectStd_ = sqrt(aVariance);
  nSubs = nSubSamples;
  nReplications = nSamples / nSubs;
#if 0
  if (nReplications <= 1)
  {
    printOutTS(PL_INFO, 
         "MainEffectAnalyzer INFO: go no further since nReps = 1.\n");
    return 1.0;
  }
#endif
  vce           = new double[nInputs];
  meanVCEVar    = new double[nInputs];
  varVCEMean    = new double[nInputs];
  varVCEVar     = new double[nInputs];
  txArray       = new double[nSamples];
  tyArray       = new double[nSamples];
  checkAllocate(tyArray, "tyArray in MainEffect::analyze");

  for (ss = 0; ss < nSamples; ss++) tyArray[ss] = Y[ss];
  for (ii = 0; ii < nInputs; ii++)
  {
    for (ss = 0; ss < nSamples; ss++) txArray[ss] = X[nInputs*ss+ii];
    sortDbleList2(nSamples, txArray, tyArray);
    nReplications = 1;
    for (ss = 1; ss < nSamples; ss++)
    {
      if (txArray[ss] == txArray[0]) nReplications++;
      else                           break;
    }
    if (nReplications <= 5)
    {
      printOutTS(PL_INFO,
        "* MainEffect INFO: nReps <= 5 for input %d.\n",ii+1);
      printOutTS(PL_INFO,
        "*     ==> probably not replicated Latin hypercube\n");
      printOutTS(PL_INFO,"*     ==> crude main effect analysis.\n");
      computeVCECrude(nInputs, nSamples, X, Y, iLowerB, iUpperB, 
                      aVariance, vce);
      if (ioPtr != NULL)
      {
        pData *pPtr = ioPtr->getAuxData();
        pPtr->nDbles_ = nInputs;
        pPtr->dbleArray_ = new double[nInputs * nInputs];
        for (ii = 0; ii < nInputs; ii++)
          pPtr->dbleArray_[ii] = vce[ii];
        pPtr->dbleData_ = aVariance;
      }
      printResults(nInputs, aVariance, vce, ioPtr);

      delete [] vce;
      delete [] meanVCEVar;
      delete [] varVCEMean;
      delete [] varVCEVar;
      delete [] Y;
      delete [] txArray;
      delete [] tyArray;
      delete constrPtr;
      return 1.0;
    }
  }

  fp = NULL;
  if (psAnaExpertMode_ == 1)
  {
    sprintf(pString,"Create main effect scatter plot ? (y or n) ");
    getString(pString, winput1);
    if (winput1[0] == 'y')
    {
      if (psPlotTool_ == 1)
        sprintf(pString,"Enter scatter plot file name (ends with .sci): ");
      else
        sprintf(pString,"Enter scatter plot file name (ends with .m): ");
      getString(pString, meFileName);
      meFileName[strlen(meFileName)-1] = '\0';
      fp = fopen(meFileName, "w");
      if (fp != NULL)
        printOutTS(PL_INFO,
             "MainEffect: main effect file = %s\n",meFileName);
    }
  }

  if (fp != NULL)
  {
    strcpy(pString," **********************************************");
    fwriteComment(fp, pString);
    strcpy(pString," *  Main Effect Analysis                     **");
    fwriteComment(fp, pString);
    strcpy(pString," *-------------------------------------------**");
    fwriteComment(fp, pString);
    strcpy(pString," file for main effect plots \n");
    fwriteComment(fp, pString);
  }

  if (fp != NULL)
  {
    for (ii = 0; ii < nInputs; ii++)
    {
      for (ss = 0; ss < nSamples; ss++)
        txArray[ss] = X[nInputs*ss+ii];
      sortDbleList2(nSamples, txArray, tyArray);
      fprintf(fp, "X%d = [ \n", ii);
      fprintf(fp,"%24.16e\n",txArray[0]);
      for (ss = 1; ss < nSamples; ss++)
        if (txArray[ss] != txArray[ss-1])
          fprintf(fp,"%24.16e\n",txArray[ss]);
      fprintf(fp, "];\n");
    }
  }

  status = computeVCE(nInputs, nSamples, nSubSamples, X, Y, whichOutput, 
                       fp, varVCEMean, meanVCEVar, varVCEVar, vce);
  if (status == -1)
  {
    computeVCECrude(nInputs,nSamples,X,Y,iLowerB,iUpperB,aVariance,vce);
    printResults(nInputs, aVariance, vce, ioPtr);
    if (ioPtr != NULL)
    {
      pData *pPtr = ioPtr->getAuxData();
      pPtr->nDbles_ = nInputs;
      pPtr->dbleArray_ = new double[nInputs * nInputs];
      for (ii = 0; ii < nInputs; ii++) pPtr->dbleArray_[ii] = vce[ii];
      pPtr->dbleData_ = aVariance;
    }

    delete [] Y;
    delete [] vce;
    delete [] varVCEMean;
    delete [] meanVCEVar;
    delete [] varVCEVar;
    delete [] txArray;
    delete [] tyArray;
    if (constrPtr != NULL) delete constrPtr;
    if(fp != NULL) fclose(fp);
    return 0.0;
  }
  if (ioPtr != NULL)
  {
    pData *pPtr = ioPtr->getAuxData();
    pPtr->nDbles_ = nInputs;
    pPtr->dbleArray_ = new double[nInputs * nInputs];
    for (ii = 0; ii < nInputs; ii++) pPtr->dbleArray_[ii] = vce[ii];
    pPtr->dbleData_ = aVariance;
  }


#if 1
  if (printLevel >= 0)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"            McKay's correlation ratio\n");
    printEquals(PL_INFO,0);
    totalVCE = 0.0;
    if (aVariance == 0) printOutTS(PL_INFO, "Total VCE = %9.2e\n",totalVCE);
    else
    {
      for (ii = 0; ii < nInputs; ii++)
      {
        totalVCE += vce[ii] / aVariance;
        printOutTS(PL_INFO, 
           "Input %4d, normalized 1st-order effect = %9.2e (raw = %9.2e)\n",
           ii+1, vce[ii]/aVariance, vce[ii]);
      }
      printOutTS(PL_INFO, "Total VCE = %9.2e\n", totalVCE);
    }
  }
#endif

#if 1
  if (printLevel > 2)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,
         "     McKay's biased correlation ratio (stdVCEMean)\n");
    printDashes(PL_INFO, 0);
    totalVCE = 0.0;
    for (ii = 0; ii < nInputs; ii++)
    {
      printOutTS(PL_INFO, "Input %2d = %9.2e >>? %9.2e\n", ii+1,
             varVCEMean[ii]/aVariance,1.0/(double)(nSubs*nSubs));
      totalVCE += varVCEMean[ii] / aVariance;
    }
    printOutTS(PL_INFO, "Total VCE = %9.2e\n", totalVCE);
  }
#endif

  if (printLevel > 2)
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO, "           Strength of interaction (varVCEVar)\n");
    printDashes(PL_INFO, 0);
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO,"    Input %2d = %9.2e (meanVCEVar = %9.2e)\n",
                 ii+1, varVCEVar[ii], meanVCEVar[ii]);
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,"            Hora and Iman sensitivity index\n");
    printOutTS(PL_INFO,
         "   (may not be valid in the presence of constraints)\n");
    printDashes(PL_INFO, 0);
    totalVCE = 0.0;
    for (ii = 0; ii < nInputs; ii++) 
      totalVCE += sqrt(aVariance-meanVCEVar[ii]);
    for (ii = 0; ii < nInputs; ii++)
      printOutTS(PL_INFO, "    Input %2d = %9.2e\n", ii+1,
             sqrt(aVariance-meanVCEVar[ii])/totalVCE);
    printAsterisks(PL_INFO, 0);
  }

#ifdef HAVE_PYTHON
  PyObject *temp, *Xlist, *XIlist, *Ylist, *VCElist, *vVCEMlist;

  PyDict_SetItemString( AnalysisDataDict, "nInputs",
                 temp=PyInt_FromLong(nInputs) ); Py_DECREF(temp);
  PyDict_SetItemString( AnalysisDataDict, "outputID",
                 temp=PyInt_FromLong(whichOutput) ); Py_DECREF(temp);
  PyDict_SetItemString( AnalysisDataDict, "nSamples",
                 temp=PyInt_FromLong(nSamples) ); Py_DECREF(temp);
  PyDict_SetItemString( AnalysisDataDict, "nSubs",
                 temp=PyInt_FromLong(nSubs) ); Py_DECREF(temp);
  PyDict_SetItemString( AnalysisDataDict, "aVariance",
                 temp=PyFloat_FromDouble(aVariance) ); Py_DECREF(temp);

  Xlist     = PyList_New(0);
  VCElist   = PyList_New(0);
  vVCEMlist = PyList_New(0);
  for (ii = 0; ii < nInputs; ii++)
  {
    PyList_Append( VCElist, temp=PyFloat_FromDouble(vce[ii]) );
    Py_DECREF(temp);
    PyList_Append( vVCEMlist, temp=PyFloat_FromDouble(varVCEMean[ii]) );
    Py_DECREF(temp);

    XIlist = PyList_New(0);
    for (ss = 0; ss < nSamples; ss++)
    {
      PyList_Append( XIlist, temp=PyFloat_FromDouble(X[nInputs*ss+ii]) );
      Py_DECREF(temp);
    }
    PyList_Append( Xlist, XIlist ); Py_DECREF(XIlist);
  }
  PyDict_SetItemString( AnalysisDataDict, "X", Xlist ); Py_DECREF(Xlist);
  PyDict_SetItemString( AnalysisDataDict, "vce", VCElist ); Py_DECREF(VCElist);
  PyDict_SetItemString( AnalysisDataDict, "varVCEMean", vVCEMlist );
  Py_DECREF(vVCEMlist);

  Ylist = PyList_New(0);
  for (ss = 0; ss < nSamples; ss++)
  {
    PyList_Append( Ylist, temp=PyFloat_FromDouble(Y[ss]) ); Py_DECREF(temp);
  }
  PyDict_SetItemString( AnalysisDataDict, "Y",	Ylist ); Py_DECREF(Ylist);
#endif

  if (fp != NULL)
  {
    fwriteHold(fp, 0);
    fprintf(fp,"nx = ceil(sqrt(%d));", nInputs);
    for (ii = 0; ii < nInputs; ii++)
    {
      fprintf(fp,"subplot(nx, nx, %d)\n", ii+1);
      if (psPlotTool_ == 1)
        fprintf(fp,"Y = matrix(Y%d_%d,nReps%d,nSymbols%d);\n",
                whichOutput,ii,ii,ii);
      else
        fprintf(fp,"Y = reshape(Y%d_%d,nReps%d,nSymbols%d);\n",
                whichOutput,ii,ii,ii);
      fprintf(fp,"ymax = max(Y%d_%d);\n", whichOutput, ii);
      fprintf(fp,"ymin = min(Y%d_%d);\n", whichOutput, ii);
      fprintf(fp,"xmax = max(X%d);\n", ii);
      fprintf(fp,"xmin = min(X%d);\n", ii);
      if (psPlotTool_ == 1)
      {
        fprintf(fp,"a = get(\"current_axes\");\n");
        fprintf(fp,"a.data_bounds=[xmin,ymin;xmax,ymax];\n");
      }
      else
      {
        fprintf(fp,"axis([xmin xmax ymin ymax])\n");
      }
      fprintf(fp,"for k = 1 : nSymbols%d \n",ii);
      fprintf(fp,"   yt = Y(:, k);\n");
      fprintf(fp,"   mean = sum(yt) / nReps%d;\n",ii);
      fprintf(fp,"   x = X%d(k) * ones(nReps%d, 1);\n",ii,ii);
      fprintf(fp,"   plot(x, yt, '.', 'MarkerSize', 1)\n");
      fprintf(fp,"   if ( k == 1 )\n");
      fwriteHold(fp, 1);
      fprintf(fp,"   end;\n");
      fprintf(fp,"   meanArray(k) = mean;\n");
      fprintf(fp,"   tmaxArray(k) = max(yt);\n");
      fprintf(fp,"   tminArray(k) = min(yt);\n");
      fprintf(fp,"end;\n");
      fwritePlotAxes(fp);
      fwritePlotXLabel(fp, "Input");
      fwritePlotYLabel(fp, "Output");
      fprintf(fp, "meanMean = sum(meanArray)/nSymbols%d;\n",ii);
      fprintf(fp, "adjMean  = (meanArray - meanMean);\n");
      fprintf(fp,"varMean = sum(adjMean.^2)/nSymbols%d\n",ii);
      fprintf(fp,"[xx,I] = sort(X%d);\n", ii);
      fprintf(fp,"yy1 = meanArray(I);\n");
      fprintf(fp,"yy2 = tmaxArray(I);\n");
      fprintf(fp,"yy3 = tminArray(I);\n");
      fprintf(fp,"plot(xx, yy1, 'r-', 'lineWidth', 1)\n");
      fprintf(fp,"plot(xx, yy2, 'r-', 'lineWidth', 1)\n");
      fprintf(fp,"plot(xx, yy3, 'r-', 'lineWidth', 1)\n");
      fprintf(fp,"title('Output %d vs Input %d')\n",whichOutput+1,ii+1);
      fprintf(fp,"disp('The center red lines are the conditional means.')\n");
      fprintf(fp,"disp('Conditional means show trends wrt each input.')\n");
      fprintf(fp,"disp('The upper/lower red lines are upper/lower envelopes.')\n");
      fwriteHold(fp, 0);
    }
    fclose(fp);
    printOutTS(PL_INFO, 
         "Main effect matlab plot %s has been generated.\n", meFileName);
  }

  if (psAnaExpertMode_ == 1)
  {
    printOutTS(PL_INFO, 
         "Bootstrap analysis takes the sample, replicates it n times,\n");
    printOutTS(PL_INFO, 
         "and assess whether the sensitivity indices have converged.\n");
    printOutTS(PL_INFO, 
         "If you are performed iterative analysis with refinements,\n");
    printOutTS(PL_INFO, 
         "you will need to enter 'no index reuse' below at the first\n");
    printOutTS(PL_INFO, 
         "iteration and 'yes' afterward until the final refinement.\n");
    sprintf(pString,"Perform bootstrap main effect analysis? (y or n) ");
    getString(pString, winput1);
    ncount = 0;
    if (winput1[0] == 'y')
    {
      sprintf(pString,"Number of bootstrap samples to use (>=100): ");
      ncount = getInt(100, 2000, pString);
      bsVCEs = new double*[ncount];
      for (ii = 0; ii < ncount; ii++) bsVCEs[ii] = new double[nInputs];
      XX = new double[nInputs*nSamples];
      YY = new double[nSamples];
      checkAllocate(YY, "YY in MainEffect::analyze");
      nReplications = nSamples / nSubSamples;
      winput1[0] = 'n';
      fp1 = fopen(".ME_bootstrap_indset", "r");
      if (fp1 != NULL)
      {
        printOutTS(PL_INFO, ".ME_bootstrap_indset file found.\n");
        sprintf(pString,"Re-use file? (y or n) ");
        getString(pString, winput1);
        if (winput1[0] == 'y')
        {
          fscanf(fp1, "%d", &ii);
          if (ii != nReplications*ncount)
          {
            printOutTS(PL_ERROR,"ERROR: expect the first line to be %d.\n",
                   nReplications*ncount);
            printOutTS(PL_ERROR, 
                 "       Instead found the first line to be %d.n",ii);
            exit(1);
          }
        }
        else
        {
          fclose(fp1);
          fp1 = fopen(".ME_bootstrap_indset", "w");
          if (fp1 == NULL)
          {
            printOutTS(PL_ERROR, 
                 "ERROR: cannot open ME_bootstrap_indset file.\n");
            exit(1);
          }
          fprintf(fp1, "%d\n", nReplications*ncount);
        }
      }
      else
      {
        fp1 = fopen(".ME_bootstrap_indset", "w");
        if (fp1 == NULL)
        {
          printOutTS(PL_ERROR, 
               "ERROR: cannot open .ME_bootstrap_indset file.\n");
          exit(1);
        }
        fprintf(fp1, "%d\n", nReplications*ncount);
      }
      fp = NULL;
      for (ii = 0; ii < ncount; ii++)
      {
        for (ir = 0; ir < nReplications; ir++)
        {
          if (fp1 != NULL && winput1[0] == 'y')
          {
            fscanf(fp1, "%d\n", &index);
            if (index < 0 || index >= nReplications)
            {
              printOutTS(PL_ERROR, 
                   "ERROR: reading index from file .ME_bootstrap_indset\n");
              printOutTS(PL_ERROR,"       index read = %d\n", index);
              printOutTS(PL_ERROR,"       expected   = [0,%d]\n", 
                         nReplications-1);
            }
          }
          else index = PSUADE_rand() % nReplications;

          if (fp1 != NULL && winput1[0] != 'y')
            fprintf(fp1, "%d\n", index);
	  // Bill Oliver range check
	  if((index*nSubSamples + nSubSamples - 1) >= nSamples)
          {
	    printOutTS(PL_ERROR, "Buffer overflow in file %s line %d\n",
                      __FILE__, __LINE__);
	    exit(1);
	  }
          for (ss = 0; ss < nSubSamples; ss++)
          { 
            for (jj = 0; jj < nInputs; jj++)
              XX[(ir*nSubSamples+ss)*nInputs+jj] =
                      X[(index*nSubSamples+ss)*nInputs+jj];
            YY[ir*nSubSamples+ss] = Y[index*nSubSamples+ss];
          }
        }
        computeMeanVariance(nInputs,1,nSamples,YY,&aMean,&aVariance,0);
        computeVCE(nInputs, nSamples, nSubSamples, XX, YY, 
                   iZero, fp, varVCEMean, meanVCEVar, 
                   varVCEVar, vce);
        for (jj = 0; jj < nInputs; jj++) 
          bsVCEs[ii][jj] = vce[jj]/aVariance;
      }
      if (fp1 != NULL) fclose(fp1);
      printf("Enter name of matlab file to store bootstrap info: ");
      scanf("%s", winput1);
      fgets(winput2,500,stdin);
      fp = fopen(winput1, "w");
      if (fp != NULL)
      {
        strcpy(winput2,"bootstrap sample of main effects");
        fwriteComment(fp, winput2);
        strcpy(winput2,"VCE(i,j) = VCE for input i, bootstrapped sample j");
        fwriteComment(fp, winput2);
        fprintf(fp, "VCE = zeros(%d,%d);\n",nInputs,ncount);
        for (ii = 0; ii < ncount; ii++)
           for (jj = 0; jj < nInputs; jj++)
              fprintf(fp, "VCE(%d,%d) = %e;\n",jj+1,ii+1,
                      bsVCEs[ii][jj]);
        fprintf(fp, "nn = %d;\n", nInputs);
        if (ioPtr != NULL) ioPtr->getParameter("input_names",pdata);
        if (pdata.strArray_ != NULL) inputNames = pdata.strArray_;
        else                         inputNames = NULL;
        if (inputNames == NULL)
        {
          fprintf(fp, "Str = {");
          for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
          fprintf(fp,"'X%d'};\n",nInputs);
        }
        else
        {
          fprintf(fp, "Str = {");
          for (ii = 0; ii < nInputs-1; ii++)
          {
            if (inputNames[ii] != NULL) fprintf(fp,"'%s',",inputNames[ii]);
            else                        fprintf(fp,"'X%d',",ii+1);
          }
          if (inputNames[nInputs-1] != NULL)
               fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
          else fprintf(fp,"'X%d'};\n",nInputs);
        }
        fwriteHold(fp,0);
        fprintf(fp,"ymin = 0.0;\n");
        fprintf(fp,"ymax = max(max(VCE));\n");
        fprintf(fp,"hh = 0.05 * (ymax - ymin);\n");
        fprintf(fp,"ymax = ymax + hh;\n");
        fprintf(fp,"VMM = sum(VCE')/%d;\n",ncount);
        fprintf(fp,"VMA = max(VCE');\n");
        fprintf(fp,"VMI = min(VCE');\n");
        fprintf(fp,"bar(VMM,0.8);\n");
        fwriteHold(fp,1);
        fprintf(fp,"for ii = 1 : nn\n");
        fprintf(fp,"  XX = [ii ii];\n");
        fprintf(fp,"  YY = [VMI(ii)  VMA(ii)];\n");
        fprintf(fp,"  plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
        fprintf(fp,"'k','MarkerFaceColor','g','MarkerSize',13)\n");
        fprintf(fp,"end;\n");
        fwritePlotAxes(fp);
        if (psPlotTool_ == 1)
        {
          fprintf(fp,"a=gca();\n");
          fprintf(fp,"a.data_bounds=[0, ymin; nn+1, ymax];\n");
          fprintf(fp,"newtick = a.x_ticks;\n");
          fprintf(fp,"newtick(2) = [1:nn]';\n");
          fprintf(fp,"newtick(3) = Str';\n");
          fprintf(fp,"a.x_ticks = newtick;\n");
          fprintf(fp,"a.x_label.font_size = 3;\n");
          fprintf(fp,"a.x_label.font_style = 4;\n");
        }
        else
        {
          fprintf(fp,"axis([0  nn+1 ymin ymax])\n");
          fprintf(fp,"set(gca,'XTickLabel',[]);\n");
          fprintf(fp,"th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),");
          fprintf(fp,"Str,'HorizontalAlignment','left','rotation',90);\n");
          fprintf(fp,"set(th, 'fontsize', 12)\n");
          fprintf(fp,"set(th, 'fontweight', 'bold')\n");
        }
        fwritePlotTitle(fp,"Bootstrapped Sobol first  Order Indices");
        fwritePlotYLabel(fp,"Sobol Indices");
        fclose(fp);
      }
      else
      {
        printOutTS(PL_ERROR, "ERROR: cannot open file %s\n", winput1);
      }
      for (ii = 0; ii < ncount; ii++) delete [] bsVCEs[ii];
      delete [] bsVCEs;
      delete [] XX;
      delete [] YY;
    }
  }

  printResults(nInputs, aVariance, vce, ioPtr);

  delete [] Y;
  delete [] vce;
  delete [] varVCEMean;
  delete [] meanVCEVar;
  delete [] varVCEVar;
  delete [] txArray;
  delete [] tyArray;
  if (constrPtr != NULL) delete constrPtr;

  // return 1.0 to facilitate continuous refinement
  return 1.0;
}

// *************************************************************************
// Compute agglomerated mean and variance
// -------------------------------------------------------------------------
int MainEffectAnalyzer::computeMeanVariance(int nInputs, int nOutputs, 
                        int nSamples, double *Y, double *aMean, 
                        double *aVariance, int outputID)
{
   int    ss, count;
   double mean, variance;

   count = 0;
   mean = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      if (Y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         mean += Y[nOutputs*ss+outputID];
         count++;
      }
   }
   if (count <= 0)
   {
      printOutTS(PL_ERROR, "MainEffect ERROR: no valid data.\n");
      exit(1);
   }
   mean /= (double) count;
   variance = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      if (Y[nOutputs*ss+outputID] < 0.9*PSUADE_UNDEFINED)
      {
         variance += ((Y[nOutputs*ss+outputID] - mean) *
                      (Y[nOutputs*ss+outputID] - mean));
      }
   }
   variance /= (double) count;
   (*aMean) = mean;
   (*aVariance) = variance;
   return 0;
}

// *************************************************************************
// Compute VCE 
// -------------------------------------------------------------------------
int MainEffectAnalyzer::computeVCE(int nInputs, int nSamples, int nSubSamples, 
                                   double *X, double *Y, int whichOutput, 
                                   FILE *fp, double *varVCEMean,
                                   double *meanVCEVar, double *varVCEVar,
                                   double *vce)
{
   int    ii, ss, nReplications, nSubs, ncount, symIndex, repID, subID;
   int    validnReps, validnSubs, totalCnt, *bins;
   double *txArray, *tyArray, *vceMean, *vceVariance, fixedInput, aMean;

   bins    = new int[nSubSamples];
   txArray = new double[nSamples];
   tyArray = new double[nSamples];
   vceMean     = new double[nSubSamples];
   vceVariance = new double[nSubSamples];
   checkAllocate(vceVariance, "vceVariance in MainEffect::computeVCE");

   for (ii = 0; ii < nInputs; ii++)
   {
      if (fp != NULL) fprintf(fp, "Y%d_%d = [\n",whichOutput,ii);

      for (ss = 0; ss < nSamples; ss++)
      {
         txArray[ss] = X[nInputs*ss+ii];
         tyArray[ss] = Y[ss];
      }
      sortDbleList2(nSamples, txArray, tyArray);
      nReplications = 1;
      for (ss = 1; ss < nSamples; ss++)
      {
         if (txArray[ss] == txArray[0]) nReplications++;
         else                           break;
      }
      nSubs = nSamples / nReplications;

      ncount = 0;
      for (ss = 0; ss < nSamples; ss+=nReplications)
      {
         symIndex = ss / nReplications;
         fixedInput = txArray[ss];
         vceMean[symIndex] = 0.0;
         validnReps = 0;
         for (repID = 0; repID < nReplications; repID++)
         {
            if (PABS(txArray[ss+repID] - fixedInput) < 1.0E-12)
            {
               ncount++;
               if (tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
               {
                  validnReps++;
                  vceMean[symIndex] += tyArray[ss+repID];
                  if (fp != NULL)
                     fprintf(fp,"%24.16e\n",tyArray[ss+repID]);
               }
            }
         }
         bins[symIndex] = validnReps;
         if (validnReps > 0) vceMean[symIndex] /= (double) validnReps;
         else                vceMean[symIndex] = PSUADE_UNDEFINED;
         vceVariance[symIndex] = 0.0;
         for (repID = 0; repID < nReplications; repID++)
         {
            if (vceMean[symIndex] != PSUADE_UNDEFINED &&
                tyArray[ss+repID] < 0.9*PSUADE_UNDEFINED)
            {
               vceVariance[symIndex] += 
                    ((tyArray[ss+repID]-vceMean[symIndex])*
                     (tyArray[ss+repID]-vceMean[symIndex]));
            }
         }
         if (validnReps > 0)
              vceVariance[symIndex] /= (double) validnReps;
         else vceVariance[symIndex] = PSUADE_UNDEFINED;
      }
      if (fp != NULL)
      { 
         fprintf(fp,"];\n");
         fprintf(fp, "nSymbols%d = %d;\n", ii, nSubs);
         fprintf(fp, "nReps%d    = %d;\n", ii, nReplications);
      }

      if (ncount != nSamples)
      {
         printf("MainEffect ERROR: not rLHS, input %d\n", ii+1);
         printf("       error data = %d (%d)\n",ncount, nSamples);
         printf("       Did you use rLHS ?\n");
         delete [] bins;
         delete [] txArray;
         delete [] tyArray;
         delete [] vceMean;
         delete [] vceVariance;
         return -1;
      }

      totalCnt = 0;
      for (subID = 0; subID < nSubs; subID++) totalCnt += bins[subID];
      varVCEMean[ii] = 0.0;
      meanVCEVar[ii] = 0.0;
      aMean = 0.0;
      validnSubs = 0;
      for (subID = 0; subID < nSubs; subID++)
      {
         if (vceVariance[subID] != PSUADE_UNDEFINED)
         {
            aMean += (vceMean[subID] / totalCnt * bins[subID]);
            validnSubs++;
         }
      }
      for (subID = 0; subID < nSubs; subID++)
      {
         if (vceVariance[subID] != PSUADE_UNDEFINED)
         {
            varVCEMean[ii] += ((vceMean[subID] - aMean) *
                               (vceMean[subID] - aMean) * bins[subID]/totalCnt);
            meanVCEVar[ii] += vceVariance[subID] * bins[subID] / totalCnt;
         }
      }
      varVCEVar[ii] = 0.0;
      for (subID = 0; subID < nSubs; subID++)
      {
         if (vceVariance[subID] != PSUADE_UNDEFINED)
            varVCEVar[ii] += ((vceVariance[subID] - meanVCEVar[ii]) *
                              (vceVariance[subID] - meanVCEVar[ii]) * 
                              bins[subID] / totalCnt);
      }
      if (validnSubs > 0)
         vce[ii] = varVCEMean[ii];
      else
         vce[ii] = 0.0;
   }
   delete [] txArray;
   delete [] tyArray;
   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   return 0;
}

// *************************************************************************
// Compute VCE 
// -------------------------------------------------------------------------
int MainEffectAnalyzer::computeVCECrude(int nInputs, int nSamples, 
                          double *X, double *Y, double *iLowerB, 
                          double *iUpperB, double aVariance, double *vce)
{
   int    ii, ss, nSize, index;
   int    totalCnt, *bins, nIntervals, *tags;
   double *vceMean, *vceVariance, aMean, ddata, hstep;
   char   pString[500];

   nInputs_ = nInputs;
   if (nSize < 10) nSize = 10;
   if (isScreenDumpModeOn())
   {
     printAsterisks(PL_INFO, 0);
     printOutTS(PL_INFO,"*                Crude Main Effect\n");
     printEquals(PL_INFO, 0);
     printOutTS(PL_INFO, 
       "* For small to moderate sample sizes, this method gives rough\n");
     printOutTS(PL_INFO, 
       "* estimates of main effect (first order sensitivity). These\n");
     printOutTS(PL_INFO, 
       "* estimates can vary with different choices of internal settings.\n");
     printOutTS(PL_INFO, 
       "* For example, try different number of levels to assess sensitivity\n");
     printOutTS(PL_INFO, 
       "* of the computed measures with respect to it.\n");
     printOutTS(PL_INFO, 
       "* Turn on analysis expert mode to change the settings.\n");
   }
   nIntervals = (int) sqrt(1.0 * nSamples);
   nSize = nSamples / nIntervals;
   if (isScreenDumpModeOn())
   {
     printOutTS(PL_INFO,"* MainEffect: number of levels   = %d\n", nIntervals);
     printOutTS(PL_INFO,"* MainEffect: sample size/levels = %d\n", nSize);
   }
   if (psAnaExpertMode_ == 1 && isInteractiveModeOn())
   {
      sprintf(pString,"number of levels (>5, default = %d): ", nIntervals);
      nIntervals = getInt(5, nSamples/2, pString);
      nSize = nSamples / nIntervals;
   }
   if (isScreenDumpModeOn()) printEquals(PL_INFO, 0);
   bins = new int[nIntervals];
   tags = new int[nSamples];
   vceMean     = new double[nIntervals];
   vceVariance = new double[nIntervals];

   inputVCE_ = new double[nInputs_];
   checkAllocate(inputVCE_, "inputVCE in MainEffect::computeVCECrude");
   for (ii = 0; ii < nInputs; ii++)
   {
      hstep = (iUpperB[ii] - iLowerB[ii]) / nIntervals;
      for (ss = 0; ss < nIntervals; ss++)
      {
         vceMean[ss] = 0.0;
         vceVariance[ss] = 0.0;
         bins[ss] = 0;
      }
      for (ss = 0; ss < nSamples; ss++)
      {
         ddata = (X[nInputs*ss+ii] - iLowerB[ii]) / hstep;
         index = (int) ddata;
         if (index <  0) index = 0;
         if (index >= nIntervals) index = nIntervals - 1;
         tags[ss] = -1;
         if (Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceMean[index] += Y[ss];
            bins[index]++;
            tags[ss] = index;
         }
      }
      for (ss = 0; ss < nIntervals; ss++)
         if (bins[ss] > 0) vceMean[ss] /= (double) bins[ss];
      for (ss = 0; ss < nSamples; ss++)
      {
         index = tags[ss];
         if (index >= 0 && Y[ss] < 0.9*PSUADE_UNDEFINED)
         {
            vceVariance[index] += ((Y[ss]-vceMean[index])*
                                   (Y[ss]-vceMean[index]));
         }
      }
      for (ss = 0; ss < nIntervals; ss++)
      {
         if (bins[ss] > 0)
              vceVariance[ss] /= (double) bins[ss];
         else vceVariance[ss] = PSUADE_UNDEFINED;
      }
         
      totalCnt = 0;
      for (ss = 0; ss < nIntervals; ss++) totalCnt += bins[ss];
      vce[ii] = 0.0;
      aMean = 0.0;
      for (ss = 0; ss < nIntervals; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
            aMean += (vceMean[ss] / totalCnt * bins[ss]);
      }
      for (ss = 0; ss < nIntervals; ss++)
      {
         if (vceVariance[ss] != PSUADE_UNDEFINED)
         {
            vce[ii] += ((vceMean[ss] - aMean) *
                        (vceMean[ss] - aMean) * bins[ss]/totalCnt);
         }
      }
      inputVCE_[ii] = vce[ii];
   }

   ddata = 0.0;
   if (isScreenDumpModeOn()) 
   {
     if (aVariance == 0) printOutTS(PL_INFO, "Total VCE = %9.2e\n", ddata);
     else
     {
       for (ii = 0; ii < nInputs; ii++)
       {
         ddata += vce[ii] / aVariance;
         printOutTS(PL_INFO, 
              "Input %4d, normalized 1st-order effect = %9.2e (raw = %9.2e)\n",
              ii+1, vce[ii]/aVariance, vce[ii]);
       }
       printOutTS(PL_INFO, "Total VCE = %9.2e\n", ddata);
     }
     printAsterisks(PL_INFO, 0);
   }
   else
   {
     if (aVariance != 0) 
       for (ii = 0; ii < nInputs; ii++) ddata += vce[ii] / aVariance;
   }
   totalInputVCE_ = ddata;


   delete [] vceMean;
   delete [] vceVariance;
   delete [] bins;
   delete [] tags;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MainEffectAnalyzer& MainEffectAnalyzer::operator=(const MainEffectAnalyzer &)
{
   printOutTS(PL_ERROR,"MainEffect operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// print result
// ------------------------------------------------------------------------
int MainEffectAnalyzer::printResults(int nInputs, double variance,
                                     double *mEffect, PsuadeData *ioPtr)
{
   int   ii;
   FILE  *fp;
   char  **iNames, pString[500];
   pData qData;

   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);
   if (qData.strArray_ != NULL) iNames = qData.strArray_;
   else                         iNames = NULL;
   printEquals(PL_INFO, 0);
   if (variance == 0.0)
   {
      printOutTS(PL_WARN, "Total variance = 0. Hence, no main effect plot.\n");
      return 0;
   }
   if (psPlotTool_ == 1) fp = fopen("scilabme.sci", "w");
   else                  fp = fopen("matlabme.m", "w");
   if (fp != NULL)
   {
      strcpy(pString," This file contains Sobol' first order indices");
      fwriteComment(fp, pString);
      strcpy(pString," set sortFlag = 1 and set nn to be the number");
      fwriteComment(fp, pString);
      strcpy(pString," of inputs to display.");
      fwriteComment(fp, pString);
      fprintf(fp, "sortFlag = 0;\n");
      fprintf(fp, "nn = %d;\n", nInputs);
      fprintf(fp, "Mids = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp,"%24.16e\n", mEffect[ii]/variance);
      fprintf(fp, "];\n");
      if (iNames == NULL)
      {
         fprintf(fp, "Str = {");
         for (ii = 0; ii < nInputs-1; ii++) fprintf(fp,"'X%d',",ii+1);
         fprintf(fp,"'X%d'};\n",nInputs);
      }
      else
      {
         fprintf(fp, "Str = {");
         for (ii = 0; ii < nInputs-1; ii++)
         {
            if (iNames[ii] != NULL) fprintf(fp,"'%s',",iNames[ii]);
            else                    fprintf(fp,"'X%d',",ii+1);
         }
         if (iNames[nInputs-1] != NULL)
              fprintf(fp,"'%s'};\n",iNames[nInputs-1]);
         else fprintf(fp,"'X%d'};\n",nInputs);
      }
      fwritePlotCLF(fp);
      fprintf(fp, "if (sortFlag == 1)\n");
      if (psPlotTool_ == 1) 
           fprintf(fp, "  [Mids, I2] = gsort(Mids);\n");
      else fprintf(fp, "  [Mids, I2] = sort(Mids,'descend');\n");
      fprintf(fp, "  Str  = Str(I2);\n");
      fprintf(fp, "  I2 = I2(1:nn);\n");
      fprintf(fp, "  Mids = Mids(1:nn);\n");
      fprintf(fp, "  Str  = Str(1:nn);\n");
      fprintf(fp, "end\n");
      fprintf(fp, "ymin = min(Mids);\n");
      fprintf(fp, "ymin = 0.0;\n");
      fprintf(fp, "ymax = max(Mids);\n");
      fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
      fprintf(fp, "bar(Mids,0.8);\n");
      fwritePlotAxes(fp);
      if (psPlotTool_ == 1)
      {
         fprintf(fp,"a=gca();\n");
         fprintf(fp,"a.data_bounds=[0, ymin; nn+1, ymax];\n");
         fprintf(fp,"newtick = a.x_ticks;\n");
         fprintf(fp,"newtick(2) = [1:nn]';\n");
         fprintf(fp,"newtick(3) = Str';\n");
         fprintf(fp,"a.x_ticks = newtick;\n");
         fprintf(fp,"a.x_label.font_size = 3;\n");
         fprintf(fp,"a.x_label.font_style = 4;\n");
      }
      else
      {
         fprintf(fp,"axis([0  nn+1 ymin ymax])\n");
         fprintf(fp,"set(gca,'XTickLabel',[]);\n");
         fprintf(fp,"th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str,");
         fprintf(fp,"'HorizontalAlignment','left','rotation',90);\n");
         fprintf(fp,"set(th, 'fontsize', 12)\n");
         fprintf(fp,"set(th, 'fontweight', 'bold')\n");
      }
      fwritePlotTitle(fp,"Sobol First Order Indices");
      fwritePlotYLabel(fp,"Sobol Indices");
      fprintf(fp,"disp('Switch sortFlag to display ranked Sobol indices')\n");
      fwriteHold(fp, 0);
      fclose(fp);
      if (psPlotTool_ == 1) 
           printOutTS(PL_INFO, "MainEffect plot matlab file = scilabme.sci\n");
      else printOutTS(PL_INFO, "MainEffect plot matlab file = matlabme.m\n");
      return 0;
   }
   else
   {
      printOutTS(PL_ERROR,"MainEffect ERROR: cannot create matlabme.m file.\n");
      return 0;
   }
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int MainEffectAnalyzer::get_nInputs()
{
   return nInputs_;
}
int MainEffectAnalyzer::get_outputID()
{
   return outputID_;
}
double *MainEffectAnalyzer::get_inputVCE()
{
   double* retVal = NULL;
   if (inputVCE_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in MainEffect::get_inputVCE");
      std::copy(inputVCE_, inputVCE_+nInputs_, retVal);
   }
   return retVal;
}
double MainEffectAnalyzer::get_totalInputVCE()
{
   return totalInputVCE_;
}
double MainEffectAnalyzer::get_mainEffectMean()
{
   return mainEffectMean_;
}
double MainEffectAnalyzer::get_mainEffectStd()
{
   return mainEffectStd_;
}

