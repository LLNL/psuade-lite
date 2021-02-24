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
// Functions for the class RSFuncApproxAnalyzer  
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "FuncApprox.h"
#include "Ann.h"
#include "Regression.h"
#include "UserRegression.h"
#include "RSFuncApproxAnalyzer.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
RSFuncApproxAnalyzer::RSFuncApproxAnalyzer() : Analyzer()
{
  setName("RSFA");
  rsType_ = PSUADE_RS_MARS;
  useCV_ = 0;
  numCVGroups_ = 10;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
RSFuncApproxAnalyzer::~RSFuncApproxAnalyzer() 
{ 
} 

// ************************************************************************ 
// perform analysis 
// ------------------------------------------------------------------------
double RSFuncApproxAnalyzer::analyze(aData &adata)
{
  int        nInputs, nOutputs, nSamples, outputID, nLevels, *levelSeps;
  int        ss, ss2, iI, iL, wgtID, nLast, nPtsPerDim=64, status, rsState;
  int        count, nSubSamples, printLevel, testFlag, iOne=1, *arrayS;
  double     ddata, ymax, ymin, retdata=0, sdata, ymean, yvar, ydiff, ssum;
  double     cvErr1, cvErr1s, cvErr2, cvErr2s,cvMax, cvMaxs, cvMaxBase;
  double     cvMaxBases, CVMaxBases, maxBase, maxBases, sumErr11, sumErr11s;
  double     CVErr2, CVErr2s, CVErr1, CVErr1s, CVMax, CVMaxs, CVMaxBase;
  double     sumErr1, sumErr1s, sumErr2, sumErr2s, maxErr, maxErrs;
  double     *lower, *upper, *arrayX, *arrayY, *arrayXX;
  char       winput1[500], winput2[500], dataFile[500], errFile[500];
  char       pString[500], *targv[2], *cString;
  FILE       *fpData, *fpErr;
  psVector   VecYLocal,VecYT,VecWgts,VecX2,VecY2,VecSigmas,VecS,VecE,VecW;
  psVector   VecS2, VecXX, VecYY, VecD, VecDD, VecD2;
  psIVector  IVec1, IVec2;
  FuncApprox *faPtr;

  printLevel = adata.printLevel_;
  nInputs   = adata.nInputs_;
  nOutputs  = adata.nOutputs_;
  nSamples  = adata.nSamples_;
  lower     = adata.iLowerB_;
  upper     = adata.iUpperB_;
  arrayX    = adata.sampleInputs_;
  arrayY    = adata.sampleOutputs_;
  arrayS    = adata.sampleStates_;
  outputID  = adata.outputID_;
  nLevels   = adata.currRefineLevel_;
  levelSeps = adata.refineSeparators_;
  wgtID     = adata.regWgtID_;
  if (adata.inputPDFs_ != NULL)
  {
    count = 0;
    for (iI = 0; iI < nInputs; iI++) count += adata.inputPDFs_[iI];
    if (count > 0)
    {
      printOutTS(PL_INFO,
           "RSAnalysis INFO: some inputs have non-uniform PDFs,\n");
      printOutTS(PL_INFO,
           "    but they are not relevant in this analysis.\n");
    }
  }

  if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
  {
     printOutTS(PL_ERROR, "RSAnalyzer ERROR: invalid arguments.\n");
     printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
     printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
     printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
     return PSUADE_UNDEFINED;
  } 
  status = 0;
  for (ss = 0; ss < nSamples; ss++)
     if (arrayY[nOutputs*ss+outputID] > 0.9*PSUADE_UNDEFINED) status = 1;
  if (status == 1)
  {
    printOutTS(PL_ERROR,"RSAnalysis ERROR: Some outputs are undefined.\n");
    printOutTS(PL_ERROR,"    Prune the undefined sample points first.\n");
    return PSUADE_UNDEFINED;
  }
  if (arrayS != NULL)
  {
    status = 0;
    for (ss = 0; ss < nSamples; ss++)
      if (arrayS[ss] != 1) status = 1;
    if (status == 1)
    {
      printOutTS(PL_ERROR,
         "RSAnalysis WARNING: Some sample statuses are not valid (=1).\n");
      return PSUADE_UNDEFINED;
    }
  }
  if (printLevel > 0)
  {
    printAsterisks(PL_INFO, 0);
    printThisFA(rsType_);
    printEquals(PL_INFO, 0);
  }
  if (rsType_ == PSUADE_RS_REGRGL) nLevels = 1;
  if (rsType_ == PSUADE_RS_REGRGL && nOutputs != nInputs+1) 
  {
    printOutTS(PL_ERROR,
       "RSAnalysis ERROR: for gradient-based RS, nOutputs=nInputs+1\n");
    printOutTS(PL_ERROR,
       "  Your loaded sample has nInputs/nOutputs = %d/%d\n",nInputs,nOutputs);
    exit(1);
  }
  if (rsType_ == PSUADE_RS_REGRGL && outputID != 0) 
  {
    printOutTS(PL_ERROR,
       "RSAnalysis ERROR: for gradient-based RS, the output ID to be\n");
    printOutTS(PL_ERROR,
       "  analyzed must be 1. (You set it to %d)\n",outputID+1);
    exit(1);
  }
  
  VecYLocal.setLength(nSamples);
  for (ss = 0; ss < nSamples; ss++) 
    VecYLocal[ss] = arrayY[ss*nOutputs+outputID];
  ymax = 0.0;
  ymin = PSUADE_UNDEFINED;
  for (ss = 0; ss < nSamples; ss++)
  {
    if (VecYLocal[ss] > ymax) ymax = VecYLocal[ss];
    if (VecYLocal[ss] < ymin) ymin = VecYLocal[ss];
  }
  printOutTS(PL_INFO,"RSA: Output ID = %d\n", outputID+1);
  printOutTS(PL_INFO,
     "RSA: Output Maximum/Minimum = %14.6e %14.6e\n",ymax,ymin);
  printOutTS(PL_INFO,
     "INFO: Set printlevel higher (1-4) to display more information.\n");
  printOutTS(PL_INFO,
     "INFO: Set print level to 4 for interpolation error graphics file.\n");
  if (ymax == PSUADE_UNDEFINED)
  {
    printOutTS(PL_ERROR,"RSAnalyzer ERROR: some outputs are undefined.\n");
    printOutTS(PL_ERROR,"           Prune them first before analyze.\n");
    return PSUADE_UNDEFINED;
  }

  VecYT.setLength(nSamples);
  for (iL = 0; iL < nLevels; iL++)
  {
    nLast = nSamples;
    if (levelSeps != NULL) nLast = levelSeps[iL];

    faPtr = genFA(rsType_, nInputs, iOne, nLast);
    if (faPtr == NULL)
    {
      printOutTS(PL_INFO,
           "RSFAnalyzer INFO: cannot create response surface.\n");
      delete faPtr;
      return 1.0;
    }
    faPtr->setNPtsPerDim(nPtsPerDim);
    faPtr->setBounds(lower, upper);
    if (iL == nLevels-1) faPtr->setOutputLevel(printLevel);
    if (wgtID >= 0 && wgtID < nOutputs)
    {
      VecWgts.setLength(nLast);
      for (ss = 0; ss < nLast; ss++) 
        VecWgts[ss] = arrayY[ss*nOutputs+wgtID];
      faPtr->loadWeights(nLast, VecWgts.getDVector());
    }
    if (rsType_ == PSUADE_RS_REGRGL)
    {
      VecD.setLength(nSamples*nInputs);
      for (ss = 0; ss < nSamples; ss++)
        for (iI = 0; iI < nInputs; iI++)
          VecD[ss*nInputs+iI] = arrayY[ss*nOutputs+iI+1];
      strcpy(pString, "deriv_sample");
      targv[0] = (char *) pString;
      targv[1] = (char *) VecD.getDVector();
      count = 2;
      faPtr->setParams(count, targv);
    }
    if (iL == nLevels - 1 && psRSCodeGen_ == 2) psRSCodeGen_ = 1;
    status = faPtr->initialize(arrayX, VecYLocal.getDVector());
    if (psRSCodeGen_ == 1) psRSCodeGen_ = 2;
    if (status != 0)
    {
      printOutTS(PL_ERROR,
           "RSAnalysis ERROR: something wrong in FA initialize.\n");
      delete faPtr;
      return PSUADE_UNDEFINED;
    }

    sumErr1  = sumErr2 = maxErr = sumErr1s = sumErr2s = maxErrs = 0.0;
    sumErr11 = sumErr11s = ydiff = 0.0;
    maxBase = maxBases = ymean = yvar = 0.0;
    faPtr->evaluatePoint(nLast, arrayX, VecYT.getDVector());
    for (ss = 0; ss < nLast; ss++)
    {
      ddata = PABS(VecYT[ss] - VecYLocal[ss]);
      if (VecYLocal[ss] != 0.0) sdata = ddata / PABS(VecYLocal[ss]);
      else                      sdata = ddata;
      sumErr1   += ddata;
      sumErr1s  += sdata;
      sumErr11  += (VecYT[ss] - VecYLocal[ss]);
      if (VecYLocal[ss] != 0.0)
           sumErr11s += (VecYT[ss] - VecYLocal[ss])/PABS(VecYLocal[ss]);
      else sumErr11s += (VecYT[ss] - VecYLocal[ss]);
      if (ddata > maxErr ) 
      {
        maxErr = ddata;
        maxBase = PABS(VecYLocal[ss]);
      }
      if (sdata > maxErrs)
      {
        maxErrs = sdata; 
        maxBases = PABS(VecYLocal[ss]);
      }
      sumErr2  += (ddata * ddata);
      sumErr2s += (sdata * sdata);
      ymean += VecYLocal[ss]; 
    }
    ymean /= (double) nLast;
    for (ss = 0; ss < nLast; ss++)
      yvar += (VecYLocal[ss] - ymean) * (VecYLocal[ss] - ymean);
    sumErr1   = sumErr1 / (double) nLast;
    sumErr1s  = sumErr1s / (double) nLast;
    sumErr11  = sumErr11 / (double) nLast;
    sumErr11s = sumErr11s / (double) nLast;
    sumErr2  = sqrt(sumErr2 / (double) nLast);
    sumErr2s = sqrt(sumErr2s / (double) nLast);
    if (printLevel > 2 || iL == nLevels-1)
    {
      printOutTS(PL_INFO,
         "RSAnalysis: L %2d: interpolation error on training set \n", iL);
      printOutTS(PL_INFO,
         "             avg error far from 0 ==> systematic bias.\n");
      printOutTS(PL_INFO,
         "             rms error large      ==> average   error large.\n");
      printOutTS(PL_INFO,
         "             max error large      ==> pointwise error large.\n");
      printOutTS(PL_INFO,
         "             R-square may not always be a reliable measure.\n");
      printOutTS(PL_INFO,"  avg error   = %11.3e (unscaled)\n", sumErr11);
      printOutTS(PL_INFO,"  avg error   = %11.3e (scaled)\n", sumErr11s);
      printOutTS(PL_INFO,"  rms error   = %11.3e (unscaled)\n", sumErr2);
      printOutTS(PL_INFO,"  rms error   = %11.3e (scaled)\n", sumErr2s);
      printOutTS(PL_INFO,"  max error   = %11.3e (unscaled, BASE=%9.3e)\n",
             maxErr, maxBase);
      printOutTS(PL_INFO,
           "  max error   = %11.3e (  scaled, BASE=%9.3e)\n",
           maxErrs, maxBases);
      printOutTS(PL_INFO,
           "  R-square    = %16.8e\n",1.0-sumErr2*sumErr2*nLast/yvar);
      printOutTS(PL_INFO,
           "Based on %d training points (total=%d).\n",nLast,nSamples);
    }

    if (nLevels > 1 && iL < nLevels-1)
    {
      sumErr1 = sumErr2 = maxErr = sumErr1s = sumErr2s = maxErrs = 0.0;
      sumErr11 = sumErr11s;
      faPtr->evaluatePoint(nSamples, arrayX, VecYT.getDVector());
      for (ss = nLast; ss < nSamples; ss++)
      {
        ddata = PABS(VecYT[ss] - VecYLocal[ss]);
        if (VecYLocal[ss] != 0.0) sdata = ddata / PABS(VecYLocal[ss]);
        else                      sdata = ddata;
        sumErr1  += ddata;
        sumErr1s += sdata;
        sumErr11  += (VecYT[ss] - VecYLocal[ss]);
        if (VecYLocal[ss] != 0.0)
             sumErr11s += (VecYT[ss]-VecYLocal[ss])/PABS(VecYLocal[ss]);
        else sumErr11s += (VecYT[ss]-VecYLocal[ss]);
        if (ddata > maxErr )
        {
          maxErr  = ddata; 
          maxBase = PABS(VecYLocal[ss]);
        }
        if (sdata > maxErrs) 
        {
          maxErrs = sdata; 
          maxBases = PABS(VecYLocal[ss]);
        }
        sumErr2  += (ddata * ddata);
        sumErr2s += (sdata * sdata);
      }
      sumErr1   = sumErr1 / (double) (nSamples-nLast);
      sumErr1s  = sumErr1s / (double) (nSamples-nLast);
      sumErr11  = sumErr11 / (double) (nSamples-nLast);
      sumErr11s = sumErr11s / (double) (nSamples-nLast);
      sumErr2   = sqrt(sumErr2 / (double) (nSamples-nLast));
      sumErr2s  = sqrt(sumErr2s / (double) (nSamples-nLast));

      printOutTS(PL_INFO,
         "RSAnalysis: L %2d: Prediction error on the remaining data:\n",iL);
      printOutTS(PL_INFO,"  avg error = %11.3e (unscaled)\n", sumErr11);
      printOutTS(PL_INFO,"  avg error = %11.3e (scaled)\n", sumErr11s);
      printOutTS(PL_INFO,"  rms error = %11.3e (unscaled)\n", sumErr2);
      printOutTS(PL_INFO,"  rms error = %11.3e (scaled)\n", sumErr2s);
      printOutTS(PL_INFO,
         "  max error = %11.3e (unscaled, BASE=%11.3e)\n",maxErr,maxBase);
      printOutTS(PL_INFO,
         "  max error = %11.3e (  scaled, BASE=%11.3e)\n",maxErrs,maxBases);
      printOutTS(PL_INFO, 
         "Based on %d training points (rest=%d).\n",nLast,nSamples-nLast);
    }

    if (psAnaExpertMode_ == 1 && iL == nLevels-1) 
    {
      if (psPlotTool_ == 1)
      {
        fpErr = fopen("RSFA_training_err.sci", "w");
        if (fpErr == NULL)
          printOutTS(PL_INFO,
               "INFO: cannot open file RSFA_training_err.sci.\n");
      }
      else
      {
        fpErr = fopen("RSFA_training_err.m", "w");
        if (fpErr == NULL)
          printOutTS(PL_INFO,
               "INFO: cannot open file RSFA_training_err.m.\n");
      }
      if (fpErr != NULL)
      {
        strcpy(pString, "Surface Fitting Error Histogram");
        fwriteComment(fpErr, pString);
        strcpy(pString, "Interpolation errors on all points");
        fwriteComment(fpErr, pString);
        strcpy(pString, "col 1: interpolated data");
        fwriteComment(fpErr, pString);
        strcpy(pString, "col 2: training data");
        fwriteComment(fpErr, pString);
        strcpy(pString, "col 3: col 1 - col2");
        fwriteComment(fpErr, pString);
        strcpy(pString, "col 4 on: input data");
        fwriteComment(fpErr, pString);
        fprintf(fpErr, "E = [\n");
      }
      sumErr1 = sumErr2 = maxErr = sumErr1s = sumErr2s = maxErrs = 0.0;
      maxBase = sumErr11 = sumErr11s = 0.0;
      faPtr->evaluatePoint(nSamples, arrayX, VecYT.getDVector());
      for (ss = 0; ss < nSamples; ss++)
      {
        ddata = VecYT[ss];
        if (fpErr != NULL)
        {
          fprintf(fpErr, "%24.16e %24.16e %24.16e ", ddata, 
                  VecYLocal[ss], ddata-VecYLocal[ss]);
          for (iI = 0; iI < nInputs; iI++)
            fprintf(fpErr, "%24.16e ", arrayX[ss*nInputs+iI]); 
          fprintf(fpErr, "\n");
        }
        ddata = ddata - VecYLocal[ss];
        sumErr11 += ddata;
        if (VecYLocal[ss] != 0) sumErr11s += ddata/PABS(VecYLocal[ss]);
        else                    sumErr11s += ddata;
        ddata = PABS(ddata);
        if (VecYLocal[ss] != 0) sdata = ddata / PABS(VecYLocal[ss]);
        else                    sdata = ddata;
        sumErr1  += ddata;
        sumErr1s += sdata;
        if (ddata > maxErr )
        {
          maxErr  = ddata; 
          maxBase = PABS(VecYLocal[ss]);
        }
        if (sdata > maxErrs)
        {
          maxErrs = sdata; 
          maxBases = PABS(VecYLocal[ss]);
        }
        sumErr2  += (ddata * ddata);
        sumErr2s += (sdata * sdata);
      }
      if (fpErr != NULL)
      {
        fprintf(fpErr, "];\n");
        fwritePlotCLF(fpErr);
        fwritePlotFigure(fpErr, 1);
        fprintf(fpErr, "subplot(1,2,1)\n");
        fprintf(fpErr, "plot(E(:,3),'x')\n");
        fwritePlotAxes(fpErr);
        fwritePlotTitle(fpErr, "Interpolation Error Plot");
        fwritePlotXLabel(fpErr, "Sample Number");
        fwritePlotYLabel(fpErr, "Interpolation Error");
        fprintf(fpErr, "subplot(1,2,2)\n");
        fprintf(fpErr, "xmax = max(E(:,2));\n");
        fprintf(fpErr, "xmin = min(E(:,2));\n");
        fprintf(fpErr, "ymax = max(E(:,1));\n");
        fprintf(fpErr, "ymin = min(E(:,1));\n");
        fprintf(fpErr, "xmin = min(xmin, ymin);\n");
        fprintf(fpErr, "xmax = max(xmax, ymax);\n");
        fprintf(fpErr, "XX   = xmin : xmax-xmin : xmax;\n");
        fprintf(fpErr, "plot(E(:,2), E(:,1),'x', XX, XX)\n");
        if (psPlotTool_ == 1)
        {
          fprintf(fpErr, "a = get(\"current_axes\");\n");
          fprintf(fpErr, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
        }
        else
        {
          fprintf(fpErr, "axis([xmin xmax xmin xmax])\n");
        }
        fwritePlotAxes(fpErr);
        fwritePlotTitle(fpErr, "Interpolated vs actual data");
        fwritePlotXLabel(fpErr, "Actual data");
        fwritePlotYLabel(fpErr, "Interpolated data");

        fwritePlotFigure(fpErr, 2);
        fwritePlotCLF(fpErr);
        fprintf(fpErr, "subplot(1,2,1)\n");
        if (psPlotTool_ == 1)
        {
          fprintf(fpErr, "ymin = min(E(:,3));\n");
          fprintf(fpErr, "ymax = max(E(:,3));\n");
          fprintf(fpErr, "ywid = 0.1 * (ymax - ymin);\n");
          fprintf(fpErr, "if (ywid < 1.0e-12)\n");
          fprintf(fpErr, "   disp('range too small.')\n");
          fprintf(fpErr, "   halt\n");
          fprintf(fpErr, "end;\n");
          fprintf(fpErr, "histplot(10, E(:,3), style=2);\n");
          fprintf(fpErr, "a = gce();\n");
          fprintf(fpErr, "a.children.fill_mode = \"on\";\n");
          fprintf(fpErr, "a.children.thickness = 2;\n");
          fprintf(fpErr, "a.children.foreground = 0;\n");
          fprintf(fpErr, "a.children.background = 2;\n");
        }
        else
        {
          fprintf(fpErr, "[nk,xk]=hist(E(:,3),10);\n");
          fprintf(fpErr, "bar(xk,nk/%d,1.0)\n",nSamples);
        }
        fwritePlotAxes(fpErr);
        fwritePlotTitle(fpErr, "Interpolation Errors Histogram");
        fwritePlotXLabel(fpErr, "Error");
        fwritePlotYLabel(fpErr, "Probabilities");
        fprintf(fpErr, "subplot(1,2,2)\n");
        if (psPlotTool_ == 1)
        {
          fprintf(fpErr, "ymin = min(E(:,3)/E(:,2));\n");
          fprintf(fpErr, "ymax = max(E(:,3)/E(:,2));\n");
          fprintf(fpErr, "ywid = 0.1 * (ymax - ymin);\n");
          fprintf(fpErr, "if (ywid < 1.0e-12)\n");
          fprintf(fpErr, "   disp('range too small.')\n");
          fprintf(fpErr, "   halt\n");
          fprintf(fpErr, "end;\n");
          fprintf(fpErr, "histplot(10, E(:,3)/E(:,2), style=2);\n");
          fprintf(fpErr, "a = gce();\n");
          fprintf(fpErr, "a.children.fill_mode = \"on\";\n");
          fprintf(fpErr, "a.children.thickness = 2;\n");
          fprintf(fpErr, "a.children.foreground = 0;\n");
          fprintf(fpErr, "a.children.background = 2;\n");
        }
        else
        {
          fprintf(fpErr, "[nk,xk]=hist(E(:,3)./E(:,2),10);\n");
          fprintf(fpErr, "bar(xk,nk/%d,1.0)\n",nSamples);
        }
        fwritePlotAxes(fpErr);
        fwritePlotTitle(fpErr, 
              "Interpolation Errors Histogram (normalized)");
        fwritePlotXLabel(fpErr, "Error");
        fwritePlotYLabel(fpErr, "Probabilities");

        fwritePlotFigure(fpErr, 3);
        fprintf(fpErr, "nn = %d;\n", nInputs);
        iL = (int) pow(1.0*nInputs, 0.5);
        if (iL*iL < nInputs) iL++;
        if (psPlotTool_ == 1) fprintf(fpErr, "drawlater\n");
        fprintf(fpErr, "for ii = 1 : nn\n");
        fprintf(fpErr, "   subplot(%d,%d,ii)\n", iL, iL);
        fprintf(fpErr, "   plot(E(:,ii+3),E(:,3),'x')\n");
        fwritePlotAxes(fpErr);
        if (psPlotTool_ == 1)
        {
          fwritePlotTitle(fpErr,"Error Plot for Input");
          sprintf(winput1, 
               "a.title.text = \"Error Plot for Input\" + string(ii);\n");
          fprintf(fpErr, "%s", winput1);
          fwritePlotXLabel(fpErr, "Input Values");
        }
        else
        {
          fprintf(fpErr, "title(['Error Plot for Input ',int2str(ii)])\n");
          fwritePlotXLabel(fpErr, "Input Values");
        }
        fwritePlotYLabel(fpErr, "Interpolation Error");
        if (psPlotTool_ == 1) fprintf(fpErr, "drawnow\n");
        fprintf(fpErr, "end\n");

        if (nInputs > 2)
        {
          fwritePlotFigure(fpErr, 4);
          if (psPlotTool_ == 1)
          {
            fprintf(fpErr, "f = gcf();\n");
            fprintf(fpErr, "f.color_map = jetcolormap(3);\n");
            fprintf(fpErr, "drawlater\n");
            fprintf(fpErr, "param3d1([E(:,4)' ; E(:,4)'],");
            fprintf(fpErr, "[E(:,5)' ; E(:,5)'],[E(:,3)' ; E(:,3)'])\n");
            fprintf(fpErr, "e = gce();\n");
            fprintf(fpErr, "e.children.mark_mode = \"on\";\n");
            fprintf(fpErr, "e.children.mark_size_unit = \"point\";\n");
            fprintf(fpErr, "e.children.mark_style = 10;\n");
            fprintf(fpErr, "e.children.mark_size = 6;\n");
            fprintf(fpErr, "for i = 1:length(e.children)\n");
            fprintf(fpErr, "   e.children(i).mark_foreground = 1;\n");
            fprintf(fpErr, "end\n");
            fprintf(fpErr, "set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fpErr, "drawnow\n");
          }
          else
          {
            fprintf(fpErr, "   plot3(E(:,4),E(:,5),E(:,3),'bp')\n");
          }
          fwritePlotXLabel(fpErr, "Input 1");
          fwritePlotYLabel(fpErr, "Input 2");
          fwritePlotTitle(fpErr, "Output Error Scatter Plot");
          fwritePlotAxes(fpErr);
        }
        fclose(fpErr);
        if (psPlotTool_ == 1)
           printOutTS(PL_INFO,
                "Interpolation error info are in RSFA_training_err.sci\n");
        else
           printOutTS(PL_INFO,
                "Interpolation error info are in RSFA_training_err.m\n");
      }
      sumErr1   = sumErr1 / (double) nSamples;
      sumErr1s  = sumErr1s / (double) nSamples;
      sumErr11  = sumErr11 / (double) nSamples;
      sumErr11s = sumErr11s / (double) nSamples;
      sumErr2   = sqrt(sumErr2);
      sumErr2s  = sqrt(sumErr2s);
      retdata = sumErr1 / ymax;
    }
    delete faPtr;
  }

  nSubSamples = nSamples / numCVGroups_;
#if 0
  if (psConfig_ != NULL)
  {
    cString = psConfig_->getParameter("RSFA_cross_validation");
    if (cString != NULL)
    {
      useCV_ = 1;
      printOutTS(PL_INFO, 
           "RSFA: turn on cross validation (from config file)\n");
    }
    cString = psConfig_->getParameter("RSFA_cv_ngroups");
    if (cString != NULL)
    {
      sscanf(cString, "%s %s %d",winput1,winput2,&ss);
      if (ss > 0)
      {
        printOutTS(PL_INFO, "RSFA: number of CV groups = %d\n",ss);
        numCVGroups_ = ss;
        nSubSamples = nSamples / numCVGroups_;
      }
      else
      {
        printOutTS(PL_INFO, "RSFA: invalid number of CV groups = %d\n",ss);
        useCV_ = 0;
      }
    }
  }
  else
#endif
  if (rsType_ == PSUADE_RS_SPLINES || rsType_ == PSUADE_RS_REGSG)
  {
    printOutTS(PL_INFO,"Cross validation (CV) cannot be performed on\n");
    printOutTS(PL_INFO,"    splines and sparse grid response surface.\n");
  }
  else
  {
    printAsterisks(PL_INFO, 0);
    printOutTS(PL_INFO,
       "Next you will be asked whether to do cross validation or not.\n");
    printOutTS(PL_INFO,
       "Since cross validation iterates as many times as the number\n");
    printOutTS(PL_INFO,
       "of groups. The rs_expert mode will be turned off. To change\n");
    printOutTS(PL_INFO,
       "the default parameters for different response surface, you\n");
    printOutTS(PL_INFO,
       "will need to exit, create a config file (use genconfigfile\n");
    printOutTS(PL_INFO,
       "in command line mode), and set config option in your data file.\n");
    printDashes(PL_INFO,0);
    sprintf(pString, "Perform cross validation ? (y or n) ");
    getString(pString, winput1);
    if (winput1[0] == 'y')
    {
      useCV_ = 1;
      sprintf(pString, "Enter the number of groups to validate : (2 - %d) ",
               nSamples);
      ss = getInt(1, nSamples, pString);
      printOutTS(PL_INFO, "RSFA: number of CV groups = %d\n",ss);
      numCVGroups_ = ss;
      nSubSamples = nSamples / numCVGroups_;
      if (nSubSamples * numCVGroups_ < nSamples)
      {
         nSubSamples++;
         //numCVGroups_++;
         //printOutTS(PL_INFO,"INFO: number of CV groups adjusted to %d.\n",
         //           numCVGroups_);
         printOutTS(PL_INFO,"      Each CV group has <= %d sample points\n",
                    nSubSamples);
      }
    }
  }

  if (useCV_ == 1)
  {
    printOutTS(PL_INFO,"RSAnalysis: L %2d:cross validation (CV) begins...\n",
               nLevels);

    adata.sampleErrors_ = new double[nSamples];
    for (ss = 0; ss < nSamples; ss++) adata.sampleErrors_[ss] = 0.0;

    faPtr = genFA(rsType_, nInputs, iOne, nSamples-nSubSamples);
    if (faPtr == NULL)
    {
      printOutTS(PL_INFO,
           "RSAnalysis: genFA returned NULL in file %s line %d\n",
           __FILE__, __LINE__ );
      exit(1);
    }
    faPtr->setNPtsPerDim(nPtsPerDim);
    faPtr->setBounds(lower, upper);
    faPtr->setOutputLevel(printLevel);
    rsState = psRSExpertMode_;
    psRSExpertMode_ = 0;

    VecXX.setLength(nSamples*nInputs);
    arrayXX = VecXX.getDVector();
    VecYY.setLength(nSamples);
    IVec1.setLength(nSamples);
    IVec2.setLength(nSamples);
    VecS2.setLength(nSamples);
    VecW.setLength(nSamples);
    VecE.setLength(nSamples);
    VecS.setLength(nSamples);
    VecSigmas.setLength(nSamples);
    if (rsType_ == PSUADE_RS_REGRGL) VecDD.setLength(nSamples*nInputs);
    sprintf(pString, "Random selection of leave-out groups ? (y or n) ");
    getString(pString, winput1);
    if (winput1[0] == 'y')
    {
      generateRandomIvector(nSamples, IVec1.getIVector());
    }
    else
    {
      for (ss = 0; ss < nSamples; ss++) IVec1[ss] = ss;
    }
    for (iI = 0; iI < nInputs; iI++)
    {
      for (ss = 0; ss < nSamples; ss++)
        arrayXX[IVec1[ss]*nInputs+iI] = arrayX[ss*nInputs+iI];
    }
    for (ss = 0; ss < nSamples; ss++) 
      VecYY[IVec1[ss]] = arrayY[ss*nOutputs+outputID];
    if (rsType_ == PSUADE_RS_REGRGL)
    {
      VecDD.setLength(nSamples*nInputs);
      for (ss = 0; ss < nSamples; ss++) 
        for (iI = 0; iI < nInputs; iI++) 
          VecDD[IVec1[ss]*nInputs+iI] = arrayY[ss*nOutputs+iI+1];
    }
    if (wgtID >= 0 && wgtID < nOutputs)
    {
      for (ss = 0; ss < nSamples; ss++)
        VecW[IVec1[ss]] = arrayY[ss*nOutputs+wgtID];
    }
    for (ss = 0; ss < nSamples; ss++) IVec2[IVec1[ss]] = ss;

    VecX2.setLength(nSamples*nInputs);
    VecY2.setLength(nSamples);
    VecD2.setLength(nSamples*nInputs);
    VecWgts.setLength(nSamples);
    CVErr1 = CVErr1s = CVErr2 = CVErr2s = CVMax = CVMaxs = 0.0;
    cvMaxBase = cvMaxBases = 0.0;
    for (ss = 0; ss < nSamples; ss+=nSubSamples)
    {
      printOutTS(PL_INFO,
           "RSAnalysis:: L %2d:CV processes %d out of %d.\n",
           nLevels,ss/nSubSamples+1, numCVGroups_);

      count = 0;
      for (ss2 = 0; ss2 < nSamples; ss2++)
      {
        if (ss2 < ss || ss2 >= (ss+nSubSamples))
        {
          for (iI = 0; iI < nInputs; iI++)
            VecX2[count*nInputs+iI] = VecXX[ss2*nInputs+iI];
          if (wgtID >= 0 && wgtID < nOutputs)
            VecWgts[count] = VecW[ss2*nOutputs+wgtID];
          VecY2[count] = VecYY[ss2];
          if (rsType_ == PSUADE_RS_REGRGL)
          {
            for (iI = 0; iI < nInputs; iI++)
              VecD2[count*nInputs+iI] = VecDD[ss2*nInputs+iI];
          }
          count++;
        }
      }
      if (wgtID >= 0 && wgtID < nOutputs)
        faPtr->loadWeights(nSamples-nSubSamples, VecWgts.getDVector());
      if (rsType_ == PSUADE_RS_REGRGL)
      {
        strcpy(pString, "deriv_sample");
        targv[0] = (char *) pString;
        targv[1] = (char *) VecD2.getDVector();
        count = 2;
        faPtr->setParams(count, targv);
      }
      status = faPtr->initialize(VecX2.getDVector(), VecY2.getDVector());
      if (status == -1) break;
      count = nSubSamples;
      if ((ss + nSubSamples) > nSamples) count = nSamples - ss;
      faPtr->evaluatePointFuzzy(count, &(arrayXX[ss*nInputs]), 
                                VecYT.getDVector(), VecS2.getDVector());
      //faPtr->evaluatePoint(count, &(arrayXX[ss*nInputs]),
      //VecYT.getDVector());

      cvErr1 = cvErr2 = cvErr1s = cvErr2s = cvMax = cvMaxs = 0.0;
      for (ss2 = 0; ss2 < count; ss2++)
      {
        ddata = VecYT[ss2] - VecYY[ss+ss2];
        VecE[IVec2[ss+ss2]] = ddata;
        VecS[IVec2[ss+ss2]] = VecYT[ss2];
        VecSigmas[IVec2[ss+ss2]] = VecS2[ss2];
        cvErr1  += ddata;
        cvErr2  += (ddata * ddata);
        if (PABS(ddata) > cvMax)
        {
          cvMax  = PABS(ddata);
          cvMaxBase = PABS(VecYY[ss+ss2]);
        }
        if (VecYY[ss+ss2] != 0.0) sdata = ddata / PABS(VecYY[ss+ss2]);
        else                   sdata = ddata;
        cvErr1s += sdata;
        cvErr2s += (sdata * sdata);
        if (PABS(sdata) > cvMaxs)
        {
          cvMaxs = PABS(sdata);
          cvMaxBases = PABS(VecYY[ss+ss2]);
        }
        if (printLevel > 4) 
          printOutTS(PL_INFO, 
               "Sample %6d: predicted =  %e, actual =  %e\n",
               IVec2[IVec1[ss+ss2]], VecYT[ss2], VecYY[ss+ss2]);
      }
      adata.sampleErrors_[ss/nSubSamples] = cvErr2s;
      CVErr1  += cvErr1;
      CVErr1s += cvErr1s;
      CVErr2  += cvErr2;
      CVErr2s += cvErr2s;
      cvErr1  = cvErr1 / (double) nSubSamples;
      cvErr1s = cvErr1s / (double) nSubSamples;
      cvErr2  = sqrt(cvErr2 / nSubSamples);
      cvErr2s = sqrt(cvErr2s / nSubSamples);
      if (cvMax > CVMax ) {CVMax = cvMax; CVMaxBase = cvMaxBase;}
      if (cvMaxs > CVMaxs ) {CVMaxs = cvMaxs; CVMaxBases = cvMaxBases;}

      printOutTS(PL_INFO,"RSA: first member of sample group %5d = %d (%d)\n",
             ss/nSubSamples+1, IVec1[ss]+1, nSamples-nSubSamples);
      printOutTS(PL_INFO,
           "RSA: CV error for sample group %5d = %11.3e (avg unscaled)\n",
           ss/nSubSamples+1, cvErr1);
      printOutTS(PL_INFO,
           "RSA: CV error for sample group %5d = %11.3e (avg scaled)\n",
           ss/nSubSamples+1, cvErr1s);
      printOutTS(PL_INFO,
           "RSA: CV error for sample group %5d = %11.3e (rms unscaled)\n",
           ss/nSubSamples+1, cvErr2);
      printOutTS(PL_INFO,
           "RSA: CV error for sample group %5d = %11.3e (rms scaled)\n",
           ss/nSubSamples+1, cvErr2s);
      printOutTS(PL_INFO,
           "RSA: CV error for sample group %5d = %11.3e (max",
           ss/nSubSamples+1, cvMax);
      printOutTS(PL_INFO," unscaled,BASE=%9.3e)\n", cvMaxBase);
      printOutTS(PL_INFO,
           "RSA: CV error for sample group %5d = %11.3e (max",
           ss/nSubSamples+1, cvMaxs);
      printOutTS(PL_INFO,"   scaled,BASE=%9.3e)\n", cvMaxBases);
    }
    if (status >= 0)
    {
      CVErr1  = CVErr1 / (double) nSamples;
      CVErr1s = CVErr1s / (double) nSamples;
      CVErr2  = sqrt(CVErr2 / nSamples);
      CVErr2s = sqrt(CVErr2s / nSamples);
      printOutTS(PL_INFO,"RSA: final CV error  = %11.3e (avg unscaled)\n",
                 CVErr1);
      printOutTS(PL_INFO,"RSA: final CV error  = %11.3e (avg   scaled)\n",
                 CVErr1s);
      printOutTS(PL_INFO,"RSA: final CV error  = %11.3e (rms unscaled)\n",
                 CVErr2);
      printOutTS(PL_INFO,"RSA: final CV error  = %11.3e (rms   scaled)\n",
                 CVErr2s);
      printOutTS(PL_INFO,
           "RSA: final CV error  = %11.3e (max unscaled, BASE=%9.3e)\n",
           CVMax, CVMaxBase);
      printOutTS(PL_INFO,
           "RSA: final CV error  = %11.3e (max   scaled, BASE=%9.3e)\n",
           CVMaxs, CVMaxBases);
      printOutTS(PL_INFO,
           "RSA: L %2d:cross validation (CV) completed.\n",nLevels);
      if (psPlotTool_ == 1)
      {
        fpData = fopen("RSFA_CV_err.sci", "w");
        if (fpData != NULL)
          printOutTS(PL_INFO,"INFO: cannot open file RSFA_CV_err.sci.\n");
      }
      else
      {
        fpData = fopen("RSFA_CV_err.m", "w");
        if (fpData == NULL)
          printOutTS(PL_INFO,"INFO: cannot open file RSFA_CV_err.m.\n");
      }
      if (fpData != NULL)
      {
        strcpy(pString,"This file stores CV error for each point");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 1: error (true - predicted)");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 2: true");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 3: predicted");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 4: standard deviations"); 
        fwriteComment(fpData, pString);
        strcpy(pString,"Set morePlots=1 for normalized residual error plot");
        fwriteComment(fpData, pString);
        ssum = 0.0;
        fprintf(fpData, "morePlots = 0;\n");
        fprintf(fpData, "A = [\n");
        for (ss = 0; ss < nSamples; ss++)
        {
          fprintf(fpData, "  %e %e %e %e\n", VecE[ss], VecYLocal[ss],
                  VecS[ss], VecSigmas[ss]);
          ssum += PABS(VecSigmas[ss]);
        }
        fprintf(fpData, "];\n");
        fwriteHold(fpData, 0);
        fwritePlotFigure(fpData, 1);
        fprintf(fpData, "subplot(1,2,1)\n");
        if (psPlotTool_ == 1)
        {
          fprintf(fpData, "ymin = min(A(:,1));\n");
          fprintf(fpData, "ymax = max(A(:,1));\n");
          fprintf(fpData, "ywid = 0.1 * (ymax - ymin);\n");
          fprintf(fpData, "if (ywid < 1.0e-12)\n");
          fprintf(fpData, "   disp('range too small.')\n");
          fprintf(fpData, "   halt\n");
          fprintf(fpData, "end;\n");
          fprintf(fpData, "histplot(10, A(:,1), style=2);\n");
          fprintf(fpData, "a = gce();\n");
          fprintf(fpData, "a.children.fill_mode = \"on\";\n");
          fprintf(fpData, "a.children.thickness = 2;\n");
          fprintf(fpData, "a.children.foreground = 0;\n");
          fprintf(fpData, "a.children.background = 2;\n");
        }
        else
        {
          fprintf(fpData, "[nk, xk] = hist(A(:,1), 10);\n");
          fprintf(fpData, "bar(xk,nk/sum(nk))\n");
        }
        fwritePlotAxes(fpData);
        fwritePlotTitle(fpData, "CV Error Histogram");
        fwritePlotXLabel(fpData, "Error (unnormalized)");
        fwritePlotYLabel(fpData, "Probabilities");
        if (psPlotTool_ != 1)
        {
          fprintf(fpData,
              "disp(['Error Mean  = ' num2str(mean(A(:,1)))])\n");
          fprintf(fpData,
              "disp(['Error stdev = ' num2str(std(A(:,1)))])\n");
        }
        fprintf(fpData, "subplot(1,2,2)\n");
        fprintf(fpData, "xmax = max(A(:,2));\n");
        fprintf(fpData, "xmin = min(A(:,2));\n");
        fprintf(fpData, "if (xmax == xmin) \n");
        fprintf(fpData, "   xmin = 0.9 * xmin; \n");
        fprintf(fpData, "   xmax = 1.1 * xmax; \n");
        fprintf(fpData, "end;\n");
        fprintf(fpData, "if (xmax == xmin) \n");
        fprintf(fpData, "   xmin = -0.1; \n");
        fprintf(fpData, "   xmax = 0.1; \n");
        fprintf(fpData, "end;\n");
        fprintf(fpData, "ymax = max(A(:,3)+A(:,4));\n");
        fprintf(fpData, "ymin = min(A(:,3)-A(:,4));\n");
        fprintf(fpData, "if (ymax == ymin) \n");
        fprintf(fpData, "   ymin = 0.9 * ymin; \n");
        fprintf(fpData, "   ymax = 1.1 * ymax; \n");
        fprintf(fpData, "end;\n");
        fprintf(fpData, "if (ymax == ymin) \n");
        fprintf(fpData, "   ymin = -0.1; \n");
        fprintf(fpData, "   ymax = 0.1; \n");
        fprintf(fpData, "end;\n");
        fprintf(fpData, "xmin = min(xmin, ymin);\n");
        fprintf(fpData, "xmax = max(xmax, ymax);\n");
        fprintf(fpData, "XX = xmin : xmax-xmin : xmax;\n");
        fprintf(fpData, "plot(A(:,2), A(:,3),'*','MarkerSize',12)\n");
        fwriteHold(fpData, 1);
        fprintf(fpData, "plot(XX, XX)\n");
        if (psPlotTool_ == 1)
        {
          fprintf(fpData, "a = get(\"current_axes\");\n");
          fprintf(fpData, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
        }
        else fprintf(fpData, "axis([xmin xmax xmin xmax])\n");
        if (ssum > 0.0)
        {
          if (psPlotTool_ == 1) fprintf(fpData,"drawlater\n");
          fprintf(fpData,"cnt1 = 0;\n");
          fprintf(fpData,"cnt2 = 0;\n");
          fprintf(fpData,"cnt3 = 0;\n");
          fprintf(fpData,"cnt4 = 0;\n");
          fprintf(fpData,"for ii = 1 : %d\n", nSamples);
          fprintf(fpData,"  xx = [A(ii,2) A(ii,2)];\n");
          fprintf(fpData,"  d1 = A(ii,3)-A(ii,4);\n");
          fprintf(fpData,"  d2 = A(ii,3)+A(ii,4);\n");
          fprintf(fpData,"  yy = [d1 d2];\n");
          fprintf(fpData,"  if (xx(1) < d1 | xx(1) > d2)\n");
          fprintf(fpData,"    plot(xx, yy, 'r-', 'lineWidth', 1)\n");
          fprintf(fpData,"    d3 = A(ii,3)-2*A(ii,4);\n");
          fprintf(fpData,"    d4 = A(ii,3)+2*A(ii,4);\n");
          fprintf(fpData,"    if (xx(1) < d3 | xx(1) > d4)\n");
          fprintf(fpData,"      d5 = A(ii,3)-3*A(ii,4);\n");
          fprintf(fpData,"      d6 = A(ii,3)+3*A(ii,4);\n");
          fprintf(fpData,"      if (xx(1) > d5 & xx(1) < d6)\n");
          fprintf(fpData,"        cnt3 = cnt3 + 1;\n");
          fprintf(fpData,"      else\n");
          fprintf(fpData,"        d7 = A(ii,3)-4*A(ii,4);\n");
          fprintf(fpData,"        d8 = A(ii,3)+4*A(ii,4);\n");
          fprintf(fpData,"        if (xx(1) > d7 & xx(1) < d8)\n");
          fprintf(fpData,"          cnt4 = cnt4 + 1;\n");
          fprintf(fpData,"        else\n");
          if (psPlotTool_ == 1)
            fprintf(fpData,
             "          disp('Point outside 4-sigma = ',ii)\n");
          else
            fprintf(fpData,
             "          disp(['Point outside 4-sigma = ',int2str(ii)])\n");
          fprintf(fpData,"        end;\n");
          fprintf(fpData,"      end;\n");
          fprintf(fpData,"    else\n");
          fprintf(fpData,"      cnt2 = cnt2 + 1;\n");
          fprintf(fpData,"    end;\n");
          fprintf(fpData,"  else\n");
          fprintf(fpData,"    plot(xx, yy, 'g-', 'lineWidth', 1);\n");
          fprintf(fpData,"    cnt1 = cnt1 + 1;\n");
          fprintf(fpData,"  end;\n");
          if (psPlotTool_ == 1)
          {
            fprintf(fpData,"//plot(xx(1), yy(1),'rv','markerSize',10);\n");
            fprintf(fpData,"//plot(xx(2), yy(2),'r^','markerSize',10);\n");
          }
          else
          {
            fprintf(fpData,"%%plot(xx(1), yy(1),'rv','markerSize',10);\n");
            fprintf(fpData,"%%plot(xx(2), yy(2),'r^','markerSize',10);\n");
          }
          fprintf(fpData,"end;\n");
          if (psPlotTool_ == 1) fprintf(fpData,"drawnow\n");
          else
          {
            fprintf(fpData,"cnt4 = cnt1 + cnt2 + cnt3 + cnt4;\n");
            fprintf(fpData,"cnt3 = cnt1 + cnt2 + cnt3;\n");
            fprintf(fpData,"cnt2 = cnt1 + cnt2;\n");
            fprintf(fpData,"disp('Total sample size = %d')\n",nSamples);
            fprintf(fpData,
               "disp(['Num points inside 1 sigma = ',int2str(cnt1)])\n");
            fprintf(fpData,
               "disp(['Num points inside 2 sigma = ',int2str(cnt2)])\n");
            fprintf(fpData,
               "disp(['Num points inside 3 sigma = ',int2str(cnt3)])\n");
            fprintf(fpData,
               "disp(['Num points inside 4 sigma = ',int2str(cnt4)])\n");
            fprintf(fpData,
                "text(0.1,0.9,'RED: prediction outside +/- 1 std ");
            fprintf(fpData,
                "dev','sc','fontSize',11,'fontweight','bold')\n");
          }
        }
        fwriteHold(fpData, 0);
        fwritePlotAxes(fpData);
        sprintf(pString,
           "CV Predictions (scaled rmse=%12.4e, R2=%12.4e)",
           CVErr2s,1.0-CVErr2*CVErr2*nSamples/yvar);
        fwritePlotTitle(fpData, pString);
        fwritePlotXLabel(fpData, "actual data");
        fwritePlotYLabel(fpData, "predicted data");
        fprintf(fpData, "rsme_scaled = %e;\n", CVErr2s);
        fprintf(fpData, "R2 = %e;\n",1.0-CVErr2*CVErr2*nSamples/yvar);
        fprintf(fpData,"if morePlots == 1\n");
        strcpy(pString," For the following B matrix");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 1: true values");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 2: normalized residual");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 3: predicted values");
        fwriteComment(fpData, pString);
        strcpy(pString,"Column 4-(m+3): inputs");
        fwriteComment(fpData, pString);
        fwritePlotFigure(fpData, 2);
        fprintf(fpData, "B = [\n");
        for (ss = 0; ss < nSamples; ss++)
        {
          if (VecYLocal[ss] == 0) 
            fprintf(fpData, " %e 0 ",VecYLocal[ss]);
          else 
            fprintf(fpData," %e %e ",VecYLocal[ss],VecE[ss]/VecYLocal[ss]);
          fprintf(fpData," %e ", VecS[ss]);
          for (ss2 = 0; ss2 < nInputs; ss2++)
            fprintf(fpData," %e ",arrayXX[ss*nInputs+ss2]);
          fprintf(fpData,"\n");
        }
        fprintf(fpData,"];\n");
        fprintf(fpData,"AA = B(:,1);\n");
        fprintf(fpData,"BB = B(:,2);\n");
        fprintf(fpData,"CC = B(:,3);\n");
        fprintf(fpData,
                "subplot(1,2,1),plot(AA, BB, '*','markerSize',12);\n");
        fwritePlotAxes(fpData);
        fwritePlotTitle(fpData, "Normalized Error Analysis");
        fwritePlotXLabel(fpData, "Actual Data");
        fwritePlotYLabel(fpData, "Normalized Error");
        fprintf(fpData,
                "subplot(1,2,2),plot(AA,BB.*AA, '*','markerSize',12);\n");
        fwritePlotAxes(fpData);
        fwritePlotTitle(fpData, "Unnormalized Error Analysis");
        fwritePlotXLabel(fpData, "Actual Data");
        fwritePlotYLabel(fpData, "Error");
        fprintf(fpData,"figure(3)\n");
        fprintf(fpData,"for ii = 1 : %d\n", nInputs);
        fprintf(fpData,"  XX = B(:,3+ii);\n");
        fprintf(fpData,"  plot(XX, BB.*AA, '*','markerSize',12);\n");
        fprintf(fpData,"  disp(['Error vs Input ' int2str(ii)])\n");
        fprintf(fpData,"  disp('Press enter to continue')\n");
        fwritePlotAxes(fpData);
        fprintf(fpData,"  xlabel(['Input ' int2str(ii)])\n");
        fwritePlotYLabel(fpData,"Error");
        fprintf(fpData,"  pause\n");
        fprintf(fpData,"end;\n");
        fprintf(fpData,"end;\n");
        fclose(fpData);
        if (psPlotTool_ == 1)
             printOutTS(PL_INFO, "CV error file is RSFA_CV_err.sci\n");
        else printOutTS(PL_INFO, "CV error file is RSFA_CV_err.m\n");
      }
    }
    psRSExpertMode_ = rsState;
    delete faPtr;
  }

  testFlag = 0;
  if (psConfig_ != NULL)
  {
    cString = psConfig_->getParameter("RSFA_test_datafile");
    if (cString != NULL)
    {
      sscanf(cString, "%s %s %s",winput1,winput2,dataFile);
      fpData = fopen(dataFile, "r");
      if (fpData != NULL)
      {
        testFlag = 1;
        printOutTS(PL_INFO, "RSFA:: test data file = %s\n", dataFile);
        fclose(fpData);
        fpData = NULL;
      }
    }
    if (testFlag == 1)
    {
      cString = psConfig_->getParameter("RSFA_test_errorfile");
      if (cString != NULL)
      {
        sscanf(cString, "%s %s %s",winput1,winput2,errFile);
        fpErr = fopen(errFile, "w");
        if (fpErr != NULL)
        {
          printOutTS(PL_INFO, "RSFA:: error file = %s\n", errFile);
          testFlag += 2;
          fclose(fpErr);
        }
      }
    }
  }
  if (testFlag == 1) retdata = validate(adata, dataFile, NULL);
  if (testFlag == 3) retdata = validate(adata, dataFile, errFile);
  printAsterisks(PL_INFO, 0);

  return retdata;
}

// ************************************************************************ 
// set parameters
// ------------------------------------------------------------------------
int RSFuncApproxAnalyzer::setParams(int argc, char **argv)
{
  int   idata;
  char  *request, *dataFile, *errFile;
  aData *adata;

  request = (char *) argv[0]; 
  if (!strcmp(request, "validate"))
  {
    if (argc != 4) printOutTS(PL_WARN, "RSAnalysis WARNING: setParams.\n");
    adata    = (aData *) argv[1];
    dataFile = (char *) argv[2];
    errFile  = (char *) argv[3];
    validate(*adata, dataFile, errFile);
  }
  else if (!strcmp(request, "rstype"))
  {
    if (argc != 2) printOutTS(PL_WARN, "RSAnalysis WARNING: setParams.\n");
    rsType_ = *(int *) argv[1];
    if (rsType_ < 0 || rsType_ > PSUADE_NUM_RS)
    {
      printOutTS(PL_ERROR, 
           "RSAnalysis ERROR: INVALID rstype, set to REGR2.\n");
      rsType_ = PSUADE_RS_REGR2;
    }
  }
  else if (!strcmp(request, "usecv"))
  {
    useCV_ = 1;
  }
  else if (!strcmp(request, "numcvgroups"))
  {
    idata = *(int *) argv[1];
    if (idata > 1) numCVGroups_ = idata;
  }
  else
  {
    printOutTS(PL_ERROR, "RSAnalysis ERROR: setParams - not valid.\n");
    exit(1);
  }
  return 0;
}

// ************************************************************************ 
// validate
// ------------------------------------------------------------------------
double RSFuncApproxAnalyzer::validate(aData &adata, char *dataFile, 
                                      char *errFile)
{
  int        nInputs, nOutputs, outputID, wgtID, ii, iOne=1,fatype;
  int        nSamples, status, nTestOut, nTestIn, nPtsPerDim=64, nTestSam;
  double     sumErr1, sumErr2, maxErr, sumErr1s, sumErr2s, maxErrs;
  double     dataMax, dataMin, maxBase, maxBases, sdata, ddata, ssum;
  double     *arrayX, *arrayY, *lower, *upper, *arrayXX, *arrayYY;
  char       pString[501];
  pData      pPtr, pInputs, pOutputs;
  PsuadeData *ioPtr;
  FuncApprox *fa;
  psVector   VecYLocal, VecWgts, VecYT, VecStd;
  FILE       *fpErr;

  printOutTS(PL_INFO, "RSAnalysis: validating against a test set ...\n");

  nInputs   = adata.nInputs_;
  nOutputs  = adata.nOutputs_;
  nSamples  = adata.nSamples_;
  lower     = adata.iLowerB_;
  upper     = adata.iUpperB_;
  arrayX    = adata.sampleInputs_;
  arrayY    = adata.sampleOutputs_;
  outputID  = adata.outputID_;
  wgtID     = adata.regWgtID_;
  fatype    = adata.faType_;
  VecYLocal.setLength(nSamples);
  for (ii = 0; ii < nSamples; ii++) 
    VecYLocal[ii] = arrayY[ii*nOutputs+outputID];

  ioPtr = new PsuadeData();
  status = ioPtr->readPsuadeFile(dataFile);
  if (status != 0)
  {
    printf("ERROR: cannot read file %s in PSUADE format.\n",dataFile);
    exit(1);
  } 
  ioPtr->getParameter("output_noutputs", pPtr);
  nTestOut = pPtr.intData_;
  ioPtr->getParameter("method_nsamples", pPtr);
  nTestSam = pPtr.intData_;
  ioPtr->getParameter("input_ninputs", pPtr);
  nTestIn = pPtr.intData_;
  if (nTestSam <= 0)
  {
    printOutTS(PL_WARN, "RSAnalysis: file %s has no data.\n", dataFile);
    delete ioPtr;
    return 1.0e12;
  }
  if (nTestIn != nInputs)
  {
    printOutTS(PL_WARN, 
         "RSAnalysis: test data has different number of inputs (%d)\n",
         nTestIn);
    printOutTS(PL_WARN,
         "            than that of the sample (%d)\n", nInputs);
    delete ioPtr;
    return 1.0e12;
  }
  ioPtr->getParameter("input_sample",pInputs);
  arrayXX = pInputs.dbleArray_;
  ioPtr->getParameter("output_sample",pOutputs);
  arrayYY = pOutputs.dbleArray_;

  fa = genFA(fatype, nInputs, iOne, nSamples);
  if (fa == NULL)
  {
    printOutTS(PL_INFO, 
         "RSAnalysis INFO: cannot create function approximator.\n");
    delete ioPtr;
    return 1.0e12;
  }
  fa->setNPtsPerDim(nPtsPerDim);
  fa->setBounds(adata.iLowerB_, adata.iUpperB_);
  fa->setOutputLevel(adata.printLevel_);
  if (wgtID >= 0 && wgtID < nOutputs)
  {
    VecWgts.setLength(nSamples);
    for (ii = 0; ii < nSamples; ii++) 
      VecWgts[ii] = arrayY[ii*nOutputs+wgtID];
    fa->loadWeights(nSamples, VecWgts.getDVector());
  }
  status = fa->initialize(arrayX, VecYLocal.getDVector());

  if (errFile != NULL) 
  {
    fpErr = fopen(errFile, "w");
    if (fpErr != NULL)
    {
      strcpy(pString, "Surface Fitting Error Histogram");
      fwriteComment(fpErr, pString);
      strcpy(pString, "Interpolation errors on all points");
      fwriteComment(fpErr, pString);
      strcpy(pString, "Column 1: predicted data");
      fwriteComment(fpErr, pString);
      strcpy(pString, "Column 2: actual data");
      fwriteComment(fpErr, pString);
      strcpy(pString, "Column 3: error = predicted - actual data");
      fwriteComment(fpErr, pString);
      strcpy(pString, "Column 4: normalized error");
      fwriteComment(fpErr, pString);
      strcpy(pString, "Column 5: predicted standard deviation");
      fwriteComment(fpErr, pString);
      fprintf(fpErr, "E = [\n");
    }
  }
  else fpErr = NULL;

  dataMax = -PSUADE_UNDEFINED;
  dataMin = PSUADE_UNDEFINED;
  for (ii = 0; ii < nTestSam; ii++)
  {
    ddata = PABS(arrayYY[ii*nTestOut+outputID]);
    if (ddata > dataMax) dataMax = ddata;
    if (ddata < dataMin) dataMin = ddata;
  }
  VecYT.setLength(nTestSam);
  VecStd.setLength(nTestSam);
  fa->evaluatePointFuzzy(nTestSam,arrayXX,VecYT.getDVector(),
                         VecStd.getDVector());
  sumErr1 = sumErr1s = sumErr2 = sumErr2s = maxErr = maxErrs = 0.0;
  for (ii = 0; ii < nTestSam; ii++)
  {
    ddata = VecYT[ii];
    sdata = arrayYY[ii*nTestOut+outputID];
    if (fpErr != NULL) 
    {
      if (sdata != 0)
        fprintf(fpErr, "%16.8e %16.8e %16.8e %16.8e %16.8e\n", ddata, 
           sdata, ddata-sdata, PABS((ddata-sdata)/sdata), VecStd[ii]);
      else
        fprintf(fpErr, "%16.8e %16.8e %16.8e %16.8e %16.8e\n", ddata, 
           sdata, ddata-sdata, PABS(ddata-sdata), VecStd[ii]);
    }
    ddata = ddata - sdata;
    if (sdata != 0.0) sdata = ddata / PABS(sdata);
    else              sdata = ddata;
    sumErr1  += ddata;
    sumErr1s += sdata;
    sumErr2  += (ddata * ddata);
    sumErr2s += (sdata * sdata);
    ddata = PABS(ddata);
    sdata = PABS(sdata);
    if (ddata > maxErr ) 
    {
      maxErr = ddata;
      maxBase = PABS(arrayYY[ii*nTestOut+outputID]);
    }
    if (sdata > maxErrs)
    {
      maxErrs  = sdata;
      maxBases = PABS(arrayYY[ii*nTestOut+outputID]);
    }
  }

  sumErr1  = sumErr1 / (double) nTestSam;
  sumErr1s = sumErr1s / (double) nTestSam;
  sumErr2  = sqrt(sumErr2 / (double) nTestSam);
  sumErr2s = sqrt(sumErr2s / (double) nTestSam);
  printOutTS(PL_INFO,"RSA: Test data maximum/minimum      = %9.3e %9.3e\n",
         dataMax, dataMin);
  printOutTS(PL_INFO, 
       "RSA: Prediction errors = %11.3e (avg unscaled)\n",sumErr1);
  printOutTS(PL_INFO, 
       "RSA: Prediction errors = %11.3e (avg   scaled)\n",sumErr1s);
  printOutTS(PL_INFO, 
       "RSA: Prediction errors = %11.3e (rms unscaled)\n",sumErr2);
  printOutTS(PL_INFO, 
       "RSA: Prediction errors = %11.3e (rms   scaled)\n",sumErr2s);
  printOutTS(PL_INFO, 
       "RSA: Prediction errors = %11.3e (max unscaled, BASE=%9.3e)\n",
        maxErr, maxBase);
  printOutTS(PL_INFO, 
       "RSA: Prediction errors = %11.3e (max   scaled, BASE=%9.3e)\n",
        maxErrs, maxBases);
  ssum = 0.0;
  for (ii = 0; ii < nTestSam; ii++) ssum += VecStd[ii];
  if (ssum != 0)
    printOutTS(PL_INFO, 
         "RSA: Average std. dev. = %11.3e (sum of all points)\n",
         ssum/nSamples);

  if (fpErr != NULL) 
  {
#if 0
    fprintf(fpErr, "];\n");
    fprintf(fpErr, "Emax = max(E(:,3));\n");
    fprintf(fpErr, "Emin = min(E(:,3));\n");
    fprintf(fpErr, "Estp = 0.1 * (Emax - Emin);\n");
    fprintf(fpErr, "X = zeros(11,1);\n");
    fprintf(fpErr, "Y = zeros(10,1);\n");
    fprintf(fpErr, "for i = 1 : 11\n");
    fprintf(fpErr, "   X(i) = Emin + (i - 1) * Estp;\n");
    fprintf(fpErr, "end;\n");
    fprintf(fpErr, "for i = 1 : 10\n");
    fprintf(fpErr, "   inds = find(E(:,3) >= X(i));\n");
    fprintf(fpErr, "   vals = E(inds,3);\n");
    fprintf(fpErr, "   inds = find(vals < X(i+1));\n");
    fprintf(fpErr, "   Y(i) = length(inds);\n");
    fprintf(fpErr, "end;\n");
    fprintf(fpErr, "Y(10) = Y(10) + 1;\n");
    fprintf(fpErr, "bar(X(1:10),Y)\n");
    fprintf("Error file in %s.\n", errFile);
#else
    fprintf(fpErr, "];\n");
    fwritePlotCLF(fpErr);
    fwritePlotFigure(fpErr, 1);
    fprintf(fpErr, "subplot(1,2,1)\n");
    fprintf(fpErr, "plot(E(:,3),'x')\n");
    fwritePlotAxes(fpErr);
    fwritePlotTitle(fpErr, "Error Plot");
    fwritePlotXLabel(fpErr, "Sample Number");
    fwritePlotYLabel(fpErr, "Error");
    fprintf(fpErr, "subplot(1,2,2)\n");
    if (psPlotTool_ == 1)
    {
      fprintf(fpErr, "ymin = min(E(:,3));\n");
      fprintf(fpErr, "ymax = max(E(:,3));\n");
      fprintf(fpErr, "ywid = 0.1 * (ymax - ymin);\n");
      fprintf(fpErr, "if (ywid < 1.0e-12)\n");
      fprintf(fpErr, "   disp('range too small.')\n");
      fprintf(fpErr, "   halt\n");
      fprintf(fpErr, "end;\n");
      fprintf(fpErr, "histplot(10, E(:,3), style=2);\n");
      fprintf(fpErr, "a = gce();\n");
      fprintf(fpErr, "a.children.fill_mode = \"on\";\n");
      fprintf(fpErr, "a.children.thickness = 2;\n");
      fprintf(fpErr, "a.children.foreground = 0;\n");
      fprintf(fpErr, "a.children.background = 2;\n");
    }
    else
    {
      fprintf(fpErr, "[nk,xk]=hist(E(:,3),10);\n");
      fprintf(fpErr, "bar(xk,nk/%d,1.0)\n",nTestSam);
    }
    fwritePlotAxes(fpErr);
    fwritePlotTitle(fpErr, "Histogram of Errors");
    fwritePlotXLabel(fpErr, "Error");
    fwritePlotYLabel(fpErr, "Probabilities");
    fwritePlotFigure(fpErr, 2);
    fprintf(fpErr, "xmax = max(E(:,2));\n");
    fprintf(fpErr, "xmin = min(E(:,2));\n");
    fprintf(fpErr, "if (xmax == xmin) \n");
    fprintf(fpErr, "   xmin = 0.9 * xmin; \n");
    fprintf(fpErr, "   xmax = 1.1 * xmax; \n");
    fprintf(fpErr, "end;\n");
    fprintf(fpErr, "if (xmax == xmin) \n");
    fprintf(fpErr, "   xmin = -0.1; \n");
    fprintf(fpErr, "   xmax = 0.1; \n");
    fprintf(fpErr, "end;\n");
    fprintf(fpErr, "ymax = max(E(:,1)+E(:,5));\n");
    fprintf(fpErr, "ymin = min(E(:,1)-E(:,5));\n");
    fprintf(fpErr, "if (ymax == ymin) \n");
    fprintf(fpErr, "   ymin = 0.9 * ymin; \n");
    fprintf(fpErr, "   ymax = 1.1 * ymax; \n");
    fprintf(fpErr, "end;\n");
    fprintf(fpErr, "if (ymax == ymin) \n");
    fprintf(fpErr, "   ymin = -0.1; \n");
    fprintf(fpErr, "   ymax = 0.1; \n");
    fprintf(fpErr, "end;\n");
    fprintf(fpErr, "xmin = min(xmin, ymin);\n");
    fprintf(fpErr, "xmax = max(xmax, ymax);\n");
    fprintf(fpErr, "plot(E(:,2), E(:,1),'*','MarkerSize',12)\n");
    fwriteHold(fpErr, 1);
    fprintf(fpErr, "XX = xmin : xmax-xmin : xmax;\n");
    fprintf(fpErr, "plot(XX, XX)\n");
    if (psPlotTool_ == 1)
    {
      fprintf(fpErr, "a = get(\"current_axes\");\n");
      fprintf(fpErr, "a.data_bounds=[xmin,xmin;xmax,xmax];\n");
    }
    else
    {
      fprintf(fpErr, "axis([xmin xmax xmin xmax])\n");
    }
    if (ssum > 0.0)
    {
      fprintf(fpErr,"for ii = 1 : %d\n", nTestSam);
      fprintf(fpErr,"  xx = [E(ii,2) E(ii,2)];\n");
      fprintf(fpErr,"  d1 = E(ii,1)-E(ii,5);\n");
      fprintf(fpErr,"  d2 = E(ii,1)+E(ii,5);\n");
      fprintf(fpErr,"  yy = [d1 d2];\n");
      fprintf(fpErr,"  if (xx(1) < d1 | xx(1) > d2);\n");
      fprintf(fpErr,"    plot(xx, yy, 'r-', 'lineWidth', 1);\n");
      fprintf(fpErr,"   else\n");
      fprintf(fpErr,"     plot(xx, yy, 'g-', 'lineWidth', 1);\n");
      fprintf(fpErr,"  end;\n");
      fprintf(fpErr,"%% plot(xx(1), yy(1), 'rv', 'markerSize', 10);\n");
      fprintf(fpErr,"%% plot(xx(2), yy(2), '^', 'markerSize', 10);\n");
      fprintf(fpErr,"end;\n");
    }
    fwriteHold(fpErr, 0);
    if (psPlotTool_ != 1)
    {
      if (ssum > 0.0)
      {
        fprintf(fpErr,"text(0.1,0.9,'RED: prediction outside +/- 1 std");
        fprintf(fpErr," dev','sc','fontSize',11,'fontweight','bold')\n");
      }
    }
    fwritePlotAxes(fpErr);
    fwritePlotTitle(fpErr, "Interpolated versus actual data");
    fwritePlotXLabel(fpErr, "Actual data");
    fwritePlotYLabel(fpErr, "Interpolated data");
#endif
    fclose(fpErr);
    fpErr = NULL;
    printOutTS(PL_INFO, 
         "RSAnalysis: individual prediction errors can be found in %s.\n",
         errFile);
  }

  delete fa;
  delete ioPtr;
  return sumErr1;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
RSFuncApproxAnalyzer& 
RSFuncApproxAnalyzer::operator=(const RSFuncApproxAnalyzer &)
{
  printOutTS(PL_ERROR, 
       "RSAnalysis operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

