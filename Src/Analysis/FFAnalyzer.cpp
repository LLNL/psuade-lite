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
// functions for the class FFAnalyzer (Fractional Factorial analysis)  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include "FFAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
FFAnalyzer::FFAnalyzer() : Analyzer(), nInputs_(0), nSamples_(0), nReps_(0),
                           means_(0), stds_(0), mEffects_(0)
{
   setName("FF");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
FFAnalyzer::~FFAnalyzer()
{
   if (means_) delete [] means_;
   if (stds_) delete [] stds_;
   if (mEffects_) delete [] mEffects_;
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double FFAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, outputID, checkSample, printLevel;
   int    ii, ii2, rr, ss,whichOutput,ncount,nplots, *iArray,nReps, offset;
   double *X, *Y,*txArray,*twArray,*tyArray, **mEffects,accum, ***iEffects;
   double ***lolos, ***lohis, ***hilos, ***hihis, stdev, mean;
   char   pString[500], winput[500], **iNames;
   FILE   *fp=NULL;
   PsuadeData *ioPtr=NULL;
   pData  qData;

   printLevel = adata.printLevel_;
   nInputs    = adata.nInputs_;
   nInputs_   = nInputs;
   nOutputs   = adata.nOutputs_;
   nSamples   = adata.nSamples_;
   nSamples_  = nSamples;
   outputID   = adata.outputID_;
   X          = adata.sampleInputs_;
   Y          = adata.sampleOutputs_;
   ioPtr      = adata.ioPtr_;
   if (adata.inputPDFs_ != NULL)
   {
      ncount = 0;
      for (ii = 0; ii < nInputs; ii++) ncount += adata.inputPDFs_[ii];
      if (ncount > 0)
      {
         printOutTS(PL_INFO, 
              "FFAnalysis INFO: some inputs have non-uniform PDFs, but\n");
         printOutTS(PL_INFO, 
              "           they are not relevant in this analysis.\n");
      }
   }
   whichOutput = outputID;
   if (whichOutput >= nOutputs || whichOutput < 0) whichOutput = 0;
   if (ioPtr != NULL) ioPtr->getParameter("input_names", qData);

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR, "FFAnalysis ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "   nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "   nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nSamples / 2 * 2 != nSamples)
   {
      printOutTS(PL_ERROR, "FFAnalysis ERROR: nSamples has to be even.\n");
      printOutTS(PL_ERROR, "   nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 

   if (means_)    delete [] means_;
   if (stds_)     delete [] stds_;
   if (mEffects_) delete [] mEffects_;
   means_    = NULL;
   stds_     = NULL;
   mEffects_ = NULL;

   txArray = new double[nSamples];
   tyArray = new double[nSamples];
   twArray = new double[nSamples];
   checkAllocate(twArray, "twArray in FF::analyze");

   printOutTS(PL_INFO, "\n");
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"* Fractional Factorial Main Effect Analysis\n");
   printOutTS(PL_INFO,"* This analysis works for FF4, FF5 and PBD designs\n");
   printOutTS(PL_INFO,"* This analysis is good for models that are linear\n");
   printOutTS(PL_INFO,"* with respect to each input with some pairwise\n");
   printOutTS(PL_INFO,"* interactions, e.g. Y = X_1 + X2 + X1 * X2.\n");
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO, "* total number of samples = %10d \n",nSamples);
   printOutTS(PL_INFO, "* number of Inputs        = %10d \n",nInputs);
   printOutTS(PL_INFO, "* Output number           = %d\n", whichOutput+1);
   printDashes(PL_INFO, 0);
   if (psPlotTool_ == 0) fp = fopen("matlabmeff.m", "w");
   else                  fp = fopen("scilabmeff.sci", "w");
   if (fp != NULL)
   {
      if (psPlotTool_ == 1)
      {
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// ********************************************** \n");
         fprintf(fp,"// * Fractional Factorial Main Effect Analysis ** \n");
         fprintf(fp,"// * (2 level fractional factorial only)       ** \n");
         fprintf(fp,"// *-------------------------------------------** \n");
         fprintf(fp,"// Output %d\n", whichOutput+1);
         fprintf(fp,"// *-------------------------------------------** \n");
      }
      else
      {
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% ********************************************** \n");
         fprintf(fp,"%% * Fractional Factorial Main Effect Analysis ** \n");
         fprintf(fp,"%% * (2 level fractional factorial only)       ** \n");
         fprintf(fp,"%% *-------------------------------------------** \n");
         fprintf(fp,"%% Output %d\n", whichOutput+1);
         fprintf(fp,"%% *-------------------------------------------** \n");
      }
   }

   nReps = 1;
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ss = 0; ss < nSamples; ss++)
      {
         txArray[ss] = X[nInputs*ss+ii];
         tyArray[ss] = Y[nOutputs*ss+whichOutput];
      }
      sortDbleList2(nSamples, txArray, tyArray);
      checkSample = 1;
      for (ss = 1; ss < nSamples/2; ss++)
         if (txArray[ss] != txArray[0]) checkSample = 0;
      for (ss = nSamples/2+1; ss < nSamples; ss++)
         if (txArray[ss] != txArray[nSamples/2]) checkSample = 0;
      if (checkSample == 0)
      {
         printOutTS(PL_INFO, 
              "FFAnalysis ERROR: sample not fractional factorial.\n");
         printOutTS(PL_INFO, 
              "If you are using replicated Fractional Factorial\n");
         printOutTS(PL_INFO, "enter the number of replications.\n");
         sprintf(pString, "Number of replications = (2 - %d) ",nSamples/2);
         nReps = getInt(2, nSamples/2, pString);
         break;
      }
   }
   nReps_ = nReps;
   if (nReps > 1)
   {
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = 0; ii2 < nReps; ii2++)
         {
            for (ss = 0; ss < nSamples/nReps; ss++)
            {
               txArray[ss] = X[(ii2*nSamples/nReps+ss)*nInputs+ii];
               tyArray[ss] = Y[(ii2*nSamples/nReps+ss)*nOutputs+whichOutput];
            }
            sortDbleList2(nSamples/nReps, txArray, tyArray); 
            checkSample = 1;
            for (ss = 1; ss < nSamples/nReps/2; ss++)
               if (txArray[ss] != txArray[0]) checkSample = 0;
            for (ss = nSamples/nReps/2+1; ss < nSamples/nReps; ss++)
               if (txArray[ss] != txArray[nSamples/nReps-1])
                  checkSample = 0;
            if (checkSample == 0)
            {
               printOutTS(PL_ERROR, 
                  "FFAnalysis ERROR: sample not fractional factorial.\n");
               printOutTS(PL_ERROR, 
                  "                  nor replicated fractional factorial.\n");
               delete [] twArray;
               delete [] txArray;
               delete [] tyArray;
               if(fp != NULL) fclose(fp);
               return 0.0;
            }
         }
      }
   }

   mEffects = new double*[nInputs];
   checkAllocate(mEffects, "mEffects in FF::analyze");
   for (ii = 0; ii < nInputs; ii++) mEffects[ii] = new double[nReps+1];

   mEffects_ = new double*[nInputs_];
   for (ii = 0; ii < nInputs_; ii++) mEffects_[ii] = new double[nReps_+1];
   means_ = new double[nInputs_];
   stds_  = new double[nInputs_];
   checkAllocate(stds_, "stds_ in FF::analyze");

   iArray  = new int[nInputs];
   for (ii = 0; ii < nInputs; ii++) iArray[ii] = ii;

   for (rr = 0; rr < nReps; rr++)
   {
      offset = rr * nSamples / nReps;
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ss = 0; ss < nSamples/nReps; ss++)
         {
            txArray[ss] = X[(offset+ss)*nInputs+ii];
            tyArray[ss] = Y[(offset+ss)*nOutputs+whichOutput];
         }
         sortDbleList2(nSamples/nReps, txArray, tyArray); 
         accum = 0.0;
         for (ss = nSamples/nReps/2; ss < nSamples/nReps; ss++)
            accum += tyArray[ss];
         for (ss = 0; ss < nSamples/nReps/2; ss++) accum -= tyArray[ss];
         mEffects[ii][rr+1] = accum * 2 / (double) (nSamples/nReps);

         //save main effects
         mEffects_[ii][rr] = mEffects[ii][rr+1];
      }
   }
   printOutTS(PL_INFO,"* Fractional Factorial Main Effect (normalized)\n");
   printOutTS(PL_INFO,"* Note: std err is the standard error or mean.\n");
   printDashes(PL_INFO, 0);
   for (ii = 0; ii < nInputs; ii++)
   {
      mean = 0.0;
      for (rr = 0; rr < nReps; rr++) mean += mEffects[ii][rr+1];
      mean /= (double) nReps;
      mEffects[ii][0] = mean;
      stdev = 0.0;
      if (nReps > 1)
      {
         for (rr = 0; rr < nReps; rr++)
            stdev += pow(mEffects[ii][rr+1]-mean, 2.0);
         stdev = sqrt(stdev / (double) (nReps - 1.0));
         printOutTS(PL_INFO, 
              "* Input %3d  =  %12.4e (std err = %12.4e)\n", ii+1,mean,
              stdev/sqrt(1.0*nReps));
      }  
      else printOutTS(PL_INFO, "* Input %3d  =  %12.4e\n", ii+1, mean);

      //save means & stds
      means_[ii] = mean;
      stds_[ii]  = stdev/sqrt(1.0*nReps);
   }
   printDashes(PL_INFO, 0);
   printOutTS(PL_INFO,"* Fractional Factorial Main Effect (ordered)\n");
   printDashes(PL_INFO, 0);
   for (ii = 0; ii < nInputs; ii++) twArray[ii] = PABS(mEffects[ii][0]);
   sortDbleList2a(nInputs, twArray, iArray);
   for (ii = nInputs-1; ii >= 0; ii--)
      printOutTS(PL_INFO, 
             "* Rank %3d : Input %3d (measure = %12.4e)\n", nInputs-ii,
             iArray[ii]+1, twArray[ii]);

   if (fp != NULL)
   {
      fprintf(fp, "Y = [\n");
      for (ii = 0; ii < nInputs; ii++)
         fprintf(fp, "%24.16e \n", mEffects[ii][0]);
      fprintf(fp, "]; \n");
      fwritePlotCLF(fp);
      fprintf(fp, "bar(Y,0.8); \n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Main Effect Plot for Output");
      fwritePlotXLabel(fp, "Input Number");
      fwritePlotYLabel(fp, "Average Gradient");
      fprintf(fp, "disp('Press enter to continue to the next plot')\n");
   }

   if (nInputs < 2 || adata.samplingMethod_ == PSUADE_SAMP_PBD)
   {
      delete [] txArray;
      delete [] twArray;
      delete [] tyArray;
      delete [] iArray;
      for(int i = 0; i < nInputs; i++) delete [] mEffects[i];
      delete [] mEffects;
      if (fp != NULL) fclose(fp);
      return 0.0;
   }

   iEffects = new double**[nInputs];
   hihis = new double**[nInputs];
   lolos = new double**[nInputs];
   lohis = new double**[nInputs];
   hilos = new double**[nInputs];
   checkAllocate(hilos, "hilos in FF::analyze");
   for (ii = 0; ii < nInputs; ii++)
   {
      iEffects[ii] = new double*[nInputs];
      hihis[ii] = new double*[nInputs];
      lolos[ii] = new double*[nInputs];
      lohis[ii] = new double*[nInputs];
      hilos[ii] = new double*[nInputs];
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         hihis[ii][ii2] = new double[nReps+1];
         hilos[ii][ii2] = new double[nReps+1];
         lohis[ii][ii2] = new double[nReps+1];
         lolos[ii][ii2] = new double[nReps+1];
         iEffects[ii][ii2] = new double[nReps+2];
         for (rr = 0; rr <= nReps; rr++)
         {
            iEffects[ii][ii2][rr] = 0.0;
            hihis[ii][ii2][rr] = 0.0;
            lolos[ii][ii2][rr] = 0.0;
            hilos[ii][ii2][rr] = 0.0;
            lohis[ii][ii2][rr] = 0.0;
         }
      }
   }

   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, 
        "* Fractional Factorial (2-level) Interaction Analysis\n");
   if (adata.samplingMethod_ == PSUADE_SAMP_FF4 ||
       adata.samplingMethod_ == PSUADE_SAMP_RFF4)
   {
      printOutTS(PL_INFO, 
         "* Note: Since Fractional Factorial Resolution 4 is used,\n");
      printOutTS(PL_INFO, 
         "* the first and second order effects are confounded.\n");
   }
   printDashes(PL_INFO, 0);

   for (rr = 0; rr < nReps; rr++)
   {
      offset = rr * nSamples / nReps;
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
         {
   
            for (ss = 0; ss < nSamples/nReps; ss++)
            {
               txArray[ss] = X[(offset+ss)*nInputs+ii];
               twArray[ss] = X[(offset+ss)*nInputs+ii2];
               tyArray[ss] = Y[(offset+ss)*nOutputs+whichOutput];
            }


            sortDbleList3(nSamples/nReps, txArray, twArray, tyArray);


            for (ss = 0; ss < nSamples/nReps; ss+=nSamples/nReps/2)
               sortDbleList2(nSamples/nReps/2,&twArray[ss],&tyArray[ss]);

            accum = 0.0;
            for (ss = 0; ss < nSamples/nReps/4; ss++) accum += tyArray[ss];
            lolos[ii][ii2][rr+1] += accum;
            accum = 0.0;
            for (ss = nSamples*3/nReps/4; ss < nSamples/nReps; ss++) 
               accum += tyArray[ss];
            hihis[ii][ii2][rr+1] += accum;
            accum = 0.0;
            for (ss = nSamples/nReps/4; ss < nSamples/nReps/2; ss++) 
               accum += tyArray[ss];
            lohis[ii][ii2][rr+1] += accum;
            accum = 0.0;
            for (ss = nSamples/nReps/2; ss < nSamples*3/nReps/4; ss++)
               accum += tyArray[ss];
            hilos[ii][ii2][rr+1] += accum;
            iEffects[ii][ii2][rr+1] = 0.5 * (lolos[ii][ii2][rr+1] +
               hihis[ii][ii2][rr+1]-lohis[ii][ii2][rr+1]-
               hilos[ii][ii2][rr+1]);
         }
      }
      accum = 0.0;
      for (ii = 0; ii < nInputs; ii++)
         for (ii2 = ii+1; ii2 < nInputs; ii2++) 
            accum += PABS(iEffects[ii][ii2][rr+1]);
      if (accum == 0.0) accum = 1.0;
      for (ii = 0; ii < nInputs; ii++)
         for (ii2 = ii+1; ii2 < nInputs; ii2++) 
            iEffects[ii][ii2][rr+1] = iEffects[ii][ii2][rr+1]/accum;
   }
   for (ii = 0; ii < nInputs; ii++)
   {
      for (ii2 = ii+1; ii2 < nInputs; ii2++)
      {
         mean = 0.0;
         for (rr = 0; rr < nReps; rr++) mean += iEffects[ii][ii2][rr+1];
         mean = mean / (double) nReps;
         iEffects[ii][ii2][0] = mean;
         stdev = 0.0;
         if (nReps > 1)
         {
            for (rr = 0; rr < nReps; rr++)
               stdev += pow(iEffects[ii][ii2][rr+1]-mean, 2.0);
            stdev = sqrt(stdev / (double) (nReps - 1.0));
         }
         iEffects[ii][ii2][nReps+1] = stdev/sqrt(1.0*nReps);
         if (nReps == 1)
            printOutTS(PL_INFO,
                 "* Input %3d %3d =  %12.4e\n",ii+1,ii2+1, mean);
         else
            printOutTS(PL_INFO,
                 "* Input %3d %3d =  %12.4e (std err = %12.4e)\n",
                 ii+1,ii2+1, mean, stdev/sqrt(1.0*nReps));
      }
   }
   printAsterisks(PL_INFO, 0);

   if (fp != NULL)
   {
      if (psPlotTool_ == 1)
           fprintf(fp,"// matrix of pairwise effect: ind1,ind2,effect\n");
      else fprintf(fp,"%% matrix of pairwise effect: ind1,ind2,effect\n");
      fprintf(fp,"figure(2)\n");
      fprintf(fp,"clf\n");
      fprintf(fp,"nInputs = %d;\n", nInputs);
      fprintf(fp,"inputSwitches = [\n");
      if (nInputs <= 5)
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"1\n");
      else
      {
         for (ii = 0; ii < 5; ii++) fprintf(fp,"1\n");
         for (ii = 5; ii < nInputs; ii++) fprintf(fp,"0\n");
      }
      fprintf(fp,"];\n");
      fprintf(fp,"nPlots = sum(inputSwitches);\n");
      iNames = qData.strArray_;
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
      fprintf(fp, "A2 = [\n");
      for (ii = 0; ii < nInputs; ii++)
         for (ii2 = ii+1; ii2 < nInputs; ii2++)
            fprintf(fp, "%3d %3d %12.4e %12.4e\n",ii+1,ii2+1, 
                   iEffects[ii][ii2][0], iEffects[ii][ii2][nReps+1]);
      fprintf(fp, "];\n");
      fprintf(fp, "Ya = [];\n");
      fprintf(fp, "Yb = [];\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (ii2 = 0; ii2 < nInputs; ii2++)
         {
            fprintf(fp, "Y%d%da = [\n",ii+1,ii2+1);
            for (rr = 0; rr < nReps; rr++)
            {
               fprintf(fp, "   %24.16e\n", lolos[ii][ii2][rr+1]);
               fprintf(fp, "   %24.16e\n", hilos[ii][ii2][rr+1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "Y%d%db = [\n",ii+1,ii2+1);
            for (rr = 0; rr < nReps; rr++)
            {
               fprintf(fp, "   %24.16e\n", lohis[ii][ii2][rr+1]);
               fprintf(fp, "   %24.16e\n", hihis[ii][ii2][rr+1]);
            }
            fprintf(fp, "];\n");
            fprintf(fp, "Ya = [Ya Y%d%da];\n", ii+1,ii2+1);
            fprintf(fp, "Yb = [Yb Y%d%db];\n", ii+1,ii2+1);
         }
      }
      fprintf(fp,"cnt = 0;\n");
      fprintf(fp,"for ii = 1 : nInputs\n");
      fprintf(fp,"  if inputSwitches(ii) == 1\n");
      fprintf(fp,"    for jj = 1 : nInputs\n");
      fprintf(fp,"      if inputSwitches(jj) == 1 & ii ~= jj\n");
      fprintf(fp,"        cnt = cnt + 1;\n");
      fprintf(fp,"      end;\n");
      fprintf(fp,"      if inputSwitches(jj) == 1 & jj > ii\n");
      fprintf(fp,"        subplot(nPlots-1,nPlots-1,cnt)\n");
      fprintf(fp,"        Y1 = Ya(:,(ii-1)*nInputs+jj);\n");
      fprintf(fp,"        plot(Y1,'k')\n");
      if (psPlotTool_ == 1) 
           fprintf(fp,"        set(gca(),\"auto_clear\",\"off\")\n");
      else fprintf(fp,"        hold on\n");
      fprintf(fp,"        Y2 = Yb(:,(ii-1)*nInputs+jj);\n");
      fprintf(fp,"        plot(Y2,'b')\n");
      fprintf(fp,"        xlabel(strcat(Str(ii),'&',Str(jj)),'FontWeight',");
      fprintf(fp,"'bold','FontSize',12)\n");
      fprintf(fp,"        ylabel(Str(ii),'FontWeight','bold',");
      fprintf(fp,"'FontSize',12)\n");
      if (psPlotTool_ == 1) 
         fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
      else
      {
         fprintf(fp,"        hold off\n");
         fprintf(fp,"text(0.2,0.2,'black: P1 (lo to hi), P2 = lo','sc')\n");
         fprintf(fp,"text(0.2,0.3,'blue:  P1 (lo to hi), P2 = hi','sc')\n");
      }
      fwritePlotAxes(fp);
      fprintf(fp,"      end;\n");
      fprintf(fp,"    end;\n");
      fprintf(fp,"  end;\n");
      fprintf(fp,"end;\n");
      fprintf(fp,"subplot(nPlots-1,nPlots-1,1)\n");
      fwritePlotYLabel(fp,"Output Change");
      fprintf(fp,"set(gcf,'NextPlot','add');\n");
      fprintf(fp,"axes;\n");
      fprintf(fp,"h=title('Pairwise Interaction','fontSize',12,'fontWeight'");
      fprintf(fp,",'bold');\n");
      fprintf(fp,"set(gca,'Visible','off');\n");
      fprintf(fp,"set(h,'Visible','on');\n");
      fclose(fp);
      if (psPlotTool_ == 1)
         printOutTS(PL_INFO, 
              "The main effect plot has been generated in scilabmeff.sci.\n");
      else
         printOutTS(PL_INFO, 
              "The main effect plot has been generated in matlabmeff.m.\n");
   }
   printAsterisks(PL_INFO, 0);

   delete [] txArray;
   delete [] twArray;
   delete [] tyArray;
   for (ii = 0; ii < nInputs; ii++)
   {
      delete [] mEffects[ii];
      for (ii2 = 0; ii2 < nInputs; ii2++)
      {
         delete [] iEffects[ii][ii2];
         delete [] lolos[ii][ii2];
         delete [] lohis[ii][ii2];
         delete [] hilos[ii][ii2];
         delete [] hihis[ii][ii2];
      }
      delete [] iEffects[ii];
      delete [] lolos[ii];
      delete [] lohis[ii];
      delete [] hilos[ii];
      delete [] hihis[ii];
   }
   delete [] iArray;
   delete [] iEffects;
   delete [] lolos;
   delete [] lohis;
   delete [] hilos;
   delete [] hihis;
   delete [] mEffects;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
FFAnalyzer& FFAnalyzer::operator=(const FFAnalyzer &)
{
   printOutTS(PL_ERROR,"FFAnalysis operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int FFAnalyzer::get_nInputs()
{
   return nInputs_;
}

int FFAnalyzer::get_nSamples()
{
   return nSamples_;
}

int FFAnalyzer::get_nReps()
{
   return nReps_;
}

double *FFAnalyzer::get_means()
{
   double* retVal = NULL;
   if (means_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in FF::get_means");
      std::copy(means_, means_+nInputs_, retVal);
   }
   return retVal;
}
double *FFAnalyzer::get_stds()
{
   double* retVal = NULL;
   if (means_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in FF::get_stds");
      std::copy(stds_, stds_+nInputs_, retVal);
   }
   return retVal;
}
double **FFAnalyzer::get_mEffects()
{
   double** retVal = NULL;
   if (mEffects_)
   {
      retVal = new double*[nInputs_];
      checkAllocate(retVal, "retVal in FF::get_mEffects");
      for (int i=0; i<nInputs_; i++)
      {
    	  retVal[i] = new double[nInputs_];
    	  std::copy(mEffects_[i], mEffects_[i]+nReps_, retVal[i]);
      }
   }
   return retVal;
}

