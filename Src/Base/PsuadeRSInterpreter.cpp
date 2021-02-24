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
// Functions for the class PsuadeBase
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
//
// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <PsuadeCmakeConfig.h>

#ifdef WINDOWS
#include <windows.h>
#endif

// ------------------------------------------------------------------------
// local includes : class definition and utilities
// ------------------------------------------------------------------------
#include "Psuade.h"
#include "PsuadeBase.h"
#include "dtype.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PsuadeConfig.h"
#include "Vector.h"
#include "MainEffectAnalyzer.h"
#include "TwoParamAnalyzer.h"

// ------------------------------------------------------------------------
// local includes : function approximator and others
// ------------------------------------------------------------------------
#include "FuncApprox.h"
#include "pData.h"
#include "PDFManager.h"
#include "PDFNormal.h"

// ------------------------------------------------------------------------
// local defines 
// ------------------------------------------------------------------------
#define PABS(x)  ((x) > 0 ? x : -(x))

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
int PsuadeBase::RSBasedAnalysis(char *lineIn, PsuadeSession *session)
{
   int    nInputs, nOutputs, outputID, flag, faType;
   int    nSamples, *sampleStates, outputLevel, ss, jj, ii, kk, status;
   double *sampleInputs, *sampleOutputs, *iLowerB, *iUpperB;
   double ddata;
   char   command[1001], winput[1001], pString[1001], **inputNames;
   char   **outputNames, lineIn2[1001], dataFile[1001];
   FILE   *fp;
   PsuadeData *psuadeIO;
   pData      pPtr;

   // read in command and data
   sscanf(lineIn,"%s", command);
   nSamples = session->nSamples_;
   outputLevel = session->outputLevel_;
   nInputs = session->nInputs_;
   nOutputs = session->nOutputs_;
   sampleInputs = session->sampleInputs_;
   sampleOutputs = session->sampleOutputs_;
   sampleStates  = session->sampleStates_;
   psuadeIO = (PsuadeData *) session->psuadeIO_;
   inputNames = session->inputNames_;
   outputNames = session->outputNames_;
   iLowerB = session->inputLBounds_;
   iUpperB = session->inputUBounds_;
   
   // +++ rsadd 
   // ----------------------------------------------------------------
   if (!strcmp(command, "rsadd"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rsadd: add the outputs of two sample files.\n");
       printf("syntax: rsadd (no argument needed)\n");
       printf("This command merges 2 sample files by creating a new file\n");
       printf("with its outputs the sum of the outputs of the samples.\n");
       printf("The 2 samples do not need to have the same set of sample\n");
       printf("inputs. However, they need to have the same number of\n");
       printf("outputs. If they have different number of inputs, \n");
       printf("rs_index_file has to be provided inside the second sample\n");
       printf("to match. This command may be useful to combine the \n");
       printf("simulation sample and discrepancy sample produced by MCMC.\n");
       return 0;
     }

     int    nInps1, nOuts1, nInps2, nOuts2, nSamp1;
     char   filename1[1001], filename2[1001];
     PsuadeData *pd1=NULL, *pd2=NULL;
      
     printf("Please provide two sample files. Make sure the second sample\n");
     printf("has the target response surface type in it. If the second\n");
     printf("sample has different number of inputs from the first (primary)\n");
     printf("sample, make sure to use the rs_index_file feature to set\n");
     printf("things up properly.\n");
     printf("The final sample will have the same set of sample inputs as\n");
     printf("the primary sample but different output values.\n");
     printf("Enter the name of the primary sample (in PSUADE data format): ");
     scanf("%s", filename1);
     printf("Enter the name of the second sample (in PSUADE data format): ");
     scanf("%s", filename2);
     fgets(lineIn2, 500, stdin);
     pd1 = new PsuadeData();
     status = pd1->readPsuadeFile(filename1);
     if (status != 0)
     {
       printf("ERROR: cannot read sample file %s.\n", filename1);
       delete pd1;
       return -1;
     }
     pd2 = new PsuadeData();
     status = pd2->readPsuadeFile(filename2);
     if (status != 0)
     {
       printf("ERROR: cannot read sample file %s.\n", filename2);
       delete pd1;
       delete pd2;
       return -1;
     }
     pd1->getParameter("input_ninputs", pPtr);
     nInps1 = pPtr.intData_;
     pd2->getParameter("input_ninputs", pPtr);
     nInps2 = pPtr.intData_;
     if (nInps1 != nInps2)
     {
       printf("WARNING: The two sample files have different nInputs.\n");
       printf("         Make sure to use rs_index_file.\n"); 
     }
     pd1->getParameter("output_noutputs", pPtr);
     nOuts1 = pPtr.intData_;
     pd2->getParameter("output_noutputs", pPtr);
     nOuts2 = pPtr.intData_;
     if (nOuts1 != 1 || nOuts2 != 1)
     {
       printf("ERROR: the two sample files have nOutputs != 1.\n");
       delete pd1;
       delete pd2;
       return -1;
     }
     printf("** Creating response surface for second sample.\n");
     pd1->updateApplicationSection(filename2,NULL,NULL,NULL,NULL,-1);
     FunctionInterface *funcIO = createFunctionInterface(pd1);
     if (funcIO == NULL)
     {
       printf("ERROR: fail to create response surface.\n");
       delete pd1;
       delete pd2;
       return -1;
     }

     pd1->getParameter("method_nsamples", pPtr);
     nSamp1 = pPtr.intData_;
     pd1->getParameter("input_sample", pPtr);
     double *samInputs = pPtr.dbleArray_;
     pPtr.dbleArray_ = NULL;
     pd1->getParameter("output_sample", pPtr);
     double *samOuts1 = pPtr.dbleArray_;
     pPtr.dbleArray_ = NULL;
     pd1->getParameter("output_states", pPtr);
     int *samStates = pPtr.intArray_;
     pPtr.intArray_ = NULL;
     double *samOuts2 = new double[nSamp1];
     for (ss = 0; ss < nSamp1; ss++) 
        funcIO->evaluate(ss+1,nInps1,&samInputs[ss*nInps1],nOuts1,
                         &samOuts2[ss],0);
     for (ss = 0; ss < nSamp1; ss++) samOuts1[ss] += samOuts2[ss];
     pd1->updateOutputSection(nSamp1,nOuts1,samOuts1,samStates,NULL);
     strcpy(filename1, "psuade_rsadd_sample");
     pd1->writePsuadeFile(filename1,0);
     delete pd1;
     delete pd2;
     delete funcIO;
     delete [] samInputs;
     delete [] samOuts2;
     delete [] samStates;
     delete [] samOuts1;
     printf("New sample has been stored in psuade_rsadd_sample.\n");
   }

   // +++ rsua 
   // ----------------------------------------------------------------
   else if (!strcmp(command, "rsua"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rsua: uncertainty analysis on response surface\n");
       printf("syntax: rsua (no argument needed)\n");
       printf("This command perform uncertainty analysis on the response\n");
       printf("surface built from the LOADED sample. Uncertainty analysis\n");
       printf("is performed using a user-specified sample in PSUADE data\n");
       printf("format (created by running psuade on an input file). If you\n");
       printf("select a stochastic response surface type (Kriging, MARSB,\n");
       printf("or polynomial regression), the effect of response surface\n");
       printf("uncertainty will be shown on the PDF and CDF plots.\n");
       printf("NOTE: This analysis is intended to replace 'rs_ua'.\n");
       printf("      Turn on master mode to select between average case\n");
       printf("      and worst case analysis.\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }

     int   discFile=0, nInps, nOuts, dnInps;
     char  discFileName[1001], uaFileName[1001];
     PsuadeData *discIO=NULL, *sampleIO=NULL;
     FuncApprox **faPtrsRsEval=NULL;
      
     printAsterisks(PL_INFO, 0);
     printf("* Response surface-based Uncertainty Analysis\n");
     printDashes(PL_INFO, 0);
     printf("* To include response surface uncertainties, use stochastic\n");
     printf("* response surface such as polynomial regression, MARSB,\n");
     printf("* Kriging, ... (specified in your loaded data file).\n");
     printf("* This command computes worst case RS uncertainties. Turn\n");
     printf("* on master mode to select average case RS uncertainties.\n");
     printAsterisks(PL_INFO, 0);
     sscanf(lineIn,"%s %s", command, winput);
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;
     faPtrsRsEval = new FuncApprox*[2];
     faPtrsRsEval[0] = NULL;
     faPtrsRsEval[1] = NULL;
     printf("Use discrepancy model (in PSUADE data format)? ('y' or 'n') ");
     scanf("%s", winput);
     fgets(lineIn2, 500, stdin);
     if (winput[0] == 'y')
     {
       discFile = 1;
       printf("Enter discrepancy model file (in PSUADE data format): ");
       scanf("%s", discFileName);
       fgets(lineIn2, 500, stdin);
       discIO = new PsuadeData();
       status = discIO->readPsuadeFile(discFileName);
       if (status != 0)
       {
         printf("ERROR: cannot read discrepancy model file.\n");
         delete [] faPtrsRsEval;
         delete discIO;
         return -1;
       }
       discIO->getParameter("input_ninputs", pPtr);
       dnInps = pPtr.intData_;
       if (dnInps < nInputs)
       {
         printf("Discrepancy model has %d inputs. So the first\n", dnInps);
         printf("%d inputs in the model file will be assumed to\n", dnInps);
         printf("be associated with the inputs of the discrepancy model.\n");
       }
       discIO->getParameter("output_noutputs", pPtr);
       nOuts = pPtr.intData_;
       if (nOuts > 1)
       {
         printf("The discrepancy model has nOutputs > 1.\n");
         printf("This is currently not supported.\n");
         printf("Use 'odelete' to modify your discrepancy model file.\n");
         delete [] faPtrsRsEval;
         delete discIO;
         return -1;
       }
       printf("** CREATING RESPONSE SURFACE FOR DISCREPANCY MODEL\n");
       faPtrsRsEval[1] = genFAInteractive(discIO, 3);
       delete discIO;
       discIO = NULL;
     }

     printf("A sample is needed from you to propagate through the RS\n");
     printf("Enter UA sample file name (in PSUADE data format): ");
     scanf("%s", uaFileName);
     fgets(lineIn2, 500, stdin);
     sampleIO = new PsuadeData();
     status = sampleIO->readPsuadeFile(uaFileName);
     if (status != 0)
     {
       printf("ERROR: cannot read UA sample file.\n");
       delete [] faPtrsRsEval;
       delete sampleIO;
       return -1;
     }
     sampleIO->getParameter("input_ninputs", pPtr);
     nInps = pPtr.intData_;
     if (nInps != nInputs)
     {
       printf("ERROR: UA nInputs not equal to nInputs in local memory.\n");
       printf(":      input size in local memory = %d.\n",nInputs);
       printf(":      input size from UA file    = %d.\n",nInps);
       if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
       delete [] faPtrsRsEval;
       delete sampleIO;
       return -1;
     }
     sampleIO->getParameter("method_nsamples", pPtr);
     kk = pPtr.intData_;
     if (kk < 1000)
     {
       printf("ERROR: Your sample for UA should be at least 1000\n");
       printf("       to give any reasonable UA results.\n");
       if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
       delete [] faPtrsRsEval;
       delete sampleIO;
       return -1;
     }

     // perform UA
     int    userNSams, uaMethod=0;
     double *userSamInps, *userSamOuts, *userSamStds, *samOuts;
     sampleIO->getParameter("method_nsamples", pPtr);
     userNSams = pPtr.intData_;
     sampleIO->getParameter("input_sample", pPtr);
     userSamInps = pPtr.dbleArray_;
     pPtr.dbleArray_ = NULL;
     userSamOuts = new double[userNSams];
     userSamStds = new double[userNSams];
     if (psMasterMode_ == 1)
     {
       printf("The default is to perform the average case analysis (1): \n");
       printf("  - that is, for each sample point, generate a sub-sample\n");
       printf("    based on its uncertainty returned from the response\n"); 
       printf("    surface; and afterward, create output probability\n");
       printf("    distribution based on this enlarged sample.\n");
       printf("However, you can also perform a worst case analysis (2): \n");
       printf("  - that is, for each sample point, take the max and min\n");
       printf("    as its upper and lower bounds; and compute statistics\n");
       printf("    using all upper bounds and lower bounds, respectively.\n");
       sprintf(pString,"Enter 1 (average case) or 2 (worst case) analysis : ");
       uaMethod = getInt(1,2,pString);
       uaMethod--;
     }

     // stochastic RS with average case analysis
     if (uaMethod == 0)
     {
       printf("** CREATING RESPONSE SURFACE FOR PRIMARY MODEL\n");
       faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples);
       if (faPtrsRsEval[0] == NULL)
       {
          printf("ERROR: cannot generate response surface.\n");
          if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
          delete [] faPtrsRsEval;
          delete [] userSamInps;
          delete [] userSamOuts;
          delete [] userSamStds;
          delete sampleIO;
          return -1;
       }
       faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
       faPtrsRsEval[0]->setOutputLevel(0);
       samOuts = new double[nSamples];
       for (ss = 0; ss < nSamples; ss++)
          samOuts[ss] = sampleOutputs[ss*nOutputs+outputID];
       status = faPtrsRsEval[0]->initialize(sampleInputs,samOuts);
       delete [] samOuts;
       if (status != 0)
       {
          printf("ERROR: cannot initialize response surface.\n");
          if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
          delete [] faPtrsRsEval;
          delete [] userSamInps;
          delete [] userSamOuts;
          delete [] userSamStds;
          delete sampleIO;
          return -1;
       }
       faPtrsRsEval[0]->evaluatePointFuzzy(userNSams,userSamInps,
                                           userSamOuts,userSamStds);
       double discSamStd;
       if (discFile == 1)
       {
         for (ss = 0; ss < userNSams; ss++)
         {
           ddata = 
             faPtrsRsEval[1]->evaluatePointFuzzy(&userSamInps[ss*nInputs],
                                                  discSamStd);
           userSamOuts[ss] += ddata;
           ddata = pow(userSamStds[ss],2.0) + discSamStd * discSamStd;
           userSamStds[ss] = sqrt(ddata);
         }
       }

       fp = fopen("rsua_sample","w");
       if (fp != NULL)
       {
         fprintf(fp,"%% This file is primarily for diagnostics and \n");
         fprintf(fp,"%% expert analysis\n");
         fprintf(fp,"%% First line: nSamples nInputs\n");
         fprintf(fp,"%% All inputs, output, output-3*sigma, output+3*sigma\n");
         fprintf(fp,"%d %d 3\n", userNSams, nInputs);
         for (ss = 0; ss < userNSams; ss++)
         {
           for (ii = 0; ii < nInputs; ii++)
             fprintf(fp, "%e ", userSamInps[ss*nInputs+ii]);
           fprintf(fp, "%e ", userSamOuts[ss]);
           fprintf(fp, "%e ", userSamOuts[ss]-3*userSamStds[ss]);
           fprintf(fp, "%e\n", userSamOuts[ss]+3*userSamStds[ss]);
         }
         fclose(fp);
         printf("The outputs and stds of your sample has been written to\n");
         printf(" 'rsua_sample'.\n");
       }

       double mean=0, stdev=0;
       for (ss = 0; ss < userNSams; ss++) mean += userSamOuts[ss];
       mean /= (double) userNSams;
       for (ss = 0; ss < userNSams; ss++)
         stdev += pow(userSamOuts[ss]-mean, 2.0);
       stdev = sqrt(stdev/(double) userNSams);
       printAsterisks(PL_INFO, 0);
       printf("Sample mean    = %e (RS uncertainties not included)\n", mean);
       printf("Sample std dev = %e (RS uncertainties not included)\n", stdev);
       printEquals(PL_INFO, 0);

       int    nbins = 100, ntimes=20;
       int    **Fcounts = new int*[ntimes+1];
       double Fmax=-PSUADE_UNDEFINED;
       double Fmin=PSUADE_UNDEFINED;
       PDFNormal *rsPDF=NULL;
       for (ss = 0; ss < userNSams; ss++)
       {
         if (userSamOuts[ss]+3*userSamStds[ss] > Fmax)
           Fmax = userSamOuts[ss] + 3 * userSamStds[ss];
         if (userSamOuts[ss]-3*userSamStds[ss] < Fmin)
           Fmin = userSamOuts[ss] - 3 * userSamStds[ss];
       }
       Fmax = Fmax + 0.1 * (Fmax - Fmin);
       Fmin = Fmin - 0.1 * (Fmax - Fmin);
       if (Fmax == Fmin)
       {
         Fmax = Fmax + 0.1 * PABS(Fmax);
         Fmin = Fmin - 0.1 * PABS(Fmin);
       }
       for (ii = 0; ii <= ntimes; ii++)
       {
         Fcounts[ii] = new int[nbins];
         for (kk = 0; kk < nbins; kk++) Fcounts[ii][kk] = 0;
       }

       double *samOutTime = new double[ntimes*nInputs];
       double *samOutSave = new double[userNSams*ntimes];
       double d1, d2;
       for (ss = 0; ss < userNSams; ss++)
       {
          if (userSamStds[ss] == 0)
             for (ii = 0; ii < ntimes; ii++) samOutTime[ii] = userSamOuts[ss];
          else
          {
             rsPDF = new PDFNormal(userSamOuts[ss],userSamStds[ss]);
             d1 = userSamOuts[ss] - 3 * userSamStds[ss];
             d2 = userSamOuts[ss] + 3 * userSamStds[ss];
             rsPDF->genSample(ntimes,samOutTime,&d1,&d2);
             delete rsPDF;
          }
          for (ii = 0; ii < ntimes; ii++) 
             samOutSave[ss*ntimes+ii] = samOutTime[ii];

          ddata = userSamOuts[ss] - Fmin;
          if (Fmax > Fmin) ddata = ddata / ((Fmax - Fmin) / nbins);
          else             ddata = nbins / 2;
          kk = (int) ddata;
          if (kk < 0)      kk = 0;
          if (kk >= nbins) kk = nbins - 1;
          Fcounts[ntimes][kk]++;

          for (ii = 0; ii < ntimes; ii++)
          {
             ddata = samOutTime[ii] - Fmin;
             if (Fmax > Fmin)
                  ddata = ddata / ((Fmax - Fmin) / nbins);
             else ddata = nbins / 2;
             kk = (int) ddata;
             if (kk < 0)      kk = 0;
             if (kk >= nbins) kk = nbins - 1;
             Fcounts[ii][kk]++;
          }
       }
       delete [] samOutTime;
       double mean2=0, stdev2=0;
       for (ss = 0; ss < userNSams*ntimes; ss++) mean2 += samOutSave[ss];
       mean2 /= (double) (userNSams*ntimes);
       stdev2 = 0.0;
       for (ss = 0; ss < userNSams*ntimes; ss++)
          stdev2 += pow(samOutSave[ss] - mean2, 2.0);
       stdev2 = sqrt(stdev2/(double) (userNSams*ntimes));
       printf("Sample mean    = %e (RS uncertainties included)\n", mean2);
       printf("Sample std dev = %e (RS uncertainties included)\n", stdev2);
       printAsterisks(PL_INFO, 0);
       delete [] samOutSave;

       double dsum = 0.0;
       for (ss = 0; ss < userNSams; ss++) dsum += userSamStds[ss];
       fp = fopen("matlabrsua.m", "w");
       if (fp == NULL)
       {
          printf("INFO: cannot write the PDFs/CDFs to matlab file.\n");
       }
       else
       {
          fwriteHold(fp, 0);
          fprintf(fp, "subplot(2,2,1)\n");
          fprintf(fp, "X = [\n");
          for (kk = 0; kk < nbins; kk++)
             fprintf(fp, "%e\n",(Fmax-Fmin)/nbins*(0.5+kk)+Fmin);
          fprintf(fp, "];\n");
          for (ii = 0; ii <= ntimes; ii++)
          {
             fprintf(fp, "N%d = [\n", ii+1);
             for (kk = 0; kk < nbins; kk++)
                fprintf(fp, "%d\n",  Fcounts[ii][kk]);
             fprintf(fp, "];\n");
          }
          fprintf(fp, "N = [");
          for (ii = 0; ii <= ntimes; ii++)
             fprintf(fp, "N%d/sum(N%d) ", ii+1, ii+1);
          fprintf(fp, "];\n");
          fprintf(fp, "NA = N(:,%d+1);\n",ntimes);
          fprintf(fp, "NA = NA / sum(NA);\n");
          fprintf(fp, "NB = sum(N(:,1:%d)');\n",ntimes);
          fprintf(fp, "NB = NB' / sum(NB);\n");
          fprintf(fp, "NN = [NA NB];\n");
          fprintf(fp, "bar(X,NA,1.0)\n");
          fprintf(fp, "ymin = min(min(NA),min(NB));\n");
          fprintf(fp, "ymax = max(max(NA),max(NB));\n");
          fprintf(fp, "axis([min(X) max(X) ymin ymax])\n");
          fwritePlotAxes(fp);
          fwritePlotTitle(fp, "Prob. Dist. (means of RS)");
          fwritePlotXLabel(fp, "Output Value");
          fwritePlotYLabel(fp, "Probabilities)");
          fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",
                  mean);
          fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                  stdev);
          if (dsum == 0.0)
          {
             printf("Deterministic RS used ==> no RS uncertainties.\n");
             fprintf(fp,"subplot(2,2,3)\n");
             fwritePlotAxes(fp);
             fwritePlotTitle(fp,"Prob. Dist. (RS with uncertainties)");
             fprintf(fp,"text(0.01,0.5,'No RS uncertainty info.','sc',");
             fprintf(fp,"'FontSize',11)\n");
             fprintf(fp,"text(0.01,0.4,'Deterministic RS used.','sc',");
             fprintf(fp,"'FontSize',11)\n");
          }
          else
          {
             fprintf(fp,"subplot(2,2,3)\n");
             fprintf(fp,"bar(X,NB,1.0)\n");
             fprintf(fp,"axis([min(X) max(X) ymin ymax])\n");
             fwritePlotAxes(fp);
             fwritePlotTitle(fp,"Prob. Dist. (RS with uncertainties)");
             fwritePlotXLabel(fp,"Output Value");
             fwritePlotYLabel(fp,"Probabilities");
             fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",
                  mean2);
             fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                  stdev2);
          }
          for (ii = 0; ii <= ntimes; ii++)
          {
             fprintf(fp,"for ii = 2 : %d\n", nbins);
             fprintf(fp,"  N%d(ii) = N%d(ii) + N%d(ii-1);\n",ii+1,ii+1,ii+1);
             fprintf(fp,"end;\n");
          }
          fprintf(fp, "N = [");
          for (ii = 0; ii <= ntimes; ii++)
             fprintf(fp,"N%d/N%d(%d) ", ii+1, ii+1, nbins);
          fprintf(fp, "];\n");
          fprintf(fp, "subplot(2,2,[2 4])\n");
          fprintf(fp, "NA = N(:,%d+1);\n",ntimes);
          fprintf(fp, "NA = NA / NA(%d);\n",nbins);
          fprintf(fp, "NB = sum(N(:,1:%d)');\n",ntimes);
          fprintf(fp, "NB = NB' / NB(%d);\n", nbins);
          fprintf(fp, "NN = [NA NB];\n");
          if (dsum == 0.0)
          {
             fprintf(fp, "plot(X,NA,'linewidth',3)\n");
             fwritePlotTitle(fp,"Cum. Dist.: (b) mean; (g) with uncertainties");
          }
          else
          {
             fprintf(fp, "plot(X,NN,'linewidth',3)\n");
             fwritePlotTitle(fp,"Cum. Dist.: (*) uncertainties unavailable");
          }
          fwritePlotAxes(fp);
          fwritePlotXLabel(fp, "Output Value");
          fwritePlotYLabel(fp, "Probabilities");
          fclose(fp);
          printf("Output distribution plots are in matlabrsua.m.\n");
       }
       for (ii = 0; ii <= ntimes; ii++) delete [] Fcounts[ii];
       delete [] Fcounts;
     }

     // stochastic RS with worst case analysis
     else if (uaMethod == 1)
     {
       printf("** CREATING RESPONSE SURFACE FOR PRIMARY MODEL\n");
       faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples);
       if (faPtrsRsEval[0] == NULL)
       {
         printf("ERROR: cannot generate response surface.\n");
         if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
         delete [] faPtrsRsEval;
         delete [] userSamInps;
         delete [] userSamOuts;
         delete [] userSamStds;
         delete sampleIO;
         return -1;
       }
       faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
       faPtrsRsEval[0]->setOutputLevel(0);
       samOuts = new double[nSamples];
       for (ss = 0; ss < nSamples; ss++)
         samOuts[ss] = sampleOutputs[ss*nOutputs+outputID];
       status = faPtrsRsEval[0]->initialize(sampleInputs,samOuts);
       delete [] samOuts;
       if (status != 0)
       {
         printf("ERROR: cannot initialize response surface.\n");
         if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
         delete [] faPtrsRsEval;
         delete [] userSamInps;
         delete [] userSamOuts;
         delete [] userSamStds;
         delete sampleIO;
         return -1;
       }
       
       faPtrsRsEval[0]->evaluatePointFuzzy(userNSams,userSamInps,
                                           userSamOuts,userSamStds);
       
       double discSamStd;
       if (discFile == 1)
       {
         for (ss = 0; ss < userNSams; ss++)
         {
           userSamOuts[ss] += 
              faPtrsRsEval[1]->evaluatePointFuzzy(&userSamInps[ss*nInputs],
                                                  discSamStd);
           ddata = pow(userSamStds[ss],2.0) + discSamStd * discSamStd;
           userSamStds[ss] = sqrt(ddata);
         }
       }
       
       double mean=0, stdev=0;
       for (ss = 0; ss < userNSams; ss++) mean += userSamOuts[ss];
       mean /= (double) userNSams;
       for (ss = 0; ss < userNSams; ss++)
          stdev += pow(userSamOuts[ss]-mean, 2.0);
       stdev = sqrt(stdev/(double) userNSams);
       printAsterisks(PL_INFO, 0);
       printf("Sample mean    = %e (RS uncertainties not included)\n",mean);
       printf("Sample std dev = %e (RS uncertainties not included)\n",stdev);
       printEquals(PL_INFO, 0);

       fp = fopen("rsua_sample","w");
       fprintf(fp,"%% This file is primarily for diagnostics and \n");
       fprintf(fp,"%% expert analysis\n");
       fprintf(fp,"%% First line: nSamples nInputs\n");
       fprintf(fp,"%% All inputs, output, output-3*sigma, output+3*sigma\n");
       fprintf(fp,"%d %d 3\n", userNSams, nInputs);
       for (ss = 0; ss < userNSams; ss++)
       {
         for (ii = 0; ii < nInputs; ii++)
           fprintf(fp, "%e ", userSamInps[ss*nInputs+ii]);
         fprintf(fp,"%e ", userSamOuts[ss]);
         fprintf(fp,"%e ", userSamOuts[ss]-3*userSamStds[ss]);
         fprintf(fp,"%e\n", userSamOuts[ss]+3*userSamStds[ss]);
       }
       fclose(fp);
       
       int    nbins = 100, ntimes=7;
       int    **Fcounts = new int*[ntimes+1];
       double Fmax=-PSUADE_UNDEFINED;
       double Fmin=PSUADE_UNDEFINED;
       PDFNormal *rsPDF=NULL;
       for (ss = 0; ss < userNSams; ss++)
       {
         if (userSamOuts[ss]+3*userSamStds[ss] > Fmax)
            Fmax = userSamOuts[ss] + 3 * userSamStds[ss];
         if (userSamOuts[ss]-3*userSamStds[ss] < Fmin)
            Fmin = userSamOuts[ss] - 3 * userSamStds[ss];
       }
       Fmax = Fmax + 0.1 * (Fmax - Fmin);
       Fmin = Fmin - 0.1 * (Fmax - Fmin);
       if (Fmax == Fmin)
       {
         Fmax = Fmax + 0.1 * PABS(Fmax);
         Fmin = Fmin - 0.1 * PABS(Fmin);
       }
       for (ii = 0; ii <= ntimes; ii++)
       {
         Fcounts[ii] = new int[nbins];
         for (kk = 0; kk < nbins; kk++) Fcounts[ii][kk] = 0;
       }

       double dsum = 0.0;
       for (ss = 0; ss < userNSams; ss++)
       {
         for (ii = 0; ii < ntimes; ii++)
         {
           ddata = userSamOuts[ss]+userSamStds[ss]*(ii-3) - Fmin;
           if (Fmax > Fmin)
                ddata = ddata / ((Fmax - Fmin) / nbins);
           else ddata = nbins / 2;
           kk = (int) ddata;
           if (kk < 0)      kk = 0;
           if (kk >= nbins) kk = nbins - 1;
           Fcounts[ii][kk]++;
         }
         dsum += userSamStds[ss];
       }

       fp = fopen("matlabrsua.m", "w");
       if (fp == NULL)
       {
         printf("INFO: cannot write the PDFs/CDFs to matlab file.\n");
       }
       else
       {
         fwriteHold(fp, 0);
         fprintf(fp, "%% worst case analysis\n");
         fprintf(fp, "X = [\n");
         for (kk = 0; kk < nbins; kk++)
           fprintf(fp, "%e\n", (Fmax-Fmin)/nbins*(0.5+kk)+Fmin);
         fprintf(fp, "];\n");
         for (ii = 0; ii < ntimes; ii++)
         {
           fprintf(fp, "E%d = [\n", ii+1);
           for (kk = 0; kk < nbins; kk++) fprintf(fp, "%d\n",  Fcounts[ii][kk]);
           fprintf(fp, "];\n");
         }
         fprintf(fp, "EE = [");
         for (ii = 0; ii < ntimes; ii++)
           fprintf(fp, "E%d/sum(E%d) ", ii+1, ii+1);
         fprintf(fp, "];\n");
         fprintf(fp, "subplot(2,2,1)\n");
         fprintf(fp, "bar(X,EE(:,4),1.0)\n");
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Prob. Dist. (means of RS)");
         fwritePlotXLabel(fp, "Output Value");
         fwritePlotYLabel(fp, "Probabilities)");
         fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",
                 mean);
         fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                 stdev);
         fprintf(fp, "axis([min(X) max(X) min(min(EE)) max(max(EE)))\n");
         if (dsum == 0.0)
         {
           printf("Deterministic RS used ==> no RS uncertainties.\n");
           fprintf(fp,"subplot(2,2,3)\n");
           fwritePlotAxes(fp);
           fwritePlotTitle(fp,"Prob. Dist. (RS with uncertainties)");
           fprintf(fp,"text(0.01,0.5,'No RS uncertainty info.','sc',");
           fprintf(fp,"'FontSize',11)\n");
           fprintf(fp,"text(0.01,0.4,'Deterministic RS used.','sc',");
           fprintf(fp,"'FontSize',11)\n");
         }
         else
         {
           fprintf(fp, "subplot(2,2,3)\n");
           fprintf(fp, "plot(X,EE,'lineWidth',2)\n");
           fprintf(fp, "axis([min(X) max(X) min(min(EE)) max(max(EE)))\n");
           fwritePlotAxes(fp);
           fwritePlotTitle(fp, "Prob. Dist. (-3,2,1,0,1,2,3 std.)");
           fwritePlotXLabel(fp, "Output Value");
           fwritePlotYLabel(fp, "Probabilities");
         }
         fprintf(fp, "subplot(2,2,[2 4])\n");
         for (ii = 0; ii < ntimes; ii++)
         {
           fprintf(fp, "for ii = 2 : %d\n", nbins);
           fprintf(fp, "   E%d(ii) = E%d(ii) + E%d(ii-1);\n",ii+1,ii+1,ii+1);
           fprintf(fp, "end;\n");
         }
         fprintf(fp, "EE = [");
         for (ii = 0; ii < ntimes; ii++)
           fprintf(fp, "E%d/E%d(%d) ", ii+1, ii+1, nbins);
         fprintf(fp, "];\n");
         fprintf(fp, "plot(X,EE,'linewidth',2)\n");
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Cum. Dist. (-3,2,1,0,1,2,3 std.)");
         fwritePlotXLabel(fp, "Output Value");
         fwritePlotYLabel(fp, "Probabilities");
         fclose(fp);
         printf("Output distribution plots has been created in matlabrsua.m\n");
         for (ii = 0; ii < ntimes; ii++) delete [] Fcounts[ii];
         delete [] Fcounts;
       }
     }
     if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
     if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
     delete [] faPtrsRsEval;
     if (sampleIO != NULL) delete sampleIO;
     delete [] userSamInps;
     delete [] userSamOuts;
     delete [] userSamStds;
   }

   // +++ rsuab 
   else if (!strcmp(command, "rsuab"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rsuab: uncertainty analysis on bootstrapped response surfaces\n");
       printf("syntax: rsuab (no argument needed)\n");
       printf("This command perform uncertainty analysis on the response\n");
       printf("surface built from the loaded sample. The effect of response\n");
       printf("surface uncertainties will be created by bootstrapping.\n");
       printf("NOTE: Turn on master mode to finetune bootstrap selections.\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }

     int    discFile=0, nInps, nOuts, dnInps, numBS;
     double bsPC=0;
     char   discFileName[1001], uaFileName[1001];
     PsuadeData *discIO=NULL, *sampleIO=NULL;
     FuncApprox **faPtrsRsEval=NULL;
     
     sscanf(lineIn,"%s %s", command, winput);
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;

     faPtrsRsEval = new FuncApprox*[2];
     faPtrsRsEval[0] = NULL;
     faPtrsRsEval[1] = NULL;
     printf("You have the option to add a correction (discrepancy) to\n");
     printf("your sample. The discrepancy model will be constructed\n");
     printf("from another sample (in PSUADE data format).\n");
     printf("Use discrepancy model? ('y' or 'n') ");
     scanf("%s", winput);
     fgets(lineIn2, 500, stdin);
     if (winput[0] == 'y')
     {
       discFile = 1;
       printf("Enter discrepancy model file (in PSUADE data format): ");
       scanf("%s", discFileName);
       fgets(lineIn2, 500, stdin);
       discIO = new PsuadeData();
       status = discIO->readPsuadeFile(discFileName);
       if (status != 0)
       {
         printf("ERROR: cannot read discrepancy model file.\n");
         delete [] faPtrsRsEval;
         delete discIO;
         return -1;
       }
       discIO->getParameter("input_ninputs", pPtr);
       dnInps = pPtr.intData_;
       if (dnInps < nInputs)
       {
         printf("Discrepancy model has %d inputs. So the first\n", dnInps);
         printf("%d inputs in the model file will be assumed to\n", dnInps);
         printf("be associated with the inputs of the discrepancy model.\n");
       }
       discIO->getParameter("output_noutputs", pPtr);
       nOuts = pPtr.intData_;
       if (nOuts > 1)
       {
         printf("The discrepancy model has nOutputs > 1.\n");
         printf("This is currently not supported.\n");
         printf("Use 'odelete' to modify your discrepancy model file.\n");
         delete [] faPtrsRsEval;
         delete discIO;
         return -1;
       }
       printf("** CREATING RESPONSE SURFACE FOR DISCREPANCY MODEL\n");
       faPtrsRsEval[1] = genFAInteractive(discIO, 3);
       delete discIO;
       discIO = NULL;
     }

     printf("A sample is needed to be propagated through the response\n");
     printf("surface. This sample can be created via the gensample\n");
     printf("command.\n");
     printf("Enter UA sample file name (in PSUADE data format): ");
     scanf("%s", uaFileName);
     fgets(lineIn2, 500, stdin);
     sampleIO = new PsuadeData();
     status = sampleIO->readPsuadeFile(uaFileName);
     if (status != 0)
     {
       printf("ERROR: cannot read UA sample file.\n");
       delete [] faPtrsRsEval;
       delete sampleIO;
       return -1;
     }
     sampleIO->getParameter("input_ninputs", pPtr);
     nInps = pPtr.intData_;
     if (nInps != nInputs)
     {
       printf("ERROR: UA nInputs not equal to nInputs in local memory.\n");
       printf(":      input size in local memory = %d.\n",nInputs);
       printf(":      input size from UA file    = %d.\n",nInps);
       if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
       delete [] faPtrsRsEval;
       delete sampleIO;
       return -1;
     }
     sampleIO->getParameter("method_nsamples", pPtr);
     kk = pPtr.intData_;
     if (kk < 1000)
     {
       printf("ERROR: Your sample for UA should be at least 1000\n");
       printf("       to give any reasonable UA results.\n");
       if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
       delete [] faPtrsRsEval;
       delete sampleIO;
       return -1;
     }
     sprintf(pString, "How many bootstrapped samples to use (10 - 300) : ");
     numBS = getInt(10, 300, pString);
     if (psMasterMode_ == 1)
     {
       printf("Bootstrapped samples will be created from randomly drawing\n");
       printf("from your response surface sample. Normally a random draw\n");
       printf("may include around 60%% of the original sample. You can\n");
       printf("increase this percentage by entering your desired percentage.\n");
       kk = 0;
       while ((bsPC < 60 || bsPC > 90) && kk < 10)
       {
          printf("Enter percentage (60-90) : ");
          scanf("%lg", &bsPC);
          kk++;
       }
       if (bsPC < 60 || bsPC > 90 || kk >= 10) bsPC = 0; 
       else                                    bsPC *= 0.01;
     }

     int    userNSams;
     double *userSamInps, *userSamOuts, *userSamStds;
     sampleIO->getParameter("method_nsamples", pPtr);
     userNSams = pPtr.intData_;
     sampleIO->getParameter("input_sample", pPtr);
     userSamInps = pPtr.dbleArray_;
     pPtr.dbleArray_ = NULL;
     userSamOuts = new double[userNSams];
     userSamStds = new double[userNSams];

     // bootstrapped method
     int bsnSams, rsMethod;
     printf("** CREATING RESPONSE SURFACE FOR PRIMARY MODEL\n");
     faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples);
     if (faPtrsRsEval[0] == NULL)
     {
       printf("ERROR: cannot generate primary response surface.\n");
       if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
       delete [] faPtrsRsEval;
       delete sampleIO;
       return -1;
     }
     rsMethod = faPtrsRsEval[0]->getID(); 
     delete faPtrsRsEval[0];
     faPtrsRsEval[0] = NULL;
     
     int    its, *useFlags = new int[nSamples];
     double *bsSamInps = new double[nSamples*nInputs];
     double *bsSamOuts = new double[nSamples];
     double *bsMeans = new double[numBS];
     double *bsStds  = new double[numBS];

     fp = NULL;
     fp = fopen("matlabrsuab.m", "w");
     for (its = 0; its < numBS; its++)
     {
       for (ss = 0; ss < nSamples; ss++) useFlags[ss] = 0;
       bsnSams = 0;
       kk = 0;
       while (kk < nSamples || (1.0*bsnSams/nSamples) < bsPC)
       {
         jj = PSUADE_rand() % nSamples;
         if (useFlags[jj] == 0)
         {
           for (ii = 0; ii < nInputs; ii++)
              bsSamInps[bsnSams*nInputs+ii] = sampleInputs[jj*nInputs+ii];
           bsSamOuts[bsnSams] = sampleOutputs[jj*nOutputs+outputID];
           useFlags[jj] = 1;
           bsnSams++;
         }
         kk++;
       }
       printf("Bootstrap %d has sample size = %d (%d)\n",its+1,bsnSams,nSamples);
       faPtrsRsEval[0] = genFA(rsMethod, nInputs, -1, bsnSams);
       faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
       faPtrsRsEval[0]->setOutputLevel(0);
       status = faPtrsRsEval[0]->initialize(bsSamInps,bsSamOuts);
       if (status != 0)
       {
         printf("ERROR: in initializing response surface (1).\n");
         if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
         if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
         delete [] faPtrsRsEval;
         delete [] bsSamInps;
         delete [] bsSamOuts;
         delete [] userSamInps;
         delete [] userSamOuts;
         delete [] userSamStds;
         delete [] useFlags;
         delete [] bsMeans;
         delete [] bsStds;
         delete sampleIO;
         return -1;
       } 
       faPtrsRsEval[0]->evaluatePoint(userNSams,userSamInps,userSamOuts);
       delete faPtrsRsEval[0];
       faPtrsRsEval[0] = NULL;
       if (discFile == 1)
       {
         for (ss = 0; ss < userNSams; ss++)
         {
           ddata = faPtrsRsEval[1]->evaluatePoint(&userSamInps[ss*nInputs]);
           userSamOuts[ss] += ddata;
         }
       }
       bsMeans[its] = bsStds[its] = 0.0;
       for (ss = 0; ss < userNSams; ss++) bsMeans[its] += userSamOuts[ss];
       bsMeans[its] /= (double) userNSams;
       for (ss = 0; ss < userNSams; ss++)
         bsStds[its] += pow(userSamOuts[ss] - bsMeans[its], 2.0);
       bsStds[its] = sqrt(bsStds[its] / userNSams);
       if (fp != NULL)
       {
         fprintf(fp, "Y = [\n");
         for (ss = 0; ss < userNSams; ss++) fprintf(fp,"%e\n",userSamOuts[ss]);
         fprintf(fp, "];\n");
         fprintf(fp, "Y%d = sort(Y);\n",its+1);
         fprintf(fp, "X%d = (1 : %d)';\n", its+1, userNSams);
         fprintf(fp, "X%d = X%d / %d;\n", its+1, its+1, userNSams);
         if (its == 0)
         {
           fprintf(fp, "YY = Y%d;\n", its+1);
           fprintf(fp, "XX = X%d;\n", its+1);
         }
         else
         {
           fprintf(fp, "YY = [YY Y%d];\n", its+1);
           fprintf(fp, "XX = [XX X%d];\n", its+1);
         }
       }
     }
     faPtrsRsEval[0] = genFA(rsMethod, nInputs, -1, nSamples);
     faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
     faPtrsRsEval[0]->setOutputLevel(0);
     for (ss = 0; ss < nSamples; ss++)
       bsSamOuts[ss] = sampleOutputs[ss*nOutputs+outputID];
     status = faPtrsRsEval[0]->initialize(sampleInputs,bsSamOuts);
     faPtrsRsEval[0]->evaluatePoint(userNSams,userSamInps,userSamOuts);
     delete faPtrsRsEval[0];
     faPtrsRsEval[0] = NULL;
     if (discFile == 1)
     {
       for (ss = 0; ss < userNSams; ss++)
       {
         ddata = faPtrsRsEval[1]->evaluatePoint(&userSamInps[ss*nInputs]);
         userSamOuts[ss] += ddata;
       }
     }
     printAsterisks(PL_INFO, 0);
     double mean, stdev;
     mean = stdev = 0.0;
     for (its = 0; its < numBS; its++) mean += bsMeans[its];
     mean /= (double) numBS;
     for (ss = 0; ss < numBS; ss++) stdev += pow(bsMeans[ss]-mean, 2.0);
     stdev = sqrt(stdev/(double) numBS);
     printf("Sample mean    = %e (std = %e)\n", mean, stdev);
     mean = stdev = 0.0;
     for (its = 0; its < numBS; its++) mean += bsStds[its];
     mean /= (double) numBS;
     for (ss = 0; ss < numBS; ss++) stdev += pow(bsStds[ss]-mean, 2.0);
     stdev = sqrt(stdev/(double) numBS);
     printf("Sample std dev = %e (std = %e)\n", mean, stdev);
     printEquals(PL_INFO, 0);
     if (fp != NULL)
     {
       fprintf(fp, "Y0 = [\n");
       for (ss = 0; ss < userNSams; ss++) fprintf(fp,"%e\n",userSamOuts[ss]);
       fprintf(fp, "];\n");
       fwriteHold(fp, 0);
       fprintf(fp,"subplot(2,2,3);\n");
       fprintf(fp,"[nk,xk] = hist(YY,50);\n");
       fprintf(fp,"nk = nk / %d;\n",userNSams);
       fprintf(fp,"plot(xk,nk, 'lineWidth',2)\n");
       fwritePlotAxes(fp);
       fwritePlotTitle(fp, "Prob. Dist. (RS with uncertainties)");
       fwritePlotXLabel(fp,"Output Value");
       fwritePlotYLabel(fp,"Probabilities");
       fprintf(fp,"subplot(2,2,1);\n");
       fprintf(fp,"[nk0,xk0] = hist(Y0,xk,50);\n");
       fprintf(fp,"nk0 = nk0 / %d;\n",userNSams);
       fprintf(fp,"bar(xk0,nk0,1.0)\n");
       fwritePlotAxes(fp);
       fwritePlotTitle(fp, "Prob. Dist. (Original Sample)");
       fwritePlotXLabel(fp,"Output Value");
       fwritePlotYLabel(fp,"Probabilities");
       mean = stdev = 0.0;
       for (ss = 0; ss < userNSams; ss++) mean += userSamOuts[ss];
       mean /= (double) userNSams;
       for (ss = 0; ss < userNSams; ss++)
          stdev += pow(userSamOuts[ss] - mean, 2.0);
       stdev = sqrt(stdev / userNSams);
       fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",mean);
       fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",stdev);
       fprintf(fp,"subplot(2,2,[2 4]);\n");
       fprintf(fp,"plot(YY, XX, 'lineWidth',3)\n");
       fwritePlotAxes(fp);
       fwritePlotTitle(fp, "Cumulative Distributions");
       fwritePlotXLabel(fp, "Output Value");
       fwritePlotYLabel(fp, "Probabilities");
       fclose(fp);
       printf("Output distribution plots has been created in matlabrsuab.m.\n");
     }
     delete [] bsSamInps;
     delete [] bsSamOuts;
     delete [] useFlags;
     delete [] bsMeans;
     delete [] bsStds;
     if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
     if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
     delete [] faPtrsRsEval;
     if (sampleIO != NULL) delete sampleIO;
     delete [] userSamInps;
     delete [] userSamOuts;
     delete [] userSamStds;
   }

   // +++ rssobol1 
   else if (!strcmp(command, "rssobol1"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rssobol1: RS-based Sobol' first order indices\n");
       printf("syntax: rssobol1 (no argument needed)\n");
       printf("note: to compute error bars for indices, turn on\n");
       printf("      ana_expert mode and set ntimes to 100.\n");
       return -1;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;

     int   analysisMethod = PSUADE_ANA_RSSOBOL1;
     AnalysisManager *anaManager = new AnalysisManager();
     anaManager->setup(analysisMethod, 0);
     psuadeIO->getParameter("ana_diagnostics",pPtr);
     int diagKeep = pPtr.intData_;
     psuadeIO->updateAnalysisSection(-1,-1,-1,outputLevel,-1,-1);
     anaManager->analyze(psuadeIO, 0, NULL, outputID);
     psuadeIO->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
     pData *pdata = psuadeIO->getAuxData();
     if (pdata->nDbles_ >= nInputs)
     {
       printEquals(PL_INFO, 0);
       printf("Main Effect Statistics: \n");
       if (pdata->dbleData_ > 0)
       {
         for (ii = 0; ii < nInputs; ii++)
         {
           printf("Input %4d: Sobol' main effect = %12.4e",
                  ii+1,pdata->dbleArray_[ii]);
           if (pdata->nDbles_ == 3*nInputs)
                printf(", bounds = [%12.4e, %12.4e]\n",
                       pdata->dbleArray_[ii+nInputs],
                       pdata->dbleArray_[ii+2*nInputs]);
           else printf(" (unnormalized)\n");
         } 
         if (psPlotTool_ == 1)
              fp = fopen("scilabrssobol1.sci", "w");
         else fp = fopen("matlabrssobol1.m", "w");
         if (fp == NULL) 
           printf("rssobol1 ERROR: cannot open file to save data\n");
         else
         {
           if (psPlotTool_ == 1)
           {
             fprintf(fp,"// This file contains Sobol' indices\n");
             fprintf(fp,"// set sortFlag = 1 and set nn to be the number\n");
             fprintf(fp,"// of inputs to display.\n");
           }
           else
           {
             fprintf(fp,"%% This file contains Sobol' indices\n");
             fprintf(fp,"%% set sortFlag = 1 and set nn to be the number\n");
             fprintf(fp,"%% of inputs to display.\n");
           }
           fprintf(fp, "sortFlag = 0;\n");
           fprintf(fp, "nn = %d;\n", nInputs);
           fprintf(fp, "Mids = [\n");
           for (ii = 0; ii < nInputs; ii++) 
             fprintf(fp,"%24.16e\n", pdata->dbleArray_[ii]/pdata->dbleData_);
           fprintf(fp, "];\n");
           if (pdata->nDbles_ == 3*nInputs)
           {
             fprintf(fp, "Mins = [\n");
             for (ii = 0; ii < nInputs; ii++) 
               fprintf(fp,"%24.16e\n", 
                        pdata->dbleArray_[nInputs+ii]/pdata->dbleData_);
             fprintf(fp, "];\n");
             fprintf(fp, "Maxs = [\n");
             for (ii = 0; ii < nInputs; ii++) 
               fprintf(fp,"%24.16e\n", 
                        pdata->dbleArray_[2*nInputs+ii]/pdata->dbleData_);
             fprintf(fp, "];\n");
           }
           if (inputNames == NULL)
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
               if (inputNames[ii] != NULL) 
                    fprintf(fp,"'%s',",inputNames[ii]);
               else fprintf(fp,"'X%d',",ii+1);
             }
             if (inputNames[nInputs-1] != NULL)
                  fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
             else fprintf(fp,"'X%d'};\n",nInputs);
           }
           fwriteHold(fp, 0);
           fprintf(fp, "if (sortFlag == 1)\n");
           if (psPlotTool_ == 1)
                fprintf(fp, "  [Mids, I2] = gsort(Mids);\n");
           else fprintf(fp, "  [Mids, I2] = sort(Mids,'descend');\n");
           if (pdata->nDbles_ == 3*nInputs)
           {
             fprintf(fp, "  Maxs = Maxs(I2);\n");
             fprintf(fp, "  Mins = Mins(I2);\n");
           }
           fprintf(fp, "  Str  = Str(I2);\n");
           fprintf(fp, "  I2 = I2(1:nn);\n");
           fprintf(fp, "  Mids = Mids(1:nn);\n");
           if (pdata->nDbles_ == 3*nInputs)
           {
             fprintf(fp, "  Maxs = Maxs(1:nn);\n");
             fprintf(fp, "  Mins = Mins(1:nn);\n");
           }
           fprintf(fp, "  Str  = Str(1:nn);\n");
           fprintf(fp, "end\n");
           if (pdata->nDbles_ == 3*nInputs)
           {
             fprintf(fp, "ymin = min(Mins);\n");
             fprintf(fp, "ymax = max(Maxs);\n");
           }
           else
           {
             fprintf(fp, "ymin = min(Mids);\n");
             fprintf(fp, "ymax = max(Mids);\n");
           }
           fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
           if (psPlotTool_ == 1) fprintf(fp, "drawlater\n");
           fprintf(fp, "bar(Mids,0.8);\n");
           if (pdata->nDbles_ == 3*nInputs)
           {
             fprintf(fp,"for ii = 1:nn\n");
             if (psPlotTool_ == 1)
                fprintf(fp,
                   "// h = plot(ii,Means(ii),'r*','MarkerSize',13);\n");
             else 
                fprintf(fp,
                   "%% h = plot(ii,Means(ii),'r*','MarkerSize',13);\n");
             fprintf(fp,"   if (ii == 1)\n");
             fwriteHold(fp, 1);
             fprintf(fp,"   end;\n");
             fprintf(fp,"   XX = [ii ii];\n");
             fprintf(fp,"   YY = [Mins(ii) Maxs(ii)];\n");
             fprintf(fp,
                "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
             fprintf(fp,"'k','MarkerFaceColor','g','MarkerSize',13)\n");
             fprintf(fp,"end;\n");
           }
           fwritePlotAxes(fp);
           fprintf(fp,"ymin=0;\n");
           if (psPlotTool_ == 1)
           {
             fprintf(fp, "a=gca();\n");
             fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
             fprintf(fp, "newtick = a.x_ticks;\n");
             fprintf(fp, "newtick(2) = [1:nn]';\n");
             fprintf(fp, "newtick(3) = Str';\n");
             fprintf(fp, "a.x_ticks = newtick;\n");
             fprintf(fp, "a.x_label.font_size = 3;\n");
             fprintf(fp, "a.x_label.font_style = 4;\n");
           }
           else
           {
             fprintf(fp,"axis([0 nn+1 ymin ymax])\n");
             fprintf(fp,"set(gca,'XTickLabel',[]);\n");
             fprintf(fp,
               "th=text(1:nn, repmat(ymin-0.07*(ymax-ymin),nn,1),Str,");
             fprintf(fp,"'HorizontalAlignment','left','rotation',90);\n");
             fprintf(fp,"set(th, 'fontsize', 12)\n");
             fprintf(fp,"set(th, 'fontweight', 'bold')\n");
           }
           fwritePlotTitle(fp,"Sobol First Order Indices");
           fwritePlotYLabel(fp, "Sobol Indices");
           fwriteHold(fp, 0);
           if (psPlotTool_ == 1)
           {
             fprintf(fp, "drawnow\n");
             printf("rssobol1 plot file = scilabrssobol1.sci\n");
           }
           else printf("rssobol1 plot file = matlabrssobol1.m\n");
           fclose(fp);
         }
       }
       else
       {
         printf("Total variance = 0. Hence, no main effect plot.\n");
       }
       pdata->clean();
     }
     psuadeIO->updateAnalysisSection(-1,-1,-1,diagKeep,-1,-1);
     delete anaManager;
   }

   // +++ rssobol2 
   else if (!strcmp(command, "rssobol2"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rssobol2: RS-based 2-input Sobol' decomposition\n");
       printf("syntax: rssobol2 (no argument needed)\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }
     if (nInputs <= 2)
     {
       printf("INFO: no point doing this for nInputs <= 2.\n");
       return -1;
     }
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;

     int analysisMethod = PSUADE_ANA_RSSOBOL2;
     AnalysisManager *anaManager = new AnalysisManager();
     anaManager->setup(analysisMethod, 0);
     psuadeIO->getParameter("ana_diagnostics",pPtr);
     int diagKeep = pPtr.intData_;
     psuadeIO->updateAnalysisSection(-1,-1,-1,outputLevel,-1,-1);
     anaManager->analyze(psuadeIO, 0, NULL, outputID);
     psuadeIO->updateAnalysisSection(-1,-1,-1,ii,-1,-1);
     pData *pdata = psuadeIO->getAuxData();
     if (pdata->nDbles_ >= nInputs)
     {
       printEquals(PL_INFO, 0);
       printf("Two-way Interaction Effect Statistics: \n");
       if (pdata->dbleData_ > 0)
       {
         for (ii = 0; ii < nInputs; ii++)
         {
           for (jj = ii+1; jj < nInputs; jj++)
           {
             printf("Inputs %4d %4d: Sobol' interaction effect = ",
                    ii+1,jj+1);
             printf("%12.4e (unnormalized)\n",
                    pdata->dbleArray_[ii*nInputs+jj]);
           }
         }
         if (psPlotTool_ == 1)
         {
           fp = fopen("scilabrssobol2.sci", "w");
           if (fp == NULL)
           {
             printf("ERROR : cannot open file scilabrssobol2.sci\n");
           }
           else
           {
             fprintf(fp,"// This file contains Sobol' 2nd order indices\n");
             fprintf(fp,"// set sortFlag = 1 and set nn to be the number\n");
             fprintf(fp,"// of inputs to display.\n");
           }
         }
         else
         {
           fp = fopen("matlabrssobol2.m", "w");
           if (fp == NULL)
           {
             printf("ERROR : cannot open file matlabrssobol2.m\n");
           }
           else
           {
             fprintf(fp,"%% This file contains Sobol' 2nd order indices\n");
             fprintf(fp,"%% set sortFlag = 1 and set nn to be the number\n");
             fprintf(fp,"%% of inputs to display.\n");
           }
         }
         if (fp != NULL)
         {
           fprintf(fp, "sortFlag = 0;\n");
           fprintf(fp, "nn = %d;\n", nInputs);
           fprintf(fp, "Mids = [\n");
           for (ii = 0; ii < nInputs; ii++) 
           {
             for (jj = 0; jj <= ii; jj++) fprintf(fp,"0.0\n");
             for (jj = ii+1; jj < nInputs; jj++) 
                fprintf(fp,"%24.16e\n",
                       pdata->dbleArray_[ii*nInputs+jj]/pdata->dbleData_);
           }
           fprintf(fp, "];\n");
           fprintf(fp, "Mids = Mids';\n");
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
               if (inputNames[ii] != NULL) 
                    fprintf(fp,"'%s',",inputNames[ii]);
               else fprintf(fp,"'X%d',",ii+1);
             }
             if (inputNames[nInputs-1] != NULL)
                  fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
             else fprintf(fp,"'X%d'};\n",nInputs);
           }
           fwriteHold(fp, 0);
           fprintf(fp, "ymin = min(Mids);\n");
           fprintf(fp, "ymax = max(Mids);\n");
           fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
           if (psPlotTool_ == 1)
           {
             fprintf(fp, "Mids = matrix(Mids, %d, %d);\n",nInputs,nInputs);
             fprintf(fp, "Mids = Mids';\n");
             fprintf(fp, "drawlater\n");
             fprintf(fp, "hist3d(Mids);\n");
             fprintf(fp, "a=gca();\n");
             fprintf(fp, "a.data_bounds=[0, 0, 0; %d+1, %d+1, ymax];\n",
                     nInputs,nInputs);
             fprintf(fp, "newtick = a.x_ticks;\n");
             fprintf(fp, "newtick(2) = [1:%d]';\n",nInputs);
             fprintf(fp, "newtick(3) = Str';\n");
             fprintf(fp, "a.x_ticks = newtick;\n");
             fprintf(fp, "a.x_label.font_size = 3;\n");
             fprintf(fp, "a.x_label.font_style = 4;\n");
             fprintf(fp, "a.y_ticks = newtick;\n");
             fprintf(fp, "a.y_label.font_size = 3;\n");
             fprintf(fp, "a.y_label.font_style = 4;\n");
             fprintf(fp, "a.rotation_angles = [5 -70];\n");
             fprintf(fp, "drawnow\n");
           }
           else
           {
             fprintf(fp, "Mids = reshape(Mids, %d, %d);\n",nInputs,nInputs);
             fprintf(fp, "Mids = Mids';\n");
             fprintf(fp, "ymin = min(min(Mids));\n");
             fprintf(fp, "ymax = max(max(Mids));\n");
             fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
             fprintf(fp, "hh = bar3(Mids,0.8);\n");
             fprintf(fp, "alpha = 0.2;\n");
             fprintf(fp, "set(hh,'FaceColor','b','facea',alpha);\n");
             fprintf(fp, "axis([0.5 %d+0.5 0.5 %d+0.5 0 ymax])\n",
                     nInputs,nInputs);
             //fprintf(fp, "bar3(Mids,0.8);\n");
             //fprintf(fp, "axis([0 %d+1 0 %d+1 0 ymax])\n",nInputs,nInputs);
             fprintf(fp, "set(gca,'XTickLabel',Str);\n");
             fprintf(fp, "set(gca,'YTickLabel',Str);\n");
             fwritePlotAxesNoGrid(fp);
           }
           fwritePlotAxes(fp);
           fwritePlotTitle(fp,"Sobol Second Order Indices (+ first order)");
           fwritePlotZLabel(fp, "Sobol Indices");
           fwritePlotXLabel(fp, "Inputs");
           fwritePlotYLabel(fp, "Inputs");
           fclose(fp);
           if (psPlotTool_ == 1)
                printf("rssobol2 plot file = scilabrssobol2.sci\n");
           else printf("rssobol2 plot file = matlabrssobol2.m\n");
         }
       }
       else
       {
         printf("Total variance = 0. Hence, no interaction effect plot.\n");
       }
       pdata->clean();
     }
     psuadeIO->updateAnalysisSection(-1,-1,-1,diagKeep,-1,-1);
     delete anaManager;
   }

   // +++ rssoboltsi 
   else if (!strcmp(command, "rssoboltsi"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rssoboltsi: RS-based Sobol' total sensitivity indices\n");
       printf("syntax: rssoboltsi (no argument needed)\n");
       printf("note: to compute error bars for indices, turn on\n");
       printf("      ana_expert mode and set ntimes to 100.\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }
     if (nInputs < 2)
     {
       printf("INFO: no point doing this for nInputs < 2.\n");
       return -1;
     }
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;

     int analysisMethod = PSUADE_ANA_RSSOBOLTSI;
     AnalysisManager *anaManager = new AnalysisManager();
     anaManager->setup(analysisMethod, 0);
     psuadeIO->getParameter("ana_diagnostics",pPtr);
     int diagKeep = pPtr.intData_;
     psuadeIO->updateAnalysisSection(-1,-1,-1,outputLevel,-1,-1);
     anaManager->analyze(psuadeIO, 0, NULL, outputID);
     psuadeIO->updateAnalysisSection(-1,-1,-1,diagKeep,-1,-1);
     delete anaManager;
   }

   // +++ rssobolg 
   else if (!strcmp(command, "rssobolg"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rssobolg: RS-based group Sobol' decomposition\n");
       printf("syntax: rssobolg (no argument needed)\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;

     int analysisMethod = PSUADE_ANA_RSSOBOLG;
     AnalysisManager *anaManager = new AnalysisManager();
     anaManager->setup(analysisMethod, 0);
     psuadeIO->getParameter("ana_diagnostics",pPtr);
     int diagKeep = pPtr.intData_;
     psuadeIO->updateAnalysisSection(-1,-1,-1,outputLevel,-1,-1);
     anaManager->analyze(psuadeIO, 0, NULL, outputID);
     psuadeIO->updateAnalysisSection(-1,-1,-1,diagKeep,-1,-1);
     fgets(lineIn,500,stdin); 
     delete anaManager;
   }

   // +++ rssobol1b 
   else if (!strcmp(command, "rssobol1b"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rssobol1b: RS-based Sobol' sensitivity analysis\n");
         printf("syntax: rssobol1b (no argument needed)\n");
         printf("Note: This command computes the first order Sobol'\n");
         printf("      sensitivity indices using response surface\n");
         printf("      constructed from the loaded sampling. It also\n");
         printf("      uses bootstrapping to estimate response surface\n");
         printf("      errors.\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data (load sample first).\n");
         return -1;
      }
      if (nSamples < 5)
      {
         printf("INFO: This command is not suitable for small samples.\n");
         printf("      nSamples needs to be at least 5.\n");
         return -1;
      }
      printf("rssobol1b INFO: MAKE SURE YOU SET YOUR DESIRED RESPONSE\n");
      printf("                SURFACE TYPE IN YOUR DATA FILE.\n");
      outputID = 0;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
      outputID = getInt(1, nOutputs, pString);
      outputID--;

      int analysisMethod = PSUADE_ANA_RSSOBOL1;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
 
      psuadeIO->getParameter("ana_diagnostics",pPtr);
      int saveDiag = pPtr.intData_;
      psuadeIO->updateAnalysisSection(-1,-1,-1,-2,-1,-1);
 
      psuadeIO->getParameter("ana_rstype", pPtr);
      int saveRS = pPtr.intData_;
      if (saveRS == PSUADE_RS_MARSB)
      {
         printf("rssobol1b INFO: MarsBagg response surface selected but\n");
         printf("                it is redundant - use MARS instead.\n");
         psuadeIO->updateAnalysisSection(-1,-1,PSUADE_RS_MARS,-3,-1,-1);
      }

      int saveMode = psAnaExpertMode_;
      psAnaExpertMode_ = 0;
      if (psRSExpertMode_ == 1)
      {
         printf("rssobol1b INFO: since RS expert mode has been enabled,\n");
         printf("     you should expect to be asked to set RS parameters\n");
         printf("     many times. If you would not prefer to be asked so\n");
         printf("     many times, use config file to set RS parameters.\n");
      }

      double *tempX = new double[nSamples*nInputs];
      double *tempY = new double[nSamples];
      int    *tempI = new int[nSamples];
      int    *states = new int[nSamples];

      sprintf(pString, "How many bootstrapped samples to use (10 - 300) : ");
      int count = getInt(10, 300, pString);
      double *tempT = new double[count*nInputs];

      int nSamples2, ind;
      printf("NOTE: the response surface type will be obtained from your data\n");
      printf("      file. If you would like to finetune your response surface,\n");
      printf("      do not turn on rs_expert mode, set the parameters in the\n");
      printf("      config file.\n");
      for (kk = 0; kk < count; kk++)
      {
         printf("rssobol1b: iteration %d\n", kk+1);
         for (jj = 0; jj < nSamples; jj++) tempI[jj] = 0;
         ss = nSamples2 = 0;
         while (ss < nSamples)
         {
            ind = PSUADE_rand() % nSamples;
            if (tempI[ind] == 0)
            {
               for (ii = 0; ii < nInputs; ii++)
                  tempX[nSamples2*nInputs+ii] = sampleInputs[ind*nInputs+ii];
               tempY[nSamples2] = sampleOutputs[ind*nOutputs+outputID];
               states[nSamples2] = sampleStates[ind];
               tempI[ind] = 1;
               nSamples2++;
            }
            ss++;
         }
         psuadeIO->updateInputSection(nSamples2,nInputs,NULL,NULL,
                                      NULL,tempX,NULL,NULL,NULL,NULL,NULL); 
         psuadeIO->updateOutputSection(nSamples2,1,tempY,states,outputNames);
         psuadeIO->updateMethodSection(PSUADE_SAMP_MC,nSamples2,-1,-1,-1);

         anaManager->analyze(psuadeIO, 0, NULL, 0);
         pData *pdata = psuadeIO->getAuxData(); 
         if (pdata->nDbles_ != nInputs)
         {
            printf("ERROR: nInputs do not match (%d, %d).\n",
                   pdata->nDbles_, nInputs);
            printf("       Consult PSUADE developers.\n");
            delete [] tempT;
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            delete [] states;
            return -1;
         }

         if (pdata->dbleData_ > 0)
            for (ii = 0; ii < nInputs; ii++)
               tempT[kk*nInputs+ii] = 
                    pdata->dbleArray_[ii]/pdata->dbleData_;
         else
            for (ii = 0; ii < nInputs; ii++)
               tempT[kk*nInputs+ii] = pdata->dbleArray_[ii];

         pdata->clean();
      }
      delete [] tempX;
      delete [] tempY;
      delete [] tempI;
      delete [] states;
      double *tempM = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         tempM[ii] = tempT[ii];
         for (jj = 1; jj < count; jj++) tempM[ii] += tempT[jj*nInputs+ii];
         tempM[ii] /= (double) count;
      }
      double *tempV = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         tempV[ii] = pow(tempT[ii]-tempM[ii], 2.0);
         for (jj = 1; jj < count; jj++) 
            tempV[ii] += pow(tempT[jj*nInputs+ii]-tempM[ii],2.0);
         tempV[ii] /= (double) (count - 1);
         tempV[ii] = sqrt(tempV[ii]);
      }
      printEquals(PL_INFO, 0);
      printf("Sobol1 Statistics (based on %d replications): \n",count);
      for (ii = 0; ii < nInputs; ii++)
         printf("Input %4d: mean = %16.8e, std = %16.8e\n",ii+1,
                tempM[ii],tempV[ii]);

      if (psPlotTool_ == 1) fp = fopen("scilabrssobol1b.sci","w");
      else                  fp = fopen("matlabrssobol1b.m","w");
      if (fp == NULL) printf("ERROR: cannot open plot file.\n");
      else
      {
         strcpy(pString," This file contains first order Sobol' indices");
         fwriteComment(fp, pString);
         strcpy(pString," with error bars coming from bootstrapping.");
         fwriteComment(fp, pString);
         strcpy(pString," to select the most important ones to display,");
         fwriteComment(fp, pString);
         strcpy(pString," set sortFlag = 1 and set nn to be the number");
         fwriteComment(fp, pString);
         strcpy(pString," of inputs to display.\n");
         fwriteComment(fp, pString);
         fprintf(fp, "sortFlag = 0;\n");
         fprintf(fp, "nn = %d;\n", nInputs);
         fprintf(fp, "Means = [\n");
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempM[ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "Stds = [\n");
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempV[ii]);
         fprintf(fp, "];\n");
         if (inputNames == NULL)
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
               if (inputNames[ii] != NULL) 
                    fprintf(fp,"'%s',",inputNames[ii]);
               else fprintf(fp,"'X%d',",ii+1);
            }
            if (inputNames[nInputs-1] != NULL) 
                 fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
            else fprintf(fp,"'X%d'};\n",nInputs);
         }
         fwriteHold(fp, 0);
         fprintf(fp, "if (sortFlag == 1)\n");
         if (psPlotTool_ == 1)
              fprintf(fp, "  [Means, I2] = gsort(Means);\n");
         else fprintf(fp, "  [Means, I2] = sort(Means,'descend');\n");
         fprintf(fp, "  Stds = Stds(I2);\n");
         fprintf(fp, "  I2 = I2(1:nn);\n");
         fprintf(fp, "  Means = Means(1:nn);\n");
         fprintf(fp, "  Stds = Stds(1:nn);\n");
         fprintf(fp, "  Str  = Str(I2);\n");
         fprintf(fp, "end\n");
         fprintf(fp, "ymin = min(Means-Stds);\n");
         fprintf(fp, "if ymin < 0 \n");
         fprintf(fp, "   ymin = 0;\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "ymax = max(Means+Stds);\n");
         fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
         if (psPlotTool_ == 1) fprintf(fp, "drawlater\n");
         fprintf(fp, "bar(Means,0.8);\n");
         fprintf(fp, "for ii = 1:nn\n");
         fprintf(fp, "   if (ii == 1)\n");
         fwriteHold(fp, 1);
         fprintf(fp, "   end;\n");
         fprintf(fp, "   XX = [ii ii];\n");
         fprintf(fp, "   d1 = Means(ii)-Stds(ii);\n");
         fprintf(fp, "   d2 = Means(ii)+Stds(ii);\n");
         fprintf(fp, "   if (d1 < 0)\n");
         fprintf(fp, "      d1 = 0.0;\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "   YY = [d1 d2];\n");
         fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
         fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',13)\n");
         fprintf(fp, "end;\n");
         fwritePlotAxes(fp);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
            fprintf(fp, "newtick = a.x_ticks;\n");
            fprintf(fp, "newtick(2) = [1:nn]';\n");
            fprintf(fp, "newtick(3) = Str';\n");
            fprintf(fp, "a.x_ticks = newtick;\n");
            fprintf(fp, "a.x_label.font_size = 3;\n");
            fprintf(fp, "a.x_label.font_style = 4;\n");
         }
         else
         {
            fprintf(fp, "axis([0  nn+1 ymin ymax])\n");
            fprintf(fp, "set(gca,'XTickLabel',[]);\n");
            fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str,");
            fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
            fprintf(fp, "set(th, 'fontsize', 12)\n");
            fprintf(fp, "set(th, 'fontweight', 'bold')\n");
         }
         fwritePlotTitle(fp,"First Order Sobol Indices (with bootstrap)");
         fwritePlotYLabel(fp, "First Order Sobol Index (Normalized)");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "drawnow\n");
            printf("rssobol1b plot file = scilabrssobol1b.sci\n");
         }
         else printf("rssobol1b plot file = matlabrssobol1b.m\n");
         fclose(fp);
      }

      delete anaManager;
      delete [] tempT;
      delete [] tempM;
      delete [] tempV;

      psAnaExpertMode_ = saveMode;
      if (saveRS == PSUADE_RS_MARSB)
         psuadeIO->updateAnalysisSection(-1,-1,saveRS,-3,-1,-1);
      psuadeIO->updateAnalysisSection(-1,-1,-1,saveDiag,-1,-1);
   }

   // +++ rssobol2b 
   else if (!strcmp(command, "rssobol2b"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rssobol2b: RS-based 2-input Sobol' method with bootstrapping\n");
         printf("syntax: rssobol2b (no argument needed)\n");
         printf("Note: This command computes the second order Sobol'\n");
         printf("      sensitivity indices using response surface\n");
         printf("      constructed from the loaded sampling. It also\n");
         printf("      uses bootstrapping to estimate response surface\n");
         printf("      errors.\n");
         return 0;
      }
      if (nInputs <= 0 || psuadeIO == NULL)
      {
         printf("ERROR: data file not loaded.\n");
         return -1;
      }
      if (nInputs <= 2)
      {
         printf("INFO: no point doing this for nInputs <= 2.\n");
         return -1;
      }
      if (nSamples < (nInputs+1)*5/3+1)
      {
         printf("INFO: This command is not suitable for small samples.\n");
         printf("      nSamples needs to be at least %d.\n",(nInputs+1)/3+1);
         printf("      nSamples needs to be at least %d.\n",(nInputs+1)/3+1);
         return -1;
      }
      printf("rssobol2b INFO: MAKE SURE YOU SET YOUR DESIRED RESPONSE\n");
      printf("          SURFACE TYPE IN YOUR DATA FILE.\n");
      outputID = 0;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
      outputID = getInt(1, nOutputs, pString);
      outputID--;

      int analysisMethod = PSUADE_ANA_RSSOBOL2;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
 
      psuadeIO->getParameter("ana_diagnostics",pPtr);
      int saveDiag = pPtr.intData_;
      psuadeIO->updateAnalysisSection(-1,-1,-1,-2,-1,-1);
      int saveMode = psAnaExpertMode_;
      psAnaExpertMode_ = 0;
      if (psRSExpertMode_ == 1)
      {
         printf("rssobol2b INFO: since RS expert mode has been enabled,\n");
         printf("     you should expect to be asked to set RS parameters\n");
         printf("     many times. If you would not prefer to be asked so\n");
         printf("     many times, use config file to set RS parameters.\n");
      }

      psuadeIO->getParameter("ana_rstype", pPtr);
      int saveRS = pPtr.intData_;
      if (saveRS == PSUADE_RS_MARSB)
      {
         printf("rssobol2b INFO: MarsBagg response surface selected but\n");
         printf("                it is redundant - use MARS instead.\n");
         psuadeIO->updateAnalysisSection(-1,-1,PSUADE_RS_MARS,-3,-1,-1);
      }

      int    ind, nSamples2;
      double *tempX = new double[nSamples*nInputs];
      double *tempY = new double[nSamples];
      int    *tempI = new int[nSamples];
      int    *states = new int[nSamples];

      sprintf(pString, "How many bootstrapped samples to use (10 - 300) : ");
      int count = getInt(10, 300, pString);
      double *tempT = new double[(count+3)*nInputs*nInputs];
      
      printf("NOTE: the response surface type will be obtained from your\n");
      printf("      data file. If you desire to fine-tune your response\n");
      printf("      surface, do not turn on rs_expert mode. Instead, set\n");
      printf("      the parameters in the config file.\n");
      for (kk = 0; kk < count; kk++)
      {
         printf(">>>>>>>>>>>>>>>> rssobol2b: ITERATION %d\n", kk+1);
         for (jj = 0; jj < nSamples; jj++) tempI[jj] = 0;
         ss = nSamples2 = 0;
         while (ss < nSamples)
         {
            ind = PSUADE_rand() % nSamples;
            if (tempI[ind] == 0)
            {
               for (ii = 0; ii < nInputs; ii++)
                  tempX[nSamples2*nInputs+ii] = sampleInputs[ind*nInputs+ii];
               tempY[nSamples2] = sampleOutputs[ind*nOutputs+outputID];
               states[nSamples2] = sampleStates[ind];
               tempI[ind] = 1;
               nSamples2++;
            }
            ss++;
         }
         psuadeIO->updateInputSection(nSamples2,nInputs,NULL,NULL,
                                      NULL,tempX,NULL,NULL,NULL,NULL,NULL);
         psuadeIO->updateOutputSection(nSamples2,1,tempY,states,
                                        outputNames);
         psuadeIO->updateMethodSection(PSUADE_SAMP_MC,nSamples2,
                                       -1,-1,-1);
         status = anaManager->analyze(psuadeIO, 0, NULL, 0);
         if (status != 0)
         {
            printf("rssobol2b ERROR: not successful.\n");
            delete [] tempT;
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            delete [] states;
            delete anaManager;
            return -1;
         }
         pData *pdata = psuadeIO->getAuxData();
         if (pdata->nDbles_ < nInputs)
         {
            printf("rssobol2b SYSTEM ERROR: nInputs do not match (%d, %d).\n",
                   pdata->nDbles_, nInputs);
            printf("       Consult PSUADE developers.\n");
            delete [] tempT;
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            delete [] states;
            delete anaManager;
            return -1;
         }
         if (pdata->dbleData_ > 0)
         {
            for (ii = 0; ii < nInputs*nInputs; ii++)
               tempT[kk*nInputs*nInputs+ii] =
                  pdata->dbleArray_[ii]/pdata->dbleData_;
         }
         else
         {
            for (ii = 0; ii < nInputs*nInputs; ii++)
               tempT[kk*nInputs*nInputs+ii] = pdata->dbleArray_[ii];
         }
         pdata->clean();
      }
      delete [] tempX;
      delete [] tempY;
      delete [] states;
      delete [] tempI;
      delete anaManager;
      
      for (ii = 0; ii < nInputs; ii++)
      {
         for (jj = 0; jj <= ii; jj++) 
            tempT[count*nInputs*nInputs+ii*nInputs+jj] = 0.0;
         for (jj = ii+1; jj < nInputs; jj++)
         {
            ddata = 0.0;
            for (kk = 0; kk < count; kk++)
               ddata += tempT[kk*nInputs*nInputs+ii*nInputs+jj]; 
            tempT[count*nInputs*nInputs+ii*nInputs+jj] = ddata/(double) count;
            tempT[(count+1)*nInputs*nInputs+ii*nInputs+jj] =  
                                tempT[ii*nInputs+jj]; 
            tempT[(count+2)*nInputs*nInputs+ii*nInputs+jj] =  
                                      tempT[ii*nInputs+jj]; 
            for (kk = 1; kk < count; kk++)
            {
               ddata = tempT[kk*nInputs*nInputs+ii*nInputs+jj]; 
               if (ddata < tempT[(count+1)*nInputs*nInputs+ii*nInputs+jj]) 
                  tempT[(count+1)*nInputs*nInputs+ii*nInputs+jj] = ddata; 
               if (ddata > tempT[(count+2)*nInputs*nInputs+ii*nInputs+jj]) 
                  tempT[(count+2)*nInputs*nInputs+ii*nInputs+jj] = ddata; 
            }
         }
      }
      printDashes(PL_INFO, 0);
      printf("Sobol2b Statistics (based on %d replications): \n", count);
      printf("Quantities are normalized.\n");
      printf("Multiply by variance to compute actual statistics.\n");
      printEquals(PL_INFO, 0);
      for (ii = 0; ii < nInputs; ii++)
      {
         for (jj = 0; jj <= ii; jj++) tempT[ii*nInputs+jj] = 0.0;
         for (jj = ii+1; jj < nInputs; jj++)
         {
            ddata = 0.0;
            for (kk = 0; kk < count; kk++)
            {
               tempT[kk*nInputs*nInputs+ii*nInputs+jj] -=  
                  tempT[count*nInputs*nInputs+ii*nInputs+jj];
               ddata += pow(tempT[kk*nInputs*nInputs+ii*nInputs+jj],2.0);
            }
            ddata /= (double) (count - 1);
            tempT[ii*nInputs+jj] = ddata;
            printf("Input (%4d %4d): mean = %16.8e, std = %16.8e\n",ii+1,
                  jj+1,tempT[count*nInputs*nInputs+ii*nInputs+jj],ddata);
            printf("Input (%4d %4d): min  = %16.8e, max = %16.8e\n",ii+1,
                  jj+1,tempT[(count+1)*nInputs*nInputs+ii*nInputs+jj],
                  tempT[(count+2)*nInputs*nInputs+ii*nInputs+jj]);
         }
      }
      printEquals(PL_INFO, 0);
      if (psPlotTool_ == 1)
      {
         fp = fopen("scilabrssobol2b.sci", "w");
         if (fp == NULL) printf("ERROR : cannot open file scilabrssobol2b.sci\n");
         else
         {
            fprintf(fp,"// This file contains Sobol' 2nd order indices\n");
            fprintf(fp,"// set sortFlag = 1 and set nn to be the number\n");
            fprintf(fp,"// of inputs to display.\n");
         }
      }
      else
      {
         fp = fopen("matlabrssobol2b.m", "w");
         if (fp == NULL) printf("ERROR : cannot open file matlabrssobol2b.sci\n");
         else
         {
            fprintf(fp, "%% This file contains Sobol' 2nd order indices\n");
            fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
            fprintf(fp, "%% of inputs to display.\n");
         }
      }
      if (fp != NULL) 
      {
         fprintf(fp, "sortFlag = 0;\n");
         fprintf(fp, "nn = %d;\n", nInputs);
         fprintf(fp, "Means = [\n");
         for (ii = 0; ii < nInputs*nInputs; ii++) 
            fprintf(fp,"%24.16e\n", tempT[count*nInputs*nInputs+ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "Stds = [\n");
         for (ii = 0; ii < nInputs*nInputs; ii++) 
            fprintf(fp,"%24.16e\n", tempT[ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "Lows = [\n");
         for (ii = 0; ii < nInputs*nInputs; ii++) 
            fprintf(fp,"%24.16e\n", tempT[(count+1)*nInputs*nInputs+ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "Highs = [\n");
         for (ii = 0; ii < nInputs*nInputs; ii++) 
            fprintf(fp,"%24.16e\n", tempT[(count+2)*nInputs*nInputs+ii]);
         fprintf(fp, "];\n");
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
               if (inputNames[ii] != NULL) 
                  fprintf(fp,"'%s',",inputNames[ii]);
               else fprintf(fp,"'X%d',",ii+1);
            }
            if (inputNames[nInputs-1] != NULL)
               fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
            else fprintf(fp,"'X%d'};\n",nInputs);
         }
         fwriteHold(fp, 0);
         fprintf(fp, "ymin = min(Means-Stds);\n");
         fprintf(fp, "ymax = max(Means+Stds);\n");
         fprintf(fp, "ymin = min(Lows);\n");
         fprintf(fp, "ymax = max(Highs);\n");
         fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "nn    = %d;\n",nInputs);
            fprintf(fp, "Means = matrix(Means, nn, nn);\n");
            fprintf(fp, "Means = Means';\n");
            fprintf(fp, "Stds  = matrix(Stds, nn, nn);\n");
            fprintf(fp, "Stds  = Stds';\n");
            fprintf(fp, "Lows  = matrix(Lows, nn, nn);\n");
            fprintf(fp, "Lows  = Lows';\n");
            fprintf(fp, "Highs = matrix(Highs, nn, nn);\n");
            fprintf(fp, "Highs = Highs';\n");
            fprintf(fp, "drawlater\n");
            fprintf(fp, "hist3d(Means);\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "//for ii = 1:nn\n");
            fprintf(fp, "//  for jj = ii+1:nn\n");
            fprintf(fp, "//    XX = [ii ii];\n");
            fprintf(fp, "//    YY = [jj jj];\n");
            fprintf(fp, "//    MM = Means(ii,jj);\n");
            fprintf(fp, "//    SS = Stds(ii,jj);\n");
            fprintf(fp, "//    ZZ = [MM-SS MM+SS];\n");
            fprintf(fp, "//    plot3d(XX,YY,ZZ,'-ko','LineWidth',3.0,");
            fprintf(fp, "//      'MarkerEdgeColor','k','MarkerFaceColor',");
            fprintf(fp, "//      'g','MarkerSize',13)\n");
            fprintf(fp, "//  end;\n");
            fprintf(fp, "//end;\n");
            fprintf(fp, "//a=gca();\n");
            fprintf(fp, "//a.data_bounds=[0, 0, 0; nn, nn+1, ymax];\n");
            fprintf(fp, "//newtick = a.x_ticks;\n");
            fprintf(fp, "//newtick(2) = [1:nn]';\n");
            fprintf(fp, "//drawlater\n");
            fprintf(fp, "//hist3d(Means);\n");
            fprintf(fp, "//set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "//for ii = 1:nn\n");
            fprintf(fp, "//  for jj = ii+1:nn\n");
            fprintf(fp, "//    XX = [ii ii];\n");
            fprintf(fp, "//    YY = [jj jj];\n");
            fprintf(fp, "//    MM = Means(ii,jj);\n");
            fprintf(fp, "//    SS = Stds(ii,jj);\n");
            fprintf(fp, "//    ZZ = [MM-SS MM+SS];\n");
            fprintf(fp, "//    plot3d(XX,YY,ZZ,'-ko','LineWidth',3.0,");
            fprintf(fp, "//      'MarkerEdgeColor','k','MarkerFaceColor',");
            fprintf(fp, "//      'g','MarkerSize',13)\n");
            fprintf(fp, "//  end;\n");
            fprintf(fp, "//end;\n");
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[0, 0, 0; nn, nn+1, ymax];\n");
            fprintf(fp, "newtick = a.x_ticks;\n");
            fprintf(fp, "newtick(2) = [1:nn]';\n");
            fprintf(fp, "newtick(3) = Str';\n");
            fprintf(fp, "a.x_ticks = newtick;\n");
            fprintf(fp, "a.x_label.font_size = 3;\n");
            fprintf(fp, "a.x_label.font_style = 4;\n");
            fprintf(fp, "a.y_ticks = newtick;\n");
            fprintf(fp, "a.y_label.font_size = 3;\n");
            fprintf(fp, "a.y_label.font_style = 4;\n");
            fprintf(fp, "a.rotation_angles = [5 -70];\n");
            fprintf(fp, "drawnow\n");
         }
         else
         {
            fprintf(fp, "nn    = %d;\n",nInputs);
            fprintf(fp, "Means = reshape(Means, nn, nn);\n");
            fprintf(fp, "Means = Means';\n");
            fprintf(fp, "Stds  = reshape(Stds, nn, nn);\n");
            fprintf(fp, "Stds  = Stds';\n");
            fprintf(fp, "Lows  = reshape(Lows, nn, nn);\n");
            fprintf(fp, "Lows  = Lows';\n");
            fprintf(fp, "Highs = reshape(Highs, nn, nn);\n");
            fprintf(fp, "Highs = Highs';\n");
            fprintf(fp, "hh = bar3(Means,0.8);\n");
            fprintf(fp, "alpha = 0.2;\n");
            fprintf(fp, "set(hh,'FaceColor','b','facea',alpha);\n");
            fprintf(fp, "Lstds = Means - Stds;\n");
            fprintf(fp, "Ustds = Means + Stds;\n");
            fprintf(fp, "Lstds = Lows;\n");
            fprintf(fp, "Ustds = Highs;\n");
            fprintf(fp, "[X,Y] = meshgrid(1:nn,1:nn);\n");
            fwriteHold(fp, 1);
            fprintf(fp, "for k = 1:nn\n");
            fprintf(fp, "  for l = k:nn\n");
            fprintf(fp, "    mkl = Means(k,l);\n");
            fprintf(fp, "    ukl = Ustds(k,l);\n");
            fprintf(fp, "    lkl = Lstds(k,l);\n");
            fprintf(fp, "    if (mkl > .02 & (ukl-lkl)/mkl > .02)\n");
            fprintf(fp, "      xkl = [X(k,l), X(k,l)];\n");
            fprintf(fp, "      ykl = [Y(k,l), Y(k,l)];\n");
            fprintf(fp, "      zkl = [lkl, ukl];\n");
            fprintf(fp, "      plot3(xkl,ykl,zkl,'-mo',...\n");
            fprintf(fp, "        'LineWidth',5,'MarkerEdgeColor','k',...\n");
            fprintf(fp, "        'MarkerFaceColor','k','MarkerSize',10);\n");
            fprintf(fp, "    end\n");
            fprintf(fp, "  end\n");
            fprintf(fp, "end\n");
            fwriteHold(fp, 0);
            fprintf(fp, "axis([0.5 nn+0.5 0.5 nn+0.5 0 ymax])\n");
            fprintf(fp, "set(gca,'XTickLabel',Str);\n");
            fprintf(fp, "set(gca,'YTickLabel',Str);\n");
            fprintf(fp, "set(gca, 'fontsize', 12)\n");
            fprintf(fp, "set(gca, 'fontweight', 'bold')\n");
            fprintf(fp, "set(gca, 'linewidth', 2)\n");
         }
         fwritePlotAxes(fp);
         fwritePlotTitle(fp,"Sobol 1st+2nd Order Indices (with bootstrap)");
         fwritePlotZLabel(fp, "Sobol Indices (Normalized)");
         fwritePlotXLabel(fp, "Inputs");
         fwritePlotYLabel(fp, "Inputs");
         fclose(fp);
         if (psPlotTool_ == 1)
              printf("rssobol2b plot file = scilabrssobol2b.sci\n");
         else printf("rssobol2b plot file = matlabrssobol2b.m\n");
      }
      psAnaExpertMode_ = saveMode;
      if (saveRS == PSUADE_RS_MARSB)
         psuadeIO->updateAnalysisSection(-1,-1,saveRS,-3,-1,-1);
      psuadeIO->updateAnalysisSection(-1,-1,-1,saveDiag,-1,-1);
      delete [] tempT;
   }

   // +++ rssoboltsib 
   else if (!strcmp(command, "rssoboltsib"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rssoboltsib: RS-based Sobol' total sensitivity analysis\n");
         printf("syntax: rssoboltsib (no argument needed)\n");
         printf("Note: This command computes the total order Sobol'\n");
         printf("      sensitivity indices using response surface\n");
         printf("      constructed from the loaded sampling. It also\n");
         printf("      uses bootstrapping to estimate response surface\n");
         printf("      errors.\n");
         return 0;
      }
      if (nInputs <= 0 || psuadeIO == NULL)
      {
         printf("ERROR: data file not loaded.\n");
         return -1;
      }
      if (nInputs < 2)
      {
         printf("INFO: no point doing this for nInputs < 2.\n");
         return -1;
      }
      if (nSamples < 10)
      {
         printf("WARNING: your sample size is quite small.\n");
         printf("         Bootstrapped samples will be smaller.\n");
         return -1;
      }
      printf("rssoboltsi INFO: MAKE SURE TO SET THE DESIRED RESPONSE\n");
      printf("                 SURFACE TYPE IN YOUR DATA FILE.\n");
      outputID = 0;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
      outputID = getInt(1, nOutputs, pString);
      outputID--;

      int analysisMethod = PSUADE_ANA_RSSOBOLTSI;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
 
      psuadeIO->getParameter("ana_diagnostics",pPtr);
      int saveDiag = pPtr.intData_;
      psuadeIO->updateAnalysisSection(-1,-1,-1,-2,-1,-1);
      int saveMode = psAnaExpertMode_;
      psAnaExpertMode_ = 0;
      if (psRSExpertMode_ == 1)
      {
         printf("rssoboltsib INFO: since RS expert mode has been enabled,\n");
         printf("     you should expect to be asked to set RS parameters\n");
         printf("     many times. If you would not prefer to be asked so\n");
         printf("     many times, use config file to set RS parameters.\n");
      }

      psuadeIO->getParameter("ana_rstype", pPtr);
      int saveRS = pPtr.intData_;
      if (saveRS == PSUADE_RS_MARSB)
      {
         printf("rssoboltsib INFO: MarsBagg response surface selected but\n");
         printf("                  it is redundant - use MARS instead.\n");
         psuadeIO->updateAnalysisSection(-1,-1,PSUADE_RS_MARS,-3,-1,-1);
      }

      int    ind, nSamples2;
      pData  *pdata;
      double *tempX = new double[nSamples*nInputs];
      double *tempY = new double[nSamples];
      int    *tempI = new int[nSamples];
      int    *states = new int[nSamples];

      sprintf(pString,"How many bootstrapped samples (10 - 300) ? ");
      int count = getInt(10, 1000, pString);
      double *tempT = new double[count*nInputs];

      printf("NOTE: the response surface type will be obtained from\n");
      printf("      your data file. If you would like to fine-tune\n");
      printf("      your response surface, do not turn on rs_expert\n");
      printf("      mode. Instead, set the parameters in the\n");
      printf("      config file.\n");
      for (kk = 0; kk < count; kk++)
      {
         printf("rssoboltsib: iteration %d\n", kk+1);
         for (jj = 0; jj < nSamples; jj++) tempI[jj] = 0;
         ss = nSamples2 = 0;
         while (ss < nSamples)
         {
            ind = PSUADE_rand() % nSamples;
            if (tempI[ind] == 0)
            {
               for (ii = 0; ii < nInputs; ii++)
                  tempX[nSamples2*nInputs+ii] = 
                                         sampleInputs[ind*nInputs+ii];
               tempY[nSamples2] = sampleOutputs[ind*nOutputs+outputID];
               states[nSamples2] = sampleStates[ind];
               tempI[ind] = 1;
               nSamples2++;
            }
            ss++;
         }
         psuadeIO->updateInputSection(nSamples2,nInputs,NULL,NULL,
                                      NULL,tempX,NULL,NULL,NULL,NULL,NULL); 
         psuadeIO->updateOutputSection(nSamples2,1,tempY,states,
                                        outputNames);
         psuadeIO->updateMethodSection(PSUADE_SAMP_MC,nSamples2,
                                        -1,-1,-1);

         anaManager->analyze(psuadeIO, 0, NULL, 0);
         pdata = psuadeIO->getAuxData(); 
         if (pdata->nDbles_ != nInputs)
         {
            printf("ERROR: nInputs do not match (%d, %d).\n",
                   pdata->nDbles_, nInputs);
            printf("       Consult PSUADE developers.\n");
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            delete [] tempT;
            delete [] states;
            return -1;
         }

         if (pdata->dbleData_ > 0)
            for (ii = 0; ii < nInputs; ii++)
               tempT[kk*nInputs+ii] = 
                    pdata->dbleArray_[ii]/pdata->dbleData_;
         else
            for (ii = 0; ii < nInputs; ii++)
               tempT[kk*nInputs+ii] = pdata->dbleArray_[ii];

         pdata->clean();
      }
      delete [] tempX;
      delete [] tempY;
      delete [] tempI;
      delete [] states;
      double *tempM = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         tempM[ii] = tempT[ii];
         for (jj = 1; jj < count; jj++) tempM[ii] += tempT[jj*nInputs+ii];
         tempM[ii] /= (double) count;
      }
      double *tempV = new double[nInputs];
      for (ii = 0; ii < nInputs; ii++)
      {
         tempV[ii] = pow(tempT[ii]-tempM[ii], 2.0);
         for (jj = 1; jj < count; jj++) 
            tempV[ii] += pow(tempT[jj*nInputs+ii]-tempM[ii],2.0);
         tempV[ii] /= (double) (count - 1);
         tempV[ii] = sqrt(tempV[ii]);
      }
      printEquals(PL_INFO, 0);
      printf("Sobol' Total Senstivities (based on %d replications): \n",
             count);
      printf("Quantities are normalized.\n");
      printf("Multiply by variance to compute actual statistics.\n");
      for (ii = 0; ii < nInputs; ii++)
         printf("Input %4d: mean = %16.8e, std = %16.8e\n",ii+1,
                tempM[ii],tempV[ii]);

      if (psPlotTool_ == 1)
      {
         fp = fopen("scilabrssoboltsib.sci","w");
         if (fp == NULL) 
            printf("ERROR : cannot open file scilabrssoboltsib.sci\n");
         else
         {
            fprintf(fp,"// This file contains total order Sobol' indices\n");
            fprintf(fp,"// with error bars coming from bootstrapping.\n");
            fprintf(fp,"// to select the most important ones to display,\n");
            fprintf(fp,"// set sortFlag = 1 and set nn to be the number\n");
            fprintf(fp,"// of inputs to display.\n");
         }
      }
      else
      {
         fp = fopen("matlabrssoboltsib.m","w");
         if (fp == NULL) 
            printf("ERROR : cannot open file matlabrssoboltsib.sci\n");
         else
         {
            fprintf(fp,"%% This file contains total order Sobol' indices\n");
            fprintf(fp,"%% with error bars coming from bootstrapping.\n");
            fprintf(fp,"%% to select the most important ones to display,\n");
            fprintf(fp,"%% set sortFlag = 1 and set nn to be the number\n");
            fprintf(fp,"%% of inputs to display.\n");
         }
      }
      if (fp != NULL)
      {
         fprintf(fp, "sortFlag = 0;\n");
         fprintf(fp, "nn = %d;\n", nInputs);
         fprintf(fp, "Means = [\n");
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempM[ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "Stds = [\n");
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%24.16e\n",tempV[ii]);
         fprintf(fp, "];\n");
         if (inputNames == NULL)
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
               if (inputNames[ii] != NULL) 
                    fprintf(fp,"'%s',",inputNames[ii]);
               else fprintf(fp,"'X%d',",ii+1);
            }
            if (inputNames[nInputs-1] != NULL) 
                 fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
            else fprintf(fp,"'X%d'};\n",nInputs);
         }
         fwriteHold(fp, 0);
         fprintf(fp, "if (sortFlag == 1)\n");
         if (psPlotTool_ == 1)
              fprintf(fp, "  [Means, I2] = gsort(Means);\n");
         else fprintf(fp, "  [Means, I2] = sort(Means,'descend');\n");
         fprintf(fp, "  Stds = Stds(I2);\n");
         fprintf(fp, "  I2 = I2(1:nn);\n");
         fprintf(fp, "  Means = Means(1:nn);\n");
         fprintf(fp, "  Stds = Stds(1:nn);\n");
         fprintf(fp, "  Str  = Str(I2);\n");
         fprintf(fp, "end\n");
         fprintf(fp, "ymin = min(Means-Stds);\n");
         fprintf(fp, "if ymin < 0 \n");
         fprintf(fp, "    ymin = 0;\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "ymax = max(Means+Stds);\n");
         fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
         if (psPlotTool_ == 1) fprintf(fp, "drawlater\n");
         fprintf(fp, "bar(Means,0.8);\n");
         fprintf(fp, "for ii = 1:nn\n");
         fprintf(fp, "   if (ii == 1)\n");
         fwriteHold(fp, 1);
         fprintf(fp, "   end;\n");
         fprintf(fp, "   XX = [ii ii];\n");
         fprintf(fp, "   YY = [Means(ii)-Stds(ii) Means(ii)+Stds(ii)];\n");
         fprintf(fp, "   if YY(1) < 0 \n");
         fprintf(fp, "      YY(1) = 0;\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
         fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',12)\n");
         fprintf(fp, "end;\n");
         fwritePlotAxes(fp);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
            fprintf(fp, "newtick = a.x_ticks;\n");
            fprintf(fp, "newtick(2) = [1:nn]';\n");
            fprintf(fp, "newtick(3) = Str';\n");
            fprintf(fp, "a.x_ticks = newtick;\n");
            fprintf(fp, "a.x_label.font_size = 3;\n");
            fprintf(fp, "a.x_label.font_style = 4;\n");
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
         fwritePlotTitle(fp,"Total Order Sobol Indices (with bootstrap)");
         fwritePlotYLabel(fp, "Total Order Sobol Index (Normalized)");
         fwriteHold(fp, 0);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "drawnow\n");
            printf("rssoboltsib plot file = scilabrssoboltsib.sci\n");
         }
         else printf("rssoboltsib plot file = matlabrssoboltsib.m\n");
         fclose(fp);
      }
 
      delete anaManager;
      delete [] tempT;
      delete [] tempM;
      delete [] tempV;

      psAnaExpertMode_ = saveMode;
      if (saveRS == PSUADE_RS_MARSB)
         psuadeIO->updateAnalysisSection(-1,-1,saveRS,-3,-1,-1);
      psuadeIO->updateAnalysisSection(-1,-1,-1,saveDiag,-1,-1);
   }

   // +++ rsua2
   else if (!strcmp(command, "rsua2") || !strcmp(command, "rs_ua"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rsua2: uncertainty analysis on response surface\n");
         printf("syntax: rsua2 (no argument needed)\n");
         printf("This command perform uncertainty analysis on the response\n");
         printf("surface built from the loaded sample. If you select a\n");
         printf("stochastic response surface type (Kriging, MARSB, or\n");
         printf("polynomial regression, the effect of response surface\n");
         printf("uncertainty (in the average sense) will be shown on the \n");
         printf("PDF and CDF plots.\n");
         printf("NOTE: This analysis supports non-uniform distributions for\n");
         printf("      the inputs. Simply prescribe PDF in the data file\n");
         printf("      and turn on use_input_pdfs in ANALYSIS.\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data (load sample first).\n");
         return -1;
      }
      
      Sampling   *samPtr;
      FuncApprox *faPtr;
      PDFManager *pdfman;
      psVector   vecOut, vecLower, vecUpper;
      psuadeIO->getParameter("ana_rstype", pPtr);
      faType = pPtr.intData_;
      outputID = 0;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
      outputID = getInt(1, nOutputs, pString);
      outputID--;
      sprintf(pString,
          "Sample size for generating distribution? (10000 - 100000) ");
      int nSamp = getInt(10000, 100000, pString);
      flag = 0;
      printf("Save the generated sample in a file? (y or n) ");
      fgets(winput,10,stdin); 
      if (winput[0] == 'y') flag = 1; 

      printf("Phase 1 of 4: create response surface\n");
      printf("NOTE: the response surface type is taken from your data file.\n");
      faPtr = genFA(faType, nInputs, -1, nSamples);
      faPtr->setNPtsPerDim(32);
      faPtr->setBounds(iLowerB, iUpperB);
      faPtr->setOutputLevel(outputLevel);
      double *tempY;
      if (nOutputs > 1)
      {
         tempY = new double[nSamples];
         for (ss = 0; ss < nSamples; ss++) 
            tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
      }
      else tempY = sampleOutputs;
      faPtr->initialize(sampleInputs,tempY);
      if (nOutputs > 1) delete [] tempY;
      tempY = NULL;
            
      printEquals(PL_INFO, 0);
      printf("Phase 2 of 4: create MC sample\n");
      double *samInputs  = new double[nInputs*nSamp];
      double *samOutputs = new double[nSamp];
      double *samStds    = new double[nSamp];
      psuadeIO->getParameter("ana_use_input_pdfs", pPtr);
      int usePDFs = pPtr.intData_;
      if (usePDFs == 1)
      {
         printf("NOTE: Some inputs have non-uniform PDFs.\n");
         printf("      A MC sample will be created with these PDFs.\n");
         psuadeIO->getParameter("method_sampling", pPtr);
         kk = pPtr.intData_;
         psuadeIO->updateMethodSection(PSUADE_SAMP_MC,-1,-1,-1,-1);
         pdfman = new PDFManager();
         pdfman->initialize(psuadeIO);
         vecOut.setLength(nSamp*nInputs);
         vecUpper.load(nInputs, iUpperB);
         vecLower.load(nInputs, iLowerB);
         pdfman->genSample(nSamp, vecOut, vecLower, vecUpper);
         for (ii = 0; ii < nSamp*nInputs; ii++) samInputs[ii] = vecOut[ii];
         psuadeIO->updateMethodSection(kk,-1,-1,-1,-1);
         delete pdfman;
      }
      else
      {
         printf("NOTE: Uniform distributions will be used for all inputs.\n");
         printf("      To use other than uniform distributions, prescribe\n");
         printf("      them in the data file and set use_input_pdfs in the\n");
         printf("      ANALYSIS section.\n");
         int *samStates  = new int[nSamp];
         if (nInputs < 51)
              samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
         else samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
         samPtr->setPrintLevel(0);
         samPtr->setInputBounds(nInputs, iLowerB, iUpperB);
         samPtr->setOutputParams(1);
         samPtr->setSamplingParams(nSamp, -1, -1);
         samPtr->initialize(0);
         samPtr->getSamples(nSamp,nInputs,1,samInputs,samOutputs,samStates);
         delete samPtr;
         samPtr = NULL;
         delete [] samStates;
         samStates = NULL;
      }

      printf("Phase 3 of 4: evaluate sample\n");
      faPtr->evaluatePointFuzzy(nSamp,samInputs,samOutputs,samStds); 
      delete faPtr;
      faPtr = NULL;
      if (flag == 1)
      {
         if (!strcmp(command, "rsua2")) fp = fopen("rsua2_sample","w");
         else                           fp = fopen("rsua_sample","w");
         fprintf(fp, "%% inputs, output, output-3 sigma, output+3sigma\n");
         fprintf(fp, "%d %d 3\n", nSamp, nInputs);
         for (ss = 0; ss < nSamp; ss++)
         {
            for (ii = 0; ii < nInputs; ii++) 
               fprintf(fp, "%e ", samInputs[ss*nInputs+ii]);
            fprintf(fp, "%e ", samOutputs[ss]);
            fprintf(fp, "%e ", samOutputs[ss]-3*samStds[ss]);
            fprintf(fp, "%e\n", samOutputs[ss]+3*samStds[ss]);
         }
         fclose(fp);
         if (!strcmp(command, "rsua2"))
            printf("A MC sample has been written to the file 'rsua2_sample'\n");
         else
            printf("A MC sample has been written to the file 'rsua_sample'\n");
      }

      int    nbins = 100, ntimes=20;
      int    **Fcounts = new int*[ntimes+1];
      double Fmax=-PSUADE_UNDEFINED;
      double Fmin=PSUADE_UNDEFINED;
      PDFNormal *rsPDF;

      printf("Phase 4 of 4: binning\n");
      for (ss = 0; ss < nSamp; ss++)
      {
         if (samOutputs[ss] > Fmax) Fmax = samOutputs[ss];
         if (samOutputs[ss] < Fmin) Fmin = samOutputs[ss];
         if (samOutputs[ss]+3*samStds[ss] > Fmax) 
            Fmax = samOutputs[ss] + 3 * samStds[ss];
         if (samOutputs[ss]-3*samStds[ss] < Fmin) 
            Fmin = samOutputs[ss] - 3 * samStds[ss];
      }
      Fmax = Fmax + 0.1 * (Fmax - Fmin);
      Fmin = Fmin - 0.1 * (Fmax - Fmin);
      if (Fmax == Fmin)
      {
         Fmax = Fmax + 0.1 * PABS(Fmax);
         Fmin = Fmin - 0.1 * PABS(Fmin);
      }
      for (ii = 0; ii <= ntimes; ii++)
      {
         Fcounts[ii] = new int[nbins];
         for (kk = 0; kk < nbins; kk++) Fcounts[ii][kk] = 0;
      }

      double mean=0, stdev=0, mean2=0, stdev2=0, d1, d2;
      double *samFuzzy = new double[ntimes*nInputs];
      double *samOutSave = new double[nSamp*ntimes];
      for (ss = 0; ss < nSamp; ss++) mean += samOutputs[ss];
      mean /= (double) nSamp;
      for (ss = 0; ss < nSamp; ss++) 
         stdev += pow(samOutputs[ss]-mean, 2.0);
      stdev = sqrt(stdev/(double) nSamp);
      printAsterisks(PL_INFO, 0);
      printf("Sample mean    = %e (RS uncertainties not included)\n",mean);
      printf("Sample std dev = %e (RS uncertainties not included)\n",stdev);
      printEquals(PL_INFO, 0);
      
      for (ss = 0; ss < nSamp; ss++)
      {
         if (samStds[ss] == 0)
            for (ii = 0; ii < ntimes; ii++) samFuzzy[ii] = samOutputs[ss];
         else
         {
            rsPDF = new PDFNormal(samOutputs[ss],samStds[ss]);
            d1 = samOutputs[ss] - 3 * samStds[ss];
            d2 = samOutputs[ss] + 3 * samStds[ss];
            rsPDF->genSample(ntimes,samFuzzy,&d1,&d2);
            delete rsPDF;
         }
         for (ii = 0; ii < ntimes; ii++) samOutSave[ss*ntimes+ii] = samFuzzy[ii];

         ddata = samOutputs[ss] - Fmin;
         if (Fmax > Fmin) ddata = ddata / ((Fmax - Fmin) / nbins);
         else             ddata = nbins / 2;
         kk = (int) ddata;
         if (kk < 0)      kk = 0;
         if (kk >= nbins) kk = nbins - 1;
         Fcounts[ntimes][kk]++;
         
         for (ii = 0; ii < ntimes; ii++) 
         {
            ddata = samFuzzy[ii] - Fmin;
            if (Fmax > Fmin)
                 ddata = ddata / ((Fmax - Fmin) / nbins);
            else ddata = nbins / 2;
            kk = (int) ddata;
            if (kk < 0)      kk = 0;
            if (kk >= nbins) kk = nbins - 1;
            Fcounts[ii][kk]++;
         }
         if (ss % (nSamp / 8) == 0)
         {
            printf(".");
            fflush(stdout);
         }
      }
      mean2 = 0.0;
      for (ss = 0; ss < nSamp*ntimes; ss++) mean2 += samOutSave[ss];
      mean2 /= (double) (nSamp*ntimes);
      stdev2 = 0.0;
      for (ss = 0; ss < nSamp*ntimes; ss++) 
         stdev2 += pow(samOutSave[ss] - mean2, 2.0);
      stdev2 = sqrt(stdev2/(double) (nSamp*ntimes));
      printf("Sample mean    = %e (RS uncertainties included)\n", mean2);
      printf("Sample std dev = %e (RS uncertainties included)\n", stdev2);
      printAsterisks(PL_INFO, 0);

      delete [] samStds;
      delete [] samInputs;
      delete [] samOutputs;
      delete [] samFuzzy;
      delete [] samOutSave;

      if (psPlotTool_ == 1)
      {
         if (!strcmp(command, "rsua2")) fp = fopen("scilabrsua2.sci", "w");
         else                           fp = fopen("scilabrsua.sci", "w");
         if (fp == NULL)
         {
            printf("rsua2 ERROR: cannot open scilab file.\n");
            for (ii = 0; ii <= ntimes; ii++) delete [] Fcounts[ii];
            delete [] Fcounts;
            return -1;
         }
      }
      else
      {
         if (!strcmp(command, "rsua2")) fp = fopen("matlabrsua2.m", "w");
         else                           fp = fopen("matlabrsua.m", "w");
         if (fp == NULL)
         {
            printf("rsua2 ERROR: cannot open matlab file.\n");
            for (ii = 0; ii <= ntimes; ii++) delete [] Fcounts[ii];
            delete [] Fcounts;
            return -1;
         }
      }
      fwriteHold(fp, 0);
      fprintf(fp, "subplot(2,2,1)\n");
      fprintf(fp, "X = [\n");
      for (kk = 0; kk < nbins; kk++)
         fprintf(fp, "%e\n", (Fmax-Fmin)/nbins*(0.5+kk)+Fmin);
      fprintf(fp, "];\n");
      for (ii = 0; ii <= ntimes; ii++)
      {
         fprintf(fp, "N%d = [\n", ii+1);
         for (kk = 0; kk < nbins; kk++)
            fprintf(fp, "%d\n",  Fcounts[ii][kk]);
         fprintf(fp, "];\n");
      }
      fprintf(fp, "N = [");
      for (ii = 0; ii <= ntimes; ii++)
         fprintf(fp, "N%d/sum(N%d) ", ii+1, ii+1);
      fprintf(fp, "];\n");
      fprintf(fp, "NA = N(:,%d+1);\n",ntimes);
      fprintf(fp, "NA = NA / sum(NA);\n");
      fprintf(fp, "NB = sum(N(:,1:%d)');\n",ntimes);
      fprintf(fp, "NB = NB' / sum(NB);\n");
      fprintf(fp, "NN = [NA NB];\n");
      fprintf(fp, "bar(X,NA,1.0)\n");
      fprintf(fp, "ymin = min(min(NA),min(NB));\n");
      fprintf(fp, "ymax = max(max(NA),max(NB));\n");
      fprintf(fp, "axis([min(X) max(X) ymin ymax])\n");
      fwritePlotAxes(fp);
      fwritePlotTitle(fp, "Prob. Dist. (means of RS)");
      fwritePlotXLabel(fp, "Output Value");
      fwritePlotYLabel(fp, "Probabilities)");
      if (psPlotTool_ == 0)
      {
         fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",mean);
         fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                 stdev);
      }
      if (faType == PSUADE_RS_MARS)
      {
         printf("Deterministic RS used ==> no RS uncertainties.\n");
      }
      else
      {
         fprintf(fp, "subplot(2,2,3)\n");
         fprintf(fp, "bar(X,NB,1.0)\n");
         fprintf(fp, "axis([min(X) max(X) ymin ymax])\n");
         fwritePlotAxes(fp);
         fwritePlotTitle(fp, "Prob. Dist. (RS with uncertainties)");
         fwritePlotXLabel(fp, "Output Value");
         fwritePlotYLabel(fp, "Probabilities)");
         if (psPlotTool_ == 0)
         {
            fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",
                    mean2);
            fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                    stdev2);
         }
      }
      for (ii = 0; ii <= ntimes; ii++)
      {
         fprintf(fp, "for ii = 2 : %d\n", nbins);
         fprintf(fp, "   N%d(ii) = N%d(ii) + N%d(ii-1);\n",ii+1,ii+1,ii+1);
         fprintf(fp, "end;\n");
      }
      fprintf(fp, "N = [");
      for (ii = 0; ii <= ntimes; ii++)
         fprintf(fp, "N%d/N%d(%d) ", ii+1, ii+1, nbins);
      fprintf(fp, "];\n");
      fprintf(fp, "subplot(2,2,[2 4])\n");
      fprintf(fp, "NA = N(:,%d+1);\n",ntimes);
      fprintf(fp, "NA = NA / NA(%d);\n",nbins);
      fprintf(fp, "NB = sum(N(:,1:%d)');\n",ntimes);
      fprintf(fp, "NB = NB' / NB(%d);\n", nbins);
      fprintf(fp, "NN = [NA NB];\n");
      if (faType == PSUADE_RS_MARS)
      {
         fprintf(fp, "plot(X,NA,'linewidth',3)\n");
         fwritePlotTitle(fp,"Cum. Dist.: (b) mean; (g) with uncertainties");
      }
      else 
      {
         fprintf(fp, "plot(X,NN,'linewidth',3)\n");
         fwritePlotTitle(fp,"Cum. Dist.: (*) uncertainties unavailable");
      }
      fwritePlotAxes(fp);
      fwritePlotXLabel(fp, "Output Value");
      fwritePlotYLabel(fp, "Probabilities");
      fclose(fp);
      if (!strcmp(command, "rsua2")) 
      {
         if (psPlotTool_ == 1)
              printf("Output distribution plots is in scilabrsua2.sci.\n");
         else printf("Output distribution plots is in matlabrsua2.m.\n");
      }
      else
      {
         if (psPlotTool_ == 1)
              printf("Output distribution plots is in scilabrsua.sci.\n");
         else printf("Output distribution plots is in matlabrsua.m.\n");
      }
      for (ii = 0; ii < ntimes; ii++) delete [] Fcounts[ii];
      delete [] Fcounts;
   }

   // +++ rsb_ua 
   else if (!strcmp(command, "rsb_ua"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rsb_ua: uncertainty analysis on response surface\n");
         printf("syntax: rsb_ua (no argument needed)\n");
         printf("This command perform uncertainty analysis on the response\n");
         printf("surface built from the LOADED sample. Uncertainty analysis\n");
         printf("is performed using a user-specified sample in PSUADE data\n");
         printf("format (created by running psuade on an input file). If you\n");
         printf("select a stochastic response surface type (Kriging, MARSB,\n");
         printf("or polynomial regression, the effect of response surface\n");
         printf("uncertainty will be shown on the PDF and CDF plots.)\n");
         printf("NOTE: This command is more general by allowing users to\n");
         printf("      provide the UA sample instead of generating it \n");
         printf("      internally.\n");
         printf("NOTE: This analysis will be replaced by rsua + rsuab.\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data (load sample first).\n");
         return -1;
      }
      int   discFile=0, nInps, nOuts, dnInps;
      char  discFileName[1001], uaFileName[1001];
      PsuadeData *discIO=NULL, *sampleIO=NULL;
      FuncApprox **faPtrsRsEval=NULL;
      
      sscanf(lineIn,"%s %s", command, winput);
      outputID = 0;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
      outputID = getInt(1, nOutputs, pString);
      outputID--;

      faPtrsRsEval = new FuncApprox*[2];
      faPtrsRsEval[0] = NULL;
      faPtrsRsEval[1] = NULL;
      printf("Use discrepancy model (in PSUADE data format)? ('y' or 'n') ");
      scanf("%s", winput);
      fgets(lineIn2, 500, stdin);
      if (winput[0] == 'y')
      {
         discFile = 1;
         printf("Enter discrepancy model file (in PSUADE data format): ");
         scanf("%s", discFileName);
         fgets(lineIn2, 500, stdin);
         discIO = new PsuadeData();
         status = discIO->readPsuadeFile(discFileName);
         if (status != 0)
         {
            printf("ERROR: cannot read discrepancy model file.\n");
            delete [] faPtrsRsEval;
            delete discIO;
            return -1;
         }
         discIO->getParameter("input_ninputs", pPtr);
         dnInps = pPtr.intData_;
         if (dnInps < nInputs)
         {
            printf("Discrepancy model has %d inputs. So the first\n", dnInps);
            printf("%d inputs in the model file will be assumed to\n", dnInps);
            printf("be associated with the inputs of the discrepancy model.\n");
         }
         discIO->getParameter("output_noutputs", pPtr);
         nOuts = pPtr.intData_;
         if (nOuts > 1)
         {
            printf("The discrepancy model has nOutputs > 1.\n");
            printf("This is currently not supported.\n");
            printf("Use 'odelete' to modify your discrepancy model file.\n");
            delete [] faPtrsRsEval;
            delete discIO;
            return -1;
         }
         printf("** CREATING RESPONSE SURFACE FOR DISCREPANCY MODEL\n");
         faPtrsRsEval[1] = genFAInteractive(discIO, 3);
         delete discIO;
         discIO = NULL;
      }

      printf("Enter UA sample file name (in PSUADE data format): ");
      scanf("%s", uaFileName);
      fgets(lineIn2, 500, stdin);
      sampleIO = new PsuadeData();
      status = sampleIO->readPsuadeFile(uaFileName);
      if (status != 0)
      {
         printf("ERROR: cannot read UA sample file.\n");
         delete [] faPtrsRsEval;
         delete sampleIO;
         return -1;
      }
      sampleIO->getParameter("input_ninputs", pPtr);
      nInps = pPtr.intData_;
      if (nInps != nInputs)
      {
         printf("ERROR: UA nInputs not equal to nInputs in local memory.\n");
         printf(":      input size in local memory = %d.\n",nInputs);
         printf(":      input size from UA file    = %d.\n",nInps);
         if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
         delete [] faPtrsRsEval;
         delete sampleIO;
         return -1;
      }
      sampleIO->getParameter("method_nsamples", pPtr);
      kk = pPtr.intData_;
      if (kk < 1000)
      {
         printf("ERROR: Your sample for UA should be at least 1000\n");
         printf("       to give any reasonable UA results.\n");
         if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
         delete [] faPtrsRsEval;
         delete sampleIO;
         return -1;
      }

      int uaMethod=0;
      int includeRSErr=0, numBS=1;
      printf("Include response surface uncertainties in UA? (y or n) ");
      scanf("%s", winput);
      fgets(lineIn2, 500, stdin);
      if (winput[0] == 'y')
      {
         includeRSErr = 1;
         printf("Three options are available for including RS uncertainties:\n");
         printf("1. use bootstrapping + RS (or deterministic RS, e.g. MARS)\n");
         printf("2. use stochastic RS (Kriging, MARS+Bootstrap, regression)\n");
         printf("3. use (2) but perform worst-case analysis (2 - average case)\n");
         sprintf(pString, "Select 1, 2, or 3 : ");
         uaMethod = getInt(1, 3, pString);
         if (uaMethod == 1)
         {
            sprintf(pString, "How many bootstrapped samples to use (10 - 300) : ");
            numBS = getInt(10, 300, pString);
         }
      }

      // perform UA
      int    userNSams;
      double *userSamInps, *userSamOuts, *userSamStds;
      sampleIO->getParameter("method_nsamples", pPtr);
      userNSams = pPtr.intData_;
      sampleIO->getParameter("input_sample", pPtr);
      userSamInps = pPtr.dbleArray_;
      pPtr.dbleArray_ = NULL;
      userSamOuts = new double[userNSams];
      userSamStds = new double[userNSams];

      // use deterministic 
      if (uaMethod == 0)
      {
         printf("** CREATING RESPONSE SURFACE FOR PRIMARY MODEL\n");
         faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples);
         if (faPtrsRsEval[0] == NULL)
         {
            printf("ERROR: cannot generate response surface.\n");
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            delete [] userSamInps;
            delete [] userSamOuts;
            delete [] userSamStds;
            delete sampleIO;
            return -1;
         }
         faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
         faPtrsRsEval[0]->setOutputLevel(0);
         status = faPtrsRsEval[0]->initialize(sampleInputs,sampleOutputs);
         if (status != 0)
         {
            printf("ERROR: cannot initialize response surface.\n");
            if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            delete [] userSamInps;
            delete [] userSamOuts;
            delete [] userSamStds;
            delete sampleIO;
            return -1;
         }
         faPtrsRsEval[0]->evaluatePoint(userNSams,userSamInps,userSamOuts);
         if (discFile == 1)
         {
            for (ss = 0; ss < userNSams; ss++)
            {
               ddata = faPtrsRsEval[1]->evaluatePoint(&userSamInps[ss*nInputs]);
               userSamOuts[ss] += ddata;
            }
         }
         
         double mean=0, stdev=0;
         for (ss = 0; ss < userNSams; ss++) mean += userSamOuts[ss];
         mean /= (double) userNSams;
         for (ss = 0; ss < userNSams; ss++)
            stdev += pow(userSamOuts[ss]-mean, 2.0);
         stdev = sqrt(stdev/(double) userNSams);
         printAsterisks(PL_INFO, 0);
         printf("Sample mean    = %e (without RS uncertainties)\n", mean);
         printf("Sample std dev = %e (without RS uncertainties)\n", stdev);
         printEquals(PL_INFO, 0);

         fp = NULL;
         fp = fopen("matlabrsbua.m", "w");
         if (fp != NULL)
         {
            fprintf(fp,"Y = [\n");
            for (ss = 0; ss < userNSams; ss++) 
               fprintf(fp,"%e\n",userSamOuts[ss]);
            fprintf(fp, "];\n");
            fwriteHold(fp, 0);
            fprintf(fp,"subplot(1,2,1);\n");
            fprintf(fp,"[nk,xk] = hist(Y,50);\n");
            fprintf(fp,"nk = nk / %d;\n", userNSams);
            fprintf(fp,"bar(xk,nk,1.0);\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Probability Distribution");
            fwritePlotXLabel(fp,"Output Value");
            fwritePlotYLabel(fp,"Probabilities");
            fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",
                    mean);
            fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                    stdev);
            fprintf(fp,"subplot(1,2,2);\n");
            fprintf(fp, "Y = sort(Y);\n");
            fprintf(fp, "X = (1 : %d)' / %d;\n", userNSams, userNSams);
            fprintf(fp,"plot(Y, X, 'lineWidth',3)\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Cumulative Distribution");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Cum. Prob.");
            fclose(fp);
            printf("Output distribution plots are in matlabrsbua.m.\n");
         }
      }

      // bootstrapped method
      if (uaMethod == 1)
      {
         int bsnSams, rsMethod;
         printf("** CREATING RESPONSE SURFACE FOR PRIMARY MODEL\n");
         faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples);
         if (faPtrsRsEval[0] == NULL)
         {
            printf("ERROR: cannot generate primary response surface.\n");
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            delete sampleIO;
            return -1;
         }
         rsMethod = faPtrsRsEval[0]->getID(); 
         delete faPtrsRsEval[0];
         faPtrsRsEval[0] = NULL;
         
         int    its, *useFlags = new int[nSamples];
         double *bsSamInps = new double[nSamples*nInputs];
         double *bsSamOuts = new double[nSamples];
         double *bsMeans = new double[numBS];
         double *bsStds  = new double[numBS];

         fp = NULL;
         fp = fopen("matlabrsbua.m", "w");
         for (its = 0; its < numBS; its++)
         {
            for (ss = 0; ss < nSamples; ss++) useFlags[ss] = 0;
            bsnSams = 0;
            for (ss = 0; ss < nSamples; ss++) 
            {
               jj = PSUADE_rand() % nSamples;
               if (useFlags[jj] == 0)
               {
                  for (ii = 0; ii < nInputs; ii++)
                     bsSamInps[bsnSams*nInputs+ii] = sampleInputs[jj*nInputs+ii];
                  bsSamOuts[bsnSams] = sampleOutputs[jj*nOutputs+outputID];
                  useFlags[jj] = 1;
                  bsnSams++;
               }
            }
            printf("Bootstrap %d has sample size = %d (%d)\n",its+1,bsnSams,
                   nSamples);
            faPtrsRsEval[0] = genFA(rsMethod, nInputs, -1, bsnSams);
            faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
            faPtrsRsEval[0]->setOutputLevel(0);
            status = faPtrsRsEval[0]->initialize(bsSamInps,bsSamOuts);
            if (status != 0)
            {
               printf("ERROR: in initializing response surface (1).\n");
               if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
               if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
               delete [] faPtrsRsEval;
               delete [] bsSamInps;
               delete [] bsSamOuts;
               delete [] userSamInps;
               delete [] userSamOuts;
               delete [] userSamStds;
               delete [] useFlags;
               delete [] bsMeans;
               delete [] bsStds;
               delete sampleIO;
               return -1;
            } 
            faPtrsRsEval[0]->evaluatePoint(userNSams,userSamInps,userSamOuts);
            delete faPtrsRsEval[0];
            faPtrsRsEval[0] = NULL;
            if (discFile == 1)
            {
               for (ss = 0; ss < userNSams; ss++)
               {
                  ddata = faPtrsRsEval[1]->evaluatePoint(&userSamInps[ss*nInputs]);
                  userSamOuts[ss] += ddata;
               }
            }
            bsMeans[its] = bsStds[its] = 0.0;
            for (ss = 0; ss < userNSams; ss++) bsMeans[its] += userSamOuts[ss];
            bsMeans[its] /= (double) userNSams;
            for (ss = 0; ss < userNSams; ss++)
               bsStds[its] += pow(userSamOuts[ss] - bsMeans[its], 2.0);
            bsStds[its] = sqrt(bsStds[its] / userNSams);
            if (fp != NULL)
            {
               fprintf(fp, "%% bootstrapped samples\n");
               fprintf(fp, "Y = [\n");
               for (ss = 0; ss < userNSams; ss++) 
                  fprintf(fp,"%e\n",userSamOuts[ss]);
               fprintf(fp, "];\n");
               fprintf(fp, "Y%d = sort(Y);\n",its+1);
               fprintf(fp, "X%d = (1 : %d)';\n", its+1, userNSams);
               fprintf(fp, "X%d = X%d / %d;\n", its+1, its+1, userNSams);
               if (its == 0)
               {
                  fprintf(fp, "YY = Y%d;\n", its+1);
                  fprintf(fp, "XX = X%d;\n", its+1);
               }
               else
               {
                  fprintf(fp, "YY = [YY Y%d];\n", its+1);
                  fprintf(fp, "XX = [XX X%d];\n", its+1);
               }
            }
         }
         printAsterisks(PL_INFO, 0);
         double mean, stdev;
         mean = stdev = 0.0;
         for (its = 0; its < numBS; its++) mean += bsMeans[its];
         mean /= (double) numBS;
         for (ss = 0; ss < numBS; ss++) stdev += pow(bsMeans[ss]-mean, 2.0);
         stdev = sqrt(stdev/(double) numBS);
         printf("Sample mean    = %e (std = %e)\n", mean, stdev);
         mean = stdev = 0.0;
         for (its = 0; its < numBS; its++) mean += bsStds[its];
         mean /= (double) numBS;
         for (ss = 0; ss < numBS; ss++) stdev += pow(bsStds[ss]-mean, 2.0);
         stdev = sqrt(stdev/(double) numBS);
         printf("Sample std dev = %e (std = %e)\n", mean, stdev);
         printEquals(PL_INFO, 0);
         if (fp != NULL)
         {
            fwriteHold(fp, 0);
            fprintf(fp,"subplot(1,2,1);\n");
            fprintf(fp,"[nk,xk] = hist(YY,50);\n");
            fprintf(fp,"plot(xk,nk, 'lineWidth',2)\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Probability Distribution");
            fwritePlotXLabel(fp,"Output Value");
            fwritePlotYLabel(fp,"Probabilities");
            fprintf(fp,"subplot(1,2,2);\n");
            fprintf(fp,"plot(YY, XX, 'lineWidth',3)\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Cumulative Distribution");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
            fclose(fp);
            printf("Output distribution plots are in matlabrsbua.m.\n");
         }
         delete [] bsSamInps;
         delete [] bsSamOuts;
         delete [] useFlags;
         delete [] bsMeans;
         delete [] bsStds;
      }
      
      // use stochastic response surface with average case analysis
      if (uaMethod == 2)
      {
         printf("** CREATING RESPONSE SURFACE FOR PRIMARY MODEL\n");
         faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples);
         if (faPtrsRsEval[0] == NULL)
         {
            printf("ERROR: cannot generate response surface.\n");
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            delete [] userSamInps;
            delete [] userSamOuts;
            delete [] userSamStds;
            delete sampleIO;
            return -1;
         }
         faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
         faPtrsRsEval[0]->setOutputLevel(0);
         status = faPtrsRsEval[0]->initialize(sampleInputs,sampleOutputs);
         if (status != 0)
         {
            printf("ERROR: cannot initialize response surface.\n");
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            delete [] userSamInps;
            delete [] userSamOuts;
            delete [] userSamStds;
            delete sampleIO;
            return -1;
         }
         faPtrsRsEval[0]->evaluatePointFuzzy(userNSams,userSamInps,
                                             userSamOuts,userSamStds);
         double discSamStd;
         if (discFile == 1)
         {
            for (ss = 0; ss < userNSams; ss++)
            {
               userSamOuts[ss] += 
                 faPtrsRsEval[1]->evaluatePointFuzzy(&userSamInps[ss*nInputs],
                                                     discSamStd);
               ddata = pow(userSamStds[ss],2.0) + discSamStd * discSamStd;
               userSamStds[ss] = sqrt(ddata);
            }
         }

         fp = fopen("rsbua_sample","w");
         if (fp != NULL)
         {
            fprintf(fp,"%% This file is primarily for diagnostics and\n");
            fprintf(fp,"%% expert analysis\n");
            fprintf(fp,"%% First line: nSamples nInputs\n");
            fprintf(fp,"%% All inputs, output, output-3*sigma, output+3*sigma\n");
            fprintf(fp,"%d %d 3\n", userNSams, nInputs);
            for (ss = 0; ss < userNSams; ss++)
            {
               for (ii = 0; ii < nInputs; ii++)
                  fprintf(fp, "%e ", userSamInps[ss*nInputs+ii]);
               fprintf(fp, "%e ", userSamOuts[ss]);
               fprintf(fp, "%e ", userSamOuts[ss]-3*userSamStds[ss]);
               fprintf(fp, "%e\n", userSamOuts[ss]+3*userSamStds[ss]);
            }
            fclose(fp);
            printf("The outputs and stds of your sample has been written");
            printf(" to 'rsbua_sample'.\n");
         }

         double mean=0, stdev=0;
         for (ss = 0; ss < userNSams; ss++) mean += userSamOuts[ss];
         mean /= (double) userNSams;
         for (ss = 0; ss < userNSams; ss++)
            stdev += pow(userSamOuts[ss]-mean, 2.0);
         stdev = sqrt(stdev/(double) userNSams);
         printAsterisks(PL_INFO, 0);
         printf("Sample mean    = %e (RS uncertainties not included)\n", mean);
         printf("Sample std dev = %e (RS uncertainties not included)\n", stdev);
         printEquals(PL_INFO, 0);

         int    nbins = 100, ntimes=20;
         int    **Fcounts = new int*[ntimes+1];
         double Fmax=-PSUADE_UNDEFINED;
         double Fmin=PSUADE_UNDEFINED;
         PDFNormal *rsPDF=NULL;
         for (ss = 0; ss < userNSams; ss++)
         {
            if (userSamOuts[ss]+3*userSamStds[ss] > Fmax)
               Fmax = userSamOuts[ss] + 3 * userSamStds[ss];
            if (userSamOuts[ss]-3*userSamStds[ss] < Fmin)
               Fmin = userSamOuts[ss] - 3 * userSamStds[ss];
         }
         Fmax = Fmax + 0.1 * (Fmax - Fmin);
         Fmin = Fmin - 0.1 * (Fmax - Fmin);
         if (Fmax == Fmin)
         {
            Fmax = Fmax + 0.1 * PABS(Fmax);
            Fmin = Fmin - 0.1 * PABS(Fmin);
         }
         for (ii = 0; ii <= ntimes; ii++)
         {
            Fcounts[ii] = new int[nbins];
            for (kk = 0; kk < nbins; kk++) Fcounts[ii][kk] = 0;
         }

         double *samOutTime = new double[ntimes*nInputs];
         double *samOutSave = new double[userNSams*ntimes], d1, d2;
         for (ss = 0; ss < userNSams; ss++)
         {
            if (userSamStds[ss] == 0)
               for (ii = 0; ii < ntimes; ii++) samOutTime[ii] = userSamOuts[ss];
            else
            {
               rsPDF = new PDFNormal(userSamOuts[ss],userSamStds[ss]);
               d1 = userSamOuts[ss] - 3.0 * userSamStds[ss];
               d2 = userSamOuts[ss] + 3.0 * userSamStds[ss];
               rsPDF->genSample(ntimes,samOutTime,&d1,&d2);
               delete rsPDF;
            }
            for (ii = 0; ii < ntimes; ii++) 
               samOutSave[ss*ntimes+ii] = samOutTime[ii];

            ddata = userSamOuts[ss] - Fmin;
            if (Fmax > Fmin) ddata = ddata / ((Fmax - Fmin) / nbins);
            else             ddata = nbins / 2;
            kk = (int) ddata;
            if (kk < 0)      kk = 0;
            if (kk >= nbins) kk = nbins - 1;
            Fcounts[ntimes][kk]++;

            for (ii = 0; ii < ntimes; ii++)
            {
               ddata = samOutTime[ii] - Fmin;
               if (Fmax > Fmin)
                    ddata = ddata / ((Fmax - Fmin) / nbins);
               else ddata = nbins / 2;
               kk = (int) ddata;
               if (kk < 0)      kk = 0;
               if (kk >= nbins) kk = nbins - 1;
               Fcounts[ii][kk]++;
            }
         }
         delete [] samOutTime;
         double mean2=0, stdev2=0;
         for (ss = 0; ss < userNSams*ntimes; ss++) mean2 += samOutSave[ss];
         mean2 /= (double) (userNSams*ntimes);
         stdev2 = 0.0;
         for (ss = 0; ss < userNSams*ntimes; ss++)
            stdev2 += pow(samOutSave[ss] - mean2, 2.0);
         stdev2 = sqrt(stdev2/(double) (userNSams*ntimes));
         printf("Sample mean    = %e (RS uncertainties included)\n", mean2);
         printf("Sample std dev = %e (RS uncertainties included)\n", stdev2);
         printAsterisks(PL_INFO, 0);
         delete [] samOutSave;

         fp = fopen("matlabrsbua.m", "w");
         if (fp == NULL)
         {
            printf("INFO: cannot write the PDFs/CDFs to matlab file.\n");
         }
         else
         {
            fwriteHold(fp, 0);
            fprintf(fp, "subplot(2,2,1)\n");
            fprintf(fp, "X = [\n");
            for (kk = 0; kk < nbins; kk++)
               fprintf(fp, "%e\n",(Fmax-Fmin)/nbins*(0.5+kk)+Fmin);
            fprintf(fp, "];\n");
            for (ii = 0; ii <= ntimes; ii++)
            {
               fprintf(fp, "N%d = [\n", ii+1);
               for (kk = 0; kk < nbins; kk++)
                  fprintf(fp, "%d\n",  Fcounts[ii][kk]);
               fprintf(fp, "];\n");
            }
            fprintf(fp, "N = [");
            for (ii = 0; ii <= ntimes; ii++)
               fprintf(fp, "N%d/sum(N%d) ", ii+1, ii+1);
            fprintf(fp, "];\n");
            fprintf(fp, "NA = N(:,%d+1);\n",ntimes);
            fprintf(fp, "NA = NA / sum(NA);\n");
            fprintf(fp, "NB = sum(N(:,1:%d)');\n",ntimes);
            fprintf(fp, "NB = NB' / sum(NB);\n");
            fprintf(fp, "NN = [NA NB];\n");
            fprintf(fp, "bar(X,NA,1.0)\n");
            fprintf(fp, "ymin = min(min(NA),min(NB));\n");
            fprintf(fp, "ymax = max(max(NA),max(NB));\n");
            fprintf(fp, "axis([min(X) max(X) ymin ymax])\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Prob. Dist. (means of RS)");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities)");
            fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",
                    mean);
            fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                    stdev);
            if (faType == PSUADE_RS_MARS)
            {
               printf("Deterministic RS used ==> no RS uncertainties.\n");
            }
            else
            {
               fprintf(fp,"subplot(2,2,3)\n");
               fprintf(fp,"bar(X,NB,1.0)\n");
               fprintf(fp,"axis([min(X) max(X) ymin ymax])\n");
               fwritePlotAxes(fp);
               fwritePlotTitle(fp,"Prob. Dist. (RS with uncertainties)");
               fwritePlotXLabel(fp,"Output Value");
               fwritePlotYLabel(fp,"Probabilities)");
               fprintf(fp,"text(0.05,0.9,'Mean = %12.4e','sc','FontSize',11)\n",
                       mean2);
               fprintf(fp,"text(0.05,0.85,'Std  = %12.4e','sc','FontSize',11)\n",
                       stdev2);
            }
            for (ii = 0; ii <= ntimes; ii++)
            {
               fprintf(fp,"for ii = 2 : %d\n", nbins);
               fprintf(fp,"  N%d(ii) = N%d(ii) + N%d(ii-1);\n",ii+1,ii+1,ii+1);
               fprintf(fp,"end;\n");
            }
            fprintf(fp, "N = [");
            for (ii = 0; ii <= ntimes; ii++)
               fprintf(fp,"N%d/N%d(%d) ", ii+1, ii+1, nbins);
            fprintf(fp, "];\n");
            fprintf(fp, "subplot(2,2,[2 4])\n");
            fprintf(fp, "NA = N(:,%d+1);\n",ntimes);
            fprintf(fp, "NA = NA / NA(%d);\n",nbins);
            fprintf(fp, "NB = sum(N(:,1:%d)');\n",ntimes);
            fprintf(fp, "NB = NB' / NB(%d);\n", nbins);
            fprintf(fp, "NN = [NA NB];\n");
            if (faType == PSUADE_RS_MARS)
            {
               fprintf(fp, "plot(X,NA,'linewidth',3)\n");
               fwritePlotTitle(fp,"Cum. Dist.: (b) mean; (g) with uncertainties");
            }
            else
            {
               fprintf(fp, "plot(X,NN,'linewidth',3)\n");
               fwritePlotTitle(fp,"Cum. Dist.: (*) uncertainties unavailable");
            }
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
            fclose(fp);
            printf("Output distribution plots are in matlabrsbua.m.\n");
         }
         for (ii = 0; ii <= ntimes; ii++) delete [] Fcounts[ii];
         delete [] Fcounts;
      }

      // use stochastic response surface with worst case analysis
      if (uaMethod == 3)
      {
         printf("** CREATING RESPONSE SURFACE FOR PRIMARY MODEL\n");
         faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples);
         if (faPtrsRsEval[0] == NULL)
         {
            printf("ERROR: cannot generate response surface.\n");
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            delete [] userSamInps;
            delete [] userSamOuts;
            delete [] userSamStds;
            delete sampleIO;
            return -1;
         }
         faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
         faPtrsRsEval[0]->setOutputLevel(0);
         status = faPtrsRsEval[0]->initialize(sampleInputs,sampleOutputs);
         if (status != 0)
         {
            printf("ERROR: cannot initialize response surface.\n");
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            delete [] userSamInps;
            delete [] userSamOuts;
            delete [] userSamStds;
            delete sampleIO;
            return -1;
         }
         
         faPtrsRsEval[0]->evaluatePointFuzzy(userNSams,userSamInps,
                                             userSamOuts,userSamStds);
         
         double discSamStd;
         if (discFile == 1)
         {
            for (ss = 0; ss < userNSams; ss++)
            {
               userSamOuts[ss] += 
                 faPtrsRsEval[1]->evaluatePointFuzzy(&userSamInps[ss*nInputs],
                                                     discSamStd);
               ddata = pow(userSamStds[ss],2.0) + discSamStd * discSamStd;
               userSamStds[ss] = sqrt(ddata);
            }
         }
         
         double mean=0, stdev=0;
         for (ss = 0; ss < userNSams; ss++) mean += userSamOuts[ss];
         mean /= (double) userNSams;
         for (ss = 0; ss < userNSams; ss++)
            stdev += pow(userSamOuts[ss]-mean, 2.0);
         stdev = sqrt(stdev/(double) userNSams);
         printAsterisks(PL_INFO, 0);
         printf("Sample mean    = %e (RS uncertainties not included)\n",mean);
         printf("Sample std dev = %e (RS uncertainties not included)\n",stdev);
         printEquals(PL_INFO, 0);

         fp = fopen("rsbua_sample","w");
         fprintf(fp, "%% This file is primarily for diagnostics and \n");
         fprintf(fp, "%% expert analysis\n");
         fprintf(fp, "%% First line: nSamples nInputs\n");
         fprintf(fp, "%% All inputs, output, output-3*sigma, output+3*sigma\n");
         fprintf(fp, "%d %d 3\n", userNSams, nInputs);
         for (ss = 0; ss < userNSams; ss++)
         {
            for (ii = 0; ii < nInputs; ii++)
               fprintf(fp, "%e ", userSamInps[ss*nInputs+ii]);
            fprintf(fp, "%e ", userSamOuts[ss]);
            fprintf(fp, "%e ", userSamOuts[ss]-3*userSamStds[ss]);
            fprintf(fp, "%e\n", userSamOuts[ss]+3*userSamStds[ss]);
         }
         fclose(fp);
         
         int    nbins = 100, ntimes=7;
         int    **Fcounts = new int*[ntimes+1];
         double Fmax=-PSUADE_UNDEFINED;
         double Fmin=PSUADE_UNDEFINED;
         PDFNormal *rsPDF=NULL;
         for (ss = 0; ss < userNSams; ss++)
         {
            if (userSamOuts[ss]+3*userSamStds[ss] > Fmax)
               Fmax = userSamOuts[ss] + 3 * userSamStds[ss];
            if (userSamOuts[ss]-3*userSamStds[ss] < Fmin)
               Fmin = userSamOuts[ss] - 3 * userSamStds[ss];
         }
         Fmax = Fmax + 0.1 * (Fmax - Fmin);
         Fmin = Fmin - 0.1 * (Fmax - Fmin);
         if (Fmax == Fmin)
         {
            Fmax = Fmax + 0.1 * PABS(Fmax);
            Fmin = Fmin - 0.1 * PABS(Fmin);
         }
         for (ii = 0; ii <= ntimes; ii++)
         {
            Fcounts[ii] = new int[nbins];
            for (kk = 0; kk < nbins; kk++) Fcounts[ii][kk] = 0;
         }

         for (ss = 0; ss < userNSams; ss++)
         {
            for (ii = 0; ii < ntimes; ii++)
            {
               ddata = userSamOuts[ss]+userSamStds[ss]*(ii-3) - Fmin;
               if (Fmax > Fmin)
                    ddata = ddata / ((Fmax - Fmin) / nbins);
               else ddata = nbins / 2;
               kk = (int) ddata;
               if (kk < 0)      kk = 0;
               if (kk >= nbins) kk = nbins - 1;
               Fcounts[ii][kk]++;
            }
         }

         fp = fopen("matlabrsbua.m", "w");
         if (fp == NULL)
         {
            printf("INFO: cannot write the PDFs/CDFs to matlab file.\n");
         }
         else
         {
            fwriteHold(fp, 0);
            fprintf(fp, "%% worst case analysis\n");
            fprintf(fp, "X = [\n");
            for (kk = 0; kk < nbins; kk++)
               fprintf(fp, "%e\n", (Fmax-Fmin)/nbins*(0.5+kk)+Fmin);
            fprintf(fp, "];\n");
            for (ii = 0; ii < ntimes; ii++)
            {
               fprintf(fp, "E%d = [\n", ii+1);
               for (kk = 0; kk < nbins; kk++) 
                  fprintf(fp, "%d\n",Fcounts[ii][kk]);
               fprintf(fp, "];\n");
            }
            fprintf(fp, "EE = [");
            for (ii = 0; ii < ntimes; ii++)
               fprintf(fp, "E%d/sum(E%d) ", ii+1, ii+1);
            fprintf(fp, "];\n");
            fprintf(fp, "subplot(1,2,1)\n");
            fprintf(fp, "plot(X,EE,'lineWidth',2)\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Prob. Dist. (-3,2,1,0,1,2,3 std.)");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
            fprintf(fp, "subplot(1,2,2)\n");
            for (ii = 0; ii < ntimes; ii++)
            {
               fprintf(fp, "for ii = 2 : %d\n", nbins);
               fprintf(fp, "   E%d(ii) = E%d(ii) + E%d(ii-1);\n",ii+1,ii+1,ii+1);
               fprintf(fp, "end;\n");
            }
            fprintf(fp, "EE = [");
            for (ii = 0; ii < ntimes; ii++)
               fprintf(fp, "E%d/E%d(%d) ", ii+1, ii+1, nbins);
            fprintf(fp, "];\n");
            fprintf(fp, "plot(X,EE,'linewidth',2)\n");
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Cum. Dist. (-3,2,1,0,1,2,3 std.)");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
            fclose(fp);
            printf("Output distribution plots are in matlabrsbua.m.\n");
            for (ii = 0; ii < ntimes; ii++) delete [] Fcounts[ii];
            delete [] Fcounts;
         }
      }
      if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
      if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
      delete [] faPtrsRsEval;
      if (sampleIO != NULL) delete sampleIO;
      delete [] userSamInps;
      delete [] userSamOuts;
      delete [] userSamStds;
   }

   // +++ rs_ua2 
   else if (!strcmp(command, "rs_ua2"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs_ua2: uncertainty analysis on response surface (worst case)\n");
       printf("syntax: rs_ua2 (no argument needed)\n");
       printf("This command perform uncertainty analysis on the response\n");
       printf("surface built from the loaded sample. If you select a\n");
       printf("stochastic response surface type (Kriging, MARSB, or\n");
       printf("polynomial regression, the effect of response surface\n");
       printf("uncertainty will be shown on the PDF and CDF plots.\n");
       printf("This is a worst case analysis in the sense that the each\n");
       printf("histogram is constructed from perturbing each sample point\n");
       printf("with the same fraction of its standard deviation.\n");
       printf("NOTE: This analysis supports non-uniform distributions\n");
       printf("      for the inputs. Simply prescribe the distributions in\n");
       printf("      the data file and turn on use_input_pdfs in ANALYSIS.\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }
     Sampling   *samPtr;
     FuncApprox *faPtr;
     PDFManager *pdfman;
     psVector   vecOut, vecLower, vecUpper;
     psuadeIO->getParameter("ana_rstype", pPtr);
     faType = pPtr.intData_;
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;
     sprintf(pString,
            "Sample size for generating distribution? (10000 - 100000) ");
     int nSamp = getInt(10000, 100000, pString);
     flag = 0;
     printf("Save the generated sample in a file? (y or n) ");
     fgets(winput,10,stdin); 
     if (winput[0] == 'y') flag = 1; 

     printf("Phase 1 of 4: create response surface\n");
     printf("NOTE: the response surface type is taken from your data file.\n");
     faPtr = genFA(faType, nInputs, -1, nSamples);
     faPtr->setNPtsPerDim(32);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);
     double *tempY;
     if (nOutputs > 1)
     {
       tempY = new double[nSamples];
       for (ss = 0; ss < nSamples; ss++) 
         tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
     }
     else tempY = sampleOutputs;
     faPtr->initialize(sampleInputs,tempY);
     if (nOutputs > 1) delete [] tempY;
     tempY = NULL;
            
     printEquals(PL_INFO, 0);
     printf("Phase 2 of 4: create MC sample\n");
     double *samInputs  = new double[nInputs*nSamp];
     double *samOutputs = new double[nSamp];
     double *samStds    = new double[nSamp];
     psuadeIO->getParameter("ana_use_input_pdfs", pPtr);
     int usePDFs = pPtr.intData_;
     if (usePDFs == 1)
     {
       printf("NOTE: Some inputs have non-uniform PDFs.\n");
       printf("      A MC sample will be created with these PDFs.\n");
       psuadeIO->getParameter("method_sampling", pPtr);
       kk = pPtr.intData_;
       psuadeIO->updateMethodSection(PSUADE_SAMP_MC,-1,-1,-1,-1);
       pdfman = new PDFManager();
       pdfman->initialize(psuadeIO);
       vecOut.setLength(nSamp*nInputs);
       vecUpper.load(nInputs, iUpperB);
       vecLower.load(nInputs, iLowerB);
       pdfman->genSample(nSamp, vecOut, vecLower, vecUpper);
       for (ii = 0; ii < nSamp*nInputs; ii++) samInputs[ii] = vecOut[ii];
       psuadeIO->updateMethodSection(kk,-1,-1,-1,-1);
     }
     else
     {
       printf("NOTE: Uniform distributions will be used for all inputs.\n");
       printf("      To use other than uniform distributions, prescribe\n");
       printf("      them in the data file and set use_input_pdfs in the\n");
       printf("      ANALYSIS section.\n");
       int *samStates  = new int[nSamp];
       if (nInputs < 51)
            samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
       else samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
       samPtr->setPrintLevel(0);
       samPtr->setInputBounds(nInputs, iLowerB, iUpperB);
       samPtr->setOutputParams(1);
       samPtr->setSamplingParams(nSamp, -1, -1);
       samPtr->initialize(0);
       samPtr->getSamples(nSamp,nInputs,1,samInputs,samOutputs,samStates);
       delete samPtr;
       samPtr = NULL;
       delete [] samStates;
       samStates = NULL;
     }

     printf("Phase 3 of 4: evaluate sample\n");
     faPtr->evaluatePointFuzzy(nSamp,samInputs,samOutputs,samStds); 
     delete faPtr;
     faPtr = NULL;
     if (flag == 1)
     {
       fp = fopen("rsua2_sample","w");
       fprintf(fp, "%% inputs, output, output-3 sigma, output+3sigma\n");
       fprintf(fp, "%d %d 3\n", nSamp, nInputs);
       for (ss = 0; ss < nSamp; ss++)
       {
         for (ii = 0; ii < nInputs; ii++) 
           fprintf(fp, "%e ", samInputs[ss*nInputs+ii]);
         fprintf(fp, "%e ", samOutputs[ss]);
         fprintf(fp, "%e ", samOutputs[ss]-3*samStds[ss]);
         fprintf(fp, "%e\n", samOutputs[ss]+3*samStds[ss]);
       }
       fclose(fp);
       printf("A MC sample has been written to the file 'rsua2_sample'.\n");
     }

     double mean=0, stdev=0;
     for (ss = 0; ss < nSamp; ss++) mean += samOutputs[ss];
     mean /= (double) nSamp;
     for (ss = 0; ss < nSamp; ss++) 
       stdev += pow(samOutputs[ss]-mean, 2.0);
     stdev = sqrt(stdev/(double) nSamp);
     printAsterisks(PL_INFO, 0);
     printf("Sample mean    = %e (RS uncertainties not included)\n", mean);
     printf("Sample std dev = %e (RS uncertainties not included)\n", stdev);
     printEquals(PL_INFO, 0);
     printf("Phase 4 of 4: binning\n");

     int    nbins = 100, ntimes=7;
     int    **Fcounts = new int*[ntimes];
     double Fmax=-PSUADE_UNDEFINED;
     double Fmin=PSUADE_UNDEFINED;

     for (ss = 0; ss < nSamp; ss++)
     {
       if (samOutputs[ss] > Fmax) Fmax = samOutputs[ss];
       if (samOutputs[ss] < Fmin) Fmin = samOutputs[ss];
       if (samOutputs[ss]+3*samStds[ss] > Fmax) 
         Fmax = samOutputs[ss] + 3 * samStds[ss];
       if (samOutputs[ss]-3*samStds[ss] < Fmin) 
         Fmin = samOutputs[ss] - 3 * samStds[ss];
     }
     Fmax = Fmax + 0.1 * (Fmax - Fmin);
     Fmin = Fmin - 0.1 * (Fmax - Fmin);
     if (Fmax == Fmin)
     {
       Fmax = Fmax + 0.1 * PABS(Fmax);
       Fmin = Fmin - 0.1 * PABS(Fmin);
     }
     for (ii = 0; ii < ntimes; ii++)
     {
       Fcounts[ii] = new int[nbins];
       for (kk = 0; kk < nbins; kk++) Fcounts[ii][kk] = 0;
     }

     for (ss = 0; ss < nSamp; ss++)
     {
       for (ii = 0; ii < ntimes; ii++) 
       {
         ddata = samOutputs[ss]+samStds[ss]*(ii-3) - Fmin;
         if (Fmax > Fmin)
              ddata = ddata / ((Fmax - Fmin) / nbins);
         else ddata = nbins / 2;
         kk = (int) ddata;
         if (kk < 0)      kk = 0;
         if (kk >= nbins) kk = nbins - 1;
         Fcounts[ii][kk]++;
       }
     }

     delete [] samStds;
     delete [] samInputs;
     delete [] samOutputs;

     if (psPlotTool_ == 1)
     {
       fp = fopen("scilabrsua2.sci", "w");
       if (fp == NULL)
       {
          printf("rs_ua2 ERROR: cannot open scilab file.\n");
          for (ii = 0; ii < ntimes; ii++) delete [] Fcounts[ii];
          delete [] Fcounts;
          return -1;
       }
     }
     else
     {
       fp = fopen("matlabrsua2.m", "w");
       if (fp == NULL)
       {
         printf("rs_ua2 ERROR: cannot open matlab file.\n");
         for (ii = 0; ii < ntimes; ii++) delete [] Fcounts[ii];
         delete [] Fcounts;
         return -1;
       }
     }
     fprintf(fp, "X = [\n");
     for (kk = 0; kk < nbins; kk++)
       fprintf(fp, "%e\n", (Fmax-Fmin)/nbins*(0.5+kk)+Fmin);
     fprintf(fp, "];\n");
     fwriteHold(fp,0);
     for (ii = 0; ii < ntimes; ii++)
     {
       fprintf(fp, "E%d = [\n", ii+1);
       for (kk = 0; kk < nbins; kk++) fprintf(fp, "%d\n",  Fcounts[ii][kk]);
       fprintf(fp, "];\n");
     }
     fprintf(fp, "EE = [");
     for (ii = 0; ii < ntimes; ii++)
       fprintf(fp, "E%d/sum(E%d) ", ii+1, ii+1);
     fprintf(fp, "];\n");
     fprintf(fp, "subplot(1,2,1)\n");
     fprintf(fp, "plot(X,EE,'lineWidth',2)\n");
     fwritePlotAxes(fp);
     fwritePlotTitle(fp, "Prob. Dist. (-3,2,1,0,1,2,3 std.)");
     fwritePlotXLabel(fp, "Output Value");
     fwritePlotYLabel(fp, "Probabilities");
     fprintf(fp, "subplot(1,2,2)\n");
     for (ii = 0; ii < ntimes; ii++)
     {
       fprintf(fp, "for ii = 2 : %d\n", nbins);
       fprintf(fp, "   E%d(ii) = E%d(ii) + E%d(ii-1);\n",ii+1,ii+1,ii+1);
       fprintf(fp, "end;\n");
     }
     fprintf(fp, "EE = [");
     for (ii = 0; ii < ntimes; ii++)
       fprintf(fp, "E%d/E%d(%d) ", ii+1, ii+1, nbins);
     fprintf(fp, "];\n");
     fprintf(fp, "plot(X,EE,'linewidth',2)\n");
     fwritePlotAxes(fp);
     fwritePlotTitle(fp, "Cum. Dist. (-3,2,1,0,1,2,3 std.)");
     fwritePlotXLabel(fp, "Output Value");
     fwritePlotYLabel(fp, "Probabilities");
     fclose(fp);
     if (psPlotTool_ == 1)
          printf("Output distribution plots is in scilabrsua2.sci.\n");
     else printf("Output distribution plots is in matlabrsua2.m.\n");
     for (ii = 0; ii < ntimes; ii++) delete [] Fcounts[ii];
     delete [] Fcounts;
   }

   // +++ rs_uab 
   else if (!strcmp(command, "rs_uab"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rs_uab: this is a generic RS-based command that can\n");
         printf("        accommodate a discrepancy model, a pre-generated\n");
         printf("        sample (in a file), and bootstrapping.\n");
         printf("        The sample file should have the following format: \n");
         printf("PSUADE_BEGIN \n");
         printf("<nPts> <nInputs> \n");
         printf("1 <input 1> <input 2> ... \n");
         printf("2 <input 1> <input 2> ... \n");
         printf("...... \n");
         printf("PSUADE_END \n");
         printf("syntax: rs_uab (no argument needed)\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data (load sample first).\n");
         return -1;
      }
      if (nOutputs > 1)
      {
         printf("Currently this command does not support nOutputs > 1.\n");
         printf("Use 'write' to generate a one-output data file first.\n");
         return -1;
      }
      else
      {
         int    discFile=1,nInps,nOuts,dnSamp,it,ind,nSamples2,*tempI, nbs;
         double mean=0.0, stdev=0.0, dtemp;
         double *outVals, *tempX, *tempY, *tempV, *tempW, *inputVals=NULL;
         PsuadeData *localIO = NULL;
         FuncApprox **faPtrsRsEval=NULL;
         outputID = 0;
         faPtrsRsEval = new FuncApprox*[2];
         faPtrsRsEval[0] = NULL;
         faPtrsRsEval[1] = NULL;

         printf("Enter discrepancy model PSUADE file (if none, just 'n'): ");
         scanf("%s", winput);
         fgets(lineIn2,500,stdin); 
         if (winput[0] == 'n') discFile = 0; 
         else
         {
            localIO = new PsuadeData();
            status = localIO->readPsuadeFile(winput);
            if (status == 0)
            {
               localIO->getParameter("input_ninputs", pPtr);
               nInps = pPtr.intData_;
               if (nInps < nInputs)
               {
                  printf("Discrepancy model has %d inputs.\n", nInps);
                  printf("So the first %d inputs in the model file ",nInps);
                  printf("are assumed to associate with the inputs of\n");
                  printf("the discrepancy model.\n");
               }
               localIO->getParameter("output_noutputs", pPtr);
               nOuts = pPtr.intData_;
               if (nOuts > 1)
               {
                  printf("The discrepancy model has nOutputs > 1.\n");
                  printf("This is currently not supported.\n");
                  delete [] faPtrsRsEval;
                  delete localIO;
                  return -1;
               }
               printf("** CREATING RESPONSE SURFACE FOR DISCREPANCY MODEL\n");
               faPtrsRsEval[1] = genFAInteractive(localIO, 3);
               delete localIO;
               localIO = NULL;
            }
            else
            {
               printf("ERROR: in reading the discrepancy model file %s.\n",
                      winput);
               discFile = 0;
               delete [] faPtrsRsEval;
               faPtrsRsEval = NULL;
               delete localIO;
               localIO = NULL;
               return -1;
            }
         }

         printf("Enter sample file (in some standard format): ");
         scanf("%s", dataFile);
         fgets(lineIn2,500,stdin); 
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: sample data file %s not found.\n", dataFile);
            if (discFile == 1) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            faPtrsRsEval = NULL;
            if (localIO != NULL) delete localIO;
            localIO = NULL;
            return -1;
         }
         else
         {
            fscanf(fp, "%s", winput);
            if (strcmp(winput, "PSUADE_BEGIN"))
            {
               printf("ERROR: file must begin with PSUADE_BEGIN\n");
               fclose(fp);
               printf("File format: \n");
               printf("PSUADE_BEGIN \n");
               printf("<nPts> <nInputs> \n");
               printf("1 <input 1> <input 2> ... \n");
               printf("2 <input 1> <input 2> ... \n");
               printf("...... \n");
               printf("PSUADE_END \n");
               delete [] faPtrsRsEval;
               faPtrsRsEval = NULL;
               if (localIO != NULL) delete localIO;
               localIO = NULL;
               return -1;
            }
            else
            {
               fscanf(fp, "%d %d", &dnSamp, &kk);
               if (dnSamp <= 0)
               {
                  printf("ERROR: invalid sample size\n");
                  fclose(fp);
                  delete [] faPtrsRsEval;
                  faPtrsRsEval = NULL;
                  if (localIO != NULL) delete localIO;
                  localIO = NULL;
                  return -1;
               }
               if (kk != nInputs)
               {
                  printf("ERROR: input size does not match nInputs.\n");
                  printf(":      input size in local memory = %d.\n", 
                         nInputs);
                  printf(":      input size from file       = %d.\n",kk);
                  fclose(fp);
                  delete [] faPtrsRsEval;
                  if (localIO != NULL) delete localIO;
                  faPtrsRsEval = NULL;
                  localIO = NULL;
                  return -1;
               }
               inputVals = new double[dnSamp*nInputs];
               for (jj = 0; jj < dnSamp; jj++)
               {
                  fscanf(fp, "%d", &ind);
                  if (ind != (jj+1))
                  {
                     printf("ERROR: input index mismatch (%d,%d)\n",jj+1,ind);
                     printf("       read     index = %d\n", ind);
                     printf("       expected index = %d\n", jj+1);
                     printf("File format: \n");
                     printf("PSUADE_BEGIN \n");
                     printf("<nPts> <nInputs> \n");
                     printf("1 <input 1> <input 2> ... \n");
                     printf("2 <input 1> <input 2> ... \n");
                     printf("...... \n");
                     printf("PSUADE_END \n");
                     delete [] faPtrsRsEval;
                     delete [] inputVals;
                     if (localIO != NULL) delete localIO;
                     faPtrsRsEval = NULL;
                     localIO = NULL;
                     fclose(fp);
                     return -1;
                  }
                  for (ii = 0; ii < nInputs; ii++)
                     fscanf(fp, "%lg", &(inputVals[jj*nInputs+ii]));
               }
               if (jj != dnSamp)
               {
                  return -1;
                  fscanf(fp, "%s", winput);
                  fscanf(fp, "%s", winput);
                  if (strcmp(winput, "PSUADE_END"))
                  {
                     fclose(fp);
                     printf("ERROR: file must end with PSUADE_END\n");
                     delete [] faPtrsRsEval;
                     delete [] inputVals;
                     if (localIO != NULL) delete localIO;
                     faPtrsRsEval = NULL;
                     localIO = NULL;
                     return -1;
                  }
               }
               fclose(fp);
            }
            
            sprintf(pString, "How many bootstrapped samples to use (10 - 300) : ");
            nbs = getInt(1, 300, pString);
            printf("Write the CDFs to a matlab/scilab file? (y or n) ");
            scanf("%s", winput);
            fgets(lineIn,500,stdin); 
            flag = 0;
            fp = NULL;
            if (winput[0] == 'y')
            {
               if (dnSamp > 50000)
               {
                  printf("INFO: sample size %d too large (>50000) for matlab\n",
                         dnSamp);
                  printf("      plot. CDF plots not to be generated.\n");
               }
               else
               {
                  flag = 1; 
                  if (psPlotTool_ == 1)
                  {
                     fp = fopen("scilabrsuab_cdf.sci", "w");
                     if (fp == NULL)
                     {
                        printf("ERROR: cannot open file.\n");
                        flag = 0;
                     }
                     else
                     {
                        fprintf(fp, "// CDFs for rs_uab\n");
                        fwritePlotCLF(fp);
                     }
                  }
                  else
                  {
                     fp = fopen("matlabrsuab_cdf.m", "w");
                     if (fp == NULL)
                     {
                        printf("ERROR: cannot open file.\n");
                        flag = 0;
                     }
                     else
                     {
                        fprintf(fp, "%% CDFs for rs_uab\n");
                        fwritePlotCLF(fp);
                     }
                  }
               }
            }
            if (nbs == 1) nSamples2 = nSamples;
            else
            {
               nSamples2 = (int) (0.9 * nSamples);
               if ((double) nSamples2 / (double) nSamples < 0.9) nSamples2++;
            }
            faPtrsRsEval[0] = genFA(-1, nInputs, -1, nSamples2);
            if (faPtrsRsEval[0] == NULL)
            {
               printf("ERROR: cannot generate response surface.\n");
               delete [] faPtrsRsEval;
               faPtrsRsEval = NULL;
               if (inputVals != NULL) delete [] inputVals;
               inputVals = NULL;
               if (localIO != NULL) delete localIO;
               localIO = NULL;
               return -1;
            }
            faPtrsRsEval[0]->setBounds(iLowerB, iUpperB);
            faPtrsRsEval[0]->setOutputLevel(0);
            outVals = new double[dnSamp];
            tempX = new double[nSamples*nInputs];
            tempY = new double[nSamples];
            tempI = new int[nSamples];
            tempV = new double[nbs];
            tempW = new double[nbs];
            for (it = 0; it < nbs; it++)
            {
               printf("rs_uab: iteration %d\n", it+1);
               if (nbs == 1)
               {
                  for (jj = 0; jj < nSamples*nInputs; jj++)
                     tempX[jj] = sampleInputs[jj];
                  for (jj = 0; jj < nSamples; jj++)
                     tempY[jj] = sampleOutputs[jj*nOutputs+outputID];
               }
               else
               {   
                  for (jj = 0; jj < nSamples; jj++) tempI[jj] = 0;
                  kk = 0;
                  while (kk < nSamples2)
                  {
                     ind = PSUADE_rand() % nSamples;
                     if (tempI[ind] == 0)
                     {
                        for (ii = 0; ii < nInputs; ii++)
                           tempX[kk*nInputs+ii] = sampleInputs[ind*nInputs+ii];
                        tempY[kk] = sampleOutputs[ind*nOutputs+outputID];
                        tempI[ind] = 1;
                        kk++;
                     }
                  }
               }
               if (discFile == 1)
               {
                  for (jj = 0; jj < nSamples2; jj++)
                  {
                     dtemp = faPtrsRsEval[1]->evaluatePoint(&tempX[jj*nInputs]);
                     tempY[jj] += dtemp;
                  }
               }
               status = faPtrsRsEval[0]->initialize(tempX,tempY);
               faPtrsRsEval[0]->evaluatePoint(dnSamp, inputVals, outVals);
               mean = stdev = 0.0;
               for (jj = 0; jj < dnSamp; jj++) mean += outVals[jj];
               mean /= (double) dnSamp;
               for (jj = 0; jj < dnSamp; jj++)
                  stdev += pow(outVals[jj] - mean, 2.0);
               stdev = sqrt(stdev / dnSamp);
               tempV[it] = mean;
               tempW[it] = stdev;
               if (fp != NULL && flag == 1)
               {
                  fprintf(fp, "Y = [\n");
                  for (jj = 0; jj < dnSamp; jj++) 
                     fprintf(fp,"%e\n",outVals[jj]);
                  fprintf(fp, "];\n");
                  fprintf(fp, "Y%d = sort(Y);\n",it+1);
                  fprintf(fp, "X%d = (1 : %d)';\n", it+1, dnSamp);
                  fprintf(fp, "X%d = X%d / %d;\n", it+1, it+1, dnSamp);
                  if (it == 0)
                  {
                     fprintf(fp, "YY = Y%d;\n", it+1);
                     fprintf(fp, "XX = X%d;\n", it+1);
                  }
                  else
                  {
                     fprintf(fp, "YY = [YY Y%d];\n", it+1);
                     fprintf(fp, "XX = [XX X%d];\n", it+1);
                  }
               }
            }
            if (fp != NULL)
            {
               fprintf(fp, "plot(YY, XX, 'lineWidth',3)\n");
               fwritePlotTitle(fp, "Cumulative Distribution");
               fwritePlotAxes(fp);
               fwritePlotXLabel(fp, "Output Value");
               fwritePlotYLabel(fp, "Probabilities");
               fclose(fp);
               if (psPlotTool_ == 1)
                    printf("rs_uab: scilabrsuab_cdf.sci has the CDF plots.\n");
               else printf("rs_uab: matlabrsuab_cdf.m has the CDF plots.\n");
            }
            delete faPtrsRsEval[0];
            faPtrsRsEval[0] = NULL;
            if (discFile == 1) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            faPtrsRsEval = NULL;
            delete [] inputVals;
            delete [] outVals;
            delete [] tempX;
            delete [] tempY;
            delete [] tempI;
            if (it == nbs && nbs > 1)
            {
               mean = 0.0;
               for (jj = 0; jj < nbs; jj++) mean += tempV[jj];
               mean /= (double) nbs;
               for (jj = 0; jj < nbs; jj++)
                  stdev += pow(tempV[jj] - mean, 2.0);
               stdev = sqrt(stdev / nbs);
               printf("rs_uab: Sample Mean  = %e (%e)\n", mean, stdev);
               mean = 0.0;
               for (jj = 0; jj < nbs; jj++) mean += tempW[jj];
               mean /= (double) nbs;
               for (jj = 0; jj < nbs; jj++)
                  stdev += pow(tempW[jj] - mean, 2.0);
               stdev = sqrt(stdev / nbs);
               printf("rs_uab: Sample Stdev = %e (%e)\n", mean, stdev);
            }
            else if (kk == nbs && nbs == 1)
            {
               printf("rs_uab: Sample Mean  = %e\n", tempV[0]);
               printf("rs_uab: Sample Stdev = %e\n", tempW[0]);
            }
            delete [] tempW;
            delete [] tempV;
         }
      }
   }

   // +++ rs_uap 
   else if (!strcmp(command, "rs_uap"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rs_uap: rs_ua when the response surface is appended\n");
         printf("        with a discrepancy model and the sample is.\n");
         printf("        provided by users.\n");
         printf("syntax: rs_uap (no argument needed)\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data (load sample first).\n");
         return -1;
      }
      if (nOutputs > 1)
      {
         printf("Currently this command does not support nOutputs > 1.\n");
         printf("Use 'write' to choose one output for processing.\n");
         return -1;
      }
      else
      {
         int    discFile=1, nInps, nOuts, faFlag;
         FuncApprox **faPtrsRsEval=NULL;
         PsuadeData *localIO=NULL;
         faFlag = 3;
         psuadeIO->getParameter("ana_outputid", pPtr);
         kk = pPtr.intData_;
         outputID = 0;
         faPtrsRsEval = new FuncApprox*[2];
         psuadeIO->updateAnalysisSection(-1, -1, -1, -1, outputID, -1);
         faPtrsRsEval[0] = genFAInteractive(psuadeIO, faFlag);
         faPtrsRsEval[1] = NULL;
         psuadeIO->updateAnalysisSection(-1, -1, -1, -1, kk, -1);
         printf("Enter discrepancy model PSUADE file (if none, enter NONE): ");
         scanf("%s", winput);
         fgets(lineIn2,500,stdin); 
         if (!strcmp(winput, "NONE")) discFile = 0; 
         else
         {
            localIO = new PsuadeData();
            status = localIO->readPsuadeFile(winput);
            if (status == 0)
            {
               localIO->getParameter("input_ninputs", pPtr);
               nInps = pPtr.intData_;
               if (nInps < nInputs)
               {
                  printf("Discrepancy model has %d inputs.\n", nInps);
                  printf("So the first %d inputs in the model file ",nInps);
                  printf("are assumed to associate with the inputs of\n");
                  printf("the discrepancy model.\n");
               }
               localIO->getParameter("output_noutputs", pPtr);
               nOuts = pPtr.intData_;
               if (nOuts > 1)
               {
                  printf("The discrepancy model has nOutputs > 1.\n");
                  printf("This is currently not supported.\n");
                  delete localIO;
                  if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
                  delete [] faPtrsRsEval;
                  return -1;
               }
               printf("** CREATING RESPONSE SURFACE FOR DISCREPANCY MODEL\n");
               faPtrsRsEval[1] = genFAInteractive(localIO, 3);
               delete localIO;
            }
            else
            {
               printf("ERROR: in reading the discrepancy model file %s.\n",
                      winput);
               discFile = 0;
               delete localIO;
               delete faPtrsRsEval[0];
               delete [] faPtrsRsEval;
               return -1;
            }
            localIO = NULL;
         }

         int    dnInps, dnSamp, ind;
         double *inputVals=NULL, *outVals;
         char   dataFile[1001];
         printf("Enter sample file: ");
         scanf("%s", dataFile);
         fgets(lineIn2,500,stdin); 
         fp = fopen(dataFile, "r");
         if (fp == NULL)
         {
            printf("ERROR: sample data file %s not found.\n", dataFile);
            if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
            if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            return -1;
         }
         else
         {
            fscanf(fp, "%s", winput);
            if (strcmp(winput, "PSUADE_BEGIN"))
            {
               printf("ERROR: file must begin with PSUADE_BEGIN\n");
               fclose(fp);
               printf("File format: \n");
               printf("PSUADE_BEGIN \n");
               printf("<nPts> <nInputs> \n");
               printf("1 <input 1> <input 2> ... \n");
               printf("2 <input 1> <input 2> ... \n");
               printf("...... \n");
               printf("PSUADE_END \n");
               if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
               if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
               delete [] faPtrsRsEval;
               return -1;
            }
            else
            {
               while (1)
               {
                  kk = getc(fp);
                  if (kk != '#')
                  {
                     ungetc(kk, fp);
                     break;
                  }
               }
               fscanf(fp, "%d %d", &dnSamp, &dnInps);
               if (dnSamp <= 0)
               {
                  printf("ERROR: invalid sample size\n");
                  fclose(fp);
                  if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
                  if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
                  delete [] faPtrsRsEval;
                  return -1;
               }
               printf("Sample size read = %d\n", dnSamp);
               if (dnInps != nInputs)
               {
                  printf("ERROR: input size does not match nInputs.\n");
                  printf(":      input size in local memory = %d.\n", 
                         nInputs);
                  printf(":      input size from file       = %d.\n",dnInps);
                  fclose(fp);
                  if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
                  if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
                  delete [] faPtrsRsEval;
                  return -1;
               }
               fgets(lineIn2, 500, fp);
               while (1)
               {
                  kk = getc(fp);
                  if (kk != '#')
                  {
                     ungetc(kk, fp);
                     break;
                  }
                  else fgets(lineIn2, 500, fp);
               }
               inputVals = new double[dnSamp*dnInps];
               outVals = new double[dnSamp];
               for (jj = 0; jj < dnSamp; jj++)
               {
                  fscanf(fp, "%d", &ind);
                  if (ind != (jj+1))
                  {
                     printf("ERROR: input index mismatch (%d,%d)\n",
                            jj+1,ind);
                     printf("       read     index = %d\n", ind);
                     printf("       expected index = %d\n", jj+1);
                     delete [] inputVals;
                     delete [] outVals;
                     if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
                     if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
                     delete [] faPtrsRsEval;
                     return -1;
                  }
                  for (ii = 0; ii < nInputs; ii++)
                     fscanf(fp, "%lg", &(inputVals[jj*dnInps+ii]));
               }
               fgets(lineIn2, 500, fp);
               fscanf(fp, "%s", winput);
               if (strcmp(winput, "PSUADE_END"))
               {
                  fclose(fp);
                  printf("ERROR: file must end with PSUADE_END (%s)\n",winput);
                  delete [] inputVals;
                  delete [] outVals;
                  if (faPtrsRsEval[0] != NULL) delete faPtrsRsEval[0];
                  if (faPtrsRsEval[1] != NULL) delete faPtrsRsEval[1];
                  delete [] faPtrsRsEval;
                  return -1;
               }
            }
            fclose(fp);
         }
         faPtrsRsEval[0]->evaluatePoint(dnSamp, inputVals, outVals);
         if (discFile == 1)
         {
            for (jj = 0; jj < dnSamp; jj++)
            {
               ddata = faPtrsRsEval[1]->evaluatePoint(&inputVals[jj*nInputs]);
               outVals[jj] += ddata;
            }
         }
         double mean=0.0, stdev=0.0;
         for (ss = 0; ss < dnSamp; ss++) mean += outVals[ss];
         mean /= (double) dnSamp;
         for (ss = 0; ss < dnSamp; ss++) stdev += pow(outVals[ss] - mean, 2.0);
         stdev = sqrt(stdev / dnSamp);
         if (psPlotTool_ == 1)
         {
            fp = fopen("scilabrsuap.sci", "w");
            if (fp == NULL)
               printf("rs_uap ERROR: cannot open scilabrsuap.sci file.\n");
         }
         else
         {
            fp = fopen("matlabrsuap.m", "w");
            if (fp == NULL)
               printf("rs_uap ERROR: cannot open matlabrsuap.m file.\n");
         }
         if (fp != NULL)
         {
            fprintf(fp, "Y = [ \n");
            for (jj = 0; jj < dnSamp; jj++)
               fprintf(fp, "%16.8e\n", outVals[jj]);
            fprintf(fp, "];\n");
            if (psPlotTool_ == 1)
            {
               fprintf(fp, "histplot(10, Y, style=2);\n");
               fprintf(fp, "a = gce();\n");
               fprintf(fp, "a.children.fill_mode = \"on\";\n");
               fprintf(fp, "a.children.thickness = 2;\n");
               fprintf(fp, "a.children.foreground = 0;\n");
               fprintf(fp, "a.children.background = 2;\n");
            }
            else
            {
               fprintf(fp, "[nk,xk]=hist(Y,10);\n");
               fprintf(fp, "bar(xk,nk/%d,1.0)\n",dnSamp);
            }
            fwritePlotAxes(fp);
            fwritePlotTitle(fp, "Probability Distribution");
            fwritePlotXLabel(fp, "Output Value");
            fwritePlotYLabel(fp, "Probabilities");
            fclose(fp);
            if (psPlotTool_ == 1)
                 printf("rs_uap: distribution in scilabrsuap.scin");
            else printf("rs_uap: distribution in matlabrsuap.m.\n");
            printf("Sample mean  = %e\n", mean);
            printf("Sample stdev = %e\n", stdev);
            delete [] inputVals;
            delete [] outVals;
            delete faPtrsRsEval[0];
            if (discFile == 1) delete faPtrsRsEval[1];
            delete [] faPtrsRsEval;
            faPtrsRsEval = NULL;
         }
      }
   }

   // several qsa methods 
   else if (!strcmp(command, "rs_qsa2"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs_qsa: RS-based quantitative sensitivity analysis\n");
       printf("syntax: rs_qsa (no argument needed)\n");
       printf("Note: to facilitate processing, all expert modes have\n");
       printf("      been suppressed.\n");
       printf("Note: This command differs from rssobol1, rssoboltsi,\n");
       printf("      and the command 'me' in that it uses bootstrapped\n");
       printf("      samples multiple times to get the errors in Sobol'\n");
       printf("      indices due to response surface errors.\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data (load sample first).\n");
       return -1;
     }
     outputID = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     outputID = getInt(1, nOutputs, pString);
     outputID--;

     printf("Which main or total effect analyzer ? \n");
     printf("1. Sobol' main effect using McKay's method (replicated LH).\n");
     printf("2. Sobol' main effect method using on numerical integration. \n");
     printf("3. Sobol' total sensitivity using numerical integration.\n");
     sprintf(pString, "Which  method (1, 2, or 3) ? ");
     int method = getInt(1, 3, pString);
     int analysisMethod, iOne=1;
     if      (method == 1) analysisMethod = PSUADE_ANA_ME;
     else if (method == 2) analysisMethod = PSUADE_ANA_RSSOBOL1;
     else if (method == 3) analysisMethod = PSUADE_ANA_RSSOBOLTSI;
     AnalysisManager *anaManager = new AnalysisManager();
     anaManager->setup(analysisMethod, 0);

     FuncApprox *faPtr;
     if (method == 1)
     {
       faType = -1;
       faPtr = genFA(faType, nInputs, iOne, nSamples);
       faPtr->setBounds(iLowerB, iUpperB);
       faPtr->setOutputLevel(0);
     }

     psuadeIO->getParameter("ana_diagnostics",pPtr);
     int saveDiag = pPtr.intData_;
     psuadeIO->updateAnalysisSection(-1,-1,-1,-1,-1,-1);
     int saveMode = psAnaExpertMode_;
     psAnaExpertMode_ = 0;
     psuadeIO->getParameter("method_sampling", pPtr);
     int saveMethod = pPtr.intData_;

     Sampling *sampPtr;
     psVector vecUpper, vecLower, vecIn, vecOut;
     int    usePDFs, ind;
     int    count2 = 100000;
     int    nReps  = 200;
     double *tempM = new double[nSamples*nInputs];
     double *tempV = new double[nSamples];
     int    *states = new int[nSamples];
     double *tempX  = new double[count2*nInputs];
     double *tempY  = new double[count2];
     int    *tempI  = new int[count2];
     for (ii = 0; ii < count2*nInputs; ii++) tempX[ii] = 0.0;
     for (ii = 0; ii < count2; ii++) tempY[ii] = 0.0;
     for (ii = 0; ii < count2; ii++) tempI[ii] = 1;

     sprintf(pString, "How many times to run it (10 - 1000) : ");
     int count = getInt(10, 1000, pString);
     double *tempT  = new double[count*nInputs];
     psuadeIO->getParameter("ana_use_input_pdfs", pPtr);
     usePDFs = pPtr.intData_;
     PDFManager *pdfman = NULL;
     if (usePDFs == 1 && method == 1)
     {
       printf("NOTE: Some inputs have non-uniform distributions.\n");
       pdfman = new PDFManager();
       pdfman->initialize(psuadeIO);
       vecOut.setLength(count2*nInputs);
       vecUpper.load(nInputs, iUpperB);
       vecLower.load(nInputs, iLowerB);
     }

     for (kk = 0; kk < count; kk++)
     {
       printf("rq_qsa: iteration %d\n", kk+1);
       for (ss = 0; ss < nSamples; ss++)
       {
         ind = PSUADE_rand() % nSamples;
         for (ii = 0; ii < nInputs; ii++)
           tempM[ss*nInputs+ii] = sampleInputs[ind*nInputs+ii];
         tempV[ss] = sampleOutputs[ind*nOutputs+outputID];
         states[ss] = sampleStates[ss];
       }
       if (method == 1)
       {
         status = faPtr->initialize(tempM,tempV);

         sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
         sampPtr->setPrintLevel(0);
         sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
         sampPtr->setOutputParams(1);
         sampPtr->setSamplingParams(count2, nReps, 1);
         sampPtr->initialize(0);
         sampPtr->getSamples(count2, nInputs, 1, tempX, tempY, tempI);
         if (usePDFs == 1)
         {
           vecIn.load(count2*nInputs, tempX);
           pdfman->invCDF(count2, vecIn, vecOut, vecLower, vecUpper);
           for (ii = 0; ii < count2*nInputs; ii++)
             tempX[ii] = vecOut[ii];
         }

         faPtr->evaluatePoint(count2, tempX, tempY);

         psuadeIO->updateInputSection(count2,nInputs,NULL,NULL,NULL,
                                      tempX,NULL,NULL,NULL,NULL,NULL);
         psuadeIO->updateOutputSection(count2,1,tempY,tempI,NULL);
         psuadeIO->updateMethodSection(PSUADE_SAMP_LHS,count2,
                                       nReps,-1,-1);
       }
       else
       {
         psuadeIO->updateInputSection(nSamples,nInputs,NULL,NULL,
                                      NULL,tempM,NULL,NULL,NULL,NULL,NULL);
         psuadeIO->updateOutputSection(nSamples,1,tempV,states,
                                       outputNames);
         psuadeIO->updateMethodSection(PSUADE_SAMP_MC,nSamples,
                                       -1,-1,-1);
       }
       anaManager->analyze(psuadeIO, 0, NULL, 0);
       pData *pdata = psuadeIO->getAuxData();
       if (pdata->nDbles_ != nInputs)
       {
         printf("ERROR: nInputs do not match (%d, %d).\n",
                pdata->nDbles_, nInputs);
         printf("       Consult PSUADE developers.\n");
         delete [] tempX;
         delete [] tempY;
         delete [] tempI;
         delete [] tempT;
         delete [] tempM;
         delete [] tempV;
         delete [] states;
         if (method == 1) delete sampPtr;
         return -1;
       }

       if (pdata->dbleData_ > 0)
         for (ii = 0; ii < nInputs; ii++)
           tempT[kk*nInputs+ii] =
                 pdata->dbleArray_[ii]/pdata->dbleData_;
       else
         for (ii = 0; ii < nInputs; ii++)
           tempT[kk*nInputs+ii] = pdata->dbleArray_[ii];

       pdata->clean();
       if (method == 1) delete sampPtr;
     }
     if (usePDFs == 1 && method == 1) delete pdfman;
     delete [] tempM;
     delete [] tempV;
     delete [] states;
     tempM = new double[nInputs];
     for (ii = 0; ii < nInputs; ii++)
     {
       tempM[ii] = tempT[ii];
       for (jj = 1; jj < count; jj++) tempM[ii] += tempT[jj*nInputs+ii];
       tempM[ii] /= (double) count;
     }
     tempV = new double[nInputs];
     for (ii = 0; ii < nInputs; ii++)
     {
       tempV[ii] = pow(tempT[ii]-tempM[ii], 2.0);
       for (jj = 1; jj < count; jj++)
         tempV[ii] += pow(tempT[jj*nInputs+ii]-tempM[ii],2.0);
       tempV[ii] /= (double) (count - 1);
       tempV[ii] = sqrt(tempV[ii]);
     }
     printEquals(PL_INFO, 0);
     printf("Statistics (based on %d replications): \n", count);
     for (ii = 0; ii < nInputs; ii++)
       printf("Input %4d: mean = %16.8e, std = %16.8e\n",ii+1,
              tempM[ii],tempV[ii]);
     delete anaManager;
     delete [] tempX;
     delete [] tempY;
     delete [] tempI;
     delete [] tempT;
     delete [] tempM;
     delete [] tempV;
     if (faPtr != NULL) delete faPtr;

     psuadeIO->updateInputSection(nSamples,nInputs,NULL,NULL,NULL,
                                  sampleInputs,NULL,NULL,NULL,NULL,NULL);
     psuadeIO->updateOutputSection(nSamples,nOutputs,sampleOutputs,
                                    sampleStates,outputNames);
     psuadeIO->updateMethodSection(saveMethod,nSamples,-1,-1,-1);
     psuadeIO->updateAnalysisSection(-1,-1,-1,saveDiag,-1,-1);
     psAnaExpertMode_ = saveMode;
   }

   // +++ rs1
   else if (!strcmp(command, "rs1"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs1: response surface plot in one parameter\n");
       printf("syntax: rs1 (no argument needed)\n");
       return -1;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data to analyze (load sample first).\n");
       return -1;
     }
     int nPtsPerDim = 256;
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);
     double *inputSettings = new double[nInputs];
     int    iplot1, iInd1, jplot, sInd;
     if (nInputs > 1)
     {
       sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
       iplot1 = getInt(1, nInputs, pString);
       iplot1--;
     }
     else iplot1 = 0;
     if (nInputs > 1)
     {
       sprintf(pString,
               "Set other nominal values automatically ? (y or n) ");
       getString(pString, winput);
     }
     if (winput[0] == 'y')
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
         else inputSettings[iInd1] = 1.0;
       }
     }
     else
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1)
         {
           sprintf(pString,
                   "Enter nominal value for input %d (%e - %e): ", 
                   iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
           inputSettings[iInd1] = getDouble(pString);
         }
         else inputSettings[iInd1] = 1.0;
       }
     }
     jplot = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     jplot = getInt(1, nOutputs, pString);
     jplot--;

     int    faLeng=0;
     double *faXOut=NULL, *faYOut=NULL;
     double *faYIn = new double[nSamples];
     for (sInd = 0; sInd < nSamples; sInd++)
       faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

     faPtr->gen1DGridData(sampleInputs,faYIn,iplot1,
                   inputSettings, &faLeng, &faXOut,&faYOut);

     if (psPlotTool_ == 1)
     {
       fp = fopen("scilabrs1.sci", "w");
       if (fp == NULL)
       {
         printf("ERROR: cannot open file scilabrs1.sci.\n");
         delete [] inputSettings;
         delete [] faYIn;
         delete [] faXOut;
         delete [] faYOut;
         delete faPtr;
         return -1;
       }
       fwritePlotCLF(fp);
       fprintf(fp, "A = [\n");
       for (sInd = 0; sInd < faLeng; sInd++)
         fprintf(fp, "%e\n", faYOut[sInd]);
       fprintf(fp, "];\n");
       fprintf(fp, "X = [\n");
       for (sInd = 0; sInd < faLeng; sInd++)
         fprintf(fp, "%e\n", faXOut[sInd]);
       fprintf(fp, "];\n");
       fprintf(fp, "plot(X,A,'-')\n");
       fprintf(fp, "a = gca();\n");
       fprintf(fp, "a.children.children.thickness = 4;\n");
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, outputNames[jplot]);
       sprintf(winput, "Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fclose(fp);
       printf("scilabrs1.sci is now available.\n");
     }
     else
     {
       fp = fopen("matlabrs1.m", "w");
       if (fp == NULL)
       {
         printf("ERROR: cannot open file matlabrs1.m.\n");
         delete [] inputSettings;
         delete [] faYIn;
         delete [] faXOut;
         delete [] faYOut;
         delete faPtr;
         return -1;
       }
       fwritePlotCLF(fp);
       fprintf(fp, "A = [\n");
       for (sInd = 0; sInd < faLeng; sInd++)
         fprintf(fp, "%e\n", faYOut[sInd]);
       fprintf(fp, "];\n");
       fprintf(fp, "X = [\n");
       for (sInd = 0; sInd < faLeng; sInd++)
         fprintf(fp, "%e\n", faXOut[sInd]);
       fprintf(fp, "];\n");
       fprintf(fp, "plot(X,A,'-','lineWidth',4)\n");
       fprintf(fp, "hold on\n");
       double Ymin = faYOut[0];
       for (sInd = 1; sInd < faLeng; sInd++)
         if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
       double Ymax = faYOut[0];
       for (sInd = 1; sInd < faLeng; sInd++)
         if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];
       printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
       sprintf(pString,"Set lower threshold ? (y or n) : ");
       getString(pString, winput);
       fprintf(fp, "yminFlag = 0;\n");
       double thresh;
       if (winput[0] == 'y')
       {
         sprintf(pString,"Enter the lower threshold (min = %e) : ", 
                 Ymin);
         thresh = getDouble(pString);
         fprintf(fp, "ymin = %e;\n", thresh);
         fprintf(fp, "plot(X,ones(%d,1)*ymin,'r-')\n",faLeng);
       }
       sprintf(pString,"Set upper threshold ? (y or n) : ");
       getString(pString, winput);
       if (winput[0] == 'y')
       {
         sprintf(pString,"Enter the upper threshold (max = %e) : ", 
                    Ymax);
         thresh = getDouble(pString);
         fprintf(fp, "ymax = %e;\n", thresh);
         fprintf(fp, "plot(X,ones(%d,1)*ymax,'r-')\n",faLeng);
       }
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, outputNames[jplot]);
       sprintf(winput, "Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fclose(fp);
       printf("matlabrs1.m is now available.\n");
     }
     delete [] faXOut;
     delete [] faYIn;
     delete [] faYOut;
     delete [] inputSettings;
     delete faPtr;
   }

   // +++ rs1s
   else if (!strcmp(command, "rs1s"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs1s: 1-parameter response surface (with uncertainty) plot\n");
       printf("syntax: rs1s (no argument needed)\n");
       return -1;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data to analyze (load sample first).\n");
       return -1;
     }
     int nPtsPerDim = 256;
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);
     double *inputSettings = new double[nInputs];
     int    iplot1, iInd1, jplot, sInd;
     if (nInputs > 1)
     {
       sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
       iplot1 = getInt(1, nInputs, pString);
       iplot1--;
     }
     else iplot1 = 0;
     if (nInputs > 1)
     {
       sprintf(pString,
               "Set other nominal values automatically ? (y or n) ");
       getString(pString, winput);
     }
     if (winput[0] == 'y')
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
         else inputSettings[iInd1] = 1.0;
       }
     }
     else
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1)
         {
           sprintf(pString,
                   "Enter nominal value for input %d (%e - %e): ", 
                   iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
           inputSettings[iInd1] = getDouble(pString);
         }
         else inputSettings[iInd1] = 1.0;
       }
     }
     jplot = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     jplot = getInt(1, nOutputs, pString);
     jplot--;

     double *faYIn = new double[nSamples];
     for (sInd = 0; sInd < nSamples; sInd++)
       faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];
     faPtr->initialize(sampleInputs,faYIn);
     delete [] faYIn;

     double *faXOut=NULL, *faYOut=NULL, *faYStds=NULL, hx;
     hx = 1.0 / (double) (nPtsPerDim - 1.0);
     faXOut = new double[nPtsPerDim*nInputs];
     for (ii = 0; ii < nPtsPerDim; ii++)
     {
       for (jj = 0; jj < nInputs; jj++)
         faXOut[ii*nInputs+jj] = inputSettings[jj];
       faXOut[ii*nInputs+iplot1] = hx * ii + iLowerB[iplot1];
     }
     faYOut = new double[nPtsPerDim];
     faYStds = new double[nPtsPerDim];

     faPtr->evaluatePointFuzzy(nPtsPerDim, faXOut, faYOut,faYStds);

     if (psPlotTool_ == 1)
     {
       fp = fopen("scilabrs1s.sci", "w");
       if (fp == NULL)
       {
         printf("ERROR: cannot open file scilabrs1s.sci.\n");
         delete [] inputSettings;
         delete [] faYStds;
         delete [] faXOut;
         delete [] faYOut;
         delete faPtr;
         return -1;
       }
       fwritePlotCLF(fp);
       fprintf(fp, "A = [\n");
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
         fprintf(fp, "%16.8e %16.8e\n", faYOut[sInd], faYStds[sInd]);
       fprintf(fp, "];\n");
       fprintf(fp, "X = [\n");
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
         fprintf(fp, "%e\n", faXOut[sInd*nInputs+iplot1]);
       fprintf(fp, "];\n");
       fprintf(fp, "plot(X,A(:,1),'-')\n");
       fprintf(fp, "a = gca();\n");
       fprintf(fp, "a.children.children.thickness = 4;\n");
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, outputNames[jplot]);
       sprintf(winput, "Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fclose(fp);
       printf("scilabrs1s.sci is now available.\n");
     }
     else
     {
       fp = fopen("matlabrs1s.m", "w");
       if (fp == NULL)
       {
         printf("ERROR: cannot open file matlabrs1s.m.\n");
         delete [] inputSettings;
         delete [] faYStds;
         delete [] faXOut;
         delete [] faYOut;
         delete faPtr;
         return -1;
       }
       fwritePlotCLF(fp);
       fprintf(fp, "%% 1D plot with +/- 2 std dev\n");
       fprintf(fp, "A = [\n");
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
         fprintf(fp, "%16.8e %16.8e\n", faYOut[sInd], faYStds[sInd]);
       fprintf(fp, "];\n");
       fprintf(fp, "X = [\n");
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
         fprintf(fp, "%e\n", faXOut[sInd*nInputs+iplot1]);
       fprintf(fp, "];\n");
       fprintf(fp, "plot(X,A(:,1),'-','lineWidth',4)\n");
       fprintf(fp, "hold on\n");
       fprintf(fp, "for ii = 1 : %d\n", nPtsPerDim);
       fprintf(fp, "  xx = [X(ii) X(ii)];\n");
       fprintf(fp, "  yy = [A(ii,1)-2*A(ii,2) A(ii,1)+2*A(ii,2)];\n");
       fprintf(fp, "  plot(xx,yy,'b-','lineWidth',2)\n");
       fprintf(fp, "end\n");
       
       double Ymin = faYOut[0];
       for (sInd = 1; sInd < nPtsPerDim; sInd++)
         if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
       double Ymax = faYOut[0];
       for (sInd = 1; sInd < nPtsPerDim; sInd++)
         if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];
       printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
       sprintf(pString,"Set lower threshold ? (y or n) : ");
       getString(pString, winput);
       fprintf(fp, "yminFlag = 0;\n");
       double thresh;
       if (winput[0] == 'y')
       {
         sprintf(pString,"Enter the lower threshold (min = %e) : ", 
                 Ymin);
         thresh = getDouble(pString);
         fprintf(fp, "ymin = %e;\n", thresh);
         fprintf(fp, "plot(X,ones(%d,1)*ymin,'r-')\n",nPtsPerDim);
       }
       sprintf(pString,"Set upper threshold ? (y or n) : ");
       getString(pString, winput);
       if (winput[0] == 'y')
       {
         sprintf(pString,"Enter the upper threshold (max = %e) : ", 
                    Ymax);
         thresh = getDouble(pString);
         fprintf(fp, "ymax = %e;\n", thresh);
         fprintf(fp, "plot(X,ones(%d,1)*ymax,'r-')\n",nPtsPerDim);
       }
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, outputNames[jplot]);
       sprintf(winput, "Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fclose(fp);
       printf("matlabrs1s.m is now available.\n");
     }
     delete [] faXOut;
     delete [] faYStds;
     delete [] faYOut;
     delete [] inputSettings;
     delete faPtr;
   }

   // +++ rs2
   else if (!strcmp(command, "rs2"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs2: response surface plot in two parameters\n");
       printf("syntax: rs2 (no argument needed)\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data to analyze (load sample first).\n");
       return -1;
     }
     if (nInputs < 2)
     {
       printf("ERROR: rs2 requires 2 or more inputs.\n");
       return -1;
     }
     int nPtsPerDim = 64;
     sprintf(pString, "Grid resolution ? (32 - 256) ");
     nPtsPerDim = getInt(32, 256, pString);
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);
     double *inputSettings = new double[nInputs];
     int    iplot1=-1, iplot2=-1;
     sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
     iplot1 = getInt(1, nInputs, pString);
     iplot1--;
     if (nInputs == 2)
     {
       if (iplot1 == 0) iplot2 = 1;
       else             iplot2 = 0;
     }
     while (iplot2 < 0 || iplot2 >= nInputs || iplot1 == iplot2)
     {
       sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
               nInputs, iplot1+1);
       iplot2 = getInt(1, nInputs, pString);
       iplot2--;
       if (iplot2 == iplot1)
       {
         printf("ERROR: cannot have same index for x and y axes.\n");
         iplot2 = -1;
       }
     }
     if (nInputs > 2)
     {
       sprintf(pString,
               "Set other nominal values automatically ? (y or n) ");
       getString(pString, winput);
     }
     int iInd1, sInd, jplot;
     if (winput[0] == 'y')
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
         else inputSettings[iInd1] = 1.0;
       }
     }
     else
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2)
         {
           sprintf(pString,
                   "Enter nominal value for input %d (%e - %e): ", 
                   iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
           inputSettings[iInd1] = getDouble(pString);
         }
         else inputSettings[iInd1] = 1.0;
       }
     }
     jplot = 0;
     sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
     jplot = getInt(1, nOutputs, pString);
     jplot--;

     int    faLeng=0;
     double *faXOut=NULL, *faYOut=NULL;
     double *faYIn = new double[nSamples];
     for (sInd = 0; sInd < nSamples; sInd++)
        faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

     if (outputLevel_ > 1)
        printf("Please wait while generating RS ....\n");
     faPtr->gen2DGridData(sampleInputs,faYIn,iplot1,iplot2,
                   inputSettings, &faLeng, &faXOut,&faYOut);

     if (psPlotTool_ == 1)
     {
       fp = fopen("scilabrs2.sci", "w");
       if (fp == NULL)
       {
         printf("ERROR: cannot open file scilabrs2.sci.\n");
         delete [] inputSettings;
         delete [] faYIn;
         delete [] faXOut;
         delete [] faYOut;
         delete faPtr;
         return -1;
       }
       fprintf(fp,"twoPlots = 1;\n");
       fprintf(fp,"A = [\n");
       for (sInd = 0; sInd < faLeng; sInd++)
         fprintf(fp, "%e\n", faYOut[sInd]);
       fprintf(fp,"];\n");
       fprintf(fp,"A = matrix(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
       fprintf(fp,"x = [\n");
       for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
         fprintf(fp, "%e\n", faXOut[sInd*2]);
       fprintf(fp,"];\n");
       fprintf(fp,"y = [\n");
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
         fprintf(fp, "%e\n", faXOut[sInd*2+1]);
       fprintf(fp,"];\n");
       fwritePlotCLF(fp);
       fprintf(fp,"if twoPlots == 1\n");
       fprintf(fp,"drawlater\n");
       fprintf(fp,"subplot(1,2,1)\n");
       fprintf(fp,"mesh(x,y,A)\n");
       fprintf(fp,"h = get(\"hdl\");\n");
       fprintf(fp,"h.color_flag=1;\n");
       fprintf(fp,"h.color_mode=-2;\n");
       fprintf(fp,"bmin = min(min(A)); bmax = max(max(A));\n");
       fprintf(fp,"xset(\"colormap\",jetcolormap(64));\n");
       fprintf(fp,"colorbar(bmin,bmax);\n");
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, inputNames[iplot2]);
       fwritePlotZLabel(fp, outputNames[jplot]);
       sprintf(winput, "Mesh Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fprintf(fp,"a=gca();\n");
       fprintf(fp,"a.data_bounds=[%e,%e;%e,%e];\n",iLowerB[iplot1],
               iLowerB[iplot2], iUpperB[iplot1], iUpperB[iplot2]);
       fprintf(fp,"a.axes_visible=\"on\";\n");
       fprintf(fp,"drawnow\n");
       fprintf(fp,"subplot(1,2,2)\n");
       fprintf(fp,"end;\n");
       fprintf(fp,"drawlater\n");
       fprintf(fp,"B = A;\n");
       fprintf(fp,"nX = length(x);\n");
       fprintf(fp,"nY = length(y);\n");
       fprintf(fp,"for ii = 1 : nX\n");
       fprintf(fp,"for jj = 1 : nY\n");
       fprintf(fp,"B(ii,jj) = A(nX-ii+1,jj);\n");
       fprintf(fp,"end;\n");
       fprintf(fp,"end;\n");
       fprintf(fp,"a=gca();\n");
       fprintf(fp,"a.data_bounds=[%e,%e;%e,%e];\n",iLowerB[iplot1],
               iLowerB[iplot2], iUpperB[iplot1], iUpperB[iplot2]);
       fprintf(fp,"bmin = min(min(B)); bmax = max(max(B));\n");
       fprintf(fp,"Matplot1((B-bmin)/(bmax-bmin)*64,[%e,%e,%e,%e])\n",
               iLowerB[iplot1],iLowerB[iplot2],iUpperB[iplot1], 
               iUpperB[iplot2]);
       fprintf(fp,"set(gca(),\"auto_clear\",\"off\")\n");
       fprintf(fp,"contour2d(x,y,flipdim(B',1),6);\n");
       fprintf(fp,"xset(\"colormap\",jetcolormap(64));\n");
       fprintf(fp,"colorbar(bmin,bmax);\n");
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, inputNames[iplot2]);
       sprintf(winput, "Contour Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fprintf(fp,"drawnow\n");
       fclose(fp);
       printf("scilabrs2.sci is now available for response surface and ");
       printf("contour plots\n");
     }
     else
     {
       fp = fopen("matlabrs2.m", "w");
       if (fp == NULL)
       {
         printf("ERROR: cannot open file matlabrs2.m.\n");
         delete [] inputSettings;
         delete [] faYIn;
         delete [] faXOut;
         delete [] faYOut;
         delete faPtr;
         return -1;
       }
       fwritePlotCLF(fp);
       fprintf(fp, "twoPlots = 1;\n");
       fprintf(fp, "A = [\n");
       for (sInd = 0; sInd < faLeng; sInd++)
         fprintf(fp, "%e\n", faYOut[sInd]);
       fprintf(fp, "];\n");
       fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
       fprintf(fp, "x = [\n");
       for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
         fprintf(fp, "%e\n", faXOut[sInd*2]);
       fprintf(fp, "];\n");
       fprintf(fp, "y = [\n");
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
         fprintf(fp, "%e\n", faXOut[sInd*2+1]);
       fprintf(fp, "];\n");
       if (nInputs == 2)
       {
         fprintf(fp, "xx = [\n");
         for (sInd = 0; sInd < nSamples; sInd++)
           fprintf(fp, "%e\n", sampleInputs[sInd*2]);
         fprintf(fp, "];\n");
         fprintf(fp, "yy = [\n");
         for (sInd = 0; sInd < nSamples; sInd++)
           fprintf(fp, "%e\n", sampleInputs[sInd*2+1]);
         fprintf(fp, "];\n");
         fprintf(fp, "zz = [\n");
         for (sInd = 0; sInd < nSamples; sInd++)
           fprintf(fp, "%e\n", faYIn[sInd]);
         fprintf(fp, "];\n");
       }
       fprintf(fp, "B = A;\n");
       double Ymin = faYOut[0];
       for (sInd = 1; sInd < faLeng; sInd++)
         if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
       double Ymax = faYOut[0];
       for (sInd = 1; sInd < faLeng; sInd++)
         if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];
       printf("Ymin and Ymax found = %e %e.\n", Ymin, Ymax);
       sprintf(pString,"Set lower threshold ? (y or n) : ");
       getString(pString, winput);
       fprintf(fp, "n1 = 0;\n");
       fprintf(fp, "n2 = 0;\n");
       double thresh;
       if (winput[0] == 'y')
       {
         sprintf(pString,"Enter the lower threshold (min = %e) : ",Ymin);
         thresh = getDouble(pString);
         fprintf(fp, "ymin = %e;\n", thresh);
         fprintf(fp, "[ia,ja,aa] = find(A<ymin);\n");
         fprintf(fp, "for ii = 1 : length(ia)\n");
         //fprintf(fp, "   B(ia(ii),ja(ii)) = %e;\n",Ymin-PABS(Ymin)*0.9);
         fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "n1 = length(ia);\n");
       }
       sprintf(pString,"Set upper threshold ? (y or n) : ");
       getString(pString, winput);
       if (winput[0] == 'y')
       {
         sprintf(pString,"Enter the upper threshold (max = %e) : ",Ymax);
         thresh = getDouble(pString);
         fprintf(fp, "ymax = %e;\n", thresh);
         fprintf(fp, "[ia,ja,aa] = find(A>ymax);\n");
         fprintf(fp, "for ii = 1 : length(ia)\n");
         //fprintf(fp, "   B(ia(ii),ja(ii)) = %e;\n",Ymin-PABS(Ymin)*0.9);
         fprintf(fp, "   B(ia(ii),ja(ii)) = NaN;\n");
         fprintf(fp, "end;\n");
         fprintf(fp, "n2 = length(ia);\n");
       }
       fprintf(fp, "nB = size(B,1);\n");
       fprintf(fp, "if (n1 + n2 == nB * nB)\n");
       fprintf(fp, "   B(1,1) = 0;\n");
       fprintf(fp, "   B(%d,%d) = 1;\n",nPtsPerDim,nPtsPerDim);
       fprintf(fp, "end\n");
       fprintf(fp, "if twoPlots == 1\n");
       fprintf(fp, "subplot(1,2,1), mesh(x,y,A)\n");
       if (nInputs == 2)
       {
         fwriteHold(fp,1);
         fprintf(fp,"plot3(xx,yy,zz,'k*','markersize',13);\n");
       } 
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, inputNames[iplot2]);
       fwritePlotZLabel(fp, outputNames[jplot]);
       sprintf(winput, "Mesh Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fprintf(fp,"colorbar\n");
       fprintf(fp,"subplot(1,2,2)\n");
       fprintf(fp,"end\n");
       fprintf(fp,"contourf(x,y,B)\n");
       if (nInputs == 2)
       {
         fwriteHold(fp,1);
         fprintf(fp,"plot(xx,yy,'k*','markersize',13);\n");
       } 
       fwritePlotAxes(fp);
       fwritePlotXLabel(fp, inputNames[iplot1]);
       fwritePlotYLabel(fp, inputNames[iplot2]);
       fprintf(fp,"colorbar\n");
       fprintf(fp,"colormap(jet)\n");
       sprintf(winput,"Contour Plot for %s", outputNames[jplot]);
       fwritePlotTitle(fp, winput);
       fclose(fp);
       printf("matlabrs2.m is now available for response surface and ");
       printf("contour plots\n");
     }
     delete [] faXOut;
     delete [] faYIn;
     delete [] faYOut;
     delete [] inputSettings;
     delete faPtr;
   }

   // +++ rs3 
   else if (!strcmp(command, "rs3"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs3: response surface plot in three parameters\n");
       printf("syntax: rs3 (no argument needed)\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data to analyze (load sample first).\n");
       return -1;
     }
     if (nInputs < 3)
     {
       printf("ERROR: rs3 requires 3 or more inputs.\n");
       return -1;
     }
     if (psPlotTool_ == 1)
     {
       printf("INFO: rs3 is currently not available in scilab.\n");
       return -1;
     }
     int nPtsPerDim = 16;
     sprintf(pString, "Grid resolution ? (16 - 32) ");
     nPtsPerDim = getInt(16, 32, pString);
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);

     int    iplot1, iplot2, iplot3;
     double *inputSettings = new double[nInputs];
     iplot1 = iplot2 = iplot3 = -1;
     sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
     iplot1 = getInt(1, nInputs, pString);
     iplot1--;
     iplot2 = iplot1;
     while (iplot1 == iplot2)
     {
       sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
               nInputs, iplot1+1);
       iplot2 = getInt(1, nInputs, pString);
       iplot2--;
       if (iplot1 == iplot2)
         printf("ERROR: duplicate input number %d.\n",iplot2+1);
     }
     if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
     while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
     {
       sprintf(pString,
               "Enter the input for z axis (1 - %d), not %d nor %d: ",
               nInputs, iplot1+1, iplot2+1);
       iplot3 = getInt(1, nInputs, pString);
       iplot3--;
       if (iplot3 == iplot1 || iplot3 == iplot2)
         printf("ERROR: duplicate input number %d.\n",iplot3+1);
     }
     strcpy(winput, "y");
     if (nInputs > 3)
     {
       sprintf(pString,"Set other nominal values automatically? (y/n) ");
       getString(pString, winput);
     }
     int iInd1, sInd, jplot;
     if (winput[0] == 'y')
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
         else inputSettings[iInd1] = 1.0;
       }
     }
     else
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
         {
           inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
           sprintf(pString,
                   "Enter nominal value for input %d (%e - %e): ", 
                   iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
           while (inputSettings[iInd1] < iLowerB[iInd1] ||
                  inputSettings[iInd1] > iUpperB[iInd1])
              inputSettings[iInd1] = getDouble(pString);
         }
         else inputSettings[iInd1] = 1.0;
       }
     }
     jplot = 0;
     sprintf(pString, "Enter the output number (1 - %d) : ", nOutputs);
     jplot = getInt(1, nOutputs, pString);
     jplot--;

     int    faLeng=0;
     double *faYIn = new double[nSamples];
     double *faXOut=NULL, *faYOut=NULL;
     for (sInd = 0; sInd < nSamples; sInd++)
       faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

     printf("Please wait while generating the RS data \n");
     faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                          inputSettings, &faLeng, &faXOut,&faYOut);

     double GYmin = faYOut[0];
     for (sInd = 1; sInd < faLeng; sInd++)
       if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
     double GYmax = faYOut[0];
     for (sInd = 1; sInd < faLeng; sInd++)
       if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];
     printf("\nYmin and Ymax found = %e %e.\n", GYmin, GYmax);
     double threshL = GYmin - 0.2 * PABS(GYmax-GYmin);
     double gamma = threshL;
     sprintf(pString,"Set lower threshold? (y or n) ");
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       sprintf(pString,"Enter the lower threshold (min = %e): ",GYmin);
       threshL = getDouble(pString);
       if (threshL < GYmin)
       {
         threshL = GYmin;
         printf("rs3 INFO: lower threshold set to %e.\n", threshL);
       }
     }
     int    ind;
     double threshU = GYmax + 0.2 * PABS(GYmax-GYmin);
     sprintf(pString,"Set upper threshold? (y or n) ");
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       sprintf(pString,"Enter the upper threshold (max = %e): ",GYmax);
       threshU = getDouble(pString);
       if (threshU > GYmax)
       {
         threshU = GYmax;
         printf("rs3 INFO: upper threshold set to %e.\n", threshU);
       }
     }
     if (threshL >= threshU)
     {
       printf("rs3 ERROR: lower threshold (%e) >= upper threshold (%e)\n",
              threshL, threshU);
       delete [] inputSettings;
       delete [] faYIn;
       delete [] faXOut;
       delete [] faYOut;
       delete faPtr;
       return -1;
     }
     fp = fopen("matlabrs3.m", "w");
     if (fp == NULL)
     {
       printf("ERROR: cannot open file matlabrs3.m.\n");
       delete [] inputSettings;
       delete [] faYIn;
       delete [] faXOut;
       delete [] faYOut;
       delete faPtr;
       return -1;
     }
     fwritePlotCLF(fp);
     fprintf(fp,"xlo = %e; \n", iLowerB[iplot2]);
     fprintf(fp,"xhi = %e; \n", iUpperB[iplot2]);
     fprintf(fp,"ylo = %e; \n", iLowerB[iplot1]);
     fprintf(fp,"yhi = %e; \n", iUpperB[iplot1]);
     fprintf(fp,"zlo = %e; \n", iLowerB[iplot3]);
     fprintf(fp,"zhi = %e; \n", iUpperB[iplot3]);
     fprintf(fp,"X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp,"Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp,"Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp,"V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     for (jj = 0; jj < nPtsPerDim; jj++)
     {
       fprintf(fp,"Y(:,:,%d) = [\n", jj + 1);
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
       {
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
           ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
           fprintf(fp,"%e ", faXOut[ind*3]);
         }
         fprintf(fp,"\n");
       }
       fprintf(fp, "];\n");
       fprintf(fp, "X(:,:,%d) = [\n", jj + 1);
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
       {
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
           ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
           fprintf(fp, "%e ", faXOut[ind*3+1]);
         }
         fprintf(fp, "\n");
       }
       fprintf(fp, "];\n");
       fprintf(fp, "Z(:,:,%d) = [\n", jj + 1);
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
       {
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
           ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
           fprintf(fp, "%e ", faXOut[ind*3+2]);
         }
         fprintf(fp, "\n");
       }
       fprintf(fp, "];\n");
     }
     int count=0;
     for (jj = 0; jj < nPtsPerDim; jj++)
     {
       fprintf(fp, "V(:,:,%d) = [\n", jj + 1);
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
       {
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
           ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
           if (faYOut[ind] < threshL)
           {
             fprintf(fp, "%e ", gamma);
             count++;
           }
           else if (faYOut[ind] > threshU)
           {
             fprintf(fp, "%e ", gamma);
             count++;
           }
           else fprintf(fp, "%e ", faYOut[ind]);
         }
         fprintf(fp, "\n");
       }
       fprintf(fp, "];\n");
     }
     if (count == nPtsPerDim*nPtsPerDim*nPtsPerDim)
     {
       fprintf(fp, "V(1,1,1)=0;\n");
       fprintf(fp, "V(%d,%d,%d)=1;\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     }
     fprintf(fp,"xt = [%e:%e:%e];\n", iLowerB[iplot2],
             (iUpperB[iplot2]-iLowerB[iplot2])*0.01, iUpperB[iplot2]);
     fprintf(fp,"yt = [%e:%e:%e];\n", iLowerB[iplot1],
             (iUpperB[iplot1]-iLowerB[iplot1])*0.01, iUpperB[iplot1]);
     fprintf(fp,"zt = [%e:%e:%e];\n", iLowerB[iplot3],
             (iUpperB[iplot3]-iLowerB[iplot3])*0.01, iUpperB[iplot3]);
     fprintf(fp,"isoval = %e;\n", gamma);
     fprintf(fp,"h = patch(isosurface(X,Y,Z,V,isoval),... \n");
     fprintf(fp,"          'FaceColor', 'blue', ... \n");
     fprintf(fp,"          'EdgeColor', 'none', ... \n");
     fprintf(fp,"          'AmbientStrength', 0.2, ... \n");
     fprintf(fp,"          'SpecularStrength', 0.7, ... \n");
     fprintf(fp,"          'DiffuseStrength', 0.4);\n");
     fprintf(fp,"isonormals(X,Y,Z,V,h);\n");
     fprintf(fp,"patch(isocaps(X,Y,Z,V,isoval), ...\n");
     fprintf(fp,"      'FaceColor', 'interp', ... \n");
     fprintf(fp,"      'EdgeColor', 'none'); \n");
     fprintf(fp,"axis([xlo xhi ylo yhi zlo zhi])\n");
     fprintf(fp,"daspect([%e,%e,%e])\n",iUpperB[iplot2]-iLowerB[iplot2],
             iUpperB[iplot1]-iLowerB[iplot1],
             iUpperB[iplot3]-iLowerB[iplot3]);
     fprintf(fp,"   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
             inputNames[iplot2]);
     fprintf(fp,"   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot1]);
     fprintf(fp,"   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot3]);
     fprintf(fp,"   title('%s','Fontsize',12,'FontWeight','bold')\n", 
             outputNames[jplot]);
     fwritePlotAxes(fp);
     fprintf(fp,"colormap('default'); colorbar\n");
     fprintf(fp,"%%axis tight\n");
     fprintf(fp,"view(3) \n");
     fprintf(fp,"set(gcf,'Renderer','zbuffer')\n");
     fprintf(fp,"lighting phong\n");
     fprintf(fp,"cin = input('generate slices ? (y or n) ','s');\n");
     fprintf(fp,"if (cin == 'y')\n");
     fprintf(fp,"xin = input('axis to slide through ? (x,y,z) ','s');\n");
     fprintf(fp,"for i = 1 : 101\n");
     fprintf(fp,"   display(['displaying ' int2str(i) ' of 100'])\n");
     fprintf(fp,"   if (xin == 'y')\n");
     fprintf(fp,"      h = slice(X,Y,Z,V,xt(i),[],[]);\n");
     fprintf(fp,"   elseif (xin == 'x')\n");
     fprintf(fp,"      h = slice(X,Y,Z,V,[],yt(i),[]);\n");
     fprintf(fp,"   elseif (xin == 'z')\n");
     fprintf(fp,"      h = slice(X,Y,Z,V,[],[],zt(i));\n");
     fprintf(fp,"   end\n");
     fprintf(fp,"   axis([%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e ",
             iLowerB[iplot2], iUpperB[iplot2], iLowerB[iplot1],
             iUpperB[iplot1], iLowerB[iplot3], iUpperB[iplot3]);
     fprintf(fp,"%11.4e %11.4e])\n",
             threshL-0.2*(threshU-threshL),threshU+0.2*(threshU-threshL));
     fwritePlotAxes(fp);
     fprintf(fp,"   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
             inputNames[iplot2]);
     fprintf(fp,"   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot1]);
     fprintf(fp,"   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot3]);
     fprintf(fp,"   title('3D Contour Plot',");
     fprintf(fp,"'FontWeight','bold','FontSize',12)\n");
     fprintf(fp,"   view(3)\n");
     fprintf(fp,"   colorbar\n");
     fprintf(fp,"   pause(1)\n");
     fprintf(fp,"   if (i < 101)\n");
     fprintf(fp,"      clf\n");
     fprintf(fp,"   end\n");
     fprintf(fp,"end\n");
     fprintf(fp,"end\n");
     fclose(fp);
     printf("\nmatlabrs3.m is now available.\n");
     delete [] inputSettings;
     delete [] faYIn;
     delete [] faXOut;
     delete [] faYOut;
     delete faPtr;
   }

   // +++ rs3m 
   else if (!strcmp(command, "rs3m"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs3m: response surface plot in 3 parameters using movie mode\n");
       printf("syntax: rs3m (no argument needed)\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data to analyze (load sample first).\n");
       return -1;
     }
     if (nInputs < 3)
     {
       printf("ERROR: rs3m requires 3 or more inputs.\n");
       return -1;
     }
     if (psPlotTool_ == 1)
     {
       printf("INFO: rs3m is currently not available in scilab.\n");
       return -1;
     }
     int nPtsPerDim = 16;
     sprintf(pString, "Grid resolution ? (16 - 32) ");
     nPtsPerDim = getInt(16, 32, pString);
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);

     int    iplot1, iplot2, iplot3, jplot, iInd1, sInd;
     double *inputSettings = new double[nInputs];
     iplot1 = iplot2 = iplot3 = -1;
     sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
     iplot1 = getInt(1, nInputs, pString);
     iplot1--;
     iplot2 = iplot1;
     while (iplot1 == iplot2)
     {
       sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
               nInputs, iplot1+1);
       iplot2 = getInt(1, nInputs, pString);
       iplot2--;
       if (iplot1 == iplot2)
          printf("ERROR: duplicate input number.\n");
     }
     if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
     while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
     {
       sprintf(pString,
               "Enter the input for t axis (1 - %d), not %d nor %d: ",
               nInputs, iplot1+1, iplot2+1);
       iplot3 = getInt(1, nInputs, pString);
       iplot3--;
       if (iplot3 == iplot1 || iplot3 == iplot2)
         printf("ERROR: duplicate input number %d.\n", iplot3+1);
     }
     if (nInputs > 3)
     {
       sprintf(pString,
               "Set other nominal values at mid point ? (y or n) ");
       getString(pString, winput);
       if (winput[0] == 'y')
       {
         for (iInd1 = 0; iInd1 < nInputs; iInd1++)
         {
           if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
           else inputSettings[iInd1] = 1.0;
         }
       }
       else
       {
         for (iInd1 = 0; iInd1 < nInputs; iInd1++)
         {
           if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
           {
             sprintf(pString,"Enter setting for input %d (%e - %e): ", 
                     iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
             inputSettings[iInd1] = getDouble(pString);
           }
           else inputSettings[iInd1] = 1.0;
         }
       }
     }
     jplot = 0;
     sprintf(pString,"Enter the output number (1 - %d) : ", nOutputs);
     jplot = getInt(1, nOutputs, pString);
     jplot--;

     double *faYIn = new double[nSamples];
     for (sInd = 0; sInd < nSamples; sInd++)
       faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

     fp = fopen("matlabrs3m.m", "w");
     if (fp == NULL)
     {
       printf("ERROR: cannot open file matlabrs3m.m.\n");
       delete faPtr;
       delete [] inputSettings;
       delete [] faYIn;
       return -1;
     }

     int faLeng=0;
     double *faXOut=NULL, *faYOut=NULL;
     faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                          inputSettings, &faLeng, &faXOut,&faYOut);
     double GYmin = faYOut[0];
     for (sInd = 1; sInd < faLeng; sInd++)
       if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
     double GYmax = faYOut[0];
     for (sInd = 1; sInd < faLeng; sInd++)
       if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];
     printf("\nYmin and Ymax found = %e %e.\n", GYmin, GYmax);
     double threshL = GYmin - 0.2 * PABS(GYmin);
     sprintf(pString, "Set lower threshold? (y or n) ");
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       sprintf(pString,"Enter the lower threshold (min = %e): ",GYmin);
       threshL = getDouble(pString);
     }
     double threshU = GYmax + 0.2 * PABS(GYmax);
     sprintf(pString, "Set upper threshold? (y or n) ");
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       sprintf(pString,"Enter the upper threshold (min = %e): ",GYmax);
       threshU = getDouble(pString);
     }
     fprintf(fp,"twoPlots = 1;\n");
     fprintf(fp,"disp(\'Please wait while loading data.\')\n");
     fprintf(fp,"hold off\n");
     fwritePlotCLF(fp);

     for (ii = 0; ii < nPtsPerDim; ii++)
     {
       inputSettings[iplot3] = (iUpperB[iplot3] - iLowerB[iplot3]) / 
                               (nPtsPerDim - 1.0) * ii + iLowerB[iplot3];
       fprintf(fp,"x = [\n");
       for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim*nPtsPerDim)
         fprintf(fp, "%e\n", faXOut[sInd*3]);
       fprintf(fp,"];\n");
       fprintf(fp,"y = [\n");
       for (sInd = 0; sInd < nPtsPerDim*nPtsPerDim; sInd+=nPtsPerDim)
         fprintf(fp, "%e\n", faXOut[sInd*3+1]);
       fprintf(fp,"];\n");

       fprintf(fp,"A%d = [\n", ii + 1);
       for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
         fprintf(fp, "%e\n", faYOut[sInd+ii]);
       fprintf(fp,"];\n");
       fprintf(fp,"A%d = reshape(A%d,%d,%d);\n", ii+1, ii+1,
               nPtsPerDim, nPtsPerDim);
       fprintf(fp,"disp(\'Plotting frame %d of %d\')\n",ii+1,nPtsPerDim);
       fprintf(fp,"B%d = A%d;\n", ii+1, ii+1);
       fprintf(fp,"yLo = %e;\n", threshL);
       fprintf(fp,"yHi = %e;\n", threshU);
       fprintf(fp,"nA  = size(A%d,1);\n", ii+1);
       fprintf(fp,"[ia,ja,aa] = find(A%d<yLo);\n", ii+1);
       fprintf(fp,"for ii = 1 : length(ia)\n");
       fprintf(fp,"   B%d(ia(ii),ja(ii)) = NaN;\n", ii+1);
       fprintf(fp,"end;\n");
       fprintf(fp,"n1 = length(ia);\n");
       fprintf(fp,"[ia,ja,aa] = find(A%d>yHi);\n", ii+1);
       fprintf(fp,"for ii = 1 : length(ia)\n");
       fprintf(fp,"   B%d(ia(ii),ja(ii)) = NaN;\n", ii+1);
       fprintf(fp,"end;\n");
       fprintf(fp,"n2 = length(ia);\n");
       fprintf(fp,"if (n1 + n2 == nA*nA)\n");
       fprintf(fp,"   B%d(1,1) = 0;\n",ii+1);
       fprintf(fp,"   B%d(%d,%d) = 1;\n",ii+1,nPtsPerDim,nPtsPerDim);
       fprintf(fp,"end;\n");
       fprintf(fp,"if twoPlots == 1\n");
       fprintf(fp,"subplot(1,2,1), surf(x,y,A%d)\n", ii+1);
       fprintf(fp,"axis([%e %e %e %e %e %e])\n",iLowerB[iplot1],
               iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2],GYmin, GYmax); 
       fwritePlotAxes(fp);
       fprintf(fp,"xlabel('%s','FontSize',12,'FontWeight','bold')\n",
               inputNames[iplot1]);
       fprintf(fp,"ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
               inputNames[iplot2]);
       fprintf(fp,"zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
               outputNames[jplot]);
       fprintf(fp,"colorbar\n");
       fprintf(fp,"title(\'%s Mesh plot, val(3) = %14.7e\',",
               outputNames[jplot], inputSettings[iplot3]);
       fprintf(fp,"'FontWeight','bold','FontSize',12)\n");
       fprintf(fp,"subplot(1,2,2)\n");
       fprintf(fp,"end\n");
       fprintf(fp,"contourf(x,y,B%d)\n",ii+1);
       fprintf(fp,"axis([%e %e %e %e])\n",iLowerB[iplot1],
               iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
       fwritePlotAxes(fp);
       fprintf(fp,"colorbar\n");
       fprintf(fp,"colormap(jet)\n");
       fprintf(fp,"caxis([%e %e])\n",GYmin, GYmax);
       fprintf(fp,"title(\'%s contour plot, val(3) = %14.7e\',",
               outputNames[jplot], inputSettings[iplot3]);
       fprintf(fp,"'FontWeight','bold','FontSize',12)\n");
       fprintf(fp,"pause(1)\n");
     }
     fprintf(fp,"rotate3d on\n");
     fclose(fp);
     printf("matlabrs3m.m is now available for response surface and ");
     printf("contour plots\n");
     delete [] faXOut;
     delete [] faYOut;
     delete [] faYIn;
     delete [] inputSettings;
     delete faPtr;
   }

   // +++ rs4 
   else if (!strcmp(command, "rs4"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rs4: response surface plot in 4 parameters\n");
       printf("syntax: rs4 (no argument needed)\n");
       return 0;
     }
     if (nInputs <= 0 || psuadeIO == NULL)
     {
       printf("ERROR: data not loaded yet.\n");
       return -1;
     }
     if (nInputs < 4)
     {
       printf("ERROR: rs4 requires 4 or more inputs.\n");
       return -1;
     }
     if (psPlotTool_ == 1)
     {
       printf("INFO: rs4 is currently not available in scilab.\n");
       return -1;
     }
     int nPtsPerDim = 16;
     printf("NOTE: if matlab crashes, it may be due to high grid resolution\n");
     sprintf(pString, "Grid resolution ? (16 - 32) ");
     nPtsPerDim = getInt(16, 32, pString);
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);

     int    iplot1, iplot2, iplot3, iplot4, iInd1;
     double *inputSettings = new double[nInputs];
     iplot1 = iplot2 = iplot3 = iplot4 = -1;
     sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
     iplot1 = getInt(1, nInputs, pString);
     iplot1--;
     iplot2 = iplot1;
     while (iplot1 == iplot2)
     {
       sprintf(pString, "Enter the input for y axis (1 - %d), not %d : ",
               nInputs, iplot1+1);
       iplot2 = getInt(1, nInputs, pString);
       iplot2--;
       if (iplot1 == iplot2)
          printf("ERROR: duplicate input number %d.\n",iplot2+1);
     }
     iplot3 = iplot1;
     while (iplot3 == iplot1 || iplot3 == iplot2)
     {
       sprintf(pString, "Enter the input for z axis (1 - %d), not %d,%d: ",
               nInputs, iplot1+1, iplot2+1);
       iplot3 = getInt(1, nInputs, pString);
       iplot3--;
       if (iplot3 == iplot1 || iplot3 == iplot2)
         printf("ERROR: duplicate input number %d.\n",iplot3+1);
     }
     if (nInputs == 4) iplot4 = 6 - iplot1 - iplot2 - iplot3;
     while (iplot4 < 0 || iplot4 == iplot1 || iplot4 == iplot2 || 
            iplot4 == iplot3)
     {
       sprintf(pString,
               "Enter the input for t axis (1 - %d), not %d nor %d,%d: ",
               nInputs, iplot1+1, iplot2+1, iplot3+1);
       iplot4 = getInt(1, nInputs, pString);
       iplot4--;
       if (iplot4 == iplot1 || iplot4 == iplot2 || iplot4 == iplot3)
         printf("ERROR: duplicate input number %d.\n",iplot4+1);
     }
     strcpy(winput, "y");
     if (nInputs > 4)
     {
       sprintf(pString,"Set other nominal values at mid point ? (y/n) ");
       getString(pString, winput);
     }
     if (winput[0] == 'y')
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
             iInd1 != iplot4)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
         else inputSettings[iInd1] = 1.0;
       }
     }
     else
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
             iInd1 != iplot4)
         {
            inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
            sprintf(pString,
                    "Enter nominal value for input %d (%e - %e): ", 
                    iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
            while (inputSettings[iInd1] < iLowerB[iInd1] ||
                   inputSettings[iInd1] > iUpperB[iInd1])
              inputSettings[iInd1] = getDouble(pString);
         }
         else inputSettings[iInd1] = 1.0;
       }
     }
     int jplot=0, sInd;
     sprintf(pString, "Enter the output number (1 - %d) : ", nOutputs);
     jplot = getInt(1, nOutputs, pString);
     jplot--;

     double *faYIn = new double[nSamples];
     for (sInd = 0; sInd < nSamples; sInd++)
       faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

     int    faLeng=0;
     double *faXOut=NULL, *faYOut=NULL;
     faPtr->gen4DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                          iplot4, inputSettings, &faLeng, &faXOut,&faYOut);
     double GYmin =   PSUADE_UNDEFINED;
     double GYmax = - PSUADE_UNDEFINED;
     for (sInd = 0; sInd < faLeng; sInd++)
       if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
     for (sInd = 0; sInd < faLeng; sInd++)
     if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];
     printf("\nYmin and Ymax found = %e %e.\n", GYmin, GYmax);
     double threshL = GYmin - 0.2 * PABS(GYmax - GYmin);
     sprintf(pString,"Set lower threshold? (y or n) ");
     double gamma = threshL;
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       sprintf(pString,"Enter the lower threshold (min = %e): ",GYmin);
       threshL = getDouble(pString);
     }
     double threshU = GYmax + 0.2 * PABS(GYmax - GYmin);
     sprintf(pString,"Set upper threshold? (y or n) ");
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       sprintf(pString,"Enter the upper threshold (max = %e): ",GYmax);
       threshU = getDouble(pString);
     }

     fp = fopen("matlabrs4.m", "w");
     if (fp == NULL)
     {
       printf("ERROR: cannot open file matlabrs4.m.\n");
       if (inputSettings != NULL) delete [] inputSettings;
       if (faYIn  != NULL) delete [] faYIn;
       if (faXOut != NULL) delete [] faXOut;
       if (faYOut != NULL) delete [] faYOut;
       if (faPtr  != NULL) delete faPtr;
       return -1;
     }

     fprintf(fp,"%% user adjustable parameter section begins *****\n");
     fprintf(fp,"%% use nSubplots, nSubNx and nSubNy to spread \n");
     fprintf(fp,"%% the movie frames into a number of subplots.\n");
     fprintf(fp,"nSubplots = 1;\n");
     fprintf(fp,"nSubNx = 1;\n");
     fprintf(fp,"nSubNy = 1;\n");
     fprintf(fp,"%% user adjustable parameter section ends *****\n");
     fwritePlotCLF(fp);
     fprintf(fp,"nFrames = %d;\n", nPtsPerDim);
     fprintf(fp,"nSubCnt = 0;\n");
     fprintf(fp,"isoval = %e;\n", threshL);
     fprintf(fp,"xlo = %e; \n", iLowerB[iplot2]);
     fprintf(fp,"xhi = %e; \n", iUpperB[iplot2]);
     fprintf(fp,"ylo = %e; \n", iLowerB[iplot1]);
     fprintf(fp,"yhi = %e; \n", iUpperB[iplot1]);
     fprintf(fp,"zlo = %e; \n", iLowerB[iplot3]);
     fprintf(fp,"zhi = %e; \n", iUpperB[iplot3]);
     fprintf(fp,"X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp,"Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp,"Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp,"V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     int ind;
     for (ii = 0; ii < nPtsPerDim; ii++)
     {
       fprintf(fp,"Y(:,:,%d) = [\n", ii + 1);
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
       {
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
           ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*nPtsPerDim;
           fprintf(fp, "%e ", faXOut[ind*4]);
         }
         fprintf(fp, "\n");
       }
       fprintf(fp,"];\n");
       fprintf(fp,"X(:,:,%d) = [\n", ii + 1);
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
       {
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
           ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*nPtsPerDim;
           fprintf(fp, "%e ", faXOut[ind*4+1]);
         }
         fprintf(fp, "\n");
       }
       fprintf(fp,"];\n");
       fprintf(fp,"Z(:,:,%d) = [\n", ii + 1);
       for (sInd = 0; sInd < nPtsPerDim; sInd++)
       {
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
           ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*nPtsPerDim;
           fprintf(fp,"%e ", faXOut[ind*4+2]);
         }
         fprintf(fp, "\n");
       }
       fprintf(fp, "];\n");
     }
     int count, ll;
     for (ll = 0; ll < nPtsPerDim; ll++)
     {
       for (ii = 0; ii < nPtsPerDim; ii++)
       {
         count = 0;
         fprintf(fp,"V(:,:,%d) = [\n", ii + 1);
         for (sInd = 0; sInd < nPtsPerDim; sInd++)
         {
           for (jj = 0; jj < nPtsPerDim; jj++)
           {
             ind = ((sInd*nPtsPerDim+jj)*nPtsPerDim+ii)*nPtsPerDim+ll;
             if (faYOut[ind] < threshL)
             {
               fprintf(fp, "%e ", threshL);
               count++;
             }
             else if (faYOut[ind] > threshU)
             {
               fprintf(fp, "%e ", threshL);
               count++;
             }
             else fprintf(fp, "%e ", faYOut[ind]);
           }
           fprintf(fp, "\n");
         }
         fprintf(fp,"];\n");
         if (count == nPtsPerDim*nPtsPerDim)
         {
           if (threshL-0.2*(threshU-threshL) > gamma)
              fprintf(fp,"V(:,:,%d) = %e * ones(%d,%d);\n",ii+1,gamma,
                      nPtsPerDim, nPtsPerDim);
           else
              fprintf(fp,"V(:,:,%d) = %e * ones(%d,%d);\n",ii+1,
                      threshL-0.2*(threshU-threshL),
                      nPtsPerDim,nPtsPerDim);
           printf("Frame %d, slice %d nonfeasible -> set to ground.\n",
                  ll+1, ii+1);
         }
       }
       fprintf(fp,"frame = %d;\n", ll+1);
       fprintf(fp,"if nSubplots > 1\n");
       fprintf(fp,"   if frame <= 2\n");
       fprintf(fp,"      nSubCnt = nSubCnt + 1;\n");
       fprintf(fp,"      subplot(nSubNx, nSubNy, nSubCnt)\n");
       fprintf(fp,"   elseif frame == nFrames\n");
       fprintf(fp,"      subplot(nSubNx, nSubNy, nSubplots)\n");
       fprintf(fp,"   else\n");
       fprintf(fp,"      ft1 = (nFrames-1) / (nSubplots-1);\n");
       fprintf(fp,"      ft2 = round(ft1 * (nSubCnt-1)) + 2;\n");
       fprintf(fp,"      if frame == ft2\n");
       fprintf(fp,"         nSubCnt = nSubCnt + 1;\n");
       fprintf(fp,"         subplot(nSubNx, nSubNy, nSubCnt)\n");
       fprintf(fp,"      end\n");
       fprintf(fp,"   end\n");
       fprintf(fp,"else\n");
       fprintf(fp,"   clf\n");
       fprintf(fp,"end\n");
       fprintf(fp,"disp('Frame %d of %d')\n", ll+1, nPtsPerDim);
       fprintf(fp,"h = patch(isosurface(X,Y,Z,V,isoval),... \n");
       fprintf(fp,"          'FaceColor', 'blue', ... \n");
       fprintf(fp,"          'EdgeColor', 'none', ... \n");
       fprintf(fp,"          'AmbientStrength', 0.2, ... \n");
       fprintf(fp,"          'SpecularStrength', 0.7, ... \n");
       fprintf(fp,"          'DiffuseStrength', 0.4);\n");
       fprintf(fp,"isonormals(X,Y,Z,V,h);\n");
       fprintf(fp,"patch(isocaps(X,Y,Z,V,isoval), ...\n");
       fprintf(fp,"      'FaceColor', 'interp', ... \n");
       fprintf(fp,"      'EdgeColor', 'none'); \n");
       fprintf(fp,"axis([xlo xhi ylo yhi zlo zhi])\n");
       fprintf(fp,"daspect([xhi-xlo, yhi-ylo, zhi-zlo])\n");
       fprintf(fp,"colormap('default')\n");
       fprintf(fp,"if nSubplots == 1\n");
       fprintf(fp,"   colorbar\n");
       fprintf(fp,"end\n");
       fprintf(fp,"%%axis tight\n");
       fprintf(fp,"view(3) \n");
       fprintf(fp,"set(gcf,'Renderer','zbuffer')\n");
       fprintf(fp,"box on\n");
       fprintf(fp,"grid on\n");
       fprintf(fp,"lighting phong\n");
       fwritePlotAxes(fp);
       if (ll == 0)
       {
         fprintf(fp,"xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot2]);
         fprintf(fp,"ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp,"zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
                 inputNames[iplot3]);
       }
       fprintf(fp,"title('%s=%12.4e',",inputNames[iplot4],faXOut[ll*4+3]);
       fprintf(fp,"'FontWeight','bold','FontSize',12)\n");
       fprintf(fp,"pause(1)\n");
     }
     fclose(fp);
     printf("\nmatlabrs4.m is now available.\n");
     if (faXOut != NULL) delete [] faXOut;
     if (faYOut != NULL) delete [] faYOut;
     if (faYIn  != NULL) delete [] faYIn;
     if (inputSettings != NULL) delete [] inputSettings;
     if (faPtr != NULL) delete faPtr;
   }

   // +++ rssd 
   else if (!strcmp(command, "rssd"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rssd: response surface plots for the std. deviations.\n");
         printf("INFO: rssd not available for >2 inputs for scilab.\n");
         printf("syntax: rssd (no argument needed.\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data to analyze (load sample first).\n");
         return -1;
      }
      int    iplot1, iplot2, iplot3, iplot4, count, iInd1, sInd, ind, ll;
      double *inputSettings = new double[nInputs];
      iplot1 = iplot2 = iplot3 = iplot4 = -1;
      sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
      iplot1 = getInt(1, nInputs, pString);
      iplot1--;
      count = 1;
      if (nInputs == 2) iplot2 = nInputs - iplot1 - 1;
      else if (nInputs > 2)
      {
         iplot2 = iplot1;
         while (iplot1 == iplot2)
         {
            sprintf(pString, "Y-axis input ? (1-%d, 0 if not used, not %d) ",
                    nInputs, iplot1+1);
            iplot2 = getInt(0, nInputs, pString);
            iplot2--;
            if (iplot2 == -1) break;
            if (iplot1 == iplot2)
               printf("ERROR: duplicate input number %d.\n",iplot2+1);
         }
      }
      if (iplot2 != -1) count++;
      if (psPlotTool_ == 0 && iplot2 != -1)
      {
         if (nInputs == 3) iplot3 = nInputs - iplot1 - iplot2;
         else if (nInputs > 3)
         {
            iplot3 = iplot1;
            while (iplot3 == iplot1 || iplot3 == iplot2)
            {
               sprintf(pString,
                  "Z axis input ? (1-%d, 0 if not used, not %d nor %d) ",
                  nInputs, iplot1+1, iplot2+1);
               iplot3 = getInt(0, nInputs, pString);
               iplot3--;
               if (iplot3 == -1) break;
               if (iplot3 == iplot1 || iplot3 == iplot2)
                  printf("ERROR: duplicate input number %d.\n",iplot3+1);
            }
         }
         if (iplot3 != -1) count++;
         if (nInputs >= 4 && iplot3 != -1)
         {
            while (iplot4 < 0 || iplot4 == iplot1 || iplot4 == iplot2 || 
                   iplot4 == iplot3)
            {
               sprintf(pString,
                  "Enter the input for t axis (1 - %d), not %d nor %d,%d: ",
                       nInputs, iplot1+1, iplot2+1, iplot3+1);
               iplot4 = getInt(1, nInputs, pString);
               iplot4--;
               if (iplot4 == iplot1 || iplot4 == iplot2 || iplot4 == iplot3)
                  printf("ERROR: duplicate input number %d.\n",iplot4+1);
            }
         }
         if (iplot4 != -1) count++;
      }
      strcpy(winput, "y");
      if (nInputs > count)
      {
         sprintf(pString,"Set other nominal values at mid point ? (y/n) ");
         getString(pString, winput);
      }
      if (winput[0] == 'y')
      {
         for (iInd1 = 0; iInd1 < nInputs; iInd1++)
         {
            if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
                iInd1 != iplot4)
                 inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
            else inputSettings[iInd1] = 1.0;
         }
      }
      else
      {
         for (iInd1 = 0; iInd1 < nInputs; iInd1++)
         {
            if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3 &&
                iInd1 != iplot4)
            {
               inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
               sprintf(pString,
                       "Enter nominal value for input %d (%e - %e): ", 
                       iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
               while (inputSettings[iInd1] < iLowerB[iInd1] ||
                      inputSettings[iInd1] > iUpperB[iInd1])
                  inputSettings[iInd1] = getDouble(pString);
            }
            else inputSettings[iInd1] = 1.0;
         }
      }
      int nPtsPerDim;
      if      (iplot2 == -1) nPtsPerDim = 1024;
      else if (iplot3 == -1) nPtsPerDim = 128;
      else if (iplot4 == -1) nPtsPerDim = 24;
      else                   nPtsPerDim = 10;
      printf("This command works with the following response surfaces:\n");
      printf("1. Linear    regression\n");
      printf("2. Quadratic regression\n");
      printf("3. cubic     regression\n");
      printf("4. quartic   regression\n");
      printf("5. GP1\n");
      printf("6. GP2\n");
      printf("7. MarsBagg\n");
      printf("8. Tree GP\n");
      printf("9. Kriging\n");
      sprintf(pString, "Enter your choice: (1, 2, ..., 9) ");
      int faType = getInt(1, 9, pString);
      if      (faType == 1) faType = PSUADE_RS_REGR1;
      else if (faType == 2) faType = PSUADE_RS_REGR2;
      else if (faType == 3) faType = PSUADE_RS_REGR3;
      else if (faType == 4) faType = PSUADE_RS_REGR4;
      else if (faType == 5) faType = PSUADE_RS_GP1;
      else if (faType == 6) faType = PSUADE_RS_GP2;
      else if (faType == 7) faType = PSUADE_RS_MARSB;
      else if (faType == 8) faType = PSUADE_RS_TGP;
      else if (faType == 9) faType = PSUADE_RS_KR;
      int faFlag = 1, iOne=1;
      FuncApprox *faPtr = genFA(faType, nInputs, iOne, nSamples);
      if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(iLowerB, iUpperB);
      faPtr->setOutputLevel(outputLevel);

      int jplot = 0;
      sprintf(pString, "Enter the output number (1 - %d) : ",nOutputs);
      jplot = getInt(1, nOutputs, pString);
      jplot--;

      double *faYIn = new double[nSamples];
      for (sInd = 0; sInd < nSamples; sInd++)
         faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

      int faLeng = 0;
      double *faXOut=NULL, *faYOut=NULL;
      if (iplot2 == -1)
         faPtr->gen1DGridData(sampleInputs,faYIn,iplot1,
                              inputSettings, &faLeng, &faXOut,&faYOut);
      else if (iplot3 == -1)
         faPtr->gen2DGridData(sampleInputs,faYIn,iplot1,iplot2,
                             inputSettings, &faLeng, &faXOut,&faYOut);
      else if (iplot4 == -1)
         faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                              inputSettings, &faLeng, &faXOut,&faYOut);
      else
         faPtr->gen4DGridData(sampleInputs,faYIn, iplot1, iplot2, 
                              iplot3, iplot4, inputSettings, &faLeng, 
                              &faXOut,&faYOut);

      double *tempW = new double[faLeng];
      double *tempX = new double[faLeng*nInputs];
      for (sInd = 0; sInd < faLeng; sInd++)
         for (jj = 0; jj < nInputs; jj++)
            tempX[sInd*nInputs+jj] = inputSettings[jj];
      for (sInd = 0; sInd < faLeng; sInd++)
      {
         tempX[sInd*nInputs+iplot1] = faXOut[sInd*count];
         if (iplot2 != -1)
            tempX[sInd*nInputs+iplot2] = faXOut[sInd*count+1];
         if (iplot3 != -1)
            tempX[sInd*nInputs+iplot3] = faXOut[sInd*count+2];
         if (iplot4 != -1)
            tempX[sInd*nInputs+iplot4] = faXOut[sInd*count+3];
      }
      faPtr->evaluatePointFuzzy(faLeng, tempX, faYOut, tempW);
      double gamma = PSUADE_UNDEFINED;
      for (sInd = 0; sInd < faLeng; sInd++)
         if (tempW[sInd] < gamma) gamma = tempW[sInd];
         
      if (psPlotTool_ == 1)
      {
         fp = fopen("scilabrssd.sci", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file scilabrssd.sci.\n");
            delete [] faXOut;
            delete [] faYOut;
            delete [] faYIn;
            delete [] tempW;
            delete [] tempX;
            delete faPtr;
            delete [] inputSettings;
            return -1;
         }
      }
      else
      {
         fp = fopen("matlabrssd.m", "w");
         if (fp == NULL)
         {
            printf("ERROR: cannot open file matlabrssd.m.\n");
            delete [] faXOut;
            delete [] faYOut;
            delete [] faYIn;
            delete [] tempW;
            delete [] tempX;
            delete faPtr;
            delete [] inputSettings;
            return -1;
         }
      }
      fwritePlotCLF(fp);
      if (count == 1)
      {
         fprintf(fp, "A = [\n");
         for (sInd = 0; sInd < faLeng; sInd++)
            fprintf(fp, "%e\n", tempW[sInd]);
         fprintf(fp, "];\n");
         fprintf(fp, "X = [\n");
         for (sInd = 0; sInd < faLeng; sInd++)
            fprintf(fp, "%e\n", faXOut[sInd]);
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "plot(X,A);");
            fprintf(fp, "a = gca();\n");
            fprintf(fp, "a.children.children.thickness = 4;\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
         }
         else 
         {
            fprintf(fp, "plot(X,A,'lineWidth',4)\n");
            fprintf(fp, "hold on\n");
         }
         fwritePlotAxes(fp);
         fwritePlotXLabel(fp, inputNames[iplot1]);
         fwritePlotYLabel(fp, outputNames[jplot]);
         sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
         fwritePlotTitle(fp, winput);
      }
      else if (count == 2)
      {
         if (psPlotTool_ == 0) fprintf(fp, "twoPlots = 1;\n");
         fprintf(fp, "A = [\n");
         for (sInd = 0; sInd < faLeng; sInd++)
            fprintf(fp, "%e\n", tempW[sInd]);
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1) 
            fprintf(fp, "A = matrix(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
         else
            fprintf(fp, "A = reshape(A,%d,%d);\n", nPtsPerDim, nPtsPerDim);
         fprintf(fp, "X = [\n");
         for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
            fprintf(fp, "%e\n", faXOut[sInd*2]);
         fprintf(fp, "];\n");
         fprintf(fp, "Y = [\n");
         for (sInd = 0; sInd < nPtsPerDim; sInd++)
            fprintf(fp, "%e\n", faXOut[sInd*2+1]);
         fprintf(fp, "];\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "mesh(X,Y,A)\n");
            fprintf(fp, "h = get(\"hdl\");\n");
            fprintf(fp, "h.color_flag=1;\n");
            fprintf(fp, "h.color_mode=-2;\n");
            fprintf(fp, "bmin = min(min(A)); bmax = max(max(A));\n");
            fprintf(fp, "xset(\"colormap\",jetcolormap(64));\n");
            fprintf(fp, "colorbar(bmin,bmax);\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fwritePlotZLabel(fp, outputNames[jplot]);
            sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp, "scf(2);\n");
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[%e,%e;%e,%e];\n",
                    iLowerB[iplot1], iLowerB[iplot2],
                    iUpperB[iplot1], iUpperB[iplot2]);
            fprintf(fp, "a.axes_visible=\"on\";\n");
            fprintf(fp, "B = A;\n");
            fprintf(fp, "nX = length(X);\n");
            fprintf(fp, "nY = length(Y);\n");
            fprintf(fp, "for ii = 1 : nX\n");
            fprintf(fp, "for jj = 1 : nY\n");
            fprintf(fp, "B(ii,jj) = A(nX-ii+1,jj);\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "end;\n");
            fprintf(fp, "Matplot1((B-bmin)/(bmax-bmin)*64,[%e,%e,%e,%e])\n",
                    iLowerB[iplot1],iLowerB[iplot2],
                    iUpperB[iplot1], iUpperB[iplot2]);
            fprintf(fp, "xset(\"colormap\",jetcolormap(64));\n");
            fprintf(fp, "colorbar(bmin,bmax);\n");
            fprintf(fp, "a.thickness = 2;\n");
            fprintf(fp, "a.font_size = 3;\n");
            fprintf(fp, "a.font_style = 4;\n");
            fprintf(fp, "a.box = \"on\";\n");
            fprintf(fp, "a.grid = [1 1];\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
         }
         else
         { 
            fprintf(fp, "if twoPlots == 1\n");
            fprintf(fp, "subplot(1,2,1), surf(X,Y,A)\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            fwritePlotZLabel(fp, outputNames[jplot]);
            sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "subplot(1,2,2)\n");
            fprintf(fp, "end\n");
            fprintf(fp, "contourf(X,Y,A)\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot1]);
            fwritePlotYLabel(fp, inputNames[iplot2]);
            sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
            fwritePlotTitle(fp, winput);
            fprintf(fp, "colorbar\n");
            fprintf(fp, "colormap(jet)\n");
         }
      }
      else if (count == 3)
      {
         fprintf(fp,"xlo = %e; \n", iLowerB[iplot2]);
         fprintf(fp,"xhi = %e; \n", iUpperB[iplot2]);
         fprintf(fp,"ylo = %e; \n", iLowerB[iplot1]);
         fprintf(fp,"yhi = %e; \n", iUpperB[iplot1]);
         fprintf(fp,"zlo = %e; \n", iLowerB[iplot3]);
         fprintf(fp,"zhi = %e; \n", iUpperB[iplot3]);
         fprintf(fp,"X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp,"Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp,"Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp,"V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
            fprintf(fp, "Y(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", faXOut[ind*3]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "X(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", faXOut[ind*3+1]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "Z(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", faXOut[ind*3+2]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
         }
         double GYmax = - PSUADE_UNDEFINED;
         double GYmin =   PSUADE_UNDEFINED;
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
            fprintf(fp, "V(:,:,%d) = [\n", jj + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (ii = 0; ii < nPtsPerDim; ii++)
               {
                  ind = sInd*nPtsPerDim*nPtsPerDim+ii*nPtsPerDim+jj;
                  fprintf(fp, "%e ", tempW[ind]);
                  if (tempW[ind] > GYmax) GYmax = tempW[ind];
                  if (tempW[ind] < GYmin) GYmin = tempW[ind];
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
         }
         fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[iplot2],
                 (iUpperB[iplot2]-iLowerB[iplot2])*0.01, iUpperB[iplot2]);
         fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
                 (iUpperB[iplot1]-iLowerB[iplot1])*0.01, iUpperB[iplot1]);
         fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[iplot3],
                 (iUpperB[iplot3]-iLowerB[iplot3])*0.01, iUpperB[iplot3]);
         fwritePlotCLF(fp);
         fprintf(fp, "isoval = %e;\n", gamma);
         fprintf(fp, "h = patch(isosurface(X,Y,Z,V,isoval),... \n");
         fprintf(fp, "          'FaceColor', 'blue', ... \n");
         fprintf(fp, "          'EdgeColor', 'none', ... \n");
         fprintf(fp, "          'AmbientStrength', 0.2, ... \n");
         fprintf(fp, "          'SpecularStrength', 0.7, ... \n");
         fprintf(fp, "          'DiffuseStrength', 0.4);\n");
         fprintf(fp, "isonormals(X,Y,Z,V,h);\n");
         fprintf(fp, "patch(isocaps(X,Y,Z,V,isoval), ...\n");
         fprintf(fp, "      'FaceColor', 'interp', ... \n");
         fprintf(fp, "      'EdgeColor', 'none'); \n");
         fprintf(fp, "axis([xlo xhi ylo yhi zlo zhi])\n");
         fprintf(fp, "daspect([xhi-xlo, yhi-ylo, zhi-zlo])\n");
         fprintf(fp, "colormap('default'); colorbar\n");
         fprintf(fp, "%%axis tight\n");
         fprintf(fp, "view(3) \n");
         fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
         fprintf(fp, "box on\n");
         fprintf(fp, "grid on\n");
         fprintf(fp, "lighting phong\n");
         fwritePlotAxes(fp);
         fwritePlotXLabel(fp, inputNames[iplot2]);
         fwritePlotYLabel(fp, inputNames[iplot1]);
         fwritePlotZLabel(fp, inputNames[iplot3]);
         sprintf(winput, "Std. Dev. Plot for %s", outputNames[jplot]);
         fwritePlotTitle(fp, winput);
         fprintf(fp,"cin = input('generate slices ? (y or n) ','s');\n");
         fprintf(fp,"if (cin == 'y')\n");
         fprintf(fp,"xin = input('axis to slide through? (x,y,z) ','s');\n");
         fprintf(fp,"N = 101;\n");
         fprintf(fp,"for i = 1 : N\n");
         fprintf(fp,"   display(['displaying ' int2str(i) ' of 101'])\n");
         fprintf(fp,"   if (xin == 'y')\n");
         fprintf(fp,"      h = slice(X,Y,Z,V,xt(i),[],[]);\n");
         fprintf(fp,"   elseif (xin == 'x')\n");
         fprintf(fp,"      h = slice(X,Y,Z,V,[],yt(i),[]);\n");
         fprintf(fp,"   elseif (xin == 'z')\n");
         fprintf(fp,"      h = slice(X,Y,Z,V,[],[],zt(i));\n");
         fprintf(fp,"   end\n");
         fprintf(fp,"   axis([%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e ",
                 iLowerB[iplot2], iUpperB[iplot2], iLowerB[iplot1],
                 iUpperB[iplot1], iLowerB[iplot3], iUpperB[iplot3]);
         fprintf(fp, "%11.4e %11.4e])\n",
                 GYmin-0.1*(GYmax-GYmin),GYmax+0.1*(GYmax-GYmin));
         fwritePlotAxes(fp);
         fwritePlotXLabel(fp, inputNames[iplot2]);
         fwritePlotYLabel(fp, inputNames[iplot1]);
         fwritePlotZLabel(fp, inputNames[iplot3]);
         sprintf(winput, "Std. Dev. Slice Plot for %s", outputNames[jplot]);
         fwritePlotTitle(fp, winput);
         fprintf(fp, "   view(3)\n");
         fprintf(fp, "   colorbar\n");
         fprintf(fp, "   pause(1)\n");
         fprintf(fp, "   if (i < 101)\n");
         fprintf(fp, "      clf\n");
         fprintf(fp, "   end\n");
         fprintf(fp, "end\n");
         fprintf(fp, "end\n");
      }
      else if (count == 4)
      {
         fprintf(fp,"xlo = %e; \n", iLowerB[iplot2]);
         fprintf(fp,"xhi = %e; \n", iUpperB[iplot2]);
         fprintf(fp,"ylo = %e; \n", iLowerB[iplot1]);
         fprintf(fp,"yhi = %e; \n", iUpperB[iplot1]);
         fprintf(fp,"zlo = %e; \n", iLowerB[iplot3]);
         fprintf(fp,"zhi = %e; \n", iUpperB[iplot3]);
         fprintf(fp,"X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp,"Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp,"Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         fprintf(fp,"V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
         for (ii = 0; ii < nPtsPerDim; ii++)
         {
            fprintf(fp, "Y(:,:,%d) = [\n", ii + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*
                        nPtsPerDim;
                  fprintf(fp, "%e ", faXOut[ind*4]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "X(:,:,%d) = [\n", ii + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*
                        nPtsPerDim;
                  fprintf(fp, "%e ", faXOut[ind*4+1]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
            fprintf(fp, "Z(:,:,%d) = [\n", ii + 1);
            for (sInd = 0; sInd < nPtsPerDim; sInd++)
            {
               for (jj = 0; jj < nPtsPerDim; jj++)
               {
                  ind = (sInd*nPtsPerDim*nPtsPerDim+jj*nPtsPerDim+ii)*
                        nPtsPerDim;
                  fprintf(fp, "%e ", faXOut[ind*4+2]);
               }
               fprintf(fp, "\n");
            }
            fprintf(fp, "];\n");
         }
         fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[iplot2],
                 (iUpperB[iplot2]-iLowerB[iplot2])*0.05, iUpperB[iplot2]);
         fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
                 (iUpperB[iplot1]-iLowerB[iplot1])*0.05, iUpperB[iplot1]);
         fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[iplot3],
                 (iUpperB[iplot3]-iLowerB[iplot3])*0.05, iUpperB[iplot3]);
         for (ll = 0; ll < nPtsPerDim; ll++)
         {
            for (ii = 0; ii < nPtsPerDim; ii++)
            {
               fprintf(fp, "V(:,:,%d) = [\n", ii + 1);
               for (sInd = 0; sInd < nPtsPerDim; sInd++)
               {
                  for (jj = 0; jj < nPtsPerDim; jj++)
                  {
                     ind=((sInd*nPtsPerDim+jj)*nPtsPerDim+ii)*nPtsPerDim+ll;
                     fprintf(fp, "%e ", tempW[ind]);
                  }
                  fprintf(fp, "\n");
               }
               fprintf(fp, "];\n");
            }
            fprintf(fp, "disp('Frame %d of %d')\n", ll+1, nPtsPerDim);
            fwritePlotCLF(fp);
            fprintf(fp, "isoval = %e;\n", gamma);
            fprintf(fp, "h = patch(isosurface(X,Y,Z,V,isoval),... \n");
            fprintf(fp, "          'FaceColor', 'blue', ... \n");
            fprintf(fp, "          'EdgeColor', 'none', ... \n");
            fprintf(fp, "          'AmbientStrength', 0.2, ... \n");
            fprintf(fp, "          'SpecularStrength', 0.7, ... \n");
            fprintf(fp, "          'DiffuseStrength', 0.4);\n");
            fprintf(fp, "isonormals(X,Y,Z,V,h);\n");
            fprintf(fp, "patch(isocaps(X,Y,Z,V,isoval), ...\n");
            fprintf(fp, "      'FaceColor', 'interp', ... \n");
            fprintf(fp, "      'EdgeColor', 'none'); \n");
            fprintf(fp, "axis([xlo xhi ylo yhi zlo zhi])\n");
            fprintf(fp, "daspect([xhi-xlo, yhi-ylo, zhi-zlo])\n");
            fprintf(fp, "colormap('default'); colorbar\n");
            fprintf(fp, "%%axis tight\n");
            fprintf(fp, "view(3) \n");
            fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
            fprintf(fp, "lighting phong\n");
            fwritePlotAxes(fp);
            fwritePlotXLabel(fp, inputNames[iplot2]);
            fwritePlotYLabel(fp, inputNames[iplot1]);
            fwritePlotZLabel(fp, inputNames[iplot3]);
            fprintf(fp, "title('3D Std Dev Isosurface Plot at %s=%e',",
                    inputNames[iplot4],faXOut[ll*4+3]);
            fprintf(fp, "'FontWeight','bold','FontSize',12)\n");
            fprintf(fp, "pause(1)\n");
         }
      }
      fclose(fp);
      if (psPlotTool_ == 1)
           printf("\nscilabrssd.sci is now available.\n");
      else printf("\nmatlabrssd.m is now available.\n");
      delete [] faXOut;
      delete [] faYOut;
      delete [] faYIn;
      delete [] tempW;
      delete [] tempX;
      delete faPtr;
      delete [] inputSettings;
   }

   // +++ rsi2
   else if (!strcmp(command, "rsi2"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rsi2: generate intersection surfaces for >1 outputs.\n");
       printf("syntax: rsi2 (no argument needed).\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data to analyze (load sample first).\n");
       return -1;
     }
     if (nInputs < 2)
     {
       printf("ERROR: rsi2 requires 2 or more inputs.\n");
       return -1;
     }
     if (nOutputs < 2)
     {
       printf("ERROR: rsi2 requires 2 or more outputs.\n");
       return -1;
     }
     if (psPlotTool_ == 1)
     {
       printf("INFO: rsi2 is currently not available for scilab.\n");
       return -1;
     }
     int nPtsPerDim = 32;
     sprintf(pString, "Grid resolution ? (32 - 256) ");
     nPtsPerDim = getInt(32, 256, pString);
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);

     int iplot1, iplot2, iInd1;
     double *inputSettings = new double[nInputs];
     iplot1 = iplot2 = -1;
     sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
     iplot1 = getInt(1, nInputs, pString);
     iplot1--;
     iplot2 = iplot1;
     while (iplot1 == iplot2)
     {
       sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
              nInputs, iplot1+1);
       iplot2 = getInt(1, nInputs, pString);
       iplot2--;
       if (iplot1 == iplot2)
         printf("ERROR: duplicate input number %d.\n", iplot2+1);
     }
     sprintf(pString,"Set other nominal values automatically ? (y or n) ");
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
         else inputSettings[iInd1] = 1.0;
       }
     }
     else
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2)
         {
           inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
           while (inputSettings[iInd1] < iLowerB[iInd1] ||
                  inputSettings[iInd1] > iUpperB[iInd1])
           {
             sprintf(pString,
                     "Enter nominal value for input %d (%e - %e):",
                     iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
             inputSettings[iInd1] = getDouble(pString);
           }
         }
         else inputSettings[iInd1] = 1.0;
       }
     }
     fp = fopen("matlabrsi2.m", "w");
     if (fp == NULL)
     {
       printf("ERROR: cannot open file matlabrsi2.m.\n");
       delete [] inputSettings;
       delete faPtr;
       return -1;
     }
     fprintf(fp, "twoPlots = 1;\n");
     fprintf(fp, "fs = 10;\n");
     fwritePlotCLF(fp);
     
     int rsiNOutputs = 2;
     sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
     rsiNOutputs = getInt(2, nOutputs, pString);

     int *rsiSet = new int[rsiNOutputs];
     if (rsiNOutputs == nOutputs)
     {
       for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
     }
     else
     {
       for (ii = 0; ii < rsiNOutputs; ii++)
       {
         sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                 ii+1, nOutputs);
         rsiSet[ii] = getInt(1, nOutputs, pString);
         rsiSet[ii]--;
       }
     }
     if (rsiNOutputs > 5)
     {
       printf("INFO: rsi2 only shows the constrained response surfaces\n");
       printf("      for the first 5 outputs and then the aggregate.\n");
     }
      
     double *faYIn = new double[nSamples];
     int    **rsiMatrix = new int*[nPtsPerDim];
     for (ii = 0; ii < nPtsPerDim; ii++)
     {
       rsiMatrix[ii] = new int[nPtsPerDim];
       for (jj = 0; jj < nPtsPerDim; jj++)
         rsiMatrix[ii][jj] = rsiNOutputs;
     }

     int    jplot, ind, ind2, sInd, faLeng=0, count;
     double Ymin, Ymax, threshU, threshL, *faXOut=NULL, *faYOut=NULL;
     for (ii = 0; ii < rsiNOutputs; ii++)
     {
        jplot = rsiSet[ii];
        for (sInd = 0; sInd < nSamples; sInd++)
          faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

        faPtr->gen2DGridData(sampleInputs,faYIn, iplot1, iplot2, 
                             inputSettings, &faLeng, &faXOut,&faYOut);

        Ymin = faYOut[0];
        for (sInd = 1; sInd < faLeng; sInd++)
          if (faYOut[sInd] < Ymin) Ymin = faYOut[sInd];
        Ymax = faYOut[0];
        for (sInd = 1; sInd < faLeng; sInd++)
          if (faYOut[sInd] > Ymax) Ymax = faYOut[sInd];

        printf("Ymin and Ymax = %e %e\n", Ymin, Ymax);
        sprintf(pString,
                "Enter the lower threshold for output %d (min = %16.8e) : ",
                jplot+1, Ymin);
        threshL = getDouble(pString);
        sprintf(pString,
                "Enter the upper threshold for output %d (max = %16.8e) : ",
                jplot+1, Ymax);
        threshU = getDouble(pString);

        if (ii == 0)
        {
          fprintf(fp, "x = [\n");
          for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
            fprintf(fp, "%e\n", faXOut[sInd*2]);
          fprintf(fp, "];\n");
          fprintf(fp, "y = [\n");
          for (sInd = 0; sInd < nPtsPerDim; sInd++)
            fprintf(fp, "%e\n", faXOut[sInd*2+1]);
          fprintf(fp, "];\n");
        }
        if (ii < 5)
        {
          fprintf(fp, "A%d = [\n", ii+1);
          for (sInd = 0; sInd < faLeng; sInd++)
            fprintf(fp, "%e\n", faYOut[sInd]);
          fprintf(fp, "];\n");
          fprintf(fp, "A%d = reshape(A%d,%d,%d);\n",ii+1,ii+1,
                  nPtsPerDim,nPtsPerDim);
          fprintf(fp, "yLo = %e;\n", threshL);
          fprintf(fp, "yHi = %e;\n", threshU);
          fprintf(fp, "nA  = size(A%d,1);\n", ii+1);
          fprintf(fp, "[ia,ja,aa] = find(A%d<yLo);\n", ii+1);
          fprintf(fp, "for ii = 1 : length(ia)\n");
          fprintf(fp, "   A%d(ia(ii),ja(ii)) = NaN;\n", ii+1); 
          fprintf(fp, "end;\n");
          fprintf(fp, "n1 = length(ia);\n");
          fprintf(fp, "[ia,ja,aa] = find(A%d>yHi);\n", ii+1);
          fprintf(fp, "for ii = 1 : length(ia)\n");
          fprintf(fp, "   A%d(ia(ii),ja(ii)) = NaN;\n", ii+1); 
          fprintf(fp, "end;\n");
          fprintf(fp, "n2 = length(ia);\n");
          fprintf(fp, "if (n1 + n2 == nA*nA)\n");
          fprintf(fp, "   A%d(1,1) = 0;\n",ii+1);
          fprintf(fp, "   A%d(%d,%d) = 1;\n",ii+1,nPtsPerDim,
                  nPtsPerDim);
          fprintf(fp, "end;\n");
          if (ii == 0) fprintf(fp, "subplot(2,3,1)\n");
          if (ii == 1) fprintf(fp, "subplot(2,3,2)\n");
          if (ii == 2) fprintf(fp, "subplot(2,3,3)\n");
          if (ii == 3) fprintf(fp, "subplot(2,3,4)\n");
          if (ii == 4) fprintf(fp, "subplot(2,3,5)\n");
          fprintf(fp, "contourf(x,y,A%d)\n", ii+1);
          fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
                  iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
          fwritePlotAxes(fp);
          fprintf(fp, "xlabel('%s','FontSize',fs,'FontWeight','bold')\n",
                  inputNames[iplot1]);
          fprintf(fp, "ylabel('%s','Fontsize',fs,'FontWeight','bold')\n",
                  inputNames[iplot2]);
          fprintf(fp, "title('%s',",outputNames[jplot]);
          fprintf(fp, "'FontWeight','bold','FontSize',fs)\n");
          fprintf(fp, "colorbar\n");
        }

        for (sInd = 0; sInd < faLeng; sInd++)
        {
           ind  = sInd % nPtsPerDim;
           ind2 = sInd / nPtsPerDim;
           if (faYOut[sInd] < threshL) rsiMatrix[ind][ind2]--;
           if (faYOut[sInd] > threshU) rsiMatrix[ind][ind2]--;
        }
        delete [] faXOut;
        delete [] faYOut;
     }


     fprintf(fp, "A = [\n");
     count = 0;
     for (ii = 0;  ii < nPtsPerDim; ii++)
       for (jj = 0;  jj < nPtsPerDim; jj++)
         if (rsiMatrix[jj][ii] == 0) count++;
     if (count == nPtsPerDim*nPtsPerDim)
     {
       for (ii = 0;  ii < nPtsPerDim; ii++)
         for (jj = 0;  jj < nPtsPerDim; jj++) fprintf(fp, "0\n");
     }
     else
     {
       for (ii = 0;  ii < nPtsPerDim; ii++)
       {
         for (jj = 0;  jj < nPtsPerDim; jj++)
         {
           if (rsiMatrix[jj][ii] == 0) fprintf(fp, "NaN\n");
           else fprintf(fp, "%d\n", rsiMatrix[jj][ii]);
         }
       }
     }
     fprintf(fp, "];\n");
     fprintf(fp, "A = reshape(A,%d,%d);\n",nPtsPerDim, nPtsPerDim);
     fprintf(fp, "A(%d,%d) = %e;\n", nPtsPerDim, nPtsPerDim, 
             (double) rsiNOutputs);
     //fprintf(fp, "if twoPlots == 1\n");
     //fprintf(fp, "subplot(2,3,5), mesh(y,x,A)\n");
     //fprintf(fp, "axis([%e %e %e %e])\n",iLowerB[iplot1],
     //        iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
     //fwritePlotAxes(fp);
     //fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
     //        inputNames[iplot1]);
     //fprintf(fp, "ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
     //        inputNames[iplot2]);
     //fprintf(fp, "title('Intersection Plot','FontWeight',");
     //fprintf(fp, "'bold','FontSize',12)\n");
     //fprintf(fp, "colorbar\n");
     //fprintf(fp, "colormap(cool)\n");
     //fprintf(fp, "end\n");
     fprintf(fp,"subplot(2,3,6), contourf(x,y,A)\n");
     fprintf(fp,"axis([%e %e %e %e])\n",iLowerB[iplot1],
             iUpperB[iplot1],iLowerB[iplot2],iUpperB[iplot2]);
     fwritePlotAxes(fp);
     fprintf(fp,"xlabel('%s','FontSize',fs,'FontWeight','bold')\n",
             inputNames[iplot1]);
     fprintf(fp,"ylabel('%s','Fontsize',fs,'FontWeight','bold')\n",
             inputNames[iplot2]);
     fprintf(fp,"title('Intersection (color=deg of overlap)','FontWeight',");
     fprintf(fp,"'bold','FontSize',fs)\n");
     fprintf(fp,"colorbar\n");
     fprintf(fp,"colormap(cool)\n");
     fprintf(fp,"disp('On intersection plot, if a region has a color value");
     fprintf(fp," of 2, it means it is feasible for 2 outputs.')\n");
     fclose(fp);
     printf("matlabrsi2.m is now available for plotting.\n");

     delete [] inputSettings;
     delete [] rsiSet;
     delete [] faYIn;
     delete faPtr;
     for (ii = 0; ii < nPtsPerDim; ii++) delete [] rsiMatrix[ii];
      delete [] rsiMatrix;
   }

   // +++ rsi3 
   else if (!strcmp(command, "rsi3"))
   {
     sscanf(lineIn,"%s %s",command,winput);
     if (!strcmp(winput, "-h"))
     {
       printf("rsi3: generate intersection surfaces for >1 outputs\n");
       printf("syntax: rsi3 (no argument needed).\n");
       return 0;
     }
     if (psuadeIO == NULL || sampleOutputs == NULL)
     {
       printf("ERROR: no data to analyze (load sample first).\n");
       return -1;
     }
     if (nInputs < 3)
     {
       printf("ERROR: rsi3 requires 3 or more inputs.\n");
       return -1;
     }
     if (nOutputs < 2)
     {
       printf("ERROR: rsi3 requires 2 or more outputs.\n");
       return -1;
     }
     if (psPlotTool_ == 1)
     {
       printf("INFO: rsi3 is currently not available for scilab.\n");
       return -1;
     }
     int nPtsPerDim = 24;
     sprintf(pString, "Grid resolution ? (16 - 32) ");
     nPtsPerDim = getInt(16, 32, pString);
     int faFlag = 1;
     FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
     if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
     faPtr->setNPtsPerDim(nPtsPerDim);
     faPtr->setBounds(iLowerB, iUpperB);
     faPtr->setOutputLevel(outputLevel);

     int    iplot1, iplot2, iplot3, iInd1, sInd;
     double *inputSettings = new double[nInputs];
     iplot1 = iplot2 = iplot3 = -1;
     sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
     iplot1 = getInt(1, nInputs, pString);
     iplot1--;
     iplot2 = iplot1;
     while (iplot1 == iplot2)
     {
       sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
               nInputs, iplot1+1);
       iplot2 = getInt(1, nInputs, pString);
       iplot2--;
       if (iplot1 == iplot2)
         printf("ERROR: duplicate input number %d.\n",iplot2+1);
     }
     if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
     while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
     {
       sprintf(pString,
               "Enter the input for z axis (1 - %d), not %d nor %d: ",
               nInputs, iplot1+1, iplot2+1);
       iplot3 = getInt(1, nInputs, pString);
       iplot3--;
       if (iplot3 == iplot1 || iplot3 == iplot2)
         printf("ERROR: duplicate input number %d.\n",iplot3+1);
     }
     sprintf(pString,"Set other nominal values automatically ? (y or n) ");
     getString(pString, winput);
     if (winput[0] == 'y')
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
              inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
         else inputSettings[iInd1] = 1.0;
       }
     }
     else
     {
       for (iInd1 = 0; iInd1 < nInputs; iInd1++)
       {
         if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
         {
           inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
           sprintf(pString,
                   "Enter nominal value for input %d (%e - %e): ", 
                   iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
           while (inputSettings[iInd1] < iLowerB[iInd1] ||
                  inputSettings[iInd1] > iUpperB[iInd1])
             inputSettings[iInd1] = getDouble(pString);
         }
         else inputSettings[iInd1] = 1.0;
       }
     }
     fp = fopen("matlabrsi3.m", "w");
     if (fp == NULL)
     {
       printf("ERROR: cannot open file matlabrsi3.m.\n");
       delete [] inputSettings;
       delete faPtr;
     }
     fwritePlotCLF(fp);

     int rsiNOutputs = 1;
     sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
     rsiNOutputs = getInt(2, nOutputs, pString);

     double *threshLs = new double[rsiNOutputs];
     double *threshUs = new double[rsiNOutputs];

     int *rsiSet = new int[rsiNOutputs];
     if (rsiNOutputs == nOutputs)
     {
       for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
     }
     else
     {
       for (ii = 0; ii < rsiNOutputs; ii++)
       {
         sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                 ii+1, nOutputs);
         rsiSet[ii] = getInt(1, nOutputs, pString);
         rsiSet[ii]--;
       }
     }
     double *faYIn = new double[nSamples];

     printf("Please wait while generating the RS data \n");
     fprintf(fp, "xlo = %e; \n", iLowerB[iplot2]);
     fprintf(fp, "xhi = %e; \n", iUpperB[iplot2]);
     fprintf(fp, "ylo = %e; \n", iLowerB[iplot1]);
     fprintf(fp, "yhi = %e; \n", iUpperB[iplot1]);
     fprintf(fp, "zlo = %e; \n", iLowerB[iplot3]);
     fprintf(fp, "zhi = %e; \n", iUpperB[iplot3]);
     fprintf(fp, "X=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp, "Y=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp, "Z=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);
     fprintf(fp, "V=zeros(%d,%d,%d);\n",nPtsPerDim,nPtsPerDim,nPtsPerDim);

     int    faLeng, jplot, ind, count;
     double GYmax, GYmin, gamma, *faXOut=NULL, *faYOut=NULL;
     for (ii = 0; ii < rsiNOutputs; ii++)
     {
       jplot = rsiSet[ii];
       for (sInd = 0; sInd < nSamples; sInd++)
          faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

       faLeng = 0;
       faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                            inputSettings, &faLeng, &faXOut,&faYOut);
       GYmin = faYOut[0];
       for (sInd = 1; sInd < faLeng; sInd++)
         if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
       GYmax = faYOut[0];
       for (sInd = 1; sInd < faLeng; sInd++)
         if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];

       printf("\nOutput %d : Ymin and Ymax found = %e %e.\n", jplot+1,
              GYmin, GYmax);
       sprintf(pString,"Enter the lower threshold (min = %e) : ", GYmin);
       threshLs[ii] = getDouble(pString);
       sprintf(pString,"Enter the upper threshold (max = %e) : ", GYmax);
       threshUs[ii] = getDouble(pString);
       if (ii == 0) gamma = threshLs[ii];  
       else         gamma = (gamma < threshLs[ii]) ? gamma : threshLs[ii];

       if (ii == 0)
       {
         for (jj = 0; jj < nPtsPerDim; jj++)
         {
           fprintf(fp, "Y(:,:,%d) = [\n", jj + 1);
           for (sInd = 0; sInd < nPtsPerDim; sInd++)
           {
             for (kk = 0; kk < nPtsPerDim; kk++)
             {
               ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
               fprintf(fp, "%e ", faXOut[ind*3]);
             }
             fprintf(fp, "\n");
           }
           fprintf(fp, "];\n");
           fprintf(fp, "X(:,:,%d) = [\n", jj + 1);
           for (sInd = 0; sInd < nPtsPerDim; sInd++)
           {
             for (kk = 0; kk < nPtsPerDim; kk++)
             {
               ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
               fprintf(fp, "%e ", faXOut[ind*3+1]);
             }
             fprintf(fp, "\n");
           }
           fprintf(fp, "];\n");
           fprintf(fp, "Z(:,:,%d) = [\n", jj + 1);
           for (sInd = 0; sInd < nPtsPerDim; sInd++)
           {
             for (kk = 0; kk < nPtsPerDim; kk++)
             {
               ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
               fprintf(fp, "%e ", faXOut[ind*3+2]);
             }
             fprintf(fp, "\n");
           }
           fprintf(fp, "];\n");
         }
       }
       delete [] faXOut;
       delete [] faYOut;
     }

     for (ii = 0; ii < rsiNOutputs; ii++)
     {
       jplot = rsiSet[ii];
       for (sInd = 0; sInd < nSamples; sInd++)
         faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

       faLeng = 0;
       faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                            inputSettings, &faLeng, &faXOut,&faYOut);
       for (jj = 0; jj < nPtsPerDim; jj++)
       {
         fprintf(fp, "V%d(:,:,%d) = [\n", ii+1, jj+1);
         for (sInd = 0; sInd < nPtsPerDim; sInd++)
         {
           for (kk = 0; kk < nPtsPerDim; kk++)
           {
             ind = sInd*nPtsPerDim*nPtsPerDim+kk*nPtsPerDim+jj;
             if (faYOut[ind] < threshLs[ii])
             {
               fprintf(fp, "%e ", gamma);
               count++;
             }
             else if (faYOut[ind] > threshUs[ii])
             {
               fprintf(fp, "%e ", gamma);
               count++;
             }
             else fprintf(fp, "%e ", faYOut[ind]);
           }
           fprintf(fp, "\n");
         }
         fprintf(fp, "];\n");
       }
       delete [] faXOut;
       delete [] faYOut;
       if (ii == 0) fprintf(fp, "V = V%d;\n", ii+1);
       else         fprintf(fp, "V = min(V, V%d);\n", ii+1);
     }

     double threshL, threshU;
     threshL = threshLs[0];
     for (ii = 1; ii < rsiNOutputs; ii++)
       if (threshLs[ii] < threshL) threshL = threshLs[ii];
     threshU = threshUs[0];
     for (ii = 1; ii < rsiNOutputs; ii++)
       if (threshUs[ii] > threshU) threshU = threshUs[ii];
     fprintf(fp, "xt = [%e:%e:%e];\n", iLowerB[iplot2],
             (iUpperB[iplot2]-iLowerB[iplot2])*0.01, iUpperB[iplot2]);
     fprintf(fp, "yt = [%e:%e:%e];\n", iLowerB[iplot1],
             (iUpperB[iplot1]-iLowerB[iplot1])*0.01, iUpperB[iplot1]);
     fprintf(fp, "zt = [%e:%e:%e];\n", iLowerB[iplot3],
             (iUpperB[iplot3]-iLowerB[iplot3])*0.01, iUpperB[iplot3]);
     fprintf(fp, "isoval = %e;\n", gamma);
     fprintf(fp, "h = patch(isosurface(X,Y,Z,V,isoval),... \n");
     fprintf(fp, "          'FaceColor', 'blue', ... \n");
     fprintf(fp, "          'EdgeColor', 'none', ... \n");
     fprintf(fp, "          'AmbientStrength', 0.2, ... \n");
     fprintf(fp, "          'SpecularStrength', 0.7, ... \n");
     fprintf(fp, "          'DiffuseStrength', 0.4);\n");
     fprintf(fp, "isonormals(X,Y,Z,V,h);\n");
     fprintf(fp, "patch(isocaps(X,Y,Z,V,isoval), ...\n");
     fprintf(fp, "      'FaceColor', 'interp', ... \n");
     fprintf(fp, "      'EdgeColor', 'none'); \n");
     fprintf(fp, "axis([xlo xhi ylo yhi zlo zhi])\n");
     fprintf(fp, "daspect([%e,%e,%e])\n",iUpperB[iplot2]-iLowerB[iplot2],
             iUpperB[iplot1]-iLowerB[iplot1],
             iUpperB[iplot3]-iLowerB[iplot3]);
     fprintf(fp, "   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
             inputNames[iplot2]);
     fprintf(fp, "   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot1]);
     fprintf(fp, "   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot3]);
     fwritePlotAxes(fp);
     fprintf(fp, "%%colormap('default'); colorbar\n");
     fprintf(fp, "%%axis tight\n");
     fprintf(fp, "view(3) \n");
     fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
     fprintf(fp, "lighting phong\n");
     fprintf(fp, "cin = input('generate slices ? (y or n) ','s');\n");
     fprintf(fp, "if (cin == 'y')\n");
     fprintf(fp, "xin = input('axis to slide through ? (x,y,z) ','s');\n");
     fprintf(fp, "for i = 1 : 101\n");
     fprintf(fp, "   if (xin == 'y')\n");
     fprintf(fp, "      h = contourslice(X,Y,Z,V,xt(i),[],[],101);\n");
     fprintf(fp, "   elseif (xin == 'x')\n");
     fprintf(fp, "      h = contourslice(X,Y,Z,V,[],yt(i),[],101);\n");
     fprintf(fp, "   elseif (xin == 'z')\n");
     fprintf(fp, "      h = contourslice(X,Y,Z,V,[],[],zt(i),101);\n");
     fprintf(fp, "   end\n");
     fprintf(fp, "   axis([%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e ",
             iLowerB[iplot2], iUpperB[iplot2], iLowerB[iplot1],
             iUpperB[iplot1], iLowerB[iplot3], iUpperB[iplot3]);
     fprintf(fp, "%11.4e %11.4e])\n",
             threshL-0.2*(threshU-threshL),threshU+0.2*(threshU-threshL));
     fwritePlotAxes(fp);
     fprintf(fp, "   xlabel('%s','FontSize',12,'FontWeight','bold')\n",
             inputNames[iplot2]);
     fprintf(fp, "   ylabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot1]);
     fprintf(fp, "   zlabel('%s','Fontsize',12,'FontWeight','bold')\n",
             inputNames[iplot3]);
     fprintf(fp, "colormap('default'); colorbar\n");
     fprintf(fp, "view(3) \n");
     fprintf(fp, "set(gcf,'Renderer','zbuffer')\n");
     fprintf(fp, "lighting phong\n");
     fprintf(fp, "pause(1)\n");
     fprintf(fp," if (i < 101)\n");
     fprintf(fp,"   clf\n");
     fprintf(fp," end\n");
     fprintf(fp, "end\n");
     fprintf(fp, "end\n");
     fclose(fp);
     printf("matlabrsi3.m is now available for response surface and ");
     printf("contour plots\n");

     delete [] rsiSet;
     delete [] faYIn;
     delete faPtr;
     delete [] inputSettings;
     delete [] threshLs;
     delete [] threshUs;
   }

   // +++ rsi3m
   else if (!strcmp(command, "rsi3m"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rsi3m: generate intersection surfaces for >1 outputs\n");
         printf("syntax: rsi3m (no argument needed).\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data to analyze (load sample first).\n");
         return -1;
      }
      if (nInputs < 3)
      {
         printf("ERROR: rsi3m requires 3 or more inputs.\n");
         return -1;
      }
      if (nOutputs < 2)
      {
         printf("ERROR: rsi3m requires 2 or more outputs.\n");
         return -1;
      }
      if (psPlotTool_ == 1)
      {
         printf("INFO: rsi3m is currently not available for scilab.\n");
         return -1;
      }
      int nPtsPerDim = 24;
      sprintf(pString, "Grid resolution ? (16 - 32) ");
      nPtsPerDim = getInt(16, 32, pString);
      int faFlag = 1;
      FuncApprox *faPtr = genFAInteractive(psuadeIO, faFlag);
      if (faPtr == NULL) {printf("ERROR detected.\n"); return -1;}
      faPtr->setNPtsPerDim(nPtsPerDim);
      faPtr->setBounds(iLowerB, iUpperB);
      faPtr->setOutputLevel(outputLevel);

      int iplot1, iplot2, iplot3, iInd1, sInd;
      double *inputSettings = new double[nInputs];
      iplot1 = iplot2 = iplot3 = -1;
      sprintf(pString, "Enter the input for x axis (1 - %d) : ", nInputs);
      iplot1 = getInt(1, nInputs, pString);
      iplot1--;
      iplot2 = iplot1;
      while (iplot1 == iplot2)
      {
         sprintf(pString,"Enter the input for y axis (1 - %d), not %d : ",
                 nInputs, iplot1+1);
         iplot2 = getInt(1, nInputs, pString);
         iplot2--;
         if (iplot1 == iplot2)
            printf("ERROR: duplicate input number %d.\n",iplot2+1);
      }
      if (nInputs == 3) iplot3 = 3 - iplot1 - iplot2;
      while (iplot3 < 0 || iplot3 == iplot1 || iplot3 == iplot2)
      {
         sprintf(pString,
                 "Enter the input for t axis (1 - %d), not %d nor %d: ",
                 nInputs, iplot1+1, iplot2+1);
         iplot3 = getInt(1, nInputs, pString);
         iplot3--;
         if (iplot3 == iplot1 || iplot3 == iplot2)
            printf("ERROR: duplicate input number %d.\n",iplot3+1);
      }
      sprintf(pString,"Set other nominal values automatically ? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y')
      {
         for (iInd1 = 0; iInd1 < nInputs; iInd1++)
         {
            if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
                 inputSettings[iInd1] = 0.5*(iLowerB[iInd1]+iUpperB[iInd1]);
            else inputSettings[iInd1] = 1.0;
         }
      }
      else
      {
         for (iInd1 = 0; iInd1 < nInputs; iInd1++)
         {
            if (iInd1 != iplot1 && iInd1 != iplot2 && iInd1 != iplot3)
            {
               inputSettings[iInd1] = iLowerB[iInd1] - 1.0;
               sprintf(pString,
                       "Enter nominal value for input %d (%e - %e): ", 
                       iInd1+1, iLowerB[iInd1], iUpperB[iInd1]);
               while (inputSettings[iInd1] < iLowerB[iInd1] ||
                      inputSettings[iInd1] > iUpperB[iInd1])
                  inputSettings[iInd1] = getDouble(pString);
            }
            else inputSettings[iInd1] = 1.0;
         }
      }
      fp = fopen("matlabrsi3m.m", "w");
      if (fp == NULL)
      {
         printf("ERROR: cannot open file matlabrsi3m.m.\n");
         delete [] inputSettings;
         delete faPtr;
         return -1;
      }
      fprintf(fp, "hold off\n");
      fwritePlotCLF(fp);
      fprintf(fp, "disp(\'Please wait while loading.\')\n");
      fprintf(fp, "pause(1)\n");
      
      int rsiNOutputs = 2;
      sprintf(pString,"How many outputs to use ? (2 - %d) ",nOutputs);
      rsiNOutputs = getInt(2, nOutputs, pString);

      int *rsiSet = new int[rsiNOutputs];
      if (rsiNOutputs == nOutputs)
      {
         for (ii = 0; ii < rsiNOutputs; ii++) rsiSet[ii] = ii;
      }
      else
      {
         for (ii = 0; ii < rsiNOutputs; ii++)
         {
            sprintf(pString,"Enter the %d-th output index (1 - %d) : ",
                    ii+1, nOutputs);
            rsiSet[ii] = getInt(1, nOutputs, pString);
            rsiSet[ii]--;
         }
      }
      double *faYIn = new double[nSamples];
      int    jplot, faLeng;
      double *faXOut=NULL, *faYOut=NULL, GYmax, GYmin, threshL, threshU;

      printf("Please wait while generating the RS data \n");
      for (jj = 0; jj < nPtsPerDim; jj++)
         fprintf(fp, "M%d = %e * ones(%d);\n", jj+1, 1.0*rsiNOutputs, 
                 nPtsPerDim);
      for (ii = 0; ii < rsiNOutputs; ii++)
      {
         jplot = rsiSet[ii];
         for (sInd = 0; sInd < nSamples; sInd++)
            faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];

         faLeng = 0;
         faPtr->gen3DGridData(sampleInputs,faYIn, iplot1, iplot2, iplot3, 
                                 inputSettings, &faLeng, &faXOut,&faYOut);
         GYmin = faYOut[0];
         for (sInd = 1; sInd < faLeng; sInd++)
            if (faYOut[sInd] < GYmin) GYmin = faYOut[sInd];
         GYmax = faYOut[0];
         for (sInd = 1; sInd < faLeng; sInd++)
            if (faYOut[sInd] > GYmax) GYmax = faYOut[sInd];

         for (jj = 0; jj < nPtsPerDim; jj++)
         {
            printf(".");
            fflush(stdout);

            fprintf(fp, "A%d_%d = [\n", ii+1, jj+1);
            for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim)
               fprintf(fp, "%e\n", faYOut[sInd+jj]);
            fprintf(fp, "];\n");
            fprintf(fp, "A%d_%d = reshape(A%d_%d,%d,%d);\n", ii+1, jj+1,
                    ii+1, jj+1, nPtsPerDim, nPtsPerDim);

            if (ii == 0 && jj == 0)
            {
               fprintf(fp, "x = [\n");
               for (sInd = 0; sInd < faLeng; sInd+=nPtsPerDim*nPtsPerDim)
                  fprintf(fp, "%e\n", faXOut[sInd*3]);
               fprintf(fp, "];\n");
               fprintf(fp, "y = [\n");
               for (sInd = 0; sInd < nPtsPerDim*nPtsPerDim; sInd+=nPtsPerDim)
                  fprintf(fp, "%e\n", faXOut[sInd*3+1]);
               fprintf(fp, "];\n");
            }
         }
         delete [] faXOut;
         delete [] faYOut;
         printf("\nOutput %d : Ymin and Ymax found = %e %e.\n", jplot+1,
                GYmin, GYmax);
         sprintf(pString,"Enter the lower threshold (min = %e) : ", GYmin);
         threshL = getDouble(pString);
         sprintf(pString,"Enter the upper threshold (max = %e) : ", GYmax);
         threshU = getDouble(pString);

         for (jj = 0; jj < nPtsPerDim; jj++)
         {
            fprintf(fp, "B%d_%d = A%d_%d;\n",ii+1,jj+1,ii+1,jj+1);
            fprintf(fp, "nA  = size(A%d_%d,1);\n", ii+1, jj+1);
            fprintf(fp, "n1 = 0;\n");
            fprintf(fp, "n2 = 0;\n");
            if (threshL > GYmin)
            { 
               fprintf(fp, "yLo = %e;\n", threshL);
               fprintf(fp, "[ia,ja,aa] = find(A%d_%d<yLo);\n",ii+1,jj+1);
               fprintf(fp, "for ii = 1 : length(ia)\n");
               fprintf(fp, "   B%d_%d(ia(ii),ja(ii))=NaN;\n",ii+1,jj+1);
               fprintf(fp, "   M%d(ia(ii),ja(ii))=M%d(ia(ii),ja(ii))-1;\n", 
                       jj+1,jj+1);
               fprintf(fp, "end;\n");
               fprintf(fp, "n1 = length(ia);\n");
            }
            if (threshU < GYmax)
            { 
               fprintf(fp, "yHi = %e;\n", threshU);
               fprintf(fp, "[ia,ja,aa] = find(A%d_%d>yHi);\n",ii+1,jj+1);
               fprintf(fp, "for ii = 1 : length(ia)\n");
               fprintf(fp, "   B%d_%d(ia(ii),ja(ii))=NaN;\n",ii+1,jj+1);
               fprintf(fp, "   M%d(ia(ii),ja(ii))=M%d(ia(ii),ja(ii))-1;\n", 
                       jj+1,jj+1);
               fprintf(fp, "end;\n");
               fprintf(fp, "n1 = length(ia);\n");
            }
            fprintf(fp, "if (n1+n2 == nA*nA)\n");
            fprintf(fp, "   B%d_%d(1,1)=0;\n",ii+1,jj+1);
            fprintf(fp, "   B%d_%d(%d,%d)=1;\n",ii+1,jj+1,
                    nPtsPerDim,nPtsPerDim);
            fprintf(fp, "end;\n");
         }
      }
      for (jj  = 0; jj < nPtsPerDim; jj++)
      {
         fprintf(fp, "[ia,ja,aa] = find(M%d == 0);\n", jj+1);
         fprintf(fp, "nM  = size(M%d,1);\n", jj+1);
         fprintf(fp, "for ii = 1 : length(ia)\n");
         fprintf(fp, "   M%d(ia(ii),ja(ii)) = NaN;\n", jj+1);
         fprintf(fp, "end;\n");
         fprintf(fp, "if (length(ia) == nM*nM)\n");
         fprintf(fp, "   M%d(1,1) = 0;\n", jj+1);
         fprintf(fp, "   M%d(nM,nM) = %e;\n", jj+1, 1.0*rsiNOutputs);
         fprintf(fp, "end;\n");
         fprintf(fp, "Mmax = max(max(M%d));\n", jj+1);
         fprintf(fp, "if (Mmax ~= %d)\n", rsiNOutputs);
         fprintf(fp, "   M%d(%d,%d) = %d;\n", jj+1, nPtsPerDim,
                 nPtsPerDim, rsiNOutputs);
         fprintf(fp, "end;\n");
         fprintf(fp, "Mmin = min(min(M%d));\n", jj+1);
         fprintf(fp, "if (Mmin ~= 0)\n");
         fprintf(fp, "   M%d(1,1) = 0;\n", jj+1);
         fprintf(fp, "end;\n");
      }

      for (jj = 0; jj < nPtsPerDim; jj++)
      {
         inputSettings[iplot3] = (iUpperB[iplot3] - iLowerB[iplot3]) *
                                 jj / (nPtsPerDim - 1.0) + iLowerB[iplot3];
         fprintf(fp, "disp(\'Plotting frame %d of %d\')\n",jj+1,nPtsPerDim);
         fprintf(fp, "subplot(2,3,1), contourf(x,y,B1_%d)\n", jj+1);
         fwritePlotAxes(fp);
         fprintf(fp,"title(\'Contour Plot for %s\',",outputNames[rsiSet[0]]);
         fprintf(fp, "'FontSize',12,'FontWeight','bold')\n"); 
         fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                 inputNames[iplot2]);
         fprintf(fp, "subplot(2,3,2), contourf(x,y,B2_%d)\n", jj+1);
         fwritePlotAxes(fp);
         fprintf(fp,"title(\'Contour Plot for %s\',",outputNames[rsiSet[1]]);
         fprintf(fp, "'FontSize',12,'FontWeight','bold')\n"); 
         fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                 inputNames[iplot1]);
         fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                 inputNames[iplot2]);
         if (rsiNOutputs > 2)
         {
            fprintf(fp, "subplot(2,3,3), contourf(x,y,B3_%d)\n", jj+1);
            fwritePlotAxes(fp);
            fprintf(fp,"title(\'Contour Plot for %s\',",
                    outputNames[rsiSet[2]]);
            fprintf(fp, "'FontSize',12,'FontWeight','bold')\n"); 
            fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                    inputNames[iplot1]);
            fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                    inputNames[iplot2]);
         }
         if (rsiNOutputs > 3)
         {
            fprintf(fp, "subplot(2,3,4), contourf(x,y,B4_%d)\n", jj+1);
            fwritePlotAxes(fp);
            fprintf(fp, "title(\'Contour Plot for ");
            fprintf(fp, "%s\','FontSize',12,'FontWeight','bold')\n", 
                    outputNames[rsiSet[3]]);
            fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                    inputNames[iplot1]);
            fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                    inputNames[iplot2]);
         }
         if (rsiNOutputs > 4)
         {
            fprintf(fp, "subplot(2,3,5), contourf(x,y,B5_%d)\n", jj+1);
            fwritePlotAxes(fp);
            fprintf(fp, "title(\'Contour Plot for ");
            fprintf(fp, "%s\','FontSize',12,'FontWeight','bold')\n", 
                    outputNames[rsiSet[4]]);
            fprintf(fp, "xlabel('%s','FontSize',12,'FontWeight','bold')\n",
                    inputNames[iplot1]);
            fprintf(fp, "ylabel('%s','FontSize',12,'FontWeight','bold');\n",
                    inputNames[iplot2]);
         }
         fprintf(fp, "subplot(2,3,6), contourf(x,y,M%d)\n",jj+1);
         fwritePlotAxes(fp);
         fwritePlotXLabel(fp, inputNames[iplot1]);
         fwritePlotYLabel(fp, inputNames[iplot2]);
         fprintf(fp, "title('Intersection: Input %s = %11.4e',",
                 inputNames[iplot3], inputSettings[iplot3]);
         fprintf(fp, "'FontSize',12,'FontWeight','bold')\n");
         fprintf(fp, "colorbar\n");
         fprintf(fp, "colormap(jet)\n");
         fprintf(fp,"pause(1)\n");
      }
      fclose(fp);
      printf("matlabrsi3m.m is now available for response surface and ");
      printf("contour plots\n");

      delete [] rsiSet;
      delete [] faYIn;
      delete faPtr;
      delete [] inputSettings;
   }

   // +++ rssd_ua 
   else if (!strcmp(command, "rssd_ua"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rssd_ua: generate pdf for std. deviations from RS fit\n");
         printf("syntax: rssd_ua (no argument needed).\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data to analyze (load sample first).\n");
         return -1;
      }
      printf("This command works with the following response surfaces:\n");
      printf("1. Linear    regression\n");
      printf("2. Quadratic regression\n");
      printf("3. cubic     regression\n");
      printf("4. quartic   regression\n");
      printf("5. GP1\n");
      printf("6. GP2\n");
      printf("7. MarsBagg\n");
      printf("8. Tree GP\n");
      printf("9. Kriging\n");
      sprintf(pString, "Enter your choice: (1, 2, ..., 9) ");
      int faType = getInt(1, 9, pString);
      if      (faType == 1) faType = PSUADE_RS_REGR1;
      else if (faType == 2) faType = PSUADE_RS_REGR2;
      else if (faType == 3) faType = PSUADE_RS_REGR3;
      else if (faType == 4) faType = PSUADE_RS_REGR4;
      else if (faType == 5) faType = PSUADE_RS_GP1;
      else if (faType == 6) faType = PSUADE_RS_GP2;
      else if (faType == 7) faType = PSUADE_RS_MARSB;
      else if (faType == 8) faType = PSUADE_RS_TGP;
      else if (faType == 9) faType = PSUADE_RS_KR;

      int jplot = 0;
      sprintf(pString, "Enter the output number (1 - %d) : ",nOutputs);
      jplot = getInt(1, nOutputs, pString);
      jplot--;

      printf("rssd_ua: setting up function approximator\n");
      int iOne=1, sInd;
      FuncApprox *faPtr = genFA(faType, nInputs, iOne, nSamples);
      if (faPtr == NULL) {printf("ERROR detected in RS.\n"); return -1;}
      double *faYIn = new double[nSamples];
      for (sInd = 0; sInd < nSamples; sInd++)
         faYIn[sInd] = sampleOutputs[sInd*nOutputs+jplot];
      faPtr->setBounds(iLowerB, iUpperB);
      faPtr->setOutputLevel(outputLevel);
      faPtr->initialize(sampleInputs,faYIn);

      printf("rssd_ua: creating a large sample for constructing PDF\n");
      Sampling *sampPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LPTAU);
      sampPtr->setPrintLevel(0);
      sampPtr->setInputBounds(nInputs, iLowerB, iUpperB);
      sampPtr->setOutputParams(1);
      int count = 100000;
      sampPtr->setSamplingParams(count, -1, 1);
      sampPtr->initialize(0);
      double *tempX  = new double[count*nInputs];
      double *tempY  = new double[count];
      double *tempW  = new double[count];
      int    *states = new int[count];
      sampPtr->getSamples(count, nInputs, 1, tempX, tempY, states);
      faPtr->evaluatePointFuzzy(count, tempX, tempY, tempW);

      int analysisMethod = PSUADE_ANA_MOMENT;
      AnalysisManager *anaManager = new AnalysisManager();
      anaManager->setup(analysisMethod, 0);
      psVector vecUpper, vecLower, vecIn, vecOut;
      vecUpper.load(nInputs, iUpperB);
      vecLower.load(nInputs, iLowerB);
      vecIn.load(count*nInputs, tempX);
      vecOut.load(count, tempW);
      printf("rssd_ua: analyzing the distribution of std deviation\n");
      anaManager->analyze(analysisMethod, count, vecLower, vecUpper,
                          vecIn, vecOut, 0);
      int flag = 1;
      for (ii = 0; ii < count; ii++)
      {
         if (tempY[ii] == 0.0) flag = 0; 
         else                  tempW[ii] /= tempY[ii];
      }
      if (flag == 1)
      { 
         sprintf(pString,"analyze std dev with normalized output (y or n)? ");
         getString(pString, winput);
         if (winput[0] == 'y' ) 
         {
            vecOut.load(count, tempW);
            anaManager->analyze(analysisMethod, count, vecLower, vecUpper,
                                vecIn, vecOut, 0);
         }
      }
      delete anaManager;
      delete sampPtr;
      delete faPtr;
      delete [] tempX;
      delete [] tempY;
      delete [] tempW;
      delete [] states;
      delete [] faYIn;
   }
 
   // +++ rsmeb 
   else if (!strcmp(command, "rsmeb"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rsmeb: main effect analysis on response surface\n");
         printf("syntax: rsmeb (no argument needed)\n");
         printf("This command perform main effect analysis on the response\n");
         printf("surface built from the loaded sample.\n");
         printf("NOTE: This analysis supports other than uniform distributions\n");
         printf("      for the inputs. Simply prescribe the distributions in\n");
         printf("      the data file and turn on use_input_pdfs in ANALYSIS.\n");
         printf("NOTE: This analysis is equivalent to rssobol1 but using\n");
         printf("      a different algorithm.\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data (load sample first).\n");
         return -1;
      }

      int nLHS=200000, nReps=100;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
      outputID = getInt(1, nOutputs, pString);
      outputID--;
      if (psRSExpertMode_ == 1)
      {
         sprintf(pString,
           "Sample size for generating distribution? (10000 - 500000) ");
         nLHS = getInt(10000, 500000, pString);
      }
      sprintf(pString, "How many bootstrapped samples to use (10 - 300) : ");
      int numBS = getInt(10, 300, pString);

      int        iOne=1, iZero=0;
      Sampling   *samPtr;
      printEquals(PL_INFO, 0);
      printf("Phase 1 of 2: create a replicated LH sample\n");
      double *LHSInputs  = new double[nInputs*nLHS];
      double *LHSOutputs = new double[nLHS];
      int    *LHSStates  = new int[nLHS];
      samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_LHS);
      samPtr->setPrintLevel(0);
      samPtr->setInputBounds(nInputs, iLowerB, iUpperB);
      samPtr->setInputParams(nInputs, NULL, NULL, NULL);
      samPtr->setOutputParams(iOne);
      samPtr->setSamplingParams(nLHS, nReps, iZero);
      samPtr->initialize(0);
      samPtr->getSamples(nLHS, nInputs, iOne, LHSInputs,
                         LHSOutputs, LHSStates);
      delete samPtr;
      delete [] LHSStates;

      PDFManager *pdfman;
      psVector   vecIn, vecOut, vecLower, vecUpper;
      psuadeIO->getParameter("ana_use_input_pdfs", pPtr);
      int usePDFs = pPtr.intData_;
      if (usePDFs == 1)
      {
         pdfman = new PDFManager();
         psuadeIO->updateAnalysisSection(-1,-1,-1,0,-1, -1);
         psuadeIO->updateMethodSection(PSUADE_SAMP_MC, -1, -1, -1, -1);
         pdfman->initialize(psuadeIO);
         vecIn.load(nLHS*nInputs, LHSInputs);
         vecOut.setLength(nLHS*nInputs);
         vecUpper.load(nInputs, iUpperB);
         vecLower.load(nInputs, iLowerB);
         pdfman->invCDF(nLHS, vecIn, vecOut, vecLower, vecUpper);
         for (ii = 0; ii < nLHS*nInputs; ii++) LHSInputs[ii] = vecOut[ii];
         delete pdfman;
      }

      FuncApprox *faPtr;
      printEquals(PL_INFO, 0);
      printf("Phase 2 of 2: create response surfaces and run main effects\n");
      printf("NOTE: the response surface type is taken from your data file.\n");
      psuadeIO->getParameter("ana_rstype", pPtr);
      faType = pPtr.intData_;
      if (faType == PSUADE_RS_MARSB) faType = PSUADE_RS_MARS;
      double *tempY;
      if (nOutputs > 1)
      {
         tempY = new double[nSamples];
         for (ss = 0; ss < nSamples; ss++)
            tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
      }
      else tempY = sampleOutputs;

      double *tmpIns  = new double[nSamples*nInputs];
      double *tmpOuts = new double[nSamples];
      int    *mebInds = new int[nSamples];
      int    ind, nSamples2;
      double *meStore = new double[(numBS+2)*nInputs];
      MainEffectAnalyzer *meAnalyzer = new MainEffectAnalyzer();
      pData *pd = NULL;
      aData adata;
      adata.nInputs_ = nInputs;
      adata.nOutputs_ = 1;
      adata.nSamples_ = nLHS;
      adata.outputID_ = 0;
      adata.sampleInputs_ = LHSInputs;
      adata.sampleOutputs_ = LHSOutputs;
      adata.nSubSamples_ = nLHS / nReps;
      adata.iLowerB_ = iLowerB;
      adata.iUpperB_ = iUpperB;
      adata.printLevel_ = 0;
      adata.ioPtr_ = psuadeIO;
      int saveAnaMode = psAnaExpertMode_;
      psAnaExpertMode_ = 0;
      if (psRSExpertMode_ == 1)
      {
         printf("rsmeb INFO: since RS expert mode has been enabled,\n");
         printf("     you should expect to be asked to set RS parameters\n");
         printf("     many times. If you would not prefer to be asked so\n");
         printf("     many times, use config file to set RS parameters.\n");
      }
      for (kk = 0; kk < numBS; kk++)
      {
         printf("rsmeb: iteration %d (of %d)\n", kk+1, numBS);
         for (ss = 0; ss < nSamples; ss++) mebInds[ss] = 0;
         ss = nSamples2 = 0;
         while (ss < nSamples)
         {
            ind = PSUADE_rand() % nSamples;
            if (mebInds[ind] == 0)
            {
               for (ii = 0; ii < nInputs; ii++)
                  tmpIns[nSamples2*nInputs+ii] = sampleInputs[ind*nInputs+ii];
               tmpOuts[nSamples2] = tempY[ind];
               mebInds[ind] = 1;
               nSamples2++;
            }
            ss++;
         }
         faPtr = genFA(faType, nInputs, -1, nSamples2);
         faPtr->setNPtsPerDim(32);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(0);
         faPtr->initialize(tmpIns,tmpOuts);
         faPtr->evaluatePoint(nLHS,LHSInputs,LHSOutputs);
         meAnalyzer->analyze(adata);
         pd = psuadeIO->getAuxData();
         for (ii = 0; ii < nInputs; ii++)
         {
            if (pd->dbleData_ > 0)
                 meStore[kk*nInputs+ii] = pd->dbleArray_[ii]/pd->dbleData_;
            else meStore[kk*nInputs+ii] = pd->dbleArray_[ii];
         }
         pd->clean();
         delete faPtr;
      }
 
      double mean, stdev;
      printf("rsmeb main effect analysis (not scaled by total variance)\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         mean = 0.0;
         for (kk = 0; kk < numBS; kk++)
            mean += meStore[kk*nInputs+ii];
         mean /= numBS;
         meStore[numBS*nInputs+ii] = mean;
         stdev = 0.0;
         for (kk = 0; kk < numBS; kk++)
            stdev += pow(meStore[kk*nInputs+ii]-mean, 2.0);
         stdev = sqrt(stdev/(numBS-1));
         meStore[(numBS+1)*nInputs+ii] = stdev;
         printf("Input %6d = %12.4e (stdev = %12.4e)\n",ii+1,mean,stdev);
      }

      fp = NULL;
      if (psPlotTool_ == 1) fp = fopen("scilabrsmeb.sci","w");
      else                  fp = fopen("matlabrsmeb.m","w");
      if (fp == NULL) printf("ERROR: cannot open plot file.\n");
      else
      {
         strcpy(pString," This file contains main effect ");
         fwriteComment(fp, pString);
         strcpy(pString," with error bars coming from bootstrapping.");
         fwriteComment(fp, pString);
         strcpy(pString," to select the most important ones to display,");
         fwriteComment(fp, pString);
         strcpy(pString," set sortFlag = 1 and set nn to be the number");
         fwriteComment(fp, pString);
         strcpy(pString," of inputs to display.\n");
         fwriteComment(fp, pString);
         fprintf(fp, "sortFlag = 0;\n");
         fprintf(fp, "nn = %d;\n", nInputs);
         fprintf(fp, "Means = [\n");
         for (ii = 0; ii < nInputs; ii++)
            fprintf(fp,"%24.16e\n",meStore[numBS*nInputs+ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "Stds = [\n");
         for (ii = 0; ii < nInputs; ii++)
            fprintf(fp,"%24.16e\n",meStore[(numBS+1)*nInputs+ii]);
         fprintf(fp, "];\n");
         if (inputNames == NULL)
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
               if (inputNames[ii] != NULL) 
                    fprintf(fp,"'%s',",inputNames[ii]);
               else fprintf(fp,"'X%d',",ii+1);
            }
            if (inputNames[nInputs-1] != NULL) 
                 fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
            else fprintf(fp,"'X%d'};\n",nInputs);
         }
         fwriteHold(fp, 0);
         fprintf(fp, "if (sortFlag == 1)\n");
         if (psPlotTool_ == 1)
              fprintf(fp, "  [Means, I2] = gsort(Means);\n");
         else fprintf(fp, "  [Means, I2] = sort(Means,'descend');\n");
         fprintf(fp, "  Stds = Stds(I2);\n");
         fprintf(fp, "  I2 = I2(1:nn);\n");
         fprintf(fp, "  Means = Means(1:nn);\n");
         fprintf(fp, "  Stds = Stds(1:nn);\n");
         fprintf(fp, "  Str  = Str(I2);\n");
         fprintf(fp, "end\n");
         fprintf(fp, "ymin = min(Means-Stds);\n");
         fprintf(fp, "ymax = max(Means+Stds);\n");
         fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
         if (psPlotTool_ == 1) fprintf(fp, "drawlater\n");
         fprintf(fp, "bar(Means,0.8);\n");
         fprintf(fp, "for ii = 1:nn\n");
         fprintf(fp, "   if (ii == 1)\n");
         fwriteHold(fp, 1);
         fprintf(fp, "   end;\n");
         fprintf(fp, "   XX = [ii ii];\n");
         fprintf(fp, "   d1 = Means(ii)-Stds(ii);\n");
         fprintf(fp, "   d2 = Means(ii)+Stds(ii);\n");
         fprintf(fp, "   if (d1 < 0)\n");
         fprintf(fp, "      d1 = 0.0;\n");
         fprintf(fp, "   end;\n");
         fprintf(fp, "   YY = [d1 d2];\n");
         fprintf(fp, "   plot(XX,YY,'-ko','LineWidth',3.0,'MarkerEdgeColor',");
         fprintf(fp, "'k','MarkerFaceColor','g','MarkerSize',13)\n");
         fprintf(fp, "end;\n");
         fwritePlotAxes(fp);
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[0, ymin; nn+1, ymax];\n");
            fprintf(fp, "newtick = a.x_ticks;\n");
            fprintf(fp, "newtick(2) = [1:nn]';\n");
            fprintf(fp, "newtick(3) = Str';\n");
            fprintf(fp, "a.x_ticks = newtick;\n");
            fprintf(fp, "a.x_label.font_size = 3;\n");
            fprintf(fp, "a.x_label.font_style = 4;\n");
         }
         else
         {
            fprintf(fp, "axis([0  nn+1 ymin ymax])\n");
            fprintf(fp, "set(gca,'XTickLabel',[]);\n");
            fprintf(fp, "th=text(1:nn, repmat(ymin-0.05*(ymax-ymin),nn,1),Str,");
            fprintf(fp, "'HorizontalAlignment','left','rotation',90);\n");
            fprintf(fp, "set(th, 'fontsize', 12)\n");
            fprintf(fp, "set(th, 'fontweight', 'bold')\n");
         }
         fwritePlotTitle(fp,"Main Effects (with bootstrap)");
         fwritePlotYLabel(fp, "Main Effects (Normalized)");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "drawnow\n");
            printf("Scilab plot for main effect is in scilabrsmebb.sci.\n");
         }
         else printf("Matlab plot for main effect is in matlabrsmeb.m.\n");
         fclose(fp);
      }
      delete meAnalyzer;
      delete [] tmpIns;
      delete [] tmpOuts;
      delete [] mebInds;
      if (nOutputs > 1) delete [] tempY;
      tempY = NULL;
      faPtr = NULL;
      delete [] LHSInputs;
      delete [] LHSOutputs;
      delete [] meStore;
      psAnaExpertMode_ = saveAnaMode;
   }
 
   // +++ rsieb 
   else if (!strcmp(command, "rsieb"))
   {
      sscanf(lineIn,"%s %s",command,winput);
      if (!strcmp(winput, "-h"))
      {
         printf("rsieb: two-parameter effect analysis on response surface\n");
         printf("syntax: rsieb (no argument needed)\n");
         printf("This command perform two parameter effect analysis on the\n");
         printf("response surface built from the loaded sample.\n");
         printf("NOTE: This analysis supports other than uniform distributions\n");
         printf("      for the inputs. Simply prescribe the distributions in\n");
         printf("      the data file and turn on use_input_pdfs in ANALYSIS.\n");
         printf("NOTE: This analysis is equivalent to rssobol2 but using\n");
         printf("      a different algorithm.\n");
         return 0;
      }
      if (psuadeIO == NULL || sampleOutputs == NULL)
      {
         printf("ERROR: no data (load sample first).\n");
         return -1;
      }
      if (nInputs <= 2)
      {
         printf("INFO: nInputs <=2 -> no point of performing this analysis.\n");
         return -1;
      }

      int nOA, nOA1;
      sprintf(pString, "Enter output number (1 - %d) : ", nOutputs);
      outputID = getInt(1, nOutputs, pString);
      outputID--;
      sprintf(pString, "How many bootstrapped samples to use (1 - 100) : ");
      int numBS = getInt(1, 100, pString);

      int        iOne=1, iZero=0, nReps;
      Sampling   *samPtr;
      printEquals(PL_INFO, 0);
      printf("Phase 1 of 2: create a replicated OA sample\n");
      nOA1  = 293;
      nOA   = nOA1 * nOA1;
      nReps = 100;
      nOA   = nOA * nReps;
      double *OAInputs  = new double[nInputs*nOA];
      double *OAOutputs = new double[nOA];
      int    *OAStates  = new int[nOA];
      samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_OA);
      samPtr->setPrintLevel(outputLevel);
      samPtr->setInputBounds(nInputs, iLowerB, iUpperB);
      samPtr->setInputParams(nInputs, NULL, NULL, NULL);
      samPtr->setOutputParams(iOne);
      samPtr->setSamplingParams(nOA, nReps, iZero);
      samPtr->initialize(0);
      samPtr->getSamples(nOA, nInputs, iOne, OAInputs,
                         OAOutputs, OAStates);
      delete samPtr;
      delete [] OAStates;

      PDFManager *pdfman;
      psVector   vecIn, vecOut, vecLower, vecUpper;
      psuadeIO->getParameter("ana_use_input_pdfs", pPtr);
      int usePDFs = pPtr.intData_;
      if (usePDFs == 1)
      {
         pdfman = new PDFManager();
         psuadeIO->updateAnalysisSection(-1,-1,-1,0,-1, -1);
         psuadeIO->updateMethodSection(PSUADE_SAMP_MC, -1, -1, -1, -1);
         pdfman->initialize(psuadeIO);
         vecIn.load(nOA*nInputs, OAInputs);
         vecOut.setLength(nOA*nInputs);
         vecUpper.load(nInputs, iUpperB);
         vecLower.load(nInputs, iLowerB);
         pdfman->invCDF(nOA, vecIn, vecOut, vecLower, vecUpper);
         for (ii = 0; ii < nOA*nInputs; ii++) OAInputs[ii] = vecOut[ii];
         delete pdfman;
      }

      FuncApprox *faPtr;
      printEquals(PL_INFO, 0);
      printf("Phase 2 of 2: create response surfaces and run main effects\n");
      printf("NOTE: the response surface type is taken from your data file.\n");
      psuadeIO->getParameter("ana_rstype", pPtr);
      faType = pPtr.intData_;
      if (faType == PSUADE_RS_MARSB) faType = PSUADE_RS_MARS;
      double *tempY;
      if (nOutputs > 1)
      {
         tempY = new double[nSamples];
         for (ss = 0; ss < nSamples; ss++)
            tempY[ss] = sampleOutputs[ss*nOutputs+outputID];
      }
      else tempY = sampleOutputs;

      int    ind, nSamples2, ii1, ii2, jj1, jj2, kk2, bin1, bin2;
      int    *iebInds  = new int[nSamples];
      int    *iecounts = new int[nOA1*nOA1];
      double width1, width2, iemean, ievar;
      double *tmpIns  = new double[nSamples*nInputs];
      double *tmpOuts = new double[nSamples];
      double *iemeans = new double[nOA1*nOA1];
      double *ievars  = new double[nOA1*nOA1];
      double *ieStore = new double[(numBS+2)*nInputs*nInputs];
      for (kk = 0; kk < numBS; kk++)
      {
         printf("rsieb: iteration %d (of %d)\n", kk+1, numBS);
         if (numBS > 1)
         {
            for (ss = 0; ss < nSamples; ss++) iebInds[ss] = 0;
            ss = nSamples2 = 0;
            while (ss < nSamples)
            {
               ind = PSUADE_rand() % nSamples;
               if (iebInds[ind] == 0)
               {
                  for (ii = 0; ii < nInputs; ii++)
                     tmpIns[nSamples2*nInputs+ii] = sampleInputs[ind*nInputs+ii];
                  tmpOuts[nSamples2] = tempY[ind];
                  iebInds[ind] = 1;
                  nSamples2++;
               }
               ss++;
            }
         }
         else
         {
            nSamples2 = nSamples;
            for (ss = 0; ss < nSamples; ss++) 
            {
               for (ii = 0; ii < nInputs; ii++) 
                  tmpIns[ss*nInputs+ii] = sampleInputs[ss*nInputs+ii];
               tmpOuts[ss] = tempY[ss];
            }
         }
         faPtr = genFA(faType, nInputs, -1, nSamples2);
         faPtr->setNPtsPerDim(32);
         faPtr->setBounds(iLowerB, iUpperB);
         faPtr->setOutputLevel(0);
         faPtr->initialize(tmpIns,tmpOuts);
         faPtr->evaluatePoint(nOA,OAInputs,OAOutputs);
         for (ii1 = 0; ii1 < nInputs; ii1++)
         {
            ieStore[kk*nInputs*nInputs+ii1*nInputs+ii1] = 0.0;
            width1 = (iUpperB[ii1] - iLowerB[ii1]) / (nOA1 - 1);
            for (ii2 = ii1+1; ii2 < nInputs; ii2++)
            {
               width2 = (iUpperB[ii2] - iLowerB[ii2]) / (nOA1 - 1);
               for (kk2 = 0; kk2 < nOA1*nOA1; kk2++) 
               {
                  iecounts[kk2] = 0;
                  iemeans[kk2] = 0.0;
                  ievars[kk2] = 0.0;
               }
               for (kk2 = 0; kk2 < nOA; kk2++)
               {
                  bin1 = (int) ((OAInputs[kk2*nInputs+ii1]-
                                 iLowerB[ii1]+1.0e-12)/width1);
                  bin2 = (int) ((OAInputs[kk2*nInputs+ii2]-
                                 iLowerB[ii2]+1.0e-12)/width2);
                  iemeans[bin1*nOA1+bin2] += OAOutputs[kk2];
                  iecounts[bin1*nOA1+bin2]++;
               }
               for (kk2 = 0; kk2 < nOA1*nOA1; kk2++)
                  if (iecounts[kk2] > 0) iemeans[kk2] /= iecounts[kk2];
               for (kk2 = 0; kk2 < nOA; kk2++)
               {
                  bin1 = (int) ((OAInputs[kk2*nInputs+ii1]-
                                 iLowerB[ii1]+1.0e-12)/width1);
                  bin2 = (int) ((OAInputs[kk2*nInputs+ii2]-
                                 iLowerB[ii2]+1.0e-12)/width2);
                  ievars[bin1*nOA1+bin2] += 
                          pow(OAOutputs[kk2]-iemeans[bin1*nOA1+bin2],2.0);
               }
               for (kk2 = 0; kk2 < nOA1*nOA1; kk2++)
                  if (iecounts[kk2] > 0) ievars[kk2] /= iecounts[kk2];
               iemean = 0.0;
               for (kk2 = 0; kk2 < nOA1*nOA1; kk2++) iemean += iemeans[kk2];
               iemean /= (double) (nOA1 * nOA1);
               ievar = 0.0;
               for (kk2 = 0; kk2 < nOA1*nOA1; kk2++) 
                  ievar += pow(iemeans[kk2]-iemean,2.0);
               ievar /= (double) (nOA1 * nOA1);
               iemean = 0.0;
               for (kk2 = 0; kk2 < nOA1*nOA1; kk2++) iemean += ievars[kk2];
               iemean /= (double) (nOA1 * nOA1);
               ievar -= iemean / nReps;
               if (ievar < 0) ievar = 0;
               ieStore[kk*nInputs*nInputs+ii1*nInputs+ii2] = ievar;
               ieStore[kk*nInputs*nInputs+ii2*nInputs+ii1] = ievar;
            }
         }
         delete faPtr;
      }
      delete [] iecounts;
      delete [] iemeans;
      delete [] ievars;
      delete [] tmpIns;
      delete [] tmpOuts;
      delete [] iebInds;
  
      double mean, stdev;
      printf("rsieb main effect analysis (not scaled by total variance)\n");
      for (ii = 0; ii < nInputs; ii++)
      {
         for (jj = ii+1; jj < nInputs; jj++)
         {
            mean = 0.0;
            for (kk = 0; kk < numBS; kk++)
               mean += ieStore[kk*nInputs*nInputs+ii*nInputs+jj];
            mean /= numBS;
            ieStore[numBS*nInputs*nInputs+ii*nInputs+jj] = mean;
            ieStore[numBS*nInputs*nInputs+jj*nInputs+ii] = 0.0;
            stdev = 0.0;
            if (numBS > 1)
            {
               for (kk = 0; kk < numBS; kk++)
                  stdev += pow(ieStore[kk*nInputs*nInputs+ii*nInputs+jj]-mean,2.0);
               stdev = sqrt(stdev/(numBS-1));
            }
            ieStore[(numBS+1)*nInputs*nInputs+ii*nInputs+jj] = stdev;
            ieStore[(numBS+1)*nInputs*nInputs+jj*nInputs+ii] = 0.0;
            printf("Input %6d %d = %12.4e (stdev = %12.4e)\n",ii+1,
                   jj+1,mean,stdev);
         }
      }

      if (psPlotTool_ == 1)
      {
         fp = fopen("scilabrsieb.sci", "w");
         if (fp == NULL) printf("ERROR : cannot open file scilabrsieb.sci\n");
         else
         {
            fprintf(fp,"// This file contains 2nd order sensitivity indices\n");
            fprintf(fp,"// set sortFlag = 1 and set nn to be the number\n");
            fprintf(fp,"// of inputs to display.\n");
         }
      }
      else
      {
         fp = fopen("matlabrsieb.m", "w");
         if (fp == NULL) printf("ERROR : cannot open file matlabrsieb.sci\n");
         else
         {
            fprintf(fp, "%% This file contains 2nd order sensitivity indices\n");
            fprintf(fp, "%% set sortFlag = 1 and set nn to be the number\n");
            fprintf(fp, "%% of inputs to display.\n");
         }
      }
      if (fp != NULL) 
      {
         fprintf(fp, "sortFlag = 0;\n");
         fprintf(fp, "nn = %d;\n", nInputs);
         fprintf(fp, "Means = [\n");
         for (ii = 0; ii < nInputs*nInputs; ii++) 
            fprintf(fp,"%24.16e\n", ieStore[numBS*nInputs*nInputs+ii]);
         fprintf(fp, "];\n");
         fprintf(fp, "Stds = [\n");
         for (ii = 0; ii < nInputs*nInputs; ii++) 
            fprintf(fp,"%24.16e\n", ieStore[(numBS+1)*nInputs*nInputs+ii]);
         fprintf(fp, "];\n");
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
               if (inputNames[ii] != NULL) 
                  fprintf(fp,"'%s',",inputNames[ii]);
               else fprintf(fp,"'X%d',",ii+1);
            }
            if (inputNames[nInputs-1] != NULL)
               fprintf(fp,"'%s'};\n",inputNames[nInputs-1]);
            else fprintf(fp,"'X%d'};\n",nInputs);
         }
         fwriteHold(fp, 0);
         fprintf(fp, "ymin = min(Means-Stds);\n");
         fprintf(fp, "ymax = max(Means+Stds);\n");
         fprintf(fp, "h2 = 0.05 * (ymax - ymin);\n");
         if (psPlotTool_ == 1)
         {
            fprintf(fp, "nn    = %d;\n",nInputs);
            fprintf(fp, "Means = matrix(Means, nn, nn);\n");
            fprintf(fp, "Means = Means';\n");
            fprintf(fp, "Stds  = matrix(Stds, nn, nn);\n");
            fprintf(fp, "Stds  = Stds';\n");
            fprintf(fp, "drawlater\n");
            fprintf(fp, "hist3d(Means);\n");
            fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
            fprintf(fp, "a=gca();\n");
            fprintf(fp, "a.data_bounds=[0, 0, 0; nn, nn+1, ymax];\n");
            fprintf(fp, "newtick = a.x_ticks;\n");
            fprintf(fp, "newtick(2) = [1:nn]';\n");
            fprintf(fp, "newtick(3) = Str';\n");
            fprintf(fp, "a.x_ticks = newtick;\n");
            fprintf(fp, "a.x_label.font_size = 3;\n");
            fprintf(fp, "a.x_label.font_style = 4;\n");
            fprintf(fp, "a.y_ticks = newtick;\n");
            fprintf(fp, "a.y_label.font_size = 3;\n");
            fprintf(fp, "a.y_label.font_style = 4;\n");
            fprintf(fp, "a.rotation_angles = [5 -70];\n");
            fprintf(fp, "drawnow\n");
         }
         else
         {
            fprintf(fp, "nn    = %d;\n",nInputs);
            fprintf(fp, "Means = reshape(Means, nn, nn);\n");
            fprintf(fp, "Means = Means';\n");
            fprintf(fp, "Stds  = reshape(Stds, nn, nn);\n");
            fprintf(fp, "Stds  = Stds';\n");
            fprintf(fp, "hh = bar3(Means,0.8);\n");
            fprintf(fp, "alpha = 0.2;\n");
            fprintf(fp, "set(hh,'FaceColor','b','facea',alpha);\n");
            fprintf(fp, "Lstds = Means - Stds;\n");
            fprintf(fp, "Ustds = Means + Stds;\n");
            fprintf(fp, "[X,Y] = meshgrid(1:nn,1:nn);\n");
            fwriteHold(fp, 1);
            fprintf(fp, "for k = 1:nn\n");
            fprintf(fp, "  for l = k:nn\n");
            fprintf(fp, "    mkl = Means(k,l);\n");
            fprintf(fp, "    ukl = Ustds(k,l);\n");
            fprintf(fp, "    lkl = Lstds(k,l);\n");
            fprintf(fp, "    if (mkl > .02 & (ukl-lkl)/mkl > .02)\n");
            fprintf(fp, "      xkl = [X(k,l), X(k,l)];\n");
            fprintf(fp, "      ykl = [Y(k,l), Y(k,l)];\n");
            fprintf(fp, "      zkl = [lkl, ukl];\n");
            fprintf(fp, "      plot3(xkl,ykl,zkl,'-mo',...\n");
            fprintf(fp, "        'LineWidth',5,'MarkerEdgeColor','k',...\n");
            fprintf(fp, "        'MarkerFaceColor','k','MarkerSize',10);\n");
            fprintf(fp, "    end\n");
            fprintf(fp, "  end\n");
            fprintf(fp, "end\n");
            fwriteHold(fp, 0);
            fprintf(fp, "axis([0.5 nn+0.5 0.5 nn+0.5 0 ymax])\n");
            fprintf(fp, "set(gca,'XTickLabel',Str);\n");
            fprintf(fp, "set(gca,'YTickLabel',Str);\n");
            fprintf(fp, "set(gca, 'fontsize', 12)\n");
            fprintf(fp, "set(gca, 'fontweight', 'bold')\n");
            fprintf(fp, "set(gca, 'linewidth', 2)\n");
         }
         fwritePlotAxes(fp);
         fwritePlotTitle(fp,"1st+2nd Order Sensitivity Indices (with bootstrap)");
         fwritePlotZLabel(fp, "2nd Order Sensitivity Indices (Normalized)");
         fwritePlotXLabel(fp, "Inputs");
         fwritePlotYLabel(fp, "Inputs");
         fclose(fp);
         if (psPlotTool_ == 1)
              printf("rsieb plot file = scilabrsieb.sci\n");
         else printf("rsieb plot file = matlabrsieb.m\n");
      }
      if (nOutputs > 1) delete [] tempY;
      tempY = NULL;
      faPtr = NULL;
      delete [] OAInputs;
      delete [] OAOutputs;
      delete [] ieStore;
   }
   return 0;
}

