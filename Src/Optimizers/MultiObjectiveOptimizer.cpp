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
// Functions for the class BobyqaOptimizer
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "MultiObjectiveOptimizer.h"
#include "Sampling.h"
#include "PsuadeData.h"
#include "PsuadeUtil.h"
#include "Psuade.h"

extern "C" void bobyqa_(int *,int *, double *, double *, double *, double *,
                        double *, int *, int *, double*);

void   *psMOOObj_=NULL;
int    MO_NumVars_=0;
double *MO_Variables_=NULL;
double *MO_Outputs_=NULL;
char   MO_PythonFile_[1001];

#define PABS(x)  ((x) > 0 ? x : -(x))
// ************************************************************************
// ************************************************************************
// resident function to perform evaluation 
// ------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" 
{
#endif
   void *moobobyqaevalfunc_(int *nInps, double *XValues, double *YValue)
   {
      int    ii, funcID, nInputs, nOutputs;
      double *localY, ddata, dsum;
      char   fName[500], runLine[500];
      oData  *odata;
      FILE   *fp = NULL;

      nInputs = (*nInps);
      odata    = (oData *) psMOOObj_;
      nOutputs = odata->nOutputs_;
      localY   = (double *) malloc(nOutputs * sizeof(double));

      funcID = odata->numFuncEvals_;
      odata->funcIO_->evaluate(funcID,nInputs,XValues,nOutputs,localY,0);
      funcID = odata->numFuncEvals_++;
      if (strcmp(MO_PythonFile_, "NULL"))
      {
         sprintf(fName, "MOO_params.in.%d", funcID);
         fp = fopen(fName, "w");
         if (fp == NULL)
         {
            printf("MultiObjectOptimizer ERROR: cannot open file %s\n",fName);
            exit(1);
         }
         fprintf(fp, "%d\n", nOutputs+MO_NumVars_);
         for (ii = 0; ii < nInputs; ii++) fprintf(fp,"%e\n", MO_Variables_[ii]);
         for (ii = 0; ii < MO_NumVars_; ii++) fprintf(fp,"%e\n", localY[ii]);
         fprintf(fp, "# line 1: number of design variables + number of outputs\n");
         fprintf(fp, "# line 2-: design variables\n");
         fprintf(fp, "# line --: output variables\n");
         fclose(fp);
         sprintf(runLine, "%s MOO_params.in.%d MOO_params.out.%d",
                 MO_PythonFile_, funcID, funcID);
         system(runLine);
         sprintf(fName, "MOO_params.out.%d", funcID);
         fp = fopen(fName, "r");
         if (fp == NULL)
         {
            printf("MultiObjectOptimizer ERROR: cannot open file %s\n",fName);
            exit(1);
         }
         fscanf(fp, "%lg", &ddata);
      }
      else
      {
         if (nOutputs != MO_NumVars_+1)
         {
            printf("MultiObjectOptimizer ERROR: nOutputs != numVars+1\n");
            exit(1);
         }
         if (MO_NumVars_ == nOutputs-1)
         {
            ddata = 0.0;
            dsum  = 0.0;
            for (ii = 1; ii <= MO_NumVars_; ii++)
            {
               ddata += MO_Variables_[ii-1] * localY[ii];
               dsum  += MO_Variables_[ii-1];
            }
            ddata += (1.0 - dsum) * localY[0];
         }
         else
         {
            ddata = 0.0;
            for (ii = 0; ii < MO_NumVars_; ii++)
               ddata += MO_Variables_[ii] * localY[ii];
         }
      }
      (*YValue) = ddata;

      if ((*YValue) < odata->optimalY_)
      {
         odata->optimalY_ = (*YValue);
         for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = XValues[ii];
         for (ii = 0; ii < nOutputs; ii++) MO_Outputs_[ii] = localY[ii];
      }
      free(localY);
      if(fp != NULL) fclose(fp);
      return NULL;
   }
#ifdef __cplusplus
}
#endif

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
MultiObjectiveOptimizer::MultiObjectiveOptimizer()
{
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
MultiObjectiveOptimizer::~MultiObjectiveOptimizer()
{
}

// ************************************************************************
// optimize
// ------------------------------------------------------------------------
void MultiObjectiveOptimizer::optimize(oData *odata)
{
   int    nInputs, printLevel=0, ii, maxfun, currDriver, nPts=0, nn;
   int    nOutputs, iOne=1, length, numVars, nSamples, *sampleStat, pLevel;
   int    resolution;
   double *XValues, rhobeg=1.0, rhoend=1.0e-4, dtemp, *workArray;
   double *sampleIns, *sampleOut, *lbounds, *ubounds, *MO_OptX, MO_OptY;
   char   filename[501], lineIn[501], pString[501];
   string   sfname, iline;
   size_t   compFlag;
   ifstream ifile, ifile2;
   Sampling *samPtr;

   printLevel = odata->outputLevel_;
   printAsterisks(PL_INFO, 0);
   printf("Surrogate-based Multi-objective optimization: \n");
   nInputs  = odata->nInputs_;
   nOutputs = odata->nOutputs_;
   printEquals(PL_INFO, 0);
   printf("The multi-objective function is built from the sample\n");
   printf("outputs (currently there are %d outputs).\n",nOutputs);
   printf("For this optimizer to build a general multi-objective\n");
   printf("function, a configuration file is required from users.\n");
   printf("The configuration file should be in the following format:\n");
   printf("line 1: PSUADE_BEGIN\n");
   printf("line 2: number of variables in the multi-objective function\n");
   printf("line 3: 1  lbound ubound <lower and upper bounds of variable 1>\n");
   printf("line 4: 2  lbound ubound <lower and upper bounds of variable 2>\n");
   printf("....\n");
   printf("line n: name of a python file to evaluate the objective function\n");
   printf("line n+1: PSUADE_END\n");
   printf("Note: If the objective function is just a linear combination\n");
   printf("      of the outputs, the python function line should be \n");
   printf("      replaced by a 'NULL', and the number of design variables\n");
   printf("      should be nOutputs-1 (since sum of weights=1).\n");
   printf("An Example: \n");
   printDashes(PL_INFO, 0);
   printf("PSUADE_BEGIN\n");
   printf("2\n");
   printf("1 0 1\n");
   printf("2 0 1\n");
   printf("objfcn.py (NULL if the objective function is a linear combination)\n");
   printf("PSUADE_END\n");
   printDashes(PL_INFO, 0);
   printf("Note: the optimizer will evaluate the multi-objective function\n");
   printf("      by using the calling sequence:\n");
   printf("          <pythonFile> <paramFile> <objFile>\n");
   printf("where:\n");
   printf("  <paramFile> contains a sample point to evaluate the function\n");
   printf("  <objFile>   the function value written by the python file\n");
   printf("NOTE: MAKE SURE the <pythonFile> HAS EXECUTE PERMISSION.\n");
   printEquals(PL_INFO, 0);
   printf("Enter the name of the configuration file: ");
   cin >> sfname;
   fgets(lineIn, 500, stdin);
   length = sfname.size();

   if (length < 500)
   {
      sfname.copy(filename, length, 0);
      filename[length] = '\0';
      ifile.open(filename);
      if (! ifile.is_open())
      {
         printf("MOO ERROR: cannot open configuration file = %s\n",filename);
         return;
      }
   }
   else
   {
      printf("MOO ERROR: configuration file name too long.\n");
      return;
   }
   getline (ifile, iline);
   compFlag = iline.compare("PSUADE_BEGIN");
   if (compFlag == 0)
   {
      ifile >> numVars;
      if (numVars <= 0)
      {
         printf("MOO configuration file ERROR: numVars <= 0\n");
         ifile.close();
         return;
      }
      lbounds = new double[numVars];
      ubounds = new double[numVars];
      for (ii = 0; ii < numVars; ii++)
      {
         ifile >> nn;
         if (nn != (ii+1))
         {
            printf("MOO config file ERROR: variable index mismatch (%d != %d).\n",
                   nn, ii+1);
            ifile.close();
            return;
         }
         ifile >> lbounds[ii];
         ifile >> ubounds[ii];
         if (lbounds[ii] > ubounds[ii])
         {
            printf("MOO config file ERROR: lbound > ubound (input %d).\n",ii+1);
            ifile.close();
            return;
         }
      } 
      getline (ifile, iline);
      ifile >> MO_PythonFile_;
      if (strcmp(MO_PythonFile_, "NULL"))
      {
         ifile2.open(MO_PythonFile_);
         if (! ifile2.is_open())
         {
            printf("MOO config file ERROR: python file not found.\n");
            ifile.close();
            return;
         }
         ifile2.close();
      }
      else
      {
         strcpy(MO_PythonFile_, "NULL");
         if (numVars != nOutputs-1)
         {
            printf("MOO config file ERROR: if no python file, the number of\n");
            printf("   of variables in the configuration file should be equal\n");
            printf("   to nOutputs-1, which is %d\n", nOutputs-1);
            return;
         }
      }
   }
   else
   {
      printf("MOO config file ERROR: PSUADE_BEGIN not found.\n");
      ifile.close();
      return;
   }
   getline (ifile, iline);
   getline (ifile, iline);
   compFlag = iline.compare("PSUADE_END");
   if (compFlag != 0)
   {
      printf("MOO config file ERROR: PSUADE_END not found.\n");
      ifile.close();
      return;
   }
   ifile.close();
   if (numVars > 4) 
   {
      printf("MOO currently cannot handle numVars>4.\n");
      return;
   }

   printf("Next, a full factorial design will be generated to explore\n");
   printf("the design variable space (of dimension %d). Please enter\n",
          numVars);
   printf("the sample resolution n. For example, if the resolution is\n");
   printf("10 for 3 design variables, the total number of optimizations\n");
   printf("to be performed will be 10x10x10=1000.\n");
   if (numVars == 1)
   {
      sprintf(pString,
              "Enter the desired resolution (>2,<100,suggested: 11): ");
      resolution = getInt(3, 99, pString);
   }
   else if (numVars == 2)
   {
      sprintf(pString,
              "Enter the desired resolution (>2,<50,suggested: 11): ");
      resolution = getInt(3, 49, pString);
   }
   else
   {
      sprintf(pString,
              "Enter the desired resolution (>2,<20,suggested: 11): ");
      resolution = getInt(3, 19, pString);
   }
   samPtr = (Sampling *) SamplingCreateFromID(PSUADE_SAMP_FACT);
   samPtr->setPrintLevel(0);
   samPtr->setInputBounds(numVars, lbounds, ubounds);
   samPtr->setOutputParams(iOne);
   nSamples = 1;
   for (ii = 0; ii < numVars; ii++) nSamples *= resolution;
   samPtr->setSamplingParams(nSamples, -1, 0);
   samPtr->initialize(0);
   nSamples = samPtr->getNumSamples();
   sampleIns  = new double[nSamples * numVars];
   sampleOut  = new double[nSamples * (nInputs+nOutputs+1)];
   sampleStat = new int[nSamples];
   samPtr->getSamples(nSamples,numVars,iOne,sampleIns,sampleOut,sampleStat);
   delete samPtr;

   for (ii = 0; ii < nInputs; ii++) odata->optimalX_[ii] = 0.0;
   maxfun = odata->maxFEval_;
   XValues = new double[nInputs+1];
   for (ii = 1; ii < nInputs; ii++) 
   {
      dtemp = odata->upperBounds_[ii] - odata->lowerBounds_[ii];
      if (dtemp < rhobeg) rhobeg = dtemp;
   }
   rhobeg = odata->upperBounds_[0] - odata->lowerBounds_[0];
   rhobeg *= 0.5;
   rhoend = rhobeg * odata->tolerance_;
   if (rhobeg < rhoend)
   {
      printf("MOO WARNING: tolerance too large.\n");
      printf("             tolerance reset to 1.0e-6.\n");
      rhoend = rhobeg * 1.0e-6;
   }
   currDriver = odata->funcIO_->getDriver();
   odata->funcIO_->setDriver(1);
   psMOOObj_= (void *) odata;
   if (printLevel > 0)
   {
      printf("MultiObj optimizer: max fevals = %d\n", odata->maxFEval_);
      printf("MultiObj optimizer: tolerance  = %e\n", odata->tolerance_);
   }
   nPts = (nInputs + 1) * (nInputs + 2) / 2;
   workArray = new double[(nPts+5)*(nPts+nInputs)+3*nInputs*(nInputs+5)/2+1];
   MO_NumVars_ = numVars;
   MO_Variables_ = new double[MO_NumVars_];
   MO_Outputs_ = new double[nOutputs];
   odata->numFuncEvals_ = 0;
   MO_OptX = new double[numVars];
   MO_OptY = PSUADE_UNDEFINED;

   MO_OptY = PSUADE_UNDEFINED;
#if 1
   for (nn = 0; nn < nSamples; nn++)
   {
      for (ii = 0; ii < numVars; ii++)
         MO_Variables_[ii] = sampleIns[nn*numVars+ii];
      for (ii = 0; ii < nInputs; ii++) XValues[ii] = odata->initialX_[ii];
      odata->optimalY_ = 1.0e50;
#ifdef HAVE_BOBYQA
      pLevel = 9999;
      bobyqa_(&nInputs, &nPts, XValues, odata->lowerBounds_,
              odata->upperBounds_, &rhobeg, &rhoend, &pLevel, &maxfun, 
              workArray);
#else
      printf("ERROR: Bobyqa optimizer not installed.\n");
      exit(1);
#endif
      sampleOut[nn*(nInputs+nOutputs+1)] = odata->optimalY_;
      for (ii = 0; ii < nInputs; ii++)
         sampleOut[nn*(nInputs+nOutputs+1)+ii+1] = odata->optimalX_[ii];
      for (ii = 0; ii < nOutputs; ii++)
         sampleOut[nn*(nInputs+nOutputs+1)+nInputs+ii+1] = MO_Outputs_[ii];
      if (printLevel > 2)
      {
         printf("Iteration %5d (%5d) : inputs =", nn+1,nSamples);
         for (ii = 0; ii < numVars; ii++)
            printf(" %e", sampleIns[nn*numVars+ii]);
         printf(", min = %e\n",odata->optimalY_);
      }
      if (odata->optimalY_ < MO_OptY)
      {
         MO_OptY = odata->optimalY_;
         for (ii = 0; ii < numVars; ii++)
            MO_OptX[ii] = sampleIns[nn*numVars+ii];
      }
   }
#else
   int kk, numPCEs, **PCEs, pOrder, flag, flag, ind, isum, itmp, ind2;
   pOrder = (resolution - 1) * numVars;
   GenPermutations(numVars, pOrder, &numPCEs, &PCEs); 
   ind = 1;
   for (nn = 1; nn < 3*(resolution-1); nn++)
   {
      isum = nn;
      ind2 = ind;
      while (isum == nn && ind2 < numPCEs-1)
      {
         ind2++;
         isum = PCEs[ind2][0];
         for (ii = 1; ii < numVars; ii++) isum += PCEs[ind2][ii];
         if (isum != nn) break;
      }
      if (ind2 == (numPCEs-1)) ind2++;
      for (kk = ind; kk < (ind+ind2)/2; kk++)
      {
         for (ii = 0; ii < numVars; ii++)
         {
            itmp = PCEs[kk][ii];
            PCEs[kk][ii] = PCEs[ind2+ind-kk-1][ii];
            PCEs[ind2+ind-kk-1][ii] = itmp;
         }
      }
      ind = ind2;
   }
   for (nn = 0; nn < numPCEs; nn++)
   {
      flag = 1;
      for (ii = 0; ii < numVars; ii++)
         if (PCEs[nn][ii] >= resolution) flag = 0;
      if (flag == 1)
      {
         ind = PCEs[nn][numVars-1];
         for (ii = numVars-2; ii >= 0; ii--)
            ind = ind * resolution + PCEs[nn][ii]; 
         for (ii = 0; ii < numVars; ii++)
            MO_Variables_[ii] = sampleIns[ind*numVars+ii];
         if (nn == 0)
         {
            for (ii = 0; ii < nInputs; ii++)
               XValues[ii] = odata->initialX_[ii];
         }
         else
         {
            for (ii = 0; ii < nInputs; ii++)
               XValues[ii] = odata->optimalX_[ii];
         }
         odata->optimalY_ = 1.0e50;
#ifdef HAVE_BOBYQA
         bobyqa_(&nInputs, &nPts, XValues, odata->lowerBounds_,
                 odata->upperBounds_, &rhobeg, &rhoend, &pLevel, 
                 &maxfun, workArray);
#else
         printf("ERROR: Bobyqa optimizer not installed.\n");
         exit(1);
#endif
         sampleOut[ind*(nInputs+nOutputs+1)] = odata->optimalY_;
         for (ii = 0; ii < nInputs; ii++)
            sampleOut[ind*(nInputs+nOutputs+1)+ii+1] = odata->optimalX_[ii];
         for (ii = 0; ii < nOutputs; ii++)
            sampleOut[ind*(nInputs+nOutputs+1)+nInputs+ii+1] = MO_Outputs_[ii];
         if (odata->optimalY_ < MO_OptY)
         {
            MO_OptY = odata->optimalY_;
            for (ii = 0; ii < numVars; ii++)
               MO_OptX[ii] = sampleIns[ind*numVars+ii];
         }
      }
   }
   for (ii = 0; ii < numPCEs; ii++) delete [] PCEs[ii];
   delete [] PCEs;
#endif
   for (ii = 0; ii < numVars; ii++)
      printf("MOO OptimalX %2d = %e\n", ii+1, MO_OptX[ii]);
   printf("MOO OptimalY    = %e\n", MO_OptY);
   printf("MOO nFuncEval   = %d\n", odata->numFuncEvals_);

   PsuadeData *ioPtr;
   char       **iNames, **oNames;
   iNames = new char*[numVars];
   for (ii = 0; ii < numVars; ii++)
   {
      iNames[ii] = new char[100];
      sprintf(iNames[ii], "X%d", ii+1);
   }
   oNames = new char*[nInputs+nOutputs+1];
   for (ii = 0; ii < nInputs+nOutputs+1; ii++)
   {
      oNames[ii] = new char[100];
      if (ii == 0)            sprintf(oNames[ii], "Y");
      else if (ii <= nInputs) sprintf(oNames[ii], "X%d", ii);
      else                    sprintf(oNames[ii], "F%d", ii);
   }
   for (ii = 0; ii < nSamples; ii++) sampleStat[ii] = 1;
   ioPtr = new PsuadeData();
   ioPtr->updateInputSection(nSamples, numVars, NULL, lbounds,
                             ubounds, sampleIns, iNames,NULL,NULL,NULL,NULL);
   ioPtr->updateOutputSection(nSamples,nInputs+nOutputs+1,sampleOut,
                              sampleStat,oNames);
   ioPtr->updateMethodSection(PSUADE_SAMP_MC, nSamples, 1, -1, -1);
   ioPtr->writePsuadeFile("psuade_moo_sample",0);
   printf("The moo sample is in file psuade_moo_sample.\n");
   printf("The optimal values are in the first output.\n");
   printf("The rest of the outputs in the file are inputs and outputs.\n");
   printf("Use write_std (output 1) and use matlab to visualize for 1D.\n");
   printf("Use rawi2 or rawi3 (output 1) command to generate 2D/3D plots.\n");

   for (ii = 0; ii < numVars; ii++) delete [] iNames[ii];
   delete [] iNames;
   for (ii = 0; ii < nInputs+nOutputs+1; ii++) delete [] oNames[ii];
   delete [] oNames;
   odata->funcIO_->setDriver(currDriver);
   delete [] XValues;
   delete [] workArray;
   delete [] sampleIns;
   delete [] sampleOut;
   delete [] sampleStat;
   delete [] lbounds;
   delete [] ubounds;
   delete [] MO_Variables_;
   delete [] MO_OptX;
   delete ioPtr;
   MO_Variables_ = NULL;
   MO_NumVars_ = 0;
}

// *************************************************************************
// generate all combinations of a multivariate PCE expansion
// This code is a direct translation from Burkardt's matlab code)
// -------------------------------------------------------------------------
int MultiObjectiveOptimizer::GenPermutations(int nInputs, int pOrder,
                                            int *nPerms, int ***pPerms)
{
   int  ii, kk, numPerms, orderTmp, rvTmp, **pcePerms;

   numPerms = 1;
   for (ii = nInputs+pOrder; ii > pOrder; ii--) numPerms *= ii;
   for (ii = 2; ii <= nInputs; ii++) numPerms /= ii;

   pcePerms = new int*[numPerms];
   for (ii = 0; ii < numPerms; ii++) pcePerms[ii] = new int[nInputs];

   numPerms = 0;
   for (kk = 0; kk <= pOrder; kk++)
   {
      orderTmp = kk;
      rvTmp = 0;
      pcePerms[numPerms][0] = orderTmp;
      for (ii = 1; ii < nInputs; ii++) pcePerms[numPerms][ii] = 0;
      while (pcePerms[numPerms][nInputs-1] != kk)
      {
         numPerms++;
         for (ii = 0; ii < nInputs; ii++)
            pcePerms[numPerms][ii] = pcePerms[numPerms-1][ii];
         if (orderTmp > 1) rvTmp = 1;
         else              rvTmp++;
         pcePerms[numPerms][rvTmp-1] = 0;
         orderTmp = pcePerms[numPerms-1][rvTmp-1];
         pcePerms[numPerms][0] = orderTmp - 1;
         pcePerms[numPerms][rvTmp] = pcePerms[numPerms-1][rvTmp] + 1;
      }
      numPerms++;
   }
   (*nPerms) = numPerms;
   (*pPerms) = pcePerms;
   return 0;
}

// ************************************************************************
// assign operator
// ------------------------------------------------------------------------
MultiObjectiveOptimizer& MultiObjectiveOptimizer::operator=(const 
                                                 MultiObjectiveOptimizer &)
{
   printf("MultiObjectiveOptimizer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

