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
// Functions for the class MCMCAnalyzer
// ------------------------------------------------------------------------
// AUTHOR : CHARLES TONG
// DATE   : 2016
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "HMCMCAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PDFManager.h"
#include "Psuade.h"
#include "Sampling.h"
#include "PrintingTS.h"

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
HMCMCAnalyzer::HMCMCAnalyzer(): nInputs_(0),  means_(0), sigmas_(0), 
                                mostLikelyInput_(0), mostLikelyOutput_(0)
{
   setName("HMCMC");
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
HMCMCAnalyzer::~HMCMCAnalyzer()
{
   if(means_) delete[] means_;
   if(sigmas_) delete[] sigmas_;
   if(mostLikelyInput_) delete[] mostLikelyInput_;
   if(mostLikelyOutput_) delete[] mostLikelyOutput_;
}

// ************************************************************************
// perform MCMC analysis 
// ------------------------------------------------------------------------
double HMCMCAnalyzer::analyze()
{
   int    ii, jj, kk, status=0, nGroups, nSystems, *groupSizes;
   int    **groupInfo, *groupIDs, nParams;
   double mean, stdev, *sysData, dOne=1.0, ddata;
   psMatrix corMat, corMatH;
   FILE   *fp;

   if (means_) delete[] means_;
   if (sigmas_) delete[] sigmas_;
   if (mostLikelyInput_) delete[] mostLikelyInput_;
   means_ = NULL;
   sigmas_ = NULL;
   mostLikelyInput_ = NULL;

   status = readUserSpec(&nGroups, &nSystems, &sysData, &groupSizes, 
                         &groupInfo, &groupIDs);
   if (status < 0)
   {
      printf("ERROR in HMCMC - abort.\n");
      return -1.0;
   }

   nParams = nGroups * 2 + 1;
   mostLikelyInput_ = new double[nParams];
   means_  = new double[nParams];
   sigmas_ = new double[nParams];
   for (ii = 0; ii < nParams; ii++) means_[ii] = mostLikelyInput_[ii] = 0;
   mostLikelyOutput_ = new double[1];
   mostLikelyOutput_[0] = 0;

   int    *inputPDFsH  = new int[3];
   double *inputMeansH = new double[3];
   double *inputStdevH = new double[3];
   checkAllocate(inputStdevH, "inputStdevH in HMCMC::analyze");

   status = getHierarchicalPriors(nGroups,groupSizes,groupInfo,sysData, 
                   &mean, &stdev, inputPDFsH, inputMeansH, 
                   inputStdevH);
   corMatH.setDim(3,3); 
   for (ii = 0; ii < 3; ii++) corMatH.setEntry(ii,ii, dOne);
   printf("Population mean   = %e\n", mean);
   printf("Population st dev = %e\n", stdev);
   printf("Hyper prior 1 = N(%12.4e, %12.4e)\n",inputMeansH[0],
          inputStdevH[0]);
   printf("Hyper prior 2 = G(%12.4e, %12.4e)\n",inputMeansH[1],
          inputStdevH[1]);
   printf("Hyper prior 3 = U(%12.4e, %12.4e)\n",inputMeansH[2],
          inputStdevH[2]);

   psVector vecLBH, vecUBH, vecSamH;
   int    maxSample=5000000, maxIterations=10000000;
   double *upperH = new double[3];
   double *lowerH = new double[3];
   checkAllocate(lowerH, "lowerH in HMCMC::analyze");
   PDFManager *pdfmanH = new PDFManager();
   pdfmanH->initialize(3,inputPDFsH,inputMeansH,inputStdevH,
                       corMatH, NULL, NULL);
   lowerH[0] = inputMeansH[0] - 3 * inputStdevH[0];
   upperH[0] = inputMeansH[0] + 3 * inputStdevH[0];
   lowerH[1] = 1e-8;
   upperH[1] = inputStdevH[1] * 3;
   lowerH[2] = inputMeansH[2];
   upperH[2] = inputStdevH[2];
   vecLBH.load(3, lowerH);
   vecUBH.load(3, upperH);
   vecSamH.setLength(maxSample*3);

   //  perform inference
   int    ss, hh, index, passCnt=0, loopCnt;
   double exponent, mu, denom, maxLikelihood=0, dtemp;
   double *inferenceStore = new double[maxIterations];
   double gdelta[3], gsigma[3], *params;
   double *Xmax = new double[nParams];
   params = new double[maxIterations*nParams];
   checkAllocate(params, "params in HMCMC::analyze");

   printOutTS(PL_INFO, "Generating a large sample ... \n");
   pdfmanH->genSample(maxSample, vecSamH, vecLBH, vecUBH);

   printOutTS(PL_INFO, "HMCMC Inference begins ... \n");
   fflush(stdout);

   loopCnt = 0;
   for (hh = 0; hh < maxIterations*100; hh++)
   {
      for (ii = 0; ii < 3; ii++)
      {
         index = PSUADE_rand() % maxSample;
         gdelta[ii] = vecSamH[index*3];
      }
      for (ii = 0; ii < 3; ii++)
      {
         index = PSUADE_rand() % maxSample;
         gsigma[ii] = vecSamH[index*3+1];
      }
      index = PSUADE_rand() % maxSample;
      mu = vecSamH[index*3+2];
      //if (hh % maxIterations == 0) printOutTS(PL_INFO, ".");
      exponent = 0.0;
      for (ii = 0; ii < nSystems; ii++)
      {
         index = groupIDs[ii];
         ddata = sysData[ii]- mu - gdelta[index];
         ddata /= gsigma[index];
         exponent += 0.5 * pow(ddata, 2.0);
      }
      denom = 1.0;
      //for (ii = 0; ii < nSystems; ii++)
      //{
      //   index = groupIDs[ii];
      //   ddata = gsigma[index];
      //   denom *= ddata;
      //}
      ddata = 1.0 / denom * exp(-exponent);
      if (hh >= maxIterations)
      {
         index = hh % maxIterations;
         for (ii = 0; ii < 3; ii++) params[index*nParams+ii] = gdelta[ii];
         for (ii = 0; ii < 3; ii++) params[index*nParams+ii+3] = gsigma[ii];
         params[index*nParams+6] = mu;
         inferenceStore[index] = ddata;
      }
      if (ddata > maxLikelihood)
      {
         maxLikelihood = ddata;
         for (ii = 0; ii < 3; ii++) Xmax[ii] = gdelta[ii];
         for (ii = 0; ii < 3; ii++) Xmax[ii+3] = gsigma[ii];
         Xmax[6] = mu;
         printf("\nNew max likelihood = %e at\n", maxLikelihood);
         for (ii = 0; ii < nParams; ii++)
            printf("   X %3d = %16.8e\n",ii+1,Xmax[ii]);
      }
      if ((hh+1) % maxIterations == 0)
      {
         fp = fopen("psuade_stop","r");
         if (fp != NULL)
         {
            printf("psuade_stop file found ==> terminate.\n");
            fclose(fp);
            break;
         }
      }
      if (hh > maxIterations && (hh % maxIterations == 0))
      {
         denom = 0.0;
         for (ss = 0; ss < maxIterations; ss++) denom += inferenceStore[ss];
         passCnt = 0;
         printf("At iteration %d\n", hh);
         for (ii = 0; ii < nParams; ii++)
         {
            ddata = 0.0;
            for (ss = 0; ss < maxIterations; ss++) 
               ddata += params[ss*nParams+ii]*inferenceStore[ss]/denom;
            dtemp = 0.0;
            for (ss = 0; ss < maxIterations; ss++) 
               dtemp += pow(params[ss*nParams+ii]-ddata,2.0)*
                        inferenceStore[ss]/denom;
            means_[ii] = means_[ii] * loopCnt + ddata;
            loopCnt++;
            means_[ii] /= (double) loopCnt;
            printOutTS(PL_INFO,
                 "HMCMC: input %3d mean,stdev = %12.4e %12.4e\n",
                 ii+1,means_[ii],sqrt(dtemp));
            sigmas_[ii] = sqrt(dtemp);
            //status = checkConvergence(maxIterations,&params[ii],nParams);
            if (status == 1)
            {
               passCnt++;
               printf("MCMC input %3d converged.\n",ii+1);
            }
         }
         if (passCnt == nParams) break;
      }
   }
   genMatlabFile(maxIterations, nParams, params, inferenceStore);

   denom = 0.0;
   for (ii = 0; ii < nParams; ii++)
   {
      printOutTS(PL_INFO,"HMCMC: param %3d value at likelihood peak = %e\n",
                 ii+1, Xmax[ii]);
      mostLikelyInput_[ii] = Xmax[ii];
      for (ss = 0; ss < maxIterations; ss++) denom += inferenceStore[ss];
   }
   mostLikelyOutput_[0] = maxLikelihood / denom;

   delete [] groupIDs;
   delete [] groupSizes;
   for (ii = 0; ii < nGroups; ii++) delete [] groupInfo[ii];
   delete [] groupInfo;
   delete [] sysData;
   delete [] inputPDFsH;
   delete [] inputMeansH;
   delete [] inputStdevH;
   delete [] upperH;
   delete [] lowerH;
   delete [] inferenceStore;
   delete [] params;
   delete [] Xmax;
   delete pdfmanH;
   return 0.0;
}

// ************************************************************************
// read specfile
// ------------------------------------------------------------------------
int HMCMCAnalyzer::readUserSpec(int *nGroups_in, int *nSystems_in, 
                                double **sysData_in, int **groupSizes_in,
                                int ***groupInfo_in, int **groupIDs_in)
{
   int    ii, kk, leng, nGroups, nSystems, *groupIDs;
   int    **groupInfo, *groupSizes;
   double *sysData;
   char   pString[1000], sysFile[2000];
   FILE   *fp;
   // read group information 
   printf("HMCMC needs system group information to proceed. The file\n");
   printf("that contains this information should have this format:\n"); 
   printf("Line 1: <nGroups> <nSystems> (nSystems > nGroups)\n"); 
   printf("Next <nSystems> lines: \n");
   printf("<System number> <System Output>   (number should be in order)\n");
   printf(" ... \n");
   printf("Next Line: 1 <nSys>   (<nSys> is no. of system in Group 1)\n"); 
   printf("<sysID> <sysID> ...   (System numbers for Group 1)\n");
   printf("Next Line: 2 <nSys>   (<nSys> is no. of system in Group 2)\n"); 
   printf("<sysID> <sysID> ...   (System numbers for Group 2)\n");
   printf(" ... \n");
   printf("Example: (3 groups, 9 systems)\n");
   printf("3 9\n");
   printf("1 0.29\n");
   printf("2 0.61\n");
   printf("3 0.89\n");
   printf("4 0.59\n");
   printf("5 0.30\n");
   printf("6 0.91\n");
   printf("7 0.60\n");
   printf("8 0.31\n");
   printf("9 0.90\n");
   printf("1 3   (group 1 has 3 systems)\n");
   printf("1 5 8 (group 1 has systems 1, 5, and 8)\n");
   printf("2 3   (group 2 has 3 systems)\n");
   printf("2 4 7 (group 2 has systems 2, 4, and 7)\n");
   printf("3 3   (group 3 has 3 systems)\n");
   printf("3 6 9 (group 3 has systems 3, 6, and 9)\n");
   sprintf(pString,"Enter the file name for group information: ");
   getString(pString, sysFile);
   leng = strlen(sysFile);
   sysFile[leng-1] = '\0';
   fp = fopen(sysFile, "r");
   if (fp == NULL)
   {
      printf("ERROR: file %s not found\n", sysFile);
      return -1;
   }
   fscanf(fp,"%d %d", &nGroups, &nSystems);
   if (nSystems <= 0 || nGroups <= 0)
   {
      printf("ERROR: nSystems and nGroups should be > 0\n"); 
      return -1;
   }
   if (nSystems <= nGroups)
   {
      printf("ERROR: nSystems %d should be > nGroups %d\n", 
             nSystems, nGroups);
      return -1;
   }
   sysData = new double[nSystems];
   checkAllocate(sysData, "sysData in HMCMC::analyze");
   for (ii = 0; ii < nSystems; ii++)
   {
      fscanf(fp,"%d", &kk);
      if (kk != ii+1)
      {
         printf("ERROR: invalid system number %d detected.\n",kk);
         printf("       Expected number = %d.\n",ii+1);
         exit(1);
      }
      fscanf(fp,"%lg", &sysData[ii]);
      printf("System %5d has data = %e\n", ii+1, sysData[ii]);
   }
   groupInfo = new int*[nGroups];
   groupSizes = new int[nGroups];
   checkAllocate(groupSizes, "groupSizes in HMCMC::analyze");
   for (ii = 0; ii < nGroups; ii++)
   {
      fscanf(fp,"%d", &kk);
      if (kk != ii+1)
      {
         printf("ERROR: invalid group number %d detected.\n",kk);
         printf("       Expected number = %d.\n",ii+1);
         delete [] sysData;
         delete [] groupSizes;
         delete [] groupInfo;
         return -1;
      }
      fscanf(fp,"%d", &groupSizes[ii]);
      if (groupSizes[ii] <= 1)
      {
         printf("ERROR: group size should be > 1.\n");
         printf("       Size read = %d.\n",groupSizes[ii]);
         delete [] sysData;
         delete [] groupSizes;
         delete [] groupInfo;
         return -1;
      }
      groupInfo[ii] = new int[groupSizes[ii]];
      for (kk = 0; kk < groupSizes[ii]; kk++)
      {
         fscanf(fp,"%d", &(groupInfo[ii][kk]));
         if (groupInfo[ii][kk] < 1 || groupInfo[ii][kk] > nSystems)
         {
            printf("ERROR: invalid group member %d.\n",
                   groupInfo[ii][kk]);
            printf("       Expected: between 1 and %d\n",nSystems);
            delete [] sysData;
            delete [] groupSizes;
            delete [] groupInfo[ii];
            delete [] groupInfo;
            return -1;
         }
         groupInfo[ii][kk]--;
         printf("Group %4d has system %d\n", ii+1, groupInfo[ii][kk]+1);
      }
   }
   fclose(fp);
   kk = 0;
   for (ii = 0; ii < nGroups; ii++) kk += groupSizes[ii];
   if (kk != nSystems)
   {
      printf("ERROR: sum of all group sizes is %d\n",kk);
      printf("       Expected sum = %d\n",nSystems);
      delete [] sysData;
      delete [] groupSizes;
      for (ii = 0; ii < nGroups; ii++) delete [] groupInfo[ii];
      delete [] groupInfo;
      exit(1);
   }
   groupIDs = new int[nSystems];
   checkAllocate(groupIDs, "groupIDs in HMCMC::analyze");
   for (ii = 0; ii < nSystems; ii++) groupIDs[ii] = -1; 
   for (ii = 0; ii < nGroups; ii++) 
      for (kk = 0; kk < groupSizes[ii]; kk++) 
         groupIDs[groupInfo[ii][kk]] = ii;
   kk = 0;
   for (ii = 0; ii < nSystems; ii++) if (groupIDs[ii] != -1) kk++;
   if (kk != nSystems)
   {
      printf("ERROR: %d systems have been left out (of %d).\n",
             kk, nSystems);
      delete [] sysData;
      delete [] groupSizes;
      for (ii = 0; ii < nGroups; ii++) delete [] groupInfo[ii];
      delete [] groupInfo;
      delete [] groupIDs;
      exit(1);
   }
   (*nGroups_in) = nGroups;
   (*nSystems_in) = nSystems; 
   (*sysData_in) = sysData;
   (*groupSizes_in) = groupSizes;
   (*groupInfo_in) = groupInfo;
   (*groupIDs_in) = groupIDs;
   return 0;
}
 
// ************************************************************************
// get information 
// ------------------------------------------------------------------------
int HMCMCAnalyzer::getHierarchicalPriors(int nGroups, int *groupSizes,
                   int **groupInfo, double *sysData, double *mean,
                   double *stdev, int *inputPDFsH, double *inputMeansH, 
                   double *inputStdevH)
{
   int    ii, jj, nSystems, index;
   double pmean, pstdev, *gMeans, *gStdevs, meanBias, sdBias;
   double dmin, dmax;

   nSystems = 0;
   for (ii = 0; ii < nGroups; ii++) nSystems += groupSizes[ii];
   pmean = 0.0;
   for (ii = 0; ii < nSystems; ii++) pmean += sysData[ii];
   pmean /= (double) nSystems;
   (*mean) = pmean;
   pstdev = 0.0;
   for (ii = 0; ii < nSystems; ii++) 
      pstdev += pow(sysData[ii] - pmean, 2.0);
   pstdev = sqrt(pstdev/ (nSystems - 1));
   (*stdev) = pstdev;

   gMeans  = new double[nGroups];
   gStdevs = new double[nGroups];
   checkAllocate(gStdevs, "gStdevs in HMCMC::getHierarchicalPriors");
   for (ii = 0; ii < nGroups; ii++)
   {
      gMeans[ii] = 0;
      for (jj = 0; jj < groupSizes[ii]; jj++) 
      {
         index = groupInfo[ii][jj];  
         gMeans[ii] += sysData[index];
      }
      gMeans[ii] /= (double) groupSizes[ii];
      gStdevs[ii] = 0;
      for (jj = 0; jj < groupSizes[ii]; jj++) 
      {
         index = groupInfo[ii][jj];  
         gStdevs[ii] += pow(sysData[index] - gMeans[ii], 2.0);
      }
      gStdevs[ii] = sqrt(gStdevs[ii]/(double) (groupSizes[ii] - 1));
   }

   inputPDFsH[0] = PSUADE_PDF_NORMAL;
   meanBias = 0.0;
   for (ii = 0; ii < nGroups; ii++) meanBias += (gMeans[ii]-pmean);
   meanBias /= (double) nGroups;
   inputMeansH[0] = meanBias;
   sdBias = 0;
   for (ii = 0; ii < nGroups; ii++) 
      sdBias += pow(gMeans[ii]-pmean-meanBias,2.0);
   sdBias = sqrt(sdBias/(double) (nGroups-1));
   inputStdevH[0] = sdBias;
   inputPDFsH[1] = PSUADE_PDF_LOGNORMAL;
   meanBias = 0.0;
   for (ii = 0; ii < nGroups; ii++) meanBias += gStdevs[ii];
   meanBias /= (double) nGroups;
   sdBias = 0;
   for (ii = 0; ii < nGroups; ii++) 
      sdBias += pow(gStdevs[ii]-meanBias,2.0);
   if (sdBias == 0) sdBias = meanBias;
   else             sdBias = sqrt(sdBias/(double) (nGroups-1));
   inputStdevH[1] = sqrt(2.0 * log(meanBias/sdBias) + 1);
   inputMeansH[1] = log(meanBias) - 0.5 * inputStdevH[1] * inputStdevH[1];
   inputPDFsH[2] = PSUADE_PDF_UNIFORM;
   dmax = -1e35;
   dmin =  1e35;
   for (ii = 0; ii < nGroups; ii++)
   {
      if (gMeans[ii] > dmax) dmax = gMeans[ii];
      if (gMeans[ii] < dmin) dmin = gMeans[ii];
   }
   inputMeansH[2] = pmean - pstdev / sqrt(3.0);
   inputStdevH[2] = pmean + pstdev / sqrt(3.0);
   return 0;
}

// ************************************************************************
// check convergence 
// ------------------------------------------------------------------------
int HMCMCAnalyzer::checkConvergence(int num, double *values, int step)
{
   int    ii, jj, leng, part, nchains=5;
   double means[5], stds[5];
   double WStat, BStat, ddata, ddata2, thresh=1.02;

   leng = num / nchains * nchains;
   part = leng / nchains;
   for (ii = 0; ii < nchains; ii++)
   {
      ddata = 0.0;
      for (jj = 0; jj < part; jj++) ddata += values[(ii*part+jj)*step];
      ddata /= (double) part;
      ddata2 = 0.0;
      for (jj = 0; jj < part; jj++) 
         ddata2 += pow(values[(ii*part+jj)*step]-ddata,2.0);
      ddata2 /= (double) (part - 1.0);
      means[ii] = ddata;
      stds[ii] = sqrt(ddata2);
   } 

   WStat = 0.0;
   for (ii = 0; ii < nchains; ii++) WStat += stds[ii];
   WStat /= (double) nchains;
   if (WStat < 0) WStat = PSUADE_UNDEFINED;
   ddata = 0.0;
   for (ii = 0; ii < nchains; ii++) ddata += means[ii];
   ddata /= (double) nchains;
   BStat = 0.0;
   for (ii = 0; ii < nchains; ii++) 
      BStat += pow(means[ii]-ddata,2.0);
   printf("HMCMC check: std dev of chain means = %e (mean = %e)\n",
          sqrt(BStat/(nchains-1.0)), ddata);
   BStat = BStat / (nchains - 1.0) * part;
   ddata = (1 - 1.0/part) * WStat + BStat / part;
   ddata = ddata / WStat * (nchains + 1) / nchains - 
           (part - 1.0) / (double) (part * nchains); 
   if (ddata < 0) ddata2 = PSUADE_UNDEFINED;
   else           ddata2 = sqrt(ddata);
   printf("Convergence check: %e <? %e\n",ddata2,thresh);
   if (ddata2 < thresh) return 1;
   else                 return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
int HMCMCAnalyzer::genMatlabFile(int samSize, int nParams, double *sample,
                                 double *inferenceStore)
{
   int    ss, ii, jj, nbins=10, index, **bins;
   double denom, *lowers, *uppers, *ranges, ddata, *likelihoods;
   double **dbins;

   denom = 0.0;
   for (ss = 0; ss < samSize; ss++) denom += inferenceStore[ss];
   likelihoods = new double[samSize];
   for (ss = 0; ss < samSize; ss++) 
      likelihoods[ss] = inferenceStore[ss] / denom;
   lowers = new double[nParams];
   uppers = new double[nParams];
   ranges = new double[nParams];
   dbins = new double*[nbins];
   for (ii = 0; ii < nbins; ii++) 
   {
      dbins[ii] = new double[nParams];
      for (jj = 0; jj < nParams; jj++) dbins[ii][jj] = 0;
   }
   bins = new int*[nbins];
   for (ii = 0; ii < nbins; ii++) 
   {
      bins[ii] = new int[nParams];
      for (jj = 0; jj < nParams; jj++) bins[ii][jj] = 0;
   }
   for (ii = 0; ii < nParams; ii++) 
   {
      lowers[ii] = sample[ii];
      uppers[ii] = sample[ii];
      for (ss = 1; ss < samSize; ss++)
      {
         if (sample[ss*nParams+ii] < lowers[ii])
            lowers[ii] = sample[ss*nParams+ii];
         if (sample[ss*nParams+ii] > uppers[ii])
            uppers[ii] = sample[ss*nParams+ii];
      }
      ranges[ii] = uppers[ii] - lowers[ii];
   }
   for (ss = 0; ss < samSize; ss++)
   {
      for (ii = 0; ii < nParams; ii++) 
      {
         ddata = (sample[ss*nParams+ii] - lowers[ii]) / ranges[ii];
         index = (int) (ddata * nbins);
         if (index >= nbins) index = nbins - 1;
         dbins[index][ii] += likelihoods[ss];
      }
   }
   double Ymax;
   for (jj = 0; jj < nParams; jj++) 
   {
      Ymax = 0.0;
      for (ii = 0; ii < nbins; ii++)
         if (dbins[ii][jj] > Ymax) Ymax = dbins[ii][jj];
      for (ii = 0; ii < nbins; ii++)
         dbins[ii][jj] = dbins[ii][jj] / Ymax * samSize;
      for (ii = 0; ii < nbins; ii++)
         bins[ii][jj] = (int) dbins[ii][jj];
   }
   char cfname[1000], charString[1000];
   FILE *fp;
   strcpy(cfname, "matlabhmcmc.m");
   fp = fopen(cfname, "w");
   if (fp == NULL)
   {
      printOutTS(PL_ERROR, "ERROR: cannot open %s file.\n", cfname);
      return 0;
   }
   sprintf(charString,"This file shows posteriors plots");
   fwriteComment(fp, charString);
   fwritePlotCLF(fp);
   fprintf(fp, "L = [\n");
   for (ii = 0; ii < nParams; ii++) fprintf(fp, "%e ",lowers[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "U = [\n");
   for (ii = 0; ii < nParams; ii++) fprintf(fp, "%e ",uppers[ii]);
   fprintf(fp, "];\n");
   fprintf(fp, "X = zeros(%d,%d);\n", nParams, nbins);
   fprintf(fp, "D = zeros(%d,%d);\n", nParams, nbins);
   int kk, kk2, sumBins;
   for (kk = 0; kk < nParams; kk++)
   {
      for (kk2 = 0; kk2 < nParams; kk2++)
      {
         if (kk == kk2)
         {
            fprintf(fp, "X(%d,:) = [\n", kk+1);
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%e ", ranges[kk]/nbins*(jj+0.5)+lowers[kk]);
            fprintf(fp, "];\n");
            sumBins = 0;
            for (jj = 0; jj < nbins; jj++) sumBins += bins[jj][kk];
            if (sumBins == 0) sumBins = 1;
            fprintf(fp, "D(%d,:) = [\n", kk+1);
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%e ", (double) bins[jj][kk]/(double) sumBins);
            fprintf(fp, "];\n");
         }
      }
   }
   fprintf(fp,"nInps = %d;\n", nParams);
   fprintf(fp,"for ii = 1 : nInps\n");
   fprintf(fp,"   subplot(nInps,nInps,(ii-1)*nInps+ii)\n");
   fprintf(fp,"   n = length(D(ii,:));\n");
   fprintf(fp,"   DN = D(ii,:);\n");
   fprintf(fp,"   bar(X(ii,:), DN, 1.0);\n");
   fprintf(fp,"   xmin = min(X(ii,:));\n");
   fprintf(fp,"   xmax = max(X(ii,:));\n");
   fprintf(fp,"   xwid = xmax - xmin;\n");
   fprintf(fp,"   xmin = xmin - 0.5 * xwid / %d;\n", nbins);
   fprintf(fp,"   xmax = xmax + 0.5 * xwid / %d;\n", nbins);
   fprintf(fp,"   ymax = max(DN);\n");
   fprintf(fp,"   axis([xmin xmax 0 ymax])\n");
   fprintf(fp,"   set(gca,'linewidth',2)\n");
   fprintf(fp,"   set(gca,'fontweight','bold')\n");
   fprintf(fp,"   set(gca,'fontsize',12)\n");
   fprintf(fp,
      "   ylabel('Probabilities','FontWeight','bold','FontSize',12)\n");
   fprintf(fp,"   grid on\n");
   fprintf(fp,"   box on\n");
   fprintf(fp,"end;\n");
   fclose(fp);
   delete [] lowers;
   delete [] uppers;
   delete [] ranges;
   for (ii = 0; ii < nbins; ii++) delete [] dbins[ii];
   for (ii = 0; ii < nbins; ii++) delete [] bins[ii];
   delete [] dbins;
   delete [] bins;
   delete [] likelihoods;
   return 0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
HMCMCAnalyzer& HMCMCAnalyzer::operator=(const HMCMCAnalyzer &)
{
   printOutTS(PL_ERROR,
        "HMCMCAnalyzer operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

// ************************************************************************
// functions for getting results
// ------------------------------------------------------------------------
int HMCMCAnalyzer::get_nInputs()
{
   return nInputs_;
}
int HMCMCAnalyzer::get_nOutputs()
{
   return 1;
}
double *HMCMCAnalyzer::get_means()
{
   int    ii;
   double *retVal = NULL;
   if (means_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in HMCMC::get_nInputs");
      for (ii = 0; ii < nInputs_; ii++) retVal[ii] = means_[ii];
   }
   return retVal;
}
double *HMCMCAnalyzer::get_mostLikelyInput()
{
   int    ii;
   double *retVal = NULL;
   if (mostLikelyInput_)
   {
      retVal = new double[nInputs_];
      checkAllocate(retVal, "retVal in HMCMC::get_mostLikelyInput");
      for (ii = 0; ii < nInputs_; ii++) retVal[ii] = mostLikelyInput_[ii];
   }
   return retVal;
}
double *HMCMCAnalyzer::get_mostLikelyOutput()
{
   int    ii;
   double *retVal = NULL;
   if (mostLikelyOutput_)
   {
      retVal = new double[1];
      checkAllocate(retVal, "retVal in HMCMC::get_mostLikelyOutput");
      retVal[0] = mostLikelyOutput_[0];
   }
   return retVal;
}

