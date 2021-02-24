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
// Functions for the class RBF
// AUTHOR : CHARLES TONG
// DATE   : 2014
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "RBF.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"
#include "PrintingTS.h"

extern "C" {
#if 0
   void dgetrf_(int *, int *, double *, int *, int *, int *);
   void dgetrs_(char *,int *,int*,double*,int*,int*,double*,int*,int*);
#endif
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
}

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

#define PS_RBF1

// ************************************************************************
// Constructor for object class RBF
// ------------------------------------------------------------------------
RBF::RBF(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   char pString[501];

   faID_ = PSUADE_RS_RBF;
   type_ = 0;

   XNormalized_ = NULL;
   YNormalized_ = NULL;
   iRanges_ = NULL;
   regCoeffs_ = NULL;
   svdThresh_ = 1e-15;
   gaussScale_ = 1;

   // display banner and additonal information
   if (psInteractive_ == 1)
   {
      printAsterisks(PL_INFO, 0);
      printOutTS(PL_INFO,
           "*           Radial Basis Function (RBF) Analysis\n");
      printOutTS(PL_INFO,"* Set printlevel to 1-4 to see RBF details.\n");
      printOutTS(PL_INFO,"* Default kernel    = multi-quadratic \n");
      printOutTS(PL_INFO,
           "* Default threshold = 1.0e-15 (for SVD truncation)\n");
      printOutTS(PL_INFO,"* Turn on rs_expert mode to make changes.\n");
      printEquals(PL_INFO, 0);
   }
   
   if (psRSExpertMode_ == 1 && psInteractive_ == 1)
   {
      printf("In the following you have the option to select the kernel. \n");
      printf("0. multi-quadratic\n");
      printf("1. inverse multi-quadratic\n");
      printf("2. Gaussian\n");
      printf("3. thin plate spline\n");
      sprintf(pString,"Enter your choice (0 - 3) : ");
      type_ = getInt(0, 3, pString);
      if (type_ == 2)
      {
         sprintf(pString,
            "Enter scaling factor for Gaussian kernel (default=1) : ");
         gaussScale_ = getDouble(pString);
      }
      printOutTS(PL_INFO,
           "The RBF matrix to be constructed may be near-singular.\n");
      printOutTS(PL_INFO,
           "Currently, singular values < max(svd)*1e-15 are truncated.\n");
      printOutTS(PL_INFO,
           "You have the option to change this threshold (1e-15).\n");
      printOutTS(PL_INFO,
           "NOTE: truncating singular values may lead to erroneous results.\n");
      sprintf(pString, "Enter new threshold for SVD (> 0 but << 1) : ");
      svdThresh_ = getDouble(pString);
   }
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
RBF::~RBF()
{
   if (XNormalized_ != NULL) delete [] XNormalized_;
   if (YNormalized_ != NULL) delete [] YNormalized_;
   if (iRanges_     != NULL) delete [] iRanges_;
   if (regCoeffs_   != NULL) delete [] regCoeffs_;
}

// ************************************************************************
// initialize 
// ------------------------------------------------------------------------
int RBF::initialize(double *X, double *Y)
{
   int    *iArray, ii, kk, ss, ss2, nSubSamples, count, nSamp1;
   double range, *Dmat, ddata;
   char   pString[501];
   FILE   *fp;

   if (lowerBounds_ == NULL)
   {
      printOutTS(PL_ERROR,
           "RBF initialize ERROR: sample bounds not set yet.\n");
      return -1;
   }
   if (XNormalized_ != NULL) delete [] XNormalized_;
   if (YNormalized_ != NULL) delete [] YNormalized_;
   if (iRanges_     != NULL) delete [] iRanges_;
   if (regCoeffs_   != NULL) delete [] regCoeffs_;
   XNormalized_ = new double[nSamples_*nInputs_];
   checkAllocate(XNormalized_, "XNormalized_ in RBF::initialize");
   for (ii = 0; ii < nInputs_; ii++)
   {
      range = 1.0 / (upperBounds_[ii] - lowerBounds_[ii]);
      for (ss = 0; ss < nSamples_; ss++)
         XNormalized_[ss*nInputs_+ii] = 
            (X[ss*nInputs_+ii] - lowerBounds_[ii]) * range;
   }
   YNormalized_ = new double[nSamples_];
   checkAllocate(YNormalized_, "YNormalized_ in RBF::initialize");
   initOutputScaling(Y, YNormalized_);
   for (ii = 0; ii < nSamples_; ii++) YNormalized_[ii] = Y[ii] - YMean_;
   iRanges_ = new double[nInputs_];
   checkAllocate(iRanges_, "iRanges_ in RBF::initialize");
   for (ii = 0; ii < nInputs_; ii++)
      iRanges_[ii] = 1.0 / (upperBounds_[ii] - lowerBounds_[ii]);
 
#ifdef PS_RBF1
   nSamp1 = nSamples_ + 1;
#else
   nSamp1 = nSamples_;
#endif
   Dmat = new double[nSamp1*nSamp1];
   checkAllocate(Dmat, "Dmat in RBF::initialize");
   switch(type_) 
   {
      case 0: 
         if (outputLevel_ > 0) 
            printOutTS(PL_INFO,"Kernel = multi-quadratic\n");
         for (ss = 0; ss < nSamples_; ss++)
         {
            Dmat[ss*nSamp1+ss] = 1.0; 
            for (ss2 = ss+1; ss2 < nSamples_; ss2++)
            {
               ddata = 0.0;
               for (ii = 0; ii < nInputs_; ii++)
                  ddata += pow((XNormalized_[ss*nInputs_+ii]-
                                XNormalized_[ss2*nInputs_+ii]),2.0);
               Dmat[ss*nSamp1+ss2] = 
                    Dmat[ss2*nSamp1+ss] = sqrt(ddata+1.0);
            }
         }
         break;

      case 1: 
         if (outputLevel_ > 0) 
            printOutTS(PL_INFO,"Kernel = inverse multi-quadratic\n");
         for (ss = 0; ss < nSamples_; ss++)
         {
            Dmat[ss*nSamp1+ss] = 1.0; 
            for (ss2 = ss+1; ss2 < nSamples_; ss2++)
            {
               ddata = 0.0;
               for (ii = 0; ii < nInputs_; ii++)
                  ddata += pow((XNormalized_[ss*nInputs_+ii]-
                                XNormalized_[ss2*nInputs_+ii]),2.0);
               Dmat[ss*nSamp1+ss2] = 
                    Dmat[ss2*nSamp1+ss] = 1.0/sqrt(ddata+1.0);
            }
         }
         break;

      case 2: 
         if (outputLevel_ > 0) 
            printOutTS(PL_INFO,"Kernel = Gaussian\n");
         for (ss = 0; ss < nSamples_; ss++)
         {
            Dmat[ss*nSamp1+ss] = 1.0; 
            for (ss2 = ss+1; ss2 < nSamples_; ss2++)
            {
               ddata = 0.0;
               for (ii = 0; ii < nInputs_; ii++)
                  ddata += pow((XNormalized_[ss*nInputs_+ii]-
                                XNormalized_[ss2*nInputs_+ii]),2.0);
               Dmat[ss*nSamp1+ss2] = 
                    Dmat[ss2*nSamp1+ss] = exp(-gaussScale_*ddata/2.0);
            }
         }
         break;

      case 3: 
         if (outputLevel_ > 0) 
            printOutTS(PL_INFO,"Kernel = thin plate spline\n");
         for (ss = 0; ss < nSamples_; ss++)
         {
            Dmat[ss*nSamp1+ss] = 0.0; 
            for (ss2 = ss+1; ss2 < nSamples_; ss2++)
            {
               ddata = 0.0;
               for (ii = 0; ii < nInputs_; ii++)
                  ddata += pow((XNormalized_[ss*nInputs_+ii]-
                                XNormalized_[ss2*nInputs_+ii]),2.0);
               Dmat[ss*nSamp1+ss2] = 
                    Dmat[ss2*nSamp1+ss] = (ddata+1.0)*log(sqrt(ddata+1.0));
            }
         }
         break;
   }
#ifdef PS_RBF1
   for (ss = 0; ss < nSamples_; ss++)
      Dmat[ss*nSamp1+nSamples_] = Dmat[nSamples_*nSamp1+ss] = 1.0; 
   Dmat[nSamples_*nSamp1+nSamples_] = 0.0;
#endif

#if 0
   int    info, iOne=1;
   char   trans='N';
   int    *ipiv = new int[nSamp1+1];
   dgetrf_(&nSamp1,&nSamp1,Dmat,&nSamp1,ipiv, &info);
   if (info != 0) printf("RBF WARNING: dgetrf returns error %d.\n",info);
   regCoeffs_ = new double[nSamp1];
   for (ii = 0; ii < nSamples_; ii++) regCoeffs_[ii] = YNormalized_[ii];
#ifdef PS_RBF1
   regCoeffs_[nSamples_] = 0.0;
#endif
   dgetrs_(&trans,&nSamp1,&iOne,Dmat,&nSamp1,ipiv,regCoeffs_,&nSamp1,&info);
   ddata = 0.0;
   for (ii = 0; ii < nSamples_; ii++) ddata += regCoeffs_[ii];
   if (info != 0) printf("RBF WARNING: dgetrs returns error %d.\n",info);
   printf("RBF check: sum of coefficients = %e\n", ddata);
   delete [] Dmat;
   delete [] ipiv;
#else
   int    info, cnt=0;
   int    wlen = 5 * nSamp1;
   char   jobu = 'A', jobvt = 'A';
   double *SS = new double[nSamp1];
   double *UU = new double[nSamp1*nSamp1];
   double *VV = new double[nSamp1*nSamp1];
   double *WW = new double[wlen];
   checkAllocate(WW, "WW in RBF::initialize");
   dgesvd_(&jobu,&jobvt,&nSamp1,&nSamp1,Dmat,&nSamp1,SS,UU,&nSamp1,VV,
           &nSamp1,WW, &wlen,&info);
   if (info != 0) 
   {
      printOutTS(PL_WARN,"RBF ERROR: dgesvd returns error %d.\n",info);
      delete [] SS;
      delete [] UU;
      delete [] VV;
      delete [] WW;
      delete [] Dmat;
      return -1;
   }
   regCoeffs_ = new double[nSamp1];
   checkAllocate(regCoeffs_, "regCoeffs_ in RBF::initialize");
   for (ii = 0; ii < nSamples_; ii++) regCoeffs_[ii] = YNormalized_[ii];
#ifdef PS_RBF1
   regCoeffs_[nSamples_] = 0.0;
#endif
   for (ss = 1; ss < nSamp1; ss++)
   {
      if (SS[ss]/SS[0] < svdThresh_)
      {
         SS[ss] = 0;
         cnt++;
      }
   }
   if (cnt > 0 && psInteractive_ == 1 && outputLevel_ > 0) 
   {
      printOutTS(PL_WARN,
           "WARNING: RBF matrix is near-singular. Small singular values\n");
      printOutTS(PL_WARN,
           "         (%d out of %d) are truncated.\n",cnt,nSamp1);
      printOutTS(PL_WARN,"         Approximation may be inaccurate.\n");
   }
   for (ss = 0; ss < nSamp1; ss++)
   {
      WW[ss] = 0.0;
      for (ss2 = 0; ss2 < nSamp1; ss2++)
         WW[ss] += UU[ss*nSamp1+ss2] * regCoeffs_[ss2];
   }
   for (ss = 0; ss < nSamp1; ss++) 
   {
      if (SS[ss] != 0) WW[ss] /= SS[ss];
      else             WW[ss] = 0;
   }
   for (ss = 0; ss < nSamp1; ss++)
   {
      regCoeffs_[ss] = 0.0;
      for (ss2 = 0; ss2 < nSamp1; ss2++) 
         regCoeffs_[ss] += VV[ss*nSamp1+ss2] * WW[ss2];
   }
   delete [] SS;
   delete [] UU;
   delete [] VV;
   delete [] WW;
   delete [] Dmat;
#endif

   if (psRSCodeGen_ == 0) return 0;
   fp = fopen("psuade_rs.info", "w");
   if (fp != NULL)
   {
      fprintf(fp,"/* *********************************************/\n");
      fprintf(fp,"/* RBF interpolator from PSUADE.       */\n");
      fprintf(fp,"/* ============================================*/\n");
      fprintf(fp,"/* This file contains information for interpolation\n");
      fprintf(fp,"   using response surface. Follow the steps below:\n");
      fprintf(fp,"   1. move this file to *.c file (e.g. main.c)\n");
      fprintf(fp,"   2. Compile main.c (cc -o main main.c -lm) \n");
      fprintf(fp,"   3. run: main input output\n");
      fprintf(fp,"          where input has the number of inputs and\n");
      fprintf(fp,"          the input values\n");
      fprintf(fp,"*/\n");
      fprintf(fp,"/* ==========================================*/\n");
      fprintf(fp,"int nSamples = %d;\n",nSamples_);
      fprintf(fp,"int nInps = %d;\n",nInputs_);
      fprintf(fp,"static double\n");
      fprintf(fp,"LBs[%d] = \n", nInputs_);
      fprintf(fp,"{\n");
      for (ii = 0; ii < nInputs_; ii++)
         fprintf(fp,"  %24.16e ,\n", lowerBounds_[ii]);
      fprintf(fp,"};\n");
      fprintf(fp,"static double\n");
      fprintf(fp,"UBs[%d] = \n", nInputs_);
      fprintf(fp,"{\n");
      for (ii = 0; ii < nInputs_; ii++)
         fprintf(fp,"  %24.16e ,\n", upperBounds_[ii]);
      fprintf(fp,"};\n");
      fprintf(fp,"static double\n");
      fprintf(fp,"Coefs[%d] = \n", nSamp1);
      fprintf(fp,"{\n");
      for (ss = 0; ss < nSamples_; ss++)
         fprintf(fp,"  %24.16e ,\n", regCoeffs_[ss]);
#ifdef PS_RBF1
      fprintf(fp,"  %24.16e };\n", regCoeffs_[nSamples_]);
#endif
      fprintf(fp,"static double\n");
      fprintf(fp,"XN[%d][%d] = \n", nSamples_, nInputs_);
      fprintf(fp,"{\n");
      for (ss = 0; ss < nSamples_; ss++)
      {
         fprintf(fp,"   { ");
         for (ii = 0; ii < nInputs_-1; ii++)
            fprintf(fp,"  %24.16e ,", XNormalized_[ss*nInputs_+ii]);
         fprintf(fp,"  %24.16e },\n", XNormalized_[ss*nInputs_+nInputs_-1]);
      }
      fprintf(fp,"};\n");
      fprintf(fp,"double YMean = %e;\n",YMean_);
      fprintf(fp,"/* *********************************************/\n");
      fprintf(fp,"/* RBF interpolator from PSUADE.       */\n");
      fprintf(fp,"/* ==========================================*/\n");
      fprintf(fp,"#include <math.h>\n");
      fprintf(fp,"#include <stdlib.h>\n");
      fprintf(fp,"#include <stdio.h>\n");
      fprintf(fp,"int interpolate(int,double*,double*);\n");
      fprintf(fp,"main(int argc, char **argv) {\n");
      fprintf(fp,"  int    i, iOne=1, nInps;\n");
      fprintf(fp,"  double X[%d], Y, S;\n",nInputs_);
      fprintf(fp,"  FILE   *fIn=NULL, *fOut=NULL;\n");
      fprintf(fp,"  if (argc != 3) {\n");
      fprintf(fp,"     printf(\"ERROR: not enough argument.\\n\");\n");
      fprintf(fp,"     exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  fIn = fopen(argv[1], \"r\");\n");
      fprintf(fp,"  if (fIn == NULL) {\n");
      fprintf(fp,"     printf(\"ERROR: cannot open input file.\\n\");\n");
      fprintf(fp,"     exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  fscanf(fIn, \"%%d\", &nInps);\n");
      fprintf(fp,"  if (nInps != %d) {\n", nInputs_);
      fprintf(fp,"    printf(\"ERROR - wrong nInputs.\\n\");\n");
      fprintf(fp,"    exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  for (i=0; i<nInps; i++) fscanf(fIn, \"%%lg\", &X[i]);\n");
      fprintf(fp,"  fclose(fIn);\n");
      fprintf(fp,"  interpolate(iOne, X, &Y);\n");
      fprintf(fp,"  fOut = fopen(argv[2], \"w\");\n");
      fprintf(fp,"  if (fOut == NULL) {\n");
      fprintf(fp,"     printf(\"ERROR: cannot open output file.\\n\");\n");
      fprintf(fp,"     exit(1);\n");
      fprintf(fp,"  }\n");
      fprintf(fp,"  fprintf(fOut,\" %%e\\n\", Y);\n");
      fprintf(fp,"  fclose(fOut);\n");
      fprintf(fp,"}\n\n");
      fprintf(fp,"/* *************************************/\n");
      fprintf(fp,"/*  interpolation function             */\n");
      fprintf(fp,"/* X[0], X[1],   .. X[m-1]   - first point\n");
      fprintf(fp," * X[m], X[m+1], .. X[2*m-1] - second point\n");
      fprintf(fp," * ... */\n");
      fprintf(fp,"/* ====================================*/\n");
      fprintf(fp,"int interpolate(int npts,double *X,double *Y){\n");
      fprintf(fp,"  int    ss, ss2, ii;\n");
      fprintf(fp,"  double dist, dd;\n");
      fprintf(fp,"  for (ss2 = 0; ss2 < npts; ss2++) {\n");
      fprintf(fp,"    Y[ss2] = 0;\n");
      fprintf(fp,"    for (ss = 0; ss < nSamples; ss++) {\n");
      fprintf(fp,"      dist = 0.0;\n");
      fprintf(fp,"      for (ii = 0; ii < nInps; ii++) {\n");
      fprintf(fp,"        dd = X[ii+ss2*nInps];\n");
      fprintf(fp,"        dd = (dd-LBs[ii])/(UBs[ii]-LBs[ii]);\n");
      fprintf(fp,"        dd = dd - XN[ss][ii];\n");
      fprintf(fp,"        dist += dd * dd;\n");
      fprintf(fp,"      }\n");
      switch (type_) 
      {
         case 0: 
            fprintf(fp,"      dist = sqrt(dist+1);\n");
            break;
         case 1: 
            fprintf(fp,"      dist = 1.0/sqrt(dist+1);\n");
            break;
         case 2: 
            fprintf(fp,"      dist = exp(-0.5*dist*%e);\n",gaussScale_);
            break;
         case 3: 
            fprintf(fp,"      dist = (dist+1)*log(sqrt(dist+1));\n");
            break;
      }
      fprintf(fp,"      Y[ss2] += dist * Coefs[ss];\n");
      fprintf(fp,"    }\n");
#ifdef PS_RBF1
      fprintf(fp,"    Y[ss2] = Y[ss2] + YMean + Coefs[nSamples];\n");
#else
      fprintf(fp,"    Y[ss2] = Y[ss2] + YMean;\n");
#endif
      fprintf(fp,"  }\n");
      fprintf(fp,"  return 0;\n");
      fprintf(fp,"}\n");
      fclose(fp);
   }
   fp = fopen("psuade_rs.py", "w");
   if (fp != NULL)
   {
      fwriteRSPythonHeader(fp);
      fprintf(fp,"#==================================================\n");
      fprintf(fp,"# RBF regression interpolation\n");
      fprintf(fp,"#==================================================\n");
      fwriteRSPythonCommon(fp);
      fprintf(fp,"nSamples = %d;\n",nSamples_);
      fprintf(fp,"nInps = %d;\n",nInputs_);
      fprintf(fp,"LBs = [\n");
      for (ii = 0; ii < nInputs_; ii++)
         fprintf(fp," %24.16e ,\n", lowerBounds_[ii]);
      fprintf(fp,"]\n");
      fprintf(fp,"UBs = [\n");
      for (ii = 0; ii < nInputs_; ii++)
         fprintf(fp," %24.16e ,\n", upperBounds_[ii]);
      fprintf(fp,"]\n");
      fprintf(fp,"Coefs = [\n");
      for (ss = 0; ss <= nSamples_; ss++)
      {
         fprintf(fp,"  %24.16e , \n", regCoeffs_[ss]);
      }
      fprintf(fp,"]\n");
      fprintf(fp,"XN = [\n");
      for (ss = 0; ss < nSamples_; ss++)
      {
         fprintf(fp,"   [ ");
         for (ii = 0; ii < nInputs_-1; ii++)
            fprintf(fp,"  %24.16e ,", XNormalized_[ss*nInputs_+ii]);
         fprintf(fp,"  %24.16e ],\n", XNormalized_[ss*nInputs_+nInputs_-1]);
      }
      fprintf(fp,"]\n");
      fprintf(fp,"YMean = %e\n", YMean_);
      fprintf(fp,"###################################################\n");
      fprintf(fp,"# interpolation function  \n");
      fprintf(fp,"# X[0], X[1],   .. X[m-1]   - first point\n");
      fprintf(fp,"# X[m], X[m+1], .. X[2*m-1] - second point\n");
      fprintf(fp,"# ... \n");
      fprintf(fp,"#==================================================\n");
      fprintf(fp,"def interpolate(XX): \n");
      fprintf(fp,"  npts = int(len(XX) / nInps + 1.0e-8)\n");
      fprintf(fp,"  Ys = (2 * npts) * [0.0]\n");
      fprintf(fp,"  X  = nInps * [0.0]\n");
      fprintf(fp,"  for nn in range(npts) : \n");
      fprintf(fp,"    for ii in range(nInps) : \n");
      fprintf(fp,"      X[ii] = XX[nn*nInps+ii]\n");
      fprintf(fp,"    for ss in range(nSamples) : \n");
      fprintf(fp,"      dist = 0.0\n");
      fprintf(fp,"      for ii in range(nInps) : \n");
      fprintf(fp,"        dd = (X[ii]-LBs[ii])/(UBs[ii]-LBs[ii])\n");
      fprintf(fp,"        dd = dd - XN[ss][ii]\n");
      fprintf(fp,"        dist += dd * dd\n");
      switch (type_) 
      {
         case 0: 
            fprintf(fp,"      dist = math.sqrt(dist+1)\n");
            break;
         case 1: 
            fprintf(fp,"      dist = 1.0/math.sqrt(dist+1)\n");
            break;
         case 2: 
            fprintf(fp,"      dist = math.exp(-0.5*dist*%e)\n",gaussScale_);
            break;
         case 3: 
            fprintf(fp,"      dist = (dist+1)*math.log(math.sqrt(dist+1))\n");
            break;
      }
      fprintf(fp,"      Ys[nn] = Ys[nn] + Coefs[ss] * dist\n");
#ifdef PS_RBF1
      fprintf(fp,"    Ys[nn] = Ys[nn] + YMean + Coefs[nSamples]\n");
#else
      fprintf(fp,"    Ys[nn] = Ys[nn] + YMean\n");
#endif
      fprintf(fp,"  return Ys\n");
      fprintf(fp,"###################################################\n");
      fprintf(fp,"# main program\n");
      fprintf(fp,"#==================================================\n");
      fprintf(fp,"infileName  = sys.argv[1]\n");
      fprintf(fp,"outfileName = sys.argv[2]\n");
      fprintf(fp,"inputs = getInputData(infileName)\n");
      fprintf(fp,"outputs = interpolate(inputs)\n");
      fprintf(fp,"genOutputFile(outfileName, outputs)\n");
      fprintf(fp,"###################################################\n");
      printf("FILE psuade_rs.py contains the final RBF functional form.\n");
      fclose(fp);
   }
   return 0;
}
  
// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int RBF::genNDGridData(double *X,double *Y,int *N2,double **X2,double **Y2)
{
   int totPts;

   initialize(X,Y);

   if ((*N2) == -999) return 0;
  
   genNDGrid(N2, X2);
   if ((*N2) == 0) return 0;
   totPts = (*N2);

   (*Y2) = new double[totPts];
   checkAllocate(*Y2, "Y2 in RBF::genNDGridData");
   evaluatePoint(totPts, *X2, *Y2);

   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int RBF::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
                       int *N2, double **X2, double **Y2)
{
   int    ii, ss, totPts;
   double *XT, *XX, *YY, HX;

   initialize(X,Y);

   totPts = nPtsPerDim_;
   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[totPts];
   (*Y2) = new double[totPts];
   (*N2) = totPts;
   XX = (*X2);
   YY = (*Y2);

   XT = new double[totPts*nInputs_];
   checkAllocate(XT, "XT in RBF::gen1DGridData");
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ss = 0; ss < totPts; ss++) 
   {
      XT[ss*nInputs_+ind1]  = HX * ss + lowerBounds_[ind1];
      XX[ss] = HX * ss + lowerBounds_[ind1];
      YY[ss] = 0.0;
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int RBF::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                       double *settings, int *N2, double **X2, double **Y2)
{
   int    ii, ss, jj, index, totPts;
   double *XT, *XX, *YY, *HX;
 
   initialize(X,Y);

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[2*totPts];
   (*Y2) = new double[totPts];
   (*N2) = totPts;
   XX = (*X2);
   YY = (*Y2);

   XT = new double[totPts*nInputs_];
   checkAllocate(XT, "XT in RBF::gen2DGridData");
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         index = ii * nPtsPerDim_ + jj;
         XT[index*nInputs_+ind1] = HX[0] * ii + lowerBounds_[ind1];
         XT[index*nInputs_+ind2] = HX[1] * jj + lowerBounds_[ind2];
         XX[index*2]   = HX[0] * ii + lowerBounds_[ind1];
         XX[index*2+1] = HX[1] * jj + lowerBounds_[ind2];
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int RBF::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3, 
                       double *settings, int *N2, double **X2, double **Y2)
{
   int    ii, ss, jj, ll, index, totPts;
   double *XT, *XX, *YY, *HX;

   initialize(X,Y);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[3*totPts];
   (*Y2) = new double[totPts];
   (*N2) = totPts;
   XX = (*X2);
   YY = (*Y2);

   XT = new double[totPts*nInputs_];
   checkAllocate(XT, "XT in RBF::gen3DGridData");
   for (ss = 0; ss < totPts; ss++)
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii];

   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            index = ii * nPtsPerDim_ * nPtsPerDim_ + jj * nPtsPerDim_ + ll;
            XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
            XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
            XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
            XX[index*3]   = HX[0] * ii + lowerBounds_[ind1];
            XX[index*3+1] = HX[1] * jj + lowerBounds_[ind2];
            XX[index*3+2] = HX[2] * ll + lowerBounds_[ind3];
         }
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int RBF::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3, 
                       int ind4, double *settings, int *N2, double **X2, 
                       double **Y2)
{
   int    ii, ss, jj, ll, mm, index, totPts;
   double *XT, *XX, *YY, *HX;

   initialize(X,Y);

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   (*X2) = new double[4*totPts];
   (*Y2) = new double[totPts];
   (*N2) = totPts;
   XX = (*X2);
   YY = (*Y2);

   XT = new double[totPts*nInputs_];
   checkAllocate(XT, "XT in RBF::gen4DGridData");
   for (ss = 0; ss < totPts; ss++) 
      for (ii = 0; ii < nInputs_; ii++) XT[ss*nInputs_+ii] = settings[ii]; 
    
   for (ii = 0; ii < nPtsPerDim_; ii++) 
   {
      for (jj = 0; jj < nPtsPerDim_; jj++)
      {
         for (ll = 0; ll < nPtsPerDim_; ll++)
         {
            for (mm = 0; mm < nPtsPerDim_; mm++)
            {
               index = ii*nPtsPerDim_*nPtsPerDim_ * nPtsPerDim_ +
                       jj*nPtsPerDim_*nPtsPerDim_ + ll*nPtsPerDim_ + mm;
               XT[index*nInputs_+ind1]  = HX[0] * ii + lowerBounds_[ind1];
               XT[index*nInputs_+ind2]  = HX[1] * jj + lowerBounds_[ind2];
               XT[index*nInputs_+ind3]  = HX[2] * ll + lowerBounds_[ind3];
               XT[index*nInputs_+ind4]  = HX[3] * mm + lowerBounds_[ind4];
               XX[index*4]   = HX[0] * ii + lowerBounds_[ind1];
               XX[index*4+1] = HX[1] * jj + lowerBounds_[ind2];
               XX[index*4+2] = HX[2] * ll + lowerBounds_[ind3];
               XX[index*4+3] = HX[3] * mm + lowerBounds_[ind4];
            }
         }
      }
   }

   evaluatePoint(totPts, XT, YY);

   delete [] XT;
   delete [] HX;
   return 0;
}

// ************************************************************************
// Evaluate a point
// ------------------------------------------------------------------------
double RBF::evaluatePoint(double *X)
{
   int    iOne=1;
   double Y;
   evaluatePoint(1, X, &Y);
   return Y;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double RBF::evaluatePoint(int nPts, double *X, double *Y)
{
   int    ss, ss2, ii;
   double dist, Yt, ddata;

   for (ss = 0; ss < nPts; ss++) 
   {
      Yt = 0.0;
      for (ss2 = 0; ss2 < nSamples_; ss2++) 
      {
         dist = 0.0;
         for (ii = 0; ii < nInputs_; ii++) 
         {
            ddata = X[ss*nInputs_+ii];
            ddata = (ddata - lowerBounds_[ii]) * iRanges_[ii];
            ddata -= XNormalized_[ss2*nInputs_+ii];
            dist += ddata * ddata;
         }
         switch (type_)
         {
            case 0: dist = sqrt(dist + 1.0); break;
            case 1: dist = 1.0/sqrt(dist + 1.0); break;
            case 2: dist = exp(-0.5*dist*gaussScale_); break;
            case 3: dist = (dist+1)*log(sqrt(dist+1)); break;
         }
         Yt += dist * regCoeffs_[ss2];;
      }
#ifdef PS_RBF1
      Y[ss] = Yt + YMean_ + regCoeffs_[nSamples_];
#else
      Y[ss] = Yt + YMean_;
#endif
   }
   return 0.0;
}

// ************************************************************************
// Evaluate a given point with standard deviation
// ------------------------------------------------------------------------
double RBF::evaluatePointFuzzy(double *X, double &std)
{
   int    iOne=1;
   double Y=0.0;
   evaluatePoint(iOne, X, &Y);
   std = 0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points with standard deviations
// ------------------------------------------------------------------------
double RBF::evaluatePointFuzzy(int npts, double *X, double *Y, double *Ystd)
{
   evaluatePoint(npts, X, Y);
   for (int ss = 0; ss < npts; ss++) Ystd[ss] = 0.0;
   return 0.0;
}

