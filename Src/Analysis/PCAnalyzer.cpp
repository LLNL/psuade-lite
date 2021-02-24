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
// Functions for the class PCAnalyzer (principal component analysis)  
// AUTHOR : CHARLES TONG
// DATE   : 2005
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "PCAnalyzer.h"
#include "PsuadeUtil.h"
#include "sysdef.h"
#include "Psuade.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

extern "C" {
   void dgels_(char *, int *, int *, int *, double *A, int *LDA,
               double *B, int *LDB, double *WORK, int *LWORK, int *INFO);
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
}
                                                                                
// ************************************************************************
// constructor
// ------------------------------------------------------------------------
PCAnalyzer::PCAnalyzer() : Analyzer()
{
   setName("PCA");
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO,"Principal component analysis (PCA) \n");
   printOutTS(PL_INFO,
        "This function performs a PCA of all sample outputs\n"); 
   printOutTS(PL_INFO,"such that YY = U * S * V\n\n");
   printOutTS(PL_INFO,"where YY is the normalized sample outputs,\n");
   printOutTS(PL_INFO,"      U are the left singular vectors,\n");
   printOutTS(PL_INFO,"      S are the singular values, and\n");
   printOutTS(PL_INFO,"      V are the right singular vectors.\n");
   printOutTS(PL_INFO,"At the end, the sample outputs will be replaced\n");
   printOutTS(PL_INFO,"by Y = Y * V^t \n");
   printOutTS(PL_INFO,"If the analysis expert mode is on, you have the\n");
   printOutTS(PL_INFO,"option to store these matrices into matlab file.\n");
   printAsterisks(PL_INFO, 0);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
PCAnalyzer::~PCAnalyzer()
{
}

// ************************************************************************
// perform analysis
// ------------------------------------------------------------------------
double PCAnalyzer::analyze(aData &adata)
{
   int    nInputs, nOutputs, nSamples, jj, ss, ii, M, N, info;
   int    wlen, pcCnt, pcaFlag=0, count;
   double *Y, *YY, *means, *stdevs, *WW, *SS, *VV, *UU, dtemp;
   char   jobu, jobvt, winput1[500], winput2[500], pcaFile[500], *cString;
   char   pString[500];
   FILE   *fp = NULL;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   Y        = adata.sampleOutputs_;
   if (adata.inputPDFs_ != NULL)
   {
      count = 0;
      for (ii = 0; ii < nInputs; ii++) count += adata.inputPDFs_[ii];
      if (count > 0)
      {
         printOutTS(PL_INFO,
              "PCA INFO: some inputs have non-uniform PDFs, but\n");
         printOutTS(PL_INFO,
              "          they are not relevant in this analysis.\n");
      }
   }

   if (nInputs <= 0 || nOutputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR, "PCA ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "    nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   if (nOutputs == 1)
   {
      printOutTS(PL_ERROR, 
           "PCA ERROR: analysis not done since nOutputs=1.\n");
      return PSUADE_UNDEFINED;
   }
   
   if (psAnaExpertMode_ == 1 && psPlotTool_ == 0)
   {
      sprintf(pString,"Write PCA information to a file? (y or n) ");
      getString(pString, winput1);
      if (winput1[0] == 'y')
      {
         sprintf(pString,"Enter the matlab file name : ");
         getString(pString, pcaFile);
         pcaFile[strlen(pcaFile)-1] = '\0';
         fp = fopen(pcaFile, "w");
         if (fp != NULL)
         {
            fclose(fp);
            pcaFlag = 1;
         }
         else printOutTS(PL_INFO,"PCA: cannot open matlab file %s\n", 
                   pcaFile);
      }
   }
   else
   {
      if (psConfig_ != NULL)
      {
         cString = psConfig_->getParameter("PCA_matlab_file");
         if (cString != NULL)
         {
            sscanf(cString, "%s %s %s", winput1, winput2, pcaFile);
            fp = fopen(pcaFile, "w");
            if (fp != NULL)
            {
               fclose(fp);
               pcaFlag = 1;
            }
            else 
               printOutTS(PL_INFO,
                    "PCA ERROR: cannot open matlab file %s\n",pcaFile);
         }
      }
   }

   YY = new double[nSamples*nOutputs];
   checkAllocate(YY, "YY in PCA::analyze");
   for (jj = 0; jj < nOutputs; jj++)
      for (ss = 0; ss < nSamples; ss++)
         YY[nSamples*jj+ss] = Y[ss*nOutputs+jj];

   means  = new double[nOutputs];
   stdevs = new double[nOutputs];
   checkAllocate(stdevs, "stdevs in PCA::analyze");
   for (jj = 0; jj < nOutputs; jj++)
   {
      means[jj] = 0.0;
      for (ss = 0; ss < nSamples; ss++) means[jj] += YY[nSamples*jj+ss];
      means[jj] /= (double) nSamples;
      for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] -= means[jj];
      stdevs[jj] = 0.0;
      for (ss = 0; ss < nSamples; ss++)
         stdevs[jj] += YY[nSamples*jj+ss] * YY[nSamples*jj+ss];
      stdevs[jj] /= (double) (nSamples - 1);
      stdevs[jj] = sqrt(stdevs[jj]);
      if (PABS(stdevs[jj]) > 1.0e-14)
         for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] /= stdevs[jj];
      printOutTS(PL_INFO,"PCA: output %5d has mean = %e, std. dev = %e\n", 
                 jj+1, means[jj], stdevs[jj]);
   }

   M = nSamples;
   N = nOutputs;
   jobu = 'S';
   jobvt = 'A';
   wlen = 5 * M;
   SS = new double[N];
   UU = new double[M*N];
   VV = new double[N*N];
   WW = new double[wlen];
   checkAllocate(WW, "WW in PCA::analyze");
   dgesvd_(&jobu, &jobvt, &M, &N, YY, &M, SS, UU, &M, VV, &N, WW,
           &wlen, &info);
   if (info != 0)
      printOutTS(PL_INFO, 
           "* PCA INFO: dgesvd returns a nonzero (%d).\n",info);
   for (ii = 0; ii < N; ii++) SS[ii] = SS[ii] * SS[ii];
   for (ii = 0; ii < N; ii++)
      printOutTS(PL_INFO, 
         "principal component %3d has variance = %16.8e\n",ii+1,SS[ii]);
   sprintf(pString,"Enter how many principal components to keep : ");
   pcCnt = getInt(1, N, pString);

   for (jj = 0; jj < nOutputs; jj++)
      for (ss = 0; ss < nSamples; ss++)
         YY[nSamples*jj+ss] = Y[ss*nOutputs+jj];
   for (jj = 0; jj < nOutputs; jj++)
   {
      means[jj] = 0.0;
      for (ss = 0; ss < nSamples; ss++) means[jj] += YY[nSamples*jj+ss];
      means[jj] /= (double) nSamples;
      for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] -= means[jj];
      stdevs[jj] = 0.0;
      for (ss = 0; ss < nSamples; ss++)
         stdevs[jj] += YY[nSamples*jj+ss] * YY[nSamples*jj+ss];
      stdevs[jj] /= (double) (nSamples - 1);
      stdevs[jj] = sqrt(stdevs[jj]);
      if (PABS(stdevs[jj]) > 1.0e-14)
         for (ss = 0; ss < nSamples; ss++) YY[nSamples*jj+ss] /= stdevs[jj];
   }


   if (pcaFlag == 1)
   {
      fp = fopen(pcaFile, "w");
      if (fp != NULL)
      {
         fprintf(fp,"%% Principal component analysis results\n");
         fprintf(fp,"%% First display the Scree diagram\n");
         fprintf(fp,"%% A  - original sample outputs (normalized)\n");
         fprintf(fp,"%% Am - mean of each column of A\n");
         fprintf(fp,"%% As - std dev of each column of A\n");
         fprintf(fp,"%% V  - eigenvector matrix\n");
         fprintf(fp,"%% U  - left singular vector matrix\n");
         fprintf(fp,"%% S  - singular values (column vector)\n");

         fprintf(fp,"A = [\n");
         for (ss = 0; ss < M; ss++) 
         {
            for (jj = 0; jj < N; jj++) fprintf(fp,"%e ", YY[M*jj+ss]);
            fprintf(fp, "\n");
         }
         fprintf(fp,"];\n");
         fprintf(fp,"Am = [\n");
         for (ii = 0; ii < N; ii++) fprintf(fp," %e\n", means[ii]);
         fprintf(fp,"];\n");
         fprintf(fp,"As = [\n");
         for (ii = 0; ii < N; ii++) fprintf(fp," %e\n", stdevs[ii]);
         fprintf(fp,"];\n");
         fprintf(fp, "S = [\n");
         for (ii = 0; ii < pcCnt; ii++) fprintf(fp," %e\n", SS[ii]);
         fprintf(fp,"];\n");
         fprintf(fp,"hold off\n");
         fprintf(fp,"plot(S)\n");
         fprintf(fp,"grid\n");
         fprintf(fp,"xlabel('Principal component number')\n");
         fprintf(fp,"ylabel('Variances')\n");
         fprintf(fp,"title('Scree Diagram')\n");
         fprintf(fp,"disp('Press enter to continue')\n");
         fprintf(fp,"pause\n");

         fprintf(fp,"V = [\n");
         for (ii = 0; ii < N; ii++)
         {
            for (jj = 0; jj < pcCnt; jj++) fprintf(fp," %e", VV[jj*N+ii]);
            fprintf(fp, "\n");
         }
         fprintf(fp,"];\n");
         fprintf(fp,"U = [\n");
         for (ii = 0; ii < nSamples; ii++)
         {
            for (jj = 0; jj < pcCnt; jj++)
               fprintf(fp," %e", UU[jj*M+ii]);
            fprintf(fp,"\n");
         }
         fprintf(fp,"];\n");
         fprintf(fp,"disp('Display each principal component')\n");
         fprintf(fp,"US = U * diag(S);\n");
         fprintf(fp,"XU = 1:%d;\n",M);
         fprintf(fp,"for ii = 1 : %d\n", pcCnt);
         fprintf(fp,"   plot(XU,US(:,ii),XU,US(:,ii),'x')\n");
         fprintf(fp,"   disp(['principal component = ' int2str(ii)])\n");
         fprintf(fp,"   disp('Press enter to continue')\n");
         fprintf(fp,"   pause\n");
         fprintf(fp,"end;\n");
         fprintf(fp,"hold off\n");
         fprintf(fp,"%%barh(abs(A)','stacked')\n");
         fprintf(fp,"%%area(XA,abs(A))\n");
         fprintf(fp,"hold off\n");
         fclose(fp);
      }
      else printOutTS(PL_WARN,"WARNING: cannot open file for plotting.\n");
   }

   for (ss = 0; ss < nSamples; ss++)
   {
      for (jj = 0; jj < pcCnt; jj++)
      {
         dtemp = 0.0;
         for (ii = 0; ii < N; ii++)
            dtemp += VV[jj*N+ii] * YY[ii*nSamples+ss];
         Y[pcCnt*ss+jj] = dtemp;
      }
   }

	   fp = fopen("psPCA.out", "w");
   if (fp == NULL)
   {
      printOutTS(PL_ERROR, 
           "PCA ERROR: failed to store principal components.\n");
   }
   else
   {
      fprintf(fp,"# Let Y = normalized sample outputs\n");
      fprintf(fp,"# Let Y = U S V^t\n");
      fprintf(fp,"# This file contains Y * V\n");
      fprintf(fp,"%d %d\n", nSamples, pcCnt);
      for (ss = 0; ss < nSamples; ss++)
      {
         for (jj = 0; jj < pcCnt; jj++)
            fprintf(fp, "%24.16e ", Y[pcCnt*ss+jj]);
         fprintf(fp, "\n");
      }
      fclose(fp);
      printOutTS(PL_INFO, 
           "PCA: principal components (Y*V) stored in psPCA.out file.\n");
   }

   adata.nOutputs_ = pcCnt;
   delete [] SS;
   delete [] WW;
   delete [] VV;
   delete [] UU;
   delete [] YY;
   delete [] means;
   delete [] stdevs;
   return 0.0;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
PCAnalyzer& PCAnalyzer::operator=(const PCAnalyzer &)
{
   printOutTS(PL_ERROR, "PCA operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

