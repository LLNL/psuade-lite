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
// Functions for the class AnovaAnalyzer (modified from my work in DDace) 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "FuncApprox.h"
#include "AnovaAnalyzer.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
AnovaAnalyzer::AnovaAnalyzer() : Analyzer(), tableLeng_(0), dofs_(0), 
                   sumSquares_(0), meanSquares_(0), fValues_(0), code_(0)
{
   setName("ANOVA");
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
AnovaAnalyzer::~AnovaAnalyzer()
{
   if (dofs_) delete [] dofs_;
   if (sumSquares_) delete [] sumSquares_;
   if (meanSquares_) delete [] meanSquares_;
   if (fValues_) delete [] fValues_;
   if (code_)
   {
      for (int jj = 0; jj < tableLeng_; jj++) delete [] code_[jj];
      delete [] code_;
   }
}

// ************************************************************************
// perform Anova
// ------------------------------------------------------------------------
double AnovaAnalyzer::analyze(aData &adata)
{
   int        nInputs, nOutputs, nSamples, outputID, ss;
   int        ii, jj, kk, nPtsPerDim=256, ncount, index, nDataPts;
   int        inCnt, inCnt2, inCnt3, nLevels=10;
   int        ncnt, cind1, cind2, cind3;
   int        dimPts;
   double     CT, dtemp, ymax, *X;
   double     *xInterpolated, *yInterpolated, dvalue, *Y;
   double     *YY, *xLower, *xUpper;
   FuncApprox *fa;

   nInputs  = adata.nInputs_;
   nOutputs = adata.nOutputs_;
   nSamples = adata.nSamples_;
   xLower   = adata.iLowerB_;
   xUpper   = adata.iUpperB_;
   X        = adata.sampleInputs_;
   Y        = adata.sampleOutputs_;
   outputID = adata.outputID_;
   if (adata.inputPDFs_ != NULL)
   {
      ncount = 0;
      for (ii = 0; ii < nInputs; ii++) ncount += adata.inputPDFs_[ii];
      if (ncount > 0)
      {
         printOutTS(PL_WARN,
              "ANOVA INFO: some inputs have non-uniform PDFs.\n");
         printOutTS(PL_WARN,
              "      However, they are not relevant in this analysis\n");
         printOutTS(PL_WARN,
              "      (the sample should have been generated with\n");
         printOutTS(PL_WARN, "       the desired distributions.)\n");
      }
   }

   if (nInputs <= 0 || nSamples <= 0)
   {
      printOutTS(PL_ERROR, "ANOVA ERROR: invalid arguments.\n");
      printOutTS(PL_ERROR, "    nInputs  = %d\n", nInputs);
      printOutTS(PL_ERROR, "    nOutputs = %d\n", nOutputs);
      printOutTS(PL_ERROR, "    nSamples = %d\n", nSamples);
      return PSUADE_UNDEFINED;
   } 
   for (ii = 0; ii < nInputs; ii++)
   {
      if (xLower[ii] >= xUpper[ii])
      {
         printOutTS(PL_ERROR, "ANOVA ERROR: invalid input bounds.\n");
         printOutTS(PL_ERROR, "             Input %3d bounds = %e %e\n", 
                    ii+1, xLower[ii], xUpper[ii]);
         return PSUADE_UNDEFINED;
      }
   }
   
   if (dofs_)        delete [] dofs_;
   if (sumSquares_)  delete [] sumSquares_;
   if (meanSquares_) delete [] meanSquares_;
   if (fValues_)     delete [] fValues_;
   if (code_)
   {
      for (int jj = 0; jj < tableLeng_; jj++) delete [] code_[jj];
      delete [] code_;
   }
   code_    = NULL;
   fValues_ = NULL;
   meanSquares_ = NULL;
   sumSquares_ = NULL;
   dofs_ = NULL;
 
   YY = new double[nSamples];
   checkAllocate(YY, "YY in AnovaAnalyzer::analyze");

   for (ss = 0; ss < nSamples; ss++) YY[ss] = Y[nOutputs*ss+outputID];
   dimPts = 1000001;
   while (dimPts > 1000000)
   {
      nPtsPerDim = nPtsPerDim / 2;
      dimPts = 1;
      for (ii = 0; ii < nInputs; ii++)
      {
         dimPts *= nPtsPerDim;
         if (dimPts > 1000000) break;
      }
   }
   fa = genFA(PSUADE_RS_MARS, nInputs, 1, nSamples);
   fa->setNPtsPerDim(nPtsPerDim);
   fa->setBounds(xLower, xUpper);
   fa->initialize(X, YY);
   
   nLevels  = nPtsPerDim;
   nDataPts = nLevels;
   for (ii = 1; ii < nInputs; ii++) nDataPts *= nLevels;
   xInterpolated = new double[nDataPts*nInputs];
   yInterpolated = new double[nDataPts];
   checkAllocate(yInterpolated, "yInterpolated in AnovaAnalyzer::analyze");
   ncount = nDataPts / nLevels;
   for (ii = 0; ii < nInputs; ii++)
   {
      for (jj = 0; jj < nDataPts; jj+=ncount)
      {
         index = (jj / ncount) % nLevels;
         dvalue = (xUpper[ii]-xLower[ii]) / (nLevels-1)*index +
                   xLower[ii]; 
         for (kk = 0; kk < ncount; kk++) 
            xInterpolated[(jj+kk)*nInputs+ii] = dvalue;
      }
      ncount /= nLevels;
   }

   ymax = -1.0E20;
   for (jj = 0; jj < nDataPts; jj++)
   {
      yInterpolated[jj] = fa->evaluatePoint(&xInterpolated[jj*nInputs]);
      if (yInterpolated[jj] > ymax) ymax = yInterpolated[jj];
   }
   if (ymax != 0.0)
   {
      for (jj = 0; jj < nDataPts; jj++) yInterpolated[jj] /= ymax;
   }

   tableLeng_ = nInputs;
   tableLeng_ += ((nInputs - 1) * nInputs / 2);
   for (ii = 0; ii < nInputs; ii++)
      for (jj = ii+1; jj < nInputs; jj++)
         for (kk = jj+1; kk < nInputs; kk++) tableLeng_++;
   tableLeng_ += 2;
   dofs_ = new int[tableLeng_];
   code_ = new int*[tableLeng_];
   for (jj = 0; jj < tableLeng_; jj++) code_[jj] = new int[3];
   checkAllocate(code_[tableLeng_-1], "code in AnovaAnalyzer::analyze");
   for (ii = 0; ii < nInputs; ii++)
      dofs_[ii] = nLevels - 1;
   ncnt = nInputs;
   for (ii = 0; ii < nInputs; ii++)
      for (jj = ii+1; jj < nInputs; jj++)
         dofs_[ncnt++] = dofs_[ii] * dofs_[jj];
   for (ii = 0; ii < nInputs; ii++)
      for (jj = ii+1; jj < nInputs; jj++)
         for (kk = jj+1; kk < nInputs; kk++)
            dofs_[ncnt++] = dofs_[ii] * dofs_[jj] * dofs_[kk];
   dofs_[tableLeng_-1] = nDataPts - 1;
   dofs_[tableLeng_-2] = nDataPts - 1;
   ncnt = nInputs;
   ncnt += ((nInputs - 1) * nInputs / 2);
   for (jj = 0; jj < ncnt; jj++) dofs_[tableLeng_-2] -= dofs_[jj];

   sumSquares_  = new double[tableLeng_];
   meanSquares_ = new double[tableLeng_];
   fValues_     = new double[tableLeng_];
   checkAllocate(fValues_, "fValues in AnovaAnalyzer::analyze");
   ncnt        = 0;
   printAsterisks(PL_INFO, 0);
   printOutTS(PL_INFO, "*                       ANOVA table\n");
   printOutTS(PL_INFO, "*              (based on MARS interpolation)\n");
   printEquals(PL_INFO, 0);

   for (inCnt = 0; inCnt < nInputs; inCnt++)
   {
      sumSquares_[inCnt] = computeSumSquares1(nDataPts, nInputs, inCnt,
                                    dofs_[inCnt], 1, 0, xInterpolated,
                                    yInterpolated);
      code_[ncnt][0] = inCnt;
      code_[ncnt][1] = 0;
      code_[ncnt][2] = 0;
      ncnt++;
   }
   for (inCnt = 0; inCnt < nInputs; inCnt++)
   {
      for (inCnt2 = inCnt+1; inCnt2 < nInputs; inCnt2++)
      {
         sumSquares_[ncnt] = computeSumSquares2(nDataPts, nInputs,
                                   inCnt, inCnt2, dofs_[inCnt],
                                   dofs_[inCnt2], 1, 0, xInterpolated,
                                   yInterpolated, sumSquares_[inCnt],
                                   sumSquares_[inCnt2]);
         code_[ncnt][0] = inCnt;
         code_[ncnt][1] = inCnt2;
         code_[ncnt][2] = 0;
         ncnt++;
      }
   }
   for (inCnt = 0; inCnt < nInputs; inCnt++)
   {
      for (inCnt2 = inCnt+1; inCnt2 < nInputs; inCnt2++)
      {
         for (inCnt3 = inCnt2+1; inCnt3 < nInputs; inCnt3++)
         {
            kk = nInputs;
            for (ii = 0; ii < nInputs; ii++)
            {
               for (jj = ii+1; jj < nInputs; jj++)
               {
                  if (ii == inCnt && jj == inCnt2)
                  {
                     cind1 = kk; break;
                  } else kk++;
               }
            }
            kk = nInputs;
            for (ii = 0; ii < nInputs; ii++)
            {
               for (jj = ii+1; jj < nInputs; jj++)
               {
                  if (ii == inCnt && jj == inCnt3)
                  {
                     cind2 = kk; break;
                  } else kk++;
               }
            }
            kk = nInputs;
            for (ii = 0; ii < nInputs; ii++)
            {
               for (jj = ii+1; jj < nInputs; jj++)
               {
                  if (ii == inCnt2 && jj == inCnt3)
                  {
                     cind3 = kk; break;
                  } else kk++;
               }
            }
            if (cind1 < 0 || cind2 < 0 || cind3 < 0 || cind1 >= tableLeng_ ||
                cind2 >= tableLeng_ || cind3 >= tableLeng_)
            {
               printOutTS(PL_ERROR, "ANOVA ERROR: unrecoverable error.\n");
               //delete [] dofs_;
               //delete [] sumSquares_;
               //delete [] meanSquares_;
               //delete [] fValues_;
               //for (jj = 0; jj < tableLeng_; jj++) delete [] code_[jj];
               //delete [] code_;
               delete [] xInterpolated;
               delete [] yInterpolated;
               delete fa;
               return -1;
            }
            sumSquares_[ncnt] = computeSumSquares3(nDataPts, nInputs,
                                  inCnt, inCnt2, inCnt3, dofs_[inCnt],
                                  dofs_[inCnt2], dofs_[inCnt3], 1, 0,
                                  xInterpolated, yInterpolated,
                                  sumSquares_[inCnt], sumSquares_[inCnt2],
                                  sumSquares_[inCnt3], sumSquares_[cind1],
                                  sumSquares_[cind2], sumSquares_[cind3]);
            code_[ncnt][0] = inCnt;
            code_[ncnt][1] = inCnt2;
            code_[ncnt][2] = inCnt3;
            ncnt++;
         }
      }
   }
   dtemp = 0.0;
   for (jj = 0; jj < nDataPts; jj++)
      dtemp += (yInterpolated[jj] * yInterpolated[jj]);
   CT = 0.0;
   for (jj = 0; jj < nDataPts; jj++) CT += yInterpolated[jj];
   CT = (CT * CT) / nDataPts;
   sumSquares_[tableLeng_-1] = dtemp - CT;
   sumSquares_[tableLeng_-2] = sumSquares_[tableLeng_-1];
   ncnt = nInputs;
   ncnt += ((nInputs - 1) * nInputs / 2);
   for (jj = 0; jj < ncnt; jj++)
      sumSquares_[tableLeng_-2] -= sumSquares_[jj];
   for (jj = 0; jj < tableLeng_; jj++)
   {
      if (dofs_[jj] > 0)
         meanSquares_[jj] = sumSquares_[jj] / (double) dofs_[jj];
      else
         meanSquares_[jj] = 0.0;
   }
   for (jj = 0; jj < tableLeng_-2; jj++)
   {
      if (meanSquares_[tableLeng_-2] != 0.0)
         fValues_[jj] = meanSquares_[jj] / meanSquares_[tableLeng_-2];
      else
         fValues_[jj] = 0.0;
   }

   printEquals(PL_INFO, 65);
   printOutTS(PL_INFO, 
       "|  source of | deg. of|   sum of    |   mean      |            |\n");
   printOutTS(PL_INFO, 
       "|  variation | freedom|   squares   |   square    |       F    |\n");
   printDashes(PL_INFO, 65);
   for (jj = 0; jj < tableLeng_-2; jj++)
   {
      if (jj < nInputs)
      {
         printOutTS(PL_WARN,"|    %7d |%7d | %11.4e | %11.4e | %11.4e|\n",
            jj+1, dofs_[jj], sumSquares_[jj], meanSquares_[jj],fValues_[jj]);
      } 
      else if (code_[jj][2] == 0)
      {
         printOutTS(PL_WARN,"|    %3d,%3d |%7d | %11.4e | %11.4e | %11.4e|\n",
              code_[jj][0]+1,code_[jj][1]+1,dofs_[jj],sumSquares_[jj],
              meanSquares_[jj], fValues_[jj]);
      } 
#if 0
      else
      {
         printOutTS(PL_INFO,
             "| %3d,%3d,%3d|   %6d | %11.4e | %11.4e | %11.4e|\n",
             code_[jj][0]+1,code_[jj][1]+1,code_[jj][2],dofs_[jj],
             sumSquares_[jj],meanSquares_[jj],fValues_[jj]);
      }
#endif
   }
   printOutTS(PL_INFO,"|   others   |%7d | %11.4e | %11.4e |    -----   |\n",
          dofs_[tableLeng_-2],sumSquares_[tableLeng_-2],
          meanSquares_[tableLeng_-2]);
   printOutTS(PL_INFO,"|   total    |%7d | %11.4e | %11.4e |    -----   |\n",
          dofs_[tableLeng_-1],sumSquares_[tableLeng_-1],
          meanSquares_[tableLeng_-1]);
   printAsterisks(PL_INFO, 0);

   delete [] xInterpolated;
   delete [] yInterpolated;
   delete fa;
   return 0.0;
}

// *************************************************************************
// Compute sum of squares for one factor
// -------------------------------------------------------------------------
double AnovaAnalyzer::computeSumSquares1(int nSamples, int nInputs, 
                       int source, int dof, int nOutputs, int outindex, 
                       double *X, double *Y)
{
   int    ss, jj, nLevels, xcnt, found;
   double mean, sumT2, CT, *gsum, *xvalues;

   CT = 0.0;
   for (ss = 0; ss < nSamples; ss++) CT += Y[nOutputs*ss+outindex];
   mean = CT / nSamples;

   nLevels = dof + 1;
   xcnt    = 0;
   xvalues = new double[nLevels];
   checkAllocate(xvalues, "xvalues in AnovaAnalyzer::computeSumSquares1");
   for (ss = 0; ss < nSamples; ss++)
   {
      found = 0;
      for (jj = 0; jj < xcnt; jj++)
      {
         if (PABS(xvalues[jj] - X[nInputs*ss+source]) < 1.0E-8)
         {
            found = 1;
            break;
         }
      }
      if (found == 0) xvalues[xcnt++] = X[nInputs*ss+source];
   }
   if (xcnt != nLevels) 
   {
      printOutTS(PL_INFO, "ANOVA ERROR (1): %d\n",xcnt);
      exit(0);
   }

   gsum = new double[nLevels];
   checkAllocate(gsum, "gsum in AnovaAnalyzer::computeSumSquares1");
   for (jj = 0; jj < nLevels; jj++) gsum[jj] = 0.0;
   for (ss = 0; ss < nSamples; ss++)
   {
      for (jj = 0; jj < xcnt; jj++)
         if (PABS(xvalues[jj]-X[nInputs*ss+source]) < 1.0E-8) break;
      gsum[jj] += Y[nOutputs*ss+outindex];
   }
   sumT2 = 0.0;
   for (jj = 0; jj < nLevels; jj++)
   {
      gsum[jj] /= (nSamples / nLevels);
      sumT2 += ((gsum[jj] - mean) * (gsum[jj] - mean));
   }
   sumT2 = sumT2 / nLevels * nSamples;
   delete [] gsum;
   delete [] xvalues;

   return sumT2;
}

// *************************************************************************
// Compute sum of squares for two factors
// -------------------------------------------------------------------------
double AnovaAnalyzer::computeSumSquares2(int nSamples, int nInputs, 
                       int input1, int input2, int dof1, int dof2, 
                       int nOutputs, int outindex, double *X, double *Y,
                       double ss1, double ss2)
{
   int    ss, jj, kk, nlevels1, nlevels2, ncnt, xcnt1, xcnt2, found;
   double mean, sumT2, CT, *gsum12, *xvalues1, *xvalues2;
   double *gsum1, *gsum2;

   CT = 0.0;
   for (ss = 0; ss < nSamples; ss++) CT += Y[nOutputs*ss+outindex];
   mean = CT / nSamples;

   nlevels1 = dof1 + 1;
   nlevels2 = dof2 + 1;
   xcnt1    = 0;
   xcnt2    = 0;
   xvalues1 = new double[nlevels1];
   xvalues2 = new double[nlevels2];
   checkAllocate(xvalues2, "xvalues2 in AnovaAnalyzer::computeSumSquares2");
   for (ss = 0; ss < nSamples; ss++)
   {
      found = 0;
      for (jj = 0; jj < xcnt1; jj++)
      {
         if (PABS(xvalues1[jj] - X[nInputs*ss+input1]) < 1.0E-8)
         {
            found = 1;
            break;
         }
      }
      if (found == 0) xvalues1[xcnt1++] = X[nInputs*ss+input1];
   }
   if (xcnt1 != nlevels1)
   {
      printOutTS(PL_INFO, "ANOVA ERROR (2a): %d\n",nlevels1);
      exit(0);
   }
   for (ss = 0; ss < nSamples; ss++)
   {
      found = 0;
      for (jj = 0; jj < xcnt2; jj++)
      {
         if (PABS(xvalues2[jj] - X[nInputs*ss+input2]) < 1.0E-8) 
         {
            found = 1;
            break;
         }
      }
      if (found == 0) xvalues2[xcnt2++] = X[nInputs*ss+input2];
   }
   if (xcnt2 != nlevels2)
   {
      printOutTS(PL_INFO, "AnovaAnalyzer ERROR (2b): %d\n",nlevels2);
      exit(0);
   }

   gsum1 = new double[nlevels1];
   checkAllocate(gsum1, "gsum1 in AnovaAnalyzer::computeSumSquares2");
   for (jj = 0; jj < nlevels1; jj++) gsum1[jj] = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      for (jj = 0; jj < xcnt1; jj++)
         if (PABS(xvalues1[jj]-X[nInputs*ss+input1]) < 1.0E-8) break;
      gsum1[jj] += Y[nOutputs*ss+outindex];
   }
   for (jj = 0; jj < nlevels1; jj++) gsum1[jj] /= (nSamples / nlevels1);

   gsum2 = new double[nlevels2];
   checkAllocate(gsum2, "gsum2 in AnovaAnalyzer::computeSumSquares2");
   for (jj = 0; jj < nlevels2; jj++) gsum2[jj] = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      for (jj = 0; jj < xcnt2; jj++)
         if (PABS(xvalues2[jj]-X[nInputs*ss+input2]) < 1.0E-8) break;
      gsum2[jj] += Y[nOutputs*ss+outindex];
   }
   for (jj = 0; jj < nlevels2; jj++) gsum2[jj] /= (nSamples / nlevels2);

   ncnt = nlevels1 * nlevels2;
   gsum12 = new double[ncnt];
   checkAllocate(gsum12, "gsum12 in AnovaAnalyzer::computeSumSquares2");
   for (jj = 0; jj < ncnt; jj++) gsum12[jj] = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      kk = -1;
      for (jj = 0; jj < xcnt1; jj++)
      {
         found = 0;
         if (PABS(xvalues1[jj] - X[nInputs*ss+input1]) < 1.0E-8)
         {
            for (kk = 0; kk < xcnt2; kk++)
            {
                if (PABS(xvalues2[kk]-X[nInputs*ss+input2]) < 1.0E-8)
                {
                    found = 1;
                    break;
                }
             }
          }
          if (found == 1) break;
      }
      if (kk >= 0) gsum12[jj*xcnt2+kk] += Y[nOutputs*ss+outindex];
   }
   for (jj = 0; jj < nlevels1*nlevels2; jj++) 
      gsum12[jj] /= (nSamples / (nlevels1 * nlevels2));

   sumT2 = 0.0;
   for (jj = 0; jj < nlevels1; jj++) 
   {
      for (kk = 0; kk < nlevels2; kk++) 
         sumT2 += ((gsum12[jj*nlevels2+kk] - gsum1[jj] - gsum2[kk] + mean)*
                   (gsum12[jj*nlevels2+kk] - gsum1[jj] - gsum2[kk] + mean));
   }
   sumT2 = sumT2 / (nSamples / (nlevels1 * nlevels2));
   delete [] gsum1;
   delete [] gsum2;
   delete [] gsum12;
   delete [] xvalues1;
   delete [] xvalues2;

   return sumT2;
}

// *************************************************************************
// Compute sum of squares for three factors
// -------------------------------------------------------------------------
double AnovaAnalyzer::computeSumSquares3(int nSamples, int nInputs, 
                          int input1, int input2, int input3, int dof1,
                          int dof2, int dof3, int nOutputs,
                          int output, double *X, double *Y,
                          double ss1, double ss2, double ss3, double ss4,
                          double ss5, double ss6)
{
   int    nlevels1, nlevels2, nlevels3, ncnt, xcnt1, xcnt2;
   int    xcnt3, found, ss, jj, kk, mm;
   double mean, sumT2, CT, *gsum, *xvalues1, *xvalues2, *xvalues3;
   double sss1, sss2, sss3, sss4, sss5, sss6, diff;

   CT = 0.0;
   for (ss = 0; ss < nSamples; ss++) CT += Y[nOutputs*ss+output];
   mean = CT / nSamples;
   CT = (CT * CT) / nSamples;

   nlevels1 = dof1 + 1;
   nlevels2 = dof2 + 1;
   nlevels3 = dof3 + 1;
   xcnt1    = 0;
   xcnt2    = 0;
   xcnt3    = 0;
   xvalues1 = new double[nlevels1];
   xvalues2 = new double[nlevels2];
   xvalues3 = new double[nlevels2];
   checkAllocate(xvalues3,"xvalues3 in AnovaAnalyzer::computeSumSquares3");
   for (ss = 0; ss < nSamples; ss++)
   {
      found = 0;
      for (jj = 0; jj < xcnt1; jj++)
      {
         if (PABS(xvalues1[jj]-X[nInputs*ss+input1]) < 1.0E-8) 
         {
            found = 1;
            break;
         }
      }
      if (found == 0) xvalues1[xcnt1++] = X[nInputs*ss+input1];
   }
   if (xcnt1 != nlevels1) 
   {
      printOutTS(PL_WARN, "ANOVA ERROR (3a): %d\n",nlevels1);
      exit(0);
   }
   for (ss = 0; ss < nSamples; ss++)
   {
      found = 0;
      for (jj = 0; jj < xcnt2; jj++)
      {
         if (PABS(xvalues2[jj]-X[nInputs*ss+input2]) < 1.0E-8) 
         {
            found = 1;
            break;
         }
      }
      if (found == 0) xvalues2[xcnt2++] = X[nInputs*ss+input2];
   }
   if (xcnt2 != nlevels2) 
   {
      printOutTS(PL_INFO, "ANOVA ERROR (3b): %d\n",nlevels2);
      exit(1);
   }
   for (ss = 0; ss < nSamples; ss++)
   {
      found = 0;
      for (jj = 0; jj < xcnt3; jj++)
      {
         if (PABS(xvalues3[jj]-X[nInputs*ss+input3]) < 1.0E-8) 
         {
            found = 1;
            break;
         }
      }
      if (found == 0) xvalues3[xcnt3++] = X[nInputs*ss+input3];
   }
   if (xcnt3 != nlevels3) 
   {
      printOutTS(PL_INFO, "ANOVA ERROR (3c): %d\n",nlevels3);
      exit(0);
   }

   ncnt = nlevels1 * nlevels2 * nlevels3;
   gsum = new double[ncnt];
   checkAllocate(gsum,"gsum in AnovaAnalyzer::computeSumSquares3");
   for (jj = 0; jj < ncnt; jj++) gsum[jj] = 0.0;
   for (ss = 0; ss < nSamples; ss++) 
   {
      kk = mm = -1;
      for (jj = 0; jj < xcnt1; jj++)
      {
         found = 0;
         if (PABS(xvalues1[jj]-X[nInputs*ss+input1]) < 1.0E-8)
         {
            for (kk = 0; kk < xcnt2; kk++)
            {
                if (PABS(xvalues2[kk]-X[nInputs*ss+input2]) < 1.0E-8)
                {
                    for (mm = 0; mm < xcnt3; mm++)
                    {
                        diff = xvalues3[mm] - X[nInputs*ss+input3];
                        if (PABS(diff) < 1.0E-8)
                        {
                           found = 1;
                           break;
                        }
                    }
                    if (found == 1) break;
                }
             }
          }
          if (found == 1) break;
      }
      if (kk < 0 || mm < 0) printOutTS(PL_INFO, "ANOVA ERROR (4).\n");
      else
         gsum[jj*xcnt2*xcnt3+kk*xcnt3+mm] += Y[nOutputs*ss+output];
   }
   sumT2 = 0.0;
   for (jj = 0; jj < ncnt; jj++) sumT2 += (gsum[jj] * gsum[jj]);
   sumT2 = sumT2 * nlevels1 * nlevels2 * nlevels3 / nSamples;

   sss1  = ss1 + CT;
   sss2  = ss2 + CT;
   sss3  = ss3 + CT;
   sss4  = ss4 + sss1 + sss2 - CT;
   sss5  = ss5 + sss1 + sss3 - CT;
   sss6  = ss6 + sss2 + sss3 - CT;
   sumT2 = sumT2 + sss1 + sss2 + sss3 - sss4 - sss5 - sss6 - CT;
   delete [] gsum;
   delete [] xvalues1;
   delete [] xvalues2;
   delete [] xvalues3;

   return sumT2;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
AnovaAnalyzer& AnovaAnalyzer::operator=(const AnovaAnalyzer &)
{
   printOutTS(PL_ERROR, "ANOVA operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}
// ************************************************************************
// Getters for analysis results
// ------------------------------------------------------------------------

int AnovaAnalyzer::get_tableLeng()
{
   return tableLeng_;
}

int *AnovaAnalyzer::get_dofs()
{
   int* retVal = NULL;
   if(dofs_)
   {
      retVal = new int[tableLeng_];
      checkAllocate(retVal,"retVal in AnovaAnalyzer::get_dofs");
      std::copy(dofs_, dofs_+tableLeng_, retVal);
   }
   return retVal;
}

double *AnovaAnalyzer::get_sumSquare()
{
   double* retVal = NULL;
   if (sumSquares_)
   {
      retVal = new double[tableLeng_];
      checkAllocate(retVal,"retVal in AnovaAnalyzer::get_sumSquare");
      std::copy(sumSquares_, sumSquares_+tableLeng_, retVal);
   }
   return retVal;
}

double *AnovaAnalyzer::get_meanSquares()
{
   double* retVal = NULL;
   if (meanSquares_)
   {
      retVal = new double[tableLeng_];
      checkAllocate(retVal,"retVal in AnovaAnalyzer::get_meanSquares");
      std::copy(meanSquares_, meanSquares_+tableLeng_, retVal);
   }
   return retVal;
}

double *AnovaAnalyzer::get_fValues()
{
   double* retVal = NULL;
   if (fValues_)
   {
      retVal = new double[tableLeng_];
      checkAllocate(retVal,"retVal in AnovaAnalyzer::get_fValues");
      std::copy(fValues_, fValues_+tableLeng_, retVal);
   }
   return retVal;
}

int **AnovaAnalyzer::get_code()
{
   int** retVal = NULL, ii;
   if (code_)
   {
      retVal = new int*[tableLeng_];
      checkAllocate(retVal,"retVal in AnovaAnalyzer::get_code");
      for (ii = 0; ii < tableLeng_; ii++)
      {
    	  retVal[ii] = new int[3];
          checkAllocate(retVal[ii],"retVal in AnovaAnalyzer::get_code");
    	  std::copy(code_[ii], code_[ii]+3, retVal[ii]);
      }
   }
   return retVal;
}

