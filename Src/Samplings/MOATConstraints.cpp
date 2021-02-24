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
// Functions for the class MOATConstraint
// AUTHOR : CHARLES TONG
// DATE   : 2009
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "PsuadeData.h"
#include "MOATConstraints.h"
#include "PsuadeUtil.h"
#define PABS(X) (((X) > 0)? X : -(X))

// ************************************************************************
// Constructor 
// ------------------------------------------------------------------------
MOATConstraints::MOATConstraints()
{
   nConstraints_ = 0;
   constraintFAs_ = NULL;
   constraintNInputs_ = NULL;
   constraintInputIndices_ = NULL;
   constraintInputValues_ = NULL;
   nInputs_ = 0;
   XLBounds_ = NULL;
   XUBounds_ = NULL;
   YLBounds_ = NULL;
   YUBounds_ = NULL;
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
MOATConstraints::MOATConstraints(const MOATConstraints & mc)
{
  nConstraints_ = mc.nConstraints_;
  constraintFAs_ = new FuncApprox*[nConstraints_];
  constraintNInputs_ = new int[nConstraints_];
  constraintInputIndices_ = new int*[nConstraints_];
  constraintInputValues_ = new double*[nConstraints_];
  nInputs_ = mc.nInputs_;
  sizeXLBounds_ = mc.sizeXLBounds_;
  sizeXUBounds_ = mc.sizeXUBounds_;
  XLBounds_ = new double[sizeXLBounds_];
  XUBounds_ = new double[sizeXUBounds_];
  YLBounds_ = new double[nConstraints_];
  YUBounds_ = new double[nConstraints_];
  for(int i = 0; i < nConstraints_; i++) {
    constraintInputIndices_[i] = mc.constraintInputIndices_[i];
    constraintInputValues_[i] = mc.constraintInputValues_[i];
    constraintFAs_[i] = mc.constraintFAs_[i];
    constraintNInputs_[i] = mc.constraintNInputs_[i];
    XLBounds_[i] = mc.XLBounds_[i];
    XUBounds_[i] = mc.XUBounds_[i];
    for(int j = 0; j < constraintNInputs_[i]; j++) {
      constraintFAs_[i][j] = mc.constraintFAs_[i][j];
      constraintInputIndices_[i][j] = mc.constraintInputIndices_[i][j];
      constraintInputValues_[i][j] = mc.constraintInputValues_[i][j];
      
    }
  }
  for(int i = 0; i < sizeXLBounds_; i++)
    XLBounds_[i] = mc.XLBounds_[i];
  for(int i = 0; i < sizeXUBounds_; i++)
    XUBounds_[i] = mc.XUBounds_[i];

} 
   
// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
MOATConstraints::~MOATConstraints()
{
   int ii;

   if (constraintFAs_ != NULL)
   {
      for (ii = 0; ii < nConstraints_; ii++)
         if (constraintFAs_[ii] != NULL) delete constraintFAs_[ii];
      delete [] constraintFAs_;
   }
   if (constraintNInputs_ != NULL) delete [] constraintNInputs_;
   if (constraintInputIndices_ != NULL)
   {
      for (ii = 0; ii < nConstraints_; ii++)
         if (constraintInputIndices_[ii] != NULL)
            delete [] constraintInputIndices_[ii];
      delete [] constraintInputIndices_;
   }
   if (constraintInputValues_ != NULL)
   {
      for (ii = 0; ii < nConstraints_; ii++)
         if (constraintInputValues_[ii] != NULL)
            delete [] constraintInputValues_[ii];
      delete [] constraintInputValues_;
   }
   if (XLBounds_ != NULL) delete [] XLBounds_;
   if (XUBounds_ != NULL) delete [] XUBounds_;
   if (YLBounds_ != NULL) delete [] YLBounds_;
   if (YUBounds_ != NULL) delete [] YUBounds_;
}

// ************************************************************************
// generate all MOAT  constraints
// ------------------------------------------------------------------------
int MOATConstraints::initialize(PsuadeData *psuadeIO)
{
   int        printLevel, ii, jj, status, *sortArray, nInputsChk;
   double     *filterLBounds, *filterUBounds;
   char       **filterIndexFiles, **filterDataFiles;
   FILE       *fp;
   pData      pPtr, flbPtr, fubPtr, fdfPtr, fifPtr, iLPtr, iUPtr;
   PsuadeData *pIO;

   psuadeIO->getParameter("input_ninputs", pPtr);
   nInputs_ = pPtr.intData_;
   psuadeIO->getParameter("input_lbounds", iLPtr);
   XLBounds_ = iLPtr.dbleArray_;
   sizeXLBounds_ = iLPtr.nDbles_;
   iLPtr.dbleArray_ = NULL;
   psuadeIO->getParameter("input_ubounds", iUPtr);
   XUBounds_ = iUPtr.dbleArray_;
   sizeXUBounds_ = iUPtr.nDbles_;
   iUPtr.dbleArray_ = NULL;
   psuadeIO->getParameter("ana_diagnostics", pPtr);
   printLevel = pPtr.intData_;
   psuadeIO->getParameter("ana_num_moatfilters", pPtr);
   nConstraints_ = pPtr.intData_;
   if (nConstraints_ > 0)
   {
      if (printLevel > 0)
         printf("MOATConstraints: number of filters = %d\n",nConstraints_);
      psuadeIO->getParameter("ana_moatfilterlbounds", flbPtr);
      filterLBounds = flbPtr.dbleArray_;
      YLBounds_ = new double[nConstraints_];
      for (ii = 0; ii < nConstraints_; ii++) YLBounds_[ii] = filterLBounds[ii];
      psuadeIO->getParameter("ana_moatfilterubounds", fubPtr);
      filterUBounds = fubPtr.dbleArray_;
      YUBounds_ = new double[nConstraints_];
      for (ii = 0; ii < nConstraints_; ii++) YUBounds_[ii] = filterUBounds[ii];
      psuadeIO->getParameter("ana_moatfilterdatafile", fdfPtr);
      filterDataFiles = fdfPtr.strArray_;
      psuadeIO->getParameter("ana_moatfilterindexfile", fifPtr);
      filterIndexFiles = fifPtr.strArray_;
      for (ii = 0; ii < nConstraints_; ii++)
      {
         if (strcmp(filterDataFiles[ii], "NULL"))
         {
            fp = fopen(filterDataFiles[ii], "r");
            if (fp == NULL)
            {
               printf("MOATConstraints ERROR: filter data file %s not found.\n",
                      filterDataFiles[ii]);
               exit(1);
            }
            else fclose(fp);
         }
      }

      constraintFAs_ = new FuncApprox*[nConstraints_];
      constraintNInputs_ = new int[nConstraints_];
      constraintInputIndices_ = new int*[nConstraints_];
      constraintInputValues_ = new double*[nConstraints_];
      for (ii = 0; ii < nConstraints_; ii++)
      {
         if (printLevel > 0)
            printf("MOATConstraints: initializing filter %d\n",ii+1);
         constraintFAs_[ii] = NULL;
         constraintNInputs_[ii] = 0;
         constraintInputIndices_[ii] = NULL;
         constraintInputValues_[ii] = NULL;
         if (strcmp(filterDataFiles[ii], "NULL"))
         {
            pIO = new PsuadeData();
            status = pIO->readPsuadeFile(filterDataFiles[ii]);
            if (status != 0) 
            {
               printf("MOATConstraints ERROR: cannot read filter file %s.\n",
                      filterDataFiles[ii]);
               exit(1);
            }
            constraintFAs_[ii] = genFAInteractive(pIO, 2);
            pIO->getParameter("input_ninputs", pPtr);
            nInputsChk = pPtr.intData_;
            delete pIO;

            if (!strcmp(filterIndexFiles[ii],"NULL"))
            {
               printf("MOATConstraints ERROR: filter index file not given.\n");
               exit(1);
            }
            if (printLevel > 0)
               printf("MOATConstraints: filter %d has index file = %s\n",
                      ii+1, filterIndexFiles[ii]);
            fp = fopen(filterIndexFiles[ii],"r");
            if (fp != NULL)
            {
               fscanf(fp, "%d", &constraintNInputs_[ii]);
               if (constraintNInputs_[ii] != nInputsChk)
               {
                  printf("MOATConstraints: filter %d must have %d(%d) inputs\n",
                         ii+1, nInputsChk, constraintNInputs_[ii]);
                  exit(1);
               }
               constraintInputIndices_[ii] = new int[constraintNInputs_[ii]];
               constraintInputValues_[ii] = new double[constraintNInputs_[ii]];
               for (jj = 0; jj < constraintNInputs_[ii]; jj++)
               {
                  fscanf(fp,"%d %lg",&constraintInputIndices_[ii][jj],
                         &(constraintInputValues_[ii][jj]));
                  if (printLevel > 0)
                     printf("MOATConstraints: filter %d has input %d = %d\n",
                               ii+1, jj+1, constraintInputIndices_[ii][jj]);
                  if (constraintInputIndices_[ii][jj] >= 0)
                     constraintInputIndices_[ii][jj]--;
                  if (constraintInputIndices_[ii][jj] >= nInputs_ ||
                      constraintInputIndices_[ii][jj] < -1)
                  {
                     printf("MOATConstraints ERROR: invalid filter input %d.\n",
                            constraintInputIndices_[ii][jj]);
                     printf("MOATConstraints input format: ");
                     printf("nInputs \n ");
                     printf("<input (or 0)> <value (nominal val if 0)> \n ");
                     printf("<input (or 0)> <value (nominal val if 0)> \n ");
                     printf("... \n ");
                     exit(1);
                  }
               }
               sortArray = new int[constraintNInputs_[ii]];
               for (jj = 0; jj < constraintNInputs_[ii]; jj++)
                  sortArray[jj] = constraintInputIndices_[ii][jj];
               sortIntList(constraintNInputs_[ii], sortArray);
               for (jj = 1; jj < constraintNInputs_[ii]; jj++)
               {
                  if (sortArray[jj] >= 0 && sortArray[jj] == sortArray[jj-1])
                  {
                     printf("MOATConstraints ERROR: filter has duplicate");
                     printf(" input %d.\n", sortArray[jj]);
                     printf("MOATConstraints input format: ");
                     printf("nInputs \n ");
                     printf("<input (or 0)> <value (nominal val if 0)> \n ");
                     printf("<input (or 0)> <value (nominal val if 0)> \n ");
                     printf("... \n ");
                     exit(1);
                  }
               }
               delete [] sortArray;
               fclose(fp);
            }
            else
            {
               printf("MOATConstraints ERROR : filter index file %s not found.\n",
                      filterIndexFiles[ii]);
               exit(1);
            }
         }
      }
   }
   return 0;
}

// ************************************************************************
// compute the MOAT scale
// ------------------------------------------------------------------------
double MOATConstraints::getScale(double *sampleInputs, int diffIndex, 
                                 int &flag)
{
   int    ii, jj, searchIndex;
   double *filterSamplePt, filterRange, *XLs, *XUs, XL, XU, scale;
   double currX1, currX2, currY1, currY2, sLow, sHi;

   if (nConstraints_ <= 0) 
   {
      flag = 1;
      return 0.0;
   }
   XLs = new double[nConstraints_];
   XUs = new double[nConstraints_];
   for (ii = 0; ii < nConstraints_; ii++)
   {
      filterSamplePt = new double[constraintNInputs_[ii]];
      for (jj = 0; jj < constraintNInputs_[ii]; jj++)
         if (constraintInputIndices_[ii][jj] == -1)
            filterSamplePt[jj] = constraintInputValues_[ii][jj];
      filterRange = YUBounds_[ii] - YLBounds_[ii];
      searchIndex = -1;
      for (jj = 0; jj < constraintNInputs_[ii]; jj++)
      {
         if (constraintInputIndices_[ii][jj] == diffIndex)
         {
            searchIndex = jj;
            break;
         }
      }
      if (searchIndex < 0) 
      {
         XLs[ii] = XLBounds_[diffIndex];
         XUs[ii] = XUBounds_[diffIndex];
      }
      else
      {
         for (jj = 0; jj < constraintNInputs_[ii]; jj++)
            if (constraintInputIndices_[ii][jj] != -1)
               filterSamplePt[jj]=sampleInputs[constraintInputIndices_[ii][jj]];
         filterSamplePt[searchIndex] = XLBounds_[diffIndex];
         currY1 = constraintFAs_[ii]->evaluatePoint(filterSamplePt);
         filterSamplePt[searchIndex] = XUBounds_[diffIndex];
         currY2 = constraintFAs_[ii]->evaluatePoint(filterSamplePt);
         currX1 = XLBounds_[diffIndex];
         currX2 = XUBounds_[diffIndex];
         if (currY2 >= YUBounds_[ii] && currY1 >= YUBounds_[ii])
            currX1 = currX2 = 0.0;
         else if (currY2 <= YLBounds_[ii] && currY1 <= YLBounds_[ii])
            currX1 = currX2 = 0.0;
         else if (currY2 > currY1)
         {
            if (currY2 <= YUBounds_[ii]) currX2 = XUBounds_[diffIndex];
            else
            {
               sLow = XLBounds_[diffIndex];
               sHi  = XUBounds_[diffIndex];
               while (PABS((currY2-YUBounds_[ii])/filterRange)>0.001)
               {
                  filterSamplePt[searchIndex] = 0.5 * (sLow + sHi);
                  currY2 = constraintFAs_[ii]->evaluatePoint(filterSamplePt);
                  if (currY2 > YUBounds_[ii])
                     sHi  = 0.5 * (sLow + sHi);
                  else
                     sLow = 0.5 * (sLow + sHi);
               }
               currX2 = filterSamplePt[searchIndex];
            }
            if (currY1 >= YLBounds_[ii]) currX1 = XLBounds_[diffIndex];
            else
            {
               sLow = XLBounds_[diffIndex];
               sHi  = XUBounds_[diffIndex];
               while (PABS((currY1-YLBounds_[ii])/filterRange)>0.001)
               {
                  filterSamplePt[searchIndex] = 0.5 * (sLow + sHi);
                  currY1 = constraintFAs_[ii]->evaluatePoint(filterSamplePt);
                  if (currY1 < YLBounds_[ii])
                     sLow = 0.5 * (sLow + sHi);
                  else
                     sHi  = 0.5 * (sLow + sHi);
               }
               currX1 = filterSamplePt[searchIndex];
            }
         }
         else
         {
            if (currY1 <= YUBounds_[ii]) currX1 = XLBounds_[diffIndex];
            else
            {
               sLow = XLBounds_[diffIndex];
               sHi  = XUBounds_[diffIndex];
               while (PABS((currY1-YUBounds_[ii])/filterRange)>0.001)
               {
                  filterSamplePt[searchIndex] = 0.5 * (sLow + sHi);
                  currY1 = constraintFAs_[ii]->evaluatePoint(filterSamplePt);
                  if (currY1 > YUBounds_[ii])
                     sLow = 0.5 * (sLow + sHi);
                  else
                     sHi  = 0.5 * (sLow + sHi);
               }
               currX1 = filterSamplePt[searchIndex];
            }
            if (currY2 >= YLBounds_[ii]) currX2 = XUBounds_[diffIndex];
            else
            {
               sLow = XLBounds_[diffIndex];
               sHi  = XUBounds_[diffIndex];
               while (PABS((currY2-YLBounds_[ii])/filterRange)>0.001)
               {
                  filterSamplePt[searchIndex] = 0.5 * (sLow + sHi);
                  currY2 = constraintFAs_[ii]->evaluatePoint(filterSamplePt);
                  if (currY2 < YLBounds_[ii])
                     sHi  = 0.5 * (sLow + sHi);
                  else
                     sLow = 0.5 * (sLow + sHi);
               }   
               currX2 = filterSamplePt[searchIndex];
            }
         }
         if (currX1 < currX2)
         {
            XLs[ii] = currX1;
            XUs[ii] = currX2;
         }
         else
         {
            XLs[ii] = currX2;
            XUs[ii] = currX1;
         }
      }
      delete [] filterSamplePt;
   }
   XU = XUs[0];
   for (ii = 1; ii < nConstraints_; ii++)
      if (XUs[ii] < XU) XU = XUs[ii];
   XL = XLs[0];
   for (ii = 1; ii < nConstraints_; ii++)
      if (XLs[ii] > XL) XL = XLs[ii];
   scale = XU - XL;
   delete [] XLs;
   delete [] XUs;
   flag = 0;
   if (scale < 0.0) scale = 0.0;
   return scale;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
MOATConstraints& MOATConstraints::operator=(const MOATConstraints &mc)
{
   int ii;
   if(this == &mc) return (*this);

   if (constraintFAs_ != NULL)
   {
      for (ii = 0; ii < nConstraints_; ii++)
         if (constraintFAs_[ii] != NULL) delete constraintFAs_[ii];
      delete [] constraintFAs_;
   }
   if (constraintNInputs_ != NULL) delete [] constraintNInputs_;
   if (constraintInputIndices_ != NULL)
   {
      for (ii = 0; ii < nConstraints_; ii++)
         if (constraintInputIndices_[ii] != NULL)
            delete [] constraintInputIndices_[ii];
      delete [] constraintInputIndices_;
   }
   if (constraintInputValues_ != NULL)
   {
      for (ii = 0; ii < nConstraints_; ii++)
         if (constraintInputValues_[ii] != NULL)
            delete [] constraintInputValues_[ii];
      delete [] constraintInputValues_;
   }
   if (XLBounds_ != NULL) delete [] XLBounds_;
   if (XUBounds_ != NULL) delete [] XUBounds_;
   if (YLBounds_ != NULL) delete [] YLBounds_;
   if (YUBounds_ != NULL) delete [] YUBounds_;

   nConstraints_ = mc.nConstraints_;
   constraintFAs_ = new FuncApprox*[nConstraints_];
   constraintNInputs_ = new int[nConstraints_];
   constraintInputIndices_ = new int*[nConstraints_];
   constraintInputValues_ = new double*[nConstraints_];
   nInputs_ = mc.nInputs_;
   sizeXLBounds_ = mc.sizeXLBounds_;
   sizeXUBounds_ = mc.sizeXUBounds_;
   XLBounds_ = new double[sizeXLBounds_];
   XUBounds_ = new double[sizeXUBounds_];
   YLBounds_ = new double[nConstraints_];
   YUBounds_ = new double[nConstraints_];
   for(int i = 0; i < nConstraints_; i++) 
   {
      constraintInputIndices_[i] = mc.constraintInputIndices_[i];
      constraintInputValues_[i] = mc.constraintInputValues_[i];
      constraintFAs_[i] = mc.constraintFAs_[i];
      constraintNInputs_[i] = mc.constraintNInputs_[i];
      XLBounds_[i] = mc.XLBounds_[i];
      XUBounds_[i] = mc.XUBounds_[i];
      for(int j = 0; j < constraintNInputs_[i]; j++)
      {
         constraintFAs_[i][j] = mc.constraintFAs_[i][j];
         constraintInputIndices_[i][j] = mc.constraintInputIndices_[i][j];
         constraintInputValues_[i][j] = mc.constraintInputValues_[i][j];
      }
   }
   for(int i = 0; i < sizeXLBounds_; i++)
      XLBounds_[i] = mc.XLBounds_[i];
   for(int i = 0; i < sizeXUBounds_; i++)
      XUBounds_[i] = mc.XUBounds_[i];
   return (*this);
}

