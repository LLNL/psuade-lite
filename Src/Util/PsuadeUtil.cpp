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
// Utility functions 
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
using namespace std;
#include "Psuade.h"
#include "PsuadeUtil.h"
#include "PrintingTS.h"
#include "Globals.h"
#include "Matrix.h"

// ------------------------------------------------------------------------
// external functions 
// ------------------------------------------------------------------------
extern "C" {
   void dgesvd_(char *, char *, int *, int *, double *, int *, double *,
               double *, int *, double *, int *, double *, int *, int *);
}

// ------------------------------------------------------------------------
// local defines
// ------------------------------------------------------------------------
static double  elapsedTime_;
static randctx rctx_; 
static int     randFlag_=0;
static int     randMask_;
static double  pgplotXmin_, pgplotXmax_, pgplotYmin_, pgplotYmax_;
static long    pgplotFlag_; 

// ************************************************************************
// random number generator initializer 
// ------------------------------------------------------------------------
void PSUADE_randInit(long seed)
{
   int    ntime, ii, nbits;
   double dtime;

   if (seed == -1)
   {
      for (ii = 0; ii < RANDSIZ; ii++) 
      {
         dtime = getClock();
         ntime = (int) (dtime * 1.0E6);
         ntime = ntime % 10000;
         dtime = getClock();
         ntime = ntime * 10000 + (int) (dtime * 1.0E6);
         rctx_.randrsl[ii] = ntime;
         dtime = getClock();
         ntime = (int) (dtime * 1.0E6);
         ntime = ntime % 10000;
         dtime = getClock();
         ntime = ntime * 10000 + (int) (dtime * 1.0E6);
         rctx_.randmem[ii] = ntime;
      }
      randinit(&rctx_, TRUE);
   }
   else
   {
      psRandomSeed_ = seed;
      randinitBySeed(&rctx_, psRandomSeed_);
   }
   nbits = sizeof(int) * 8;  
   randMask_ = (1 << (nbits-2)) + ((1 << (nbits-2)) - 1);
   randFlag_ = 1;
}

// ************************************************************************
// random number generator 
// ------------------------------------------------------------------------
int PSUADE_rand()
{
   int  irand;

   if (randFlag_ == 0) PSUADE_randInit(psRandomSeed_);

   irand = ISAAC_RAND_RAND(&rctx_);
   irand = irand & randMask_;
      
   return irand;
}

// ************************************************************************
// random number generator (double) 
// ------------------------------------------------------------------------
double PSUADE_drand()
{
   ub4 irand;
   double drand;

   irand = PSUADE_rand();
   drand = ((double) irand / (double) randMask_);
   if (drand > 1 - 1.2e-7) drand = 1.0 - 1.2e-7;
   return drand;
}

// ************************************************************************
// copies a string into a new malloc'd pointer
// ************************************************************************
char* PSUADE_strdup(const char *s)
{
  size_t len = strlen (s) + 1;
  char *result = (char*) malloc (len);
  if (result == (char*) 0) return (char*) 0;
  return (char*) memcpy (result, s, len);
}

// ************************************************************************
// copies a n characters of a string into a new malloc'd pointer
// ************************************************************************
char* PSUADE_strndup(const char *s, size_t n)
{
  char *result;
  size_t len = strlen (s);

  if (n < len) len = n;
  result = (char *) malloc (len + 1);
  if (!result) return 0;
  result[len] = '\0';
  return (char *) memcpy (result, s, len);
}

// ************************************************************************
// generate an array of random number using standard normal distribution
// ------------------------------------------------------------------------
void PSUADE_drandn(int length, double *rArray)
{    
   int    ii, newLength;
   double R, *rData, pi, cTheta;

   newLength = length;
   pi = 4.0 * atan(1.0e0);
   if (length % 2 == 1) newLength++;
   rData = new double[newLength];
   for (ii = 0; ii < newLength; ii++) rArray[ii] = PSUADE_drand();
   for (ii = 0; ii < newLength; ii=+2)
   {
      R = -2.0 * log(rArray[ii]);
      if (R < 0.0) R = 0.0;
      R = sqrt(R);
      cTheta = 2.0 * pi * rArray[ii+1];
      rArray[ii] = R * cos(cTheta);
      if ((ii+1) < length) rArray[ii+1] = R * sin(cTheta);
   }
   delete [] rData;
}

// ************************************************************************
// sortIntList : sort an integer list
// ------------------------------------------------------------------------
void sortIntList(int length, int *intList)
{
   int i, last, mid, idata;

   if ( length <= 1 ) return;
   mid = ( length - 1 ) / 2;
   idata = intList[0];
   intList[0] = intList[mid];
   intList[mid] = idata;
   last = 0;
   for (i = 1; i < length; i++)
   {
      if (intList[i] < intList[0])
      {
         last++;
         idata = intList[last];
         intList[last] = intList[i];
         intList[i] = idata;
      }
   }
   idata = intList[0];
   intList[0] = intList[last];
   intList[last] = idata;
   sortIntList(last,intList);
   sortIntList(length-last-1,&(intList[last+1]));
}

// ************************************************************************
// sort integer list with an integer and a double array 
// ------------------------------------------------------------------------
void sortIntList2a(int length, int *intList, double *valueList)
{
   int    mid, i, idata, last, isort, *intList2, *intList3;
   double ddata, *valueList2, *valueList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   idata = intList[0];
   intList[0] = intList[mid];
   intList[mid] = idata;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   last = 0;
   isort = intList[0];
   intList2 = &(intList[last]);
   valueList2 = &(valueList[last]);
   intList3 = &(intList[1]);
   valueList3 = &(valueList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *intList3 < isort )
      {
         last++;
         intList2++; valueList2++;
         idata = *intList2;
         *intList2 = *intList3;
         *intList3 = idata;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
      }
      intList3++; valueList3++;
   }
   idata = intList[0];
   intList[0] = intList[last];
   intList[last] = idata;
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   sortIntList2a(last, intList, valueList);
   sortIntList2a(length-last-1, &(intList[last+1]), &(valueList[last+1]));
}

// ************************************************************************
// sort double list 
// ------------------------------------------------------------------------
void sortDbleList(int length, double *valueList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   valueList3 = &(valueList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
      }
      valueList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   sortDbleList(last, valueList);
   sortDbleList(length-last-1,&(valueList[last+1]));
}

// ************************************************************************
// sort double list with 2 double arrays 
// ------------------------------------------------------------------------
void sortDbleList2(int length, double *valueList, double *dataList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort, *dataList2, *dataList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   ddata = dataList[0];
   dataList[0] = dataList[mid];
   dataList[mid] = ddata;
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   dataList2 = &(dataList[last]);
   valueList3 = &(valueList[1]);
   dataList3 = &(dataList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         dataList2++; valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         ddata = *dataList2;
         *dataList2 = *dataList3;
         *dataList3 = ddata;
      }
      dataList3++; valueList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   ddata = dataList[0];
   dataList[0] = dataList[last];
   dataList[last] = ddata;
   sortDbleList2(last, valueList, dataList);
   sortDbleList2(length-last-1,&(valueList[last+1]),&(dataList[last+1]));
}

// ************************************************************************
// sort double list with 1 double array and 1 int array
// ------------------------------------------------------------------------
void sortDbleList2a(int length, double *valueList, int *dataList)
{
   int    mid, i, last, idata, *dataList2, *dataList3;
   double ddata, *valueList2, *valueList3, dsort;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   idata = dataList[0];
   dataList[0] = dataList[mid];
   dataList[mid] = idata;
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   dataList2 = &(dataList[last]);
   valueList3 = &(valueList[1]);
   dataList3 = &(dataList[1]);
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         dataList2++; valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         idata = *dataList2;
         *dataList2 = *dataList3;
         *dataList3 = idata;
      }
      dataList3++; valueList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   idata = dataList[0];
   dataList[0] = dataList[last];
   dataList[last] = idata;
   sortDbleList2a(last, valueList, dataList);
   sortDbleList2a(length-last-1,&(valueList[last+1]),&(dataList[last+1]));
}

// ************************************************************************
// sort double list with 3 double arrays 
// ------------------------------------------------------------------------
void sortDbleList3(int length, double *valueList, double *dataList,
                   double *auxList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort, *dataList2, *dataList3;
   double *auxList2, *auxList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[mid];
      dataList[mid] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[mid];
      auxList[mid] = ddata;
   }
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   valueList3 = &(valueList[1]);
   if (dataList != NULL)
   {
      dataList2 = &(dataList[last]);
      dataList3 = &(dataList[1]);
   }
   if (auxList != NULL)
   {
      auxList2 = &(auxList[last]);
      auxList3 = &(auxList[1]);
   }
   for (i = 1; i < length; i++)
   {
      if ( *valueList3 < dsort )
      {
         last++;
         valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         if (dataList != NULL)
         {
            dataList2++;
            ddata = *dataList2;
            *dataList2 = *dataList3;
            *dataList3 = ddata;
         }
         if (auxList != NULL)
         {
            auxList2++;
            ddata = *auxList2;
            *auxList2 = *auxList3;
            *auxList3 = ddata;
         }
      }
      valueList3++; 
      if (dataList != NULL) dataList3++; 
      if (auxList  != NULL) auxList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[last];
      dataList[last] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[last];
      auxList[last] = ddata;
   }
   sortDbleList3(last, valueList, dataList, auxList);
   if (dataList != NULL && auxList != NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    &(auxList[last+1]));
   else if (dataList != NULL && auxList == NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    NULL);
   else if (dataList == NULL && auxList != NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),NULL,
                    &(auxList[last+1]));
   else if (dataList == NULL && auxList == NULL)
      sortDbleList3(length-last-1,&(valueList[last+1]),NULL,NULL);
}

// ************************************************************************
// sort double list with 4 double arrays 
// ------------------------------------------------------------------------
void sortDbleList4(int length, double *valueList, double *dataList,
                   double *auxList, double *moreList)
{
   int    mid, i, last;
   double ddata, *valueList2, *valueList3, dsort, *dataList2, *dataList3;
   double *auxList2, *auxList3, *moreList2, *moreList3;

   if ( length <= 1 ) return;
   mid   = (length - 1) / 2;
   ddata = valueList[0];
   valueList[0] = valueList[mid];
   valueList[mid] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[mid];
      dataList[mid] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[mid];
      auxList[mid] = ddata;
   }
   if (moreList != NULL)
   {
      ddata = moreList[0];
      moreList[0] = moreList[mid];
      moreList[mid] = ddata;
   }
   last = 0;

   dsort = valueList[0];
   valueList2 = &(valueList[last]);
   valueList3 = &(valueList[1]);
   if (dataList != NULL)
   {
      dataList2 = &(dataList[last]);
      dataList3 = &(dataList[1]);
   }
   if (auxList != NULL)
   {
      auxList2 = &(auxList[last]);
      auxList3 = &(auxList[1]);
   }
   if (moreList != NULL)
   {
      moreList2 = &(moreList[last]);
      moreList3 = &(moreList[1]);
   }
   for (i = 1; i < length; i++)
   {
      if (*valueList3 < dsort)
      {
         last++;
         valueList2++;
         ddata = *valueList2;
         *valueList2 = *valueList3;
         *valueList3 = ddata;
         if (dataList != NULL)
         {
            dataList2++;
            ddata = *dataList2;
            *dataList2 = *dataList3;
            *dataList3 = ddata;
         }
         if (auxList != NULL)
         {
            auxList2++;
            ddata = *auxList2;
            *auxList2 = *auxList3;
            *auxList3 = ddata;
         }
         if (moreList != NULL)
         {
            moreList2++;
            ddata = *moreList2;
            *moreList2 = *moreList3;
            *moreList3 = ddata;
         }
      }
      valueList3++;
      if (dataList != NULL) dataList3++;
      if (auxList  != NULL) auxList3++;
      if (moreList != NULL) moreList3++;
   }
   ddata = valueList[0];
   valueList[0] = valueList[last];
   valueList[last] = ddata;
   if (dataList != NULL)
   {
      ddata = dataList[0];
      dataList[0] = dataList[last];
      dataList[last] = ddata;
   }
   if (auxList != NULL)
   {
      ddata = auxList[0];
      auxList[0] = auxList[last];
      auxList[last] = ddata;
   }
   if (moreList != NULL)
   {
      ddata = moreList[0];
      moreList[0] = moreList[last];
      moreList[last] = ddata;
   }
   sortDbleList4(last, valueList, dataList, auxList, moreList);
   if (dataList != NULL && auxList != NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    &(auxList[last+1]), &(moreList[last+1]));
   else if (dataList != NULL && auxList != NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    &(auxList[last+1]), NULL);
   else if (dataList != NULL && auxList == NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    NULL, &(moreList[last+1]));
   else if (dataList != NULL && auxList == NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),&(dataList[last+1]),
                    NULL, NULL);
   else if (dataList == NULL && auxList != NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL,
                    &(auxList[last+1]), &(moreList[last+1]));
   else if (dataList == NULL && auxList != NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL,
                    &(auxList[last+1]), NULL);
   else if (dataList == NULL && auxList == NULL && moreList != NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL, NULL,
                    &(moreList[last+1]));
   else if (dataList == NULL && auxList == NULL && moreList == NULL)
      sortDbleList4(length-last-1,&(valueList[last+1]),NULL,NULL, NULL);
}

// ************************************************************************
// special sort and delete (very specialized for sparse grid)
// sort double list with 4 double arrays 
// ------------------------------------------------------------------------
int sortAndDelete(int nrows, int ncols, double **ddata)
{
   int    *isortList, ii, jj, kk, ll, ibegin, iclose, nrows2;
   double *dsortList;
   
   // range checking by Bill Oliver
   if (nrows <= 0)
   { 
      printf("First parameter to function sortAndDelete is <= 0 in file %s \
               and must be >0 returning -1", __FILE__);
      return -1;
   }
   if(ncols <= 0)
   { 
      printf("Second parameter to function sortAndDelete is <= 0 in file %s \
               and must be >0 returning -1", __FILE__);
      return -1;
   }
             
   dsortList = (double *) malloc(nrows * sizeof(double));
   isortList = (int *) malloc(nrows * sizeof(int));

   // Bill Oliver check for out of memory
   if(dsortList == NULL || isortList == NULL){
     printf("out of memory in file %s line %d aborting\n",__FILE__,__LINE__);
     abort();
   }

   for (ii = 0; ii < ncols-1; ii++)
   {
      if (ii == 0)
      {
         for (jj = 0; jj < nrows; jj++)
         {
            dsortList[jj] = ddata[jj][ii];
            isortList[jj] = jj;
         }
         sortDbleList2a(nrows, dsortList, isortList);
         for (jj = 0; jj < nrows; jj++)
            ddata[jj][ii] = dsortList[jj] ;
         for (kk = ii+1; kk < ncols; kk++)
         {
            for (jj = 0; jj < nrows; jj++)
               dsortList[jj] = ddata[isortList[jj]][kk];
            for (jj = 0; jj < nrows; jj++) ddata[jj][kk] = dsortList[jj];
         }
      }
      else
      {
         ibegin = 0;
         while (ibegin < nrows)
         {
           iclose = ibegin + 1;
           for (jj = iclose; jj < nrows; jj++)
           {
              if (ddata[jj][ii-1] != ddata[jj-1][ii-1]) break;
              else
              {
                 for (kk = 0; kk < ii-1; kk++)
                    if (ddata[jj][kk] != ddata[jj-1][kk]) break;
                 if (kk < ii-1) break;
              }
           }
           iclose = jj;
           if (iclose - ibegin > 1)
           {
              for (jj = ibegin; jj < iclose; jj++)
              {
                 dsortList[jj-ibegin] = ddata[jj][ii];
                 isortList[jj-ibegin] = jj;
              }
              sortDbleList2a(iclose-ibegin, dsortList, isortList);
              for (jj = ibegin; jj < iclose; jj++)
                 ddata[jj][ii] = dsortList[jj-ibegin];
              for (kk = ii+1; kk < ncols; kk++)
              {
                 for (jj = 0; jj < iclose-ibegin; jj++)
                    dsortList[jj] = ddata[isortList[jj]][kk];
                 for (jj = 0; jj < iclose-ibegin; jj++)
                     ddata[jj+ibegin][kk] = dsortList[jj];
              }
            }
            ibegin = iclose;
         }
      }
   }
   free(dsortList);
   free(isortList);
   /* delete */
   kk = 1;
   nrows2 = nrows;
   while (kk < nrows2)
   {
      for (jj = 0; jj < ncols-1; jj++)
         if (ddata[kk-1][jj] != ddata[kk][jj]) break;
      if (jj == ncols-1)
      {
         ddata[kk-1][ncols-1] += ddata[kk][ncols-1];
         for (jj = kk+1; jj < nrows2; jj++)
            for (ll = 0; ll < ncols; ll++) ddata[jj-1][ll] = ddata[jj][ll];
         nrows2--;
      }
      else kk++;
   }
   return nrows2;
}

// ************************************************************************
// sort the rows of a matrix (specially made for ProbMatris)
// ------------------------------------------------------------------------
int sortnDbleList(int nrows, int ncols, double **dbleData, int *intData)
{
  int    *isortList, ii, jj, kk, ll, ibegin, iclose, nrows2, *isortList2;
  double *dsortList;
   
  if (nrows <= 0)
  { 
    printf("First parameter to function sortnDbleList  is <= 0 in file %s \
           and must be >0 returning -1", __FILE__);
    return -1;
  }
  if (ncols <= 0)
  { 
    printf("Second parameter to function sortnDbleList  is <= 0 in file %s \
           and must be >0 returning -1", __FILE__);
    return -1;
  }
  dsortList = new double[nrows];
  isortList = new int[nrows];
  isortList2 = new int[nrows];

  if (dsortList == NULL || isortList == NULL)
  {
    printf("out of memory in file %s line %d aborting\n",__FILE__,__LINE__);
    abort();
  }

  ii = 0;
  for (jj = 0; jj < nrows; jj++)
  {
    dsortList[jj] = dbleData[jj][ii];
    isortList[jj] = jj;
  }
  sortDbleList2a(nrows, dsortList, isortList);
  for (jj = 0; jj < nrows; jj++)
    dbleData[jj][ii] = dsortList[jj] ;
  for (kk = ii+1; kk < ncols; kk++)
  {
    for (jj = 0; jj < nrows; jj++)
      dsortList[jj] = dbleData[isortList[jj]][kk];
    for (jj = 0; jj < nrows; jj++) dbleData[jj][kk] = dsortList[jj];
  }
  for (jj = 0; jj < nrows; jj++)
    isortList2[jj] = intData[isortList[jj]];
  for (jj = 0; jj < nrows; jj++)
    intData[jj] = isortList2[jj] ;

  for (ii = 1; ii < ncols-1; ii++)
  {
    ibegin = 0;
    while (ibegin < nrows)
    {
      iclose = ibegin + 1;
      for (jj = iclose; jj < nrows; jj++)
      {
        if (dbleData[jj][ii-1] != dbleData[jj-1][ii-1]) break;
        else
        {
          for (kk = 0; kk < ii-1; kk++)
            if (dbleData[jj][kk] != dbleData[jj-1][kk]) break;
          if (kk < ii-1) break;
        }
      }
      iclose = jj;
      if (iclose - ibegin > 1)
      {
        for (jj = ibegin; jj < iclose; jj++)
        {
          dsortList[jj-ibegin] = dbleData[jj][ii];
          isortList[jj-ibegin] = jj;
        }
        sortDbleList2a(iclose-ibegin, dsortList, isortList);
        for (jj = ibegin; jj < iclose; jj++)
          dbleData[jj][ii] = dsortList[jj-ibegin];
        for (kk = ii+1; kk < ncols; kk++)
        {
          for (jj = 0; jj < iclose-ibegin; jj++)
             dsortList[jj] = dbleData[isortList[jj]][kk];
          for (jj = 0; jj < iclose-ibegin; jj++)
              dbleData[jj+ibegin][kk] = dsortList[jj];
        }
        for (jj = 0; jj < iclose-ibegin; jj++)
          isortList2[jj] = intData[isortList[jj]];
        for (jj = 0; jj < iclose-ibegin; jj++)
          intData[jj+ibegin] = isortList2[jj];
      }
      ibegin = iclose;
    }
  }
  delete [] dsortList;
  delete [] isortList;
  delete [] isortList2;
  return 0;
}

// ************************************************************************
// search
// ------------------------------------------------------------------------
int binarySearchInt(int sKey, int *sData, int length)
{
   int  left, right, mid, found, index;

   if (length <= 0) return -1;
   left  = 0;
   right = length - 1;
   if (sKey > sData[right]) return -(right+1);
   if (sKey < sData[left])  return -(left+1);
   found = 0;
   while ((found == 0) && ((right-left)>1))
   {
      mid = (left + right) / 2;
      if      (sKey == sData[mid]) {index  = mid; found = 1;}
      else if (sKey >  sData[mid])  left  = mid;
      else                          right = mid;
   }
   if (found == 1)                return index;
   else if (sKey == sData[left])  return left;
   else if (sKey == sData[right]) return right;
   else                           return -(left+1);
}

// ************************************************************************
// output a line of the provided character
// ------------------------------------------------------------------------
void printCharLine(int printLevel, int length, char printchar)
{
  if (length <= 0) length = 70;
  for (int ii = 0; ii < length; ii++)
    printOutTS(printLevel, "%c", printchar);
  printOutTS(printLevel, "\n");
}

// ************************************************************************
// output a line of asterisks
// ------------------------------------------------------------------------
void printAsterisks(int printLevel, int length)
{
  printCharLine(printLevel, length, '*');
}

// ************************************************************************
// output a line of dashes
// ------------------------------------------------------------------------
void printDashes(int printLevel, int length)
{
  printCharLine(printLevel, length, '-');
}

// ************************************************************************
// output a line of equal signs
// ------------------------------------------------------------------------
void printEquals(int printLevel, int length)
{
   printCharLine(printLevel, length, '=');
}

// ************************************************************************
// check null for an allocation
// ------------------------------------------------------------------------
void checkAllocate(void *ptr, const char *mssg)
{
   if (ptr == NULL)
   {
      cerr << "Memory allocation ERROR in: " << mssg << "\n";
      exit(1);
   }
}

// ************************************************************************
// search
// ------------------------------------------------------------------------
int binarySearchDble(double sKey, double *sData, int length)
{
   int  left, right, mid, found, index;

   if (length <= 0) return -1;
   left  = 0;
   right = length - 1;
   if (sKey > sData[right]) return -(right+1);
   if (sKey < sData[left])  return -(left+1);
   found = 0;
   while ((found == 0) && ((right-left)>1))
   {
      mid = (left + right) / 2;
      if      (sKey == sData[mid]) {index  = mid; found = 1;}
      else if (sKey >  sData[mid])  left  = mid;
      else                          right = mid;
   }
   if (found == 1)                return index;
   else if (sKey == sData[left])  return left;
   else if (sKey == sData[right]) return right;
   else                           return -(left+1);
}

// ************************************************************************
// timer functions
// ************************************************************************

// ************************************************************************
// get the time of the day in seconds
// ------------------------------------------------------------------------
double getClock()
{
   double time_i;
   double time_d;
   struct timeval tp;
   struct timezone tzp;
   Gettimeofday(&tp,&tzp);
   time_i = tp.tv_sec % 10;
   time_d = (double) time_i;
   time_i = tp.tv_usec;
   time_d = time_d + (double) time_i / 1000000.0; 
   return(time_d);
}

// ************************************************************************
// startTimer 
// ------------------------------------------------------------------------
void startTimer()
{
   elapsedTime_ = getClock();
}

// ************************************************************************
// StopTimer 
// ------------------------------------------------------------------------
void StopTimer()
{
   elapsedTime_ = getClock() - elapsedTime_;
}

// ************************************************************************
// Get the elapsed time (from startTimer to stopTimer)
// ------------------------------------------------------------------------
double getElapsedTime()
{
   return elapsedTime_;
}

// ************************************************************************
// set screen dump mode
// ------------------------------------------------------------------------
void setScreenDumpMode(int mode)
{
   if (mode == 0) psScreenOutput_ = 0;
   else           psScreenOutput_ = 1;
}

// ************************************************************************
// set screen dump mode
// ------------------------------------------------------------------------
void setLibraryMode(int mode)
{
   if (mode == 0) psLibMode_ = 0;
   else           psLibMode_ = 1;
}

// ************************************************************************
// set interactive mode
// ------------------------------------------------------------------------
void setInteractiveMode(int mode)
{
   if (mode == 0) psInteractive_ = 0;
   else           psInteractive_ = 1;
}

// ************************************************************************
// check screen dump mode
// ------------------------------------------------------------------------
int isScreenDumpModeOn()
{
   if (psScreenOutput_ == 1) return 1;
   else                      return 0;
}

// ************************************************************************
// check library mode
// ------------------------------------------------------------------------
int isLibraryModeOn()
{
   if (psLibMode_ == 1) return 1;
   else                 return 0;
}

// ************************************************************************
// check interactive mode
// ------------------------------------------------------------------------
int isInteractiveModeOn()
{
   if (psInteractive_ == 1) return 1;
   else                     return 0;
}

// ************************************************************************
// plot facilities : call this before calling anything else
// ------------------------------------------------------------------------
void Plotbegin(double xmin, double xmax, double ymin, double ymax)
{
   pgplotFlag_ = 0;
   pgplotXmin_ = xmin;
   pgplotXmax_ = xmax;
   pgplotYmin_ = ymin;
   pgplotYmax_ = ymax;
}

// ************************************************************************
// plot 2D data
// ------------------------------------------------------------------------
void PlotSamples2D(long n, double *x, double *y, long *nstat)
{
#ifdef IRIX
   plot2dirix_(&n, x, y, nstat, &pgplotXmin_, &pgplotXmax_, 
               &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
#else
   plotsamples2d_(&n, x, y, nstat, &pgplotXmin_, &pgplotXmax_, 
                  &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
#endif
}

// ************************************************************************
// 2D scatter plots
// ------------------------------------------------------------------------
void PlotScatter2D(long n, double *x, double *y)
{
   plotscatter2d_(&n, x, y, &pgplotXmin_, &pgplotXmax_, 
                  &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
}

// ************************************************************************
// 2D multiple scatter plots
// ------------------------------------------------------------------------
void PlotScatterM2D(long n, double *x, double *x2, double *y)
{
   plotscatterm2d_(&n, x, x2, y, &pgplotXmin_, &pgplotXmax_, 
                  &pgplotYmin_, &pgplotYmax_, &pgplotFlag_);
}

// ************************************************************************
// plot 3D data
// ------------------------------------------------------------------------
void Plot3d(long n, double *, double *, double *z)
{
   plot3d_(&n, z, &pgplotFlag_);
}

void Plotend()
{
   plotend_();
}

// ************************************************************************
// generate random integer vector
// ------------------------------------------------------------------------
void generateRandomIvector(int leng, int *vec)
{
   int    k, j, irand, *iArray;

   iArray = new int[leng];
   for ( j = 0; j < leng; j++ ) iArray[j] = j; 
   k = leng;
   while (k > 1)
   {
      irand = PSUADE_rand() % k;
      vec[leng-k] = iArray[irand];
      for ( j = irand; j < k-1; j++ ) iArray[j] = iArray[j+1];
      k--; 
   }
   if (k == 1) vec[leng-1] = iArray[0];
   delete [] iArray;
}

// ************************************************************************
// check a number is prime
// ------------------------------------------------------------------------
int checkPrime(int idata)
{
  int  i, upper;

  upper = (int) (sqrt((double) idata) + 1.0);
  if (idata <= 5) return 1;
  else 
  {
     for ( i = 2; i < upper; i++ )
       if (((i & 1) != 0) && (idata / i * i == idata)) return 0;
  }
  return 1;
}

// ************************************************************************
// factorize a number
// ------------------------------------------------------------------------
int *factorize(int iValue)
{
   int  ind, itemp, sqrtN, *factors, nFactors;

   factors = new int[201];
   sqrtN    = (int)(sqrt((double) iValue));
   itemp    = iValue;
   ind      = 2;
   nFactors = 0;

   while (ind <= sqrtN)
   {
      if (itemp % ind == 0)
      {
         if (nFactors > 200)
         {
            printf("Factorize ERROR : factors > 200.\n");
            exit(1);
         } 
         factors[nFactors+1] = ind;
         nFactors++;
         itemp /= ind;
      }
      else ind++;
   }
   factors[0] = nFactors;
   return factors;
}

// ************************************************************************
// compare samples to see if there are duplicates
// ------------------------------------------------------------------------
int compareSamples(int sampleID, int nSamples, int nInputs,
                   double *sampleData, int *sampleStates)
{
   int  sampleBase, indexBase, found=0, ss, ii;

   sampleBase = sampleID * nInputs;
   for (ss = 0; ss < nSamples; ss++)
   {
      if ((ss != sampleID) && (sampleStates[ss] > 0))
      {
         indexBase = ss * nInputs;
         for (ii = 0; ii < nInputs; ii++)
            if (sampleData[sampleBase+ii] != sampleData[indexBase+ii]) break;
         if (ii == nInputs)
         {
            found = 1;
            break;
         }
      }
   } 
   if (found == 0) return -1;
   else            return ss;
}

// ************************************************************************
// compute the number of terms in the permutations
// ------------------------------------------------------------------------
int computeNumPCEPermutations(int nRVs, int pOrder)
{
   int ii, nTerms;

   nTerms = 1;
   for (ii = nRVs; ii < nRVs+pOrder; ii++) nTerms *= (ii + 1);
   for (ii = 0; ii < pOrder; ii++) nTerms /= (ii + 1);
   return nTerms;
}

// ************************************************************************
// compute OUU sample file formats
// ------------------------------------------------------------------------
int checkOUUFileFormat(char *fname, int which, int nInps, int printLevel)
{
   int    intparam1, intparam2, ii, kk, index;
   double value;
   char   lineIn[20001];
   FILE   *fp = fopen(fname, "r");
   if (fp == NULL)
   {
      if (printLevel > 0) printf("ERROR: file %s not readable.\n", fname);
      return -1;
   }
   lineIn[0] = '#';
   fgets(lineIn, 20000, fp);
   while (lineIn[0] == '#') fgets(lineIn, 20000, fp);
   sscanf(lineIn, "%d %d", &intparam1, &intparam2);
   if (intparam1 < nInps+1)
   {
      if (printLevel > 0)
         printf("ERROR: first parameter in line 1 %d not large enough.\n",
                intparam1);
      fclose(fp);
      return -2;
   }
   if (intparam2 != nInps)
   {
      if (printLevel > 0)
         printf("ERROR: second parameter in line 1 %d not equal to %d.\n",
                intparam2, nInps);
      fclose(fp);
      return -3;
   }
   if (which == 1)
   {
      for (ii = 0; ii < intparam1; ii++)
      {
         fgets(lineIn, 20000, fp);
         index = kk = 0;
         while (kk < nInps)
         {
            while (lineIn[index] == ' ') index++;
            if (lineIn[index] == '\0' || lineIn[index] == '\n')
            {
               if (printLevel > 0)
                  printf("ERROR: in reading sample data.\n");
               fclose(fp);
               return -4;
            }
            sscanf(&lineIn[index],"%lg",&value);
            while (lineIn[index] != ' ') index++;
            if (lineIn[index] == '\0' || lineIn[index] == '\n')
            {
               if (printLevel > 0)
                  printf("ERROR: in reading sample data.\n");
               fclose(fp);
               return -4;
            }
            kk++;
         }
         while (lineIn[index] == ' ') index++;
         if (lineIn[index] == '\0' || lineIn[index] == '\n')
         {
            if (printLevel > 0)
               printf("ERROR: in reading sample probabilities.\n");
            fclose(fp);
            return -4;
         }
         else sscanf(&lineIn[index],"%lg",&value);
      }
      fclose(fp);
   }
   else if (which == 2)
   {
      for (ii = 0; ii < intparam1; ii++)
      {
         fgets(lineIn, 20000, fp);
         index = kk = 0;
         while (kk < nInps)
         {
            while (lineIn[index] == ' ') index++;
            if (lineIn[index] == '\0' || lineIn[index] == '\n')
            {
               if (printLevel > 0)
                  printf("ERROR: in reading sample data.\n");
               fclose(fp);
               return -4;
            }
            sscanf(&lineIn[index],"%lg",&value);
            while (lineIn[index] != ' ') index++;
            if (lineIn[index] == '\0' || lineIn[index] == '\n')
            {
               if (printLevel > 0)
                  printf("ERROR: in reading sample data.\n");
               fclose(fp);
               return -4;
            }
            kk++;
         }
      }
      fclose(fp);
   }
   else return -1;
   return 0;
}

// ************************************************************************
// compute MCMC specification file format
// ------------------------------------------------------------------------
int checkMCMCFileFormat(char *fname, int option, int printLevel)
{
   int    cnt, nExps, nOuts, nDesigns, ii, kk, jj;
   double ddata, ddata2;
   char   lineIn[10001], cword[1001];
   FILE   *fp=NULL;

   if (option == 0)
   {
      cnt = strlen(fname);
      if (cnt <= 5000)
      {
         fname[cnt] = '\0';
         fp = fopen(fname, "r");
         if (fp == NULL) 
         {
            printf("ERROR: file %s not found.\n", fname);
            return -1;
         }
      }
      else return -1;
      lineIn[0] = '#';
      while (lineIn[0] == '#') fgets(lineIn, 2000, fp);
      sscanf(lineIn, "%s", cword); 
      if (!strcmp(cword, "PSUADE_BEGIN")) fgets(lineIn, 2000, fp);
      sscanf(lineIn, "%d %d %d", &nExps, &nOuts, &nDesigns);
      if (nExps <= 0)
      {
         if (printLevel > 0)
            printf("ERROR: No experimental data available.\n");
         fclose(fp);
         return -1;
      }
      if (nOuts <= 0)
      {
         if (printLevel > 0)
            printf("ERROR: number of output parameters <= 0.\n");
         fclose(fp);
         return -1;
      }
      if (nDesigns < 0)
      {
         if (printLevel > 0)
            printf("ERROR: number of design parameters < 0.\n");
         fclose(fp);
         return -1;
      }
      if (nDesigns > 0)
      {
         cnt = 0;
         for (ii = 0; ii < nDesigns; ii++)
         {
            fscanf(fp, "%d", &kk);
            if (kk <= 0)
            {
               if (printLevel > 0)
                  printf("ERROR: wrong design parameter index.\n");
               fclose(fp);
               return -1;
            }
            if (kk <= cnt)
            {
               if (printLevel > 0)
                 printf("ERROR: wrong design parameter order (ascending?)\n");
               fclose(fp);
               return -1;
            }
            cnt = kk;
         }
      }
      for (ii = 0; ii < nExps; ii++)
      {
         fscanf(fp, "%d", &kk);
         if (kk != ii+1)
         {
            if (printLevel > 0) printf("ERROR: wrong experiment index.\n");
            fclose(fp);
            return -1;
         }
         for (jj = 0; jj < nDesigns; jj++) fscanf(fp, "%lg", &ddata);
         for (jj = 0; jj < nOuts; jj++)
         {
            fscanf(fp, "%lg %lg", &ddata, &ddata2);
            if (ddata2 < 0.0)
            {
               fclose(fp);
               if (printLevel > 0)
                  printf("ERROR: std dev zero or negative.\n");
               return -1;
            }
         }
      }
      fclose(fp);
   }
   else if (option == 1)
   {
      cnt = strlen(fname);
      if (cnt <= 5000)
      {
         fname[cnt] = '\0';
         fp = fopen(fname, "r");
         if (fp == NULL) 
         {
            printf("ERROR: file %s not found.\n", fname);
            return -1;
         }
      }
      else return -1;
      fscanf(fp, "%s", lineIn);
      if (strcmp(lineIn, "PSUADE_BEGIN"))
      {
         fclose(fp);
         return -1;
      }
      fscanf(fp, "%d %d", &cnt, &kk);
      if (cnt <= 0 || kk <= 0) return -1;
      fgets(lineIn, 5000, fp);
      while (1)
      {
         ii = getc(fp);
         if (ii != '#')
         {
            ungetc(ii, fp);
            break;
         }
         else fgets(lineIn, 5000, fp);
      }
      for (ii = 0; ii < cnt; ii++)
      {
         fscanf(fp, "%d", &jj);
         if ((ii+1) != jj) 
         {
            fclose(fp);
            return -1;
         }
         for (jj = 0; jj < kk; jj++) fscanf(fp,"%lg", &ddata);
      }
      fclose(fp);
   }
   return 0;
}

// ************************************************************************
// check S PDF file format
// ------------------------------------------------------------------------
int checkSPDFFileFormat(char *fname, int printLevel)
{
   int    nSamp, nInps, ii, jj, nn;
   double ddata;
   char   pString[10001];
   FILE   *fp;

   fp = fopen(fname, "r");
   if (fp == NULL)
   {
      if (printLevel > 0)
         printf("ERROR: S PDF file %s does not exist.\n",fname);
      return -1;
   }
   fscanf(fp, "%s", pString);
   if (strcmp(pString, "PSUADE_BEGIN"))
   {
      fclose(fp);
      fp = fopen(fname, "r");
   }
   fscanf(fp, "%d %d", &nSamp, &nInps);
   if (nSamp < 1)
   {
      if (printLevel > 0)
         printf("ERROR: S PDF sample file has nSamples <= 0.\n");
      return -1;
   }
   if (nInps < 1)
   {
      if (printLevel > 0)
         printf("ERROR: S PDF sample file has nInputs <= 0.\n");
      return -1;
   }
   fgets(pString, 1000, fp);
   while (1)
   {
      nn = getc(fp);
      if (nn == '#') fgets(pString, 1000, fp);
      else
      {
         ungetc(nn, fp);
         break;
      }
   }
   for (ii = 0; ii < nSamp; ii++)
   {
      fscanf(fp, "%d", &nn);
      if (nn != (ii+1))
      {
         if (printLevel > 0)
            printf("ERROR: S PDF sample file has wrong sample index.\n");
         return -1;
      }
      for (jj = 0; jj < nInps; jj++) fscanf(fp, "%lg", &ddata);
      fgets(pString, 1000, fp);
   }
   fclose(fp);
   return 0;
}

// ************************************************************************
// read standard file data 
// ------------------------------------------------------------------------
int readStdDataFile(char *fname, psMatrix &Amat)
{
  int      ii, idata, count, nSamp, nInps;
  double   ddata;
  ifstream infile;

  infile.open(fname, ios::in);
  if (!infile.is_open())
  {
    cout << "readStdDataFile ERROR: cannot read file " << fname << endl;
    return 1;
  }
  infile >> nSamp;
  if (nSamp <= 0)
  {
    cout << "readStdDataFile ERROR: wrong nSamples <= 0" << endl;
    infile.close();
    return 1;
  }
  infile >> nInps;
  if (nInps <= 0)
  {
    cout << "readStdDataFile ERROR: nInps <= 0" << endl;
    infile.close();
    return 1;
  }
  
  Amat.setDim(nSamp, nInps);
  count = 1;
  while (!infile.eof())
  {
    infile >> idata;
    if (idata != count)
    {
      cout << "readStdDataFile ERROR: wrong format " << endl;
      cout << "    expected sample index = " << count << endl;
      cout << "    actual   sample index = " << idata << endl;
      infile.close();
      return 1;
    }
    for (ii = 0; ii < nInps; ii++)
    {
      infile >> ddata;
      Amat.setEntry(count-1, ii, ddata);
    }
    count++;
    if (count == nSamp) break;
  }
  infile.close();
  if (count != nSamp)
  {
    cout << "readStdDataFile ERROR: not enough data" << endl;
    return 1;
  }
  return 0;
}

// ************************************************************************
// read standard iread file data 
// ------------------------------------------------------------------------
int readIReadDataFile(char *fname, psMatrix &Amat)
{
  int      ii, idata, count, nSamp, nInps;
  double   ddata;
  ifstream infile;

  infile.open(fname, ios::in);
  if (!infile.is_open())
  {
    printf("readIReadDataFile ERROR: cannot read file %s.\n", fname);
    return 1;
  }
  infile >> nSamp;
  if (nSamp <= 0)
  {
    printf("readIReadDataFile ERROR: wrong nSamples <= 0\n");
    infile.close();
    return 1;
  }
  infile >> nInps;
  if (nInps <= 0)
  {
    printf("readIReadDataFile ERROR: nInps <= 0");
    infile.close();
    return 1;
  }
  
  Amat.setDim(nSamp, nInps);
  count = 1;
  while (!infile.eof())
  {
    infile >> idata;
    if (idata != count)
    {
      printf("readIReadDataFile ERROR: wrong format\n");
      printf("    expected sample index = %d\n",count);
      printf("    actual   sample index = %d\n",idata);
      infile.close();
      return 1;
    }
    for (ii = 0; ii < nInps; ii++)
    {
      infile >> ddata;
      Amat.setEntry(count-1, ii, ddata);
    }
    count++;
    if (count > nSamp) break;
  }
  infile.close();
  if (count <= nSamp)
  {
    printf("readIReadDataFile ERROR: not enough data\n");
    return 1;
  }
  return 0;
}

// ************************************************************************
// read sample input data 
// ------------------------------------------------------------------------
int readSampleInputFile(const char *fname, psIVector &rsIndices,
                        psMatrix &Amat)
{
  int      ii, jj, nn, nSamp, nInps;
  double   ddata;
  char     cString[10001], inChar;
  string   line;
  ifstream infile;
  psIVector indices;

  infile.open(fname, ios::in);
  if (!infile.is_open())
  {
    cout << "readSampleFile ERROR: cannot read file " << fname << endl;
    return -1;
  }
  infile >> cString;
  if (strcmp(cString, "PSUADE_BEGIN"))
  {
    infile.close();
    infile.open(fname, ios::in);
  }
  infile >> nSamp;
  if (nSamp <= 0)
  {
    cout << "readSampleFile ERROR: wrong nSamples <= 0" << endl;
    infile.close();
    return -1;
  }
  infile >> nInps;
  if (nInps <= 0)
  {
    cout << "readSampleFile ERROR: nInps <= 0" << endl;
    infile.close();
    return -1;
  }
  std::getline(infile, line);
  while (!infile.eof())
  {
    inChar = infile.peek();
    if (inChar != '#') break;
    std::getline(infile, line);
  }
  if (infile.eof())
  {
    cout << "readSampleFile ERROR: EOF reached.\n";
    infile.close();
    return -1;
  }
  Amat.setDim(nSamp, nInps);
  indices.setLength(nInps);
  nn = 0;
  for (ii = 0; ii < rsIndices.length(); ii++)
  {
    if (rsIndices[ii] >= 1000) indices[nn++] = rsIndices[ii] - 1000; 
  }
  for (ii = 0; ii < nSamp; ii++)
  {
    infile >> nn;
    if (nn != (ii+1))
    {
      cout << "readSampleFile ERROR: invalid sample number.\n";
      cout << "           Expected: " << ii+1 << endl;
      cout << "           Read:     " << nn << endl;
      cout << "Advice: check your data format at line " << ii+2 << endl;
      cout << "Correct Format: \n";
      cout << "line 1: (optional) PSUADE_BEGIN\n";
      cout << "line 2: <number of sample points> <number of inputs>\n";
      cout << "line 3: (optional) : '#' followed by input names\n";
      cout << "line 4: 1 sample point 1 inputs \n";
      cout << "line 5: 2 sample point 2 inputs \n";
      cout << "line 6: 3 sample point 3 inputs \n";
      cout << "...\n";
      cout << "line n: (optional) PSUADE_END\n";
      infile.close();
      return -1;
    }
    for (jj = 0; jj < nInps; jj++)
    {
      infile >> ddata;
      nn = indices[jj];
      Amat.setEntry(ii,nn,ddata);
    }
  }
  infile.close();
  return 0;
}

// ************************************************************************
// findMaxMinDistance (N (dimension of candMatrix) choose M)
// candMatrix: distance matrix for the candidate sample points
// distMatrix: distance matrix for candidate + pre-selected sample points
// vecIndMap:  mapping between candMatrix and distMatrix
//             (if kk == vecIndMap[ii], rank kk in candMatrix is
//              rank ii in distMatrix) 
// vecIndTest: passing runtime selection for recursion
// currPos:    sample point for consideration in the present stage
// Upon return, vecBestInds[ii] = 1 means candidate ii (based on rank in 
//                             the candMatrix) has been selected
// ------------------------------------------------------------------------
int findMaxMinDistance(psIVector vecIndTest, int currPos, int M,
                    int MLeft,psMatrix &candMatrix, psMatrix &distMatrix,
                    psIVector &vecIndMap, psIVector &vecBestInds,
                    double &maxMinDist, int &iter, int printLevel)
{
  int    N, doEval, ii, kk, jj, ind, ind2, debugLevel;
  double ddata, sumDist, minDist;
  FILE   *fp;
  psIVector vecLocalInds;

  debugLevel = printLevel;
  fp = fopen("ps_print", "r");
  if (fp != NULL)
  {
    debugLevel = 5;
    fclose(fp);
  }

  N = candMatrix.nrows();
  vecLocalInds = vecIndTest;
  doEval = 0;
  if (MLeft == 0)
  {
    doEval = 1;
  }
  if (MLeft != 0 && (N-currPos) <= MLeft)
  {
    for (ii = 0; ii < MLeft; ii++) vecLocalInds[currPos+ii] = 1;
    doEval = 1;
  }
  ii = vecLocalInds.sum();

  if (doEval == 1 && ii > 0)
  {
    minDist = 1e35;
    for (kk = 0; kk < N; kk++)
    {
      if (vecLocalInds[kk] != 0)
      { 
        for (ii = kk+1; ii < N; ii++)
        {
          if (vecLocalInds[ii] != 0)
          {
            ddata = candMatrix.getEntry(kk,ii);
            sumDist += ddata;
            if (ddata < minDist) minDist = ddata;
          }
        }
        for (ii = 0; ii < vecIndMap.length(); ii++)
          if (vecIndMap[ii] == kk) break;
        ind2 = ii;
        for (ii = 0; ii < vecIndMap.length(); ii++)
        {
          if (vecIndMap[ii] < 0)
          {
            ind = - vecIndMap[ii] - 1;
            ddata = distMatrix.getEntry(ind2,ind);
            if (ddata < minDist) minDist = ddata;
          }
        }
      }
    }
    if (debugLevel > 4)
    {
      kk = 0;
      printf("Enumeration: ");
      for (ii = 0; ii < N; ii++)
      {
        if (vecLocalInds[ii] == 1)
        {
          printf("%5d ", ii+1); 
          kk++;
          if (kk > 10) break;
        }
      }
      printf("- min distance = %e ", minDist);
      if (minDist > maxMinDist) printf("*\n");
      else                      printf("\n");
    }
    if (minDist > maxMinDist)
    {
      maxMinDist = minDist;
      vecBestInds = vecLocalInds;
      if (debugLevel > 0)
        printf("Iteration %d: Selected sample points (maxMinDist = %e)\n",
               iter, maxMinDist);
      kk = 0;
      for (ii = 0; ii < N; ii++)
      {
        if (vecBestInds[ii] == 1)
        {
          printf("%4d ", ii+1); 
          kk++;
          if (kk >= 10)
          {
            kk = 0;
            printf("\n");
          }
        }
      }
      printf("\n");
    }  
    iter++;
    return 0;
  }

  vecLocalInds[currPos] = 0;
  findMaxMinDistance(vecLocalInds, currPos+1, M, MLeft,
                     candMatrix,distMatrix,vecIndMap,vecBestInds,
                     maxMinDist,iter,printLevel);
  vecLocalInds[currPos] = 1;
  findMaxMinDistance(vecLocalInds, currPos+1, M, MLeft-1,
                     candMatrix,distMatrix,vecIndMap,vecBestInds,
                     maxMinDist,iter,printLevel);
  return 0;
}

// ************************************************************************
// Given a sample, compose its covariance matrix and compute eigenvalues
// ------------------------------------------------------------------------
int computeFromSampleCovMatEigen(psMatrix matSample,psVector &vecEigs,
                                 psMatrix &matCov)
{
  int      ii, jj, kk, nSamp, nInps;
  double   dmean, ddata;
  psMatrix matLocalSample, matEig;

  //**/ centering
  nInps = matSample.ncols();
  nSamp = matSample.nrows();
  matLocalSample = matSample;
  for (ii = 0; ii < nInps; ii++)
  {
    dmean = 0.0;
    for (kk = 0; kk < nSamp; kk++)
      dmean += matSample.getEntry(kk,ii);
    dmean /= (double) nSamp;
    for (kk = 0; kk < nSamp; kk++)
    {
      ddata = matSample.getEntry(kk,ii) - dmean;
      matLocalSample.setEntry(kk,ii,ddata);
    }
  }
  //**/ construct covariance matrix
  matCov.setDim(nInps, nInps);
  for (ii = 0; ii < nInps; ii++)
  {
    for (jj = ii; jj < nInps; jj++)
    {
      ddata = 0.0;
      for (kk = 0; kk < nSamp; kk++)
         ddata += matLocalSample.getEntry(kk,ii) *
                  matLocalSample.getEntry(kk,jj);
      ddata /= (double) nSamp;
      matCov.setEntry(ii,jj,ddata);;
      matCov.setEntry(jj,ii,ddata);
    }
  }
  //**/ compute eigenvalues
  matCov.eigenSolve(matEig, vecEigs, 1);
  return 0;
}

// ************************************************************************
// Python header for response surface interpolators
// ------------------------------------------------------------------------
void fwriteRSPythonHeader(FILE *fp)
{
   fprintf(fp,"#!/usr/bin/python\n");
   fprintf(fp,"###################################################\n");
   fprintf(fp,"# Response surface interpolator from PSUADE\n");
   fprintf(fp,"#==================================================\n");
   fprintf(fp,"# This file contains information for interpolation\n");
   fprintf(fp,"# using response surface. Follow the steps below:\n");
   fprintf(fp,"#  1. move this file to *.py file (e.g. interpolate.py)\n");
   fprintf(fp,"#  2. make sure the first line points to your Python\n");
   fprintf(fp,"#  3. prepare your new sample points to be interpolated\n");
   fprintf(fp,"#     in a text file (e.g. infile) having the format below:\n");
   fprintf(fp,"#    <number of sample points M> <number of inputs n>\n");
   fprintf(fp,"#    1 input1 input2 ...inputn\n");
   fprintf(fp,"#    2 input1 input2 ...inputn\n");
   fprintf(fp,"#    ....\n");
   fprintf(fp,"#    M input1 input2 ...inputn\n");
   fprintf(fp,"#  4. run: interpolate.py infile outfile\n");
   fprintf(fp,"#     where <outfile> will have the interpolated values.\n");
   fprintf(fp,"#==================================================\n");
   fprintf(fp,"import sys\n");
   fprintf(fp,"import string\n");
   fprintf(fp,"import math\n\n");
}

// ************************************************************************
// Python header for response surface interpolators
// ------------------------------------------------------------------------
void fwriteRSPythonCommon(FILE *fp)
{
   fprintf(fp,"###################################################\n");
   fprintf(fp,"# Function to get input data from PSUADE-generated\n");
   fprintf(fp,"# parameter files (standard format, do not change).\n");
   fprintf(fp,"# The return data will contain the inputs values.\n");
   fprintf(fp,"#==================================================\n");
   fprintf(fp,"format = 0\n");
   fprintf(fp,"def getInputData(inFileName):\n");
   fprintf(fp,"   inFile  = open(inFileName, 'r')\n");
   fprintf(fp,"   lineIn  = inFile.readline()\n");
   fprintf(fp,"   nCols   = lineIn.split()\n");
   fprintf(fp,"   if len(nCols) == 1:\n");
   fprintf(fp,"     format  = 1\n");
   fprintf(fp,"     nInputs = eval(nCols[0])\n");
   fprintf(fp,"     inData  = (nInputs) * [0]\n");
   fprintf(fp,"     for ind in range(nInputs):\n");
   fprintf(fp,"       lineIn = inFile.readline()\n");
   fprintf(fp,"       nCols  = lineIn.split()\n");
   fprintf(fp,"       inData[ind] = eval(nCols[0])\n");
   fprintf(fp,"   else:\n");
   fprintf(fp,"     format  = 2\n");
   fprintf(fp,"     nSamp   = eval(nCols[0])\n");
   fprintf(fp,"     nInputs = eval(nCols[1])\n");
   fprintf(fp,"     inData  = (nSamp * nInputs) * [0]\n");
   fprintf(fp,"     for cnt in range(nSamp):\n");
   fprintf(fp,"       lineIn = inFile.readline()\n");
   fprintf(fp,"       nCols  = lineIn.split()\n");
   fprintf(fp,"       for ind in range(nInputs):\n");
   fprintf(fp,"         inData[cnt*nInputs+ind] = eval(nCols[ind+1])\n");
   fprintf(fp,"   inFile.close()\n");
   fprintf(fp,"   return inData\n");
   fprintf(fp,"###################################################\n");
   fprintf(fp,"# Function to generate output file\n");
   fprintf(fp,"# This function writes the output data (which should\n");
   fprintf(fp,"# have been generated in outData) to the PSUADE-based\n");
   fprintf(fp,"# output file.\n");
   fprintf(fp,"#==================================================\n");
   fprintf(fp,"def genOutputFile(outFileName, outData):\n");
   fprintf(fp,"   if format == 1:\n");
   fprintf(fp,"      nLeng = len(outData)\n");
   fprintf(fp,"      for ind in range(nLeng):\n");
   fprintf(fp,"         outfile.write(\"%%e \" %% outData[ind])\n");
   fprintf(fp,"   else:\n");
   fprintf(fp,"      nLeng = len(outData) / 2\n");
   fprintf(fp,"      nLeng = math.ceil(nLeng)\n");
   fprintf(fp,"      outfile = open(outFileName, 'w')\n");
   fprintf(fp,"      for ind in range(nLeng):\n");
   fprintf(fp,"         outfile.write(\"%%d \" %% (ind+1))\n");
   fprintf(fp,"         outfile.write(\"%%e \" %% outData[2*ind])\n");
   fprintf(fp,"         outfile.write(\"%%e \\n\" %% outData[2*ind+1])\n");
   fprintf(fp,"   outfile.close()\n");
   fprintf(fp,"   return\n");
   fprintf(fp,"###################################################\n");
}

// ************************************************************************
// generate matrix of histogram plots
// ------------------------------------------------------------------------
int genMatlabPlotFile(int nInps,double *LB,double *UB,int nSamps, 
                      double *sample, char **inpNames, char *fname, 
                      int nbins)
{
   int    kk, kk2, ii, ii2, jj, jj2, index, index2, sumBins;
   int    **bins, ****bins2;
   double ddata, ddata2;
   char   cfname[1001], charString[1001];;
   FILE   *fp;

   bins = new int*[nbins];
   for (ii = 0; ii < nbins; ii++)
   {
      bins[ii] = new int[nInps];
      for (jj = 0; jj < nInps; jj++) bins[ii][jj] = 0;
   }
   bins2 = new int***[nbins];
   for (jj = 0; jj < nbins; jj++)
   {
      bins2[jj] = new int**[nbins];
      for (jj2 = 0; jj2 < nbins; jj2++)
      {
         bins2[jj][jj2] = new int*[nInps];
         for (ii = 0; ii < nInps; ii++)
         {
            bins2[jj][jj2][ii] = new int[nInps];
            for (ii2 = 0; ii2 < nInps; ii2++)
               bins2[jj][jj2][ii][ii2] = 0;
         }
      }
   }
   for (ii = 0; ii < nSamps; ii++)
   {
      for (ii2 = 0; ii2 < nInps; ii2++)
      {
         ddata = sample[ii*nInps+ii2];
         ddata = (ddata - LB[ii2]) / (UB[ii2] - LB[ii2]);
         index = (int) (ddata * nbins);
         if (index >= nbins) index = nbins - 1;
         if (index < 0) index = 0;
         bins[index][ii2]++;
      }
      for (ii2 = 0; ii2 < nInps; ii2++)
      {
         ddata = sample[ii*nInps+ii2];
         ddata = (ddata - LB[ii2]) / (UB[ii2] - LB[ii2]);
         index = (int) (ddata * nbins);
         if (index >= nbins) index = nbins - 1;
         if (index < 0) index = 0;
         for (jj = 0; jj < nInps; jj++)
         {
            ddata2 = sample[ii*nInps+jj];
            ddata2 = (ddata2 - LB[jj]) / (UB[jj] - LB[jj]);
            index2 = (int) (ddata2 * nbins);
            if (index2 >= nbins) index2 = nbins - 1;
            if (index2 < 0) index2 = 0;
            bins2[index][index2][ii2][jj]++;
         }
      }
   }

   strcpy(cfname, fname);
   fp = fopen(cfname, "w");
   if (fp == NULL)
   {
      printOutTS(PL_ERROR, "ERROR: cannot open %s file.\n", cfname);
      return 0;
   }
   sprintf(charString,"This file shows posteriors plots");
   fwriteComment(fp, charString);
   sprintf(charString,"ns  - set to 1 for 1-step smoothing of 2D contours");
   fwriteComment(fp, charString);
   sprintf(charString,"ns1 - set to 1 for 1-step smoothing of 1D histgrams");
   fwriteComment(fp, charString);
   fprintf(fp, "ns  = 0;\n");
   fprintf(fp, "ns1 = 0;\n");
   fwritePlotCLF(fp);
   fprintf(fp, "active = [\n");
   for (kk = 0; kk < nInps; kk++) fprintf(fp, "1\n");
   fprintf(fp, "];\n");
   fprintf(fp, "L = [\n");
   for (kk = 0; kk < nInps; kk++) fprintf(fp, "%24.16e ",LB[kk]);
   fprintf(fp, "];\n");
   fprintf(fp, "U = [\n");
   for (kk = 0; kk < nInps; kk++) fprintf(fp, "%24.16e ",UB[kk]);
   fprintf(fp, "];\n");
   fprintf(fp, "iStr = {\n");
   for (kk = 0; kk < nInps-1; kk++) fprintf(fp, "'%s',", inpNames[kk]);
   fprintf(fp, "'%s'};\n", inpNames[nInps-1]);
   fprintf(fp, "X = zeros(%d,%d);\n", nInps, nbins);
   fprintf(fp, "D = zeros(%d,%d);\n", nInps, nbins);
   fprintf(fp, "NC = zeros(%d,%d,%d,%d);\n",nInps,nInps,nbins,nbins);
   for (kk = 0; kk < nInps; kk++)
   {
      for (kk2 = 0; kk2 < nInps; kk2++)
      {
         if (kk == kk2)
         {
            fprintf(fp, "X(%d,:) = [\n", kk+1);
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%24.16e ",(UB[kk]-LB[kk])/nbins*(jj+0.5)+LB[kk]);
            fprintf(fp, "];\n");
            fprintf(fp, "D(%d,:) = [\n", kk+1);
            sumBins = 0;
            for (jj = 0; jj < nbins; jj++) sumBins += bins[jj][kk];
            if (sumBins == 0) sumBins = 1;
            for (jj = 0; jj < nbins; jj++)
               fprintf(fp, "%24.16e ",(double) bins[jj][kk]/(double) sumBins);
            fprintf(fp, "];\n");
         }
         else
         {
            fprintf(fp, "NC(%d,%d,:,:) = [\n", kk+1, kk2+1);
            for (jj = 0; jj < nbins; jj++)
            {
               for (jj2 = 0; jj2 < nbins; jj2++)
                  fprintf(fp, "%d ", bins2[jj][jj2][kk][kk2]);
               fprintf(fp, "\n");
            }
            fprintf(fp, "]';\n");
         }
      }
   }
   fprintf(fp,"fs = 12;\n");
   fprintf(fp,"nInps  = length(active);\n");
   fprintf(fp,"nPlots = 0;\n");
   fprintf(fp,"for ii = 1 : nInps\n");
   fprintf(fp,"   if (active(ii) == 1)\n");
   fprintf(fp,"      nPlots = nPlots + 1;\n");
   fprintf(fp,"      active(ii) = nPlots;\n");
   fprintf(fp,"   end;\n");
   fprintf(fp,"end;\n");
   fprintf(fp,"dzero = 0;\n");
   fprintf(fp,"for ii = 1 : nInps\n");
   fprintf(fp,"  for jj = ii : nInps\n");
   fprintf(fp,"    if (active(ii) ~= 0 & active(jj) ~= 0)\n");
   fprintf(fp,"      index = (active(ii)-1) * nPlots + active(jj);\n");
   fprintf(fp,"      subplot(nPlots,nPlots,index)\n");
   fprintf(fp,"      if (ii == jj)\n");
   fprintf(fp,"        n = length(D(ii,:));\n");
   fprintf(fp,"        DN = D(ii,:);\n");
   fprintf(fp,"        for kk = 1 : ns1\n");
   fprintf(fp,"          DN1 = DN;\n");
   fprintf(fp,"          for ll = 2 : n-1\n");
   fprintf(fp,"            DN(ll) = DN(ll) + DN1(ll+1);\n");
   fprintf(fp,"            DN(ll) = DN(ll) + DN1(ll-1);\n");
   fprintf(fp,"            DN(ll) = DN(ll) / 3;\n");
   fprintf(fp,"          end;\n");
   fprintf(fp,"        end;\n");
   fprintf(fp,"        bar(X(ii,:), DN, 1.0);\n");
   fprintf(fp,"        xmin = min(X(ii,:));\n");
   fprintf(fp,"        xmax = max(X(ii,:));\n");
   fprintf(fp,"        xwid = xmax - xmin;\n");
   fprintf(fp,"        xmin = xmin - 0.5 * xwid / %d;\n", nbins);
   fprintf(fp,"        xmax = xmax + 0.5 * xwid / %d;\n", nbins);
   fprintf(fp,"        ymax = max(DN);\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp,"        e = gce();\n");
      fprintf(fp,"        e.children.thickness = 2;\n");
      fprintf(fp,"        e.children.foreground = 0;\n");
      fprintf(fp,"        e.children.background = 2;\n");
      fprintf(fp,"        a = gca();\n");
      fprintf(fp,"        a.data_bounds=[xmin,0;xmax,ymax];\n");
      fprintf(fp,"        a.x_label.text = iStr(ii);\n");
      fprintf(fp,"        a.x_label.font_size = 3;\n");
      fprintf(fp,"        a.x_label.font_style = 4;\n");
      fprintf(fp,"//      a.grid = [1 1];\n");
      fprintf(fp,"        a.y_label.text = iStr(jj);\n");
      fprintf(fp,"        a.y_label.font_size = 3;\n");
      fprintf(fp,"        a.y_label.font_style = 4;\n");
      fprintf(fp,"        a.thickness = 2;\n");
      fprintf(fp,"        a.font_size = 3;\n");
      fprintf(fp,"        a.font_style = 4;\n");
      fprintf(fp,"        a.box = \"on\";\n");
   }
   else
   {
      fprintf(fp,"        axis([xmin xmax 0 ymax])\n");
      fprintf(fp,"        set(gca,'linewidth',2)\n");
      fprintf(fp,"        set(gca,'fontweight','bold')\n");
      fprintf(fp,"        set(gca,'fontsize',12)\n");
      fprintf(fp,
         "        xlabel(iStr(ii),'FontWeight','bold','FontSize',fs)\n");
      fprintf(fp,"        if (ii == 1)\n");
      fprintf(fp,
         "        ylabel('Probabilities','FontWeight','bold','FontSize',fs)\n");
      fprintf(fp,"        end;\n");
      fprintf(fp,"        grid on\n");
      fprintf(fp,"        box on\n");
   }
   fprintf(fp,"      else\n");
   fprintf(fp,"        n = length(X(jj,:));\n");
   fprintf(fp,"        XT = X(jj,:);\n");
   fprintf(fp,"        YT = X(ii,:);\n");
   fprintf(fp,"        HX = (XT(n) - XT(1)) / (n-1);\n");
   fprintf(fp,"        HY = (YT(n) - YT(1)) / (n-1);\n");
   fprintf(fp,"        ZZ = squeeze(NC(ii,jj,:,:));\n");
   fprintf(fp,"        for kk = 1 : ns\n");
   fprintf(fp,"          ZZ1 = ZZ;\n");
   fprintf(fp,"          for ll = 2 : n-1\n");
   fprintf(fp,"            for mm = 2 : n-1\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm+1);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm-1);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm+1);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm-1);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm+1);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm-1);\n");
   fprintf(fp,"              ZZ(ll,mm) = ZZ(ll,mm) / 9;\n");
   fprintf(fp,"            end;\n");
   fprintf(fp,"          end;\n");
   fprintf(fp,"        end;\n");
   fprintf(fp,"        ZZ = ZZ / (sum(sum(ZZ)));\n");
   if (psPlotTool_ == 1)
   {
      fprintf(fp,"        XX = [XT(1):HX:XT(n)];\n");
      fprintf(fp,"        YY = [YT(1):HY:YT(n)];\n");
      fprintf(fp,"        DD = splin2d(XX,YY,ZZ);\n");
      fprintf(fp,"        HX = 0.01 * (XT(n) - XT(1));\n");
      fprintf(fp,"        HY = 0.01 * (YT(n) - YT(1));\n");
      fprintf(fp,"        X2 = [XT(1):HX:XT(n)];\n");
      fprintf(fp,"        Y2 = [YT(1):HY:YT(n)];\n");
      fprintf(fp,"        [XI, YI] = ndgrid(X2, Y2);\n");
      fprintf(fp,"        disp('interpolation')\n");
      fprintf(fp,"        ZI =interp2d(XI, YI, XX, YY, DD, \"natural\");\n");
      fprintf(fp,"        disp('interpolation done')\n");
      fprintf(fp,"        ZB = ZI;\n");
      fprintf(fp,"        nX = length(X2);\n");
      fprintf(fp,"        nY = length(Y2);\n");
      fprintf(fp,"        for kk = 1 : nX\n");
      fprintf(fp,"          for ll = 1 : nY\n");
      fprintf(fp,"            ZI(kk,ll) = ZB(kk,nY-ll+1);\n");
      fprintf(fp,"          end;\n");
      fprintf(fp,"        end;\n");
      fprintf(fp,"        zmax = max(max(ZI));\n");
      fprintf(fp,"        zmin = min(min(ZI)) / zmax;\n");
      fprintf(fp,"        ZI   = ZI / zmax;\n");
      fprintf(fp,"        zmax = 1;\n");
      fprintf(fp,"        xset(\"colormap\",jetcolormap(64));\n");
      fprintf(fp,"        colorbar(zmin,zmax);\n");
      fprintf(fp,
           "        contour2d(X2,Y2,ZB,5,rect=[L(jj),L(ii),U(jj),U(ii)]);\n");
      fprintf(fp,"        xset(\"fpf\",\" \");\n");
      fprintf(fp,"        Matplot1((ZI-zmin)/(zmax-zmin)*64,[L(jj),L(ii),");
      fprintf(fp,"U(jj),U(ii)]);\n");
      fprintf(fp,"        a = gca();\n");
      fprintf(fp,"        a.x_label.text = iStr(jj);\n");
      fprintf(fp,"        a.x_label.font_size = 3;\n");
      fprintf(fp,"        a.x_label.font_style = 4;\n");
      fprintf(fp,"        a.y_label.text = iStr(ii);\n");
      fprintf(fp,"        a.y_label.font_size = 3;\n");
      fprintf(fp,"        a.y_label.font_style = 4;\n");
      fwritePlotAxesNoGrid(fp);
   }
   else
   {
      fprintf(fp,"        imagesc(ZZ')\n");
      fprintf(fp,"        xtick = L(ii):(U(ii)-L(ii))/4:U(ii);\n");
      fprintf(fp,"        set(gca,'XTick',0:n/4:n);\n");
      fprintf(fp,"        set(gca,'XTickLabel', xtick);\n");
      fprintf(fp,"        ytick = L(jj):(U(jj)-L(jj))/4:U(jj);\n");
      fprintf(fp,"        set(gca,'YTick',0:n/4:n);\n");
      fprintf(fp,"        set(gca,'YTickLabel', ytick);\n");
      fprintf(fp,"        set(gca,'YDir', 'normal');\n");
      fprintf(fp,"%% axis off\n");
      fprintf(fp,
         "        xlabel(iStr(jj),'FontWeight','bold','FontSize',fs)\n");
      fprintf(fp,
         "        ylabel(iStr(ii),'FontWeight','bold','FontSize',fs)\n");
      fwritePlotAxesNoGrid(fp);
   }
   fprintf(fp,"      end;\n");
   fprintf(fp,"    end;\n");
   fprintf(fp,"  end;\n");
   fprintf(fp,"end;\n");
   if (psPlotTool_ == 0)
   {
      fprintf(fp,"set(gcf,'NextPlot','add');\n");
      fprintf(fp,"axes;\n");
      fprintf(fp,"h=title('One- and Two-Input Distributions',");
      fprintf(fp,"'fontSize',fs,'fontWeight','bold');\n");
      fprintf(fp,"set(gca,'Visible','off');\n");
      fprintf(fp,"set(h,'Visible','on');\n");
   }
   fclose(fp);
   printOutTS(PL_INFO,"The %s file has been created.\n",fname);

   for (ii = 0; ii < nbins; ii++) delete [] bins[ii];
   delete [] bins;
   for (jj = 0; jj < nbins; jj++)
   {
      for (jj2 = 0; jj2 < nbins; jj2++)
      {
         for (ii = 0; ii < nInps; ii++) delete [] bins2[jj][jj2][ii];
         delete [] bins2[jj][jj2];
      }
      delete [] bins2[jj];
   }
   delete [] bins2;
   return 0;
}

// ************************************************************************
// plot stuff for scilab or matlab
// ------------------------------------------------------------------------
int fwritePlotAxes(FILE *fp)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.thickness = 2;\n");
      fprintf(fp, "a.font_size = 3;\n");
      fprintf(fp, "a.font_style = 4;\n");
      fprintf(fp, "a.box = \"on\";\n");
      fprintf(fp, "a.grid = [1 1];\n");
   }
   else
   {
      fprintf(fp, "set(gca,'linewidth',2)\n");
      fprintf(fp, "set(gca,'fontweight','bold')\n");
      fprintf(fp, "set(gca,'fontsize',12)\n");
      fprintf(fp, "grid on\n");
      fprintf(fp, "box on\n");
   }
   return 0;
}
// ************************************************************************
int fwritePlotAxesNoGrid(FILE *fp)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.thickness = 2;\n");
      fprintf(fp, "a.font_size = 3;\n");
      fprintf(fp, "a.font_style = 4;\n");
      fprintf(fp, "a.box = \"on\";\n");
   }
   else
   {
      fprintf(fp, "set(gca,'linewidth',2)\n");
      fprintf(fp, "set(gca,'fontweight','bold')\n");
      fprintf(fp, "set(gca,'fontsize',10)\n");
      fprintf(fp, "box on\n");
   }
   return 0;
}
// ************************************************************************
int fwritePlotXLabel(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.x_label.text = \"%s\";\n", label);
      fprintf(fp, "a.x_label.font_size = 3;\n");
      fprintf(fp, "a.x_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp,"xlabel('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotYLabel(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.y_label.text = \"%s\";\n", label);
      fprintf(fp, "a.y_label.font_size = 3;\n");
      fprintf(fp, "a.y_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp,"ylabel('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotZLabel(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.z_label.text = \"%s\";\n", label);
      fprintf(fp, "a.z_label.font_size = 3;\n");
      fprintf(fp, "a.z_label.font_style = 4;\n");
   }
   else
   {
      fprintf(fp,"zlabel('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotTitle(FILE *fp, const char *label)
{
   if (psPlotTool_ == 1)
   {
      fprintf(fp, "a = gca();\n");
      fprintf(fp, "a.title.text = \"%s\";\n", label);
      fprintf(fp, "a.title.font_size = 4;\n");
      fprintf(fp, "a.title.font_style = 4;\n");
   }
   else
   {
      fprintf(fp, "title('%s','FontWeight','bold','FontSize',12)\n",label);
   }
   return 0;
}
// ************************************************************************
int fwritePlotCLF(FILE *fp)
{
   if (psPlotTool_ == 1) fprintf(fp, "clf();\n");
   else
   {
      fprintf(fp, "if exist('noCLF') \n");
      fprintf(fp, "   hold off\n");
      fprintf(fp, "else\n");
      fprintf(fp, "   clf\n");
      fprintf(fp, "end;\n");
   }
   return 0;
}
// ************************************************************************
int fwritePlotFigure(FILE *fp, int num)
{
   if (psPlotTool_ == 1) fprintf(fp, "scf(%d);\n", num);
   else                  fprintf(fp, "figure(%d)\n", num);
   return 0;
}
// ************************************************************************
int fwriteHold(FILE *fp, int hswitch)
{
   if (psPlotTool_ == 1)
   {
      if (hswitch == 1) fprintf(fp, "set(gca(),\"auto_clear\",\"off\")\n");
      else              fprintf(fp, "set(gca(),\"auto_clear\",\"on\")\n");
   }
   else
   {
      if (hswitch == 1) fprintf(fp, "hold on\n");
      else              fprintf(fp, "hold off\n");
   }
   return 0;
}
// ************************************************************************
int fwriteComment(FILE *fp, char *line)
{
   if (psPlotTool_ == 1) fprintf(fp, "// %s\n", line);
   else                  fprintf(fp, "%% %s\n", line);
   return 0;
}

