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
// psStrings functions
// AUTHOR : CHARLES TONG
// DATE   : 2019
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "sysdef.h"
#include "Vector.h"
#include <string.h>
#include "psStrings.h"

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
psStrings::psStrings()
{
  nStrings_ = 0;
  strData_  = NULL;
}

// ************************************************************************
// operator=
// ------------------------------------------------------------------------
psStrings & psStrings::operator=(const psStrings & strs)
{
  int ii, leng;

  if (this == &strs) return *this;
  clean();
  nStrings_ = strs.nStrings_;
  if (nStrings_ > 0)
  {
    strData_ = new char*[nStrings_];
    for (ii = 0; ii < nStrings_; ii++)
    {
      if (strs.strData_[ii] != NULL)
      {
        leng = strlen(strs.strData_[ii]);
        strData_[ii] = new char[leng+2];
        strcpy(strData_[ii], strs.strData_[ii]);
      } 
    }
  } 
  return *this;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
psStrings::~psStrings()
{
  clean();
}

// ************************************************************************
// get string
// ------------------------------------------------------------------------
char *psStrings::operator[](int ind)
{
  if (ind < 0 || ind >= nStrings_)
  {
    printf("psStrings operator[] ERROR: index = %d (%d)\n",ind+1,nStrings_);
    exit(1);
  }
  return strData_[ind];
}

// ************************************************************************
// get strings
// ------------------------------------------------------------------------
char **psStrings::getStrings()
{
  return strData_;
}

// ************************************************************************
// get number of rows 
// ------------------------------------------------------------------------
int psStrings::numStrings()
{
  return nStrings_;
}

// ************************************************************************
// load strings
// ------------------------------------------------------------------------
int psStrings::load(psStrings &inStrs)
{
  int ii, leng;

  clean();
  assert(this != &inStrs);

  nStrings_ = inStrs.nStrings_;
  strData_ = new char*[nStrings_];
  for (ii = 0; ii < nStrings_; ii++)
  {
    leng = strlen(inStrs.strData_[ii]);
    strData_[ii] = new char[leng+2];
    strcpy(strData_[ii], inStrs.strData_[ii]);
  } 
  return 0;
}

// ************************************************************************
// load strings
// ------------------------------------------------------------------------
int psStrings::load(int numStrs, const char **inStrs)
{
  int ii, leng;

  clean();
  if (numStrs <= 0)
  {
    printf("psStrings load ERROR: numStrs <= 0\n");
    exit(1);
  }
  nStrings_ = numStrs;
  strData_ = new char*[nStrings_];
  for (ii = 0; ii < nStrings_; ii++)
  {
    leng = strlen(inStrs[ii]);
    strData_[ii] = new char[leng+2];
    strcpy(strData_[ii], inStrs[ii]);
  } 
  return 0;
}

// ************************************************************************
// load one string
// ------------------------------------------------------------------------
int psStrings::loadOneString(int row, const char *inStr)
{
  if (row < 0 || row >= nStrings_)
  {
    printf("psStrings loadOneString ERROR: invalid row (%d)\n",row);
    exit(1);
  }
  if (strData_[row] != NULL) delete [] strData_[row];
  int leng = strlen(inStr);
  strData_[row] = new char[leng+2];
  strcpy(strData_[row], inStr);
  return 0;
}

// ************************************************************************
// get one string
// ------------------------------------------------------------------------
char *psStrings::getOneString(int row)
{
  if (row < 0 || row >= nStrings_)
  {
    printf("psStrings getOneString ERROR: invalid index (%d)\n",row);
    exit(1);
  }
  return strData_[row];
}

// ************************************************************************
// set number of strings
// ------------------------------------------------------------------------
int psStrings::setNumStrings(int numStrs)
{
  clean();
  if (numStrs <= 0)
  {
    printf("psStrings setNumStrings ERROR: numStrs <= 0\n");
    exit(1);
  }
  nStrings_ = numStrs;
  strData_ = new char*[nStrings_];
  for (int ii = 0; ii < nStrings_; ii++) strData_[ii] = NULL;
  return 0;
} 

// ************************************************************************
// find out length of string
// ------------------------------------------------------------------------
int psStrings::strLength(int row)
{
  if (row < 0 || row >= nStrings_)
  {
    printf("psStrings strLength ERROR: invalid row (%d)\n",row);
    exit(1);
  }
  return strlen(strData_[row]);
}

// ************************************************************************
// find out length of string
// ------------------------------------------------------------------------
int psStrings::findString(const char *inStr)
{
  for (int ii = 0; ii < nStrings_; ii++) 
    if (!strcmp(strData_[ii], inStr)) return ii;
  return -1;
}

// ************************************************************************
// remove one string
// ------------------------------------------------------------------------
int psStrings::removeOneString(int row)
{
  int ii;
  if (row < 0 || row >= nStrings_)
  {
    printf("psStrings removeOneString ERROR: invalid row (%d)\n",row);
    exit(1);
  }
  if (nStrings_ == 1)
  {
    if (strData_[0] != NULL) delete [] strData_[0];
    delete [] strData_;
    return 0;
  }
  char **tmpStrs = new char*[nStrings_-1];
  for (ii = 0; ii < row; ii++)
    tmpStrs[ii] = strData_[ii];
  if (strData_[row] != NULL) delete [] strData_[row];
  for (ii = row; ii < nStrings_-1; ii++)
    tmpStrs[ii] = strData_[ii+1];
  delete [] strData_;
  strData_ = tmpStrs;
  nStrings_--;
  return 0;
}

// ************************************************************************
// add more string 
// ------------------------------------------------------------------------
int psStrings::addMoreStrings(int num)
{
  int ii;
  if (num <= 0)
  {
    printf("psStrings addMoreStrings ERROR: invalid entry (%d)\n",num);
    exit(1);
  }
  char **tmpStrs = new char*[nStrings_+num];
  for (ii = 0; ii < nStrings_; ii++)
    tmpStrs[ii] = strData_[ii];
  for (ii = nStrings_; ii < nStrings_+num; ii++)
    tmpStrs[ii] = NULL;
  delete [] strData_;
  strData_ = tmpStrs;
  nStrings_ += num;
  return 0;
}

// ************************************************************************
// swap two strings 
// ------------------------------------------------------------------------
int psStrings::swapStrings(int str1, int str2)
{
  if (str1 < 0 || str2 < 0)
  {
    printf("psStrings swapString ERROR: invalid entry (%d,%d)\n",
           str1,str2);
    exit(1);
  }
  if (str1 >= nStrings_ || str2 >= nStrings_)
  {
    printf("psStrings swapString ERROR: invalid entry (%d,%d)\n",
           str1,str2);
    exit(1);
  }
  if (str1 == str2)
  {
    printf("psStrings swapString INFO: same index - no change\n");
    return 0;
  }
  if (strData_[str1] == NULL || strData_[str2] == NULL)
  {
    printf("psStrings swapString INFO: empty string - no change\n");
    return 0;
  }
  char *tmpStr = strData_[str1];
  strData_[str1] = strData_[str2];
  strData_[str2] = tmpStr;
  return 0;
}

// ************************************************************************
// shuffle strings 
// ------------------------------------------------------------------------
int psStrings::shuffleStrings(psIVector & vecIndices)
{
  if (vecIndices.length() != nStrings_)
  {
    printf("psStrings shuffle ERROR: invalid vector length %d != %d\n",
           vecIndices.length(), nStrings_);
    exit(1);
  }
  char **tmpStrs = new char*[nStrings_];
  for (int ii = 0; ii < nStrings_; ii++) 
  {
    if (vecIndices[ii] < 0 || vecIndices[ii] >= nStrings_)
    {
      printf("psStrings shuffleStrings ERROR: invalid index %d\n",
             vecIndices[ii]);
      exit(1);
    }
    tmpStrs[ii] = strData_[vecIndices[ii]];
  }
  if (strData_ != NULL) delete [] strData_;
  strData_ = tmpStrs;
  return 0;
}

// ************************************************************************
// print
// ------------------------------------------------------------------------
void psStrings::print()
{
  for (int ii = 0; ii < nStrings_; ii++)
    if (strData_[ii] != NULL) printf("String %d = %s\n",ii+1,strData_[ii]);
}

// ************************************************************************
// clean
// ------------------------------------------------------------------------
void psStrings::clean()
{
  for (int ii = 0; ii < nStrings_; ii++)
    if (strData_[ii] != NULL) delete [] strData_[ii];
  if (strData_ != NULL) delete [] strData_;
  strData_ = NULL;
  nStrings_ = 0;
}

