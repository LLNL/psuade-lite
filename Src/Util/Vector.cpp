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
// psVector functions
// AUTHOR : CHARLES TONG
// DATE   : 2008
// ************************************************************************
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "sysdef.h"
#include "Vector.h"
#include "PsuadeUtil.h"

//#define PS_DEBUG

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
psVector::psVector()
{
#ifdef PS_DEBUG
  printf("psVector constructor\n");
#endif
  length_ = 0;
  Vec_ = NULL;
#ifdef PS_DEBUG
  printf("psVector constructor ends\n");
#endif
}

// ************************************************************************
// Copy Constructor by Bill Oliver
// ------------------------------------------------------------------------
psVector::psVector(const psVector & v)
{
  length_ = v.length_;
  Vec_ = NULL;
  if (length_ > 0)
  {
    Vec_ = new double[length_];
    for (int ii = 0; ii < length_; ii++) Vec_[ii] = v.Vec_[ii];
  }
}

// ************************************************************************
// overload operator= by Bill Oliver
// ------------------------------------------------------------------------
psVector & psVector::operator=(const psVector & v)
{
  if (this == &v) return *this;
  delete [] Vec_;
  Vec_ = NULL;
  length_ = v.length_;
  if (length_ > 0)
  {
    Vec_ = new double[length_];
    for (int ii = 0; ii < length_; ii++) Vec_[ii] = v.Vec_[ii];
  }
  return *this;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
psVector::~psVector()
{
#ifdef PS_DEBUG
  printf("psVector destructor\n");
#endif
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = 0;
#ifdef PS_DEBUG
  printf("psVector destructor ends\n");
#endif
}

// ************************************************************************
// get length 
// ------------------------------------------------------------------------
int psVector::length() 
{
  return length_;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int psVector::load(psVector &inVec)
{
#ifdef PS_DEBUG
  printf("psVector load\n");
#endif
  assert(this != &inVec);
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = inVec.length();
  if (length_ <= 0) return -1;
  Vec_ = new double[length_];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < length_; ii++) Vec_[ii] = inVec[ii];
#ifdef PS_DEBUG
  printf("psVector load ends\n");
#endif
  return 0;
}

// ************************************************************************
// set dimension
// ------------------------------------------------------------------------
int psVector::setLength(int leng)
{
#ifdef PS_DEBUG
  printf("psVector setLength\n");
#endif
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  assert(leng >= 0);
  if (leng == 0) return -1;
  length_ = leng;
  Vec_ = new double[leng];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < leng; ii++) Vec_[ii] = 0.0;
#ifdef PS_DEBUG
  printf("psVector setLength ends\n");
#endif
  return 0;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int psVector::load(int leng, double *data)
{
#ifdef PS_DEBUG
  printf("psVector load, length = %d\n", leng);
#endif
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  assert(leng > 0);
  assert(data != NULL);
  length_ = leng;
  Vec_ = new double[leng];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < leng; ii++) Vec_[ii] = data[ii];
#ifdef PS_DEBUG
  printf("psVector load ends\n");
#endif
  return 0;
}

// ************************************************************************
// set entry
// ------------------------------------------------------------------------
int psVector::setEntry(int ind, double ddata)
{
  if (ind < 0 || ind >= length_)
  {
    printf("psVector setEntry ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  Vec_[ind] = ddata;
  return 0;
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
double psVector::getEntry(int ind)
{
  if (ind < 0 || ind >= length_)
  {
    printf("psVector getEntry ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  return Vec_[ind];
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
double& psVector::operator[](int ind) 
{
  if (ind < 0 || ind >= length_)
  {
    printf("psVector operator[] ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  return Vec_[ind];
}

// ************************************************************************
// get vector
// ------------------------------------------------------------------------
double *psVector::getDVector() 
{
  return Vec_;
}

// ************************************************************************
// find max 
// ------------------------------------------------------------------------
double psVector::max()
{
  int    ii;
  double dmax = -PSUADE_UNDEFINED;
  for (ii = 0; ii < length_; ii++) 
    if (Vec_[ii] > dmax) dmax = Vec_[ii];
  return dmax;
}
   
// ************************************************************************
// find min 
// ------------------------------------------------------------------------
double psVector::min()
{
  int    ii;
  double dmin = PSUADE_UNDEFINED;
  for (ii = 0; ii < length_; ii++) 
    if (Vec_[ii] < dmin) dmin = Vec_[ii];
  return dmin;
}
   
// ************************************************************************
// compute vector sum 
// ------------------------------------------------------------------------
double psVector::sum()
{
  int    ii;
  double dsum=0.0;
  for (ii = 0; ii < length_; ii++) dsum += Vec_[ii];
  return dsum;
}

// ************************************************************************
// scale vector 
// ------------------------------------------------------------------------
void psVector::scale(double alpha)
{
  for (int ii = 0; ii < length_; ii++) Vec_[ii] *= alpha;
}

// ************************************************************************
// sort 
// ------------------------------------------------------------------------
void psVector::sort()
{
  sortDbleList(length_, Vec_);
}
   
// ************************************************************************
// add to vector
// ------------------------------------------------------------------------
int psVector::addElements(int leng, double *data)
{
#ifdef PS_DEBUG
  printf("psVector add, length = %d\n", leng);
#endif
  int    ii;
  double *tmpVec = Vec_;
  Vec_ = new double[leng+length_];
  for (ii = 0; ii < length_; ii++) Vec_[ii] = tmpVec[ii];
  if (data == NULL)
       for (ii = 0; ii < leng; ii++) Vec_[length_+ii] = 0.0;
  else for (ii = 0; ii < leng; ii++) Vec_[length_+ii] = data[ii];
  delete [] tmpVec;
  length_ += leng;
#ifdef PS_DEBUG
   printf("psVector add ends\n");
#endif
   return 0;
}

// ************************************************************************
// clean
// ------------------------------------------------------------------------
void psVector::clean()
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = 0;
}

void psVector::print(char *name)
{
  for (int ii = 0; ii < length_; ii++)
    printf("%s %5d = %16.8e\n", name, ii+1, Vec_[ii]);
}

// ************************************************************************
// generate a subvector
// ------------------------------------------------------------------------
void psIVector::subvector(int ibeg, int iend)
{
  int  ii, leng;
  int  *tmpVec;

  leng = iend - ibeg + 1;
  if (leng < 0 || ibeg < 0 || iend >= length_)
  {
    printf("psIVector subvector range ERROR: beg/end = %d %d\n",ibeg,iend);
    exit(1);
  }
  tmpVec = Vec_;
  length_ = leng;
  Vec_ = new int[leng];
  for (ii = ibeg; ii < iend+1; ii++) Vec_[ii-ibeg] = tmpVec[ii];
  if (tmpVec != NULL) delete [] tmpVec;
  return;
}

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
psIVector::psIVector()
{
  length_ = 0;
  Vec_ = NULL;
}

// ************************************************************************
// Copy Constructor 
// ------------------------------------------------------------------------
psIVector::psIVector(const psIVector & v)
{
  length_ = v.length_;
  Vec_ = NULL;
  if (length_ > 0)
  {
    Vec_ = new int[length_];
    for (int ii = 0; ii < length_; ii++) Vec_[ii] = v.Vec_[ii];
  }
}

// ************************************************************************
// set entry
// ------------------------------------------------------------------------
int psIVector::setEntry(int ind, int idata)
{
  if (ind < 0 || ind >= length_)
  {
    printf("psIVector setEntry ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  Vec_[ind] = idata;
  return 0;
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
int psIVector::getEntry(int ind)
{
  if (ind < 0 || ind >= length_)
  {
    printf("psIVector getEntry ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  return Vec_[ind];
}

// ************************************************************************
// compute vector sum 
// ------------------------------------------------------------------------
int psIVector::sum()
{
  int    ii, isum=0;
  for (ii = 0; ii < length_; ii++) isum += Vec_[ii];
  return isum;
}

// ************************************************************************
// overload operator= 
// ------------------------------------------------------------------------
psIVector & psIVector::operator=(const psIVector & v)
{
  if (this == &v) return *this;
  delete [] Vec_;
  Vec_ = NULL;
  length_ = v.length_;
  if (length_ > 0)
  {
    Vec_ = new int[length_];
    for (int ii = 0; ii < length_; ii++) Vec_[ii] = v.Vec_[ii];
  }
  return *this;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
psIVector::~psIVector()
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = 0;
}

// ************************************************************************
// get length 
// ------------------------------------------------------------------------
int psIVector::length() 
{
  return length_;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int psIVector::load(psIVector &inVec)
{
  assert(this != &inVec);
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = inVec.length();
  if (length_ <= 0) return -1;
  Vec_ = new int[length_];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < length_; ii++) Vec_[ii] = inVec[ii];
  return 0;
}

// ************************************************************************
// set dimension
// ------------------------------------------------------------------------
int psIVector::setLength(int leng)
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  assert(leng >= 0);
  if (leng == 0) return -1;
  length_ = leng;
  Vec_ = new int[leng];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < leng; ii++) Vec_[ii] = 0;
  return 0;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int psIVector::load(int leng, int *data)
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  assert(leng > 0);
  assert(data != NULL);
  length_ = leng;
  Vec_ = new int[leng];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < leng; ii++) Vec_[ii] = data[ii];
  return 0;
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
int& psIVector::operator[](int ind) 
{
  if (ind < 0 || ind >= length_)
  {
    printf("psIVector operator[] ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  return Vec_[ind];
}

// ************************************************************************
// get vector
// ------------------------------------------------------------------------
int *psIVector::getIVector() 
{
  return Vec_;
}

// ************************************************************************
// clean
// ------------------------------------------------------------------------
void psIVector::clean()
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = 0;
}

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
psFVector::psFVector()
{
  length_ = 0;
  Vec_ = NULL;
}

// ************************************************************************
// Copy Constructor 
// ------------------------------------------------------------------------
psFVector::psFVector(const psFVector & v)
{
  length_ = v.length_;
  Vec_ = NULL;
  if (length_ > 0)
  {
    Vec_ = new float[length_];
    for (int ii = 0; ii < length_; ii++) Vec_[ii] = v.Vec_[ii];
  }
}

// ************************************************************************
// overload operator= 
// ------------------------------------------------------------------------
psFVector & psFVector::operator=(const psFVector & v)
{
  if (this == &v) return *this;
  delete [] Vec_;
  Vec_ = NULL;
  length_ = v.length_;
  if (length_ > 0)
  {
    Vec_ = new float[length_];
    for (int ii = 0; ii < length_; ii++) Vec_[ii] = v.Vec_[ii];
  }
  return *this;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
psFVector::~psFVector()
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = 0;
}

// ************************************************************************
// get length 
// ------------------------------------------------------------------------
int psFVector::length() 
{
  return length_;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int psFVector::load(psFVector &inVec)
{
  assert(this != &inVec);
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = inVec.length();
  if (length_ <= 0) return -1;
  Vec_ = new float[length_];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < length_; ii++) Vec_[ii] = inVec[ii];
  return 0;
}

// ************************************************************************
// set dimension
// ------------------------------------------------------------------------
int psFVector::setLength(int leng)
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  assert(leng > 0);
  length_ = leng;
  Vec_ = new float[leng];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < leng; ii++) Vec_[ii] = 0;
  return 0;
}

// ************************************************************************
// load vector
// ------------------------------------------------------------------------
int psFVector::load(int leng, float *data)
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  assert(leng > 0);
  assert(data != NULL);
  length_ = leng;
  Vec_ = new float[leng];
  assert(Vec_ != NULL);
  for (int ii = 0; ii < leng; ii++) Vec_[ii] = data[ii];
  return 0;
}

// ************************************************************************
// set entry
// ------------------------------------------------------------------------
int psFVector::setEntry(int ind, float fdata)
{
  if (ind < 0 || ind >= length_)
  {
    printf("psFVector setEntry ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  Vec_[ind] = fdata;
  return 0;
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
float psFVector::getEntry(int ind)
{
  if (ind < 0 || ind >= length_)
  {
    printf("psFVector getEntry ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  return Vec_[ind];
}

// ************************************************************************
// get entry
// ------------------------------------------------------------------------
float& psFVector::operator[](int ind) 
{
  if (ind < 0 || ind >= length_)
  {
    printf("psFVector operator[] ERROR: index = %d (%d)\n",ind+1,length_);
    exit(1);
  }
  return Vec_[ind];
}

// ************************************************************************
// get vector
// ------------------------------------------------------------------------
float *psFVector::getFVector() 
{
  return Vec_;
}

// ************************************************************************
// clean
// ------------------------------------------------------------------------
void psFVector::clean()
{
  if (Vec_ != NULL) delete [] Vec_;
  Vec_ = NULL;
  length_ = 0;
}

