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
#ifndef __VECTORH__
#define __VECTORH__

#ifdef WINDOWS
#undef min
#undef max
#endif

/**
 * @name psVector class
 *
 **/
/*@{*/

// ************************************************************************
// class definition
// ************************************************************************
class psVector
{
   int    length_;
   double *Vec_;

public:

   psVector();
   // Copy Constructor by Bill Oliver
   psVector(const psVector & v);

   // Overload operator= by Bill Oliver
   psVector & operator=(const psVector & v);

   ~psVector();
   int     length();
   int     load(psVector &);
   int     setLength(int);
   int     load(int, double *);
   int     setEntry(int, double);
   double  getEntry(int);
   int     addElements(int, double *);
   double& operator[](int ind);
   double  max();
   double  min();
   double  sum();
   void    scale(double);
   void    sort();
   void    subvector(int, int);
   void    clean();
   void    print(char *);
   double  *getDVector();
};

// ------------------------------------------------------------------------
// integer vector
// ------------------------------------------------------------------------
class psIVector
{
   int length_;
   int *Vec_;

public:

   psIVector();
   psIVector(const psIVector & v);
   psIVector & operator=(const psIVector & v);

   ~psIVector();
   int length();
   int load(psIVector &);
   int setEntry(int, int);
   int getEntry(int);
   int setLength(int);
   int load(int, int *);
   int & operator[](int ind);
   int sum();
   void subvector(int, int);
   int *getIVector();
   void clean();
};

// ------------------------------------------------------------------------
// float vector
// ------------------------------------------------------------------------
class psFVector
{
   int   length_;
   float *Vec_;

public:

   psFVector();
   psFVector(const psFVector & v);
   psFVector & operator=(const psFVector & v);

   ~psFVector();
   int   length();
   int   load(psFVector &);
   int   setEntry(int, float);
   float getEntry(int);
   int   setLength(int);
   int   load(int, float *);
   float & operator[](int ind);
   float *getFVector();
   void  clean();
};

/*@}*/

#endif /* __VECTORH__ */

