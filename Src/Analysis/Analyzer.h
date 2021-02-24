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
// Definition for the base class Analyzer
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#ifndef __ANALYZERH__
#define __ANALYZERH__
                                                                                
#ifdef HAVE_PYTHON
#include <Python.h>
#endif

#include "aData.h"

// ************************************************************************
// class definition
// ************************************************************************
class Analyzer
{
public:
   char analyzerName_[100];
   int  rstype_;

public:
                                                                                
   Analyzer();
   virtual ~Analyzer();

   int     setName(const char *);
   virtual double analyze(aData &);
   virtual int    setParams(int, char **);

#ifdef HAVE_PYTHON
   PyObject *AnalysisDataDict;
#endif
};

#endif // __ANALYZERH__

