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
// Functions for the base class Analyzer
// AUTHOR : CHARLES TONG
// DATE   : 2004
// ************************************************************************
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Analyzer.h"

// ************************************************************************
// constructor
// ------------------------------------------------------------------------
Analyzer::Analyzer()
{
   strcpy(analyzerName_, "NONE");
   rstype_ = -1;
#ifdef HAVE_PYTHON
   AnalysisDataDict = PyDict_New();
#endif
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
Analyzer::~Analyzer()
{
#ifdef HAVE_PYTHON
   Py_DECREF(AnalysisDataDict);
#endif
}

// ************************************************************************
// set name 
// ------------------------------------------------------------------------
int Analyzer::setName(const char *name)
{
   strncpy(analyzerName_, name, 100);

#ifdef HAVE_PYTHON
   PyObject* s = PyString_FromString(name);
   PyDict_SetItemString(AnalysisDataDict, "name", s);
   Py_DECREF(s);
#endif

   return 0;
}

// ************************************************************************
// analyze function
// ------------------------------------------------------------------------
double Analyzer::analyze(aData &adata)
{
   (void) adata;
   return -1.0;
}

// ************************************************************************
// set parameter function 
// ------------------------------------------------------------------------
int Analyzer::setParams(int, char **)
{
   printf("Analyzer::setParams ERROR: should not be called.\n");
   exit(1);
   return 0;
}

