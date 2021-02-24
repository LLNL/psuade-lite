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
// Functions for the class Acosso
// Reference: This class uses Curt Storlie's (LANL) Acosso method (in R)
// DATE   : 2015
// ************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "Acosso.h"
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "Psuade.h"

#define PABS(x) ((x) > 0 ? (x) : (-(x)))

// ************************************************************************
// Constructor for object class Acosso
// ------------------------------------------------------------------------
Acosso::Acosso(int nInputs,int nSamples) : FuncApprox(nInputs,nSamples)
{
   printf("Acosso ERROR: not included in this release\n");
   exit(1);
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
Acosso::~Acosso()
{
}

// ************************************************************************
// initialize
// ------------------------------------------------------------------------
int Acosso::initialize(double *, double *)
{
   return 0;
}

// ************************************************************************
// Generate results for display
// ------------------------------------------------------------------------
int Acosso::genNDGridData(double *, double *, int *, double **, double **)
{
   return 0;
}

// ************************************************************************
// Generate 1D results for display
// ------------------------------------------------------------------------
int Acosso::gen1DGridData(double *, double *, int, double *, int *, 
                          double **, double **)
{
   return 0;
}

// ************************************************************************
// Generate 2D results 
// ------------------------------------------------------------------------
int Acosso::gen2DGridData(double *, double *, int, int, double *, int *, 
                          double **, double **)
{
   return 0;
}

// ************************************************************************
// Generate 3D results for display
// ------------------------------------------------------------------------
int Acosso::gen3DGridData(double *, double *, int, int, int, double *, 
                          int *, double **, double **)
{
   return 0;
}

// ************************************************************************
// Generate 4D results for display
// ------------------------------------------------------------------------
int Acosso::gen4DGridData(double *, double *, int, int, int , int , 
                          double *, int *, double **, double **)
{
   return 0;
}

// ************************************************************************
// Evaluate a given point 
// ------------------------------------------------------------------------
double Acosso::evaluatePoint(double *)
{
   double Y=0.0;
   return Y;
}

// ************************************************************************
// Evaluate a number of points 
// ------------------------------------------------------------------------
double Acosso::evaluatePoint(int, double *, double *)
{
   return 0.0;
}

// ************************************************************************
// Evaluate a given point, return also the standard deviation 
// ------------------------------------------------------------------------
double Acosso::evaluatePointFuzzy(double *, double &)
{
   return 0.0;
}

// ************************************************************************
// Evaluate a number of points, return also the standard deviation 
// ------------------------------------------------------------------------
double Acosso::evaluatePointFuzzy(int,double *, double *,double *)
{
   return 0.0;
}

// ************************************************************************
// run Storlie's Acosso code
// ------------------------------------------------------------------------
int Acosso::runAcosso(int, double *, double *, int, int *)
{
   return 0;
}

// ************************************************************************
// generate Acosso code
// ------------------------------------------------------------------------
int Acosso::genAcosso()
{
   return 0;
}

// ************************************************************************
// set parameters
// ------------------------------------------------------------------------
double Acosso::setParams(int, char **)
{
   return 0.0;
}

