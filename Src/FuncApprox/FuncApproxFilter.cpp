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
// Functions for the class FuncApproxFilter  
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "sysdef.h"
#include "PsuadeData.h"
#include "FuncApproxFilter.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------
FuncApproxFilter::FuncApproxFilter(char *filename)
{
   char       inputLine[500], keyword[500];
   FILE*      fp;
   PsuadeData *ioPtr;

   fp = fopen(filename, "r");
   if (fp == NULL)
   {
      printf("FuncApproxFilter ERROR - no data file. \n");
      assert(fp != NULL);
   }
   fgets(inputLine, 500, fp);
   sscanf(inputLine, "%s", keyword);
   while (keyword[0] == '#')
   {
      fgets(inputLine, 500, fp);
      sscanf(inputLine, "%s", keyword);
   }
   if (strcmp(keyword, "PSUADE_IO"))
   {
      printf("FuncApproxFilter ERROR - PSUADE_IO section absent.\n");
      printf("                         filename = %s.\n",filename);
      exit(1);
   }
   fclose(fp);

   printLevel_ = 0; 
   strcpy(dataFileName_, filename);
   YLBound_ = - PSUADE_UNDEFINED;
   YUBound_ =   PSUADE_UNDEFINED;

   ioPtr = new PsuadeData();
   ioPtr->readPsuadeFile(dataFileName_);
   faPtr_ = genFAInteractive(ioPtr, 2);
   delete ioPtr;
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
FuncApproxFilter::~FuncApproxFilter()
{
   delete faPtr_;
}

// ************************************************************************
// Set print level
// ------------------------------------------------------------------------
int FuncApproxFilter::setPrintLevel(int level)
{
   printLevel_ = level;
   return 0;
}

// ************************************************************************
// Set output bounds 
// ------------------------------------------------------------------------
int FuncApproxFilter::setYBounds( double& lower, double& upper )
{
   YLBound_ = lower;
   YUBound_ = upper;
   return 0;
}

// ************************************************************************
// Get output bounds 
// ------------------------------------------------------------------------
int FuncApproxFilter::getYBounds( double &lower, double &upper )
{
   lower = YLBound_;
   upper = YUBound_;
   return 0;
}

// ************************************************************************
// evaluate output given a sample point
// ------------------------------------------------------------------------
double FuncApproxFilter::evaluatePoint(double *X, int &flag)
{
   double Y = faPtr_->evaluatePoint(X);
   flag = 1;
   if (Y < YLBound_ || Y > YUBound_) flag = 0;
   return Y;
}

