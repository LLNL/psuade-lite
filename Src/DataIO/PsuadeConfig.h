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
// Definition for the class PsuadeConfig
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#ifndef __PSUADECONFIGH__
#define __PSUADECONFIGH__

// ************************************************************************
// ************************************************************************
// main class definition 
// ************************************************************************
// ************************************************************************

class PsuadeConfig 
{
   int  printLevel_;
   int  nLines_;
   int  nLinesUsed_;
   char **fileData_;

public:

   // Constructor
   // fname : configure file name
   // printLevel : diagnostics print level
   PsuadeConfig(char *fname, int printLevel);

   // Constructor
   PsuadeConfig();

   // Copy Contructor by Bill Oliver
   PsuadeConfig(const PsuadeConfig & ps);

   // Destructor
   ~PsuadeConfig();

   // Get parameter 
   char *getParameter(const char *);

   // put parameter 
   void putParameter(const char *);

   // remove parameter 
   int removeParameter(const char *);

   // write file
   // fname : configure file name
   void writeToFile(char *fname);

   // add from file
   // fname : configure file name
   void addFromFile(char *fname);

   // print 
   void print();
};

// ************************************************************************
// friend function
// ------------------------------------------------------------------------
extern "C"
{
   int genConfigFileTemplate(char *filename);
}

#endif // __PSUADECONFIGH__

