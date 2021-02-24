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
// Functions for the class PsuadeBase
// AUTHOR : CHARLES TONG
// DATE   : 2012
// ************************************************************************
#include <stddef.h>
#include <Globals.h>

// -------------------------------------------------------------------------
// global variables
// -------------------------------------------------------------------------
int          psIOExpertMode_=0;
int          psRSExpertMode_=0;
int          psRSCodeGen_=0;
int          psSamExpertMode_=0;
int          psAnaExpertMode_=0;
int          psOptExpertMode_=0;
int          psPDFDiagMode_=0;
int          psMasterMode_=0;
int          psInteractive_=0;
int          psScreenOutput_=0;
int          psLibMode_=1;
int          psGMMode_=0;
int          psPythonOverride_=0;
int          psPlotTool_=0;
long         psRandomSeed_=-1;
int          psFAMaxDataPts_=10000;
int          psConstraintSetOp_=0;
char         *psConfigFileName_=NULL;
char         *psPythonInterpreter_=NULL;
PsuadeConfig *psConfig_=NULL;
CommManager  *psCommMgr_=NULL;
const char   *psInputFilename_  = "psuadeData";
const char   *psOutputFilename_ = "psuadeData";

