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
// Definition for the Global Variables class
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __GLOBALSH__ // header guards
#define __GLOBALSH__

#include "CommManager.h"
#include "PsuadeConfig.h"

extern int          psIOExpertMode_;
extern int          psRSExpertMode_;
extern int          psRSCodeGen_;
extern int          psSamExpertMode_;
extern int          psAnaExpertMode_;
extern int          psOptExpertMode_;
extern int          psPDFDiagMode_;
extern int          psMasterMode_;
extern int          psInteractive_;
extern int          psScreenOutput_;
extern int          psLibMode_;
extern int          psGMMode_;
extern int          psPythonOverride_;
extern int          psPlotTool_;
extern long         psRandomSeed_;
extern int          psFAMaxDataPts_;
extern int          psConstraintSetOp_;
extern char         *psConfigFileName_;
extern char         *psPythonInterpreter_;
extern PsuadeConfig *psConfig_;
extern CommManager  *psCommMgr_;
extern const char   *psInputFilename_;
extern const char   *psOutputFilename_;

#endif

