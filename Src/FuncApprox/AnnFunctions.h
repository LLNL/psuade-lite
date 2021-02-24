/* ************************************************************************
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
// Public Functions for the class AnnFunctions
// Definitions for the class ANN
// AUTHOR : Christopher K. Chan
// DATE   : 2005
// ************************************************************************
*/

#ifndef __ANNFUNCTIONSH__
#define __ANNFUNCTIONSH__

#define ANNF_CC	  1
#define ANNF_NOCC 2

extern void ANNF_initialize(int nNets, int nPats);
extern void ANNF_createNets(int nInput, int nHidden);
extern void ANNF_saveNet(int iNet, char *fileName);
extern void ANNF_swapNet(int iNet1, int iNet2);

extern void ANNF_createPat(int iPat, int iStart, int iEnd, int nInput,
                           int nSample, double *aInput, double *aOutput);
extern void ANNF_createPatRand(int iPat, int iStart, int iEnd, int nInput,
			       double *aInput, double *aOutput);
extern void ANNF_savePat(int iPat, char *fileName);

extern void ANNF_trainNet(int iNet, int iPat, char *learnFuncName, 
			  float *aLearnFuncParam, int nLearnFuncParam,
			  float **aResults,int ccFlag);									
extern float ANNF_evaluatePoint(int iNet, float *xPoint);

extern void ANNF_setVerboseFlag(void);
extern void ANNF_exit(void);							  
#endif
