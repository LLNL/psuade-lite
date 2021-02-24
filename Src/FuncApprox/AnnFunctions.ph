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
   Definitions for the class ANN
// AUTHOR : Christopher K. Chan
// DATE   : 2005
// ************************************************************************
*/

#ifndef __ANNFUNCTIONSPH__
#define __ANNFUNCTIONSPH__

#include <glob_typ.h>	  /* SNNS-Kernel Global Constants, data type defns */
#include <kr_typ.h>       /* SNNS-Kernel Kernel data type definitions */
#include <enzo_mem_typ.h> /* SNNS-Kernel memNet definitions */

/* default values */
static char *DEF_ACT_FUNC = "Act_Logistic"; /* hidden unit activation function */
static int   DEF_INIT_BIAS = 0.0;           /* initial bias */
static int   DEF_INIT_ACT = 0.0;            /* initial activation */
static float DEF_INIT_WEIGHTS[2] = {0.1, -0.1};	/* range of init weights */
static char *DEF_LEARN_FUNC = "CC";         /* learning function */
static char *DEF_UPD_FUNC = "CC_Order";	    /* update function */

static const int DEF_MAX_EPOCHS = 200;	    /* max learning cycles */
static const int DEF_CC_STEP = 2;
static const int DEF_NOCC_STEP = 11;        /* cycles before validation  
                                               for CC and non-CC algorthms */

/* class variables */
static int netLoaded_ = 0;                  /* current situation of class */
static int verboseFlag_ = 0;                /* verbose option */

/* Maintain pattern set array */
static int *PAT_array;
static int totalPats_ = 0;

/* Maintain NET array */
static memNet *NET_array;
static int totalNets_ = 0;

/* Function definitions */
static void checkErr(int err); 
static void makeUnits(int nInput, int nHidden);
static void makeLinks(int nInput, int nHidden);

#endif 
