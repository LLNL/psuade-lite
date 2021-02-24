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
// Definitions for the class AnnFunctions
// AUTHOR : Christopher K. Chan
// DATE   : 2005
// ************************************************************************
*/

/*************************************************************************
  Notes  : Functions directly manipulate SNNS-kernel. Must configure 
           SNNS-kernel with "--enable-enzo" to allow for krm_getNet() and
	   krm_putNet(). 
 ************************************************************************/
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>	

#ifdef HAVE_SNNS
#include "AnnFunctions.h"  /* Public function declarations */
#include "AnnFunctions.ph" /* Private function declarations */

#include <glob_typ.h> /* SNNS-Kernel Global Constants and data type defns */
#include <kr_ui.h>    /* SNNS-Kernel User-Interface function prototypes */
#include <cc_mac.h>   /* SNNS-Kernel Cascade-correlation parameter macros */
#include <kr_mem.h>	  /* SNNS-Kernel Memory manager */
#include <enzo_mem_typ.h> /* SNNS-Kernel Network data type */
#include <kr_typ.h>	  /* SNNS-Kernel Kernel data type definitions */
#include <kr_const.h>	  /* SNNS-Kernel Kernel constants */
#include <cc_glob.h>      /* SNNS-Kernel Cascade-correlation pointers */
#endif

/*****************************************************************************
  FUNCTION : checkErr  
  PURPOSE  : checks and prints errors with krui_ functions
******************************************************************************/
#ifdef HAVE_SNNS
static void checkErr(int err) 
{
   if (err != KRERR_NO_ERROR) 
   { 
      fprintf(stderr,"%s\n", krui_error(err));
      exit(1);
   }
}
#endif

/*****************************************************************************
  FUNCTION : makeUnits  
  PURPOSE  : makes all units in SNNS-kernel memory
******************************************************************************/

#ifdef HAVE_SNNS
static void makeUnits(int nInput, int nHidden) 
{
   struct PosType unitPos;
   int i, unitNo;
	
   /* Allocate SNNS-kernel memory for units */
   checkErr(krui_allocateUnits(nInput+nHidden+1));
	
   /*** Create input units ***/
   for (i = 0; i < nInput; i++)
   {
      unitNo = krui_createUnit("input", "Out_Identity", "Act_Identity",
                                 DEF_INIT_ACT,DEF_INIT_BIAS);
      if(unitNo < 0)
      {
         fprintf(stderr,"%s\n", krui_error(unitNo));
         exit(1);
      }
      krui_setUnitTType(unitNo, INPUT);		

      unitPos.x = 2;
      unitPos.y = i+2;
      unitPos.z = 0;
      krui_setUnitPosition(unitNo, &unitPos);
   }

   /*** Create hidden units ***/
   for (i = 0; i < nHidden; i++)
   {
      unitNo = krui_createUnit("hidden", "Out_Identity", DEF_ACT_FUNC,
                               DEF_INIT_ACT,DEF_INIT_BIAS);
      if(unitNo < 0)
      {
         fprintf(stderr,"%s\n", krui_error(unitNo));
         exit(1);
      }
      krui_setUnitTType(unitNo, HIDDEN);

      unitPos.x = 5;
      unitPos.y = i+2;
      unitPos.z = 0;		 
      krui_setUnitPosition(unitNo, &unitPos);		
   }
	
   /*** Create output unit ***/
   unitNo = krui_createUnit("output", "Out_Identity", "Act_Identity",
                            DEF_INIT_ACT,DEF_INIT_BIAS);
	
   if(unitNo < 0)
   {
      fprintf(stderr, "%s\n", krui_error(unitNo));
      exit(1);
   }
	
   krui_setUnitTType(unitNo, OUTPUT);
	
   if (nHidden == 0) { unitPos.x = 5; }
   else              { unitPos.x = 8; }
	
   unitPos.y = 2;
   unitPos.z = 0;	 
   krui_setUnitPosition(unitNo, &unitPos);
}
#endif

/*****************************************************************************
  FUNCTION : makeLinks 
  PURPOSE  : makes "full connections" between adjacent layers of units  
  NOTE	  : must call makeUnits() before makeLinks()
******************************************************************************/

#ifdef HAVE_SNNS
static void makeLinks(int nInput, int nHidden) 
{
   int inputUnitNo, hiddenUnitNo;
	
   if (nHidden == 0)
   {
      /*** Create links between input and output layers ***/
      krui_setCurrentUnit(nInput+1);
		
      for (inputUnitNo = 1; inputUnitNo <= nInput; inputUnitNo++)
         krui_createLink(inputUnitNo, 0.0);
   }
   else 
   {
      /*** Create links between input and hidden layers ***/
      for (hiddenUnitNo = 1; hiddenUnitNo <= nHidden; hiddenUnitNo++) 
      {
         /* Set target link to hidden unit */
         krui_setCurrentUnit(hiddenUnitNo+nInput);				

         for (inputUnitNo = 1; inputUnitNo <= nInput; inputUnitNo++)
            krui_createLink(inputUnitNo, 0.0);			
      }

      /*** Create links between hidden and output layers ***/
      krui_setCurrentUnit(nInput+nHidden+1);

      for (hiddenUnitNo = 1; hiddenUnitNo <= nHidden; hiddenUnitNo++) 
         krui_createLink(hiddenUnitNo+nInput, 0.0);
   }
}
#endif

/*****************************************************************************
  FUNCTION : ANNF_initialize  
  PURPOSE  : Initializes ANN class for <nNets> networks and <nPats> pattern 
  				 sets. 
  NOTES    : Must initialize before using ANNF_ functions.
******************************************************************************/

extern void ANNF_initialize(int nNets, int nPats)
{
#ifdef HAVE_SNNS
   if (netLoaded_)
   {
      fprintf(stderr,"ANNF_initialize ERROR - class already initialized. \n");
      exit(1);
   }
   if (nNets <= 0 || nPats < 0)
   {
      fprintf(stderr,"ANNF::ANNF_initialize ERROR - invalid inputs. \n");
      exit(1);		
   }

   /* Initialize class variables */	
   NET_array = (memNet *)malloc(sizeof(memNet) * nNets);	
   PAT_array = (int *)malloc(sizeof(int) * nPats);

   totalPats_ = nPats;
   totalNets_ = nNets;
   netLoaded_ = 1;
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_createNets  
  PURPOSE  : create feedforward network given number of input, hidden 
             and output units. Sets the network's update, learning and 
             activation functions according to the default values. 
******************************************************************************/

extern void ANNF_createNets(int nInput, int nHidden)
{
#ifdef HAVE_SNNS
   int i;

   if (!netLoaded_)
   {
      fprintf(stderr,"ANNF::ANNF_createNets ERROR - class not initialized. \n");
      exit(1);
   }
   if (nInput <= 0 || nHidden < 0)
   {
      fprintf(stderr,"ANNF::ANNF_createNets ERROR - invalid inputs. \n");
      exit(1);
   }

   /* Create network */
   for (i = 0; i < totalNets_; i++) 
   {
      makeUnits(nInput, nHidden);
      makeLinks(nInput, nHidden);

      /* Default learning and update functions */			
      checkErr(krui_setUpdateFunc(DEF_UPD_FUNC));
      checkErr(krui_setLearnFunc(DEF_LEARN_FUNC));
      krm_getNet(&NET_array[i]);
   }
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_saveNet  
  PURPOSE  : Saves network in NET_array[iNet] as fileName.  
******************************************************************************/

extern void ANNF_saveNet(int iNet, char *fileName) 
{
#ifdef HAVE_SNNS
   if (!netLoaded_)
   {
      fprintf(stderr,"ANNF::ANNF_saveNet ERROR - class not initialized. \n");
      exit(1);
   }
   if (iNet < 0 || iNet >= totalNets_)
   {
      fprintf(stderr,"ANNF::ANNF_saveNet ERROR - invalid iNet value. \n");
      exit(1);	
   }
   krm_putNet(&NET_array[iNet]);	
   checkErr(krui_saveNet(fileName, fileName));	
   krm_getNet(&NET_array[iNet]);
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_swapNet  
  PURPOSE  : Swap NET_array[iNet1] and NET_array[iNet2]
******************************************************************************/

extern void ANNF_swapNet(int iNet1, int iNet2)
{
#ifdef HAVE_SNNS
   memNet tempNet;

   if (!netLoaded_) 
   {
      fprintf(stderr,"ANNF::ANNF_swapNet ERROR - class not initialized. \n");
      exit(1);
   }
   if (iNet1 < 0 || iNet1 >= totalNets_)
   {
      fprintf(stderr,"ANNF::ANNF_swapNet ERROR - invalid iNet1 value. \n");
      exit(1);	
   }
   if (iNet2 < 0 || iNet2 >= totalNets_)
   {
      fprintf(stderr,"ANNF::ANNF_swapNet ERROR - invalid iNet2 value. \n");
      exit(1);	
   }
   if (iNet1 != iNet2)
   {
      tempNet = NET_array[iNet1];
      NET_array[iNet1] = NET_array[iNet2];
      NET_array[iNet2] = tempNet;
   }
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_createPat  
  PURPOSE  : creates a pattern set at PAT_array[iPat] using data from locations 
            <start> to <end> in arrays <aInput>/<aOutput>. Must correctly
            specify <nInput> units for function to work.				 			  				 
  NOTES    : Network at NET_array[0] must be created before using this function.
******************************************************************************/

extern void ANNF_createPat(int iPat, int iStart, int iEnd, int nInput, 
	                   int nSample, double *aInput, double *aOutput)
{
#ifdef HAVE_SNNS
   int i, j, totalUnits;

   if (!netLoaded_)
   {
      fprintf(stderr,"ANNF::ANNF_createPat ERROR - class not initialized. \n");
      exit(1);
   }
   if (aInput == NULL || aOutput == NULL)
   {
      fprintf(stderr,"ANNF::ANNF_createPat ERROR - invalid input arrays. \n");
      exit(1);	
   }
   if (nInput <= 0 || iStart < 0)
   {
      fprintf(stderr,"ANNF::ANNF_createPat ERROR - invalid inputs. \n");
      exit(1);
   }
   if (iPat < 0 || iPat >= totalPats_)
   {
      fprintf(stderr,"ANNF::ANNF_createPat ERROR - invalid iPat value. \n");
      exit(1);	
   }
   
   /* Uses network at NET_array[0] as topology template */
   krm_putNet(&NET_array[0]);
   totalUnits = krui_getNoOfUnits();
	
   /* Create pattern set  */
   checkErr(krui_allocNewPatternSet(&PAT_array[iPat]));
   if (iStart < iEnd)
   {
      for (i = iStart; i <= iEnd; i++)
      {
         for (j = 1; j <= nInput; j++)
            krui_setUnitActivation(j, aInput[i*nInput+j-1]); 		
         krui_setUnitActivation(totalUnits, aOutput[i]);
         checkErr(krui_newPattern());
      }
   }
	
   /* When iStart > iEnd, holdout sampling is being used*/
   else
   {
      for (i = 0; i <= iEnd; i++)
      {			
         for (j = 1; j <= nInput; j++)
            krui_setUnitActivation(j, aInput[i*nInput+j-1]); 
         krui_setUnitActivation(totalUnits, aOutput[i]);
         checkErr(krui_newPattern());
      }
      for (i = iStart; i < nSample; i++)
      {			
         for (j = 1; j <= nInput; j++)
            krui_setUnitActivation(j, aInput[i*nInput+j-1]); 
         krui_setUnitActivation(totalUnits, aOutput[i]);
         checkErr(krui_newPattern());
      }
   }
   krm_getNet(&NET_array[0]);	
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_createPatRand  
  PURPOSE  : Same as ANNF_createPat but patterns are randomly selected (with 
             replacement) from data arrays.   			  				 
  NOTES    : Network at NET_array[0] must be created before using this function.
             Random generator seed set according to computer time.
******************************************************************************/

extern void ANNF_createPatRand(int iPat, int iStart, int iEnd, int nInput, 
	                       double *aInput, double *aOutput)
{
#ifdef HAVE_SNNS
   int i, j, r, totalUnits; 

   if (!netLoaded_)
   {
      fprintf(stderr,"ANNF_createPatRand ERROR - class not initialized. \n");
      exit(1);
   }
   if (aInput == NULL || aOutput == NULL)
   {
      fprintf(stderr,"ANNF_createPatRand ERROR - invalid input arrays.\n");
      exit(1);	
   }
   if (nInput <= 0 || iStart < 0)
   {
      fprintf(stderr,"ANNF::ANNF_createPatRand ERROR - invalid inputs. \n");
      exit(1);
   }
   if (iPat < 0 || iPat >= totalPats_) 
   {
      fprintf(stderr,"ANNF_createPatRand ERROR - invalid iPat value. \n");
      exit(1);	
   }

   srand(time(NULL));
	
   /* Uses network at NET_array[0] as topology template */	
   krm_putNet(&NET_array[0]);
   totalUnits = krui_getNoOfUnits();

   /* Generate bootstrap sample */
   checkErr(krui_allocNewPatternSet(&PAT_array[iPat]));
   for (i = iStart; i <= iEnd; i++)
   {
      r = rand() % (iEnd-iStart+1);
      for (j = 1; j <= nInput; j++)
         krui_setUnitActivation(j, aInput[r*nInput+j-1]); 
      krui_setUnitActivation(totalUnits, aOutput[r]);
      checkErr(krui_newPattern());
   }

   krm_getNet(&NET_array[0]);
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_savePat  
  PURPOSE  : Saves pattern set at PAT_array[iPat] as fileName.
 ****************************************************************************/

extern void ANNF_savePat(int iPat, char *fileName) 
{
#ifdef HAVE_SNNS
   if (iPat < 0 || iPat >= totalPats_)
   {
      fprintf(stderr,"ANNF::ANNF_savePat ERROR - invalid iPat value. \n");
      exit(1);
   }
   checkErr(krui_saveNewPatterns(fileName, PAT_array[iPat]));
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_trainNet  
  PURPOSE  : Trains neural network with given network location, training 
             pattern set and learning algorithm. Returns the testing set MSE 
            in <aResults[0]> and final number of hidden units in <aResults[1]>. 
            <ccFlag>==1 indicates that CasCor algorithm will be used so all
            previous special and hidden units will be deleted from network
            to be trained.				   			  				 
  NOTES    : Assumes testing pattern set is located at PAT_array[0].
******************************************************************************/

extern void ANNF_trainNet(int iNet, int iPat, char *learnFuncName,  
	                  float *aLearnFuncParam, int nLearnFuncParam,
	                  float **aResults, int ccFlag)
{
#ifdef HAVE_SNNS
   int subPatISize[5] = {0,0,0,0,0}, subPatOSize[5] = {0,0,0,0,0};
   int subPatIStep[5] = {0,0,0,0,0}, subPatOStep[5] = {0,0,0,0,0};
   int i, j, nParamOut, UP1_stop = -1, UP2_stop = -1, iStep; 
   float localMinMSE = 0, currMSE = 0, lastMSEs[2] = {1E16,1E16};	
   float *aParamOut, *aPtr = *aResults;
	
   if (!netLoaded_)
   {
      fprintf(stderr,"ANNF::ANNF_trainNet - class not initialized. \n");
      exit(1);
   }	
   if (iNet < 0 || iNet >= totalNets_)
   {
      fprintf(stderr,"ANNF::ANNF_trainNet - invalid iNet value. \n");
      exit(1);	
   }
   if (iPat < 0 || iPat >= totalPats_) 
   {
      fprintf(stderr,"ANNF::ANNF_trainNet - invalid iPat value. \n");
      exit(1);	
   }

   /* Initialize weights */
   krm_putNet(&NET_array[iNet]);		
   krui_setCurrPatSet(PAT_array[iPat]);
   krui_setSeedNo(0);	
   krui_shufflePatterns(TRUE);
   krui_DefTrainSubPat(subPatISize,subPatOSize,subPatIStep,subPatOStep,NULL);
	
   /* Initialize network according to ccFlag */
   if (ccFlag == ANNF_CC)
   {
      /* Delete previously trained network if performing RE-train. 
         Does nothing otherwise. */
      cc_deleteAllSpecialAndAllHiddenUnits();		
		
      checkErr(krui_setInitialisationFunc("CC_Weights"));
      checkErr(krui_initializeNet(DEF_INIT_WEIGHTS, 2));
      checkErr(krui_setUpdateFunc("CC_Order"));
      checkErr(krui_setLearnFunc("CC"));
      iStep = DEF_CC_STEP;
   }
   else if (ccFlag == ANNF_NOCC)
   {
      checkErr(krui_setInitialisationFunc("Randomize_Weights"));
      checkErr(krui_initializeNet(DEF_INIT_WEIGHTS, 2));
      checkErr(krui_setUpdateFunc("Topological_Order"));
      checkErr(krui_setLearnFunc(learnFuncName));			
      iStep = DEF_NOCC_STEP;
   }
   else
   {
      fprintf(stderr, "ANNF::ANNF_trainNet - invalid ccFlag");
      exit(1);	
   }

   /* Train network until stopping criteria satisfied */
   for (i = 1; i <= DEF_MAX_EPOCHS; i++)
   {
      krui_learnAllPatterns(aLearnFuncParam, nLearnFuncParam, 
                            &aParamOut, &nParamOut);

      if (i % iStep == 0)
      {
         /* Test network accuracy */
         krui_setCurrPatSet(PAT_array[0]);
         krui_shufflePatterns(FALSE);		
         krui_DefTrainSubPat(subPatISize,subPatOSize,subPatIStep,
                             subPatOStep,NULL);
         krui_testAllPatterns(aLearnFuncParam, nLearnFuncParam, &aParamOut, 
                              &nParamOut);
         currMSE = aParamOut[0]/krui_getNoOfPatterns();

         /* check if NaN error encountered */
         if (isnan(aParamOut[0]))
         { 
            if (verboseFlag_ == 1)
            {
               printf("*..Epoch %d: validation MSE: %e \n", i, currMSE);
               printf("*..Stopped training due to floating-point exception.\n");
            }
            break;
         }

         /* Stopping Criteria 1: error increases after two successive times */					
         if (currMSE > lastMSEs[1])
         {
            if (UP1_stop == i-iStep && UP2_stop == i-(2*iStep))
            { 
               if (verboseFlag_ == 1)
               {
                  printf("*...Epoch %d: validation MSE: %e \n", i, currMSE);
                  printf("*...Stopped training due to Criteria 1 ");
                  printf("(successive MSE increase). \n");
               }
               break;
            }
            else if (UP1_stop == i-iStep) { UP2_stop = UP1_stop; }
            UP1_stop = i;
         }

         /* Stopping Criteria 2: progress is less than 0.1 and over 10 epochs */
         localMinMSE = currMSE;
         for (j = 0; j < 2; j++)
         {
            if (lastMSEs[j] < localMinMSE)
            localMinMSE = lastMSEs[j];
         }

         if ((lastMSEs[0]+lastMSEs[1]+currMSE)/(3*localMinMSE)<1.01 && i>12)
         {
            if (verboseFlag_ == 1)
            {
               printf("*...Epoch %d: validation MSE: %e \n", i, currMSE);
               printf("*...Stopped training due to Criteria 2 ");
               printf("(lack of progress). \n");
            }
            break;
         }

         if (verboseFlag_ == 1) 
            printf("*...Epoch %d: validation MSE: %e \n", i, currMSE);
         /* Stores the previous two validation MSEs */
         lastMSEs[0] = lastMSEs[1];
         lastMSEs[1] = currMSE;

         /* Reset network to training pattern */
         krui_setCurrPatSet(PAT_array[iPat]);
         krui_shufflePatterns(TRUE);
         krui_DefTrainSubPat(subPatISize,subPatOSize,subPatIStep,
                             subPatOStep,NULL);
      }
   }
		
   aPtr[0] = currMSE;
   aPtr[1] = krui_getNoOfUnits();
   krm_getNet(&NET_array[iNet]);				
#else
   printf("ANN ERROR: SNN not installed.\n");
   exit(1);
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_evaluatePoint  
  PURPOSE  : Return output of network at NET_array[iNet] given input values
******************************************************************************/  
extern float ANNF_evaluatePoint(int iNet, float *xPoint)
{
#ifdef HAVE_SNNS
   register struct Unit *UnitPtr;
   register Patterns patPtr; 	
   float result;
   int i;
		
   if (!netLoaded_)
   {
      fprintf(stderr,"ANNF_evaluatePoint ERROR - no network in memory. \n");
      exit(1);
   }
   if (iNet < 0 || iNet >= totalNets_)
   {
      fprintf(stderr,"ANNF::ANNF_evaluatePoint ERROR - invalid iNet value. \n");
      exit(1);	
   }
	
   krm_putNet(&NET_array[iNet]);	
	
   /********* PROPAGATE INPUT DATA POINTS FORWARD **********/
   patPtr = (Patterns) xPoint;   
   cc_setPointers();	/* Updates CasCor algorithm's internal pointers */
	
   /* Propagate activations through input layer */
   for (UnitPtr = *FirstInputUnitPtr, i=0; UnitPtr != NULL;
        UnitPtr = FirstInputUnitPtr[++i])
   {
      UnitPtr->Out.output = ((UnitPtr->out_func == OUT_IDENTITY) ? \
            (UnitPtr->act = *patPtr++) : \
            (*UnitPtr->out_func) (UnitPtr->act = *patPtr++));			
   }
	
   /* Propagate activations through hidden layer */
   for(UnitPtr = *FirstHiddenUnitPtr, i=0; UnitPtr != NULL;
       UnitPtr = FirstHiddenUnitPtr[++i])
   {
      UnitPtr->Out.output = ((UnitPtr->out_func == OUT_IDENTITY) ? \
         (UnitPtr->act = (*UnitPtr->act_func)(UnitPtr)) : \
         (*UnitPtr->out_func) (UnitPtr->act = (*UnitPtr->act_func)(UnitPtr)));
   }

   /* Propagate activations through output layer */
   for(UnitPtr = *FirstOutputUnitPtr,i=0; UnitPtr != NULL; 
       UnitPtr = FirstOutputUnitPtr[++i])
   {		
      UnitPtr->Out.output = ((UnitPtr->out_func == OUT_IDENTITY) ? 
          (UnitPtr->act = (*UnitPtr->act_func)(UnitPtr)) : 
          (*UnitPtr->out_func) (UnitPtr->act = (*UnitPtr->act_func)(UnitPtr))); 
   }

   /* Look at resulting output */
   UnitPtr = *FirstOutputUnitPtr;
   result = (float)UnitPtr->Out.output;

   krm_getNet(&NET_array[iNet]);
   return result;
#else
   return -1;
#endif
}

/*****************************************************************************
  FUNCTION : ANNF_setVerboseFlag 
  PURPOSE  : Prints to stdout additional information.
******************************************************************************/

extern void ANNF_setVerboseFlag(void) 
{
#ifdef HAVE_SNNS
   verboseFlag_ = 1;
#endif
}


/*****************************************************************************
  FUNCTION : ANNF_exit  
  PURPOSE  : Purge SNNS-Kernel and AnnFunctions class data from memory. Does
             nothign if class not initialized.
******************************************************************************/

extern void ANNF_exit(void) 
{
#ifdef HAVE_SNNS
   int i;
	
   if (netLoaded_)
   {
      /* Free stored networks */
      for (i = 0; i < totalNets_; i++)
      {
         krm_putNet(&NET_array[i]);
         krui_deleteNet();
      }

      /* Free stored pattern sets */
      for (i = totalPats_ -1; i >= 0; i--) {
         checkErr(krui_deletePatSet(PAT_array[i]));		
      }

      /* Free network and pattern array pointers */
      free(NET_array);
      free(PAT_array);
	
      /* Reset class variables */
      verboseFlag_ = 0;
      totalPats_ = 0;
      totalNets_ = 0;
      netLoaded_ = 0;
   }
#endif
}

/* end of file */
