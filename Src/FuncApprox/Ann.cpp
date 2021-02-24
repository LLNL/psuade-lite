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
// Functions for the class ANN
// AUTHOR : Christopher K. Chan 
// DATE   : July, 2005
// ************************************************************************
// ************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sysdef.h"
#include "Ann.h"

#ifdef HAVE_SNNS
extern "C" {
#include "AnnFunctions.h"
typedef float   FlintType;
#include "cc_mac.h"
}
#endif

#define MAX_FLOAT 1E16
#define ANN_BOOTSTRAP 1
#define ANN_HOLDOUT 2
#define ANN_ALLOUT 3

#define PABS(x) (((x) > 0.0) ? (x) : -(x))

//****************************************************************************
// FUNCTION : ANN(int nInput, int nSample) 
//****************************************************************************
Ann::Ann(int nInputs, int nSamples) : FuncApprox(nInputs, nSamples)
{	
#ifdef HAVE_SNNS
   this->nInputs_ = nInputs;
   this->nSamples_ = nSamples;
   faID_ = PSUADE_RS_ANN;
   learnFuncName_ = NULL;
   inputDataMax_  = NULL;
   inputDataMin_  = NULL;
   inputData_     = NULL;


   nNet_ = 25;
   nHidden_ = 0;	

   samplingMethod_ = ANN_ALLOUT;
   normFlag_ = true;
   normBounds_[0] = -1.0;
   normBounds_[1] = 1.0;	

   learnFuncName_ = (char *)malloc(sizeof(char) * 20);
   strcpy(learnFuncName_, "CC");
   nLearnFuncParam_ = 28;
   nMaxRetrain_ = 2;	
   trainPerformed_ = false;

   maxMagnitudeOutput_ = 0;
   accuracyLevel_ = 2;

   verboseFlag_ = false;
   saveNetFlag_ = false;
   savePatFlag_ = false;		

   scalingFactor_ = 1.0;

   aLearnFuncParam_[0] = 0.00045;  // Learning function parameters
   aLearnFuncParam_[1] = 1.75;
   aLearnFuncParam_[2] = 0.00005;
   aLearnFuncParam_[3] = 0.00045;
   aLearnFuncParam_[4] = 1.75;
   aLearnFuncParam_[7] = QUICKPROP;	  

   aLearnFuncParam_[6] = 0.0001;   // Max Output Unit Error
   aLearnFuncParam_[8] = OFF;	   // Print Covariance Data On/Off 
	
   aLearnFuncParam_[9] = 0.04;	   // Min Error Change (Candidate Units)
   aLearnFuncParam_[10] = 50;	   // Patience (Candidate Units)
   aLearnFuncParam_[11] = 200;	   // Max Update Epochs (Candidate Units)
   aLearnFuncParam_[12] = 8;	   // Max Candidate Units

   aLearnFuncParam_[13] = ASYM_SIGMOID;	// Activation Function
	
   aLearnFuncParam_[14] = 0.00;	   // Min Error Change (Output Units)
   aLearnFuncParam_[15] = 50;	   // Patience (Output Units)
   aLearnFuncParam_[16] = 200;	   // Max Update Epochs (Output Units)
	
   aLearnFuncParam_[17] = OFF;	   // Pruning On/Off

   aLearnFuncParam_[18] = 0;	   // 18 and 19 not used by SNNSv4.2 
   aLearnFuncParam_[19] = 0;				

   aLearnFuncParam_[20] = SBC;	   // Model Selection Function
   aLearnFuncParam_[21] = CC_SDCC; // CC Modified Algorithm 

   aLearnFuncParam_[22] = 0; // CC Modified Algorithm Parameters 
   aLearnFuncParam_[23] = 0;
   aLearnFuncParam_[24] = 0;
   aLearnFuncParam_[25] = 0;
   aLearnFuncParam_[26] = 0;

   aLearnFuncParam_[27] = ON; // Caching Unit Activations On/Off

   readConfigFile();	
   checkClassVars();

   if (verboseFlag_ == true)
   {
      printf("\n");
      printf("* ============================================== *\n");
      printf("* Loaded the following settings from Ann.config: \n");
      printf("* Number of Nets: %d \n", nNet_);
      printf("* Input Units: %d  Hidden Units: %d \n", nInputs_, nHidden_);

      if (samplingMethod_ == ANN_HOLDOUT)
         printf("* Sampling Method: Holdout \n");
      else 
         printf("* Sampling Method: Bootstrap \n");

      if (normFlag_ == true)
      {
         printf("* Normalization: yes \n");			
	 printf("* Normalization Bounds: [%1.2f, %1.2f] \n", 
                normBounds_[0], normBounds_[1]);
      }
		
      printf("* Training Algorithm: %s \n", learnFuncName_);
      printf("* Training Parameters: %d ", nLearnFuncParam_);
   }
#else
   printf("ERROR::SNNS not installed.\n");
   exit(1);
#endif
}

//****************************************************************************
// Copy Constructor Created by Bill Oliver
//****************************************************************************
Ann::Ann(const Ann & an) : FuncApprox(an.nInputs_, an.nSamples_)  
{	
   nNet_ = an.nNet_;
   nHidden_ = an.nHidden_;	

   samplingMethod_ = an.samplingMethod_;
   normFlag_ = an.normFlag_;
   normBounds_[0] = an.normBounds_[0];
   normBounds_[1] = an.normBounds_[1];	

   learnFuncName_ = (char *)malloc(strlen(an.learnFuncName_) + 1);
   if(learnFuncName_ == NULL)
   {
      printf("Out of Heap Memory in file %s line %d aborting.\n", 
             __FILE__, __LINE__);
      exit(1);
   }
   strcpy(learnFuncName_, an.learnFuncName_);
   nLearnFuncParam_ = an.nLearnFuncParam_;
   nMaxRetrain_ =an.nMaxRetrain_;	
   trainPerformed_ = an.trainPerformed_;

   maxMagnitudeOutput_ = an.maxMagnitudeOutput_;
   accuracyLevel_ = an.accuracyLevel_;

   verboseFlag_ = an.verboseFlag_;
   saveNetFlag_ = an.saveNetFlag_;
   savePatFlag_ = an.savePatFlag_;		

   scalingFactor_ = an.scalingFactor_;

   for(int i = 0; i < 28; i++)
   {
     aLearnFuncParam_[i] = an.aLearnFuncParam_[i];
   }
   inputDataMax_ = new double[nInputs_];
   inputDataMin_ = new double[nInputs_];
   for(int i = 0; i < nInputs_; i++)
   {
      inputDataMax_[i] = an.inputDataMax_[i];
      inputDataMin_[i] = an.inputDataMin_[i];
   }
   inputData_ = new double[nInputs_*nSamples_];
   for(int i = 0; i < nInputs_*nSamples_; i++)
   {
      inputData_[i] = an.inputData_[i];
   }
}

//****************************************************************************
// FUNCTION : ~Ann()  
// PURPOSE  : Destructor for object class Ann
//****************************************************************************
Ann::~Ann()
{
   cleanUp();
}

//****************************************************************************
// FUNCTION : checkClassVars()  
// PURPOSE  : Checks class variables for errors.
//****************************************************************************
void Ann::checkClassVars(void) 
{
#ifdef HAVE_SNNS
   if (nNet_ <= 0)
   {
      fprintf(stderr, "ANN::checkClassVars ERROR - invalid num_networks. \n");
      exit(1);
   }	
   if (nHidden_ <= -1)
   {
      fprintf(stderr,"ANN::checkClassVars ERROR - invalid num_hidden_units.\n");
      exit(1);			
   }
   if (normBounds_[0] >= normBounds_[1])
   {
      fprintf(stderr, "ANN::checkClassVars ERROR - normalization bounds");
      fprintf(stderr, " [a,b] must follow a < b. \n ");
      exit(1);
   }
   if (learnFuncName_ == NULL)
   {
      fprintf(stderr,"ANN::checkClassVars ERROR - invalid training_algorithm.\n");
      exit(1);
   }
   if (nMaxRetrain_ < 0)
   {
      fprintf(stderr, "ANN::checkClassVars ERROR - invalid max_retrains. \n");
      exit(1);
   }
   if (nNet_ == 1 && samplingMethod_ == ANN_HOLDOUT)
   {
      fprintf(stderr, "ANN::checkClassVars ERROR - cannot perform holdout");
      fprintf(stderr, " sampling on 1 network. \n");
      exit(1);	
   }
		
   if (strcmp(learnFuncName_, "CC") == 0)
   {
      if (nHidden_ > 0)
      {
         fprintf(stderr, "ANN::checkClassVars ERROR - cannot use CC algorithm");
         fprintf(stderr, " when num_hidden_units > 0. \n");
         exit(1);
      }
      if (nLearnFuncParam_ != 28)
      {
         fprintf(stderr, "ANN::checkClassVars ERROR - cannot specify ");
         fprintf(stderr, "training_num_param if using CC algorithm. \n");
         exit(1);		
      }
   }
   else if (nLearnFuncParam_ != 2 && nLearnFuncParam_ != 4)
   {
      fprintf(stderr, "ANN::checkClassVars ERROR - invalid ");
      fprintf(stderr, "training_num_param. \n ");
      exit(1);
   }	
#endif
}

//****************************************************************************
// FUNCTION : readConfigFile()  
// PURPOSE  : If "Ann.config", "ann.config" or "ANN.config" exists, read and 
//            set ANN variables.
//****************************************************************************
void Ann::readConfigFile(void)
{
#ifdef HAVE_SNNS
   int lineLength = 200;
   char line[300], winput[200];
   char *keywords[] = {"ANN","TOPOLOGY","DATA","TRAINING","OUTPUT","END"};
   bool exists = true;
   FILE *fIn;
	
	
   if      ((fIn=fopen("ANN.config", "r")) != NULL) exists = true;
   else if ((fIn=fopen("ann.config", "r")) != NULL) exists = true;
   else if ((fIn=fopen("Ann.config", "r")) != NULL) exists = true;
   else                                             exists = false;
	
   if (exists == false) return;
	
   while ((fgets(line, lineLength, fIn) != NULL) && (feof(fIn) == 0))
   {
      sscanf(line, "%s", winput);
      if (strcmp(winput, keywords[0]) == 0) { break; }
   }
	
   while ((fgets(line, lineLength, fIn) != NULL) && (feof(fIn) == 0))
   {
      sscanf(line,"%s", winput);
      if (strcmp(winput, keywords[1]) == 0) readTopologySection(fIn);
      else if (strcmp(winput, keywords[2]) == 0) readDataSection(fIn);
      else if (strcmp(winput, keywords[3]) == 0) readTrainingSection(fIn);
      else if (strcmp(winput, keywords[4]) == 0) readOutputSection(fIn);
      else if (strcmp(winput, keywords[5]) == 0) break;
      else if (strcmp(winput,"#") == 0) { /* comment line */ }
      else if (winput[0] == '#') { /* comment line */ }
      else
      {
         fprintf(stderr, "ANN::readConfigFile ERROR - unrecognized line: %s\n",
                 line);
         exit(1);
      }
   }
   fclose(fIn);
   printf("ANN:: Ann.config found. Settings loaded. \n");	
#endif
}

//****************************************************************************
// FUNCTION : readTopologySection(FILE *fp)  
//----------------------------------------------------------------------------
//****************************************************************************
void Ann::readTopologySection(FILE *fp)
{
#ifdef HAVE_SNNS
   int lineLength = 200;
   char line[200], winput[200], winput2[200];
   char *keywords[] = {"num_networks","num_hidden_units","END"};
	
   if (fp == NULL)
   {
      fprintf(stderr,"ANN::readTopologySection ERROR - null file pointer. \n");
      exit(1);
   }
	
   while ((fgets(line, lineLength, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line, "%s", winput);
      if (strcmp(winput, keywords[0]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0) 
            sscanf(line, "%s %s %d", winput, winput2, &nNet_);
         else 
            sscanf(line, "%s %d", winput, &nNet_);	
      }
		
      else if (strcmp(winput, keywords[1]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);
         if (strcmp(winput2, "=") == 0)
            sscanf(line, "%s %s %d", winput, winput2, &nHidden_);
         else 
            sscanf(line, "%s %d", winput, &nHidden_);
      }
		
      else if (strcmp(winput, keywords[2]) == 0) { break; }
      else if (strcmp(winput, "#") == 0) { /* comment line */ }
      else if (winput[0] == '#') { /* comment line */ }
      else
      {
         fprintf(stderr,"ANN::readTopologySection ERROR - ");
         fprintf(stderr, "unrecognized line: %s \n", line);
         exit(1);
      }
   }
   if (feof(fp) != 0)
   {
      fprintf(stderr, "ANN::readTopologySection ERROR - END not found. \n");
      exit(1);
   }
#endif
}

//****************************************************************************
// FUNCTION : readDataSection(FILE *fp)  
//----------------------------------------------------------------------------
//****************************************************************************
void Ann::readDataSection(FILE *fp)
{
#ifdef HAVE_SNNS
   int lineLength = 200;
   char line[200], winput[200], winput2[200], winput3[200];
   char *keywords[]={"sampling_method","nonormalize","normalize_bounds","END"};
	
   if (fp == NULL)
   {
      fprintf(stderr,"ANN::readDataSection ERROR - null file pointer. \n");
      exit(1);
   }
   while ((fgets(line, lineLength, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line, "%s", winput);

      if (strcmp(winput, keywords[0]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);

         if (strcmp(winput2, "=") == 0) 
            sscanf(line, "%s %s %s", winput, winput2, winput3);
         else 
            sscanf(line, "%s %s", winput, winput3);

         if (strcmp(winput3, "Bootstrap") == 0) 
            samplingMethod_ = ANN_BOOTSTRAP;
         else if (strcmp(winput3, "Holdout") == 0)
            samplingMethod_ = ANN_HOLDOUT;
         else
         {
            fprintf(stderr, "ANN::readDataSection ERROR - invalid ");
            fprintf(stderr, "sampling_method. \n");
            exit(1);
         }
      }

      else if (strcmp(winput, keywords[1]) == 0)
      {
         normFlag_ = false;		
      }

      else if (strcmp(winput, keywords[2]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);

         if (strcmp(winput2, "=") == 0)
         {
            sscanf(line, "%s %s %f %f", winput, winput2, &normBounds_[0], 
                   &normBounds_[1]);
         }
         else
         {
            sscanf(line, "%s %f %f", winput, &normBounds_[0], &normBounds_[1]);
         }

         normFlag_ = true;
      }
      else if (strcmp(winput, keywords[3]) == 0) { break; }
      else if (strcmp(winput, "#") == 0) { /* comment line */ }
      else if (winput[0] == '#') { /* comment line */ }
      else 
      {
         fprintf(stderr, "ANN::readDataSection ERROR - unrecognized line: %s\n",
                 line);
         exit(1);
      }
   }
   if (feof(fp) != 0)
   {
      fprintf(stderr, "ANN::readDataSection ERROR - END not found. \n");
      exit(1);
   }
#endif
}

//****************************************************************************
// FUNCTION : readTrainingSection(FILE *fp)  
//****************************************************************************
//****************************************************************************
void Ann::readTrainingSection(FILE *fp)
{
#ifdef HAVE_SNNS
   int lineLength = 200;
   char line[200], winput[200], winput2[200], winput3[200];
   char *keywords[] = {"training_algorithm","training_num_param",
                       "training_param","max_retrains","accuracy_level","END"};
   char *algorithms[] = {"CC","BackpropMomentum","Quickprop","RPROP","SCG",
                         "Std_Backpropagation"};
	
   if (fp == NULL)
   {
      fprintf(stderr,"ANN::readTrainingSection ERROR - null file pointer. \n");
      exit(1);
   }
   while ((fgets(line, lineLength, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line, "%s", winput);
		
      if (strcmp(winput, keywords[0]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);

         if (strcmp(winput2, "=") == 0) 
              sscanf(line, "%s %s %s", winput, winput2, winput3);
         else sscanf(line, "%s %s", winput, winput3);
			
         /* Determine learning function */
         if (strcmp(algorithms[0], winput3) == 0)
            strcpy(learnFuncName_, algorithms[0]);
         else if (strcmp(algorithms[1], winput3) == 0)
            strcpy(learnFuncName_, algorithms[1]);
         else if (strcmp(algorithms[2], winput3) == 0)
            strcpy(learnFuncName_, algorithms[2]);
         else if (strcmp(algorithms[3], winput3) == 0)
            strcpy(learnFuncName_, algorithms[3]);
         else if (strcmp(algorithms[4], winput3) == 0)
            strcpy(learnFuncName_, algorithms[4]);
         else if (strcmp(algorithms[5], winput3) == 0)
            strcpy(learnFuncName_, algorithms[5]);			
      }
		
      else if (strcmp(winput, keywords[1]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);

         if (strcmp(winput2, "=") == 0)
              sscanf(line, "%s %s %d",winput,winput2,&nLearnFuncParam_);
         else sscanf(line, "%s %d", winput, &nLearnFuncParam_);
      }

      else if (strcmp(winput, keywords[2]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);

         if (strcmp(learnFuncName_, "CC") == 0)
         {
            fprintf(stderr, "ANN::readTrainingSection ERROR - cannot \n");
            fprintf(stderr, "specify training_param if using CC algorithm.\n"); 
            exit(1);
         }
			
         if (strcmp(winput2, "=") == 0)
         {				
            if (nLearnFuncParam_ == 2)
               sscanf(line, "%s %s %f %f", winput, winput2, 
                      &aLearnFuncParam_[0], &aLearnFuncParam_[1]);
            else 
               sscanf(line, "%s %s %f %f %f %f", winput, winput2, 
                      &aLearnFuncParam_[0], &aLearnFuncParam_[1], 
                      &aLearnFuncParam_[2], &aLearnFuncParam_[3]);
         }
         else
         {
            if (nLearnFuncParam_ == 2)
               sscanf(line, "%s %f %f", winput, &aLearnFuncParam_[0], 
                      &aLearnFuncParam_[1]);
            else 
               sscanf(line, "%s %f %f %f %f",winput,&aLearnFuncParam_[0],
                      &aLearnFuncParam_[1], &aLearnFuncParam_[2], 
                      &aLearnFuncParam_[3]);			
         }
      }

      else if (strcmp(winput, keywords[3]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);

         if (strcmp(winput2, "=") == 0)		
              sscanf(line, "%s %s %d ", winput, winput2, &nMaxRetrain_);
         else sscanf(line, "%s %d", winput, &nMaxRetrain_);			
      }

      else if (strcmp(winput, keywords[4]) == 0)
      {
         sscanf(line, "%s %s", winput, winput2);

         if (strcmp(winput2, "=") == 0) 
              sscanf(line, "%s %s %d ", winput,winput2,&accuracyLevel_);
         else sscanf(line, "%s %d", winput, &accuracyLevel_);	

         if (accuracyLevel_ <= 0 || accuracyLevel_ > 3) 
         {
            fprintf(stderr, "ANN::checkClassVars ERROR - invalid \n");
            fprintf(stderr, "accuracy_level. Must be 1 (low), 2 (medium) \n");
            fprintf(stderr, "or 3 (high). \n");
            exit(1);
         }
         else if (accuracyLevel_ == 1)
         {
            aLearnFuncParam_[10] = 50;	// Patience (Candidate Units)
            aLearnFuncParam_[11] = 60;	// Max Update Epochs (Candidate Units)
            aLearnFuncParam_[15] = 50;	// Patience (Output Units)
            aLearnFuncParam_[16] = 60;	// Max Update Epochs (Output Units)
         }
         else if (accuracyLevel_ == 3)
         {
            aLearnFuncParam_[9] = 0.08;	// Min Error Change (Candidate Units)
            aLearnFuncParam_[10] = 50;	// Patience (Candidate Units) 
            aLearnFuncParam_[11] = 400;	// Max Update Epochs (Candidate Units) 
            aLearnFuncParam_[15] = 50;	// Patience (Output Units) 
            aLearnFuncParam_[16] = 400;	// Max Update Epochs (Output Units)
            aLearnFuncParam_[12] = 15;	// Max Candidate Units 
         }
         else { /* Default value is 2 */ }
      }

      else if (strcmp(winput, keywords[5]) == 0) { break; }
      else if (strcmp(winput, "#") == 0) { /* comment line */ }
      else if (winput[0] == '#') { /* comment line */ }
      else 
      {
         fprintf(stderr,"ANN::readTrainingSection ERROR - invalid : %s\n",line);
         exit(1);
      }
   }
   if (feof(fp) != 0)
   {
      fprintf(stderr, "ANN::readTrainingSection ERROR - END not found. \n");
      exit(1);
   }
#endif
}

//****************************************************************************
//  FUNCTION : readOutputSection(FILE *fp) 
// ---------------------------------------------------------------------------
//****************************************************************************
void Ann::readOutputSection(FILE *fp)
{
#ifdef HAVE_SNNS
   int lineLength = 200;
   char line[200], winput[200];
   char *keywords[] = {"verbose","save_nets","save_patterns","END"};
	
   if (fp == NULL)
   {
      fprintf(stderr,"ANN::readOutputSection ERROR - null file pointer. \n");
      exit(1);
   }
   while ((fgets(line, lineLength, fp) != NULL) && (feof(fp) == 0))
   {
      sscanf(line, "%s", winput);
		
      /* Check for verbose option */
      if (strcmp(winput, keywords[0]) == 0) { verboseFlag_ = true; }
		
      /* Check for save_nets option */
      else if (strcmp(winput, keywords[1]) == 0) { saveNetFlag_ = true; }
		
      /* Check for save_patterns option */
      else if (strcmp(winput, keywords[2]) == 0) { savePatFlag_ = true; }
		
      else if (strcmp(winput, keywords[3]) == 0) { break; }
      else if (strcmp(winput, "#") == 0) { /* comment line */ }
      else if (winput[0] == '#') { /* comment line */ }
      else
      {
         fprintf(stderr,"ANN::readOutputSection ERROR- unrecognized line: %s\n",
                 line);
         exit(1);
      }
   }
   if (feof(fp) != 0)
   {
      fprintf(stderr, "ANN::readOutputSection ERROR - END not found. \n");
      exit(1);
   }
#endif
}

//****************************************************************************
//  FUNCTION : generateNet(double *X, double *Y)  
//****************************************************************************
void Ann::generateNet(double *X, double *Y) 
{
#ifdef HAVE_SNNS
   float *results = new float[2], *netMSE = new float[nNet_];
   float *netHid = new float[nNet_];
   float sumMSE = 0, sumsquareMSE = 0, meanMSE, stddevMSE; 
   float sumHid = 0, sumsquareHid = 0, meanHid, stddevHid;
   int   minMSE = 0, minHid = 0 , ccFlag, k;

   if (normFlag_ == true)
   {
      inputDataMax_ = new double[nInputs_]; // Maximum input per dimension
      inputDataMin_ = new double[nInputs_]; // Minimum input per dimension
      inputData_ = new double[nInputs_*nSamples_]; // Class copy of input data

      /* Initialize [inputDataMax_] and [inputDataMin_] arrays */
      for (int i = 0; i < nInputs_; i++)
      {
         inputDataMax_[i] = -MAX_FLOAT;
         inputDataMin_[i] = MAX_FLOAT;
      }

      /* Copy and find min/max values from <inputData_> */
      for (int i = 0; i < nSamples_; i++) 
      {
         for (int j = 0; j < nInputs_; j++) 
         {
            inputData_[i*nInputs_+j] = X[i*nInputs_+j];
            if (inputData_[i*nInputs_+j] > inputDataMax_[j])
               inputDataMax_[j] = X[i*nInputs_+j];
            if (inputData_[i*nInputs_+j] < inputDataMin_[j])
               inputDataMin_[j] = X[i*nInputs_+j];	
         }
		
         /* Find max output */
         if (fabs(Y[i]) > maxMagnitudeOutput_)
            maxMagnitudeOutput_ = fabs(Y[i]);
      }

      /* Normalize input */
      for (int i=0; i < nSamples_; i++) 
      {
         for (int j=0; j < nInputs_; j++)
         inputData_[i*nInputs_+j] = normBounds_[0]+
                (normBounds_[1]-normBounds_[0]) *
                (inputData_[i*nInputs_+j]-inputDataMin_[j])/
                (inputDataMax_[j]-inputDataMin_[j]);
      }
   }	
   else { inputData_ = X; }

   scalingFactor_ = 0.0;
   for (k = 0; k < nSamples_; k++) 
      if (fabs(Y[k]) > scalingFactor_) scalingFactor_ = fabs(Y[k]);
   for (k = 0; k < nSamples_; k++) Y[k] /= scalingFactor_;

   ANNF_initialize(nNet_+nMaxRetrain_, nNet_+1);
   ANNF_createNets(nInputs_, nHidden_);

   ANNF_createPat(0, 0, nSamples_-1, nInputs_, nSamples_, inputData_, Y);
	
   if (samplingMethod_ == ANN_BOOTSTRAP) 
   {
      for (int i = 1; i <= nNet_; i++)
         ANNF_createPatRand(i, 0, nSamples_-1, nInputs_, inputData_, Y);	
   }
   else if (samplingMethod_ == ANN_HOLDOUT) 
   {
      ANNF_createPat(0, 0, ((nNet_-1)*nSamples_)/nNet_-1, nInputs_, 
                     nSamples_, inputData_, Y);
      for (int i = 2; i < nNet_; i++) 
         ANNF_createPat(i, (i*nSamples_)/nNet_, ((i-1)*nSamples_)/nNet_-1, 
                        nInputs_, nSamples_, inputData_, Y);
      ANNF_createPat(nNet_, nSamples_/nNet_, nSamples_-1, nInputs_, nSamples_, 
                     inputData_, Y);			
   }
   else if (samplingMethod_ == ANN_ALLOUT) 
   {
      for (int i = 0; i <= nNet_; i++) 
         ANNF_createPat(i, 0, nSamples_-1, nInputs_, nSamples_, inputData_, Y); 
   }
   else
   {
      fprintf(stderr, "ANN::generateNet ERROR - invalid sampling method. \n");
      exit(1);	
   }

   /* Set ANNF class verbose options if necessary */
   if (verboseFlag_ == true) { ANNF_setVerboseFlag(); }
	
   /* Set ccFlag if CasCor algorithm is used */
   if (strcmp(learnFuncName_, "CC") == 0) { ccFlag = ANNF_CC; }
   else { ccFlag = ANNF_NOCC; }
	
   printf("\n");
   printf("* ============================================== * \n");
	
   if (verboseFlag_ == true)
      printf("* ANN training begins...\n");
   for (int i = 0; i < nNet_; i++)
   {
      ANNF_trainNet(i, i+1, learnFuncName_, aLearnFuncParam_, 
                    nLearnFuncParam_, &results, ccFlag);
      if (verboseFlag_ == true)
         printf("* Network %d: finished training with %s algorithm. \n", 
                i+1, learnFuncName_);

      /* Update MSEs */
      netMSE[i] = results[0];
      sumMSE += netMSE[i];
      sumsquareMSE += netMSE[i]*netMSE[i];
      if (netMSE[i] < netMSE[minMSE]) minMSE = i;

      /* Update Hiddens */
      netHid[i] = results[1] - nInputs_ - 1;
      sumHid += netHid[i];
      sumsquareHid += netHid[i]*netHid[i];
      if (netHid[i] < netHid[minHid]) minHid = i;		
   }

   meanMSE = sumMSE / nNet_;
   stddevMSE = sqrt((sumsquareMSE - nNet_*meanMSE*meanMSE)/(nNet_));

   meanHid = sumHid / nNet_;
   stddevHid = sqrt((sumsquareHid - nNet_*meanHid*meanHid)/(nNet_));

   if (nMaxRetrain_ > 0)
   {
      float *netMSERetrain = new float[nMaxRetrain_];
      float *netHidRetrain = new float[nMaxRetrain_];
      float bestMSE; 
      int bestNet;
		
      for (int i = 0; i < nNet_; i++)
      {
         bestNet = i;
         bestMSE = netMSE[i];			
			
         for (int j = 0; netHid[i] - meanHid > 1.25*stddevHid && 
              j < nMaxRetrain_; j++)
         {
            ANNF_trainNet(nNet_+j, i+1, learnFuncName_, aLearnFuncParam_, 
                          nLearnFuncParam_, &results, ccFlag);				
            if (verboseFlag_ == true)
               printf("* Network %d: finished RE-training (attempt %d of %d).\n",
                      i+1, j+1, nMaxRetrain_);				
				
            /* Store new network data */
            netMSERetrain[j] = results[0];
            netHidRetrain[j] = results[1] - nInputs_ - 1;
				
            if (bestMSE > results[0])
            {
               bestMSE = results[0];
               bestNet = nNet_+j;
            }
				
            /* Break if criteria satisfied */ 
            if (netHidRetrain[j] - meanHid <= 1.25*stddevHid) break;
         }

         ANNF_swapNet(i, bestNet);

         /* Update data if current network has been replaced */
         if (bestNet != i)
         {
            /* Remove old data */
            sumsquareMSE -= netMSE[i]*netMSE[i];
            sumMSE -= netMSE[i];
            sumsquareHid -= netHid[i]*netHid[i];
            sumHid -= netHid[i];
	 				
            /* Input new data */
            netMSE[i] = netMSERetrain[bestNet-nNet_];
            netHid[i] = netHidRetrain[bestNet-nNet_];
            sumMSE += netMSERetrain[bestNet-nNet_];
            sumsquareMSE += netMSERetrain[bestNet-nNet_]*
                            netMSERetrain[bestNet-nNet_];
            sumHid += netHidRetrain[bestNet-nNet_];
            sumsquareHid += (netHidRetrain[bestNet-nNet_])*
                            (netHidRetrain[bestNet-nNet_]);
         }
      }
      delete [] netMSERetrain;
      delete [] netHidRetrain;
   }
	
   meanMSE = sumMSE / nNet_;
   stddevMSE = sqrt((sumsquareMSE - nNet_*meanMSE*meanMSE)/(nNet_));

   meanHid = sumHid / nNet_;
   stddevHid = sqrt((sumsquareHid - nNet_*meanHid*meanHid)/(nNet_));

   if (verboseFlag_ == true)
   {	
      meanHid = sumHid / nNet_;
      stddevHid = sqrt((sumsquareHid - nNet_*meanHid*meanHid)/(nNet_));

      printf("\n");
      printf("* MSE:    MEAN: %e  MINIMUM: %e  STDDEV: %e \n", meanMSE, 
             netMSE[minMSE], stddevMSE);
      printf("* HIDDEN: MEAN: %e  MINIMUM: %e  STDDEV: %e \n", meanHid, 
             netHid[minHid], stddevHid);
   }

   printf("* ============================================== * \n");

   if (saveNetFlag_ == true) 
   {
      char buffer[30];
      for (int i = 0; i < nNet_; i++)
      {
         sprintf(buffer, "annData.net%d.net", i+1);
         ANNF_saveNet(i, buffer);
      }
   }
   if (savePatFlag_ == true) 
   {
      ANNF_savePat(0, "annData.test.pat");
      char buffer[30];
      for (int i = 1; i <= nNet_; i++) 
      {
         sprintf(buffer, "annData.train%d.pat", i);
         ANNF_savePat(i, buffer);
      }
   }

   trainPerformed_ = true;
   delete [] results;
   delete [] netMSE;
   delete [] netHid;
#endif
}

//****************************************************************************
//  FUNCTION : initialize(double *X, double *Y)
//  PURPOSE  : Creates network 
//****************************************************************************
int Ann::initialize(double *X, double *Y)
{
#ifdef HAVE_SNNS
   if (nInputs_ <= 0 || nSamples_ <= 0)
   {
      printf("ANN::initialize ERROR - invalid argument.\n");
      exit(1);
   } 
   if (nSamples_ <= nInputs_)
   {
      printf("ANN::initialize INFO - not enough points.\n");
      return 0;
   }
   
   cleanUp();
   if (outputLevel_ > 3) verboseFlag_ = true;
   if (trainPerformed_ == false) { generateNet(X, Y); }
   if (psRSCodeGen_ == 1) 
      printf("ANN INFO: response surface stand-alone code not available.\n");
#endif
   return 0;
}

//****************************************************************************
//  FUNCTION : genNDGridData(double *X, double *Y, int *N, 
//                           double **X2, double **Y2)  
//  PURPOSE  : Creates network and generates 1-D lattice points if necessary.
//****************************************************************************
int Ann::genNDGridData(double *X, double *Y, int *N, double **X2, double **Y2)
{
#ifdef HAVE_SNNS
   int totPts, mm;

   initialize(X, Y);

   if ((*N) == -999) return 0;

   genNDGrid(N, X2);
   if ((*N) == 0) return 0;
   totPts = (*N);

   (*Y2) = new double[totPts];
   (*N)  = totPts;
   for (mm = 0; mm < totPts; mm++)
      (*Y2)[mm] = evaluatePoint(&((*X2)[mm*nInputs_]));
#endif
   return 0;
}

//****************************************************************************
//  FUNCTION : gen1DGridData
//****************************************************************************
int Ann::gen1DGridData(double *X, double *Y, int ind1, double *settings, 
                       int *NN, double **XX, double **YY)
{
#ifdef HAVE_SNNS
   int    mm, nn;
   double HX, *Xloc, std;

   initialize(X, Y);

   if ((*N) == -999) return 0;

   HX = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 

   (*NN) = nPtsPerDim_;
   (*XX) = new double[nPtsPerDim_];
   (*YY) = new double[nPtsPerDim_];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 

   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      Xloc[ind1] = HX * mm + lowerBounds_[ind1];
      (*XX)[mm]   = Xloc[ind1];
      (*YY)[mm] = evaluatePointFuzzy(Xloc,std);
      if (outputLevel_ > 3)
      {
         printf("  X %4d = %12.4e\n", ind1, Xloc[ind1]);
         printf("  Y     = %12.4e (std = %12.4e)\n", (*YY)[mm],std);
      }
   }

   delete [] Xloc;
#endif
   return 0;
}

//****************************************************************************
//  FUNCTION : gen2DGridData
//****************************************************************************
int Ann::gen2DGridData(double *X, double *Y, int ind1, int ind2, 
                       double *settings, int *NN, double **XX, double **YY)
{
#ifdef HAVE_SNNS
   int    totPts, mm, nn, index;
   double *HX, *Xloc, std;

   initialize(X, Y);

   if ((*N) == -999) return 0;

   totPts = nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[2];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 2];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
	    
   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         index = mm * nPtsPerDim_ + nn;
         Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
         Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
         (*XX)[index*2]   = Xloc[ind1];
         (*XX)[index*2+1] = Xloc[ind2];
         (*YY)[index] = evaluatePointFuzzy(Xloc,std);
         if (outputLevel_ > 3)
         {
            printf("  X %4d = %12.4e\n", ind1, Xloc[ind1]);
            printf("  X %4d = %12.4e\n", ind2, Xloc[ind2]);
            printf("  Y     = %12.4e (std = %12.4e)\n", (*YY)[index],std);
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
#endif
   return 0;
}

//****************************************************************************
//  FUNCTION : gen3DGridData
//****************************************************************************
int Ann::gen3DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                      double *settings, int *NN, double **XX, double **YY)
{
#ifdef HAVE_SNNS
   int    totPts, mm, nn, pp, index;
   double *HX, *Xloc, std;

   initialize(X, Y);

   if ((*N) == -999) return 0;

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[3];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 3];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
	    
   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            index = mm * nPtsPerDim_ * nPtsPerDim_ + nn * nPtsPerDim_ + pp;
            Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
            Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
            Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
            (*XX)[index*3]   = Xloc[ind1];
            (*XX)[index*3+1] = Xloc[ind2];
            (*XX)[index*3+2] = Xloc[ind3];
            (*YY)[index] = evaluatePointFuzzy(Xloc,std);
            if (outputLevel_ > 3)
            {
               printf("  X %4d = %12.4e\n", ind1, Xloc[ind1]);
               printf("  X %4d = %12.4e\n", ind2, Xloc[ind2]);
               printf("  X %4d = %12.4e\n", ind3, Xloc[ind3]);
               printf("  Y     = %12.4e (std = %12.4e)\n", (*YY)[index],std);
            }
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
#endif
   return 0;
}

//****************************************************************************
//  FUNCTION : gen4DGridData
//****************************************************************************
int Ann::gen4DGridData(double *X, double *Y, int ind1, int ind2, int ind3,
                       int ind4, double *settings, int *NN, double **XX, 
                       double **YY)
{
#ifdef HAVE_SNNS
   int    totPts, mm, nn, pp, qq, index;
   double *HX, *Xloc, std;

   initialize(X, Y);

   if ((*N) == -999) return 0;

   totPts = nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_ * nPtsPerDim_;
   HX    = new double[4];
   HX[0] = (upperBounds_[ind1] - lowerBounds_[ind1]) / (nPtsPerDim_ - 1); 
   HX[1] = (upperBounds_[ind2] - lowerBounds_[ind2]) / (nPtsPerDim_ - 1); 
   HX[2] = (upperBounds_[ind3] - lowerBounds_[ind3]) / (nPtsPerDim_ - 1); 
   HX[3] = (upperBounds_[ind4] - lowerBounds_[ind4]) / (nPtsPerDim_ - 1); 

   (*NN) = totPts;
   (*XX) = new double[totPts * 4];
   (*YY) = new double[totPts];
   Xloc  = new double[nInputs_];
   for (nn = 0; nn < nInputs_; nn++) Xloc[nn] = settings[nn]; 
	    
   for (mm = 0; mm < nPtsPerDim_; mm++)
   {
      for (nn = 0; nn < nPtsPerDim_; nn++)
      {
         for (pp = 0; pp < nPtsPerDim_; pp++)
         {
            for (qq = 0; qq < nPtsPerDim_; qq++)
            {
               index = mm*nPtsPerDim_*nPtsPerDim_*nPtsPerDim_ + 
                       nn*nPtsPerDim_*nPtsPerDim_ + pp * nPtsPerDim_ + qq;
               Xloc[ind1] = HX[0] * mm + lowerBounds_[ind1];
               Xloc[ind2] = HX[1] * nn + lowerBounds_[ind2];
               Xloc[ind3] = HX[2] * pp + lowerBounds_[ind3];
               Xloc[ind4] = HX[3] * qq + lowerBounds_[ind4];
               (*XX)[index*4]   = Xloc[ind1];
               (*XX)[index*4+1] = Xloc[ind2];
               (*XX)[index*4+2] = Xloc[ind3];
               (*XX)[index*4+3] = Xloc[ind4];
               (*YY)[index] = evaluatePointFuzzy(Xloc,std);
               if (outputLevel_ > 3)
               {
                  printf("  X %4d = %12.4e\n", ind1, Xloc[ind1]);
                  printf("  X %4d = %12.4e\n", ind2, Xloc[ind2]);
                  printf("  X %4d = %12.4e\n", ind3, Xloc[ind3]);
                  printf("  X %4d = %12.4e\n", ind4, Xloc[ind4]);
                  printf("  Y     = %12.4e (std = %12.4e)\n", (*YY)[index],std);
               }
            }
         }
      }
   }

   delete [] Xloc;
   delete [] HX;
#endif
   return 0;
}

//****************************************************************************
// FUNCTION : evaluatePoint  
//****************************************************************************
double Ann::evaluatePoint(double *x)
{
#ifdef HAVE_SNNS
   float output = 0;
   float *xPoint = new float[nInputs_];

   if (normFlag_ == true) 
   {
      for (int i = 0; i < nInputs_; i++)
         xPoint[i] = (float)(normBounds_[0] + 
                     (normBounds_[1]-normBounds_[0])*
                     (x[i] - inputDataMin_[i])/
                     (inputDataMax_[i] - inputDataMin_[i]));
   }	
   /* Else, copy point */	
   else
   {
      for (int i = 0; i < nInputs_; i++) xPoint[i] = (float)x[i];
   }

   /* Evaluate point on all networks */
   for (int i = 0; i < nNet_; i++)
      output += ANNF_evaluatePoint(i, xPoint);
   output *= scalingFactor_;

   delete [] xPoint;
   return output/nNet_;
#else
   return -1;
#endif
}

//****************************************************************************
// FUNCTION : evaluatePoint  
//****************************************************************************
double Ann::evaluatePoint(int npts, double *x, double *y)
{
#ifdef HAVE_SNNS
   int   ii, kk;
   float output = 0;
   float *xPoint = new float[nInputs_];

   for (kk = 0; kk < npts; kk++)
   {
      if (normFlag_ == true) 
      {
         for (ii = 0; ii < nInputs_; ii++)
            xPoint[ii] = (float)(normBounds_[0] + 
                         (normBounds_[1]-normBounds_[0])*
                         (x[kk*nInputs_+ii] - inputDataMin_[ii])/
                         (inputDataMax_[ii] - inputDataMin_[ii]));
      }
      /* Else, copy point */
      else
      {
         for (ii = 0; ii < nInputs_; ii++) xPoint[ii] = (float)x[kk*nInputs_+ii];
      }

      /* Evaluate point on all networks */
      output = 0.0;
      for (ii = 0; ii < nNet_; ii++) output += ANNF_evaluatePoint(ii, xPoint);
      y[kk] = output * scalingFactor_;
      y[kk] /= nNet_;
   }
   delete [] xPoint;
   return 0.0;
#else
   return -1;
#endif
}

//****************************************************************************
// FUNCTION : evaluatePointFuzzy  
//****************************************************************************
double Ann::evaluatePointFuzzy(double *X, double &std)
{
#ifdef HAVE_SNNS
   int    ii;
   float  *xPoint = new float[nInputs_];
   double *netValues, mean;

   if (normFlag_ == true) 
   {
      for (ii = 0; ii < nInputs_; ii++)
         xPoint[ii] = (float)(normBounds_[0] + 
                     (normBounds_[1]-normBounds_[0])*
                     (X[ii] - inputDataMin_[ii])/
                     (inputDataMax_[ii] - inputDataMin_[ii]));
   }	
   /* Else, copy point */	
   else
   {
      for (ii = 0; ii < nInputs_; ii++) xPoint[ii] = (float) X[ii];
   }

   /* Evaluate point on all networks */
   netValues = new double[nNet_];
   mean = 0.0;
   if (outputLevel_ > 2) printf("ANN evaluatePointFuzzy \n");
   for (ii = 0; ii < nInputs_; ii++)
      printf("**   Input %4d = %12.4e\n", ii+1, X[ii]);
   for (ii = 1; ii < nNet_; ii++) 
   {
      netValues[ii] = (double) ANNF_evaluatePoint(ii, xPoint);
      netValues[ii] *= scalingFactor_;
      if (outputLevel_ > 2)
         printf("ANN interpolant %3d = %12.4e\n",ii,netValues[ii]);
      mean += netValues[ii];
   }
   if (nNet_ > 1) mean /= (double) (nNet_ - 1);
   std = 0.0;
   for (ii = 1; ii < nNet_; ii++) 
      std += (netValues[ii] - mean) * (netValues[ii] - mean);
   if (nNet_ > 2) std = sqrt(std/(double) (nNet_-2));
   if (outputLevel_ > 2) printf("ANN mean/std = %12.4e %12.4e\n",mean,std);
   delete [] xPoint;
   delete [] netValues;
   return mean;
#else
   return -1;
#endif
}

//****************************************************************************
// FUNCTION : evaluatePointFuzzy  
//****************************************************************************
double Ann::evaluatePointFuzzy(int npts, double *X, double *Y, double *Ystd)
{
#ifdef HAVE_SNNS
   int    ii, kk;
   float  *xPoint = new float[nInputs_];
   double *netValues, mean;

   netValues = new double[nNet_];
   for (kk = 0; kk < npts; kk++)
   {
      if (normFlag_ == true) 
      {
         for (ii = 0; ii < nInputs_; ii++)
            xPoint[ii] = (float)(normBounds_[0] + 
                        (normBounds_[1]-normBounds_[0])*
                        (X[ii] - inputDataMin_[ii])/
                        (inputDataMax_[ii] - inputDataMin_[ii]));
      }
      /* Else, copy point */	
      else
      {
         for (ii = 0; ii < nInputs_; ii++) xPoint[ii] = (float) X[ii];
      }

      /* Evaluate point on all networks */
      if (outputLevel_ > 2)
      {
         printf("ANN evaluatePointFuzzy \n");
         for (ii = 0; ii < nInputs_; ii++)
            printf("**   Input %4d = %12.4e\n", ii+1, X[ii]);
      }
      mean = 0.0;
      for (ii = 1; ii < nNet_; ii++) 
      {
         netValues[ii] = (double) ANNF_evaluatePoint(ii, xPoint);
         netValues[ii] *= scalingFactor_;
         if (outputLevel_ > 2)
            printf("ANN interpolant %3d = %12.4e\n",ii,netValues[ii]);
         mean += netValues[ii];
      }
      if (nNet_ > 1) mean /= (double) (nNet_ - 1);
      Y[kk] = mean;
      double thisStd = 0.0;
      for (ii = 1; ii < nNet_; ii++) 
         thisStd += (netValues[ii] - mean) * (netValues[ii] - mean);
      if (nNet_ > 2) thisStd = sqrt(thisStd/(double) (nNet_-2));
      Ystd[kk] = thisStd;
      if (outputLevel_ > 2) printf("ANN mean/std = %12.4e %12.4e\n",mean,thisStd);
   }
   delete [] xPoint;
   delete [] netValues;
   return 0.0;
#else
   return -1;
#endif
}

//****************************************************************************
// FUNCTION : cleanUp  
//****************************************************************************
void Ann::cleanUp()
{
   if (trainPerformed_ == true)
   {
      if (learnFuncName_ != NULL) free(learnFuncName_);
      if (inputDataMax_ != NULL) delete [] inputDataMax_;
      if (inputDataMin_ != NULL) delete [] inputDataMin_;
      if (inputData_       != NULL) delete [] inputData_;
#ifdef HAVE_SNNS
      ANNF_exit();
#endif
      trainPerformed_ = false;
   }
}

