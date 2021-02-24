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
// DATE   : 2006
// ************************************************************************

// ------------------------------------------------------------------------
// system includes
// ------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "PsuadeBase.h"
#include "PsuadeUtil.h"
#include "PDFBase.h"
#include "PsuadeData.h"
#include "PrintingTS.h"

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
int PsuadeBase::setupGuide()
{
   printOutTS(PL_INFO, 
      "Suppose you have an application and you would like to use \n");
   printOutTS(PL_INFO, 
      "PSUADE to perform UQ on it. In the following you will find a\n");
   printOutTS(PL_INFO, 
      "short guide on how to do it:\n");
   printOutTS(PL_INFO, 
      "(0) Suppose the name of your application is called FOO and\n");
   printOutTS(PL_INFO, 
      "    it takes its input from a file called FOO.in using the\n");
   printOutTS(PL_INFO,"    following command:\n");
   printOutTS(PL_INFO,"       srun -n 4 /home/me/FOO FOO.in\n");
   printOutTS(PL_INFO, 
      "    and the parameters you will vary live in FOO.in.\n");
   printOutTS(PL_INFO, 
      "(1) First create a driver (say, FOOdriver.py) using the\n");
   printOutTS(PL_INFO, 
      "    `gendriver' command. This driver will be called PSUADE\n");
   printOutTS(PL_INFO, 
      "    to set up and to run your application code. PSUADE will,\n");
   printOutTS(PL_INFO, 
      "    for each sample point, outputs a parameter file called\n");
   printOutTS(PL_INFO, 
      "    'psuadeApps_ct.in.x', then calls FOOdriver.py with\n");
   printOutTS(PL_INFO, 
      "        FOOdriver.py psuadeApps_ct.in.x psuadeApps_ct.out.x\n");
   printOutTS(PL_INFO, 
      "    and expects the outputs to be written to psuadeApps_ct.out.x\n");
   printOutTS(PL_INFO,"    So FOOdriver.py is expected to do 4 things:\n");
   printOutTS(PL_INFO,"    (a) take the inputs from 'psuadeApps_ct.in.x'\n");
   printOutTS(PL_INFO,"    (b) insert the input values into FOO.in\n");
   printOutTS(PL_INFO, 
      "    (c) run the application code and extract the desired outputs\n");
   printOutTS(PL_INFO,"    (d) write the output to 'psuadeApps_ct.out.x'\n");
   printOutTS(PL_INFO, 
      "    Again, the python template for this step can be accessed via\n");
   printOutTS(PL_INFO,"    the 'gendriver' command.\n");
   printOutTS(PL_INFO, 
      "(2) Then use 'geninputfile' to create a PSUADE input file (say \n");
   printOutTS(PL_INFO, 
      "    psuade.in).  This file specifies the number and names of \n");
   printOutTS(PL_INFO, 
      "    inputs/outputs, sets up the sampling and analysis methods, \n");
   printOutTS(PL_INFO, 
      "    and set up a link to the driver created in step (1).\n");
   printOutTS(PL_INFO, 
      "(3) Optionally, use 'genbatchfile' to create a batch file.\n");
   printOutTS(PL_INFO,"    (If you run the jobs on LLNL's machines).\n");
   printOutTS(PL_INFO, 
      "(4) Once these files have been created, run PSUADE on your\n");
   printOutTS(PL_INFO,"    application just by typing:\n");
   printOutTS(PL_INFO,"         psuade psuade.in\n");
   printOutTS(PL_INFO, 
      "    This is done when the run time is relatively short, no more\n");
   printOutTS(PL_INFO, 
      "    than a few minutes. If the run time is long, you may want to\n");
   printOutTS(PL_INFO,"    break up this step into a few tasks: \n");
   printOutTS(PL_INFO,"    (a) create all the parameter files first\n");
   printOutTS(PL_INFO, 
      "        (use 'gen_inputfile_only' in the psuade.in file)\n");
   printOutTS(PL_INFO,"    (b) launch the jobs in whatever way you desire\n");
   printOutTS(PL_INFO,
      "    (c) run postprocessing to generate all output files\n");
   printOutTS(PL_INFO,"    (d) run PSUADE to collect all output files\n");
   return 0;
}

// ************************************************************************
// create LLNL-specific batch file
// ------------------------------------------------------------------------
int PsuadeBase::genBatchFile(int genFlag)
{
  char batchName[200], dirName[200], pString[501];
  FILE *fp;

  printOutTS(PL_INFO, 
     "INFO: your specified batch file will be appended with .Tmplt.\n");
  sprintf(pString, "Enter the name of the batch file: ");
  getString(pString, batchName);
  strcpy(&batchName[strlen(batchName)-1], ".Tmplt\0");
  fp = fopen(batchName, "w");
  if (fp == NULL)
  {
    printOutTS(PL_ERROR, "ERROR: Cannot open file %s.\n", batchName);
    return 1;
  }
  sprintf(pString,
      "Enter the absolute path of the run directory (no / at end): ");
  getString(pString, dirName);

  fprintf(fp,"###Modify the following lines as deemed necesaary\n");
  fprintf(fp,"#MSUB -s /bin/csh\n");
  fprintf(fp,"#MSUB -b <bankName>\n");
  fprintf(fp,"#MSUB -c <machine>\n");
  fprintf(fp,"#MSUB -tM 1:59 (<Time>)\n");
  fprintf(fp,"#MSUB -ln <numNodes>\n");
  fprintf(fp,"#MSUB -x\n\n");

  fprintf(fp,"###The following lines are used for daisy-chain mode\n");
  fprintf(fp,"###available on LLNL Livermore Computing (LC) systems.\n");
  fprintf(fp,"###Replace the keyword PSUB_JOBID with the correct one.\n");
  fprintf(fp,"###This assumes that your work directory is workdir.<num>\n");
  fprintf(fp,"###cd %s/workdir.PSUADE_NEXT\n", dirName);
  fprintf(fp,"###/usr/bin/msub -d PSUB_JOBID %s.PSUADE_NEXT\n\n",batchName);

  fprintf(fp, "date\n\n");
  fprintf(fp,"###The following line is running your job \n");
  fprintf(fp,"cd %s/workdir.PSUADE_COUNTER\n\n", dirName);
  fprintf(fp,"###Replace the following line with the calling sequence\n");
  fprintf(fp,"###srun -n 1 <executable code> <input decks> <options>\n");
  fprintf(fp, "\ndate\n");
  fclose(fp);
  return 0;
}

// ************************************************************************
// generate psuade driver python or C script
// ------------------------------------------------------------------------
int PsuadeBase::genDriver(int genFlag)
{
  int  ii, nInputs, nOutputs, nFiles, driverID;
  int  haveBatch=0, nSupportFiles;
  char **appInFiles, **targetFiles, response[201], pythonDir[200];
  char driverName[200], winput[501], pString[501];
  FILE *fp, *dfp;

  sprintf(pString,"(1) C or (2) Python driver (1 or 2) ? ");
  driverID = getInt(1, 2, pString);
  if (genFlag == 1)
  {
    if (driverID == 1) strcpy(driverName, "simulator.c");
    else               strcpy(driverName, "simulator.py");
  }
  else
  {
    if (driverID == 1)
       sprintf(pString,
            "Enter the name of the driver (with .c at the end) : ");
    else
       sprintf(pString,
            "Enter the name of the driver (with .py at the end) : ");
    getString(pString, driverName);
    driverName[strlen(driverName)-1] = '\0';
  }

  dfp = fopen(driverName, "r");
  if (dfp != NULL)
  {
    fclose(dfp);
    if (genFlag == 1)
    {
      dfp = fopen(driverName, "w");
      if (dfp == NULL)
      {
        printOutTS(PL_ERROR,"ERRPR: Cannot open file %s. \n", driverName);
        return 1;
      }
    }
    else
    {
      sprintf(pString, "File exists. Override ? (y or n) ");
      getString(pString, response);
      if (response[0] != 'y') return 1;
      else
      {
        dfp = fopen(driverName, "w");
        if (dfp == NULL) 
        {
          printOutTS(PL_ERROR,"ERROR: Cannot open file %s. \n",driverName);
          return 1;
        }
      }
    }
  }
  else
  {
    dfp = fopen(driverName, "w");
    if (dfp == NULL)
    {
      printOutTS(PL_ERROR,"ERROR: Cannot open file %s.\n",driverName);
      return 1;
    }
  }

  if (driverID == 1)
  {
    fprintf(dfp,"/***************************************************\n");
    fprintf(dfp," * Main program\n");
    fprintf(dfp," *=================================================*/\n\n");
    fprintf(dfp,"#include <math.h>\n");
    fprintf(dfp,"#include <stdio.h>\n");
    fprintf(dfp,"#include <stdlib.h>\n\n");
    fprintf(dfp,"main(int argc, char **argv)\n");
    fprintf(dfp,"{\n");
    fprintf(dfp,"   int    ii, nInputs;\n");
    fprintf(dfp,"   double *X, Y=1.0e35;\n");
    fprintf(dfp,"   FILE   *infile, *outfile;\n\n");
    fprintf(dfp,"   infile  = fopen(argv[1], \"r\");\n");
    fprintf(dfp,"   if (infile == NULL)\n");
    fprintf(dfp,"   {\n");
    fprintf(dfp,"      printf(\"ERROR - cannot open input file.\\n\");\n");
    fprintf(dfp,"      exit(1);\n");
    fprintf(dfp,"   }\n");
    fprintf(dfp,"   fscanf(infile, \"%%d\", &nInputs);\n");
    fprintf(dfp,"   if (nInputs <= 0)\n");
    fprintf(dfp,"   {\n");
    fprintf(dfp,"      printf(\"ERROR - nInputs <= 0.\\n\");\n");
    fprintf(dfp,"      exit(1);\n");
    fprintf(dfp,"   }\n");
    fprintf(dfp,"   X = (double *) malloc(nInputs*sizeof(double));\n");
    fprintf(dfp,"   for (ii = 0; ii < nInputs; ii++)\n");
    fprintf(dfp,"      fscanf(infile, \"%%lg\", &X[ii]);\n");
    fprintf(dfp,"   fclose(infile);\n\n");
    if (genFlag == 1)
      fprintf(dfp,"   Y = X[0] - 19.6 * X[1] / (1.5 * X[2]);\n");
    else
      fprintf(dfp,"   /* add codes for computing the function */ \n\n");
    fprintf(dfp,"   outfile  = fopen(argv[2], \"w\");\n");
    fprintf(dfp,"   if (outfile == NULL)\n");
    fprintf(dfp,"   {\n");
    fprintf(dfp,"      printf(\"ERROR-cannot write to outfile.\\n\");\n");
    fprintf(dfp,"      exit(1);\n");
    fprintf(dfp,"   }\n");
    fprintf(dfp,"   fprintf(outfile, \" %%24.16e\\n\", Y);\n");
    fprintf(dfp,"   fclose(outfile);\n");
    fprintf(dfp,"}\n");
    fclose(dfp);
    printOutTS(PL_INFO, "PSUADE: the user driver is in %s\n", driverName);
    return 0;
  }
      
  printEquals(PL_INFO, 0);
  printOutTS(PL_INFO,"Begin creating a Python-based application driver.\n");
  strcpy(pythonDir, "/usr/local/bin/python");
  fp = fopen(pythonDir, "r");
  if (fp != NULL) fclose(fp);
  else
  {
    strcpy(pythonDir, "/usr/bin/python");
    fp = fopen(pythonDir, "r");
    if (fp != NULL) fclose(fp);
    else
    {
      printOutTS(PL_INFO, 
         "Python not found in /usr/bin nor /usr/local/bin.\n");
      printOutTS(PL_INFO, "You will have to edit the driver later.\n");
     }
  }
  fprintf(dfp, "#!%s\n", pythonDir);
  fprintf(dfp, "import os\n"); 
  fprintf(dfp, "import sys\n"); 
  fprintf(dfp, "import string\n"); 
  fprintf(dfp, "import shutil\n"); 

  fprintf(dfp,"#######################################################\n");
  fprintf(dfp,"# BEGIN USER SPECIFIC SECTION\n");
  fprintf(dfp,"#======================================================\n\n");
  if (genFlag == 1) nInputs = 3;
  else
  {
    sprintf(pString,"How many uncertain variables? ");
    nInputs = getInt(1, 10000000, pString);
  }
  fprintf(dfp,"nInputs = %d\n", nInputs);
  fprintf(dfp,"inputNames = [");
  if (genFlag == 1) fprintf(dfp,"\"H0\", \"M\", \"S\"");
  else
  {
    for (ii = 0; ii < nInputs; ii++)
    {
      sprintf(pString, "Enter the name of input variable %d : ",ii+1);
      getString(pString, response);
      response[strlen(response)-1] = '\0';
      fprintf(dfp,"\"%s\"", response);
      if (ii < nInputs-1) fprintf(dfp,", ");
    }
  }
  fprintf(dfp,"]\n\n");
  nFiles = 0;
  if (genFlag != 1)
  {
    printOutTS(PL_INFO, 
       "Your uncertain variables may scatter over several files\n");
    printOutTS(PL_INFO, 
       "needed by your application code. You need to prepare for\n");
    printOutTS(PL_INFO, 
       "automatic insertion of sample values to this file by first\n");
    printOutTS(PL_INFO, 
       "changing these file names via adding .Tmplt to their names.\n");
    printOutTS(PL_INFO, 
       "Then, you need to open up these .Tmplt files, find where the\n");
    printOutTS(PL_INFO, 
       "uncertain variables are, and replace their default values\n");
    printOutTS(PL_INFO, 
       "with a unique symbol (e.g. psThresh for some threshold).\n");
    printOutTS(PL_INFO, 
       "During preprocessing, these symbols will be replaced by the\n");
    printOutTS(PL_INFO, 
       "application driver with the actual sample values.\n");
    sprintf(pString,
       "How many application files contains uncertain input variables? ");
    nFiles = getInt(0, 1000, pString);
    if (nFiles > 0) appInFiles  = new char*[nFiles];
    if (nFiles > 0) targetFiles = new char*[nFiles];
  }
  fprintf(dfp,"# ====> files to set the input variables\n");
  if (nFiles > 0)
  {
    sprintf(pString,
       "Where can your application template be found(absolute path):");
    getString(pString, response);
    response[strlen(response)-1] = '/';
    response[strlen(response)] = '\0';
    fprintf(dfp,"appDir = \"%s\"\n", response);
  }
  for (ii = 0; ii < nFiles; ii++)
  {
    appInFiles[ii] = new char[200];
    sprintf(pString,
        "Enter the name of application input template file %d : ",ii+1);
    getString(pString, appInFiles[ii]);
    targetFiles[ii] = new char[200];
    sprintf(pString,
            "Enter the name of application input file %d : ", ii+1);
    getString(pString, targetFiles[ii]);
  }
  fprintf(dfp,"appInputTmplts = ["); 
  for (ii = 0; ii < nFiles; ii++)
  {
    appInFiles[ii][strlen(appInFiles[ii])-1] = '\0';
    fprintf(dfp,"\"%s\"", appInFiles[ii]); 
    if (ii < nFiles-1) fprintf(dfp,", ");
    delete [] appInFiles[ii];
  }
  fprintf(dfp,"]\n");
  if (nFiles > 0) delete [] appInFiles;

  fprintf(dfp,"appInputFiles = ["); 
  for (ii = 0; ii < nFiles; ii++)
  {
    targetFiles[ii][strlen(targetFiles[ii])-1] = '\0';
    fprintf(dfp,"\"%s\"", targetFiles[ii]); 
    if (ii < nFiles-1) fprintf(dfp,", ");
    delete [] targetFiles[ii];
  }
  if (nFiles > 0) delete [] targetFiles;
  fprintf(dfp,"]\n\n");

  if (genFlag == 1) strcpy(winput, "n");
  else
  {
    printOutTS(PL_INFO, 
       "You might have a template batch file for your computer.\n");
    printOutTS(PL_INFO, 
       "Or, you might have use 'genbatchfile' to generate one.\n");
    sprintf(pString, "Do you have a batch file ? (y or n) ");
    getString(pString, winput);
  }
  if (winput[0] == 'y')
  {
    fprintf(dfp,"# ====> files to set up the batch file \n");
    printOutTS(PL_INFO, 
       "INFO: batch template files are to be modified with run numbers.\n");
    sprintf(pString,"Where can your batch template be found(absolute path):");
    getString(pString, response);
    response[strlen(response)-1] = '/';
    response[strlen(response)] = '\0';
    fprintf(dfp,"batchDir = \"%s\"", response);
    sprintf(pString,"What is the name of your batch template file? ");
    getString(pString, response);
    response[strlen(response)-1] = '\0';
    fprintf(dfp,"batchTmpltFile = \"%s\"\n", response); 
    sprintf(pString,"What is the name of your target batch file ? ");
    getString(pString, response);
    response[strlen(response)-1] = '\0';
    fprintf(dfp,"batchFile = \"%s\"\n\n", response); 
    haveBatch = 1;
  }

  fprintf(dfp,"# ====> files to be copied to the working directories\n");
  if (genFlag == 1) nSupportFiles = 0;
  else
  {
    printOutTS(PL_INFO, 
       "INFO: Support files are files needed to run but not modified.\n");
    sprintf(pString,
         "How many other support files are needed to run the code ? ");
    nSupportFiles = getInt(0, 100, pString);
  }
  if (nSupportFiles > 0)
  {
    sprintf(pString,
       "Where can these support files be found(absolute path):");
    getString(pString, response);
    response[strlen(response)-1] = '/';
    response[strlen(response)] = '\0';
    fprintf(dfp,"supportDir = \"%s\"", response);
    fprintf(dfp,"supportFiles = [");
    for (ii = 0; ii < nSupportFiles; ii++)
    {
      sprintf(pString,
          "Enter enter the name (absolute path) of support file %d : ",
          ii+1);
      getString(pString, response);
      response[strlen(response)-1] = '\0';
      fprintf(dfp,"\"%s\"", response);
      if (ii < nSupportFiles-1) fprintf(dfp,", ");
    }
    fprintf(dfp,"]\n\n");
  }
  printOutTS(PL_INFO, 
        "=> Psuade assumes you use psub or msub to submit your job.\n");
  printOutTS(PL_INFO, 
        "=> If this is not true, please modify the runApplication\n");
  printOutTS(PL_INFO, "=> function after this is done.\n");
  fprintf(dfp,"#======================================================\n");
  fprintf(dfp,"# END USER SPECIFIC SECTION\n");
  fprintf(dfp,"#######################################################\n\n");

  fprintf(dfp,"#######################################################\n");
  fprintf(dfp,"# Function to get input data from PSUADE-generated\n");
  fprintf(dfp,"# parameter files (standard format, do not change). \n");
  fprintf(dfp,"# The return data will contain the inputs values.\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"def getInputData(inFileName):\n");
  fprintf(dfp,"   inFile  = open(inFileName, \"r\")\n");
  fprintf(dfp,"   lineIn  = inFile.readline()\n");
  fprintf(dfp,"   nCols   = lineIn.split()\n");
  fprintf(dfp,"   nInputs = eval(nCols[0])\n");
  fprintf(dfp,"   inputData = range(nInputs)\n");
  fprintf(dfp,"   for ind in range(nInputs):\n");
  fprintf(dfp,"      lineIn  = inFile.readline()\n");
  fprintf(dfp,"      nCols   = lineIn.split()\n");
  fprintf(dfp,"      inputData[ind] = eval(nCols[0])\n");
  fprintf(dfp,"   inFile.close()\n");
  fprintf(dfp,"   return inputData\n\n");

  fprintf(dfp,"#######################################################\n");
  fprintf(dfp,"# Function to generate input file\n");
  fprintf(dfp,"# Given an application input template file which has\n");
  fprintf(dfp,"# been modified to replace the parameter values with\n");
  fprintf(dfp,"# key words of your choice (input names you have\n");
  fprintf(dfp,"# entered), substitute those key words with the parameter\n");
  fprintf(dfp,"# values returned from the getInputData function and.\n");
  fprintf(dfp,"# put the modified file into appInputFile.\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"def genAppInputFile(inputData,appTmpltFile,appInputFile,\n");
  fprintf(dfp,"                    nInputs,inputNames):\n");
  fprintf(dfp,"   infile = open(appTmpltFile, \"r\")\n");
  fprintf(dfp,"   outfile = open(appInputFile, \"w\")\n");
  fprintf(dfp,"   while 1:\n");
  fprintf(dfp,"      lineIn  = infile.readline()\n");
  fprintf(dfp,"      if lineIn == \"\":\n");
  fprintf(dfp,"         break\n");
  fprintf(dfp,"      lineLen = len(lineIn)\n");
  fprintf(dfp,"      newLine = lineIn\n");
  fprintf(dfp,"      if nInputs > 0:\n");
  fprintf(dfp,"         for fInd in range(nInputs):\n");
  fprintf(dfp,"            strLen = len(inputNames[fInd])\n");
  fprintf(dfp,"            sInd = newLine.find(inputNames[fInd])\n");
  fprintf(dfp,"            if sInd >= 0:\n");
  fprintf(dfp,"               strdata = str(inputData[fInd])\n");
  fprintf(dfp,"               next = sInd + strLen\n");
  fprintf(dfp,"               lineTemp = newLine[0:sInd] + strdata + ");
  fprintf(dfp," \" \" + newLine[next:lineLen+1]\n");
  fprintf(dfp,"               newLine = lineTemp\n");
  fprintf(dfp,"               lineLen = len(newLine)\n");
  fprintf(dfp,"      outfile.write(newLine)\n");
  fprintf(dfp,"   infile.close()\n");
  fprintf(dfp,"   outfile.close()\n");
  fprintf(dfp,"   return\n\n");
                                                                               
  fprintf(dfp,"#######################################################\n");
  fprintf(dfp,"# Function to generate batch file\n");
  fprintf(dfp,"# If you submit a job via a batch file, this is the\n");
  fprintf(dfp,"# function to generate this file. Simply create a batch\n");
  fprintf(dfp,"# template file and substitute the counter inside used\n");
  fprintf(dfp,"# for distinguishing between different runs with the\n");
  fprintf(dfp,"# keyword PSUADE_COUNTER and the daisy-chained next run\n");
  fprintf(dfp,"# with the keyword PSUADE_NEXT. The modified batch file will\n");
  fprintf(dfp,"# be written to the file in batchFile.\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"def genBatchFile(batchTmpltFile,batchFile,fileTag,newTag):\n");
  fprintf(dfp,"   infile = open(batchTmpltFile, \"r\")\n");
  fprintf(dfp,"   outfile = open(batchFile, \"w\")\n");
  fprintf(dfp,"   searchString = \"PSUADE_COUNTER\"\n");
  fprintf(dfp,"   searchString2 = \"PSUADE_NEXT\"\n");
  fprintf(dfp,"   strLen = len(searchString)\n");
  fprintf(dfp,"   strLen2 = len(searchString2)\n");
  fprintf(dfp,"   while 1:\n");
  fprintf(dfp,"      lineIn  = infile.readline()\n");
  fprintf(dfp,"      if lineIn == \"\":\n");
  fprintf(dfp,"         break\n");
  fprintf(dfp,"      lineLen = len(lineIn)\n");
  fprintf(dfp,"      newLine = lineIn\n");
  fprintf(dfp,"      sInd    = newLine.find(searchString)\n");
  fprintf(dfp,"      if sInd > 0:\n");
  fprintf(dfp,"         next = sInd + strLen\n");
  fprintf(dfp,"         strdata = str(fileTag)\n");
  fprintf(dfp,"         lineTemp = newLine[0:sInd] + strdata + ");
  fprintf(dfp,"\" \" + newLine[next:lineLen+1]\n");
  fprintf(dfp,"         newLine = lineTemp\n");
  fprintf(dfp,"      sInd = newLine.find(searchString2)\n");
  fprintf(dfp,"      if sInd > 0:\n");
  fprintf(dfp,"         next = sInd + strLen2\n");
  fprintf(dfp,"         strdata = str(newTag)\n");
  fprintf(dfp,"         lineTemp = newLine[0:sInd] + strdata + ");
  fprintf(dfp,"\" \" + newLine[next:lineLen+1]\n");
  fprintf(dfp,"         newLine = lineTemp\n");
  fprintf(dfp,"      outfile.write(newLine)\n");
  fprintf(dfp,"   infile.close()\n");
  fprintf(dfp,"   outfile.close()\n");
  fprintf(dfp,"   return\n\n");

  fprintf(dfp,"#######################################################\n");
  fprintf(dfp,"# Function to run batch file\n");
  fprintf(dfp,"# This function submits the batch job using psub.\n");
  fprintf(dfp,"# Or, you can modify this function to submit jobs the\n");
  fprintf(dfp,"# way you want it.\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"def runApplication(executable):\n");
  fprintf(dfp,"   sysComm = \"/usr/gapps/XXXX/mach/xxxx.dev -r + \"\n");
#if 0
  fprintf(dfp,"#  sysComm = \"/usr/bin/srun -n 1 -N 1 -ppdebug + ");
  fprintf(dfp,"#executable inputdeck\n");
  fprintf(dfp,"#  tagString = os.path.splitext(executable)[1]\n");
  fprintf(dfp,"#  stringLen = len(tagString)\n");
  fprintf(dfp,"#  runNumber = eval(tagString[1:stringLen])\n");
  fprintf(dfp,"#  lastNo    = runNumber - %d\n", numJobs);
  fprintf(dfp,"#  prefix    = os.path.splitext(executable)[0];\n");
  fprintf(dfp,"#  lastJob   = prefix + \".\" + str(lastNo) + \" \"\n");
  fprintf(dfp,"#  statComm  = \"/usr/local/bin/pstat -D | /bin/grep \"\n");
  fprintf(dfp,"#  statComm  = statComm + lastJob + \" > stat\"\n");
  fprintf(dfp,"#  status = os.system(statComm); \n");
  fprintf(dfp,"#  dependString = \"\"\n");
  fprintf(dfp,"#  if status == 0:\n");
  fprintf(dfp,"#     statFile = open(\"stat\", \"r\")\n");
  fprintf(dfp,"#     lineIn = statFile.readline()\n");
  fprintf(dfp,"#     sIndex = lineIn.find(\"%s\")\n", userName);
  fprintf(dfp,"#     if sIndex > 0:\n");
  fprintf(dfp,"#        sIndex = lineIn.find(lastJob)\n");
  fprintf(dfp,"#        if sIndex > 0:\n");
  fprintf(dfp,"#           dependString = \"-d \" + lineIn[0:sIndex-1]\n");
  fprintf(dfp,"#     statFile.close()\n");
  fprintf(dfp,"#  sysComm = \"/usr/local/bin/psub \" + dependString + ");
  fprintf(dfp,"executable\n");
  fprintf(dfp, "def runApplication(executable):\n");
  fprintf(dfp, "   runComm  = \"/usr/bin/msub \" + batchFile\n");
  fprintf(dfp, "   os.system(statComm)\n");
  fprintf(dfp, "   waitFlag = 1\n");
  fprintf(dfp, "   while waitFlag == 1:\n");
  fprintf(dfp, "      statComm = \"/usr/bin/pstat | /bin/grep \" + ");
  fprintf(dfp, " batchFile + \" > stat\" \n");  
  fprintf(dfp, "      stat+ batchFile + " > stat"\n");
  fprintf(dfp, "      status = os.system(statComm)\n");
  fprintf(dfp, "      if os.path.isfile(\"stat\") != 0:\n");
  fprintf(dfp, "         statFile = open(\"stat\", \"r\")\n");
  fprintf(dfp, "         while 1:\n");
  fprintf(dfp, "            lineIn = statFile.readline()\n");
  fprintf(dfp, "            if lineIn == \"\":\n");
  fprintf(dfp, "               break\n");
  fprintf(dfp, "            sIndex = lineIn.find(batchFile)\n");
  fprintf(dfp, "            if sIndex > 0:\n");
  fprintf(dfp, "               waitFlag = 0\n");
  fprintf(dfp, "         statFile.close()\n");
  fprintf(dfp, "         os.remove(\"stat\")\n");
  fprintf(dfp, "      if (waitFlag == 1):\n");
  fprintf(dfp, "         os.system(\"sleep 100\")\n");
#endif
  fprintf(dfp,"   os.system(sysComm)\n");
  fprintf(dfp,"   return\n\n");

  fprintf(dfp,"#######################################################\n");
  fprintf(dfp,"# Function to generate output file \n");
  fprintf(dfp,"# This function writes the output data (which should\n");
  fprintf(dfp,"# have been generated in outData) to the PSUADE-based\n");
  fprintf(dfp,"# output file.\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"def genOutputFile(outFileName, outData):\n");
  fprintf(dfp,"   nLeng = len(outData)\n");
  fprintf(dfp,"   outfile = open(outFileName, \"w\")\n");
  fprintf(dfp,"   for ind in range(nLeng):\n");
  fprintf(dfp,"      outfile.write(\"%%e \\n\" %% outData[ind])\n");
  fprintf(dfp,"   outfile.close()\n");
  fprintf(dfp,"   return\n\n");
                                                                                
  fprintf(dfp,"#######################################################\n");
  fprintf(dfp,"# Main program file \n");
  fprintf(dfp,"#######################################################\n\n");
  fprintf(dfp,"inputfileName  = sys.argv[1]\n");
  fprintf(dfp,"outputfileName = sys.argv[2]\n\n");

  fprintf(dfp,"doTest        = 1\n");
  fprintf(dfp,"doAnalysis    = 0\n");
  fprintf(dfp,"doPreProcess  = 1\n");
  fprintf(dfp,"doPostProcess = 1\n\n");

  fprintf(dfp,"#======================================================\n");
  fprintf(dfp,"# check to see if output file already exists\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"if os.path.isfile(outputfileName) != 0:\n");
  fprintf(dfp,"   doTest = 0\n");
  fprintf(dfp,"   exit\n\n");

  fprintf(dfp,"#======================================================\n");
  fprintf(dfp,"# extract file tag (sample number)\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"tagString = os.path.splitext(inputfileName)[1]\n");
  fprintf(dfp,"stringLen = len(tagString)\n");
  fprintf(dfp,"fileTag   = eval(tagString[1:stringLen])\n\n");

  fprintf(dfp,"#======================================================\n");
  fprintf(dfp,"# check whether the working directory exists or not\n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"dirname = \"./workdir.\" + str(fileTag)\n");
  fprintf(dfp,"if os.path.isdir(dirname):\n");
  fprintf(dfp,"   doTest = 0\n");
  fprintf(dfp,"else:\n");
  fprintf(dfp,"   os.mkdir(dirname)\n\n");

  fprintf(dfp,"#======================================================\n");
  fprintf(dfp,"# if not, create directory and copy support files into it\n");
  fprintf(dfp,"# then set up and run the application \n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"if doTest == 1:\n");
  fprintf(dfp,"   if doPreProcess == 1:\n");
  fprintf(dfp,"      # copy files\n");
  fprintf(dfp,"      for file in appInputTmplts:\n");
  fprintf(dfp,"         file1 = appDir + file\n");
  fprintf(dfp,"         shutil.copy(file1, dirname)\n\n");
  fprintf(dfp,"      shutil.copy(inputfileName, dirname)\n\n");
  if (haveBatch != 0)
  {
    fprintf(dfp,"      file1 = appDir + batchTmpltFile\n");
    fprintf(dfp,"      shutil.copy(file1, dirname)\n\n");
  }
  if (nSupportFiles > 0)
  {
    fprintf(dfp,"      for file in supportFiles:\n");
    fprintf(dfp,"         file1 = supportDir + file\n");
    fprintf(dfp,"         shutil.copy(file1, dirname)\n\n");
  }
  fprintf(dfp,"   os.chdir(dirname)\n\n");
  fprintf(dfp,"   if doPreProcess == 1:\n");
  fprintf(dfp,"      inputData = getInputData(inputfileName)\n");
  fprintf(dfp,"      for ii in range (len(appInputTmplts)):\n");
  fprintf(dfp,"         genAppInputFile(inputData,appInputTmplts[ii],\n");
  fprintf(dfp,"                 appInputFiles[ii],nInputs,inputNames)\n");
  if (haveBatch == 1)
    fprintf(dfp,"      genBatchFile(batchTmpltFile, batchFile, fileTag)\n");
  if (genFlag == 1) 
  {
    nOutputs = 1;
    fprintf(dfp,"\n#  if doAnalysis == 1:\n");
    fprintf(dfp,"#     runApplication(batchFile)\n\n");
    fprintf(dfp,"#     runApplication(appInputFiles[0])\n\n");
  }
  else
  {
    fprintf(dfp,"\n   if doAnalysis == 1:\n");
    fprintf(dfp,"      runApplication(batchFile)\n\n");
    fprintf(dfp,"#     runApplication(appInputFiles[0])\n\n");
    sprintf(pString, "How many outputs are there ? ");
    nOutputs = getInt(1, 100000, pString);
  }
  fprintf(dfp,"   outData = range(%d)\n",nOutputs);
  fprintf(dfp,"   for ii in range(%d):\n",nOutputs);
  fprintf(dfp,"      outData[ii] = 1.0e35\n");
  if (genFlag == 1)
  {
    fprintf(dfp,"   outData[0] = inputData[0] - ");
    fprintf(dfp,"19.6 * inputData[1] / (1.5 * inputData[2])\n");
  }
  fprintf(dfp,"   genOutputFile(outputfileName, outData)\n");
  fprintf(dfp,"   shutil.copy(outputfileName, \"..\")\n\n");
  fprintf(dfp,"   os.chdir(\"..\")\n\n");

  fprintf(dfp,"#======================================================\n");
  fprintf(dfp,"# create dummy output files \n");
  fprintf(dfp,"#======================================================\n\n");
  fprintf(dfp,"else:\n");
  fprintf(dfp,"   outData = range(%d)\n",nOutputs);
  fprintf(dfp,"   for ii in range(%d):\n",nOutputs);
  fprintf(dfp,"      outData[ii] = 1.0e35\n");
  fprintf(dfp,"   genOutputFile(outputfileName, outData)\n");
  fprintf(dfp,"\n");
  if (genFlag == 1)
  {
    fprintf(dfp,"os.chdir(dirname)\n");
    fprintf(dfp,"filelist = os.listdir(\".\")\n");
    fprintf(dfp,"for file in filelist:\n");
    fprintf(dfp,"   if fnmatch.fnmatch(file, '*.core'):\n");
    fprintf(dfp,"      os.remove(file)\n");
    fprintf(dfp,"os.chdir(\"..\")\n");
    fprintf(dfp,"os.rmdir(dirname)\n\n");
  }
  else
  {
    fprintf(dfp,"#os.chdir(dirname)\n");
    fprintf(dfp,"#filelist = os.listdir(\".\")\n");
    fprintf(dfp,"#for file in filelist:\n");
    fprintf(dfp,"#   if fnmatch.fnmatch(file, '*.core'):\n");
    fprintf(dfp,"#      os.remove(file)\n");
    fprintf(dfp,"#os.chdir(\"..\")\n");
    fprintf(dfp,"#os.rmdir(dirname)\n\n");
  }
  printOutTS(PL_INFO, "PSUADE: the user driver is in %s\n", driverName);
  fclose(dfp);
  return 0;
}

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
int PsuadeBase::genSetup(int genFlag, char *filename)
{
  int    nInputs, iInd, nOutputs, samplingMethod, nSamples, nReps, ii;
  int    randomize = 0, jj, oInd;
  double *iLowerB, *iUpperB;
  char   dataFile[500], pString[500], **inputNames=NULL;
  char   **outputNames=NULL, winput[500];
  pData  pPtr, pINames, pLower, pUpper, pONames;
  PsuadeData psuadeIO;
  FILE   *fp, *fp2;

  if (genFlag == 0)
  {
    sprintf(pString,"Enter the name of the PSUADE input file : ");
    getString(pString, dataFile);
    dataFile[strlen(dataFile)-1] = '\0';
    if (!strcmp(dataFile, "\0")) 
    {
      printOutTS(PL_ERROR, "ERROR : invalid file name.\n");
      return 1;
    }
    sprintf(pString,"Enter the number of inputs (> 0) : ");
    nInputs = getInt(1, 10000000, pString);
  }
  else
  {
    printOutTS(PL_INFO,
       "PSUADE: an example PSUADE input file is in psuade.in.\n");
    strcpy(dataFile, "psuade.in");
    nInputs = 3;
  }
  iLowerB = new double[nInputs];
  iUpperB = new double[nInputs];
  inputNames = new char*[nInputs];
  for (iInd = 0; iInd < nInputs; iInd++) inputNames[iInd] = new char[200];
  if (genFlag == 0)
  {
    for (iInd = 0; iInd < nInputs; iInd++)
    {
      sprintf(pString, "name for input %d ? ", iInd+1);
      getString(pString, inputNames[iInd]);
      inputNames[iInd][strlen(inputNames[iInd])-1] = '\0';
      iLowerB[iInd] = iUpperB[iInd] = 0.0;
      while (iLowerB[iInd] >= iUpperB[iInd])
      {
        sprintf(pString, "lower bound for input %d ? ", iInd+1);
        iLowerB[iInd] = getDouble(pString);
        sprintf(pString,"upper bound for input %d (> lower bound)? ",
                iInd+1);
        iUpperB[iInd] = getDouble(pString);
      }
    }
  }
  else
  {
    strcpy(inputNames[0], "H0");
    strcpy(inputNames[1], "M");
    strcpy(inputNames[2], "S");
    iLowerB[0] = 40.0; iUpperB[0] = 60.0;
    iLowerB[1] = 67.0; iUpperB[1] = 74.0;
    iLowerB[2] = 20.0; iUpperB[2] = 40.0;
  }
  if (genFlag == 0)
  {
    printOutTS(PL_INFO, 
       "=> If you want to impose other than uniform probability \n");
    printOutTS(PL_INFO, 
       "=> density functions, you will have to modify the file\n");
    printOutTS(PL_INFO, "=> yourself after this session.\n");
    sprintf(pString,"Enter the number of outputs (> 0) : ");
    nOutputs = getInt(1, 100000, pString);
    outputNames = new char*[nOutputs];
    for (oInd = 0; oInd < nOutputs; oInd++) 
    {
       sprintf(pString, "name for output %d ? ", oInd+1);
       outputNames[oInd] = new char[200];
       getString(pString, outputNames[oInd]);
       outputNames[oInd][strlen(outputNames[oInd])-1] = '\0';
    }
  }
  else
  {
    nOutputs = 1;
    outputNames = new char*[nOutputs];
    outputNames[0] = new char[200];
    strcpy(outputNames[0], "H");
  }
  if (genFlag == 0)
  {
    samplingMethod = -1;
    while (samplingMethod < 0 || samplingMethod > 10)
    {
      printOutTS(PL_INFO,"Available sampling methods : \n");
      printOutTS(PL_INFO," MC    - Monte Carlo (random)\n");
      printOutTS(PL_INFO," FACT  - full factorial\n");
      printOutTS(PL_INFO," LH    - Latin Hypercube\n");
      printOutTS(PL_INFO," OA    - Orthogonal Array\n");
      printOutTS(PL_INFO," OALH  - OA-based Latin Hypercube\n");
      printOutTS(PL_INFO," MOAT  - Morris one at a time\n");
      printOutTS(PL_INFO," LPTAU - A Pseudo-random sequence\n");
      printOutTS(PL_INFO," METIS - A space-filling design\n");
      printOutTS(PL_INFO," FAST  - Fourier Amplitude Sampling Test\n");
      printOutTS(PL_INFO," FF4   - Fractional Factorial with Resolution IV\n");
      printOutTS(PL_INFO," FF5   - Fractional Factorial with Resolution V\n");
      sprintf(pString, "Sampling Method ? ");
      getString(pString, winput);
      winput[strlen(winput)-1] = '\0';
      if (!strcmp(winput,"MC"))         samplingMethod = 0;
      else if (!strcmp(winput,"FACT"))  samplingMethod = 1;
      else if (!strcmp(winput,"LH"))    samplingMethod = 2;
      else if (!strcmp(winput,"OA"))    samplingMethod = 3;
      else if (!strcmp(winput,"OALH"))  samplingMethod = 4;
      else if (!strcmp(winput,"MOAT"))  samplingMethod = 5;
      else if (!strcmp(winput,"LPTAU")) samplingMethod = 6;
      else if (!strcmp(winput,"METIS")) samplingMethod = 7;
      else if (!strcmp(winput,"FAST"))  samplingMethod = 8;
      else if (!strcmp(winput,"FF4"))   samplingMethod = 9;
      else if (!strcmp(winput,"FF5"))   samplingMethod = 10;
    }
    nReps = 1;
    if (samplingMethod >= 2 && samplingMethod <= 4)
    {
      printOutTS(PL_INFO, 
         "If you do not understand what number of replications\n");
      printOutTS(PL_INFO, "means, enter 1 for now.\n");
      sprintf(pString, "Number of replications (>= 1) ? ");
      nReps = getInt(1, 10000000, pString);
    }
    if (samplingMethod==0 || samplingMethod==6 || samplingMethod==7)
    {
      sprintf(pString, "Sample size (> 0) ? ");
      nSamples = getInt(1, 10000000, pString);
    }
    else if (samplingMethod == 1)
    {
      sprintf(pString,"Number of levels in each input (> 1) ? ");
      ii = getInt(1, 10000000, pString);
      nSamples = ii;
      for (jj = 1; jj < nInputs; jj++) nSamples *= ii;
    }
    else if (samplingMethod == 2)
    {
      nSamples = 0;
      while (nSamples <= 0 || nSamples%nReps != 0)
      {
        sprintf(pString,"sample size (multiple of %d) ? ", nReps);
        nSamples = getInt(1, 10000000, pString);
      }
    }
    else if (samplingMethod == 3 || samplingMethod == 4)
    {
      nSamples = 0;
      while (nSamples <= 0 || (nSamples%nReps != 0))
      {
        sprintf(pString,"number of levels per input (> 1, prime) ? ");
        ii = getInt(2, 1000000, pString);
        nSamples = ii * ii * nReps;
      }
      printOutTS(PL_INFO, "sample size = %d\n", nSamples);
    }
    else if (samplingMethod == 5)
    {
      sprintf(pString,"sample size (multiple of %d) ? ", nInputs+1);
      nSamples = getInt(nInputs+1, 10000, pString);
      nSamples = nSamples / (nInputs + 1) * (nInputs + 1);
    }
    else if (samplingMethod >= 8 && samplingMethod <=11)
    {
      nSamples = 100;
    }
    printOutTS(PL_INFO, "Sample size is set to %d.\n", nSamples);
    randomize = 0;
    if (samplingMethod >= 1 && samplingMethod <= 4)
    {
      sprintf(pString,
              "Do you want to add random noise to the sample? (y or n) ");
      getString(pString, winput);
      if (winput[0] == 'y') randomize = 1;
    }
    if (samplingMethod >= 2 && samplingMethod <= 4)
       if (nReps > 1) randomize = 0;
    if (samplingMethod == 8) randomize = 1;
  }
  else
  {
    samplingMethod = 0;
    nSamples = 1000;
    nReps = 1;
  }
  fp = fopen(dataFile, "w");
  if (fp == NULL)
  {
    printOutTS(PL_ERROR, "ERROR : cannot open file %s\n",dataFile);
    printOutTS(PL_ERROR, 
       "    Do you have write access in the current directory?\n");
  }
  else
  {
    fprintf(fp, "PSUADE\n");
    fprintf(fp, "INPUT\n");
    fprintf(fp, "   dimension = %d\n", nInputs);
    for (iInd = 0; iInd < nInputs; iInd++) 
       fprintf(fp, "   variable %3d %s = %16.8e %16.8e\n", iInd+1,
               inputNames[iInd],iLowerB[iInd],iUpperB[iInd]);
    fprintf(fp, "END\n");
    fprintf(fp, "OUTPUT\n");
    fprintf(fp, "   dimension = %d\n", nOutputs);
    for (oInd = 0; oInd < nOutputs; oInd++) 
    fprintf(fp, "   variable %3d %s\n", oInd+1, outputNames[oInd]);
    fprintf(fp, "END\n");
    fprintf(fp, "METHOD\n");
    if (samplingMethod == 0) fprintf(fp,"   sampling = MC\n");
    else if (samplingMethod == 1) fprintf(fp,"   sampling = FACT\n");
    else if (samplingMethod == 2) fprintf(fp,"   sampling = LH\n");
    else if (samplingMethod == 3) fprintf(fp,"   sampling = OA\n");
    else if (samplingMethod == 4) fprintf(fp,"   sampling = OALH\n");
    else if (samplingMethod == 5) fprintf(fp,"   sampling = MOAT\n");
    else if (samplingMethod == 6) fprintf(fp,"   sampling = LPTAU\n");
    else if (samplingMethod == 7) fprintf(fp,"   sampling = METIS\n");
    else if (samplingMethod == 8) fprintf(fp,"   sampling = FAST\n");
    else if (samplingMethod == 9) fprintf(fp,"   sampling = FF4\n");
    else if (samplingMethod == 10) fprintf(fp,"   sampling = FF5\n");
    fprintf(fp, "   num_samples = %d\n", nSamples);
    if (nReps > 1) fprintf(fp, "   num_replications = %d\n", nReps);
    if (randomize == 1) fprintf(fp, "   randomize\n");
    fprintf(fp, "#  for other options, consult manual\n");
    fprintf(fp, "END\n");
    fprintf(fp, "APPLICATION\n");
    fp2 = NULL;
    if (genFlag == 1)
    {
      fprintf(fp, "   driver = ./simulator\n");
      fprintf(fp, "#  driver = ./simulator.py\n");
    }
    else
    {
      while (fp2 == NULL)
      {
         printf("Enter the absolute path for the following.\n");
         sprintf(pString,
                 "Driver program name : (enter NONE if not needed) ");
         getString(pString, winput);
         winput[strlen(winput)-1] = '\0';
         if (!strcmp(winput, "NONE")) break;
         fp2 = fopen(winput, "r");
         if (fp2 == NULL) printOutTS(PL_ERROR,"file %s not found.\n",winput);
      } 
      if (fp2 != NULL) fclose(fp2);
      fprintf(fp, "   driver = %s\n", winput);
      fprintf(fp, "#  max_parallel_jobs = 1\n");
      fprintf(fp, "#  min_job_wait_time = 1\n");
      fprintf(fp, "#  max_job_wait_time = 1000000\n");
      fprintf(fp, "#  gen_inputfile_only\n");
      fprintf(fp, "#  launch_only\n");
      fprintf(fp, "#  launch_interval = 1\n");
      fprintf(fp, "#  for other options, consult manual\n");
    }
    fprintf(fp, "END\n");
    fprintf(fp, "ANALYSIS\n");
    fprintf(fp, "#  for analyzer and optimization options, consult manual\n");
    fprintf(fp, "#  analyzer method = Moment\n");
    fprintf(fp, "#  analyzer method = MainEffect\n");
    fprintf(fp, "#  analyzer rstype = MARS\n");
    fprintf(fp, "   diagnostics 1\n");
    fprintf(fp, "END\n");
    fprintf(fp, "END\n");
    fclose(fp); 
  }
  if (genFlag == 0)
    printOutTS(PL_INFO,"PSUADE: the PSUADE input file is in %s.\n",dataFile);
  if (inputNames != NULL)
  {
    for (ii = 0; ii < nInputs; ii++) delete [] inputNames[ii];
    delete [] inputNames;
  }
  if (outputNames != NULL)
  {
    for (ii = 0; ii < nOutputs; ii++) delete [] outputNames[ii];
    delete [] outputNames;
  }
  delete [] iLowerB;
  delete [] iUpperB;
  strcpy(filename, dataFile);
  return 0;
}

// ************************************************************************
// interpret command from interactive session
// ------------------------------------------------------------------------
int PsuadeBase::genWorkFlow()
{
   return 0;
}

