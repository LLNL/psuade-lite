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
// Functions for PsuadeConfig 
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sysdef.h"
#include "PsuadeConfig.h"
#include "PsuadeUtil.h"

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------ 
PsuadeConfig::PsuadeConfig(char *fname, int printLevel)
{
  int  ii, lineLeng=500, lineCnt;
  char lineIn[500], winput[500];
  FILE *fIn;

  fIn = fopen(fname, "r");
  if (fIn == NULL) 
  {
    printf("PsuadeConfig ERROR:: configure file %s not found.\n",fname);
    exit(1);
  }

  while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
  {
    sscanf(lineIn, "%s", winput);
    if (strcmp(winput, "PSUADE_CONFIG") == 0) break;
  }
  if (feof(fIn) != 0)
  {
    printf("PsuadeConfig ERROR:: keyword PSUADE_CONFIG not found.\n");
    fclose(fIn);
    exit(1);
  }
  lineCnt = 0;
  while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
  {
    sscanf(lineIn,"%s", winput);
    if   (strcmp(winput, "PSUADE_END") == 0) break;
    else lineCnt++;
  }
  if (strcmp(winput, "PSUADE_END") != 0)
  {
    printf("PsuadeConfig ERROR:: keyword PSUADE_END not found.\n");
    fclose(fIn);
    exit(1);
  }
  fclose(fIn);

  nLinesUsed_ = 0;
  nLines_ = 1000; 
  if (lineCnt > nLines_) nLines_ = lineCnt + 1000;
  fileData_ = new char*[nLines_];
  for (ii = 0; ii < nLines_; ii++)
  {
    fileData_[ii] = new char[1000];
    strcpy(fileData_[ii], "NONE");
  }
  if (lineCnt == 0) return;
  else
  {
    fIn = fopen(fname, "r");
    while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
    {
      sscanf(lineIn,"%s", winput);
      if (strcmp(winput, "PSUADE_CONFIG") == 0) break;
    }
    for (ii = 0; ii < lineCnt; ii++)
      fgets(fileData_[ii], lineLeng, fIn);
    nLinesUsed_ = lineCnt;
    fclose(fIn);
  }
  printLevel_ = printLevel;
}

// ************************************************************************
// constructor 
// ------------------------------------------------------------------------ 
PsuadeConfig::PsuadeConfig()
{
  printLevel_ = 0;
  nLinesUsed_ = 0;
  nLines_ = 1000; 
  fileData_ = new char*[nLines_];
  for (int ii = 0; ii < nLines_; ii++)
  {
    fileData_[ii] = new char[1000];
    strcpy(fileData_[ii], "NONE");
  }
}

// ************************************************************************
// copy constructor by Bill Oliver
// ------------------------------------------------------------------------
PsuadeConfig::PsuadeConfig(const PsuadeConfig & ps)
{
  fileData_ = NULL; 
  nLinesUsed_ = ps.nLinesUsed_;
  printLevel_ = ps.printLevel_;
  nLines_ = ps.nLines_;
  if (ps.fileData_ != NULL)
  { 
    fileData_ = new char *[nLines_];
    for (int i = 0; i < nLines_; i++)
    {
      fileData_[i] = new char[strlen(ps.fileData_[i] + 1)];
      strcpy(fileData_[i], ps.fileData_[i]); 
    }
  }
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------ 
PsuadeConfig::~PsuadeConfig()
{ 
  int ii;

  if (nLines_ > 0 && fileData_ != NULL)
  {
    for (ii = 0; ii < nLines_; ii++)
      if (fileData_[ii] != NULL) delete [] fileData_[ii];
    delete [] fileData_;
  }
}

// ************************************************************************
// request data from this object 
// ------------------------------------------------------------------------ 
char *PsuadeConfig::getParameter(const char *keyword)
{ 
  int  ii;
  char firstWord[100];

  for (ii = 0; ii < nLinesUsed_; ii++)
  {
    sscanf(fileData_[ii], "%s", firstWord);
    if (strcmp(keyword,firstWord) == 0)
    {
      if (printLevel_ > 1)
        printf("PsuadeConfig::parameter found: %s\n",fileData_[ii]);
      return fileData_[ii];
    }
  }
  return NULL;
}

// ************************************************************************
// add new data to this object 
// ------------------------------------------------------------------------ 
void PsuadeConfig::putParameter(const char *putLine)
{ 
  int leng;
  if (nLinesUsed_ >= nLines_)
  {
    printf("ERROR: cannot add lines to PsuadeConfig - FULL.\n");
    return;
  }
  leng = strlen(putLine);
  strncpy(fileData_[nLinesUsed_], putLine, leng);
  fileData_[nLinesUsed_][leng] = '\0';
  nLinesUsed_++;
  return;
}

// ************************************************************************
// remove data from this object 
// ------------------------------------------------------------------------ 
int PsuadeConfig::removeParameter(const char *keyword)
{ 
  int  ii;
  char firstWord[100];

  for (ii = 0; ii < nLinesUsed_; ii++)
  {
    sscanf(fileData_[ii], "%s", firstWord);
    if (strcmp(keyword,firstWord) == 0)
    {
       strcpy(fileData_[ii], "NONE");
       return 1;
    }
  }
  return 0;
}


// ************************************************************************
// write to file 
// ------------------------------------------------------------------------ 
void PsuadeConfig::writeToFile(char *fname)
{
  int  ii;
  FILE *fOut;

  fOut = fopen(fname, "w");
  if (fOut == NULL) 
  {
    printf("PsuadeConfig ERROR:: cannot write to configure file %s.\n",
           fname);
    exit(1);
  }

  fprintf(fOut, "PSUADE_CONFIG\n");
  for (ii = 0; ii < nLines_; ii++) fprintf(fOut, "%s\n", fileData_[ii]);
  fprintf(fOut, "PSUADE_END\n");
  fclose(fOut);
}

// ************************************************************************
// print content of the config object
// ------------------------------------------------------------------------ 
void PsuadeConfig::print()
{
  printAsterisks(PL_INFO, 0);
  printf("************* PSUADE configuration information\n");
  printEquals(PL_INFO, 0);
  for (int ii = 0; ii < nLines_; ii++)
    if (strcmp(fileData_[ii], "NONE")) printf("%s\n", fileData_[ii]);
  printAsterisks(PL_INFO, 0);
}

// ************************************************************************
// add information from a configure file
// ------------------------------------------------------------------------ 
void PsuadeConfig::addFromFile(char *fname)
{
  int  ii, lineLeng=500, lineCnt, lastLines;
  char lineIn[500], winput[500], **fdata;
  FILE *fIn;

  fIn = fopen(fname, "r");
  if (fIn == NULL) 
  {
    printf("PsuadeConfig ERROR:: configure file %s not found.\n",fname);
    exit(1);
  }

  while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
  {
    sscanf(lineIn, "%s", winput);
    if (strcmp(winput, "PSUADE_CONFIG") == 0) break;
  }
  if (feof(fIn) != 0)
  {
    printf("PsuadeConfig ERROR:: keyword PSUADE_CONFIG not found.\n");
    fclose(fIn);
    exit(1);
  }
  lineCnt = 0;
  while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
  {
    sscanf(lineIn,"%s", winput);
    if   (strcmp(winput, "PSUADE_END") == 0) break;
    else lineCnt++;
  }
  if (strcmp(winput, "PSUADE_END") != 0)
  {
    printf("PsuadeConfig ERROR:: keyword PSUADE_END not found.\n");
    fclose(fIn);
    exit(1);
  }
  fclose(fIn);

  if (lineCnt+nLinesUsed_ > nLines_) 
  {
    lastLines = nLines_;
    nLines_ = lineCnt + nLinesUsed_ + 1000;
    fdata = fileData_;
    fileData_ = new char*[nLines_];
    for (ii = 0; ii < lastLines; ii++) fileData_[ii] = fdata[ii];
    delete [] fdata;
    for (ii = lastLines; ii < nLines_; ii++)
    {
      fileData_[ii] = new char[1000];
      checkAllocate(fileData_[ii],"fileData in Config::addFromFile");
      strcpy(fileData_[ii], "NONE");
    }
  }
  if (lineCnt == 0) return;
  else
  {
    fIn = fopen(fname, "r");
    while ((fgets(lineIn, lineLeng, fIn) != NULL) && (feof(fIn) == 0))
    {
      sscanf(lineIn,"%s", winput);
      if (strcmp(winput, "PSUADE_CONFIG") == 0) break;
    }
    for (ii = 0; ii < lineCnt; ii++)
      fgets(fileData_[nLinesUsed_+ii], lineLeng, fIn);
    nLinesUsed_ += lineCnt;
    fclose(fIn);
  }
}

// ************************************************************************
// friend function (create a function approximator given a file name)
// perform PDF transformation
// check invalid sample points
// RS type from file 
// ------------------------------------------------------------------------
extern "C"
int genConfigFileTemplate(char *fname)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if (fp != NULL)
  {
    fprintf(fp, "# Use this file if you need to set options but\n");
    fprintf(fp, "# but it is cumbersome to use, e.g., rs_expert mode.\n");
    fprintf(fp, "# Current options are listed in the following\n");
    fprintf(fp, "# Uncomment and add information to activate them.\n");
    fprintf(fp, "# To use this, add this line in ANALYSIS section\n");
    fprintf(fp, "#       use_configure_file = <this file name>\n");
    fprintf(fp, "PSUADE_CONFIG\n");
    fprintf(fp, "## Normalize inputs for some response surface methods\n");
    fprintf(fp, "#normalize_input (take out the # to turn on)\n");
    fprintf(fp, "## Normalize outputs for some response surface methods\n");
    fprintf(fp, "#normalize_output (take out the # to turn on)\n");
    fprintf(fp, "## MARS parameters (take out the # to turn on)\n");
    fprintf(fp, "#MARS_num_bases = 50\n");
    fprintf(fp, "#MARS_interaction = 2\n");
    fprintf(fp, "## SVM parameters (take out the # to turn on)\n");
    fprintf(fp, "#SVM_gamma = 1.0 (1e-6 - 1.0)\n");
    fprintf(fp, "#SVM_tol = 1.0 (1e-6 - 1.0)\n");
    fprintf(fp, "#SVM_kernel = 1 (1:linear, 2:cubic, 3:RBF, 4:sigmoid)\n");
    fprintf(fp, "## Kriging parameters (take out the # to turn on)\n");
    fprintf(fp, "#KRI_mode = 2\n");
    fprintf(fp, "#KRI_tol = 1.0e-6\n");
    fprintf(fp, "#KRI_DATA_STDEV_FILE = <add a file here>\n");
    fprintf(fp, "#KRI_LENG_SCALE 1 = 0.1\n");
    fprintf(fp, "#KRI_LENG_SCALE 2 = 0.1\n");
    fprintf(fp, "## Legendre parameters (take out the # to turn on)\n");
    fprintf(fp, "#Legendre_order = 1\n");
    fprintf(fp, "## RBF,GP override\n");
    fprintf(fp, "RS_no_multi_domain\n");
    fprintf(fp, "## multi-domain response surface method parameters\n");
    fprintf(fp, "MRBF_max_samples_per_group = 1000\n");
    fprintf(fp, "MGP_max_samples_per_group = 1000\n");
    fprintf(fp, "## MOAT parameter (number of levels)\n");
    fprintf(fp, "#MOAT_P = 4\n");
    fprintf(fp, "## GMOAT parameter (number of levels)\n");
    fprintf(fp, "#GMOAT_P = 4\n");
    fprintf(fp, "## RS-based Sobol index parameters\n");
    fprintf(fp, "#RSMSobol1_nsubsamples = 1000\n");
    fprintf(fp, "#RSMSobol1_nlevels = 200\n");
    fprintf(fp, "#RSMSobol2_nsubsamples = 1000\n");
    fprintf(fp, "#RSMSobol2_nlevels = 200\n");
    fprintf(fp, "#RSMSoboltsi_nsubsamples = 1000\n");
    fprintf(fp, "#RSMSoboltsi_nlevels = 100\n");
    fprintf(fp, "#RSMSoboltG_nsubsamples_ingroup = 500\n");
    fprintf(fp, "#RSMSoboltG_nsubsamples_outgroup = 2000\n");
    fprintf(fp, "PSUADE_END\n");
    fclose(fp);
    return 0;
  }
  else return -1;
}

