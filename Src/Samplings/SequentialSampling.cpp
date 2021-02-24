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
// Functions for the SequentialSampling class 
// AUTHOR : CHARLES TONG
// DATE   : 2017
// ************************************************************************

// ------------------------------------------------------------------------
// system and local includes
// ------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sysdef.h"
#include "PsuadeUtil.h"
#include "SequentialSampling.h"
#include "GP3.h"

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
SequentialSampling::SequentialSampling() : Sampling()
{
  samplingID_ = PSUADE_SAMP_SEQ;
  printAsterisks(PL_INFO, 0);
  printf("SequentialSampling:\n");
  printf("This module performs batch sequential experimental design.\n");
  printf("To track progress and get matlab plots, set print level > 1.\n");
  printDashes(PL_INFO, 0);
}

// ************************************************************************
// destructor 
// ------------------------------------------------------------------------
SequentialSampling::~SequentialSampling()
{
}

// ************************************************************************
// initialize the sampling data
// ------------------------------------------------------------------------
int SequentialSampling::initialize(int initLevel)
{
  int    ii, jj, kk, status, nTests, nDesigns, *SS, option1, option2;
  double ddata, *SY, quality[3];
  char   pString[1000], fname[10000];
  FILE   *fp;
  GP3    *rsPtr=NULL;
  psVector gpParams;
  psMatrix testMatrix, finalDesign;
  Sampling *sampler;

  if (nInputs_ == 0)
  {
    printf("SequentialSampling::initialize ERROR - input not set up.\n");
    printf("                    Did you call setInputBounds?\n");
    exit(1);
  }

  printf("SequentialSampling::initialization.\n");
  printf("This method requires the GP hyperparameters in search of\n");
  printf("new sample points. The GP hyperparameters can either be\n");
  printf("computed by PSUADE when a sample is provided by users, or,\n");
  printf("alternatively, users can supply the hyperparameters directly.\n");
  if (vecSamInps_.length() > 0)
  {
    if (printLevel_ > 1)
    {
      printf("SequentialSampling::initialize: create response surface\n");
      printf("            RS type  = GP\n");
      printf("            nSamples = %d\n", nSamples_);
      printf("            nInputs  = %d\n", nInputs_);
    }
    if (nSamples_ < 10)
    {
      printf("SequentialSampling::ERROR - nSamples must be >= 10\n");
      printf("                    RS may not be good for small samples.\n");
      return -1;
    }
    printf("A sample has been provided ==> \n");
    printf("         PSUADE will compute GP3 hyperparameters.\n");
    rsPtr = new GP3(nInputs_, nSamples_);
    rsPtr->setBounds(vecLBs_.getDVector(), vecUBs_.getDVector());
    rsPtr->setOutputLevel(printLevel_);
    status = rsPtr->initialize(vecSamInps_.getDVector(),
                               vecSamOuts_.getDVector());
    rsPtr->getHyperparameters(gpParams);

    for (ii = 0; ii < nInputs_; ii++)
      gpParams[ii] = 0.5 / exp(2.0 * gpParams[ii]);
    gpParams[nInputs_] = 1.0;
    gpParams[nInputs_+1] = exp(gpParams[nInputs_+1]);
    gpParams[nInputs_+1] *= sqrt(nSamples_);
    gpParams[nInputs_+2] = 0;
    gpParams[nInputs_+3] = exp(gpParams[nInputs_+3]); 
    delete rsPtr;
  }
  else
  {
    printf("A sample has NOT been provided ==> \n");
    printf("Please enter the hyperparameters (scale) parameters.\n");
    printf("Scale used in exponential: exp(-sum scale_i * dist_i ^ 2).\n");
    gpParams.setLength(nInputs_+4);
    for (ii = 0; ii < nInputs_; ii++)
    {
      sprintf(pString,"Enter hyperparameter for input %d (e.g. 0.1): ",ii+1);
      gpParams[ii] = 0.0;
      while (gpParams[ii] <= 0.0)
      {
        gpParams[ii] = getDouble(pString);
        if (gpParams[ii] <= 0) printf("ERROR: value needs to be positive.\n");
      }
    }
    gpParams[nInputs_] = 1.0;
    gpParams[nInputs_+1] = 0;
    gpParams[nInputs_+2] = 0;
    gpParams[nInputs_+3] = -1.0;
    sprintf(pString,"Enter nugget (added to Cdiag) value (> 0, e.g. 0.1) : ");
    while (gpParams[nInputs_+3] < 0)
    {
      gpParams[nInputs_+3] = getDouble(pString);
      if (gpParams[nInputs_+3] < 0)
        printf("ERROR: Nugget should be >= 0.\n");
    }
  }

  printf("NEXT, a candidate design set (a set of sample points from\n");
  printf("which the final design point will be selected) needs to be\n");
  printf("created by PSUADE or provided by users.\n");
  printf("Please select whether \n");
  printf("1. PSUADE is to create a candidate design set, or\n");
  printf("2. You will provide a candidate design set\n");
  sprintf(pString, "Enter you choice : (1 or 2) ");
  kk = getInt(1,2,pString);
  if (kk == 1)
  {
    printf("PSUADE is to create a candidate design set.\n");
    printf("Select candidate design method:\n");
    printf("1. Full Factorial\n");
    printf("2. LPtau (LHS if nInputs > 51)\n");
    sprintf(pString, "Enter candidate design method (1 or 2) : ");
    kk = getInt(1,2,pString);
    if (kk == 1) sampler = SamplingCreateFromID(PSUADE_SAMP_FACT);
    else
    {
      if (nInputs_ > 51) sampler = SamplingCreateFromID(PSUADE_SAMP_LHS);
      else               sampler = SamplingCreateFromID(PSUADE_SAMP_LPTAU);
    }
    sprintf(pString, "Enter candidate design sample size (>= 25): ");
    nTests = getInt(25,1000000,pString);
    sampler->setInputBounds(nInputs_, vecLBs_.getDVector(), 
                            vecUBs_.getDVector());
    sampler->setOutputParams(1);
    sampler->setSamplingParams(nTests, 1, 0);
    sampler->initialize(0);
    psVector  vecSX, vecSY;
    psIVector vecSS;
    vecSX.setLength(nTests*nInputs_);
    vecSY.setLength(nTests);
    vecSS.setLength(nTests);
    sampler->getSamples(nTests, nInputs_, 1, vecSX.getDVector(), 
                        vecSY.getDVector(), vecSS.getIVector());
    delete sampler;
    testMatrix.setDim(nTests, nInputs_);
    for (ii = 0; ii < nTests; ii++)
      for (jj = 0; jj < nInputs_; jj++)
        testMatrix.setEntry(ii,jj,vecSX[ii*nInputs_+jj]);
  }
  else
  {
    printf("The candidate design set should be stored in a text file\n");
    printf("having the following format: \n");
    printf("<number of points> <number of inputs>\n");
    printf("1 input values \n");
    printf("2 input values \n");
    printf(".... \n");
    sprintf(pString, "Enter the file name of your candidate design set : ");
    getString(pString, fname);
    kk = strlen(fname);
    fname[kk-1] = '\0';
    status = readStdDataFile(fname, testMatrix);
    if (status != 0)
    {
      printf("SequentialSampling::initialize: ERROR in reading test data\n");
      exit(1);
    }
  }
  nTests = testMatrix.nrows();
  nDesigns = nTests - 1;
  sprintf(pString,
    "How many design points to be selected from the candidate set? (1-%d) ", 
    nDesigns);
  nDesigns = getInt(1, nDesigns, pString);
  //if (nDesigns != nDesigns / 2 * 2)
  //{
  //  nDesigns = (nDesigns + 1) / 2 * 2;
  //  printf("nDesigns rounded to even number ==> %d\n", nDesigns);
  //}

  if (nSamples_ > 0)
  {
    XMeans_.setLength(nInputs_);
    XStds_.setLength(nInputs_);
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++) 
        ddata += vecSamInps_[jj*nInputs_+ii];
      XMeans_[ii] = ddata / (double) nSamples_;
      ddata = 0.0;
      for (jj = 0; jj < nSamples_; jj++)
        ddata += pow(vecSamInps_[jj*nInputs_+ii] - XMeans_[ii], 2.0);
      XStds_[ii] = sqrt(ddata / (double) (nSamples_ - 1));
      if (XStds_[ii] == 0.0) XStds_[ii] = 1.0;
    }
    for (jj = 0; jj < nTests; jj++)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = testMatrix.getEntry(jj,ii);
        ddata = (ddata - XMeans_[ii]) / XStds_[ii];
        testMatrix.setEntry(jj,ii, ddata);
      }
    }
  }
  else
  {
    XMeans_.setLength(nInputs_);
    XStds_.setLength(nInputs_);
    for (ii = 0; ii < nInputs_; ii++) XMeans_[ii] = 0;
    for (ii = 0; ii < nInputs_; ii++) XStds_[ii] = 1;
  } 
 
  psVector markers;
  printf("The procedure goes like this:\n");
  printf("I.  An initial sample will be created by PSUADE\n");
  printf("II. An insert-and-delete algorithm will be used to finetune\n");
  printf("For Step (I), you can choose between the following 4 options:\n");
  printf("1. Pick the first %d sample points from the candidate set\n",
         nDesigns);
  printf("2. Randomly select %d sample points from the candidate set\n",
         nDesigns);
  printf("3. Use sequential sampling to select sample points\n");
  printf("4. Use pairwise sequential sampling to select sample points\n");
  printf("5. Use optimal sequential sampling (time-consuming)\n");
  sprintf(pString, "Please make your selection (1, 2, 3, 4, or 5) : ");
  option1 = getInt(1, 5, pString);
  if (option1 != 5)
  {
    printf("For Step (II), you can choose between the following 2 options:\n");
    printf("1. Use add-delete algorithm\n");
    printf("2. Use swap algorithm (more expensive, better performance)\n");
    sprintf(pString, "Please make your selection (1 or 2) : ");
    option2 = getInt(1, 2, pString);
  }
  else option2 = 2;

  markers.setLength(nTests); 
  for (ii = 0; ii < nTests; ii++) markers[ii] = -1;
  finalDesign.setDim(nDesigns, nInputs_);

  fp = fopen("matlabsdoe.m", "r");
  if (fp != NULL)
  { 
    fclose(fp);
    fp = fopen("matlabsdoe.m", "w");
    fclose(fp);
  }
  else if (printLevel_ > 0)
  {
    fp = fopen("matlabsdoe.m", "w");
    fclose(fp);
  }
  printf("Phase 1 analysis\n");
  if (option1 == 1 || option1 == 2)
  {
    if (printLevel_ > 0)
    {
      fp = fopen("matlabsdoe.m", "w");
      fprintf(fp,"X = [\n");
    }
    for (jj = 0; jj < nDesigns; jj++)
    {
      if (option1 == 2)
      {
        kk = PSUADE_rand() % nTests;
        while (markers[kk] != -1) kk = PSUADE_rand() % nTests;
      }
      else kk = jj;
      if (printLevel_ > 0) fprintf(fp,"%4d ",jj+1);
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = testMatrix.getEntry(kk,ii);
        finalDesign.setEntry(jj,ii, ddata);
        if (printLevel_ > 0) 
          fprintf(fp,"%12.4e ",ddata*XStds_[ii]+XMeans_[ii]);
      }
      markers[kk] = jj;
      if (printLevel_ > 0) fprintf(fp,"\n");
    }
    if (printLevel_ > 0)
    {
      fprintf(fp,"];\n");
      fprintf(fp,"X = X(1:%d,:);\n",nDesigns);
      fprintf(fp,
        "plot(X(:,2),X(:,3),'rh','markersize',13,'markerfacecolor','r')\n");
      fprintf(fp,"axis([%e %e %e %e])\n",vecLBs_[0],vecUBs_[0],
              vecLBs_[1],vecUBs_[1]);
      fclose(fp);
    }
  }
  else if (option1 == 3)
  {
    genInitialDesigns(gpParams,testMatrix,nDesigns,finalDesign,markers);
  }
  else if (option1 == 4)
  {
    genInitialDesigns2(gpParams,testMatrix,nDesigns,finalDesign,markers);
  }
  else if (option1 == 5)
  {
    genDesignsUltimate(gpParams,testMatrix,nDesigns,finalDesign);
  }
  if (option1 != 5 && printLevel_ > 0)
  {
    fp = fopen("matlabsdoe.m", "a");
    fprintf(fp, "disp('Initial sample generated.')\n");
    fprintf(fp, "disp('Press enter to continue')\n");
    fclose(fp);
  }

  psVector vecSamInps;
  vecSamInps.setLength(nDesigns*nInputs_);
  for (jj = 0; jj < nDesigns; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = finalDesign.getEntry(jj,ii);
      vecSamInps[jj*nInputs_+ii] = ddata;
    }
  }
  SamplingQuality(nDesigns,nInputs_,vecSamInps.getDVector(),
                  vecLBs_.getDVector(),vecUBs_.getDVector(), quality);
  printf("Sample quality after initial design\n");
  printf("Min-min distance (sample-sample ) = %e (large is good)\n",
         quality[0]);
  printf("Avg-avg distance (sample pairs  ) = %e (large is good)\n",
         quality[2]);
  printf("Avg-min distance (corner-samples) = %e (small is good)\n",
         quality[1]);

  if (option1 != 5) 
  {
    printf("Phase 2 analysis\n");
    if (option2 == 1)
      genDesigns(gpParams, testMatrix, nDesigns, finalDesign, markers);
    else if (option2 == 2)
      genDesigns2(gpParams, testMatrix, nDesigns, finalDesign, markers);

    vecSamInps.setLength(nDesigns*nInputs_);
    for (jj = 0; jj < nDesigns; jj++)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = finalDesign.getEntry(jj,ii);
        vecSamInps[jj*nInputs_+ii] = ddata;
      }
    }
    SamplingQuality(nDesigns,nInputs_,vecSamInps.getDVector(),
                    vecLBs_.getDVector(),vecUBs_.getDVector(), quality);
    printf("Sample quality after final design\n");
    printf("Min-min distance (sample-sample ) = %e (large is good)\n",
           quality[0]);
    printf("Avg-avg distance (sample pairs  ) = %e (large is good)\n",
           quality[2]);
    printf("Avg-min distance (corner-samples) = %e (small is good)\n",
           quality[1]);
  }
  if (printLevel_ > 0)
    printf("INFO: design update history is in matlabsdoe.m\n");

  if (nSamples_ > 0)
  {
    kk = finalDesign.nrows();
    for (jj = 0; jj < kk; jj++)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = finalDesign.getEntry(jj,ii);
        ddata = ddata * XStds_[ii] + XMeans_[ii];
        finalDesign.setEntry(jj,ii, ddata);
      }
    }
  }
  return 0;
}

// ************************************************************************
// construct covariance matrix given hyperparameters
// ------------------------------------------------------------------------
void SequentialSampling::constructCMatrix(psMatrix &sampleMatrix,
                                       psMatrix &CMatrix,psVector params)
{
  int    ii, jj, kk, nrows, ncols;
  double dist, dtmp1, dtmp2;

  nrows = sampleMatrix.nrows();
  ncols = sampleMatrix.ncols();

  CMatrix.setDim(nrows, nrows);
  for (jj = 0; jj < nrows; jj++)
  {
    dist = params[nInputs_] + params[nInputs_+1] + params[nInputs_+3];
    CMatrix.setEntry(jj,jj,dist);
    for (kk = jj+1; kk < nrows; kk++)
    {
      dist = 0.0;
      for (ii = 0; ii < nInputs_; ii++)
      {
        dtmp1 = sampleMatrix.getEntry(jj,ii);
        dtmp2 = sampleMatrix.getEntry(kk,ii);
        dtmp1 = dtmp1 - dtmp2;
        dist += pow(dtmp1, 2.0) * params[ii];
      }
      dist = params[nInputs_] * exp(-dist);
      if (dist < 1.0e-50) dist = 0;
      CMatrix.setEntry(jj,kk,dist);
      CMatrix.setEntry(kk,jj,dist);
    }
  }
  return;
}

// ************************************************************************
// generate design method I 
// ------------------------------------------------------------------------
int SequentialSampling::genDesigns(psVector &params, psMatrix &testMatrix, 
                                   int nDesigns, psMatrix &finalDesign,
                                   psVector &markers)
{
  int    ii, jj, kk, mm, nn, status, nTests, minInd, markerSave, iter;
  double ddata, dist, meas, minVal, maxVal;
  FILE   *fp=NULL;
  psMatrix CMatrix, samMatrix, baseDesign;
  psVector correlations, tmpVec;

  nTests = testMatrix.nrows();
  baseDesign.setDim(nDesigns+1, nInputs_);
  for (jj = 0; jj < nDesigns; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = finalDesign.getEntry(jj,ii);
      baseDesign.setEntry(jj, ii, ddata);
    }
  }

  iter = 0;
  while (iter < 100)
  {
    iter++;
    samMatrix.setDim(nDesigns+1, nInputs_);
    for (jj = 0; jj < nDesigns; jj++)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = baseDesign.getEntry(jj,ii);
        samMatrix.setEntry(jj, ii, ddata);
      }
    }

    minVal = PSUADE_UNDEFINED;
    minInd = -1;
    for (mm = 0; mm < nTests; mm++)
    {
      if (markers[mm] == -1)
      {
        for (ii = 0; ii < nInputs_; ii++)
        {
          ddata = testMatrix.getEntry(mm,ii);
          samMatrix.setEntry(nDesigns, ii, ddata);
        }
        markers[mm] = nDesigns;

        constructCMatrix(samMatrix, CMatrix, params);
        status = CMatrix.computeInverse(CMatrix);
        if (status != 0)
        {
          printf("SequentialSampling ERROR in covariance inversion\n");
          exit(1);
        }

        correlations.setLength(nDesigns+1);
        maxVal = - PSUADE_UNDEFINED;
        for (nn = 0; nn < nTests; nn++)
        {
          if (markers[nn] == -1)
          {
            for (jj = 0; jj < nDesigns+1; jj++)
            {
              dist = 0.0;
              for (ii = 0; ii < nInputs_; ii++)
              {
                ddata = samMatrix.getEntry(jj, ii);
                ddata -= testMatrix.getEntry(nn,ii);
                dist += pow(ddata, 2.0) * params[ii];
              }
              dist = exp(-dist);
              if (dist < 1.0e-50) dist = 0;
              correlations[jj] = dist;
            }
            CMatrix.matvec(correlations, tmpVec, 0);
            meas = 0.0;
            for (jj = 0; jj <= nDesigns; jj++)
              meas += correlations[jj]*tmpVec[jj];
            meas = -meas;
            if (meas > maxVal) maxVal = meas;
          }
        }
        if (maxVal < minVal)
        {
          minVal = maxVal;
          minInd = mm;
        }
        else if (maxVal == minVal)
        {
          nn = PSUADE_rand() % 2;
          if (nn == 1)
          {
            minVal = maxVal;
            minInd = mm;
          }
        }
        markers[mm] = -1;
      }
    }
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = testMatrix.getEntry(minInd,ii);
      baseDesign.setEntry(nDesigns, ii, ddata);
    } 
    markers[minInd] = nDesigns;

    minVal = PSUADE_UNDEFINED;
    minInd = -1;
    samMatrix.setDim(nDesigns, nDesigns);
    for (mm = 0; mm <= nDesigns; mm++)
    {
      nn = 0;
      for (kk = 0; kk <= nDesigns; kk++)
      {
        if (kk != mm)
        {
          for (ii = 0; ii < nInputs_; ii++)
          {
            ddata = baseDesign.getEntry(kk,ii);
            samMatrix.setEntry(nn, ii, ddata);
          }
          nn++;
        }
      }
      for (kk = 0; kk < nTests; kk++)
        if (markers[kk] == mm) break;
      markerSave = kk;
      markers[kk] = -1;

      constructCMatrix(samMatrix, CMatrix, params);
      status = CMatrix.computeInverse(CMatrix);
      if (status != 0)
      {
        printf("SequentialSampling ERROR in covariance inversion\n");
        exit(1);
      }

      correlations.setLength(nDesigns);
      maxVal = - PSUADE_UNDEFINED;
      for (nn = 0; nn < nTests; nn++)
      {
        if (markers[nn] == -1)
        {
          for (jj = 0; jj < nDesigns; jj++)
          {
            dist = 0.0;
            for (ii = 0; ii < nInputs_; ii++)
            {
              ddata = samMatrix.getEntry(jj, ii);
              ddata -= testMatrix.getEntry(nn,ii);
              dist += pow(ddata, 2.0) * params[ii];
            }
            dist = exp(-dist);
            if (dist < 1.0e-50) dist = 0;
            correlations[jj] = dist;
          }
          CMatrix.matvec(correlations, tmpVec, 0);
          meas = 0.0;
          for (jj = 0; jj < nDesigns; jj++)
            meas += correlations[jj]*tmpVec[jj];
          meas = -meas;
          if (meas > maxVal) maxVal = meas;
        }
      }
      if (maxVal < minVal)
      {
        minVal = maxVal;
        minInd = mm;
      }
      markers[markerSave] = mm;
    }
    if (printLevel_ > 1)
      printf("genDesigns iteration %d: best value = %e\n",iter,minVal);

    if (minInd != nDesigns)
    {
      for (kk = 0; kk < nTests; kk++)
        if (markers[kk] == minInd) break;
      markers[kk] = -1;
      if (printLevel_ > 1)
      {
        printf("genDesigns iteration %4d: sample point deleted = %d\n",
               iter,kk);
      }
      for (kk = 0; kk < nTests; kk++)
        if (markers[kk] == nDesigns) break;
      markers[kk] = minInd;
      if (printLevel_ > 1)
      {
        printf("                           sample point added   = %d\n",
               kk);
      }
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = baseDesign.getEntry(nDesigns,ii);
        baseDesign.setEntry(minInd, ii, ddata);
      } 
      if (printLevel_ > 0)
      {
        fp = fopen("matlabsdoe.m", "a");
        fprintf(fp,"pause\n");
        fprintf(fp,"X = [\n");
        printDesigns(baseDesign, fp);
        fprintf(fp,"];\n");
        fprintf(fp,"X = X(1:%d,:);\n",nDesigns);
        fprintf(fp,
          "plot(X(:,2),X(:,3),'rh','markersize',13,'markerfacecolor','r')\n");
        fprintf(fp,"axis([%e %e %e %e])\n",vecLBs_[0],vecUBs_[0],
                vecLBs_[1],vecUBs_[1]);
        fprintf(fp,"grid on\n");
        fprintf(fp,"box on\n");
        fprintf(fp,"set(gca, 'linewidth',2)\n");
        fprintf(fp,"set(gca, 'fontsize',12)\n");
        fprintf(fp,"xlabel('X1','fontsize',12,'fontweight','bold')\n");
        fprintf(fp,"ylabel('X2','fontsize',12,'fontweight','bold')\n");
        fprintf(fp,"disp('Press enter to continue')\n");
        fclose(fp);
      }
    }
    else
    {
      for (kk = 0; kk < nTests; kk++)
        if (markers[kk] == nDesigns) break;
      markers[kk] = -1;
      if (printLevel_ > 1)
        printf("genDesigns iteration %4d: no more point exchanged\n",iter);
      break;
    }

    printf("Iter %d: Selected test sample points so far: (1-based)\n", iter);
    jj = 0;
    for (ii = 0; ii < nTests; ii++)
    {
      if (markers[ii] != -1)
      {
        printf("%4d ", ii+1); 
        jj++;
        if (jj >= 10)
        {
          jj = 0;
          printf("\n");
        }
      }
    }
    printf("\n");
  }
  for (kk = 0; kk < nDesigns; kk++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = baseDesign.getEntry(kk, ii);
      finalDesign.setEntry(kk, ii, ddata);
    }
  }
  printf("genDesigns: total number of iterations = %d\n", iter);
  return 0;
}

// ************************************************************************
// improve design 
// ------------------------------------------------------------------------
int SequentialSampling::genDesigns2(psVector &params, psMatrix &testMatrix, 
                                   int nDesigns, psMatrix &finalDesign,
                                   psVector &markers)
{
  int    ii, jj, kk, mm, nn, status, nTests, markerSave, cnt;
  int    m1, m2, iter, minInd, maxInd, converged;
  double ddata, dist, meas, minVal, maxVal;
  FILE   *fp=NULL;
  psMatrix CMatrix, samMatrix, baseDesign, CIMatrix;
  psVector correlations, tmpVec;

  nTests = testMatrix.nrows();
  baseDesign.setDim(nDesigns, nInputs_);
  for (jj = 0; jj < nDesigns; jj++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = finalDesign.getEntry(jj,ii);
      baseDesign.setEntry(jj, ii, ddata);
    }
  }

  iter = 0;
  while (1)
  {
    iter++;
    samMatrix.setDim(nDesigns, nInputs_);

    converged = 1;
    minVal = PSUADE_UNDEFINED;
    for (m1 = 0; m1 < nDesigns; m1++)
    {
      cnt = 0;
      for (kk = 0; kk < nDesigns; kk++)
      {
        if (kk != m1)
        { 
          for (ii = 0; ii < nInputs_; ii++)
          {
            ddata = baseDesign.getEntry(kk,ii);
            samMatrix.setEntry(cnt, ii, ddata);
          }
          cnt++;
        }
      }

      kk = 0;
      for (kk = 0; kk < nTests; kk++) if (markers[kk] == m1) break;
      markerSave = kk;
      markers[kk] = -1;

      minInd = -1;
      for (m2 = 0; m2 < nTests; m2++)
      {
        if (markers[m2] == -1)
        {
          for (ii = 0; ii < nInputs_; ii++)
          {
            ddata = testMatrix.getEntry(m2,ii);
            samMatrix.setEntry(nDesigns-1, ii, ddata);
          }
          markers[m2] = nDesigns - 1;

          constructCMatrix(samMatrix, CMatrix, params);

          status = CMatrix.computeInverse(CIMatrix);
          if (status != 0)
          {
            printf("SequentialSampling ERROR in covariance inversion\n");
            exit(1);
          }

          correlations.setLength(nDesigns);
          maxVal = - PSUADE_UNDEFINED;
          for (nn = 0; nn < nTests; nn++)
          {
            if (markers[nn] == -1)
            {
              for (jj = 0; jj < nDesigns; jj++)
              {
                dist = 0.0;
                for (ii = 0; ii < nInputs_; ii++)
                {
                  ddata = samMatrix.getEntry(jj, ii);
                  ddata -= testMatrix.getEntry(nn,ii);
                  dist += pow(ddata, 2.0) * params[ii];
                }
                dist = exp(-dist);
                if (dist < 1.0e-50) dist = 0;
                correlations[jj] = dist;
              }
              CIMatrix.matvec(correlations, tmpVec, 0);
              meas = 0.0;
              for (jj = 0; jj < nDesigns; jj++)
                meas += correlations[jj]*tmpVec[jj];
              meas = -meas;
              if (meas > maxVal) 
              {
                maxVal = meas;
                maxInd = nn;
              }
            } // if markers[nn] == -1 
          } // nn loop

          if (maxVal == minVal)
          {
            kk = PSUADE_rand() % 2;
            if (kk == 0)
            {
              minVal = maxVal;
              minInd = m2;
            }
          }
          else if (maxVal < minVal)
          {
            minVal = maxVal;
            minInd = m2;
          }
          markers[m2] = -1;
        } // markers[m2] != -1 loop
      } 


      if (markerSave != minInd && minInd != -1)
      {
        for (ii = 0; ii < nInputs_; ii++)
        {
          ddata = testMatrix.getEntry(minInd,ii);
          baseDesign.setEntry(m1, ii, ddata);
        } 
        markers[minInd] = m1;
        markerSave  = -1;
        if (printLevel_ > 1)
          printf("genDesigns2 iter %d: swap rows %d (in) and %d (out)\n",
                 iter, minInd, kk);
        converged = 0;
        if (printLevel_ > 0)
        {
          fp = fopen("matlabsdoe.m", "a");
          fprintf(fp,"pause\n");
          fprintf(fp,"X = [\n");
          printDesigns(baseDesign, fp);
          fprintf(fp,"];\n");
          fprintf(fp,"X = X(1:%d,:);\n",nDesigns);
          fprintf(fp,
            "plot(X(:,2),X(:,3),'rh','markersize',13,'markerfacecolor','r')\n");
          fprintf(fp,"axis([%e %e %e %e])\n",
                  vecLBs_[0],vecUBs_[0],
                  vecLBs_[1],vecUBs_[1]);
          fprintf(fp,"grid on\n");
          fprintf(fp,"box on\n");
          fprintf(fp,"set(gca, 'linewidth',2)\n");
          fprintf(fp,"set(gca, 'fontsize',12)\n");
          fprintf(fp,"xlabel('X1','fontsize',12,'fontweight','bold')\n");
          fprintf(fp,"ylabel('X2','fontsize',12,'fontweight','bold')\n");
          fprintf(fp,"disp('Press enter to continue')\n");
          fclose(fp);
        }
        printf("Selected test sample points so far: (minVal=%e)\n",minVal);
        jj = 0;
        for (ii = 0; ii < nTests; ii++)
        {
          if (markers[ii] != -1)
          {
            printf("%4d ", ii+1); 
            jj++;
            if (jj >= 10)
            {
              jj = 0;
              printf("\n");
            }
          }
        }
        printf("\n");
      } // if loop
      if (markerSave != -1) markers[markerSave] = m1;
    } // m1 loop
    if (converged == 1) break;
  }
  for (kk = 0; kk < nDesigns; kk++)
  {
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = baseDesign.getEntry(kk, ii);
      finalDesign.setEntry(kk, ii, ddata);
    }
  }
  printf("genDesigns2: total number of iterations = %d\n", iter);
  return 0;
}

// ************************************************************************
// generate optimal design 
// ------------------------------------------------------------------------
int SequentialSampling::genDesignsUltimate(psVector &params, 
                               psMatrix &testMatrix, 
                               int nDesigns, psMatrix &finalDesign)
{
  int    ii, jj, kk, nn, nTests, iter;
  double ddata, minVal;
  FILE   *fp=NULL;
  psMatrix  CMatrix, samMatrix, CIMatrix;
  psIVector vecInds, vecBestInds;

  nTests = testMatrix.nrows();
  vecInds.setLength(nTests);
  vecBestInds.setLength(nTests);
  minVal = PSUADE_UNDEFINED;
  iter = 0;
  evaluateUltimate(nTests,vecInds,0,nDesigns,nDesigns,params,
                   testMatrix,vecBestInds,minVal,iter);
  printf("Final selected test sample points: \n");
  jj = 0;
  for (ii = 0; ii < nTests; ii++)
  {
    if (vecBestInds[ii] == 1)
    {
      printf("%4d ", ii+1); 
      jj++;
      if (jj >= 10)
      {
        jj = 0;
        printf("\n");
      }
    }
  }
  printf("\n");
  finalDesign.setDim(nDesigns, nInputs_);
  nn = 0;
  for (kk = 0; kk < nTests; kk++)
  {
    if (vecBestInds[kk] == 1)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = testMatrix.getEntry(kk, ii);
        finalDesign.setEntry(nn, ii, ddata);
      }
      nn++;
    }
  }
  if (printLevel_ > 0)
  {
    fp = fopen("matlabsdoe.m", "w");
    fprintf(fp,"X = [\n");
    printDesigns(finalDesign, fp);
    fprintf(fp,"];\n");
    fprintf(fp,"X = X(1:%d,:);\n",nDesigns);
    fprintf(fp,
      "plot(X(:,2),X(:,3),'rh','markersize',13,'markerfacecolor','r')\n");
    fprintf(fp,"axis([%e %e %e %e])\n",vecLBs_[0],vecUBs_[0],
            vecLBs_[1],vecUBs_[1]);
    fclose(fp);
  }
  return 0;
} 
  
// ************************************************************************
// cycle all combination 
// ------------------------------------------------------------------------
int SequentialSampling::evaluateUltimate(int nTests, psIVector &vecInds,
              int position, int nDesigns, int K, psVector &params, 
              psMatrix &testMatrix,psIVector &vecBestInds,double &minVal,
              int &iter)
{
  int    doEval=0, kk, ii, jj, nn, count, status, *localIndices;
  double ddata, dist, meas, maxVal;
  FILE   *fp;
  psVector  correlations, tmpVec;
  psIVector vecLocalInds;
  psMatrix  samMatrix, CMatrix, CIMatrix;

  fp = fopen("psuade_stop", "r");
  if (fp != NULL)
  {
    printf("PSUADE found the psuade_stop file ==> terminate.\n");
    fclose(fp);
    return 0;
  }

  vecLocalInds = vecInds;
  doEval = 0;
  if (K == 0)
  {
    doEval = 1;
  }
  if (K != 0 && (nTests-position) <= K)
  {
    for (ii = 0; ii < K; ii++) vecLocalInds[position+ii] = 1;
    doEval = 1;
  }

  if (doEval == 1)
  {
    samMatrix.setDim(nDesigns, nInputs_);
    count = 0;
    for (kk = 0; kk < nTests; kk++)
    {
      if (vecLocalInds[kk] != 0)
      { 
        for (ii = 0; ii < nInputs_; ii++)
        {
          ddata = testMatrix.getEntry(kk,ii);
          samMatrix.setEntry(count, ii, ddata);
        }
        count++;
      }
    }
    if (count != nDesigns)
    {
      printf("Catastrophic ERROR\n");
      exit(1);
    }
    constructCMatrix(samMatrix, CMatrix, params);
    status = CMatrix.computeInverse(CIMatrix);
    if (status != 0)
    {
      printf("SequentialSampling ERROR in covariance inversion\n");
      exit(1);
    }
    correlations.setLength(nDesigns);
    maxVal = - PSUADE_UNDEFINED;
    for (nn = 0; nn < nTests; nn++)
    {
      if (vecLocalInds[nn] == 0)
      {
        for (jj = 0; jj < nDesigns; jj++)
        {
          dist = 0.0;
          for (ii = 0; ii < nInputs_; ii++)
          {
            ddata = samMatrix.getEntry(jj, ii);
            ddata -= testMatrix.getEntry(nn,ii);
            dist += pow(ddata, 2.0) * params[ii];
          }
          dist = exp(-dist);
          if (dist < 1.0e-50) dist = 0;
          correlations[jj] = dist;
        }
        CIMatrix.matvec(correlations, tmpVec, 0);
        meas = 0.0;
        for (jj = 0; jj < nDesigns; jj++)
          meas += correlations[jj]*tmpVec[jj];
        meas = -meas;
        if (meas > maxVal) maxVal = meas;
      } // if markers[nn] == -1 
    } // nn loop
    fp = fopen("ps_print", "r");
    if (fp != NULL)
    {
      fclose(fp);
      jj = 0;
      for (ii = 0; ii < nTests; ii++)
      {
        if (vecLocalInds[ii] == 1)
        {
          printf("%4d ", ii); 
          jj++;
          if (jj >= 10)
          {
            jj = 0;
            printf("\n");
          }
        }
      }
      printf("\n");
    }
    if (maxVal < minVal)
    {
      minVal = maxVal;
      vecBestInds = vecLocalInds;
      printf("Iteration %d: Selected sample points (minVal = %e)\n",
             iter, minVal);
      jj = 0;
      for (ii = 0; ii < nTests; ii++)
      {
        if (vecBestInds[ii] == 1)
        {
          printf("%4d ", ii); 
          jj++;
          if (jj >= 10)
          {
            jj = 0;
            printf("\n");
          }
        }
      }
      printf("\n");
    }  
    iter++;
    return 0;
  }

  vecLocalInds[position] = 0;
  evaluateUltimate(nTests, vecLocalInds, position+1, nDesigns, K, params, 
                   testMatrix, vecBestInds, minVal, iter);
  vecLocalInds[position] = 1;
  evaluateUltimate(nTests, vecLocalInds, position+1, nDesigns, K-1, params, 
                   testMatrix, vecBestInds, minVal, iter);
  return 0;
}

// ************************************************************************
// generate an initial design 
// ------------------------------------------------------------------------
int SequentialSampling::genInitialDesigns(psVector params, 
                         psMatrix &testMatrix, int nDesigns, 
                         psMatrix &finalDesign, psVector &markers)
{
  int    ii, jj, mm, nn, ss, status, nTests, minInd;
  double ddata, dist, meas, minVal, maxVal;
  FILE   *fp;
  psMatrix CMatrix, samMatrix;
  psVector correlations, tmpVec;

  nTests = testMatrix.nrows();
  for (ss = 0; ss < nDesigns; ss++)
  {
    samMatrix.setDim(ss+1, nInputs_);
    for (jj = 0; jj < ss; jj++)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = finalDesign.getEntry(jj,ii);
        samMatrix.setEntry(jj, ii, ddata);
      }
    }

    minVal = PSUADE_UNDEFINED;
    minInd = -1;
    for (mm = 0; mm < nTests; mm++)
    {
      if (markers[mm] == -1)
      {
        for (ii = 0; ii < nInputs_; ii++)
        {
          ddata = testMatrix.getEntry(mm,ii);
          samMatrix.setEntry(ss, ii, ddata);
        }
        markers[mm] = ss;

        constructCMatrix(samMatrix, CMatrix, params);
        status = CMatrix.computeInverse(CMatrix);
        if (status != 0)
        {
          printf("SequentialSampling ERROR in covariance inversion\n");
          exit(1);
        }

        correlations.setLength(ss+1);
        maxVal = - PSUADE_UNDEFINED;
        for (nn = 0; nn < nTests; nn++)
        {
          if (markers[nn] == -1)
          {
            for (jj = 0; jj < ss+1; jj++)
            {
              dist = 0.0;
              for (ii = 0; ii < nInputs_; ii++)
              {
                ddata = samMatrix.getEntry(jj, ii);
                ddata -= testMatrix.getEntry(nn,ii);
                dist += pow(ddata, 2.0) * params[ii];
              }
              dist = exp(-dist);
              if (dist < 1.0e-50) dist = 0;
              correlations[jj] = dist;
            }
            CMatrix.matvec(correlations, tmpVec, 0);
            meas = 0.0;
            for (jj = 0; jj <= ss; jj++)
              meas += correlations[jj]*tmpVec[jj];
            meas = -meas;
            if (meas > maxVal) maxVal = meas;
          }
        }
        if (maxVal < minVal)
        {
          minVal = maxVal;
          minInd = mm;
        }
        else if (maxVal == minVal)
        {
          nn = PSUADE_rand() % 2;
          if (nn == 1)
          {
            minVal = maxVal;
            minInd = mm;
          }
        }
        markers[mm] = -1;
      }
    }
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = testMatrix.getEntry(minInd,ii);
      finalDesign.setEntry(ss, ii, ddata);
    } 
    markers[minInd] = ss;
    if (printLevel_ > 0)
    {
      fp = fopen("matlabsdoe.m", "a");
      fprintf(fp,"X = [\n");
      for (nn = 0; nn < ss+1; nn++)
      {
        fprintf(fp, "%d ", nn+1);
        for (ii = 0; ii < nInputs_; ii++)
          fprintf(fp,"%12.4e ",
            finalDesign.getEntry(nn,ii)*XStds_[ii]+XMeans_[ii]);
        fprintf(fp, "\n");
      }
      fprintf(fp,"];\n");
      fprintf(fp,
        "plot(X(:,2),X(:,3),'rh','markersize',13,'markerfacecolor','r')\n");
      fprintf(fp,"axis([%e %e %e %e])\n",vecLBs_[0],vecUBs_[0],
              vecLBs_[1],vecUBs_[1]);
      fprintf(fp,"disp('Press enter to continue')\n");
      if (ss < nDesigns-1) fprintf(fp,"pause\n");
      fclose(fp);
    }
    if (printLevel_ > 1)
      printf("genInitDesigns iteration %d: sample point added = %d\n",
             ss+1,minInd);
  }
  return 0;
}

// ************************************************************************
// generate an initial design 
// ------------------------------------------------------------------------
int SequentialSampling::genInitialDesigns2(psVector params, 
                         psMatrix &testMatrix, int nDesigns, 
                         psMatrix &finalDesign, psVector &markers)
{
  int    ii, jj, ll, mm, nn, ss, status, nTests, minInd, minInd2;
  double ddata, dist, meas, minVal, maxVal;
  FILE   *fp;
  psMatrix CMatrix, samMatrix;
  psVector correlations, tmpVec;

  nTests = testMatrix.nrows();
  for (ss = 0; ss < nDesigns; ss+=2)
  {
    samMatrix.setDim(ss+2, nInputs_);
    for (jj = 0; jj < ss; jj++)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = finalDesign.getEntry(jj,ii);
        samMatrix.setEntry(jj, ii, ddata);
      }
    }

    minVal = PSUADE_UNDEFINED;
    minInd = -1;
    for (mm = 0; mm < nTests; mm++)
    {
      for (ll = 0; ll < nTests; ll++)
      {
        if (markers[ll] == -1 && markers[mm] == -1 && (ll != mm))
        {
          for (ii = 0; ii < nInputs_; ii++)
          {
            ddata = testMatrix.getEntry(ll,ii);
            samMatrix.setEntry(ss, ii, ddata);
          }
          for (ii = 0; ii < nInputs_; ii++)
          {
            ddata = testMatrix.getEntry(mm,ii);
            samMatrix.setEntry(ss+1, ii, ddata);
          }
          markers[ll] = ss;
          markers[mm] = ss + 1;

          constructCMatrix(samMatrix, CMatrix, params);
          status = CMatrix.computeInverse(CMatrix);
          if (status != 0)
          {
            printf("SequentialSampling ERROR in covariance inversion\n");
            exit(1);
          }

          correlations.setLength(ss+2);
          maxVal = - PSUADE_UNDEFINED;
          for (nn = 0; nn < nTests; nn++)
          {
            if (markers[nn] == -1)
            {
              for (jj = 0; jj < ss+2; jj++)
              {
                dist = 0.0;
                for (ii = 0; ii < nInputs_; ii++)
                {
                  ddata = samMatrix.getEntry(jj, ii);
                  ddata -= testMatrix.getEntry(nn,ii);
                  dist += pow(ddata, 2.0) * params[ii];
                }
                dist = exp(-dist);
                if (dist < 1.0e-50) dist = 0;
                correlations[jj] = dist;
              }
              CMatrix.matvec(correlations, tmpVec, 0);
              meas = 0.0;
              for (jj = 0; jj < ss+2; jj++)
                meas += correlations[jj]*tmpVec[jj];
              meas = -meas;
              if (meas > maxVal) maxVal = meas;
            }
          }
          if (maxVal < minVal)
          {
            minVal = maxVal;
            minInd = ll;
            minInd2 = mm;
          }
          else if (maxVal == minVal)
          {
            nn = PSUADE_rand() % 2;
            if (nn == 1)
            {
              minVal = maxVal;
              minInd = ll;
              minInd2 = mm;
            }
          }
          markers[ll] = -1;
          markers[mm] = -1;
        }
      }
    }
    for (ii = 0; ii < nInputs_; ii++)
    {
      ddata = testMatrix.getEntry(minInd,ii);
      finalDesign.setEntry(ss, ii, ddata);
    } 
    markers[minInd] = ss;
    if (ss+1 < nDesigns)
    {
      for (ii = 0; ii < nInputs_; ii++)
      {
        ddata = testMatrix.getEntry(minInd2,ii);
        finalDesign.setEntry(ss+1, ii, ddata);
      }
      markers[minInd2] = ss + 1;
    }
    if (printLevel_ > 0)
    {
      fp = fopen("matlabsdoe.m", "a");
      fprintf(fp,"X = [\n");
      for (nn = 0; nn < ss+2; nn++)
      {
        if (ss+2 <= nDesigns)
        {
          fprintf(fp, "%d ", nn+1);
          for (ii = 0; ii < nInputs_; ii++)
            fprintf(fp,"%12.4e ",
              finalDesign.getEntry(nn,ii)*XStds_[ii]+XMeans_[ii]);
          fprintf(fp, "\n");
        }
      }
      fprintf(fp,"];\n");
      fprintf(fp,
        "plot(X(:,2),X(:,3),'rh','markersize',13,'markerfacecolor','r')\n");
      fprintf(fp,"axis([%e %e %e %e])\n",vecLBs_[0],vecUBs_[0],
              vecLBs_[1],vecUBs_[1]);
      fprintf(fp,"disp('Press enter to continue')\n");
      if (ss < nDesigns-1) fprintf(fp,"pause\n");
      fclose(fp);
    }
    if (printLevel_ > 1)
      printf("genInitDesigns2 iteration %d: sample points added   = %d %d\n",
             ss+1,minInd,minInd2);
  }
  return 0;
}

// ************************************************************************
// print designs
// ------------------------------------------------------------------------
int SequentialSampling::printDesigns(psMatrix &designs, FILE *fp)
{
  int ii, kk, nrows, ncols;
  if (fp != NULL)
  {
    nrows = designs.nrows();
    ncols = designs.ncols();
    for (kk = 0; kk < nrows; kk++)
    {
      fprintf(fp, "%4d ", kk+1);
      for (ii = 0; ii < ncols; ii++)
        fprintf(fp, "%12.4e ", 
                designs.getEntry(kk,ii)*XStds_[ii]+XMeans_[ii]);
      fprintf(fp, "\n");
    }
    return 0;
  }
  else return -1;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
SequentialSampling& SequentialSampling::operator=(const SequentialSampling &)
{
  printf("SequentialSampling operator= ERROR: operation not allowed.\n");
  exit(1);
  return (*this);
}

