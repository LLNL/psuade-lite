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
// Functions for the class GridNode
// AUTHOR : CHARLES TONG
// DATE   : 2010
// ************************************************************************
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GridNode.h"

// ************************************************************************
// Constructor for object class SumOfTrees
// ------------------------------------------------------------------------
GridNode::GridNode()
{
   nInputs_      = 0;
   leftOrRight_  = 0;
   leftNodes_    = NULL;
   rightNodes_   = NULL;
   myCoordinate_ = NULL;
   myMultiIndex_ = NULL;
}

// ************************************************************************
// destructor
// ------------------------------------------------------------------------
GridNode::~GridNode()
{
   cleanUp();
}

// ************************************************************************
// clean up
// ------------------------------------------------------------------------
void GridNode::cleanUp()
{
   int ii;
   if (myMultiIndex_ != NULL) delete [] myMultiIndex_;
   if (myCoordinate_ != NULL) delete [] myCoordinate_;
   if (leftNodes_ != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (leftNodes_[ii] != NULL) delete leftNodes_[ii];
      delete leftNodes_;
   }
   if (rightNodes_ != NULL)
   {
      for (ii = 0; ii < nInputs_; ii++)
         if (rightNodes_[ii] != NULL) delete rightNodes_[ii];
      delete rightNodes_;
   }
   myMultiIndex_  = NULL;
   myCoordinate_ = NULL;
   leftNodes_    = NULL;
   rightNodes_   = NULL;
}

// ************************************************************************
// create tree with numLevels
// ------------------------------------------------------------------------
int GridNode::createTree(int nInputs, int *multiIndex, int numLevels, 
                         double *X, int leftOrRight, psVector &sample)
{
   int    num, ii, jj, numNodes=1, level, left=-1, right=1, found, cnt;
   double gridH;

   cleanUp();

   num = sample.length() / nInputs;
   found = 0;
   for (ii = 0; ii < num; ii++) 
   {
      for (jj = 0; jj < nInputs; jj++) 
         if (sample[ii*nInputs+jj] != X[jj]) break;
      if (jj == nInputs) found = 1;
   }
   if (found == 1) return -1;
   sample.addElements(nInputs, X);

   nInputs_      = nInputs;
   leftOrRight_  = leftOrRight;
   myMultiIndex_ = new int[nInputs];
   for (ii = 0; ii < nInputs_; ii++) myMultiIndex_[ii] = multiIndex[ii];
   myCoordinate_ = new double[nInputs];
   for (ii = 0; ii < nInputs_; ii++) myCoordinate_[ii] = X[ii];
   leftNodes_  = NULL;
   rightNodes_ = NULL;

   level = 1;
   for (ii = 0; ii < nInputs_; ii++) level += myMultiIndex_[ii] - 1;

   if (level < numLevels)
   {
      leftNodes_  = new GridNode*[nInputs]; 
      rightNodes_ = new GridNode*[nInputs]; 
      for (ii = 0; ii < nInputs_; ii++)
      {
         gridH = 1.0 / pow(2e0, myMultiIndex_[ii]);
         if ((X[ii] - gridH) >= 0.0)
         {
            leftNodes_[ii] = new GridNode();
            X[ii] -= gridH;
            multiIndex[ii]++;
            cnt = leftNodes_[ii]->createTree(nInputs,multiIndex,
                                             numLevels,X,left,sample);
            multiIndex[ii]--;
            X[ii] += gridH;
            if (cnt > 0) numNodes += cnt;
            else
            {
               delete leftNodes_[ii];
               leftNodes_[ii] = NULL;
            }
         }
         else leftNodes_[ii] = NULL;

         if ((X[ii] + gridH) <= 1.0)
         {
            rightNodes_[ii] = new GridNode();
            X[ii] += gridH;
            multiIndex[ii]++;
            cnt = rightNodes_[ii]->createTree(nInputs,multiIndex,
                                              numLevels,X,right,sample);
            multiIndex[ii]--;
            X[ii] -= gridH;
            if (cnt > 0) numNodes += cnt;
            else
            {
               delete rightNodes_[ii];
               rightNodes_[ii] = NULL;
            }
         }
         else rightNodes_[ii] = NULL;
      } 
   } 
   return numNodes;
}

// ************************************************************************
// equal operator
// ------------------------------------------------------------------------
GridNode& GridNode::operator=(const GridNode &)
{
   printf("GridNode operator= ERROR: operation not allowed.\n");
   exit(1);
   return (*this);
}

