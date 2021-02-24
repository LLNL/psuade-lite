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
// interface functions to MARS (fortran)
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
*/
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_MARS

#if !defined(SUNOS)
extern void mars_(int *,int *, float *, float *, float *, int *,
	      int *, int *, float *, int *, float *, double*, int*);
extern void fmod_(int *, int *, float *, float *, int *, float *,
                  float *);
#endif

/*************************************************************************/
/* call MARS to generate regression coefficients                         */
/*-----------------------------------------------------------------------*/

#if defined(SUNOS)

void mars_process(nsamples, ninputs, x, y, wgts, nk, mi, lx, fm, im)
int    nsamples, ninputs, nk, mi, *lx, *im;
double **x, *y;
float  *wgts, *fm;

#else
void mars_process(int nsamples, int ninputs, double **x, double *y, 
               float *wgts, int nk, int mi, int *lx, float *fm, int *im)
#endif
{
   int     sampleID, inputID, length, *iws;
   float   *fws, *xlocal, *ylocal;
   double  *dws;

   if (nsamples <= 0 || ninputs <= 0)
   {
      printf("mars_process ERROR: invalid parameters.\n");
      exit(1);
   }
   length = nsamples * (nk + 4) + 6 * nsamples + 2 * ninputs + 4 * nk; 
   fws = (float*)  malloc(length * 2 * sizeof(float));
   length = 2 * ((nsamples + 1) * (nk + 1) + 10 * nk);
   length = (nsamples + 1) * (nk + 1) + (nk + 2) * (nsamples + 3) + 4 * nk;
   dws = (double*) malloc(length * sizeof(double));
   length = nsamples * ninputs + 2 * nsamples;
   iws  = (int*)   malloc(length * sizeof(int));

   xlocal = (float*) malloc(nsamples * ninputs * sizeof(float));
   ylocal = (float*) malloc(nsamples * sizeof(float));
   if (ylocal == NULL)
   {
     printf("Memory allocation ERROR in: ylocal in MarsInterface\n");
     exit(1);
   }

   for (sampleID = 0; sampleID < nsamples; sampleID++) 
   {
      ylocal[sampleID] = y[sampleID];
      for (inputID = 0; inputID < ninputs; inputID++)
         xlocal[inputID*nsamples+sampleID] = x[sampleID][inputID];
   }
   mars_(&nsamples,&ninputs,xlocal,ylocal,wgts,&nk,&mi,lx,fm,im,fws,dws,iws);

   free(fws);
   free(dws);
   free(iws);
   free(xlocal);
   free(ylocal);
}

/*************************************************************************/
/* call MARS to interpolate input data                                   */
/*-----------------------------------------------------------------------*/

#if defined(SUNOS)
void mars_fmod(tot_npts, ninputs, x, y, fm, im)
int    tot_npts, ninputs, *im;
double **x, *y;
float  *fm;
#else
void mars_fmod(int tot_npts, int ninputs, double **x, double *y,
               float *fm, int *im)
#endif
{
   int    sampleID, inputID, model_flag;
   float  *sp, *xlocal, *ylocal;

   model_flag = 1; /* 1: piecewise-linear, 2: cubic */

   if (tot_npts <= 0)
   {
      printf("mars_fmod ERROR: invalid parameters (nPts=%d).\n",tot_npts);
      exit(1);
   }
   sp = (float*) malloc(tot_npts * 4 * sizeof(float));
   xlocal = (float*) malloc(tot_npts * ninputs * sizeof(float));
   ylocal = (float*) malloc(tot_npts * sizeof(float));

   for (sampleID = 0; sampleID < tot_npts; sampleID++) 
      for (inputID = 0; inputID < ninputs; inputID++)
         xlocal[inputID*tot_npts+sampleID] = x[sampleID][inputID];

   fmod_(&model_flag, &tot_npts, xlocal, fm, im, ylocal, sp);

   for (sampleID = 0; sampleID < tot_npts; sampleID++)
      y[sampleID] = ylocal[sampleID];

   free(sp);
   free(xlocal);
   free(ylocal);
}

#endif

