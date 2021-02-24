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
// internal identifiers
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************

#define PSUADE_NO_NOISE          0
#define PSUADE_ADD_NOISE         1
#define PSUADE_RAN_NOISE         2

#define PSUADE_UNDEFINED         1.0E35

#define PSUADE_SAMP_MC           0
#define PSUADE_SAMP_FACT         1
#define PSUADE_SAMP_LHS          2
#define PSUADE_SAMP_OA           3
#define PSUADE_SAMP_OALH         4
#define PSUADE_SAMP_MOAT         5
#define PSUADE_SAMP_SOBOL        6
#define PSUADE_SAMP_LPTAU        7
#define PSUADE_SAMP_METIS        8
#define PSUADE_SAMP_FAST         9
#define PSUADE_SAMP_BBD         10
#define PSUADE_SAMP_PBD         11
#define PSUADE_SAMP_FF4         12
#define PSUADE_SAMP_FF5         13
#define PSUADE_SAMP_CCI4        14
#define PSUADE_SAMP_CCI5        15
#define PSUADE_SAMP_CCIF        16
#define PSUADE_SAMP_CCF4        17
#define PSUADE_SAMP_CCF5        18
#define PSUADE_SAMP_CCFF        19
#define PSUADE_SAMP_CCC4        20
#define PSUADE_SAMP_CCC5        21
#define PSUADE_SAMP_CCCF        22
#define PSUADE_SAMP_SFAST       23
#define PSUADE_SAMP_UMETIS      24
#define PSUADE_SAMP_GMOAT       25
#define PSUADE_SAMP_GMETIS      26
#define PSUADE_SAMP_SG          27
#define PSUADE_SAMP_DISCRETE    28
#define PSUADE_SAMP_LSA         29
#define PSUADE_SAMP_RFF4        30
#define PSUADE_SAMP_RFF5        31
#define PSUADE_SAMP_SEQ         32

#define PSUADE_ANA_MOMENT       1
#define PSUADE_ANA_ME           1<<1
#define PSUADE_ANA_IE           1<<2
#define PSUADE_ANA_ANOVA        1<<3
#define PSUADE_ANA_GLSA         1<<4
#define PSUADE_ANA_RSFA         1<<5
#define PSUADE_ANA_MOAT         1<<6
#define PSUADE_ANA_SOBOL        1<<7
#define PSUADE_ANA_CORRELATION  1<<8
#define PSUADE_ANA_INTEGRATION  1<<9
#define PSUADE_ANA_FAST         1<<10
#define PSUADE_ANA_FF           1<<11
#define PSUADE_ANA_PCA          1<<12
#define PSUADE_ANA_ONESIGMA     1<<13
#define PSUADE_ANA_FORM         1<<14
#define PSUADE_ANA_RSSOBOL1     1<<15
#define PSUADE_ANA_RSSOBOL2     1<<16
#define PSUADE_ANA_RSSOBOLTSI   1<<17
#define PSUADE_ANA_BSTRAP       1<<18
#define PSUADE_ANA_RSSOBOLG     1<<19
#define PSUADE_ANA_1SAMPLE      1<<20
#define PSUADE_ANA_2SAMPLE      1<<21
#define PSUADE_ANA_MCMC         1<<22
#define PSUADE_ANA_ARSMNN       1<<23
#define PSUADE_ANA_REL          1<<24
#define PSUADE_ANA_AOPT         1<<25
#define PSUADE_ANA_DTEST        1<<26
#define PSUADE_ANA_ARSMMB       1<<27
#define PSUADE_ANA_GOWER        1<<28
#define PSUADE_ANA_ETEST        1<<29
#define PSUADE_ANA_ARSMMBBS     1<<30
#define PSUADE_ANA_LSA          1<<31

/* do not change these settings, as they
   have been hardwired. 
   both MGP1 and PWL are not used
*/
#define PSUADE_NUM_RS           29
#define PSUADE_RS_MARS          0
#define PSUADE_RS_REGR1         1
#define PSUADE_RS_REGR2         2
#define PSUADE_RS_REGR3         3
#define PSUADE_RS_REGR4         4
#define PSUADE_RS_ANN           5
#define PSUADE_RS_REGRS         6
#define PSUADE_RS_GP1           7
#define PSUADE_RS_GP2           8
#define PSUADE_RS_SVM           9
#define PSUADE_RS_REGRGL        10
#define PSUADE_RS_TGP           11
#define PSUADE_RS_MARSB         12
#define PSUADE_RS_EARTH         13
#define PSUADE_RS_SOTS          14
#define PSUADE_RS_REGRL         15
#define PSUADE_RS_REGRU         16
#define PSUADE_RS_REGSG         17
#define PSUADE_RS_KR            18
#define PSUADE_RS_SPLINES       19
#define PSUADE_RS_KNN           20
#define PSUADE_RS_RBF           21
#define PSUADE_RS_ACOSSO        22
#define PSUADE_RS_BSSANOVA      23
#define PSUADE_RS_RBFB          24
#define PSUADE_RS_PLS           25
#define PSUADE_RS_MRBF          26
#define PSUADE_RS_MGP2          27
#define PSUADE_RS_MMARS         28
#define PSUADE_RS_LOCAL         29
#define PSUADE_RS_NPL           30
#define PSUADE_RS_MGP1          31
#define PSUADE_RS_PWL           32
#define PSUADE_RS_REGRLC        33
#define PSUADE_RS_CSREG         34

#define PSUADE_GRAPHICS         1
#define PSUADE_SAMPLE_GRAPHICS  2

#define PSUADE_PDF_UNIFORM      0
#define PSUADE_PDF_NORMAL       1
#define PSUADE_PDF_LOGNORMAL    2
#define PSUADE_PDF_TRIANGLE     3
#define PSUADE_PDF_BETA         4
#define PSUADE_PDF_WEIBULL      5
#define PSUADE_PDF_GAMMA        6
#define PSUADE_PDF_EXPONENTIAL  7
#define PSUADE_PDF_SAMPLE       8
#define PSUADE_PDF_F            9
#define PSUADE_PDF_SAMPLEHIST  10
#define PSUADE_PDF_INVGAMMA    11
#define PSUADE_PDF_CAUCHY      12
#define PSUADE_PDF_USER        13

//define parameters associated with printLevel & diagnostics
//added 4/2/14 - Jim McEnerney
#define PL_MIN -2
#define PL_MAX 5
#define PL_ERROR -2
#define PL_WARN -1
#define PL_INFO 0
#define PL_BASIC 1
#define PL_INTERACTIVE 2
#define PL_MOREINFO 3
#define PL_DETAIL 4
#define PL_DUMP 5

