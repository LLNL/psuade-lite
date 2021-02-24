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
// Definition for the class PsuadeData
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __PSUADEDATAH__
#define __PSUADEDATAH__

#include <math.h>
//#include <stream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pData.h"
#include "Matrix.h"
#include "PsuadeSession.h"

// ************************************************************************
// subclasses
// ************************************************************************

class psuadeInputSection 
{
public:
   int    nFixed_;
   double *fixedValues_;
   char   **fixedNames_;
   int    nInputs_;
   double *inputLBounds_;
   double *inputUBounds_;
   char   **inputNames_;
   int    *inputPDFs_;
   double *inputMeans_;
   double *inputStdevs_;
   double *inputAuxs_;
   int    *symbolTable_;
   double **inputSettings_;
   int    *inputNumSettings_;
   double *sampleInputs_;
   psMatrix corMatrix_;
   int    useInputPDFs_;
   int    *inputSIndices_;
   char   **sampleFileNames_;

   psuadeInputSection();
   ~psuadeInputSection();
   void reset();
};

class psuadeOutputSection 
{
public:
   int    nOutputs_;
   char   **outputNames_;
   double *sampleOutputs_;
   int    *sampleStates_;

   psuadeOutputSection();
   ~psuadeOutputSection();
   void reset();
};

class psuadeMethodSection 
{
public:
   int    samplingMethod_;
   int    nSamples_;
   int    nReplications_;
   int    nRefinements_;
   int    refNumRefinements_;
   int    refinementType_;
   int    refinementSize_;
   int    sampleRandomize_;

   psuadeMethodSection();
   ~psuadeMethodSection();
   void reset();
};

class psuadeApplicationSection
{
public:
   char   appDriver_[2000];
   char   rsDriver_[2000];
   char   optDriver_[2000];
   char   auxOptDriver_[2000];
   char   ensembleDriver_[2000];
   char   ensembleOptDriver_[2000];
   char   inputTemplate_[2000];
   char   outputTemplate_[2000];
   int    maxParallelJobs_;
   int    maxJobWaitTime_;
   int    minJobWaitTime_;
   int    runType_;
   int    useFunction_;
   int    launchInterval_;
   int    saveFrequency_;

   psuadeApplicationSection();
   ~psuadeApplicationSection();
   void reset();
};

typedef struct 
{
   double FilterLBound_;
   double FilterUBound_;
   char   *FilterDataFile_;  /* psuade data file name as RS filter    */
   char   *FilterIndexFile_; /* index file that gives input indices   */
} psuadeFilter;

// ------------------------------------------------------------------------

class psuadeAnalysisSection 
{
public:

   //  analysisIntOptions[0] = analysis method
   //  analysisIntOptions[1] = output ID
   //  analysisIntOptions[2] = rstype
   //  analysisIntOptions[3] = diagnostics
   //  analysisIntOptions[4] = graphics
   int analysisIntOptions_[10];
   
   //  analysisDbleOptions[0] = threshold
   double analysisDbleOptions_[10];

   //  optimizeIntOptions[0] = turn on/off optimization
   //  optimizeIntOptions[1] = choose optimization method
   //  optimizeIntOptions[2] = store num_local_minima
   //  optimizeIntOptions[3] = use response surface for optimization
   //  optimizeIntOptions[4] = output level
   //  optimizeIntOptions[5] = num of fmin to calculate
   //  optimizeIntOptions[6] = optimization output ID
   int    optimizeIntOptions_[10];
  
   //  optimizeDbleOptions[0] = fmin
   //  optimizeDbleOptions[1] = cutoff point for optimization
   //  optimizeDbleOptions[2] = tolerance for optimization
   double optimizeDbleOptions_[10];
   int    fileWriteFlag_;
   char   specsFile_[2000];
   int    numRSFilters_;
   int    numMOATFilters_;
   char   rsIndexFile_[2000];
   char   rsIndexSampleFile_[2000];
   int    useInputPDFs_;
   int    legendreOrder_;
   int    marsNbasis_;
   int    marsNdegrees_;
   int    marsNum_;
   int    kriMode_;
   double kriTol_;
   int    optSaveHistory_;
   int    optUseHistory_;
   psuadeFilter **RSFilters_;
   psuadeFilter **MOATFilters_;

   psuadeAnalysisSection();
   ~psuadeAnalysisSection();
   void reset();
};

// ************************************************************************
// ************************************************************************
// main class definition 
// ************************************************************************
// ************************************************************************
class PsuadeData 
{

   char                     psuadeFileName_[200];
   int                      outputLevel_;
   int                      writeCnt_;
   psuadeInputSection       pInput_;
   psuadeOutputSection      pOutput_;
   psuadeMethodSection      pMethod_;
   psuadeApplicationSection pApplication_;
   psuadeAnalysisSection    pAnalysis_;
   pData                    auxData_;

public:

   // Constructor
   PsuadeData();

   // Destructor
   ~PsuadeData();

   // Read input file
   // param fname - name of the file to read
   int readPsuadeFile(const char *fname);

   // Write to output file
   // param fname - name of the file to write
   // param flag  - to indicate duplicate save or not
   void writePsuadeFile(const char *fname, int);

   // Get parameter 
   int getParameter(const char *, pData &);

   // create the input section (no samples)
   //param nInputs  - number of input variables (for checking only)
   //param nSymbols - number of symbols for the inputs
   //param lowerB - lower bounds
   //param upperb - upper bounds
   //param names    - input names
   void createInputSection(int nInputs, int *nSymbols,
                           double *lowerB, double *upperB,
                           char **names);

   // Update information about the input section
   // param nSamples - number of sample points
   // param nInputs  - number of input variables (for checking only)
   // param nSymbols - number of symbols for the inputs
   // param lowerB - lower bounds
   // param upperb - upper bounds
   // param sampleInputs - sample input values
   // param names    - input names
   // Note : set to NULL or -1 if update not requested.
   void updateInputSection(int nSamples, int nInputs, int *nSymbols,
                           double *lowerB, double *upperB,
                           double *sampleInputs, char **names,
                           int *, double *, double *, psMatrix *);

   // Basic creation of the output section (no samples)
   // param nOutputs - number of output variables (for checking only)
   // param names         - output names
   void createOutputSection(int nOutputs, char **names);
   
   // Update information about the output section
   // param nSamples - number of sample points
   // param nOutputs - number of output variables (for checking only)
   // param sampleOutputs - sample output values
   // param sampleStates  - execution states of the samples
   // param names         - output names
   void updateOutputSection(int nSamples, int nOutputs,
                            double *sampleOutputs, int *sampleStates,
                            char **names);


   // Clears out any existing samples so they can be replaced.
   void resetSamples();

   // Update information about the method section
   // param method - sampling method
   // param nSamples - number of samples
   // param nReps - number of replications
   // param nRefines - number of refinements
   // param randomize - randomize flag
   void updateMethodSection(int method, int nSamples, int nReps,
                            int nRefines, int randomize);

   // Update information about the application section
   // @param appDriver - appDriver name
   // @param optDriver - optDriver name
   // @param optAuxDriver - optAuxDriver name
   // @param ensembleDriver - ensembleDriver name
   // @param ensembleOptDriver - ensembleOptDriver name
   // @param maxJobs - maximum number of parallel jobs
   void updateApplicationSection(char *appDriver, char *optDriver,
                  char *auxOptDriver, char *ensembleDriver,
                  char *ensembleOptDriver, int maxJobs);

   // @param methods - analysis method bitmask
   // @param transform - input/output transform
   // @param rstype - response surface type
   // @param diag - diagnostics level 
   // @param outputID - output identifier (base 1)
   // @param usePDFs - flag to indicate PDFs are used
   void updateAnalysisSection(int methods, int transform, int rstype,
			      int diag, int outputID, int usePDFs);

   // Update information about the optimization section
   // @param methods   - opt method 
   // @param nLocalMin - number of local minima
   // @param tolerance - opt tolerance 
   void updateOptimizationSection(int methods, int nLocalMin,
                                  double tolerance, int maxIter);

   // get sample information 
   // @param session   - place holder
   int getSession(PsuadeSession *session);

   // special processing of output (e.g. write to file in matlab format)
   void processOutputData();

   // set output level for error outputs
   void setOutputLevel(int);

   // read in sample data from the PSUADE_IO section of a psuadeData file
   int readPsuadeIO(const char *);

   // assignment operator override
   PsuadeData& operator=(const PsuadeData &);

   // accessor functions for SWIG to read C array data
   int    getSampleState(int sampleid);
   double getSampleInput(int inputNumber, int sampleid);
   double getSampleOutput(int outputNumber, int sampleid);
   pData  *getAuxData() {return &auxData_;}

   // accessor functions for everything.  Arrays are copied before return,
   // so the caller must call delete[] on them.
   int     getInput_nInputs();
   double* getInput_inputLBounds();
   double* getInput_inputUBounds();
   char**  getInput_inputNames();
   int*    getInput_inputPDFs();
   double* getInput_inputMeans();
   double* getInput_inputStdevs();
   double* getInput_inputAuxs();
   int*    getInput_symbolTable();
   //   double* getInput_inputSettings();  //Use getParameter("input_settings""
   //int*    getInput_inputNumSettings();  //to get settings and numsettings
   double* getInput_sampleInputs();
   psMatrix  getInput_corMatrix();
   int     getInput_useInputPDFs();

   //Output getters
   int getOutput_nOutputs();
   char** getOutput_outputNames();
   double* getOutput_sampleOutputs();
   int* getOutput_sampleStates();

   //method getters
   int getMethod_samplingMethod();
   int getMethod_nSamples();
   int getMethod_nReplications();
   int getMethod_nRefinements();
   int getMethod_refNumRefinements();
   int getMethod_refinementType();
   int getMethod_refinementSize();
   int getMethod_sampleRandomize();
   
   //application section getters
   char* getApplication_appDriver();
   char* getApplication_rsDriver();
   char* getApplication_optDriver();
   char* getApplication_auxOptDriver();
   char* getApplication_ensembleDriver();
   char* getApplication_ensembleOptDriver();
   char* getApplication_inputTemplate();
   char* getApplication_outputTemplate();
   int getApplication_maxParallelJobs();
   int getApplication_maxJobWaitTime();
   int getApplication_minJobWaitTime();
   int getApplication_runType();
   int getApplication_useFunction();
   int getApplication_launchInterval();
   int getApplication_saveFrequency();

   //analysis section getters
   int getAnalysis_fileWriteFlag();
   int getAnalysis_useInputPDFs();
   int getAnalysis_legendreOrder();
   int getAnalysis_marsNbasis();
   int getAnalysis_marsNdegrees();
   int getAnalysis_marsNum();
   int getAnalysis_kriMode();
   double getAnalysis_kriTol();
   int getAnalysis_optSaveHistory();
   int getAnalysis_optUseHistory();

   //Returns an array of length 10 with the following values:
   //  analysisIntOptions[0] = analysis method
   //  analysisIntOptions[1] = output ID
   //  analysisIntOptions[2] = rstype
   //  analysisIntOptions[3] = diagnostics
   //  analysisIntOptions[4] = graphics
   int* getAnalysis_analysisIntOptions();

   //Returns an array of length 10 with the following values:
   //  analysisDbleOptions[0] = threshold
   double* getAnalysis_analysisDbleOptions();

   //Returns an array of length 10 with the following values:
   //  optimizeIntOptions[0] = turn on/off optimization
   //  optimizeIntOptions[1] = choose optimization method
   //  optimizeIntOptions[2] = store num_local_minima
   //  optimizeIntOptions[3] = use response surface for optimization
   //  optimizeIntOptions[4] = output level
   //  optimizeIntOptions[5] = num of fmin to calculate
   //  optimizeIntOptions[6] = optimization output ID
   int* getAnalysis_optimizeIntOptions();

   //Returns an array of length 10 with the following values:
   //  optimizeDbleOptions[0] = fmin
   //  optimizeDbleOptions[1] = cutoff point for optimization
   //  optimizeDbleOptions[2] = tolerance for optimization
   double* getAnalysis_optimizeDbleOptions();
   
   char* getAnalysis_specsFile();
   char* getAnalysis_rsIndexFile();
   char* getAnalysis_rsIndexSampleFile();
   
   // Use the GetParameter interface to get information on RSFilters 
   // and MOATFilters
   //   int getAnalysis_numRSFilters();
   //int getAnalysis_numMOATFilters();
   //psuadeFilter **RSFilters_;
   //psuadeFilter **MOATFilters_;

private:
    int readInputSection(FILE *);
    int readOutputSection(FILE *);
    int readMethodSection(FILE *);
    int readApplicationSection(FILE *);
    int readAnalysisSection(FILE *);
    int getInputParameter(const char *, pData &);
    int getOutputParameter(const char *, pData &);
    int getMethodParameter(const char *, pData &);
    int getApplicationParameter(const char *, pData &);
    int getAnalysisParameter(const char *, pData &);

public:
    void writePsuadeIO(FILE *, int);
    void writeInputSection(FILE *);
    void writeOutputSection(FILE *);
    void writeMethodSection(FILE *);
    void writeSamplingSection(FILE *);
    void writeApplicationSection(FILE *);
    void writeAnalysisSection(FILE *);
};

#endif // __PSUADEDATAH__

