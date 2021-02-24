#!/usr/bin/python
import os
import sys
import string
import shutil

#######################################################
# Function to get input data from PSUADE-generated
# parameter files (standard format, do not change). 
# The return data will contain the inputs values.
#======================================================

def getInputData(inFileName):
   inFile  = open(inFileName, "r")
   lineIn  = inFile.readline()
   nCols   = string.split(lineIn)
   nInputs = eval(nCols[0])
   inputData = range(nInputs)
   for ind in range(nInputs):
      lineIn  = inFile.readline()
      nCols   = string.split(lineIn)
      inputData[ind] = eval(nCols[0])
   inFile.close()
   return inputData

#######################################################
# Function to generate output file 
# This function writes the output data (which should
# have been generated in outData) to the PSUADE-based
# output file.
#======================================================

def genOutputFile(outFileName, outData):
   nLeng = len(outData)
   outfile = open(outFileName, "w")
   for ind in range(nLeng):
      outfile.write("%e \n" % outData[ind])
   outfile.close()
   return

#######################################################
# Main program file 
#######################################################

inputfileName  = sys.argv[1]
outputfileName = sys.argv[2]

inputData = getInputData(inputfileName)
F1 = inputData[0] * inputData[0] + inputData[1] * inputData[1] 
F2 = 1 - F1
outdata = range(2)
outdata[0] = F1
outdata[1] = F2
genOutputFile(outputfileName, outdata)

