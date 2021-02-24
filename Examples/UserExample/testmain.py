#!/usr/bin/python
import string
import sys
infile  = open(sys.argv[1], "r")
lineIn  = infile.readline()
ncols   = string.split(lineIn)
n       = eval(ncols[0])
X       = range(n)
for ii in range(n):
   lineIn  = infile.readline()
   ncols   = string.split(lineIn)
   X[ii]   = eval(ncols[0])
infile.close()
Y = 0
for ii in range(n-1):
   Y=Y+ pow(1-X[ii],2)+100*pow(X[ii+1]-X[ii]*X[ii],2)

outfile = open(sys.argv[2], "w")
outfile.write("%e \n" % Y)
outfile.close()

