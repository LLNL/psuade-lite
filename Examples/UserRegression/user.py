#!/usr/bin/python
import sys
import string
infile  = open(sys.argv[1], "r")
lineIn  = infile.readline()
ncols   = string.split(lineIn)
nSamp   = eval(ncols[0])
nInps   = eval(ncols[1])
outfile = open(sys.argv[2], "w")
for ii in range(nSamp):
  lineIn  = infile.readline()
  ncols   = lineIn.split()
  X1      = eval(ncols[1])
  X2      = eval(ncols[2])
  Y = X1 * X1
  outfile.write("%d " % (ii+1))
  outfile.write("%e " % Y)
  Y = X2 * X2
  outfile.write("%e " % Y)
  Y = X1 * X2
  outfile.write("%e\n" % Y)
infile.close()
outfile.close()

