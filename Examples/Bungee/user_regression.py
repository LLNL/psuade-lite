#!/usr/bin/python
import sys
import string
infile  = open(sys.argv[1], "r")
lineIn  = infile.readline()
ncols   = lineIn.split()
nSamp   = eval(ncols[0])
nInps   = eval(ncols[1])
outfile = open(sys.argv[2], "w")
for ii in range(nSamp):
   lineIn  = infile.readline()
   ncols   = lineIn.split()
   id = eval(ncols[0])
   if (id != ii+1):
      write("error")
   X1 = eval(ncols[1])
   X2 = eval(ncols[2])
   X3 = eval(ncols[3])
   outfile.write("%d " % (ii+1))
   outfile.write("1 ")
   outfile.write("%e " % X1)
   outfile.write("%e " % X2)
   outfile.write("%e " % X3)
   Y =  X1 * X1
   outfile.write("%e " % Y)
   Y =  X1 * X2
   outfile.write("%e " % Y)
   Y =  X1 * X3
   outfile.write("%e " % Y)
   Y =  X2 * X2
   outfile.write("%e " % Y)
   Y =  X2 * X3
   outfile.write("%e " % Y)
   Y =  X3 * X3
   outfile.write("%e \n" % Y)
infile.close()
outfile.close()

