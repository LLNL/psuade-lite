#!/bin/sh
# This is a test script for PSUADE.  It runs PSADE once,
# just an input file (usually for smapling)
# The terminal output of the run is verified against
# the saved output. (First 15 lines are skipped to skip version #s)
#$1 <path to psuade>     
#$2 <input psuade file>
#$3 <The offical saved output for comparison> (psScript.out)
PSUADE="$1"  #This is a trick to allow the user of the script to auto delete files by putting them after the important args
INFILE="$2"
REFOUT="$3"
shift
shift
shift
rm -rf "$@"
$PSUADE < $INFILE > ./tmp.out
tail -n+15 tmp.out > test.out
tail -n+15 $REFOUT > reference.out
diff ./test.out ./reference.out > `basename $3`.diff
