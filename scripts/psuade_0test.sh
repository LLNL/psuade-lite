#!/bin/sh
# This is a test script for PSUADE.  It runs PSADE once,
# just an input file (usually for smapling)
# The terminal output of the run is verified against
# the saved output. (First 15 lines are skipped to skip version #s)
#$1 <path to psuade>     
#$2 <input psuade file>
#$3 <The offical saved output for comparison> (psScript.out)
$1 $2 > ./tmp.out
tail -n+15 tmp.out > test.out
tail -n+15 $3 > reference.out
diff ./test.out ./reference.out > `basename $3`.diff
