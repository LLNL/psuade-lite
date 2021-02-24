#!/bin/sh
# This is a test script for PSUADE.  It runs PSADE twice,
# 1st time with just an input file (usually for smapling)
# 2nd time by running a script, usually for analysis.
# The terminal output of the second run is verified against
# the saved output. (First 15 lines are skipped to skip version #s)
#$1 <path to psuade>     
#$2 <input psuade file>
#$3 <Name of the script file>
#$4 <The offical saved output for comparison> (psScript.out)
echo $1 $ $3 $4
$1 $2 > ./test.log
mv psuadeData psSave
$1 < $3 > ./tmp.out
tail -n+15 tmp.out > test.out
tail -n+15 $4 > reference.out
diff ./test.out ./reference.out > `basename $4`.diff
