/*  

  If more than BIGWORK comparisons are required in
an oacheck routine, then a warning is printed that
a large job is underway.  If more than MEDWORK comparisons
are required then intermediate results are printed.

  No strength checking beyond strength MAXSTR is done.
Only change it if you implement the higher strength
checks!

*/

#define BIGWORK 100000000
#define MEDWORK BIGWORK

