/*
------------------------------------------------------------------------------
Support for R functions
------------------------------------------------------------------------------
*/
#include <math.h>
#include <stdio.h>

#ifdef NAN
/*
static double R_NaN=NAN;
*/
static double R_NaN=0.0/0.0;
#else
static double R_NaN=0.0/0.0;
#endif /* NAN*/

#ifdef HUGE
static double R_NegInf=-HUGE_VAL;
static double R_PostInf=HUGE_VAL;
#else 
static double R_NegInf=((-1.0) / 0.0);
static double R_PostInf=(1.0 / 0.0);
#endif


void Rf_warning(int id, char *str)
{
   printf("R warning: code = %d, message = %s\n",id,str);
}

