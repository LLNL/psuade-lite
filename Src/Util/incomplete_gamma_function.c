////////////////////////////////////////////////////////////////////////////////
// File: incomplete_gamma_function.c (from mymathlib)                         //
// Routine(s):                                                                //
//    Incomplete_Gamma_Function                                               //
//    xIncomplete_Gamma_Function                                              //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The incomplete gamma function is defined as the integral from 0 to x   //
//     of the integrand t^(nu-1) exp(-t) dt.  The parameter nu is sometimes   //
//     referred to as the shape parameter.                                    //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                        // required for expl and logl().

//                         Externally Defined Routines                        //

extern long double xGamma_Function( long double x );
extern double Gamma_Function_Max_Arg( void );
extern long double xLn_Gamma_Function( long double x );
extern long double xEntire_Incomplete_Gamma_Function(long double x,
                                                                long double nu);

//                         Internally Defined Routines                        //

double Incomplete_Gamma_Function(double x, double nu);
long double xIncomplete_Gamma_Function(long double x, long double nu);

////////////////////////////////////////////////////////////////////////////////
// double Incomplete_Gamma_Function(double x, double nu)                      //
//                                                                            //
//  Description:                                                              //
//     The incomplete gamma function is defined as the integral from 0 to x   //
//     of the integrand t^(nu-1) exp(-t) dt.  The parameter nu is sometimes   //
//     referred to as the shape parameter.                                    //
//                                                                            //
//  Arguments:                                                                //
//     double x   Upper limit of the integral with integrand given above.     //
//     double nu  The shape parameter of the incomplete gamma function.       //
//                                                                            //
//  Return Values:                                                            //
//     The incomplete gamma function Integral[0,x] t^(nu-1) exp(-t) dt.       //
//                                                                            //
//  Example:                                                                  //
//     double x, g, nu;                                                       //
//                                                                            //
//     g = Incomplete_Gamma_Function( x, nu );                                //
////////////////////////////////////////////////////////////////////////////////
double Incomplete_Gamma_Function(double x, double nu) {
   return (double) xIncomplete_Gamma_Function((long double)x, (long double)nu);
}


////////////////////////////////////////////////////////////////////////////////
// long double xIncomplete_Gamma_Function(long double x, long double nu)      //
//                                                                            //
//  Description:                                                              //
//     The incomplete gamma function is defined as the integral from 0 to x   //
//     of the integrand t^(nu-1) exp(-t) dt.  The parameter nu is sometimes   //
//     referred to as the shape parameter.                                    //
//                                                                            //
//  Arguments:                                                                //
//     long double x   Upper limit of the integral with integrand given above.//
//     long double nu  The shape parameter of the incomplete gamma function.  //
//                                                                            //
//  Return Values:                                                            //
//     The incomplete gamma function Integral[0,x] t^(nu-1) exp(-t) dt.       //
//                                                                            //
//  Example:                                                                  //
//     long double x, g, nu;                                                  //
//                                                                            //
//     g = xIncomplete_Gamma_Function( x, nu );                               //
////////////////////////////////////////////////////////////////////////////////

long double xIncomplete_Gamma_Function(long double x, long double nu) {
   
   if ( x == 0.0L ) return 0.0L;
   if ( nu <= Gamma_Function_Max_Arg() )
      return xEntire_Incomplete_Gamma_Function(x,nu) * xGamma_Function(nu);
   else
      return expl(logl(xEntire_Incomplete_Gamma_Function(x,nu))
                                                  + xLn_Gamma_Function(nu));
}
