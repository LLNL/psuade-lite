// ************************************************************************
// Util/Exceptions.h
//
// Exceptions which can be raised during psuade's execution
// ************************************************************************

#ifndef __EXCEPTIONSH__
#define __EXCEPTIONSH__


// ************************************************************************
// Base class for all psuade exceptions. Catch this exception to cover them all
// ************************************************************************

class Psuade_Exception {};

// ************************************************************************
// Exception for halting PsuadeBase:runDace*
//   when the file psuade_stop is detected, or when
//   Python requests a stop by setting PsuadeBase::psuade_stop = 1
// ************************************************************************

class Psuade_Stop_Exception : Psuade_Exception {};


#endif // __EXCEPTIONSH__

