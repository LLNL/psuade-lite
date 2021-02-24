/*===========================================================================*/
/**
 *  \file PrintingTS.c
 *
 *  \brief Utility routines for printing output to stdout, stderr, and log
 *  files.
 */
/*  Copyright:  (c) 2004-2007 The Regents of the University of California.
 *
 *===========================================================================*/

#include "PrintingTS.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/time.h>

#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>

/** \brief The overall print level we're running at */
/*------------------------------------*/
static int s_printLevel = 1;

/** \brief Do thread safe printing? */
/*------------------------------------*/
static int doPrintTS = 1;

/** \brief Do thread safe Logging printing? */
/*------------------------------------------*/
static int doLogTS = 0;

/** \brief MPI rank of symponent process. */
/*------------------------------------------*/
static int localRank = -1;

/** \brief Log file pointer. */
/*-----------------------------*/
static FILE* fpLog;

/** \brief Lock that serializes output. */
/*----------------------------------------*/
//static pthread_mutex_t printLock = PTHREAD_MUTEX_INITIALIZER;

/*---------------------------------------------------------------------------
 *  initializePrintingTS
 */
/** \brief Create directory to hold log files if needed; open log file.
 *
 *  This routine recursively removes a logging directory if it already exists
 *  and then creates a new empty top level logging directory with the same
 *  name.  It also constructs the log file name and opens it.
 *
 *  doPrintTS gets set to 1/True.
 *
 *  \param[in] printLevel  : The level of output to display and log.
 *  \param[in] logFileName : Log file name.
 *  \param[in] myRank      : Rank that is given.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
initializePrintingTS(int printLevel, const char* logFileName, int myRank)
{
  char commandString[128];
  char logFileNameFull[256];

//  pthread_mutex_lock(&printLock);

  localRank = myRank;

  s_printLevel = printLevel;

  if(logFileName) {
    sprintf(logFileNameFull, "%s.%05d", logFileName, myRank);
    fpLog = fopen(logFileNameFull, "w");
  }

  if (fpLog != NULL) { doLogTS = 1; }

//  pthread_mutex_unlock(&printLock);
}

/*---------------------------------------------------------------------------
 *  finalizePrintingTS
 */
/** \brief Closes the log file.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
finalizePrintingTS(void)
{
//  pthread_mutex_lock(&printLock);

  doPrintTS = 0;
  
  if (fpLog != NULL) { fclose(fpLog); fpLog = NULL; }

//  pthread_mutex_unlock(&printLock);
}

/*---------------------------------------------------------------------------
 *  printErrTS
 */
/** \brief Prints error message to log file & stderr.
 *
 *  \param[in] format : Text to print.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
printErrTS(int printLevel, const char* format, ...)
{
  char buffer[4096];
  va_list args;


//  pthread_mutex_lock(&printLock);

  if(printLevel <= s_printLevel) {
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    
    if(doLogTS) {
      fputs(buffer, fpLog);
      fflush(fpLog);
    }
    
    if (doPrintTS) {    
      fputs(buffer, stderr);
      fflush(stderr);
    }    

    va_end(args);
  }

//  pthread_mutex_unlock(&printLock);

}

/*---------------------------------------------------------------------------
 *  printLogTS
 */
/** \brief Prints log message to log file.
 *
 *  \param[in] format : Text to print.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
printLogTS(int printLevel, const char* format, ...)
{
  char buffer[4096];
  va_list args;

//  pthread_mutex_lock(&printLock);

  if (printLevel <= s_printLevel && doLogTS)
  {
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);

    fputs(buffer, fpLog);
    fflush(fpLog);

    va_end(args);
  }
 
//  pthread_mutex_unlock(&printLock);
}

/*---------------------------------------------------------------------------
 *  printOutTS
 */
/** \brief Prints message to log file & stdout.
 *
 *  \param[in] format : Text to print.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
printOutTS(int printLevel, const char* format, ...)
{
  char buffer[4096];
  va_list args;

//  pthread_mutex_lock(&printLock);

  if (printLevel <= s_printLevel) {
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);

    if(doLogTS) {
      fputs(buffer, fpLog);
      fflush(fpLog);
    }

    if (doPrintTS) {
      fputs(buffer, stdout);
      fflush(stdout);
    }
    
    va_end(args);

  }
//  pthread_mutex_unlock(&printLock);
}

/*---------------------------------------------------------------------------
 *  isPrintTSOn
 */
/** \brief Is thread safe printing on?
 *
 *  \return On/Off indicator. */
/*---------------------------------------------------------------------------*/
int
isPrintTSOn(void)
{

  /* I think we can get away without locking here. */

  return doPrintTS;
}

/*---------------------------------------------------------------------------
 *  turnPrintTSOff
 */
/** \brief Turn thread safe printing off.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
turnPrintTSOff(void)
{
//  pthread_mutex_lock(&printLock);

  doPrintTS = 0;

//  pthread_mutex_unlock(&printLock);
}

/*---------------------------------------------------------------------------
 *  turnPrintTSOn
 */
/** \brief Turn thread safe printing on.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
turnPrintTSOn(void)
{
//  pthread_mutex_lock(&printLock);

  doPrintTS = 1;

//  pthread_mutex_unlock(&printLock);

}

/*---------------------------------------------------------------------------
 *  isLogTSOn
 */
/** \brief Is thread safe printing on?
 *
 *  \return On/Off indicator. */
/*---------------------------------------------------------------------------*/
int
isLogTSOn(void)
{

  /* I think we can get away without locking here. */

  return doLogTS;
}

/*---------------------------------------------------------------------------
 *  turnLogTSOff
 */
/** \brief Turn thread safe printing off.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
turnLogTSOff(void)
{
//  pthread_mutex_lock(&printLock);

  doLogTS = 0;

//  pthread_mutex_unlock(&printLock);
}

/*---------------------------------------------------------------------------
 *  turnLogTSOn
 */
/** \brief Turn thread safe printing on.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
turnLogTSOn(void)
{
//  pthread_mutex_lock(&printLock);

  doLogTS = 1;

//  pthread_mutex_unlock(&printLock);

}

/*---------------------------------------------------------------------------
 *  getPrintLevel
 */
/** \brief Return the current printLevel (See description in initilaizePrintingTS)
 *
 *  \return int: printLevel. */
/*---------------------------------------------------------------------------*/
int
getPrintLevelTS(void)
{
  /* I think we can get away without locking here. */

  return s_printLevel;
}

/*---------------------------------------------------------------------------
 *  setPrintLevel
 */
/** \brief Set the printlevel, usually done at initialization. Default is 0.
 *
 *  \return void */
/*---------------------------------------------------------------------------*/
void
setPrintLevelTS(int printLevel)
{
//  pthread_mutex_lock(&printLock);

  s_printLevel = printLevel;

//  pthread_mutex_unlock(&printLock);
}

