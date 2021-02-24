#ifndef INCLUDED_COOP_PRINTINGTS_H
#define INCLUDED_COOP_PRINTINGTS_H

/*===========================================================================*/
/**
 *  \file PrintingTS.h
 *
 *  \brief Tility routines for printing output to stdout, stderr, and log
 *  files.
 */
/*  Copyright:  (c) 2004-2007 The Regents of the University of California.
 *
 *===========================================================================*/

#ifdef __cplusplus
extern "C"
{
#endif

  /** \brief Create directory to hold log files if needed; open log file.               */
  /*  printLevel gives level of output to allow.  The higher the number, the            */
  /*  more output is printed.                                                           */
  /*  0 only prints fatal errors                                                        */
  /*  1 prints warnings and fatal errors                                                */
  /*  2 prints output of interest for later viewing (non-interactive)                   */
  /*  3 and above prints everything (for interactive users)                             */
  /*  logFileName may be null if only printing (not logging) is to be enabled.          */
  /*                                                                                    */
  /* Calling initializePrintingTS is not necessary if you do not plan to do logging.    */
  /* If initilizePrinting is not called, output to terminal will occur, and the default */
  /* print level of 0.                                                                  */
  /*------------------------------------------------------------------------------------*/
  void initializePrintingTS(int printLevel, const char* logFileName, int myRank);

  /** \brief Closes the log file. */
  /*--------------------------------*/
  void finalizePrintingTS(void);

  /** \brief Prints error message to log file & stderr. */
  /*  The printLevel gives the importance */
  /*------------------------------------------------------*/
  void printErrTS(int printLevel, const char* format, ...);

  /** \brief Prints log message to log file. */
  /*-------------------------------------------*/
  void printLogTS(int printLevel, const char* format, ...);

  /** \brief Prints message to log file & stdout. */
  /*------------------------------------------------*/
  void printOutTS(int printLevel, const char* format, ...);

  /** \brief Is thread safe printing on? */
  /*---------------------------------------*/
  int isPrintTSOn(void);

  /** \brief Turn thread safe printing off. */
  /*------------------------------------------*/
  void turnPrintTSOff(void);

  /** \brief Turn thread safe printing on. */
  /*-----------------------------------------*/
  void turnPrintTSOn(void);

  /** \brief Is thread safe logging on? */
  /*---------------------------------------*/
  int isLogTSOn(void);

  /** \brief Turn thread safe logging off. */
  /*------------------------------------------*/
  void turnLogTSOff(void);

  /** \brief Turn thread safe logging on. */
  /*-----------------------------------------*/
  void turnLogTSOn(void);

  /** \brief Return the current printLevel (See description in initilaizePrintingTS) */
  /*---------------------------------------*/
  int getPrintLevelTS(void);

  /** \brief Sets the current printLevel. */
  /*  printLevel gives level of output to allow.  The higher the number, the            */
  /*  more output is printed.                                                           */
  /*  0 only prints fatal errors                                                        */
  /*  1 prints warnings and fatal errors                                                */
  /*  2 prints output of interest for later viewing (non-interactive)                   */
  /*  3 and above prints everything (for interactive users                              */
  /*------------------------------------------*/
  void setPrintLevelTS(int localPrintLevel);



#ifdef __cplusplus
}
#endif


#endif
