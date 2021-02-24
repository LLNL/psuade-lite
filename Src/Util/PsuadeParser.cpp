// ************************************************************************
// Copyright (c) 2007   Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
// Written by the PSUADE team.
// All rights reserved.
//
// Please see the COPYRIGHT_and_LICENSE file for the copyright notice,
// disclaimer, contact information and the GNU Lesser General Public License.
//
// PSUADE is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// PSUADE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// ************************************************************************
// Parser functions
// AUTHOR : CHARLES TONG
// DATE   : 2006
// ************************************************************************

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "PsuadeUtil.h"

/* ************************************************************************
// parse a line a return the tokens separated by spaces
// ***********************************************************************/
int getTokens(char *inString, int *numTokens, char ***tokens, int nlimit)
{
   int  myNumTokens;
   char **myTokens, *tokenPtr, *myString;

   if (inString == NULL)
   {
      (*numTokens) = 0;
      (*tokens) = NULL;
      return -1;
   }
   myString = new char[strlen(inString)+1];
   strcpy(myString, inString);
   myNumTokens = 0;
   tokenPtr = strtok(myString, " ");
   while (tokenPtr != NULL)
   {
      myNumTokens++; 
      tokenPtr = strtok(NULL, " ");
      if (tokenPtr == NULL) break;
      if (tokenPtr[0] == '#') break;
      if (((int) tokenPtr[0]) == 10) break;
   } 
   delete [] myString;
   myTokens = new char*[myNumTokens];
   myNumTokens = 0;
   tokenPtr = strtok(inString, " ");
   while (tokenPtr != NULL)
   {
      myTokens[myNumTokens++] = tokenPtr;
      tokenPtr = strtok(NULL, " ");
      if (tokenPtr == NULL) break;
      if (tokenPtr[0] == '#') break;
      if (((int) tokenPtr[0]) == 10) break;
   } 
   if (nlimit > 0)
   {
      if (myNumTokens != nlimit)
      {
         printf("ERROR: mismatch in the number of parameters.\n");
         myNumTokens = 0;
         delete [] myTokens;
         myTokens = NULL;
         return -1;
      }
   }
   (*numTokens) = myNumTokens;
   (*tokens) = myTokens;
   return 0;
}

/* ************************************************************************
// parse a token and return an integer
// ***********************************************************************/
int getIntFromToken(char *token, int *outdata)
{
   int status=0, leng, ii, intdata, first, last;

   if (token == NULL) status = -1;
   else
   {
      leng = strlen(token);
      if (leng == 0) status = -2;
      first = 0;
      last = leng;
      if (token[0] == '-' || token[0] == '+') first++;
      if (token[leng-1] == 10) last--;
      for (ii = first; ii < last; ii++)
      {
         if (token[ii] < '0' || token[ii] > '9')
         {
            status = -2;
            break;
         }
      }
   }
   if (status == 0)
   {
      sscanf(token, "%d", &intdata);   
      (*outdata) = intdata;
   }
   else printf("ERROR: input is not an integer.\n");
   return status;
}

/* ************************************************************************
// parse a token and return a double
// ***********************************************************************/
int getDbleFromToken(char *token, double *outdata)
{
   int    leng, ii, first, last, index, status=0;
   double dbledata;

   if (token == NULL)
   {
      printf("ERROR: input is not a number.\n");
      return -1;
   }
   leng = strlen(token);
   first = 0;
   if (token[0] == '-' || token[0] == '+') first++;
   last = leng;
   if (token[leng-1] == 10) last--;
   for (ii = first; ii < last; ii++)
   {
      if (token[ii] < '0' || token[ii] > '9') 
      {
         if (token[ii] != 'e' && token[ii] != '.') status = -2;
         break;
      }
   }
   if (status < 0)
   {
      printf("ERROR: input is not a number.\n");
      return -1;
   }
   if (ii < last)
   {
      index = ii;
      if (token[index] == 'e' && index < last-1)
      {
         if (token[index+1] == '-' || token[index+1] == '+')
            for (ii = index+2; ii < last; ii++)
               if (token[ii] < '0' || token[ii] > '9') {status = -3; break;}
      }
      else if (token[index] == 'e' && index >= last-1) status = -4;
      else if (token[index] == '.')
      {
         for (ii = index+1; ii < last; ii++)
         {
            if (token[ii] < '0' || token[ii] > '9') 
            {
               if (token[ii] == 'e') break;
               else                  {status = -5; break;}
            }
         }
         if (ii < last && status == 0)
         {
            index = ii;
            if (index >= last-1) status = -6;
            else
            {
               if (token[index+1] == '-' || token[index+1] == '+')
                  for (ii = index+2; ii < last; ii++)
                     if (token[ii] < '0' || token[ii] > '9')
                     {
                        status = -7; break;
                     } 
            }
         }
      }
   }
   if (status < 0)
   {
      printf("ERROR: input is not a number.\n");
      return -1;
   }
   sscanf(token, "%lg", &dbledata);   
   (*outdata) = dbledata;
   return 0;
}

// ************************************************************************
// get integer
// ------------------------------------------------------------------------
int getInt(int llimit, int uLimit, char *inString)
{
   int  idata, nTokens, iOne=1, status, count;
   char **tokens, winput[501];

   idata = llimit - 1;
   count = 0;
   while (idata < llimit || idata > uLimit)
   {
      printf("%s", inString);
      fgets(winput, 500, stdin);
      status = getTokens(winput, &nTokens, &tokens, iOne);
      if (status >= 0)
      {
         status = getIntFromToken(tokens[0], &idata);
         if (status < 0) idata = llimit - 1;
         else if (idata < llimit || idata > uLimit)
            printf("ERROR: input not within range.\n");
         delete [] tokens;
      }
      else printf("ERROR: Invalid input.\n");
      count++;
      if (count > 10)
      {
         printf("getInt ERROR: bad data entered more than 10 times - abort.\n");
         exit(1);
      }
   }
   return idata;
}

// ************************************************************************
// get double
// ------------------------------------------------------------------------
double getDouble(char *inString)
{
   int    nTokens, iOne=1, status=-1, count;
   double ddata;
   char   **tokens, winput[501];

   count = 0;
   while (status < 0)
   {
      printf("%s", inString);
      fgets(winput, 500, stdin);
      status = getTokens(winput, &nTokens, &tokens, iOne);
      if (status >= 0)
      {
         status = getDbleFromToken(tokens[0], &ddata);
         delete [] tokens;
      }
      count++;
      if (count > 10)
      {
         printf("getDouble ERROR: bad data entered more than 10 times - abort.\n");
         exit(1);
      }
   }
   return ddata;
}

// ************************************************************************
// get character string
// ------------------------------------------------------------------------
int getString(char *inString, char *cstr)
{
  int    nTokens, iOne=1, status=-1, count;
  char   **tokens, winput[501];

  count = 0;
  while (status < 0)
  {
    printf("%s", inString);
    fgets(winput, 500, stdin);
    status = getTokens(winput, &nTokens, &tokens, iOne);
    if (status >= 0)
    {
      strcpy(cstr, tokens[0]);
      cstr[strlen(cstr)] = '\0';
      delete [] tokens;
    }
    count++;
    if (count > 10)
    {
      printf("getString ERROR: bad data entered more than 10 times: abort\n");
      exit(1);
    }
  }
  return 0;
}

