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
// Defintion for the class CommMPICH
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#ifndef __COMMMPICHH__
#define __COMMMPICHH__

#ifdef HAVE_MPICH
#include <mpi.h>
#else
#define MPI_Request int
#endif

#include "BaseComm.h"
#include "dtype.h"

// *******************************************************************
// Definition of the class CommMPICH 
// *******************************************************************

class CommMPICH : public BaseComm 
{
   int         mypid_;
   int         numProcs_;    
   MPI_Request *requests_;
   MPI_Group   group_;

public:

   CommMPICH(int argc, void **argv);
   ~CommMPICH();
	
   inline int getPID();
   inline int getNumProcs();

   void synchronize();
   void shutdown();

   int  send(void*, int, int, int, int);
   int  recv(void*, int, int , int, int);
   int  iRecv(void*, int, int, int, int);
   int  disableIrecv(int);
   int  iProbe(int, int);
   void wait(int);
   int  waitAny();
 
   void bcast(void *, int, int, int);
   void allReduce(void*, int, int, char);

   MPI_Datatype setDataType(int);
};

#endif // __COMMMPICHH__

