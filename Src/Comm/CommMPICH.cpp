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
// Functions for the class CommMPICH
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ------------------------------------------------------------------------
// Functions provided in this file :
//
// Point to Point 
// Communication : send         - send message
//                 recv         - receive message
//                 iRecv        - nonblocking receive
//                 disableIrecv - disable processor
//                 iProbe       - probe the processor for a message
//                 wait         - wait for a posted iRecv
//                 waitAny      - wait for a any posted iRecv
// Global 
// Communication : synchronize - synchronize all processors
//                 bcast       - broadcast to all processors
//                 allReduce   - reduce(add/max) and broadcast
//
// Others        : shutdown
// ************************************************************************
// ************************************************************************

#if defined(HAVE_MPICH) 

#include <stdio.h>
#include <stdlib.h>
#include "CommMPICH.h"
#include "dtype.h"

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
CommMPICH::CommMPICH(int argc, void **argv)
{
   // get MPI information
   ::MPI_Init(&argc, (char ***) &argv);
   ::MPI_Comm_rank(MPI_COMM_WORLD, &mypid_); 
   ::MPI_Comm_size(MPI_COMM_WORLD, &numProcs_);

   // set up for requests
   requests_ = new MPI_Request[numProcs_];
   for (int i = 0; i < numProcs_; i++) requests_[i] = MPI_REQUEST_NULL;
}

// ************************************************************************
// Destructor
// ------------------------------------------------------------------------
CommMPICH::~CommMPICH()
{
   shutdown();
   if (requests_ != NULL) delete [] requests_;
}

// ************************************************************************
// get processor ID 
// ------------------------------------------------------------------------
int CommMPICH::getPID()
{
   return mypid_;
}

// ************************************************************************
// get number of processors
// ------------------------------------------------------------------------
int CommMPICH::getNumProcs()
{
   return numProcs_;
}

// ************************************************************************
// synchronize all processors
// ------------------------------------------------------------------------
void CommMPICH::synchronize()
{
   ::MPI_Barrier(MPI_COMM_WORLD);
}

// ************************************************************************
// Send message to a given destination
// ------------------------------------------------------------------------
int CommMPICH::send(void *msg,int length,int dtype,int msgID,int dest) 
{
   MPI_Datatype  mpiDtype;

   mpiDtype = setDataType(dtype);
   ::MPI_Send(msg, length, mpiDtype, dest, msgID, MPI_COMM_WORLD);
   return 0;
}

// ************************************************************************
// receive data from a given source processor
// ------------------------------------------------------------------------
int CommMPICH::recv(void *msg,int length,int dtype,int msgID,int src)
{
   MPI_Datatype mpiDtype;
   MPI_Status   status;

   mpiDtype = setDataType(dtype);
   if (src == -1)
   {
      ::MPI_Recv(msg, length, mpiDtype, MPI_ANY_SOURCE, msgID,
                 MPI_COMM_WORLD, &status);
      return status.MPI_SOURCE;
   } 
   else 
   {
      ::MPI_Recv(msg,length,mpiDtype,src,msgID,MPI_COMM_WORLD,&status);
      return src;
   }
}

// ************************************************************************
// post a nonblocking receive from a given processor
// ------------------------------------------------------------------------
int CommMPICH::iRecv(void *msg,int length,int dtype,int msgID,int src)
{
   MPI_Datatype  mpiDtype;

   mpiDtype = setDataType(dtype);
   if ( requests_[src] != MPI_REQUEST_NULL )
   {
      printf("CommMPICH::iRecv ERROR - pending iRecv not processed.\n");
      exit(-1);
   }
   ::MPI_Irecv(msg, length, mpiDtype, src, msgID, MPI_COMM_WORLD,
               &(requests_[src]));
   return src;
}

// ************************************************************************
// disableIrecv
// ------------------------------------------------------------------------
int CommMPICH::disableIrecv(int procID)
{
   if ( procID > 0 && procID < numProcs_ )
   {
      requests_[procID] = MPI_REQUEST_NULL;
      return 0;
   }
   else return -1;
}

// ************************************************************************
// searches for an incoming message that matches src, tag, and
// comm.  If found, return flag = 1, and status will contain the
// same information as in MPI_Recv.
// ------------------------------------------------------------------------
int CommMPICH::iProbe(int src, int msgID) 
{
   int        flag;
   MPI_Status status;

   if ( src == -1 )
   {
      MPI_Iprobe(MPI_ANY_SOURCE, msgID, MPI_COMM_WORLD, &flag, &status);
      if ( flag != 0 ) return status.MPI_SOURCE;
      else             return -1;
   }
   else
   {
      MPI_Iprobe(src, msgID, MPI_COMM_WORLD, &flag, &status);
      return flag;
      //if ( flag != 0 ) return src;
      //else             return -1;
   }
}

// ************************************************************************
// wait for a request to be completed
// ------------------------------------------------------------------------
void CommMPICH::wait(int procID)
{
   MPI_Status status;

   if ( procID < 0 || procID > numProcs_ )
   {
      printf("CommMPICH::wait - handler not valid.\n");
      exit(-1);
   }
   ::MPI_Wait( &requests_[procID], &status);
   requests_[procID] = MPI_REQUEST_NULL;
} 

// ************************************************************************
// wait for any request to be completed
// ------------------------------------------------------------------------
int CommMPICH::waitAny()
{
   int        errStat, procID;
   MPI_Status status;

   errStat = MPI_Waitany(numProcs_, requests_, &procID, &status);
   if ( errStat != MPI_SUCCESS )
   {
      printf("CommMPICH::waitAny fatal ERROR\n");
      exit(-1);
   }
   return procID;
}

// ************************************************************************
// broadcast data 
// ------------------------------------------------------------------------
void CommMPICH::bcast(void *msg, int length, int dtype, int root)
{
   MPI_Datatype  mpiDtype;

   mpiDtype = setDataType(dtype);
   ::MPI_Bcast(msg, length, mpiDtype, root, MPI_COMM_WORLD);
}

// ************************************************************************
// reduce and broadcast data using MPI primitive
// ------------------------------------------------------------------------
void CommMPICH::allReduce(void *msg,int length,int dtype,char op)
{
   int           i, nb;
   MPI_Datatype  mpiDtype;
   char          *msg2, *msg3;

   mpiDtype = setDataType(dtype);
   switch(dtype) 
   {
      case INT :    nb = sizeof(int)*length;    break;
      case FLOAT :  nb = sizeof(float)*length;  break;
      case DOUBLE : nb = sizeof(double)*length; break;
      case CHAR :   nb = sizeof(char)*length;   break;
   }
   msg2 = new char[nb];
   switch(op) {
      case '+' :    
         ::MPI_Allreduce(msg,msg2,length,mpiDtype,MPI_SUM,MPI_COMM_WORLD);
         break;
      case 'm' :    
         ::MPI_Allreduce(msg,msg2,length,mpiDtype,MPI_MAX,MPI_COMM_WORLD);
         break;
   }
   msg3 = (char *) msg;
   for (i = 0; i < nb; i++) msg3[i] = msg2[i];
   delete [] msg2;
}

// ************************************************************************
// exiting from the MPI
// ------------------------------------------------------------------------
void CommMPICH::shutdown()
{
   ::MPI_Finalize ();
}

// ************************************************************************
// setting MPI data type 
// ------------------------------------------------------------------------
MPI_Datatype CommMPICH::setDataType(int intype)
{
   MPI_Datatype mpiDtype;

   switch(intype) 
   {
      case INT :    mpiDtype = MPI_INT;    break;
      case FLOAT :  mpiDtype = MPI_FLOAT;  break;
      case DOUBLE : mpiDtype = MPI_DOUBLE; break;
      case CHAR   : mpiDtype = MPI_CHAR;   break;
      default :
         printf("CommMPICH::setDataType ERROR : wrong data type %d\n",
                intype);
         exit(-1);
         break;
   }
   return mpiDtype;
}
#endif

