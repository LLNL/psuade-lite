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
// Functions for the class BaseComm
// AUTHOR : CHARLES TONG
// DATE   : 2003
// ************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "BaseComm.h"

// ************************************************************************
// Constructor
// ------------------------------------------------------------------------
BaseComm::BaseComm()
{
}

// ************************************************************************
// Destructor
// ------------------------------------------------------------------------
BaseComm::~BaseComm()
{
}

// ************************************************************************
// get processor ID
// ------------------------------------------------------------------------
int BaseComm::getPID()
{
   return 0;
}

// ************************************************************************
// get number of processors
// ------------------------------------------------------------------------
int BaseComm::getNumProcs()
{
   return 1;
}

// ************************************************************************
// synchronize all processors
// ------------------------------------------------------------------------
void BaseComm::synchronize()
{
}

// ************************************************************************
// shutdown 
// ------------------------------------------------------------------------
void BaseComm::shutdown()
{
}

// ************************************************************************
// Send message to a given destination
// ------------------------------------------------------------------------
int BaseComm::send(void *msg,int leng,int dtype,int msgid,int dest) 
{
   (void) msg;
   (void) leng;
   (void) dtype;
   (void) msgid;
   (void) dest;
   return -1;
}

// ************************************************************************
// receive data from a given source processor
// ------------------------------------------------------------------------
int BaseComm::recv(void *msg,int leng,int dtype,int msgid,int src)
{
   (void) msg;
   (void) leng;
   (void) dtype;
   (void) msgid;
   (void) src;
   return -1;
}

// ************************************************************************
// post a nonblocking receive from a given processor
// ------------------------------------------------------------------------
int BaseComm::iRecv(void *msg,int leng,int dtype,int msgid,int src)
{
   (void) msg;
   (void) leng;
   (void) dtype;
   (void) msgid;
   (void) src;
   return -1;
}

// ************************************************************************
// disableIrecv
// ------------------------------------------------------------------------
int BaseComm::disableIrecv(int proc)
{
   (void) proc;
   return -1;
}

// ************************************************************************
// iProbe
// ------------------------------------------------------------------------
int BaseComm::iProbe(int proc, int msgID)
{
   (void) proc;
   (void) msgID;
   return -1;
}

// ************************************************************************
// wait for a request to be completed
// ------------------------------------------------------------------------
void BaseComm::wait(int index)
{
   (void) index;
}

// ************************************************************************
// wait for any request to be completed
// ------------------------------------------------------------------------
int BaseComm::waitAny()
{
   return -1;
}

// ************************************************************************
// broadcast data 
// ------------------------------------------------------------------------
void BaseComm::bcast(void *msg, int leng, int dtype, int root)
{
   (void) msg;
   (void) leng;
   (void) dtype;
   (void) root;
}

// ************************************************************************
// reduce and broadcast data using MPI primitive
// ------------------------------------------------------------------------
void BaseComm::allReduce(void *msg,int leng,int dtype,char op)
{
   (void) msg;
   (void) leng;
   (void) dtype;
   (void) op;
}

