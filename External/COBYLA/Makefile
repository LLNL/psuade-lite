# Makefile for Cobyla
# $Jeannot: Makefile,v 1.8 2004/04/13 16:02:01 js Exp $

#CC = gcc
#CFLAGS = -ffloat-store

include ../../config/mh-machine

LIB = libcobyla2c.a
OBJ = cobyla.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

LIB: $(OBJ)
	ar ruv libcobyla2c.a cobyla.o

clean:
	rm -f libcobyla2c.a
	rm -f *.o 
