# Feb 17 2015, converted fortran to C
# Stephen Becker, stephen.becker@colorado.edu


FC = gfortran
CC = gcc
#CC = gcc-4.9

#FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 
#FFLAGS += -fdefault-integer-8
CFLAGS = -O3 -Wall -Wno-uninitialized 
#CFLAGS = -g -Wall -lm -DDEBUG -Wuninitialized -Wextra -Wno-unused-parameter
LDFLAGS = -lm

# my f2c.h is here:
INCLUDES = -I/usr/local/include -I./

DRIVER1_77 = driver1.c
DRIVER2_77 = driver2.c
DRIVER3_77 = driver3.c

LBFGSB  = lbfgsb.c linesearch.c subalgorithms.c print.c

LINPACK = linpack.c

BLAS 	= miniCBLAS.c
#CFLAGS += -D_USE_OPTIMIZED_BLAS -lblas

TIMER   = timer.c

SRC = $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(F2CFILES)

%all :  lbfgsb_77_1 lbfgsb_77_2 lbfgsb_77_3

OBJS = lbfgsb.o linesearch.o subalgorithms.o print.o linpack.o timer.o miniCBLAS.o
default : $(OBJS)
	$(AR) ruv liblbfgs.a $(OBJS)
	true liblbfgs.a
test1: lbfgsb_77_1

clean:
	/bin/rm -f *.a *.o

lbfgsb.o: lbfgsb.c
	$(CC) $(CFLAGS) -c lbfgsb.c
linesearch.o: linesearch.c
	$(CC) $(CFLAGS) -c linesearch.c
subalgorithms.o: subalgorithms.c
	$(CC) $(CFLAGS) -c subalgorithms.c
print.o: print.c
	$(CC) $(CFLAGS) -c print.c
linpack.o: linpack.c
	$(CC) $(CFLAGS) -c linpack.c
timer.o: timer.c
	$(CC) $(CFLAGS) -c timer.c
miniCBLAS.o: miniCBLAS.c
	$(CC) $(CFLAGS) -c miniCBLAS.c
	
lbfgsb_77_1 : $(DRIVER1_77C) $(SRC)
	$(CC) $(CFLAGS) $(DRIVER1_77) $(SRC) $(LDFLAGS) -o x.lbfgsb_77_1

lbfgsb_77_2 : $(DRIVER2_77C) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) $(CFLAGS) $(DRIVER2_77) $(SRC) $(LDFLAGS) -o x.lbfgsb_77_2

lbfgsb_77_3 : $(DRIVER2_77C) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(CC) $(CFLAGS) $(DRIVER3_77) $(SRC) $(LDFLAGS) -o x.lbfgsb_77_3

test_1 : x.lbfgsb_77_1
	./x.lbfgsb_77_1
