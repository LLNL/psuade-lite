
# The following system parameters are for a Linux host.
# ARCH   = LINUX

CC     = gcc -g -Wall
CXX    = g++ -g -Wall
F77    = f77 
AR     = ar
RANLIB = true

OBJS = dgemv.o   dlange.o  dlaset.o  dlassq.o  dormbr.o  ilaenv.o\
       dbdsqr.o  dgeqr2.o  dlapy2.o  dlasq1.o  dlasv2.o  dorml2.o\
       dcopy.o   dgeqrf.o  dlarfb.o  dlasq2.o  dnrm2.o   dormlq.o\
       dgebd2.o  dger.o    dlarf.o   dlasq3.o  dorg2r.o  dormqr.o\
       dgebrd.o  dgesvd.o  dlarfg.o  dlasq4.o  dorgbr.o  dscal.o \
       dgelq2.o  dlabad.o  dlarft.o  dlasq5.o  dorgl2.o  dtrmm.o \
       dgelqf.o  dlabrd.o  dlartg.o  dlasq6.o  dorglq.o  dtrmv.o \
       dgels.o   dlacpy.o  dlas2.o   dlasr.o   dorgqr.o  dtrsm.o \
       dgemm.o   dlamch.o  dlascl.o  dlasrt.o  dorm2r.o  ieeeck.o \
       lsame.o   xerbla.o  dswap.o   drot.o    dlae2.o   dlaev2.o \
       dlanst.o  dlansy.o  dlatrd.o  dorg2l.o  dorgql.o  dorgtr.o \
       dsteqr.o  dsterf.o  dsyev.o   dsytd2.o  dsytrd.o  blas2test.o \
       daxpy.o   ddot.o    dsymv.o   dsyr2.o   dsyr2k.o  dchdc.o \
       dgesv.o   dgetf2.o  dgetrf.o  dgetrs.o  dsyrk.o   dlaswp.o \
       dsymm.o   dposv.o   dpotf2.o  dpotrf.o  dpotrs.o  idamax.o \
       dtrsv.o   dqrsl.o   dqrdc2.o  dtrsl.o   dasum.o

.f.o:
	$(F77) $(FFLAGS) -c $*.f

./liblapack.a: $(OBJS)
	$(AR) ruv $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o *.a

realclean:
	rm -f *.o ; rm -f ./liblapack.a
