main: calfun.o getact.o lincoa.o lincob.o main.o prelim.o qmstep.o trstep.o \
      update.o 
	gfortran -o main calfun.o getact.o lincoa.o lincob.o main.o prelim.o \
              qmstep.o trstep.o update.o
main.o: main.f
	gfortran -g  -c main.f
calfun.o: calfun.f
	gfortran -g  -c calfun.f
getact.o: getact.f
	gfortran -g  -c getact.f
lincoa.o: lincoa.f
	gfortran -g  -c lincoa.f
lincob.o: lincob.f
	gfortran -g  -c lincob.f
main.o: main.f
	gfortran -g  -c main.f
prelim.o: prelim.f
	gfortran -g  -c prelim.f
qmstep.o: qmstep.f
	gfortran -g  -c qmstep.f
trstep.o: trstep.f
	gfortran -g  -c trstep.f
update.o: update.f
	gfortran -g  -c update.f

