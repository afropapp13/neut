ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk
#include config_neut.gmk

#.SILENT:

NEUTCOREVER	= 5.3.6
INEUTCOREVER= 536
NUCEFFVER   = 5.3.5
INUCEFFVER  = 535
NUCCORVER   = 5.3.5
INUCCORVER  = 535
PARTNUCKVER = 5.3.5
IPARTNUCKVER= 535
SKMCSVCVER  = 5.3.5
ISKMCSVCVER = 535
ZBSFNSVER   = 5.3.5
IZBSFNSVER  = 535
SPECFUNCVER = 5.3.5
ISPECFUNCVER= 535
RADCORRVER  = 5.3.5
IRADCORRVER = 535

# compiler
#ifndef ${FC}
##FC = g77
#FC = gfortran
##STOP
#endif

# executables
MAKETABLES = makeTables.o combineTables.o

LIBDIRS = -L. -L${CERN}/2005/lib -L../../lib/Linux_pc

LIBS = -lspecfunc_${SPECFUNCVER} -lneutcore_${NEUTCOREVER} -lskmcsvc_${SKMCSVCVER} -lnuccorrspl_${NUCCORVER} -lnuceff_${NUCEFFVER} -lpartnuck_${PARTNUCKVER} -lpawlib -lpacklib -lmathlib -lgraflib -lkernlib -lstdc++

FLAGS = -pedantic-errors -Wall -W -w -Wunused -Wuninitialized -DstrictF77 $(FCOPTIONS) # -fdefault-real-8

INCLUDES = -I../../inc

# libraries
OBJECTS = fourVector.o calcLH.o selectSfValues.o buildSf.o maxDiff.o readXsecData.o sfevent.o # bbba05.o bbba07.o 

all: makeTables

makeTables: $(MAKETABLES) $(OBJECTS) libspecfunc_${SPECFUNCVER}.a
	$(FC) $(FLAGS)-o makeTables.exe  $(MAKETABLES) $(INCLUDES) $(LIBDIRS) $(LIBS)


#bbba07.o:
#	rm -rf neutmodelC.h
#	rm bbba07.cc
#	ln -s ../neutcore/bbba07.cc bbba07.cc
#	ln -s ../neutcore/neutmodelC.h neutmodelC.h
#	gcc -D$(FC) $(FLAGS) -g -c bbba07.cc -o bbba07.o
    
%.o: %.F 
	$(FC) $(FLAGS)$ ${INCLUDES} -g -c -o $@ $<

libspecfunc_${SPECFUNCVER}.a: $(OBJECTS)
#	rm bbba05.F
#	ln -s ../neutcore/bbba05.F bbba05.F
#	ar r libspecfunc.a *.o
	ar r $@ $(OBJECTS)

clean:
	rm -rf *.o *.a *.exe

install: libspecfunc_${SPECFUNCVER}.a
	cp libspecfunc_${SPECFUNCVER}.a ../../lib/Linux_pc

