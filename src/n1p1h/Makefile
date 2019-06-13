ifdef NEUT_ROOT
	include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk
else
	include inc/N1P1HVER.inc
	INSTALL   = /usr/bin/install
	FH2H = /bin/echo
	FC   = gfortran
	CXX  = g++
	CC   = gcc
	RM   = /bin/rm

	FCOPTIONS = -fno-backslash -fno-automatic -fno-second-underscore \
		    -ffixed-line-length-132 -fdollar-ok -fpic -fbounds-check \
		    -Wall -W -ffpe-trap=invalid,zero,overflow \
		    -DstrictF77 \
		    -DgFortran -Df2cFortran
	FCOPTFLAGS= -g

	CXXOPTIONS  = `root-config --cflags`
	CXXOPTFLAGS = -g -Wno-write-strings -fpermissive

	COPTIONS    = `root-config --cflags`
	COPTFLAGS   = -g -Wno-write-strings -fpermissive
endif

LIBNAME=N1p1h_$(N1P1HVER)

####################################################################
#CFLAGS = -pthread  -m64 `root-config --cflags` -Iinc
#LDFLAGS = -Llib -lN1p1h `root-config --glibs` -lm -ldl -lgfortran 

#CFLAGS=$(COPTFLAGS) $(CXXOPTIONS) $(ROOTINCDIR) $(NEUTINCDIR) -pthread -std=c++11 -m64 -Iinc
CFLAGS=$(COPTFLAGS) $(CXXOPTIONS) $(ROOTINCDIR) $(NEUTINCDIR) -pthread -m64 -Iinc -fPIC 
FFLAGS=$(FCOPTFLAGS) $(FCFLAGS) $(ROOTINCDIR) -Iinc -fPIC
LDFLAGS= -L${NEUT_ROOT}/lib -lneutcore $(ROOTLIBS)

ifdef NEUT_ROOT
CFLAGS+=-DWITH_NEUT
endif

####################################################################

INCS = inc/nieves1p1h.fh
$(shell mkdir -p lib)

all: lib/lib${LIBNAME}.a 

lib/Nucleus.o: src/Nucleus.C 
	${CXX} $(CFLAGS) -c -Iinc/ -funroll-loops src/Nucleus.C -o lib/Nucleus.o

bin/testNuclei: tst/testNuclei.C lib/lib${LIBNAME}.a
	${CXX} $(CFLAGS) -funroll-loops tst/testNuclei.C -o bin/testNuclei -Llib -l${LIBNAME} -L/usr/local/lib  -lgfortran $(CFLAGS) $(LDFLAGS) $(CERNLIB)

lib/qe.o : src/qe.F
	${FC}  $(FFLAGS) -funroll-loops -c src/qe.F -o lib/qe.o

lib/N1p1h.o : src/N1p1h.C	
	${CXX} $(CFLAGS) -funroll-loops -c src/N1p1h.C -o lib/N1p1h.o

lib/N1p1h2d.o : src/N1p1h2d.C	
	${CXX} $(CFLAGS) -funroll-loops -c src/N1p1h2d.C -o lib/N1p1h2d.o

lib/N1p1hFortranInterface.o : src/N1p1hFortranInterface.C
	${CXX} $(CFLAGS) -funroll-loops -c src/N1p1hFortranInterface.C \
								-o lib/N1p1hFortranInterface.o

lib/lib${LIBNAME}.a: lib/qe.o lib/Nucleus.o lib/N1p1h.o lib/N1p1h2d.o lib/N1p1hFortranInterface.o
	${AR} rcs lib/lib${LIBNAME}.a lib/N1p1h.o lib/N1p1h2d.o lib/Nucleus.o lib/qe.o \
								lib/N1p1hFortranInterface.o

bin/Test: tst/Test.C lib/lib${LIBNAME}.a
	${CXX} $(CFLAGS) -funroll-loops tst/Test.C -o bin/Test  $(LDFLAGS) \
		-Llib -l${LIBNAME} -lgfortran $(CERNLIBS)

bin/SearchMaximum: tst/SearchMaximum.C lib/lib${LIBNAME}.a
	${CXX} $(CFLAGS) -funroll-loops tst/SearchMaximum.C -o bin/SearchMaximum  $(LDFLAGS) -Llib -l${LIBNAME} -lgfortran $(CERNLIBS) $(ROOTLIBS) -lMinuit

install.library: lib/lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib

install.include: ${INCS}
	for i in $? ;do ${FH2H} $$i > `echo $$i | sed "s/\.fh$$/C.h/"`; done
	${INSTALL} -m 0644 ${INCS} ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS:.fh=C.h} ${NEUT_ROOT}/include


