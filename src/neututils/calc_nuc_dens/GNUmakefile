ifndef NEUT_ROOT
      NEUT_ROOT = ../../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

NEUT_SRC = $(NEUT_ROOT)/src 
#include config_neut.gmk

ROOTINC = -I`root-config --incdir`
NEUTINC = -I${NEUT_SRC}/neutcore -I${NEUT_SRC}/nuccorspl -I${NEUT_SRC}/skmcsvc -I${NEUT_SRC}/nuceff -I${NEUT_ROOT}/include
INCDIRS = ${ROOTINC} ${NEUTINC}

FCOPTIONS = -g

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive


COREDIR =	${NEUT_SRC}/neutcore/Linux_pc

LIBDIRS =	-L${COREDIR}         -L${NEUT_SRC}/nuccorspl/Linux_pc \
			-L${NEUT_SRC}/nuceff/Linux_pc -L${NEUT_SRC}/partnuck/Linux_pc \
			-L${NEUT_SRC}/skmcsvc/Linux_pc

NEUTCOREVER	= 5.1.1
INEUTCOREVER= 511

NUCEFFVER   = 5.1.0
INUCEFFVER  = 510

NUCCORVER   = 1.1.0
INUCCORVER  = 110

PARTNUCKVER = 5.0.5
IPARTNUCKVER= 505

SKMCSVCVER  = 5.0.5
ISKMCSVCVER = 505

FORTRANDEFINES += -DNECORE_VERSION=$(INEUTCOREVER)
FORTRANDEFINES += -DNENUCE_VERSION=$(INUCEFFVER)
FORTRANDEFINES += -DNENUCC_VERSION=$(INUCCORVER)

MCLIB = ${LIBDIRS}   -lneutcore_${NEUTCOREVER} -lnuceff_${NUCEFFVER} \
		-lneutcore_${NEUTCOREVER} -lnuccorrspl_${NUCCORVER} \
		-lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER}

CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
		${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
		`cernlib jetset74 photos202 mathlib packlib kernlib`

ROOTLIBS  = `root-config --libs`
LIBS      = $(ROOTLIBS) ${MCLIB} ${CLIBS} 

LDOPTFLAGS= -g

FINCDIRS  = ${NEUTINC}

OBJS = integrate_rho.o

.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(FINCDIRS) -o $@ $<

integrate: integrate.o $(OBJS) 
	$(FC) -o $@ integrate.o -Xlinker -R`pwd`  $(OBJS) $(LIBS)
clean:
	$(RM) -f *.o *~ ${MOBJS} integrate

