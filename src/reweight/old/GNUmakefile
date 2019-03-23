ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

ROOTINC = -I`root-config --incdir`
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff -I../../include
INCDIRS = ${ROOTINC} ${NEUTINC}

FCOPTIONS = -g

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive


COREDIR =	../neutcore/Linux_pc

LIBDIRS =	-L${COREDIR}         -L../nuccorspl/Linux_pc \
			-L../nuceff/Linux_pc -L../partnuck/Linux_pc \
			-L../skmcsvc/Linux_pc

NEUTCOREVER	= 5.1.0
INEUTCOREVER= 510

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
LIBS      = $(ROOTLIBS) ${MCLIB} ${CLIBS} -lstdc++

LDOPTFLAGS= -g

FINCDIRS  = ${NEUTINC}

POBJS = ${COREDIR}/structm.o ${COREDIR}/pdfset.o ${COREDIR}/grv94di.o \
		${COREDIR}/grv98_lo.o
MOBJS = t2kneutreweight.o

OBJS  = ${MOBJS} ${ROBJS} ${SOBJS} ${POBJS}


.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
			./bin/gen_root_so.sh

.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(FINCDIRS) -o $@ $<

runT2Krw: runT2Krw.o $(OBJS)
	$(FC) -o $@ runT2Krw.o -Xlinker -R`pwd` $(OBJS) $(LIBS)

clean:
	$(RM) -f *.o *~ ${MOBJS} runT2Krw.o runT2Krw
