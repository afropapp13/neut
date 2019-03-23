ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

NEUT_SRC = $(NEUT_ROOT)/src

ROOTINC = -I`root-config --incdir`
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff -I../../include -I../neutgeom
INCDIRS = ${ROOTINC} ${NEUTINC}

FCOPTIONS = -g

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive


COREDIR =	../neutcore/Linux_pc

LIBDIRS =	-L${COREDIR}         -L../nuccorspl/Linux_pc \
			-L../nuceff/Linux_pc -L../partnuck/Linux_pc \
			-L../skmcsvc/Linux_pc

NEUTCOREVER = 5.3.2.2
INEUTCOREVER= 5322

NUCEFFVER   = 5.3.2.2
INUCEFFVER  = 5322

NUCCORVER   = 5.3.2.2
INUCCORVER  = 5322

PARTNUCKVER = 5.3.2.2
IPARTNUCKVER= 5322

SKMCSVCVER  = 5.3.2.2
ISKMCSVCVER = 5322

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

POBJS = $(NEUT_SRC)/neutsmpl/grndmq.o $(NEUT_SRC)/neutsmpl/rmarin_dum.o

MOBJS =  piscatana.o

OBJS  = ${MOBJS} ${POBJS}


.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
			./bin/gen_root_so.sh

.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(FINCDIRS) -o $@ $<

all: piscat gampi run_piscatana run_piscatplot

piscat: piscat.o $(POBJS)
	$(FC) -o $@ piscat.o -Xlinker -R`pwd` $(POBJS) $(LIBS)

scatt_test: scatt_test.o $(POBJS)
	$(FC) -o $@ scatt_test.o -Xlinker -R`pwd` $(POBJS) $(LIBS)

gampi: gampi.o $(POBJS)
	$(FC) -o $@ gampi.o -Xlinker -R`pwd` $(POBJS) $(LIBS)

run_piscatana: run_piscatana.o $(OBJS)
	$(FC) -o $@ run_piscatana.o -Xlinker -zmuldefs  $(OBJS) $(LIBS)

run_piscatplot: run_piscatplot.o $(OBJS)
	$(FC) -o $@ run_piscatplot.o -Xlinker -zmuldefs  $(OBJS) $(LIBS)
clean:
	$(RM) -f *.o *~ ${MOBJS}  run_piscatana piscat  gampi run_piscatplot
