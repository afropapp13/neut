ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk


NEUT_SRC = $(NEUT_ROOT)/src

#FC = g77
FC = gfortran

ROOTINC = -I`root-config --incdir`
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff
CLASSDIR = $(PWD)/../neutclass
CLASSINC= -I${CLASSDIR}
INCDIRS = ${ROOTINC} ${NEUTINC} ${CLASSINC}

FCOPTIONS = -g

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive


COREDIR =	../neutcore/Linux_pc

LIBDIRS =	-L${COREDIR}         -L../nuccorspl/Linux_pc \
			-L../nuceff/Linux_pc -L../partnuck/Linux_pc \
                      -L../skmcsvc/Linux_pc -L../specfunc/Linux_pc -L../specfunc

NEUTCOREVER	= 5.3.5
INEUTCOREVER= 535
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

FORTRANDEFINES += -DNECORE_VERSION=$(INEUTCOREVER)
FORTRANDEFINES += -DNENUCE_VERSION=$(INUCEFFVER)
FORTRANDEFINES += -DNENUCC_VERSION=$(INUCCORVER)

MCLIB = ${LIBDIRS}   -lneutcore_${NEUTCOREVER} -lnuceff_${NUCEFFVER} \
		     -lneutcore_${NEUTCOREVER} -lnuccorrspl_${NUCCORVER} \
                     -lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER} -lspecfunc_${SPECFUNCVER}


CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
		${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
		`cernlib jetset74 photos202 mathlib packlib kernlib`

ROOTLIBS  = `root-config --libs` -lCore -lCint -lEG -lPhysics -lRIO -lNet \
        -lTree -lGeom -lGraf -lGraf3d -lHist -lMatrix -lMinuit -lPostscript \
	-lTree -lMathCore -lGpad -lGui -lGX11 -lRint -lThread

LIBS      = $(ROOTLIBS) ${MCLIB} ${CLIBS} -lstdc++

LDOPTFLAGS= -g

FINCDIRS  = ${NEUTINC}

SOBJS = ${CLASSDIR}/neutctrl.so ${CLASSDIR}/neutpart.so \
		${CLASSDIR}/neutvect.so ${CLASSDIR}/neutvtx.so

ROBJS = ${CLASSDIR}/neutfill.o \
		${CLASSDIR}/event_ratefortwrapper.o \
		${CLASSDIR}/NeutRootHandlers.o 

POBJS = ${COREDIR}/structm.o ${COREDIR}/pdfset.o ${COREDIR}/grv94di.o \
		${COREDIR}/grv98_lo.o
MOBJS = ${NEUT_SRC}/neutsmpl/nevecgen.o ${NEUT_SRC}/neutsmpl/grndmq.o ${NEUT_SRC}/neutsmpl/rmarin_dum.o

OBJS = neutev.o neutxs.o TNuTrajectory.o TNeutOutput.o TNuFlux.o


.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
			${NEUT_SRC}/neutsmpl/bin/gen_root_so.sh

.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cxx.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(FINCDIRS) -o $@ $<

all: event_rate genev

event_rate: event_rate.o $(OBJS) 
	$(FC) -O3 -o $@ event_rate.o -Xlinker -R`pwd`  $(OBJS) ${MOBJS} $(LIBS)

genev: genev.o $(OBJS)
	$(FC) -O3 -o $@ genev.o -Xlinker -R`pwd`  $(OBJS) ${MOBJS} $(LIBS)

dumptotpau_nd280: dumptotpau_nd280.o $(OBJS)
	$(FC) -O3 -o $@ dumptotpau_nd280.o -Xlinker -R`pwd` ${OBJS} ${MOBJS} $(LIBS)

clean:
	$(RM) -f *.o *~ ${OBJS} event_rate genev dumptotpau_nd280
