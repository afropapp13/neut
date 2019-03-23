ifndef NEUT_ROOT
      NEUT_ROOT = ../../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk
#include config_neut.gmk



ROOTINC = -I`root-config --incdir`
NEUTINC = -I../../neutcore -I../../nuccorspl -I../skmcsvc -I../nuceff -I../../include
CLASSDIR = $(PWD)/../neutclass
CLASSINC= -I${CLASSDIR}
INCDIRS = ${ROOTINC} ${NEUTINC} ${CLASSINC}

FCOPTIONS += -g -O

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive


COREDIR =	../../neutcore/Linux_pc

LIBDIRS =	-L${COREDIR}         -L../../nuccorspl/Linux_pc \
			-L../../nuceff/Linux_pc -L../../partnuck/Linux_pc \
			-L../../skmcsvc/Linux_pc

NEUTCOREVER	= 5.3.1
INEUTCOREVER= 531

NUCEFFVER   = 5.3.1
INUCEFFVER  = 531

NUCCORVER   = 5.3.1
INUCCORVER  = 531

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

SOBJS = ${CLASSDIR}/neutctrl.so ${CLASSDIR}/neutpart.so \
		${CLASSDIR}/neutvect.so ${CLASSDIR}/neutvtx.so \
		${CLASSDIR}/neutfsipart.so ${CLASSDIR}/neutfsivert.so 
#${CLASSDIR}/neutfillvect.o ${CLASSDIR}/neutfillvtx.o 
ROBJS = ${CLASSDIR}/neutfill.o \
		${CLASSDIR}/neutrootfortwrapper.o \
		${CLASSDIR}/NeutRootHandlers.o 

POBJS = ${COREDIR}/structm.o ${COREDIR}/pdfset.o ${COREDIR}/grv94di.o \
		${COREDIR}/grv98_lo.o
#MOBJS = nevecgen.o grndmq.o rmarin_dum.o flx_intrfc.o ufm2th1dist.o

#OBJS  = ${MOBJS} ${ROBJS} ${SOBJS} ${POBJS}


.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
			./bin/gen_root_so.sh

.cc.o:
	$(CXX) -D$(FC) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(FINCDIRS) -o $@ $<

clean:
	$(RM) -f *.o *~ ${MOBJS} crs crs_newff

DBGOBJS = trapfpe.o
#LIBOBJS = rsdcrs.o res_lambda.o res_q2_limits.o $(DBGOBJS)
LIBOBJS = res_lambda.o res_q2_limits.o $(DBGOBJS)

crs: crs.o $(OBJS) $(LIBOBJS) 
	$(FC) -o $@ crs.o -Xlinker -R`pwd`  $(LIBOBJS) $(OBJS) $(LIBS)

crs_newff: crs_newff.o $(OBJS) $(LIBOBJS) 
	$(FC) -o $@ crs_newff.o -Xlinker -R`pwd`  $(LIBOBJS) $(OBJS) $(LIBS)

