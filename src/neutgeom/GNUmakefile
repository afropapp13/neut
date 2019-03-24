ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk


NEUT_SRC = $(NEUT_ROOT)/src

#FC = g77
#FC = gfortran

ROOTINC = -I`root-config --incdir`
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff
CLASSDIR = $(PWD)/../neutclass
CLASSINC= -I${CLASSDIR}
INCDIRS = ${ROOTINC} ${NEUTINC} ${CLASSINC}

FCOPTIONS = -g

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive


COREDIR =	../neutcore/Linux_pc

ROOTLIBS  = `root-config --libs` -lCore -lCint -lEG -lPhysics -lRIO -lNet \
        -lTree -lGeom -lGraf -lGraf3d -lHist -lMatrix -lMinuit -lPostscript \
	-lTree -lMathCore -lGpad -lGui -lGX11 -lRint -lThread

LIBS      = ${NEUTLIBS} $(ROOTLIBS) ${CERNLIBS} -lstdc++

LDOPTFLAGS= -g

MOBJS = ${NEUT_SRC}/neutsmpl/nevecgen.o ${NEUT_SRC}/neutsmpl/grndmq.o ${NEUT_SRC}/neutsmpl/rmarin_dum.o

OBJS = neutev.o neutxs.o TNuTrajectory.o TNeutOutput.o TNuFlux.o

.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
			${NEUT_SRC}/neutsmpl/bin/gen_root_so.sh

all: event_rate genev

event_rate: event_rate.o $(OBJS) 
	$(FC) -O3 -o $@ event_rate.o -Xlinker -R`pwd`  $(OBJS) ${MOBJS} $(LIBS)

genev: genev.o $(OBJS)
	$(FC) -O3 -o $@ genev.o -Xlinker -R`pwd`  $(OBJS) ${MOBJS} $(LIBS)

dumptotpau_nd280: dumptotpau_nd280.o $(OBJS)
	$(FC) -O3 -o $@ dumptotpau_nd280.o -Xlinker -R`pwd` ${OBJS} ${MOBJS} $(LIBS)

clean:
	$(RM) -f *.o *~ ${OBJS} event_rate genev dumptotpau_nd280

