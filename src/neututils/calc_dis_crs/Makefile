CXX = gcc
CXXFLAGS = -O0 -Wall -g
FCFLAGS = -O0 -Wall -g

ifndef NEUT_ROOT
	NEUT_ROOT = ../../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk


TARGETS = calcapicrs


NEUTINC = -I../../neutcore -I../../nuccorspl -I../../skmcsvc -I../../nuceff
COREDIR =	../../neutcore
LIBDIRS =     -L${COREDIR} -L../../../lib/ -L../../nuccorspl/ \
                       -L../../nuceff/ -L../../partnuck/ \
                       -L../../skmcsvc/ -L../../specfunc -L../../specfunc/ \
		       -L../radcorr

MCLIB = ${LIBDIRS} \
	  	-lneutcore_${NEUTCOREVER} -lspecfunc_${SPECFUNCVER} \
		-lnuceff_${NUCEFFVER} \
		-lneutcore_${NEUTCOREVER} -lnuccorspl_${NUCCORVER} \
		-lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER} \
		-lradcorr_${RADCORRVER} 

CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
		${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
		`cernlib jetset74 photos202 mathlib packlib kernlib`

LIBS      = $(ROOTLIBS) ${MCLIB} ${CLIBS} -lstdc++

MOBJS=${NEUT_SRC}/neututils/calc_dis_crs/tpispir.o

calcapicrs: calcapicrs.o tpispir.o
	$(FC) -o $@ calcapicrs.o $(MOBJS) $(LIBS)

calcapicrs.o:
	$(FC) -c $(FCFLAGS) $(NEUTINC) -o $@ $*.F

tpispir.o:
	$(FC) -c $(FCFLAGS) $(NEUTINC) -o $@ $*.F

all: calcapicrs


clean:
	$(RM) -f *.o  calcapicrs
