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

ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

INSTALL     = /usr/bin/install
ROOTINCDIR  = `root-config --incdir`
LOCALINCDIR = ./ 
INCDIRS     = -I${LOCALINCDIR} -I${ROOTINCDIR} -I${SKOFL_ROOT}/inc \
			  -I${NEUT_ROOT}/include

CXXFLAGS  = -g -O0 `root-config --cflags`

COPTFLAGS   = -g -O0 -I${NEUT_ROOT}/include $(ALLDEFINES)
FOPTFLAGS   = -g -O0 -I${NEUT_ROOT}/include -DFLUX_10C $(ALLDEFINES)

FFLAGS      = ${FOPTFLAGS} ${FCOPTIONS} $(FCFLAGS) 

OBJS        = radcorr.o

LIBRADCORR  = libradcorr_${RADCORRVER}.a

all: ${LIBRADCORR}

${LIBRADCORR}: ${OBJS}
	ar rv $@ ${OBJS}

install: ${LIBRADORR}
	$(INSTALL) ${LIBRADCORR} ${NEUT_ROOT}/lib/Linux_pc 

.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.F.o:
	$(FC) -c $(FFLAGS) $(INCDIRS) -o $@ $<

clean:
	$(RM) -f *.o *.a

veryclean:
	$(RM) -f *.o *.a *~ #* rflist*


