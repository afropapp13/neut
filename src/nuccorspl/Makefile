## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

FCFLAGS += -DNENUCCVER=${INUCCORVER}

LIBBASE	= nuccorspl
LIBNAME = ${LIBBASE}_${NUCCORVER}

OBJS =	nrcfd.o		nrazio.o 	nrinit.o 	nrroxy.o \
		nrllor.o 	nrcone.o 	nrkin.o 	nrprd.o 	nrstp.o \
		nrloren.o 	nrfermi.o 	nrnuc.o 	nrprton.o 	nrprbin.o \
		nr2iso.o 	nrhis.o 	nrrodis.o 	nrranve.o \
		nrcst.o 	nrctpp.o 	nrctpnl.o 	nrctpnh.o \
		wgproton.o	vcmkvcbf.o  nrintr.o    nrsetxpos.o \
		nrcard.o	nrlibver.o

INCS = nrnuclparam.h nrcard.h nucleonfsihist.h 

all: lib${LIBNAME}.a

lib${LIBNAME}.a: ${OBJS}
	${AR}     rv $@ ${OBJS}
	${RANLIB} $@

install.include: ${INCS}
	for i in $? ;do ${FH2H} $$i > `echo $$i | sed "s/\.h$$/C.h/"`; done
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS} ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS:.h=C.h} ${NEUT_ROOT}/include

install.library: lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib


