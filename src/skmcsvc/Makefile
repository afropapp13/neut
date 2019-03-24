## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

LIBBASE	= skmcsvc
LIBNAME = ${LIBBASE}_$(SKMCSVCVER)

## physics related fucntions taken from SK
POBJ =	mcvecbst.o vcphsvct.o vcwtrmbo.o vcrambo.o mcmass.o mcmassgv.o

## random number functions taken from SK 
ROBJ =	nesetpos.o rlu.o rndir.o \
		rnferm.o rnswav.o rnpwav.o vcpswave.o vcppwave.o

## card related functions taken from SK 
MOBJ =	getuni.o excard.o

## common handing functions taken from SK 
COBJ =	necardev.o vcclcm.o

OBJS = $(POBJ) $(IOBJ) $(ROBJ) $(MOBJ) $(COBJ)

INCS = vcvrtx.h vcwork.h necardev.h

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


