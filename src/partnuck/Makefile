## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

LIBBASE	= partnuck
LIBNAME = $(LIBBASE)_$(PARTNUCKVER)

OBJS  = partnuc.o ranve.o cfd.o lloren.o loren.o prbin.o cone.o \
		ranf.o

INCS = hadron.h masses.h nucnewprd.h partnucint.h ran_ranlux.h

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

