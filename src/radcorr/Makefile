## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

LIBBASE	= radcorr
LIBNAME = $(LIBBASE)_$(RADCORRVER)

OBJS  = radcorr.o

INCS = hadron.h masses.h nucnewprd.h partnucint.h ran_ranlux.h

all: lib${LIBNAME}.a

install.include:
	echo "nothing to be done."

lib${LIBNAME}.a: ${OBJS}
	${AR}     rv $@ ${OBJS}
	${RANLIB} $@

install.library: lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib

