## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

FCFLAGS += -DZBSFNSVER=${IZBSFNSVER}
NEUTINCDIR += -I${SKOFL_ROOT}/inc -I${NEUT_ROOT}/src/t2kflux_zbs

LIBBASE	= zbsfns
LIBNAME = $(LIBBASE)_$(ZBSFNSVER)

OBJS = nemknebk.o  nemkmodelbk.o vcmkvc.o vcmkvx.o \
       nemknetarg.o nerdnetarg.o nemkfsibk.o nemkcrsbk.o \
       nemknucfsibk.o nerdnucfsibk.o \
       nerdnebk.o vcrdvccm.o vcrdvc.o vcrdvx.o \
       nerdnetarg.o nerdfsibk.o nerdcrsbk.o vcmkdetgm.o \
       vcclfsicm.o 

all: lib${LIBNAME}.a

lib${LIBNAME}.a: ${OBJS}
	${AR}     rv $@ ${OBJS}
	${RANLIB} $@

install.library: lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib


