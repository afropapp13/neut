## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

FCFLAGS += -DTAUOLAVER=${ITAUOLAVER}

LIBBASE	= tauola
LIBNAME = $(LIBBASE)_$(NEUTCOREVER)

OBJ =	tauola.o formf.o choice.o dcdmas.o lunpik.o tralo4.o pkorb.o \
        tainicom.o tainidcy.o tainimas.o tainipar.o tapricom.o \
        taudecay.o tauvect.o tacard.o  \
		filhep.o taurdf.o \
        nemktaubk.o taupol.o mrst2002.o

OBJS = $(OBJ)

INCS = taucom.h taumc.h

all: lib${LIBNAME}.a

lib${LIBNAME}.a: ${OBJS}
	${AR}     rv $@ ${OBJS}
	${RANLIB} $@

install.include: ${INCS}
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS} ${NEUT_ROOT}/include

install.library: lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib




