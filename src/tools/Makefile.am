#Its not ideal, but this should be manually kept aprised of any changes in build.common

include $(top_srcdir)/build.common

INST_ROOT = $(prefix)

NEUT_VERSION = @NEUT_VERSION@

HEPMC_FLAGS =

FEATURES = reweight common-iface

if NEUT_HEPMC
HEPMC_FLAGS += -L$(HEPMC_LIBDIR) -L$(NUHEPMC_LIBDIR) -lNuHepMC -lHepMC3 -lHepMC3rootIO
FEATURES += HepMC3
endif

if NEUT_USE_NEUTGEOM
FEATURES += neutgeom
endif

if NEUT_WRITE_NUCFSI
FEATURES += neutvect-nucfsistep
endif

if NEUT_DEBUG
FEATURES += debug
endif

ROOT_FLAGS = $(HEPMC_FLAGS) $(ROOT_LINKFLAGS) $(ROOT_LIBS)
NEUT_CFLAGS = -fPIC $(ROOT_CFLAGS)
CERN_FLAGS = $(CERN_LINKFLAGS) $(CERN_LIBS)

IO_LIB_FLAGS = -lNEUTClassUtils -lNEUTClass
IO_LIB_FLAGS_PRIVATE = $(HEPMC_FLAGS) $(ROOT_FLAGS)

LIB_FLAGS = -lNEUTReWeight -lNEUT $(IO_LIB_FLAGS)
LIB_FLAGS_PRIVATE = $(IO_LIB_FLAGS_PRIVATE) $(CERN_FLAGS)

NEUT.pc: NEUT.pc.in
	cat ${top_srcdir}/tools/NEUT.pc.in \
	| sed "s|__INST_ROOT__|$(INST_ROOT)|g" \
	| sed "s|__FEATURES__|$(FEATURES)|g" \
	| sed "s|__ROOT_FLAGS__|$(ROOT_FLAGS)|g" \
	| sed "s|__NEUT_CFLAGS__|$(NEUT_CFLAGS)|g" \
	| sed "s|__CERN_FLAGS__|$(CERN_FLAGS)|g" \
	| sed "s|__IO_LIB_FLAGS__|$(IO_LIB_FLAGS)|g" \
	| sed "s|__NEUT_VERSION__|$(NEUT_VERSION)|g" \
	| sed "s|__EXE_LINKFLAGS__|$(EXE_LINKFLAGS)|g" \
	| sed "s|__LIB_FLAGS__|$(LIB_FLAGS)|g" \
	| sed "s|__LIB_FLAGS_PRIVATE__|$(LIB_FLAGS_PRIVATE)|g" > NEUT.pc

setup.sh: setup.sh.in
	cat ${top_srcdir}/tools/setup.sh.in \
	| sed "s|__ROOT_VERSION_COMPILE__|$(ROOT_VERSION_COMPILE)|g" \
	| sed "s|__ROOT_ROOT_COMPILE__|$(ROOT_ROOT_COMPILE)|g" \
	| sed "s|__USE_HEPMC__|$(USE_HEPMC)|g" \
	| sed "s|__HEPMC_LIBDIR__|$(HEPMC_LIBDIR)|g" > setup.sh

install-exec-hook: neut-config setup.sh NEUT.pc
	-cp ${top_srcdir}/tools/neut-config ${bindir}/
	-mkdir -p ${prefix}/cmake
	-cp ${top_srcdir}/tools/NEUTConfig.cmake ${prefix}/cmake
	chmod +x ${bindir}/neut-config
	-cp setup.sh NEUT.pc ${prefix}/

uninstall-hook:
	-(cd ${bindir}; rm -f neut-config)
	-(cd ${prefix}; rm -f setup.sh NEUT.pc)
