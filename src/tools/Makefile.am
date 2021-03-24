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

ROOT_FLAGS = $(HEPMC_FLAGS) $(ROOT_LINKFLAGS) $(ROOT_LIBS)
CERN_FLAGS = $(CERN_LINKFLAGS) $(CERN_LIBS)

IO_LIB_LIST = libNEUTClassUtils.a libNEUTClass.so
JUST_IO_LIB_FLAGS = -lNEUTClassUtils -lNEUTClass
IO_LIB_FLAGS = $(JUST_IO_LIB_FLAGS) $(ROOT_FLAGS)
RW_LIB = libNEUTReWeight.so
RW_FLAGS = -lNEUTReWeight
LIB_LIST = $(RW_LIB) $(IO_LIB_LIST) libNEUT.a
LIB_FLAGS = $(RW_FLAGS) $(JUST_IO_LIB_FLAGS) -lNEUT $(ROOT_FLAGS) $(CERN_FLAGS)

neut-config: neut-config.in
	cat ${top_srcdir}/tools/neut-config.in | sed "s|__FEATURES__|$(FEATURES)|g" | sed "s|__INST_ROOT__|$(INST_ROOT)|g" | sed "s|__LIB_LIST__|$(LIB_LIST)|g" | sed "s|__IO_LIB_LIST__|$(IO_LIB_LIST)|g" | sed "s|__RW_LIB__|$(RW_LIB)|g" | sed "s|__LIB_FLAGS__|$(LIB_FLAGS)|g" | sed "s|__IO_LIB_FLAGS__|$(IO_LIB_FLAGS)|g" | sed "s|__RW_FLAGS__|$(RW_FLAGS)|g" | sed "s|__NEUT_VERSION__|$(NEUT_VERSION)|g" | sed "s|__ROOT_FLAGS__|$(ROOT_FLAGS)|g" | sed "s|__CERN_FLAGS__|$(CERN_FLAGS)|g" > neut-config
	chmod +x neut-config

setup.sh: setup.sh.in
	cat ${top_srcdir}/tools/setup.sh.in | sed "s|__USE_HEPMC__|$(USE_HEPMC)|g" | sed "s|__HEPMC_LIBDIR__|$(HEPMC_LIBDIR)|g" > setup.sh

install-exec-hook: neut-config setup.sh
	-cp neut-config ${bindir}/
	-cp setup.sh ${prefix}/
	chmod +x ${prefix}/setup.sh

uninstall-hook:
	-(cd ${bindir}; rm -f neut-config)
	-(cd ${prefix}; rm -f setup.sh)
