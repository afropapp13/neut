#Its not ideal, but this should be manually kept aprised of any changes in build.common

include $(top_srcdir)/build.common

INST_ROOT = $(prefix)

NEUT_VERSION = @NEUT_VERSION@

HEPMC_FLAGS =

if NEUT_HEPMC
HEPMC_FLAGS += -L$(HEPMC_LIBDIR) -L$(NUHEPMC_LIBDIR) -lNuHepMC -lHepMC3 -lHepMC3rootIO
endif

ROOT_FLAGS = $(HEPMC_FLAGS) $(ROOT_LINKFLAGS) $(ROOT_LIBS)
CERN_FLAGS = $(CERN_LINKFLAGS) $(CERN_LIBS)

if NEUT_LIBCOMPAT

IO_LIB_LIST = neutrootTreeSingleton.so neutctrl.so neutfsipart.so neutfsivert.so neutnucfsistep.so neutnucfsivert.so neutpart.so neutvect.so neutvtx.so
IO_LIB_FLAGS = $(INST_ROOT)/lib/neutrootTreeSingleton.so $(INST_ROOT)/lib/neutctrl.so $(INST_ROOT)/lib/neutfsipart.so $(INST_ROOT)/lib/neutfsivert.so $(INST_ROOT)/lib/neutnucfsistep.so $(INST_ROOT)/lib/neutnucfsivert.so $(INST_ROOT)/lib/neutpart.so $(INST_ROOT)/lib/neutvect.so $(INST_ROOT)/lib/neutvtx.so
RW_LIB = libNReWeight.so
RW_FLAGS = -lNReWeight
LIB_LIST = $(RW_LIB) $(IO_LIB_LIST) libneutcore.a libnuccorspl.a libnuceff.a libn1p1h.a libht2p2h.a libradcorr.a libtauola.a libspecfunc.a libskmcsvc.a libzbsfns.a libpartnuck.a
LIB_FLAGS = $(RW_FLAGS) $(IO_LIB_FLAGS) -lneutcore -lnuccorspl -lnuceff -ln1p1h -lht2p2h -lradcorr -ltauola -lspecfunc -lskmcsvc -lzbsfns -lpartnuck $(ROOT_FLAGS) $(CERN_FLAGS)

else

IO_LIB_LIST = libNEUTClassUtils.a libNEUTClass.so
JUST_IO_LIB_FLAGS = -lNEUTClassUtils -lNEUTClass
IO_LIB_FLAGS = $(JUST_IO_LIB_FLAGS) $(ROOT_FLAGS)
RW_LIB = libNEUTReWeight.so
RW_FLAGS = -lNEUTReWeight
LIB_LIST = $(RW_LIB) $(IO_LIB_LIST) libNEUT.a
LIB_FLAGS = $(RW_FLAGS) $(JUST_IO_LIB_FLAGS) -lNEUT $(ROOT_FLAGS) $(CERN_FLAGS)

endif

neut-config: neut-config.in
	cat ${top_srcdir}/tools/neut-config.in | sed "s|__INST_ROOT__|$(INST_ROOT)|g" | sed "s|__LIB_LIST__|$(LIB_LIST)|g" | sed "s|__IO_LIB_LIST__|$(IO_LIB_LIST)|g" | sed "s|__RW_LIB__|$(RW_LIB)|g" | sed "s|__LIB_FLAGS__|$(LIB_FLAGS)|g" | sed "s|__IO_LIB_FLAGS__|$(IO_LIB_FLAGS)|g" | sed "s|__RW_FLAGS__|$(RW_FLAGS)|g" | sed "s|__NEUT_VERSION__|$(NEUT_VERSION)|g" | sed "s|__ROOT_FLAGS__|$(ROOT_FLAGS)|g" | sed "s|__CERN_FLAGS__|$(CERN_FLAGS)|g" > neut-config
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
