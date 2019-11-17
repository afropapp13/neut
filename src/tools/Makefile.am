#Its not ideal, but this should be manually kept aprised of any changes in build.common

include $(top_srcdir)/build.common

INST_ROOT = $(prefix)

if NEUT_LIBCOMPAT

IO_LIB_LIST = neutrootTreeSingleton.so neutctrl.so neutfsipart.so neutfsivert.so neutnucfsistep.so neutnucfsivert.so neutpart.so neutvect.so neutvtx.so
IO_LIB_FLAGS = $(INST_ROOT)/lib/neutrootTreeSingleton.so $(INST_ROOT)/lib/neutctrl.so $(INST_ROOT)/lib/neutfsipart.so $(INST_ROOT)/lib/neutfsivert.so $(INST_ROOT)/lib/neutnucfsistep.so $(INST_ROOT)/lib/neutnucfsivert.so $(INST_ROOT)/lib/neutpart.so $(INST_ROOT)/lib/neutvect.so $(INST_ROOT)/lib/neutvtx.so
LIB_LIST = $(IO_LIB_LIST) libneutcore.a libnuccorspl.a libnuceff.a libn1p1h.a libht2p2h.a libradcorr.a libtauola.a libspecfunc.a libskmcsvc.a libzbsfns.a libpartnuck.a
LIB_FLAGS = $(IO_LIB_FLAGS) -lneutcore -lnuccorspl -lnuceff -ln1p1h -lht2p2h -lradcorr -ltauola -lspecfunc -lskmcsvc -lzbsfns -lpartnuck
RW_LIB = libNReWeight.so
RW_FLAGS = -lNReWeight

else

IO_LIB_LIST = libNEUTClassUtils.a libNEUTClass.so
IO_LIB_FLAGS = -lNEUTClassUtils -lNEUTClass
LIB_LIST = libNEUT.a libNEUTn1p1h.a libNEUTht2p2h.a
LIB_FLAGS = -lNEUT -lNEUTn1p1h -lNEUTht2p2h
RW_LIB = libNEUTReWeight.so
RW_FLAGS = -lNEUTReWeight
endif

NEUT_VERSION = @NEUT_VERSION@

ROOT_FLAGS = $(ROOT_LINKFLAGS) $(ROOT_LIBS)
CERN_FLAGS = $(CERN_LINKFLAGS) $(CERN_LIBS)

neut-config: neut-config.in
	cat ${top_srcdir}/tools/neut-config.in | sed "s|__INST_ROOT__|$(INST_ROOT)|g" | sed "s|__LIB_LIST__|$(LIB_LIST)|g" | sed "s|__IO_LIB_LIST__|$(IO_LIB_LIST)|g" | sed "s|__RW_LIB__|$(RW_LIB)|g" | sed "s|__LIB_FLAGS__|$(LIB_FLAGS)|g" | sed "s|__IO_LIB_FLAGS__|$(IO_LIB_FLAGS)|g" | sed "s|__RW_FLAGS__|$(RW_FLAGS)|g" | sed "s|__NEUT_VERSION__|$(NEUT_VERSION)|g" | sed "s|__ROOT_FLAGS__|$(ROOT_FLAGS)|g" | sed "s|__CERN_FLAGS__|$(CERN_FLAGS)|g" > neut-config
	chmod +x neut-config

install-exec-hook: neut-config
	-cp neut-config ${bindir}/
	-cp ${top_srcdir}/tools/setup.sh ${prefix}/
	chmod +x ${prefix}/setup.sh

uninstall-hook:
	-(cd ${bindir}; rm -f neut-config)
	-(cd ${prefix}; rm -f setup.sh)
