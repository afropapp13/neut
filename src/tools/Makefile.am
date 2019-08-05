#Its not ideal, but this should be manually kept aprised of any changes in build.common

include $(top_srcdir)/build.common

INST_ROOT = $(prefix)

if NEUT_LIBCOMPAT
LIB_LIST = libneutrootTreeSingleton.so libneutctrl.so libneutfsipart.so libneutfsivert.so libneutnucfsistep.so libneutnucfsivert.so libneutpart.so libneutvect.so libneutvtx.so
LIB_FLAGS = -lneutrootTreeSingleton.so -lneutctrl.so -lneutfsipart.so -lneutfsivert.so -lneutnucfsistep.so -lneutnucfsivert.so -lneutpart.so -lneutvect.so -lneutvtx.so
else
LIB_LIST = lib@LIBPREFIX@neutclassUtils.a lib@LIBPREFIX@neutclass.a
LIB_FLAGS = -l@LIBPREFIX@neutclassUtils -l@LIBPREFIX@neutclass
endif

LIB_LIST += lib@LIBPREFIX@neutcore.a lib@LIBPREFIX@nuccorspl.a lib@LIBPREFIX@nuceff.a lib@LIBPREFIX@n1p1h.a lib@LIBPREFIX@ht2p2h.a lib@LIBPREFIX@radcorr.a lib@LIBPREFIX@tauola.a lib@LIBPREFIX@specfunc.a lib@LIBPREFIX@skmcsvc.a lib@LIBPREFIX@zbsfns.a lib@LIBPREFIX@partnuck.a

RW_LIB = libNEUTReWeight.so
LIB_FLAGS += -l@LIBPREFIX@neutcore -l@LIBPREFIX@nuccorspl -l@LIBPREFIX@nuceff -l@LIBPREFIX@n1p1h -l@LIBPREFIX@ht2p2h -l@LIBPREFIX@radcorr -l@LIBPREFIX@tauola -l@LIBPREFIX@specfunc -l@LIBPREFIX@skmcsvc -l@LIBPREFIX@zbsfns -l@LIBPREFIX@partnuck

RW_FLAGS = -lNEUTReWeight
NEUT_VERSION = @NEUT_VERSION@

ROOT_FLAGS = $(ROOT_LINKFLAGS) $(ROOT_LIBS)
CERN_FLAGS = $(CERN_LINKFLAGS) $(CERN_LIBS)

neut-config: neut-config.in
	cat ${top_srcdir}/tools/neut-config.in | sed "s|__INST_ROOT__|$(INST_ROOT)|g" | sed "s|__LIB_LIST__|$(LIB_LIST)|g" | sed "s|__RW_LIB__|$(RW_LIB)|g" | sed "s|__LIB_FLAGS__|$(LIB_FLAGS)|g" | sed "s|__RW_FLAGS__|$(RW_FLAGS)|g" | sed "s|__NEUT_VERSION__|$(NEUT_VERSION)|g" | sed "s|__ROOT_FLAGS__|$(ROOT_FLAGS)|g" | sed "s|__CERN_FLAGS__|$(CERN_FLAGS)|g" > neut-config
	chmod +x neut-config

install-exec-hook: neut-config
	-cp neut-config ${bindir}/

uninstall-hook:
	-(cd ${bindir}; rm -f neut-config)
