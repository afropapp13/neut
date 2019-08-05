include $(top_srcdir)/build.common

LIBNAME = lib@LIBPREFIX@ht2p2h.@NEUT_VERSION_MAJ@.@NEUT_VERSION_MIN@.@NEUT_VERSION_REV@.a

lib_LIBRARIES = $(LIBNAME)

AM_CXXFLAGS = -DHT2P2HVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

lib@LIBPREFIX@ht2p2h_@NEUT_VERSION_MAJ@_@NEUT_VERSION_MIN@_@NEUT_VERSION_REV@_a_SOURCES = src/HadronTensor.C src/HT2p2h.C src/HT2p2hFotranInterface.C src/Nucleus2p2h.C


# Unused sources that live in the directory
UNUSED_srcs =

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)
