include $(top_srcdir)/build.common

LIBNAME = lib@LIBPREFIX@ht2p2h.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

AM_CXXFLAGS = -DHT2P2HVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

lib@LIBPREFIX@ht2p2h_@NEUT_VERSION_NUM@_a_SOURCES = src/HadronTensor.C src/HT2p2h.C src/HT2p2hFotranInterface.C src/Nucleus2p2h.C

# Unused sources that live in the directory
UNUSED_srcs =
