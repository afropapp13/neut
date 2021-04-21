include $(top_srcdir)/build.common

nieves2p2hC.h: inc/nieves2p2h.h
	$(top_srcdir)/tools/fh2h.pl $< > $@

BUILT_SOURCES = nieves2p2hC.h
CLEANFILES = $(BUILT_SOURCES)

include_HEADERS = $(BUILT_SOURCES)

LIBNAME = libht2p2h.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

AM_CXXFLAGS = -DHT2P2HVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

libht2p2h_@NEUT_VERSION_NUM@_a_SOURCES = src/HadronTensor.C src/HT2p2h.C src/HT2p2hFotranInterface.C src/Nucleus2p2h.C

# Unused sources that live in the directory
UNUSED_srcs =
