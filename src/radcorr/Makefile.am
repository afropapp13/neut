include $(top_srcdir)/build.common

AM_CXXFLAGS = -DRADCORRVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

LIBNAME = lib@LIBPREFIX@radcorr.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

lib@LIBPREFIX@radcorr_@NEUT_VERSION_NUM@_a_SOURCES = radcorr.cc

# Unused sources that live in the directory
UNUSED_srcs =
