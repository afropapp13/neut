if !NEUT_LIBCOMPAT

include $(top_srcdir)/build.common

AM_FFLAGS = -DNEMEC_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@neutmec.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

lib@LIBPREFIX@neutmec_@NEUT_VERSION_NUM@_a_SOURCES = nemechad.F nemeclpv.F nemecvcp.F neranbll2.F fnmec_nieves.F fnmec_select.F

include_HEADERS = meccrs.h

# Unused sources that live in the directory
UNUSED_srcs =

endif
