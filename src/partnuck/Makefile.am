include $(top_srcdir)/build.common

AM_FFLAGS = -DPARTNUCKVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@partnuck.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

lib@LIBPREFIX@partnuck_@NEUT_VERSION_NUM@_a_SOURCES = cfd.F cone.F lloren.F loren.F partnuc.F prbin.F ranf.F ranve.F

# Unused sources that live in the directory
UNUSED_srcs =
