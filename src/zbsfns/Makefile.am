include $(top_srcdir)/build.common

AM_FFLAGS = -DZBSFNSVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = libzbsfns.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

libzbsfns_@NEUT_VERSION_NUM@_a_SOURCES = nemkcrsbk.F nemkfsibk.F nemkmodelbk.F nemknebk.F nemknetarg.F nemknucfsibk.F nemktaubk.F nerdcrsbk.F nerdfsibk.F nerdnebk.F nerdnetarg.F nerdnucfsibk.F vcclfsicm.F vcmkdetgm.F vcmkfsicm.F vcmkfsip.F vcmkfsiv.F vcmkvc.F vcmkvx.F vcrdfsicm.F vcrdfsip.F vcrdfsiv.F vcrdvccm.F vcrdvc.F vcrdvx.F

# Unused sources that live in the directory
UNUSED_srcs =