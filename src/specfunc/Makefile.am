include $(top_srcdir)/build.common

AM_FFLAGS = -DSPECFUNCVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DSPECFUNCVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

LIBNAME = libspecfunc.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

libspecfunc_@NEUT_VERSION_NUM@_a_SOURCES = buildEffSf.F buildSf.F calcLH.F combineEffSFTables.F combineTables.F effsfevent.F effSFMaxDiff.F maxDiff.F readEffSFXsecData.F readXsecData.F selectEffSFValues.F selectSfValues.F sfevent.F fourVector.F HaveSpecFuncForNuclei.cc

# Unused sources that live in the directory
UNUSED_srcs =

