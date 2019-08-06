include $(top_srcdir)/build.common

AM_FFLAGS = -DSPECFUNCVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@specfunc.@NEUT_VERSION_NUM@.a

if NEUT_LIBCOMPAT
lib_LIBRARIES = $(LIBNAME)
else
noinst_LIBRARIES = $(LIBNAME)
endif

lib@LIBPREFIX@specfunc_@NEUT_VERSION_NUM@_a_SOURCES = buildEffSf.F buildSf.F calcLH.F combineEffSFTables.F combineTables.F effsfevent.F effSFMaxDiff.F maxDiff.F readEffSFXsecData.F readXsecData.F selectEffSFValues.F selectSfValues.F sfevent.F fourVector.F

# Unused sources that live in the directory
UNUSED_srcs =

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
