include $(top_srcdir)/build.common

AM_FFLAGS = -DSPECFUNCVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

# bin_PROGRAMS = makeEffSFTables makeTables

# makeEffSFTables_SOURCES = makeEffSFTables.F
# makeEffSFTables_LDADD = lib@LIBPREFIX@specfunc.a

# makeTables_SOURCES = makeTables.F
# makeTables_LDADD = lib@LIBPREFIX@specfunc.a

LIBNAME = lib@LIBPREFIX@specfunc.@NEUT_VERSION_NUM@.a

if NEUT_LIBCOMPAT
lib_LIBRARIES = $(LIBNAME)
else
noinst_LIBRARIES = $(LIBNAME)
endif

lib@LIBPREFIX@specfunc_@NEUT_VERSION_NUM@_a_SOURCES = buildEffSf.F buildSf.F calcLH.F combineEffSFTables.F combineTables.F effsfevent.F effSFMaxDiff.F maxDiff.F readEffSFXsecData.F readXsecData.F selectEffSFValues.F selectSfValues.F sfevent.F

# Unused sources that live in the directory
UNUSED_srcs = selectSfValues.old.F tauvect_dummy.F

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
