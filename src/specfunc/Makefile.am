include $(top_srcdir)/build.common

AM_FFLAGS = -DSPECFUNCVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

# bin_PROGRAMS = makeEffSFTables makeTables

# makeEffSFTables_SOURCES = makeEffSFTables.F
# makeEffSFTables_LDADD = libNEUTspecfunc.a

# makeTables_SOURCES = makeTables.F
# makeTables_LDADD = libNEUTspecfunc.a

lib_LIBRARIES = libNEUTspecfunc.a

libNEUTspecfunc_a_SOURCES = buildEffSf.F buildSf.F calcLH.F combineEffSFTables.F combineTables.F effsfevent.F effSFMaxDiff.F fourVector.F maxDiff.F readEffSFXsecData.F readXsecData.F selectEffSFValues.F selectSfValues.F sfevent.F tauvect_dummy.F

# Unused sources that live in the directory
UNUSED_srcs = selectSfValues.old.F
