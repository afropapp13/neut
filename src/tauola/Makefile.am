include $(top_srcdir)/build.common

AM_FFLAGS = -DTAUOLAVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@tauola.@NEUT_VERSION_MAJ@.@NEUT_VERSION_MIN@.@NEUT_VERSION_REV@.a

lib_LIBRARIES = $(LIBNAME)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

lib@LIBPREFIX@tauola_@NEUT_VERSION_MAJ@_@NEUT_VERSION_MIN@_@NEUT_VERSION_REV@_a_SOURCES = choice.F dcdmas.F filhep.F formf.F lunpik.F pkorb.F tacard.F tainicom.F tainidcy.F tainimas.F tainipar.F tapricom.F taudecay.F tauola.F taurdf.F tauvect.F tralo4.F

# Unused sources that live in the directory
UNUSED_srcs =
