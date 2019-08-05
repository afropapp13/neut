include $(top_srcdir)/build.common

AM_FFLAGS = -DSKMCSVCVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@skmcsvc.@NEUT_VERSION_MAJ@.@NEUT_VERSION_MIN@.@NEUT_VERSION_REV@.a

lib_LIBRARIES = $(LIBNAME)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

lib@LIBPREFIX@skmcsvc_@NEUT_VERSION_MAJ@_@NEUT_VERSION_MIN@_@NEUT_VERSION_REV@_a_SOURCES = excard.F getuni.F mcmass.F mcmassgv.F mcvecbst.F necardev.F nesetpos.F rlu.F rndir.F rnferm.F rnpwav.F rnswav.F vcclcm.F vcphsvct.F vcppwave.F vcpswave.F vcrambo.F vcwtrmbo.F

# Unused sources that live in the directory
UNUSED_srcs =
