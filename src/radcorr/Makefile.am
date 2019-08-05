include $(top_srcdir)/build.common

AM_CXXFLAGS = -DRADCORRVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

LIBNAME = lib@LIBPREFIX@radcorr.@NEUT_VERSION_MAJ@.@NEUT_VERSION_MIN@.@NEUT_VERSION_REV@.a

lib_LIBRARIES = $(LIBNAME)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

lib@LIBPREFIX@radcorr_@NEUT_VERSION_MAJ@_@NEUT_VERSION_MIN@_@NEUT_VERSION_REV@_a_SOURCES = radcorr.cc

# Unused sources that live in the directory
UNUSED_srcs =
