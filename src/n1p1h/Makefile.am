include $(top_srcdir)/build.common

LIBNAME = lib@LIBPREFIX@n1p1h.@NEUT_VERSION_NUM@.a

lib_LIBRARIES = $(LIBNAME)

N1p1hDEFINES = -DWITH_NEUT

AM_CXXFLAGS = -DN1P1HVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) $(N1p1hDEFINES)

# src/qe.F will fail to build with the standard NEUT suite of flags so they are purposefully not used here.
AM_FFLAGS = -DN1P1HVER=@NEUT_VERSION_NUM@ $(N1p1hDEFINES)

lib@LIBPREFIX@n1p1h_@NEUT_VERSION_NUM@_a_SOURCES = src/N1p1h2d.C src/N1p1h.C src/N1p1hFortranInterface.C src/Nucleus.C src/qe.F

# Unused sources that live in the directory
UNUSED_srcs = src/Simplex.C

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
