include $(top_srcdir)/build.common

AM_CXXFLAGS = -DRADCORRVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

LIBNAME = lib@LIBPREFIX@radcorr.@NEUT_VERSION_NUM@.a

if NEUT_LIBCOMPAT
lib_LIBRARIES = $(LIBNAME)
else
noinst_LIBRARIES = $(LIBNAME)
endif

lib@LIBPREFIX@radcorr_@NEUT_VERSION_NUM@_a_SOURCES = radcorr.cc

# Unused sources that live in the directory
UNUSED_srcs =

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
