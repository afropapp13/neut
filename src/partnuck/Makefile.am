include $(top_srcdir)/build.common

AM_FFLAGS = -DPARTNUCKVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@partnuck.@NEUT_VERSION_NUM@.a

if NEUT_LIBCOMPAT
lib_LIBRARIES = $(LIBNAME)
else
noinst_LIBRARIES = $(LIBNAME)
endif

lib@LIBPREFIX@partnuck_@NEUT_VERSION_NUM@_a_SOURCES = cfd.F cone.F lloren.F loren.F partnuc.F prbin.F ranf.F ranve.F

# Unused sources that live in the directory
UNUSED_srcs =

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
