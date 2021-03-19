include $(top_srcdir)/build.common

AM_CXXFLAGS = -DRADCORRVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) -std=c++11

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
	-(mkdir -p ${bindir}/../share/neut/crsdat/radcorr; cp -r $(top_srcdir)/radcorr/*.root ${bindir}/../share/neut/crsdat/radcorr )

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

else

install-exec-hook:
	-(mkdir -p ${bindir}/../share/neut/crsdat/radcorr; cp -r $(top_srcdir)/radcorr/*.root ${bindir}/../share/neut/crsdat/radcorr )

endif
