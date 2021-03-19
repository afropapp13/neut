include $(top_srcdir)/build.common

AM_FFLAGS = -DNENUCC_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@nuccorspl.@NEUT_VERSION_NUM@.a

if NEUT_LIBCOMPAT
lib_LIBRARIES = $(LIBNAME)
else
noinst_LIBRARIES = $(LIBNAME)
endif

lib@LIBPREFIX@nuccorspl_@NEUT_VERSION_NUM@_a_SOURCES = nr2iso.F nrcst.F nrhis.F nrllor.F nrprton.F nrstp.F nrazio.F nrctpnh.F nrinit.F nrloren.F nrranve.F nrcard.F nrctpnl.F nrintr.F nrnuc.F nrrodis.F nrcfd.F nrctpp.F nrkin.F nrprbin.F nrroxy.F wgproton.F nrcone.F nrfermi.F nrlibver.F nrprd.F

if NEUT_USE_ZBS

lib@LIBPREFIX@nuccorspl_@NEUT_VERSION_NUM@_a_SOURCES += vcmkvcbf.F nrsetxpos.F

endif

# Unused sources that live in the directory
UNUSED_srcs = proton_nuceff.F

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
