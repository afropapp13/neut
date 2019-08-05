include $(top_srcdir)/build.common

AM_FFLAGS = -DNENUCC_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@nuccorspl.@NEUT_VERSION_MAJ@.@NEUT_VERSION_MIN@.@NEUT_VERSION_REV@.a

lib_LIBRARIES = $(LIBNAME)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

lib@LIBPREFIX@nuccorspl_@NEUT_VERSION_MAJ@_@NEUT_VERSION_MIN@_@NEUT_VERSION_REV@_a_SOURCES = nr2iso.F nrcst.F nrhis.F nrllor.F nrprton.F nrstp.F nrazio.F nrctpnh.F nrinit.F nrloren.F nrranve.F nrcard.F nrctpnl.F nrintr.F nrnuc.F nrrodis.F vcmkvcbf.F nrcfd.F nrctpp.F nrkin.F nrprbin.F nrroxy.F wgproton.F nrcone.F nrfermi.F nrlibver.F nrprd.F nrsetxpos.F

# Unused sources that live in the directory
UNUSED_srcs = proton_nuceff.F
