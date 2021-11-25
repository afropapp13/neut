include $(top_srcdir)/build.common

AM_FFLAGS = -DNENUCC_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DNENUCC_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

LIBNAME = libnuccorspl.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

libnuccorspl_@NEUT_VERSION_NUM@_a_SOURCES = nr2iso.F nrcst.F nrhis.F nrllor.F nrprton.F \
											nrstp.F nrazio.F nrctpnh.F nrinit.F nrloren.F \
											nrranve.F nrcard.F nrcard_toml.cc nrctpnl.F \
											nrintr.F nrnuc.F nrrodis.F nrcfd.F nrctpp.F \
											nrkin.F nrprbin.F nrroxy.F wgproton.F nrcone.F \
											nrfermi.F nrlibver.F nrprd.F vcmkvcbf.F \
											nrsetxpos.F


# Unused sources that live in the directory
UNUSED_srcs = proton_nuceff.F
