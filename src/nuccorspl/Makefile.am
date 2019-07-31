include $(top_srcdir)/build.common

lib_LIBRARIES = libNEUTnuccorspl.a

libNEUTnuccorspl_a_FFLAGS = -DNENUCC_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

libNEUTnuccorspl_a_SOURCES = nr2iso.F nrcst.F nrhis.F nrllor.F nrprton.F nrstp.F nrazio.F nrctpnh.F nrinit.F nrloren.F nrranve.F nrcard.F nrctpnl.F nrintr.F nrnuc.F nrrodis.F vcmkvcbf.F nrcfd.F nrctpp.F nrkin.F nrprbin.F nrroxy.F wgproton.F nrcone.F nrfermi.F nrlibver.F nrprd.F nrsetxpos.F

# Unused sources that live in the directory
UNUSED_srcs = proton_nuceff.F
