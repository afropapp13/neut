include $(top_srcdir)/build.common

AM_FFLAGS = -DSKMCSVCVER=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DSKMCSVC_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

LIBNAME = libskmcsvc.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

libskmcsvc_@NEUT_VERSION_NUM@_a_SOURCES = excard.F getuni.F mcmass.F mcmassgv.F mcvecbst.F necardev.F necardev_toml.cc nesetpos.F rlu.F rndir.F rnferm.F rnpwav.F rnswav.F vcclcm.F vcphsvct.F vcppwave.F vcpswave.F vcrambo.F vcwtrmbo.F

# Unused sources that live in the directory
UNUSED_srcs =
