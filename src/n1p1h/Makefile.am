include $(top_srcdir)/build.common

lib_LIBRARIES = libNEUTn1p1h.a

libNEUTn1p1h_a_CXXFLAGS = -DN1P1HVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

# src/qe.F will fail to build with the standard NEUT suite of flags so they are purposefully not used here.
libNEUTn1p1h_a_FFLAGS = -DN1P1HVER=@NEUT_VERSION_NUM@

libNEUTn1p1h_a_SOURCES = src/N1p1h2d.C src/N1p1h.C src/N1p1hFortranInterface.C src/Nucleus.C src/qe.F src/Simplex.C


# Unused sources that live in the directory
UNUSED_srcs =
