include $(top_srcdir)/build.common

LIBNAME = libn1p1h.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

N1p1hDEFINES = -DWITH_NEUT

AM_CXXFLAGS = -DN1P1HVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) $(N1p1hDEFINES)

# src/qe.F will fail to build with the standard NEUT suite of flags so they are purposefully not used here.
AM_FFLAGS = -DN1P1HVER=@NEUT_VERSION_NUM@ $(N1p1hDEFINES)

libn1p1h_@NEUT_VERSION_NUM@_a_SOURCES = src/N1p1h2d.C src/N1p1h.C src/N1p1hFortranInterface.C src/Nucleus.C src/qe.F

# Unused sources that live in the directory
UNUSED_srcs = src/Simplex.C
