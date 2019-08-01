include $(top_srcdir)/build.common

AM_FFLAGS = -DNEUTSMPL_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DNEUTSMPL_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

lib_LIBRARIES = libNEUTneutsmpl.a

libNEUTneutsmpl_a_SOURCES = nevecgen.F grndmq.F rmarin_dum.F flx_intrfc.cc ufm2th1dist.cc fill_dsdq2.cc trapfpe.c

libNEUTneutsmpl_a_CXXFLAGS = -DNEUTSMPL_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) -I$(top_srcdir)/neutclass

bin_PROGRAMS = neutroot2

neutroot2_SOURCES = neutroot.F

neutroot2_DEPENDENCIES = libNEUTneutsmpl.a

neutroot2_LDFLAGS = $(EXE_LINKFLAGS) libNEUTneutsmpl.a  $(NEUT_STATIC_LIBS) $(top_builddir)/neutclass/libNEUTneutclass.la $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS) -lstdc++
