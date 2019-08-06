include $(top_srcdir)/build.common

AM_FFLAGS = $(DEFAULT_FFLAGS)
AM_CXXFLAGS = $(DEFAULT_CXXFLAGS)

bin_PROGRAMS = calc_ccqe_crs

calc_ccqe_crs_SOURCES = clenshaw_curtis.cc CrossSection.cc Initialize.cc main.cc SetParameters.cc trapfpe.c

calc_ccqe_crs_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)

calc_ccqe_crs_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)
