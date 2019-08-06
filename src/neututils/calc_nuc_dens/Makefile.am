include $(top_srcdir)/build.common

AM_FFLAGS = $(DEFAULT_FFLAGS)
AM_CXXFLAGS = $(DEFAULT_CXXFLAGS)

bin_PROGRAMS = calc_nuc_dens

calc_nuc_dens_SOURCES = integrate.cc integrate_rho.F

calc_nuc_dens_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)

calc_nuc_dens_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)
