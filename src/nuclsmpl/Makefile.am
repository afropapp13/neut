include $(top_srcdir)/build.common

AM_FFLAGS = -DNUCLSMPL_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DNUCLSMPL_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) -I$(top_srcdir)/neutclass

bin_PROGRAMS = nuclscat run_nuclscatana run_nucltransana

nuclscat_SOURCES = nuclscat.F tauvect.F
nuclscat_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) -lstdc++
nuclscat_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

run_nuclscatana_SOURCES = run_nuclscatana.cc tauvect.F
run_nuclscatana_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) -lgfortran
run_nuclscatana_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

run_nucltransana_SOURCES = run_nucltransana.cc tauvect.F
run_nucltransana_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) -lgfortran
run_nucltransana_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)
