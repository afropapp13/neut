include $(top_srcdir)/build.common

AM_CXXFLAGS = -DNREWEIGHT_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) $(NEUT_CLASS_INCLUDES)

LIBNAME = libNEUTReWeight.@NEUT_VERSION_NUM@.a

lib_LIBRARIES = $(LIBNAME)

# Don't seem to be able to have canonical names easily build from makefile variables, this is ugly, but works
SRCLIST = NReWeight.cc NSyst.cc NuXSecCCQE.cc NuXSecRES.cc CascadePion.cc

libNEUTReWeight_@NEUT_VERSION_NUM@_a_SOURCES = $(SRCLIST)

include_HEADERS = NReWeight.h NSyst.h NReWeightEngineI.h NuXSecCCQE.h NuXSecRES.h CascadePion.h NReWeightFactory.h

# Unused sources that live in the directory
UNUSED_srcs =
UNUSED_headers =

# bin_PROGRAMS =

# ReValidWeight_SOURCES = ReValidWeight.cc
# ReValidWeight_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
# ReValidWeight_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

# if NEUT_HEPMC

# bin_PROGRAMS += ReValidWeight_HepMC

# ReValidWeight_HepMC_SOURCES = ReValidWeight_HepMC.cc
# ReValidWeight_HepMC_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
# ReValidWeight_HepMC_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

# endif

# PlotReValidWeight_SOURCES = PlotReValidWeight.cc
# PlotReValidWeight_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
# PlotReValidWeight_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

# SummaryPlots_SOURCES = SummaryPlots.cc
# SummaryPlots_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
# SummaryPlots_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)
