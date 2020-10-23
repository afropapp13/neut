include $(top_srcdir)/build.common

AM_CXXFLAGS = -DNREWEIGHT_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) $(NEUT_CLASS_INCLUDES)

if NEUT_LIBCOMPAT
LIBNAME = libNReWeight.@NEUT_VERSION_NUM@.a
else
LIBNAME = libNEUTReWeight.@NEUT_VERSION_NUM@.a
endif

lib_LIBRARIES = $(LIBNAME)

# Don't seem to be able to have canonical names easily build from makefile variables, this is ugly, but works

SRCLIST = NReWeightCasc.cc NReWeightNuXSecCCQE.cc NReWeightNuXSecRES.cc NReWeightFactory.cc NReWeight.cc NSystSet.cc NSystUncertainty.cc

if NEUT_LIBCOMPAT
libNReWeight_@NEUT_VERSION_NUM@_a_SOURCES = $(SRCLIST)
else
libNEUTReWeight_@NEUT_VERSION_NUM@_a_SOURCES = $(SRCLIST)
endif

include_HEADERS = NModeDefn.h NReWeightCasc.h  NReWeightNuXSecCCQE.h    NReWeightNuXSecRES.h NReWeightFactory.h NReWeight.h NReWeightI.h NSyst.h NSystSet.h NSystUncertainty.h PDGCodes.h EventSummaryHelper.h

# Unused sources that live in the directory
UNUSED_srcs = NReWeightNuclPiless.cc NReWeightNuXSecCOH.cc NReWeightNuXSecDIS.cc NReWeightNuXSecNCEL.cc
UNUSED_headers = NReWeightNuclPiless.h NReWeightNuXSecCOH.h NReWeightNuXSecDIS.h NReWeightNuXSecNCEL.h

bin_PROGRAMS = ReValidWeight ReValidWeight_HepMC SummaryPlots PlotReValidWeight
ReValidWeight_SOURCES = ReValidWeight.cc
ReValidWeight_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
ReValidWeight_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

ReValidWeight_HepMC_SOURCES = ReValidWeight_HepMC.cc
ReValidWeight_HepMC_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
ReValidWeight_HepMC_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)


PlotReValidWeight_SOURCES = PlotReValidWeight.cc
PlotReValidWeight_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
PlotReValidWeight_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

SummaryPlots_SOURCES = SummaryPlots.cc
SummaryPlots_LDADD = $(LIBNAME) $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
SummaryPlots_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)
