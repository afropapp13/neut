include $(top_srcdir)/build.common

AM_CXXFLAGS = -DNEUTCLASSVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)
AM_LDFLAGS = $(ROOT_LINKFLAGS) $(DEFAULT_LINKFLAGS)

NEUTCLASSUTILS_SRCS = neutrootfortwrapper.cc NeutRootHandlers.cc neutrootTreeSingleton.cc neutfill.cc fill_dsdq2.cc CommonBlockIFace.cc silence.F

######### ROOT Standalone I/O library
NEUTCLASSHEADERS = neutctrl.h neutfsipart.h neutfsivert.h neutnucfsistep.h neutnucfsivert.h neutpart.h neutvect.h neutvtx.h

if ROOT_MAJOR_IS_SIX
neutclassDict.cc neutclassDict_rdict.pcm: $(NEUTCLASSHEADERS) neutclassLinkDef.h
	rootcint -f neutclassDict.cc -c $(NEUT_DEFINES) $^
pcmdir = $(prefix)/lib
pcm_DATA = neutclassDict_rdict.pcm
else
neutclassDict.cc: $(NEUTCLASSHEADERS) neutclassLinkDef.h
	rootcint -f neutclassDict.cc -c $(NEUT_DEFINES) $^
endif

BUILT_SOURCES = neutclassDict.cc
CLEANFILES = $(BUILT_SOURCES)

include_HEADERS = neutctrl.h neutfsipart.h neutfsivert.h neutnucfsistep.h neutnucfsivert.h neutpart.h NeutRootHandlers.h neutrootTreeSingleton.h neutvect.h neutvtx.h CommonBlockIFace.h

lib_LTLIBRARIES = libNEUTClass.la

libNEUTClass_la_SOURCES = neutclassDict.cc neutctrl.cc neutfsipart.cc neutfsivert.cc neutnucfsistep.cc neutnucfsivert.cc neutpart.cc neutvect.cc neutvtx.cc

if NEUT_HEPMC
libNEUTClass_la_LDFLAGS = $(HEPMC_LINKFLAGS)
libNEUTClass_la_LIBADD = $(HEPMC_LIBS)
endif

#
LIBNAME = libNEUTClassUtils.@NEUT_VERSION_NUM@.a
#
lib_LIBRARIES = $(LIBNAME)
#
libNEUTClassUtils_@NEUT_VERSION_NUM@_a_SOURCES = $(NEUTCLASSUTILS_SRCS)

if NEUT_HEPMC
libNEUTClassUtils_@NEUT_VERSION_NUM@_a_SOURCES += HepMCWriter.cc HepMCReader.cc
endif

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

# Unused sources that live in the directory
UNUSED_srcs = chkfillneutvect.cc chkfillneutvtx.cc chkreadneutvect.cc chkreadneutvtx.cc