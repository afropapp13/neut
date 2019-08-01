include $(top_srcdir)/build.common

NEUTCLASSHEADERS = neutctrl.h neutfsipart.h neutfsivert.h neutnucfsistep.h neutnucfsivert.h neutpart.h neutvect.h neutvtx.h

neutclassDict.cc neutclassDict_rdict.pcm: $(NEUTCLASSHEADERS) neutclassLinkDef.h
	rootcint -f neutclassDict.cc -c $^

BUILT_SOURCES = neutclassDict.cc
CLEANFILES = $(BUILT_SOURCES)

AM_CXXFLAGS = -DNEUTCLASSVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)
AM_LDFLAGS = $(ROOT_LINKFLAGS)

lib_LTLIBRARIES = libNEUTneutclass.la

if ROOT_MAJOR_IS_SIX
pcmdir = $(prefix)/lib
pcm_DATA = neutclassDict_rdict.pcm
endif

libNEUTneutclass_la_SOURCES = neutclassDict.cc neutctrl.cc neutfsipart.cc neutfsivert.cc neutnucfsistep.cc neutnucfsivert.cc neutpart.cc neutvect.cc neutvtx.cc

lib_LIBRARIES = libNEUTneutclassUtils.a
libNEUTneutclassUtils_a_SOURCES = neutrootfortwrapper.cc NeutRootHandlers.cc neutrootTreeSingleton.cc neutfill.cc

# Unused sources that live in the directory
UNUSED_srcs = chkfillneutvect.cc chkfillneutvtx.cc chkreadneutvect.cc chkreadneutvtx.cc
