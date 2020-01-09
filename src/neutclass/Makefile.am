include $(top_srcdir)/build.common

AM_CXXFLAGS = -DNEUTCLASSVER=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)
AM_LDFLAGS = $(ROOT_LINKFLAGS) $(DEFAULT_LINKFLAGS)

NEUTCLASSUTILS_SRCS = neutrootfortwrapper.cc NeutRootHandlers.cc neutrootTreeSingleton.cc neutfill.cc fill_dsdq2.cc CommonBlockIFace.cc

######### ROOT Standalone I/O library

##### if compat
# if we want each class to be a separate so
if NEUT_LIBCOMPAT

neutctrlLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutCtrl;\n#endif" > $@

neutfsipartLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutFsiPart;\n#endif" > $@

neutfsivertLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutFsiVert;\n#endif" > $@

neutnucfsistepLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutNucFsiStep;\n#endif" > $@

neutnucfsivertLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutNucFsiVert;\n#endif" > $@

neutpartLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutPart;\n#endif" > $@

neutvectLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutVect;\n#endif" > $@

neutvtxLinkDef.h:
	echo -e "#ifdef __CINT__\n#pragma link off all globals;\n#pragma link off all classes;\n#pragma link off all functions;\n#pragma link C++ class NeutVtx;\n#endif" > $@


if ROOT_MAJOR_IS_SIX
%Dict.cc %Dict_rdict.pcm: %.h %LinkDef.h
	rootcint -f $@ -c $^
pcmdir = $(prefix)/lib
pcm_DATA = neutctrlDict_rdict.pcm neutfsipartDict_rdict.pcm neutfsivertDict_rdict.pcm neutnucfsistepDict_rdict.pcm neutnucfsivertDict_rdict.pcm neutpartDict_rdict.pcm neutvectDict_rdict.pcm neutvtxDict_rdict.pcm
else
%Dict.cc: %.h %LinkDef.h
	rootcint -f $@ -c $^
endif

lib_LTLIBRARIES = libneutctrl.la libneutfsipart.la libneutfsivert.la libneutnucfsistep.la libneutnucfsivert.la libneutpart.la libneutvect.la libneutvtx.la libneutrootTreeSingleton.la

libneutctrl_la_SOURCES = neutctrlDict.cc neutctrl.cc
libneutfsipart_la_SOURCES = neutfsipartDict.cc neutfsipart.cc
libneutfsivert_la_SOURCES = neutfsivertDict.cc neutfsivert.cc
libneutnucfsistep_la_SOURCES = neutnucfsistepDict.cc neutnucfsistep.cc
libneutnucfsivert_la_SOURCES = neutnucfsivertDict.cc neutnucfsivert.cc
libneutpart_la_SOURCES = neutpartDict.cc neutpart.cc
libneutvect_la_SOURCES = neutvectDict.cc neutvect.cc
libneutvtx_la_SOURCES = neutvtxDict.cc neutvtx.cc
libneutrootTreeSingleton_la_SOURCES = $(NEUTCLASSUTILS_SRCS)

#you may think this is pointless, but look up autotools renaming objects...
libneutrootTreeSingleton_la_CXXFLAGS = $(AM_CXXFLAGS)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) libneutctrl.so neutctrl.so)
	-(cd ${libdir}; $(LN_S) libneutfsipart.so neutfsipart.so)
	-(cd ${libdir}; $(LN_S) libneutfsivert.so neutfsivert.so)
	-(cd ${libdir}; $(LN_S) libneutnucfsistep.so neutnucfsistep.so)
	-(cd ${libdir}; $(LN_S) libneutnucfsivert.so neutnucfsivert.so)
	-(cd ${libdir}; $(LN_S) libneutpart.so neutpart.so)
	-(cd ${libdir}; $(LN_S) libneutvect.so neutvect.so)
	-(cd ${libdir}; $(LN_S) libneutvtx.so neutvtx.so)
	-(cd ${libdir}; $(LN_S) libneutrootTreeSingleton.so neutrootTreeSingleton.so)

uninstall-hook:
	-(cd ${libdir}; rm -f neutctrl.so)
	-(cd ${libdir}; rm -f neutfsipart.so)
	-(cd ${libdir}; rm -f neutfsivert.so)
	-(cd ${libdir}; rm -f neutnucfsistep.so)
	-(cd ${libdir}; rm -f neutnucfsivert.so)
	-(cd ${libdir}; rm -f neutpart.so)
	-(cd ${libdir}; rm -f neutvect.so)
	-(cd ${libdir}; rm -f neutvtx.so)
	-(cd ${libdir}; rm -f neutrootTreeSingleton.so)


##### if new style where you have one I/O library and one static utils library.
else

NEUTCLASSHEADERS = neutctrl.h neutfsipart.h neutfsivert.h neutnucfsistep.h neutnucfsivert.h neutpart.h neutvect.h neutvtx.h

if ROOT_MAJOR_IS_SIX
neutclassDict.cc neutclassDict_rdict.pcm: $(NEUTCLASSHEADERS) neutclassLinkDef.h
	rootcint -f neutclassDict.cc -c $^
pcmdir = $(prefix)/lib
pcm_DATA = neutclassDict_rdict.pcm
else
neutclassDict.cc: $(NEUTCLASSHEADERS) neutclassLinkDef.h
	rootcint -f neutclassDict.cc -c $^
endif

BUILT_SOURCES = neutclassDict.cc
CLEANFILES = $(BUILT_SOURCES)

include_HEADERS = neutctrl.h neutfsipart.h neutfsivert.h neutnucfsistep.h neutnucfsivert.h neutpart.h NeutRootHandlers.h neutrootTreeSingleton.h neutvect.h neutvtx.h CommonBlockIFace.h

lib_LTLIBRARIES = libNEUTClass.la

libNEUTClass_la_SOURCES = neutclassDict.cc neutctrl.cc neutfsipart.cc neutfsivert.cc neutnucfsistep.cc neutnucfsivert.cc neutpart.cc neutvect.cc neutvtx.cc

if NEUT_HEPMC
libNEUTClass_la_SOURCES += HepMCWriter.cc
libNEUTClass_la_LDFLAGS = $(HEPMC_LINKFLAGS)
libNEUTClass_la_LIBADD = $(HEPMC_LIBS)
endif

#
LIBNAME = libNEUTClassUtils.@NEUT_VERSION_NUM@.a
#
lib_LIBRARIES = $(LIBNAME)
#
libNEUTClassUtils_@NEUT_VERSION_NUM@_a_SOURCES = $(NEUTCLASSUTILS_SRCS)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

# Unused sources that live in the directory
UNUSED_srcs = chkfillneutvect.cc chkfillneutvtx.cc chkreadneutvect.cc chkreadneutvtx.cc

endif
