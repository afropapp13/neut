include $(top_srcdir)/build.common

AM_FFLAGS = -DNEUTSMPL_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DNEUTSMPL_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) -I$(top_srcdir)/neutclass

bin_PROGRAMS = neutroot2

NEUTSMPL_PROGLDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
NEUTSMPL_PROGLDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

neutroot2_SOURCES = neutroot.F flx_intrfc.cc ufm2th1dist.cc trapfpe.c
neutroot2_LDADD = $(NEUTSMPL_PROGLDADD)
neutroot2_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

if NEUT_USE_GPERF
neutroot2_SOURCES += gperfstart.cc
endif

if !NEUT_LIBCOMPAT

bin_PROGRAMS += neut neut_ntpl dumptotpau dumpelspau dumpcohcrs dumpcrs dumptotpauC basic_histos neutclass_to_tree cross_section_histos

neut_SOURCES = neut.F
neut_LDADD = $(NEUTSMPL_PROGLDADD)
neut_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

neut_ntpl_SOURCES = neut_ntpl.F
neut_ntpl_LDADD = $(NEUTSMPL_PROGLDADD)
neut_ntpl_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

dumptotpau_SOURCES = dumptotpau.F
dumptotpau_LDADD = $(NEUTSMPL_PROGLDADD)
dumptotpau_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

dumpelspau_SOURCES = dumpelspau.F
dumpelspau_LDADD = $(NEUTSMPL_PROGLDADD)
dumpelspau_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

dumpcohcrs_SOURCES = dumpcohcrs.F
dumpcohcrs_LDADD = $(NEUTSMPL_PROGLDADD)
dumpcohcrs_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

dumpcrs_SOURCES = dumpcrs.F fnrawcrs.F
dumpcrs_LDADD = $(NEUTSMPL_PROGLDADD)
dumpcrs_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

dumptotpauC_SOURCES = dumptotpauC.cc
dumptotpauC_LDADD = $(NEUTSMPL_PROGLDADD)
dumptotpauC_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

basic_histos_SOURCES = basic_histos.cc
basic_histos_LDADD = $(NEUTSMPL_PROGLDADD)
basic_histos_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

neutclass_to_tree_SOURCES = neutclass_to_tree.cc
neutclass_to_tree_LDADD = $(NEUTSMPL_PROGLDADD)
neutclass_to_tree_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

cross_section_histos_SOURCES = cross_section_histos.cc fnrawcrs.F
cross_section_histos_LDADD = $(NEUTSMPL_PROGLDADD)
cross_section_histos_LDFLAGS = $(NEUTSMPL_PROGLDFLAGS)

endif
