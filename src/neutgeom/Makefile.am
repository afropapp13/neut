include $(top_srcdir)/build.common

AM_FFLAGS = -DNEUTGEOM_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DNEUTGEOM_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) -I$(top_srcdir)/neutclass

LIBNAME = lib@LIBPREFIX@neutgeom.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

lib@LIBPREFIX@neutgeom_@NEUT_VERSION_NUM@_a_SOURCES = neutev.F neutxs.F TNuTrajectory.cc TNeutOutput.cc TNuFlux.cc

bin_PROGRAMS = event_rate genev dumptotpau_nd280

event_rate_SOURCES = event_rate.cc
event_rate_LDADD = lib@LIBPREFIX@neutgeom.@NEUT_VERSION_NUM@.a $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) -lEve -lEG -lGeom $(EXE_RTLIBS)
event_rate_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

genev_SOURCES = genev.cc
genev_LDADD = lib@LIBPREFIX@neutgeom.@NEUT_VERSION_NUM@.a $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) -lEve -lEG -lGeom $(EXE_RTLIBS)
genev_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

dumptotpau_nd280_SOURCES = dumptotpau_nd280.cc
dumptotpau_nd280_LDADD = lib@LIBPREFIX@neutgeom.@NEUT_VERSION_NUM@.a $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) -lEve -lEG -lGeom $(EXE_RTLIBS)
dumptotpau_nd280_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)
