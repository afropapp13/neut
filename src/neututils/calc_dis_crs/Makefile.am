include $(top_srcdir)/build.common

AM_FFLAGS = $(DEFAULT_FFLAGS)

bin_PROGRAMS = calc_dis_crs CalcQ2dXsec2 CalcQ2Xsec2

calc_dis_crs_SOURCES = calcapicrs.F tpispir.F

calc_dis_crs_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)

calc_dis_crs_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

CalcQ2dXsec2_SOURCES = CalcQ2dXsec2.F tpispir.F

CalcQ2dXsec2_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)

CalcQ2dXsec2_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

CalcQ2Xsec2_SOURCES = CalcQ2Xsec2.F tpispir.F

CalcQ2Xsec2_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)

CalcQ2Xsec2_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)
