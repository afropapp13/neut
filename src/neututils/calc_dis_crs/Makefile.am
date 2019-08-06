include $(top_srcdir)/build.common

AM_FFLAGS = $(DEFAULT_FFLAGS)

bin_PROGRAMS = calc_dis_crs

calc_dis_crs_SOURCES = calcapicrs.F tpispir.F

calc_dis_crs_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)

calc_dis_crs_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)
