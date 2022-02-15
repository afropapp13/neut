include $(top_srcdir)/build.common

AM_FFLAGS = $(DEFAULT_FFLAGS)
AM_CXXFLAGS = $(DEFAULT_CXXFLAGS)
AM_CFLAGS = $(DEFAULT_CFLAGS)

bin_PROGRAMS = calc_res_crs calc_res_crs_newff

calc_res_crs_SOURCES = crs.F res_lambda.F res_q2_limits.F trapfpe.c
calc_res_crs_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
calc_res_crs_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

calc_res_crs_newff_SOURCES = crs_newff.F res_lambda.F res_q2_limits.F trapfpe.c
calc_res_crs_newff_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
calc_res_crs_newff_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)
