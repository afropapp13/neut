include $(top_srcdir)/build.common

AM_FFLAGS = $(DEFAULT_FFLAGS)
AM_CXXFLAGS = $(DEFAULT_CXXFLAGS) $(NEUT_CLASS_INCLUDES)

bin_PROGRAMS = cFORTRAN_boundary neut_test_card_reader

cFORTRAN_boundary_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS) -lgfortran
cFORTRAN_boundary_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)

cFORTRAN_boundary_SOURCES = cFORTRAN_boundary.cc double_me.F

neut_test_card_reader_SOURCES = neut_test_card_reader.cc
neut_test_card_reader_LDADD = $(NEUT_LIBS) $(CERN_LIBS) $(ROOT_LIBS) $(EXE_RTLIBS)
neut_test_card_reader_LDFLAGS = $(EXE_LINKFLAGS) $(CERN_LINKFLAGS) $(ROOT_LINKFLAGS)