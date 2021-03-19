ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
all: nucleontester.exe piontester.exe nuclscat2.exe nuclscat.exe nuclscat_anal.exe
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

#NEUTVER6 = 5.3.6
#NEUTVER5 = 5.3.5
NEUTVER6 = 5.4.0
NEUTVER5 = 5.4.0
TAUVER = 2.0.2

NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff -I../../include

COMMON_BLOCK_STATIC_DEFINES = -DSTATIC_COMMON_POINTERS -DNO_EXTERN_COMMON_POINTERS

CXXFLAGS  = $(shell root-config --cflags) $(NEUTINC) -g -Wno-write-strings -fpermissive $(COMMON_BLOCK_STATIC_DEFINES) -O0

CERNLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
	${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
	$(shell cernlib jetset74 pdflib804 photos202 mathlib packlib kernlib)

ROOTLIBS = $(shell root-config --libs)

RMRARINDUM = $(NEUT_ROOT)/src/neutsmpl/rmarin_dum.o

NEUTLIBDIRS = -L../../lib
NEUTLIBS = -lnuccorspl_${NEUTVER5} -lneutcore_${NEUTVER6} -lnuceff_${NEUTVER5} -lneutcore_${NEUTVER6} -lnuccorspl_${NEUTVER5} -lpartnuck_${NEUTVER5} -lskmcsvc_${NEUTVER5} -lradcorr_${NEUTVER5} -lHT2p2h_${NEUTVER5} -lN1p1h_${NEUTVER5} -lspecfunc_${NEUTVER5}

LDFLAGS = $(NEUTLIBDIRS) $(NEUTLIBS) $(NEUTLIBS) tauvect.o $(ROOTLIBS) rmarin.o $(CERNLIBS) -lgfortran


rmarin.o : rmarin_dum.F
	gfortran $< -c -o $@
tauvect.o : tauvect.F
	gfortran $< -c -o $@

nuclscat2.exe: nuclscat2.cxx NEUTCIncludes.hxx nuclscat_output_blocks.hxx rmarin.o tauvect.o
	g++ $(CXXFLAGS) -g -o $@ $<  $(LDFLAGS)

nuclscat.exe: nuclscat.cxx NEUTCIncludes.hxx nuclscat_output_blocks.hxx rmarin.o tauvect.o
	g++ $(CXXFLAGS) -g -o $@ $<  $(LDFLAGS)


piontester.exe: piontester.cxx NEUTCIncludes.hxx nuclscat_output_blocks.hxx rmarin.o tauvect.o
	g++ $(CXXFLAGS) -g -o $@ $<  $(LDFLAGS)


nucleontester.exe: nucleontester.cxx NEUTCIncludes.hxx nuclscat_output_blocks.hxx rmarin.o tauvect.o
	g++ $(CXXFLAGS) -g -o $@ $<  $(LDFLAGS)



nuclscat_anal.exe: nuclscat_anal.cxx nuclscat_output_blocks.hxx rmarin.o tauvect.o
	g++ $(CXXFLAGS) -g -o $@ $< $(LDFLAGS)

clean:
	rm -f nuclscat.exe nuclscat_anal.exe nuclscat2.exe
