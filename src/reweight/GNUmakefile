ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

NEUT_SRC =  $(NEUT_ROOT)/src

ROOTINC = -I`root-config --incdir`
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff -I../../include -I./
NEUTINC += -I../neutclass -I../specfunc
INCDIRS = ${ROOTINC} ${NEUTINC}

#FCOPTIONS = -g
FCOPTIONS = -O2 -march=native -Wall

CXXFLAGS  = `root-config --cflags`
#COPTFLAGS = -g -Wno-write-strings -fpermissive -fPIC
COPTFLAGS = -O2 -march=native -Wall -Wno-write-strings -fpermissive -fPIC


COREDIR =	../neutcore/Linux_pc
NUCEFFDIR =     ../nuceff/Linux_pc

LIBDIRS =	-L${COREDIR}         -L../nuccorspl/Linux_pc \
			-L../nuceff/Linux_pc -L../partnuck/Linux_pc \
			-L../skmcsvc/Linux_pc -L../specfunc

NEUTCOREVER = 5.3.3
INEUTCOREVER= 533

NUCEFFVER   = 5.3.3
INUCEFFVER  = 533

NUCCORVER   = 5.3.3
INUCCORVER  = 533

PARTNUCKVER = 5.3.3
IPARTNUCKVER= 533

SKMCSVCVER  = 5.3.3
ISKMCSVCVER = 533

SPECFUNCVER  = 5.3.3
ISPECFUNCVER = 533

FORTRANDEFINES += -DNECORE_VERSION=$(INEUTCOREVER)
FORTRANDEFINES += -DNENUCE_VERSION=$(INUCEFFVER)
FORTRANDEFINES += -DNENUCC_VERSION=$(INUCCORVER)

#MCLIB = ${LIBDIRS}   -lneutcore_${NEUTCOREVER} -lnuceff_${NUCEFFVER} \
#		-lneutcore_${NEUTCOREVER} -lnuccorrspl_${NUCCORVER} \
#		-lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER}

MCLIB = ${LIBDIRS} -lnuceff_${NUCEFFVER} -lspecfunc_${SPECFUNCVER} \
        -lneutcore_${NEUTCOREVER} -lnuceff_${NUCEFFVER} \
		-lnuccorrspl_${NUCCORVER} -lskmcsvc_${SKMCSVCVER} 
		

CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
		${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
		`cernlib jetset74 photos202 mathlib packlib kernlib`

ROOTLIBS  = `root-config --libs`
LIBS      = $(ROOTLIBS) ${MCLIB} ${CLIBS} -lstdc++ -lgfortran

LDOPTFLAGS= -g

FINCDIRS  = ${NEUTINC}

POBJS = ${COREDIR}/structm.o ${COREDIR}/pdfset.o ${COREDIR}/grv94di.o \
		${COREDIR}/grv98_lo.o
MOBJS =

OBJS  = ${MOBJS} ${ROBJS} ${SOBJS} ${POBJS}

# For ROOT shared library
PACKAGE                       = ReWeight
LIBNAME                       = libNReWeight
DllSuf                        = so
PACKAGE_LIB                   = $(LIBNAME).$(DllSuf)
PACKAGE_LIB_WITH_PATH         = $(NEUT_SRC)/reweight/$(LIBNAME).$(DllSuf)

DICTIONARY                    = _ROOT_DICT_ReWeight
ObjSuf                       := o
PACKAGE_ALL_SOURCES          := $(wildcard *.cc)
DICTGEN_HEADERS              := $(addsuffix .h, $(basename $(PACKAGE_ALL_SOURCES)))
DICTGEN_OBJECTS              := $(addsuffix .$(ObjSuf), $(basename $(DICTGEN_HEADERS)))
PACKAGE_ALL_SOURCES          := $(PACKAGE_ALL_SOURCES) $(DICTIONARY).cc $(wildcard *.cpp)
PACKAGE_ALL_OBJECTS           = $(addsuffix .$(ObjSuf), $(basename $(PACKAGE_ALL_SOURCES)))
PACKAGE_SOURCES              := $(wildcard *.cc)
PACKAGE_OBJECTS               = $(addsuffix .$(ObjSuf), $(basename $(PACKAGE_SOURCES)))

FORT_OBJECTS                  = $(COREDIR)/evdifcrs.o $(COREDIR)/nesetfgparams.o $(COREDIR)/nefillmodel.o $(NUCEFFDIR)/evpiprob.o

LIB_DEPEND=
ifeq ($(strip $(DICTIONARY)),undefined)
  LIB_DEPEND = $(PACKAGE_OBJECTS) 
else
  LIB_DEPEND = $(PACKAGE_ALL_OBJECTS)
endif

# default target to build the library
all: clean lib

#
#
rootcint: $(DICTIONARY).o

$(DICTIONARY).cc: LinkDef.h $(DICTGEN_OBJECTS)
	$(RM) $(DICTIONARY).*
	@echo "[package: $(PACKAGE)] Generating ROOT dictionary ..."
	rootcint -f $(DICTIONARY).cc -c -I$(ROOTSYS)/include -I$(NEUT_ROOT)/include -I$(NEUT_SRC)/reweight $(DICTGEN_HEADERS) LinkDef.h

SOCMD     = $(FC)
SOFLAGS   = -fPIC -shared -Wl,-soname,$(LIBNAME).$(DllSuf)
OutPutOpt = -o

#############################


.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS} ${LIBS}" INCDIRS="${INCDIRS}" \
			./bin/gen_root_so.sh

.cc.o:
	$(CXX) -c -w $(COPTFLAGS) $(ALLDEFINES) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(FINCDIRS) -o $@ $<

clean:
	$(RM) *.o *~ ${MOBJS} *Dict.* *.so _ROOT_DICT_ReWeight.*

lib: $(LIB_DEPEND)
	@echo "[package: $(PACKAGE)] Generating shared library ..."
	$(SOCMD) $(SOFLAGS) $(SOMINF) $(PACKAGE_ALL_OBJECTS) $(LIBS) $(EXTRA_EXT_LIBS) $(OutPutOpt) $(PACKAGE_LIB_WITH_PATH)

