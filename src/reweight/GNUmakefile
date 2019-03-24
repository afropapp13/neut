export NEUT_ROOT = ../..

include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

# Check the env. variables for CERNLIB
ifndef CERN
 $(error CERN is not set!)
endif

ifndef CERN_LEVEL
 $(error CERN_LEVEL is not set!)
endif

# Check that we can find cernlib in shell
ifeq "$(shell cernlib)" ""
  $(error Looks like you havent set up CERNLIB, please fix this)
  $(error I rely on being able to do cernlib in shell)
endif

ifndef ROOTSYS
 $(error ROOTSYS is not set!)
endif

NEUT_SRC =  $(NEUT_ROOT)/src

ROOTINC = -I`root-config --incdir`
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff -I../../include -I./ -I../neutclass
INCDIRS = ${ROOTINC} ${NEUTINC}

FCOPTIONS = -g

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive -fPIC


LIBDIRS =	-L../../lib

#MCLIB = ${LIBDIRS}   -lneutcore_${NEUTCOREVER} -lnuceff_${NUCEFFVER} \
#		-lneutcore_${NEUTCOREVER} -lnuccorrspl_${NUCCORVER} \
#		-lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER}

NEUTLIBS := $(wildcard ${LIBDIRS}/neut*.so)
MCLIB = ${LIBDIRS} ${NEUTLIBS}
#MCLIB = ${LIBDIRS} -lnuceff_${NUCEFFVER} -lneutcore_${NEUTCOREVER} \
        -lskmcsvc_${SKMCSVCVER} -lnuccorrspl_${NUCCORVER} \
	-lradcorr_${RADCORRVER}		

# CERNLIB dependencies
CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
	${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
	`cernlib jetset74 photos202 mathlib packlib kernlib`

ROOTLIBS  = `root-config --libs`
LIBS      = $(ROOTLIBS) ${MCLIB} ${CLIBS} -lstdc++

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

#DICTIONARY                    = _ROOT_DICT_ReWeight
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
all: lib

#
#
rootcint: $(DICTIONARY).o

$(DICTIONARY).cc: LinkDef.h $(DICTGEN_OBJECTS)
	$(RM) $(DICTIONARY).*
	@echo "[package: $(PACKAGE)] Generating ROOT dictionary ..."
	rootcint -f $(DICTIONARY).cc -c ${ROOTINC} ${NEUTINC} -I$(NEUT_ROOT)/include -I$(NEUT_SRC)/reweight $(DICTGEN_HEADERS) LinkDef.h

SOCMD     = $(FC)
SOFLAGS   = -fPIC -shared -Wl,-soname,$(LIBNAME).$(DllSuf)
OutPutOpt = -o

#############################


.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
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
	$(SOCMD) $(SOFLAGS) $(SOMINF) $(PACKAGE_ALL_OBJECTS) $(MCLIB) $(EXTRA_EXT_LIBS) $(OutPutOpt) $(PACKAGE_LIB_WITH_PATH)

