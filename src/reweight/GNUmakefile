ifndef NEUT_ROOT
  NEUT_ROOT = ../..
endif

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

ROOTINC = -I$(shell root-config --incdir)
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff -I../../include -I./ -I../neutclass -I../specfunc
INCDIRS = ${ROOTINC} ${NEUTINC}

FCOPTIONS = -g -O2 -Wall

CXXFLAGS  = $(shell root-config --cflags)
COPTFLAGS = -g -O2 -Wall -Wno-write-strings -fpermissive -fPIC


NEUT_LIB_DIR = $(NEUT_ROOT)/lib
STRING := $(filter-out $(wildcard ${NEUT_LIB_DIR}/*tau*),$(wildcard ${NEUT_LIB_DIR}/lib*_*.a))
NEUTLIBS = -L$(NEUT_LIB_DIR) -Wl,--as-needed,--start-group $(patsubst ${NEUT_LIB_DIR}/lib%.a,-l%,$(STRING)) -Wl,--end-group
MCLIB=${NEUTLIBS}
LIBDIRS=${NEUT_LIB_DIR}

#MCLIB = ${LIBDIRS}   -lneutcore_${NEUTCOREVER} -lnuceff_${NUCEFFVER} \
#		-lneutcore_${NEUTCOREVER} -lnuccorrspl_${NUCCORVER} \
#		-lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER}

#$(info ${LIBDIRS})
#MCLIB = -L$(NEUT_ROOT)/lib -Wl,--start-group $(patsubst ${LIBDIRS}/lib%.a,-l%,$(wildcard ${LIBDIRS}/lib*.a)) -Wl,--end-group
#NEUTLIBS = $(patsubst ${LIBDIRS}/%,%,$(wildcard ${LIBDIRS}/lib*.a))
#NEUTLIBS = $(wildcard ${LIBDIRS}/lib*.a)
#$(info ${NEUTLIBS})
#MCLIB = -L$(LIBDIRS) ${NEUTLIBS}
MCLIB = -L${LIBDIRS} -lnuceff_${NUCEFFVER} -lneutcore_${NEUTCOREVER} \
        -lskmcsvc_${SKMCSVCVER} -lnuccorspl_${NUCCORVER} \
	-lradcorr_${RADCORRVER}		

# CERNLIB dependencies
CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
	${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
	$(shell cernlib jetset74 photos202 mathlib packlib kernlib)

ROOTLIBS  = $(shell root-config --libs)
LIBS      = ${ROOTLIBS} ${MCLIB} ${CLIBS} -lstdc++ -lgfortran
#LIBS      = $(ROOTLIBS) -lstdc++

LDOPTFLAGS= -g

#FINCDIRS  = ${NEUTINC}

#POBJS = ${COREDIR}/structm.o ${COREDIR}/pdfset.o ${COREDIR}/grv94di.o \
		#${COREDIR}/grv98_lo.o
#MOBJS =

#OBJS  = ${MOBJS} ${ROBJS} ${SOBJS} ${POBJS}

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
all: lib

#
#
rootcint: $(DICTIONARY).o

$(DICTIONARY).cc: LinkDef.h $(DICTGEN_OBJECTS)
	-$(RM) $(DICTIONARY).*
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
	-$(RM) *.o *~ ${MOBJS} *Dict.* *.so _ROOT_DICT_ReWeight.*

lib: $(LIB_DEPEND)
	@echo "[package: $(PACKAGE)] Generating shared library ..."
	$(SOCMD) $(SOFLAGS) $(SOMINF) $(PACKAGE_ALL_OBJECTS) $(LIBS) $(EXTRA_EXT_LIBS) $(OutPutOpt) $(PACKAGE_LIB_WITH_PATH)

