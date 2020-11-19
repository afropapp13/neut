# Just check that NEUT core is built
ifeq "$(wildcard ../../../lib/libneutcore_*)" ""
  $(error Are you sure you have built NEUT? The neutcore library cant be found!))
endif

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

# DONE WITH CHECKS ---------------------------

SRC=$(wildcard *.cpp)
TRGS=$(patsubst %.cpp,%,${SRC})

${TRGS} : ${SRC}

# CERNLIBs
CERNLIB=$(shell cernlib jetset74 photos202 mathlib packlib kernlib pdflib)

# NEUTLIBS
NEUTLIBS=-L../../../lib -lNReWeight -lgfortran $(patsubst ../../../lib/lib%.a,-l%,$(wildcard ../../../lib/*.a)) #../../neutsmpl/rmarin_dum.o

# Includes
INC = -I../ -I../../../include -I${CERN}/${CERN_LEVEL}/include

COREDIR = ../../neutcore

# Also need to include some objects which for some reason aren't in neutcore
POBJS = ${COREDIR}/structm.o ${COREDIR}/pdfset.o #${COREDIR}/grv94di.o

ROOTLIB = $(shell root-config --libs)
ROOTINC = $(shell root-config --cflags)

# Libraries for NEUT class
NEUTCLASS=$(wildcard ../../neutclass/neut*.so)

INC += ${ROOTINC}
LIB =  ${NEUTLIBS}  ${NEUTCLASS} ${POBJS} ${CERNLIB} ${ROOTLIB} 

%: %.cpp
	g++ -g -Wl,-rpath=../../../lib,--start-group -o $@ $< ${INC}  ${LIB} -Wl,--end-group

clean:
	rm -f ${TRGS}
