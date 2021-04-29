
SOURCES=$(wildcard *.cpp)
TARGETS=$(patsubst %.cpp,%.exe,${SOURCES})

ROOT_INC=$(shell root-config --cflags)
ROOT_LIB=$(shell root-config --libs)

NEUT_INC=-I$(shell neut-config --incdir)
NEUT_LIB=-L $(shell neut-config --libdir) $(shell neut-config --libflags)
NEUT_IO_LIB=-L $(shell neut-config --libdir) $(shell neut-config --iolibflags)
NEUT_RW_LIB=-L $(shell neut-config --libdir) $(shell neut-config --rwlibflags)
CERN_LIB=$(shell neut-config --cernflags)

all: ${TARGETS}

${TARGETS}: ${SOURCES}

%.exe: %.cpp
	g++ -g $(NEUT_INC) $(ROOT_INC) -o $@ $< $(NEUT_LIB) $(NEUT_IO_LIB) $(NEUT_RW_LIB) $(ROOT_LIB) $(CERN_LIB) -lgfortran

clean: 
	rm -rfv ${TARGETS}
