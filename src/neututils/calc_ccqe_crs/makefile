CXX = gcc
CXXFLAGS = -O0 -Wall -g
FCFLAGS = -O0 -Wall -g

ifndef NEUT_ROOT
	NEUT_ROOT = ../../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

LIBOBJS = SetParameters.o bbba07.o CrossSection.o clenshaw_curtis.o bbba05.o Initialize.o # bbba07_vec.o
LIB = librecipes.a
$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
lib: $(LIB)


main: main.o trapfpe.o librecipes.a 
	$(FC) $(CXXFLAGS) -o main main.o trapfpe.o -lm $(LIB) -lstdc++
#	g++ $(CXXFLAGS) -o main main.o -lm $(LIB) -lg2c

TO_CLEAN_FILES += $(LIB) main
TO_CLEAN_FILES += *.o 
TO_CLEAN_FILES += *~ cross.dat

clean:
	rm -f $(TO_CLEAN_FILES) 
