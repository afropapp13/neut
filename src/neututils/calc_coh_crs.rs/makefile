CXX = gcc4
CXXFLAGS = -O3
FC  = gfortran
FCFLAGS = -O3

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $*.cc
%.o : %.f
	$(FC) $(FCFLAGS) -c $*.f

LIBOBJS = Initialize.o SetParameters.o CrossSection.o clenshaw_curtis.o \
	crosin.o crosto.o getr.o dbesi0.o
LIB = librecipes.a
$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
lib: $(LIB)

main: main.o librecipes.a
	g++4 $(CXXFLAGS) -o main main.o -lm $(LIB)

TO_CLEAN_FILES += $(LIB) main
TO_CLEAN_FILES += *.o 
TO_CLEAN_FILES += *~ cross.dat

clean:
	rm -f $(TO_CLEAN_FILES) 
