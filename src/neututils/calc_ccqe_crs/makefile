CXX = gcc
CXXFLAGS = -O2 -Wall -g

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $*.cc

LIBOBJS = Initialize.o SetParameters.o CrossSection.o bbba07.o clenshaw_curtis.o
LIB = librecipes.a
$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
lib: $(LIB)

main: main.o librecipes.a
	g++ $(CXXFLAGS) -o main main.o -lm $(LIB)

TO_CLEAN_FILES += $(LIB) main
TO_CLEAN_FILES += *.o 
TO_CLEAN_FILES += *~ cross.dat

clean:
	rm -f $(TO_CLEAN_FILES) 
