CXX = gcc
CXXFLAGS = -O3 -Wall -I./Cuba-1.6
FC  = g77
FCFLAGS = -O3 -Wall -ffortran-bounds-check

%.o : %.c
	$(CXX) $(CXXFLAGS) -c $*.c
%.o : %.f
	$(FC) $(FCFLAGS) -c $*.f

OBJS = CrossSection.o crosin.o crosto.o getr.o dbesi0.o

LIB = -L./Cuba-1.6 -lcuba -lm -lg2c

main: $(OBJS)
	g++ $(CXXFLAGS) -o main $(OBJS) $(LIB)

clean:
	rm -f *.o main
