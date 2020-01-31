all: install.include install.library bin/Test 

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

LIBNAME=HT2p2h_$(HT2P2HVER)
$(shell mkdir -p lib)
$(shell mkdir -p bin)

INCS = inc/nieves2p2h.fh

#CFLAGS = -pthread -stdlib=libc++ -std=c++11 -m64 -I/Applications/root-6.04.00/include -Iinc
#LDFLAGS =-L/Applications/root-6.04.00/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,/Applications/root-6.04.00/lib -stdlib=libc++ -lm -ldl

#CFLAGS = -pthread -stdlib=libc++ -std=c++11 -m64 `root-config --cflags` -Iinc
#LDFLAGS = `root-config --glibs` -lminuit -lm -ldl -lgfortran

#CFLAGS=$(COPTFLAGS) $(CXXOPTIONS) $(ROOTINCDIR) -pthread -std=c++11 -m64 -Iinc
CFLAGS=$(COPTFLAGS) $(CXXOPTIONS) $(ROOTINCDIR) -pthread -m64 -Iinc
LDFLAGS=$(ROOTLIBS)

lib/Nucleus2p2h.o: src/Nucleus2p2h.C
	c++ -c $(CFLAGS) -funroll-loops src/Nucleus2p2h.C -o lib/Nucleus2p2h.o 	

lib/HT2p2h.o: src/HT2p2h.C 
	c++ -c $(CFLAGS) -funroll-loops src/HT2p2h.C -o lib/HT2p2h.o 

lib/HadronTensor.o: src/HadronTensor.C 
	c++ -c $(CFLAGS) -funroll-loops src/HadronTensor.C -o lib/HadronTensor.o

lib/HT2p2hFotranInterface.o: src/HT2p2hFotranInterface.C
	c++ -c $(CFLAGS) -funroll-loops src/HT2p2hFotranInterface.C -o lib/HT2p2hFotranInterface.o 

lib/lib${LIBNAME}.a: lib/HT2p2hFotranInterface.o lib/HadronTensor.o lib/HT2p2h.o lib/Nucleus2p2h.o
	ar rcs lib/lib${LIBNAME}.a lib/HT2p2hFotranInterface.o lib/HadronTensor.o lib/HT2p2h.o lib/Nucleus2p2h.o

bin/Test: tst/Test.C lib/libHT2p2h_$(HT2P2HVER).a
	c++ -funroll-loops tst/Test.C  -o bin/Test -Llib/ -l${LIBNAME}  -L/usr/local/lib -lgfortran  $(CFLAGS) $(LDFLAGS) $(CERNLIBS)

lib/Test.o: tst/Test.F
	gfortran -o lib/Test.o -c tst/Test.F -Iinc

bin/TESTFORTRAN: lib/Test.o lib/HT2p2hFotranInterface.o  lib/HT2p2h.o lib/HadronTensor.o 
	c++ -o bin/TESTFORTRAN lib/Test.o -Llib/ -l${LIBNAME} -L/usr/local/lib -lgfortran -stdlib=libc++ -lm -ldl $(CERNLIBS)

install.library: ./lib/lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib

install.include: $(INCS)
	for i in $? ;do ${FH2H} $$i > `echo $$i | sed "s/\.fh$$/C.h/"`; done
	${INSTALL} -m 0644 ${INCS} ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS:.fh=C.h} ${NEUT_ROOT}/include

