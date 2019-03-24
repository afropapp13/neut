ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

COREDIR = ${NEUT_ROOT}/src/neutcore
CLASSDIR= ${NEUT_ROOT}/src/neutclass

OBJS  = 

NEUTLIBS += -lskmcsvc_${SKMCSVCVER} -lzbsfns_${ZBSFNSVER} \
			-L${SKOFL_ROOT}/lib -lmon -liolib -lzbs -llibrary \
			$(CERNLIBS) -lUnuran -lrfa

T2KOBJS = necardbm.o  \
		  necpjndfxcm.o nemkjndfx.o nerdnufxv.o rnbmom2802.o \
		  tauvect.o neopnufxv.o nerd1nufxv.o \
          neopskfxv.o rnpos.o reset_signal_.o


NTPLTQOBJS  = neutntpl_tq.o
NTOBJS= nerdjndfx.o  
T2KSKOBJS   = t2kflux_sk.o ranlux_trandom3.o ${T2KOBJS} ${OBJS}

FCOPTIONS  += -DFLUX_13
FCOPTFLAGS = -g -O0 

CXXOPTIONS += -I./ -I${SKOFL_ROOT}/include

FLAGS      = ${FOPTFLAGS} ${FCOPTIONS} $(FCFLAGS) 

TARGETS = t2kneut t2kneut_ndall t2kneut_sk skflux_dump neutdump

all: $(TARGETS)

necardbmC.h: necardbm.h	
	rm -f necardbmC.h
	fh2h.pl necardbm.h > necardbmC.h

beamntplC.h: beamntpl.h	
	rm -f beamntplC.h
	fh2h.pl beamntpl.h > beamntplC.h

skheadC.h: ${SKOFL_ROOT}/inc/skhead.h
	rm -f skheadC.h
	fh2h.pl ${SKOFL_ROOT}/inc/skhead.h > skheadC.h

t2kneut_sk.o: necardbmC.h beamntplC.h skheadC.h

t2kneut: t2kneut.o necardbm.h $(T2KOBJS) 
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` t2kneut.o ${T2KOBJS} $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

t2kneut_ndall: t2kneut_ndall.o necardbm.h $(T2KOBJS) 
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` t2kneut_ndall.o ${T2KOBJS} $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

t2kneut_sk: necardbm.h t2kneut_sk.o $(T2KSKOBJS) 
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` t2kneut_sk.o ${T2KSKOBJS} $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

skflux_dump: skflux_dump.o necardbm.h $(T2KOBJS) 
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` skflux_dump.o ${T2KSKOBJS} $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

split_data: split_data.o fort_fopen.o atoi_.o concat_.o
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` \
		  split_data.o fort_fopen.o atoi_.o concat_.o ${LIBS} -liolib

empdist: empdist.o
	$(CXX) ${LDFLAGS} -o $@ empdist.o ${ROOTLIBS}

neutntpl: neutntpl.o $(NTOBJS)
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` \
		  neutntpl.o ${NTOBJS}  $(NEUTLIBS) $(ROOTLIBS) 

neutntpl_tq: $(NTPLTQOBJS) $(NTOBJS)
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` \
		  ${NTPLTQOBJS} ${NTOBJS}  $(NEUTLIBS) $(ROOTLIBS) 

neutdump: neutdump.o
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` \
		  neutdump.o $(NEUTLIBS) $(ROOTLIBS) 

neut_test: neut_test.o $(T2KOBJS)  ${OBJS}
	$(FC) -o $@ $(FLAGS) -Xlinker -R`pwd` \
		neut_test.o ${T2KOBJS} ${OBJS} ${LIBS}

t2kflux_sktest.o : t2kflux_sk.h
t2kflux_sk.o     : t2kflux_sk.h
