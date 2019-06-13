## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

COREDIR = ${NEUT_ROOT}/src/neutcore
CLASSDIR= ${NEUT_ROOT}/src/neutclass


SOBJS = ${CLASSDIR}/neutctrl.so ${CLASSDIR}/neutpart.so \
		${CLASSDIR}/neutvect.so ${CLASSDIR}/neutvtx.so \
		${CLASSDIR}/neutfsipart.so ${CLASSDIR}/neutfsivert.so \
		${CLASSDIR}/neutnucfsistep.so ${CLASSDIR}/neutnucfsivert.so

OBJS  = nevecgen.o grndmq.o rmarin_dum.o flx_intrfc.o \
		ufm2th1dist.o trapfpe.o

TARGETS = neut neut_ntpl dumptotpau dumpelspau dumpcohcrs dumpcrs \
	neutroot2 dumptotpauC  basic_histos cross_section_histos \
	neutclass_to_tree

FLAGS      = $(FOPTFLAGS) $(FCOPTIONS) $(FCFLAGS)

all: $(TARGETS)

neut: neut.o $(OBJS)
	$(FC) -o $@ $(FLAGS) neut.o -Xlinker -R`pwd`  $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

neut_ntpl: neut_ntpl.o $(OBJS)
	$(FC) -o $@ $(FLAGS) neut_ntpl.o -Xlinker -R`pwd`  $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

dumptotpau: dumptotpau.o $(OBJS)
	$(FC) -o $@ $(FLAGS) dumptotpau.o -Xlinker -R`pwd`  $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

dumpelspau: dumpelspau.o $(OBJS)
	$(FC) -o $@ $(FLAGS) dumpelspau.o -Xlinker -R`pwd`  $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

dumpcohcrs: dumpcohcrs.o $(OBJS)
	$(FC) -o $@ $(FLAGS) dumpcohcrs.o -Xlinker -R`pwd`  $(OBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

dumpcrs: dumpcrs.o $(OBJS) fnrawcrs.o
	$(FC) -o $@ $(FLAGS) dumpcrs.o -Xlinker -R`pwd`  $(OBJS) fnrawcrs.o \
				$(NEUTLIBS) $(ROOTLIBS) 

neutroot2: neutroot.o $(OBJS)
	$(FC) -o $@ $(FLAGS) neutroot.o -Xlinker -R`pwd` $(OBJS) \
				$(SOBJS) $(NEUTLIBS) $(ROOTLIBS) 

dumptotpauC: dumptotpauC.o $(OBJS) $(SOBJS)
	$(FC) -o $@ $(FLAGS) dumptotpauC.o -Xlinker -R`pwd` ${OBJS} \
				$(SOBJS) $(NEUTLIBS) $(ROOTLIBS) 

basic_histos: basic_histos.o $(OBJS)
	$(FC) -o $@ $(FLAGS) basic_histos.o -Xlinker -R`pwd` $(OBJS) \
				$(SOBJS) $(NEUTLIBS) $(ROOTLIBS) 

neutclass_to_tree: neutclass_to_tree.o $(OBJS)
	$(FC) -o $@ $(FLAGS) neutclass_to_tree.o -Xlinker -R`pwd` $(OBJS) \
				$(SOBJS) $(NEUTLIBS) $(ROOTLIBS) 

cross_section_histos: cross_section_histos.o $(OBJS) fnrawcrs.o
	$(FC) -o $@ $(FLAGS) cross_section_histos.o -Xlinker -R`pwd` $(OBJS) \
				fnrawcrs.o $(SOBJS) $(NEUTLIBS) $(ROOTLIBS) 

clean:
	$(RM) -f *.o *.a $(TARGETS)
