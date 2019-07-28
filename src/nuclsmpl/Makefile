ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

POBJS = $(NEUT_SRC)/neutsmpl/grndmq.o $(NEUT_SRC)/neutsmpl/rmarin_dum.o

MOBJS =  piscatana.o

OBJS  = ${MOBJS} ${POBJS} tauvect.o 

all: nuclscat run_nuclscatana run_nucltransana

nuclscat: nuclscat.o $(POBJS) tauvect.o
	$(FC) -o $@ nuclscat.o -Xlinker -R`pwd` tauvect.o $(POBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

run_nuclscatana: run_nuclscatana.o nuclscatana.o $(POBJS) tauvect.o
	$(FC) -o $@ run_nuclscatana.o nuclscatana.o \
				-Xlinker -R`pwd` tauvect.o $(POBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

run_nucltransana: run_nucltransana.o nucltransana.o $(POBJS) tauvect.o
	$(FC) -o $@ run_nucltransana.o nucltransana.o \
				-Xlinker -R`pwd` tauvect.o $(POBJS) \
				$(NEUTLIBS) $(ROOTLIBS) 

clean:
	$(RM) -f *.o *~ ${MOBJS}  nuclscat run_nuclscatana run_nucltransana

