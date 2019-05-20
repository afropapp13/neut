ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk

POBJS = $(NEUT_SRC)/neutsmpl/grndmq.o $(NEUT_SRC)/neutsmpl/rmarin_dum.o

MOBJS =  piscatana.o

OBJS  = ${MOBJS} ${POBJS} tauvect.o 

all: piscat gampi run_piscatana run_piscatplot

piscat: piscat.o $(POBJS)
	$(FC) -o $@ piscat.o -Xlinker -R`pwd` $(POBJS) $(NEUTLIBS) $(ROOTLIBS) 

gampi: gampi.o $(POBJS)
	$(FC) -o $@ gampi.o -Xlinker -R`pwd` $(POBJS)  $(NEUTLIBS) $(ROOTLIBS) 

run_piscatana: run_piscatana.o $(OBJS)
	$(FC) -o $@  $(OBJS) run_piscatana.o$(NEUTLIBS) $(ROOTLIBS) 

run_piscatplot: run_piscatplot.o $(OBJS)
	$(FC) -o $@ $(OBJS) run_piscatplot.o $(NEUTLIBS) $(ROOTLIBS) 

clean:
	$(RM) -f *.o *~ ${MOBJS}  run_piscatana piscat  gampi run_piscatplot
