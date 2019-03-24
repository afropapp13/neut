## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

## For library
SOBJS = neutctrl.so neutpart.so neutvect.so neutvtx.so \
        neutfsipart.so neutfsivert.so \
	neutnucfsistep.so neutnucfsivert.so \
	neutrootTreeSingleton.so

ROBJS = neutfill.o neutrootfortwrapper.o NeutRootHandlers.o

OBJS  = ${ROBJS} ${SOBJS} 

INCS  = NeutRootHandlers.h  neutnucfsistep.h	     neutvect.h    \
		neutctrl.h	        neutnucfsivert.h	     neutvtx.h     \
		neutfsipart.h	    neutpart.h \
		neutfsivert.h	    neutrootTreeSingleton.h

all: $(OBJS)

install.include: ${INCS}
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS} ${NEUT_ROOT}/include

install.library:
	echo "nothing to do here."
