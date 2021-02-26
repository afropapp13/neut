# Perform some checks
ifndef NEUT_ROOT
	NEUT_ROOT = ../../../
endif

# Just check that NEUT core is built
ifeq "$(wildcard ../../neutcore/Linux_pc)" ""
  $(error Are you sure you have built NEUT? The neutcore library cant be found!))
endif

# Check the env. variables for CERNLIB
ifndef CERN
 $(error CERN is not set!)
endif

ifndef CERN_LEVEL
 $(error CERN_LEVEL is not set!)
endif

# Check that we can find cernlib in shell
ifeq "$(shell cernlib)" ""
  $(error Looks like you havent set up CERNLIB, please fix this)
  $(error I rely on being able to do cernlib in shell)
endif

# Compile with gfortran
FC = gfortran
# Don't use optimiser flags
#FCOPTIONS = -g -Wall -Wno-unused-parameter -Wno-unused-variable -Wtabs -O3 -ffast-math
FCOPTIONS = -g -Wall -Wno-unused-parameter -Wno-unused-variable -Wtabs -O0

# ROOT includes, don't think we really need this...
ROOTINC = -I$(shell root-config --incdir)

# NEUT includes
NEUTINC = -I./ -I../../neutcore
INCDIRS = ${ROOTINC} ${NEUTINC} ${CLASSINC}

LIBDIRS = -L../../../lib/Linux_pc
NEUTCOREVER	= 5.3.6
INEUTCOREVER= 536
NUCEFFVER   = 5.3.5
INUCEFFVER  = 535
NUCCORVER   = 5.3.5
INUCCORVER  = 535
PARTNUCKVER = 5.3.5
IPARTNUCKVER= 535
SKMCSVCVER  = 5.3.5
ISKMCSVCVER = 535
ZBSFNSVER   = 5.3.5
IZBSFNSVER  = 535
SPECFUNCVER = 5.3.5
ISPECFUNCVER= 535
RADCORRVER  = 5.3.5
IRADCORRVER = 535

# The NEUT libraries; only need to include neutcore
MCLIB = ${LIBDIRS} -lnuceff_${NUCEFFVER} -lneutcore_${NEUTCOREVER} \
        -lskmcsvc_${SKMCSVCVER} -lnuccorrspl_${NUCCORVER} \
	-lradcorr_${RADCORRVER}

# The CERN libraries from CERNLIB
CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
	${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
	`cernlib jetset74 photos202 mathlib packlib kernlib`

LIBS      = ${MCLIB} ${CLIBS} -lstdc++

# Build all the .F files
# NOT IMPLEMENTED BUT LEFT ANYWAY IF ANYONE WANTS TO REWRITE
BUILD = $(wildcard crs*.F)
# All the files in here will need to be built (also used in the rm)
TARGET = $(BUILD:.F=)


all: crs_newff crs_min_adap crs_newff_W_adap crs_newff_Q2_adap crs_newff_adap crs_min_single crs_min_W_adap crs_min_Q2_adap

crs_newff: crs_newff.F
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS)

crs_min_adap: crs_min_adap.F quadpack.o
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS) quadpack.o

crs_min_W_adap: crs_min_W_adap.F quadpack.o
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS) quadpack.o

crs_min_Q2_adap: crs_min_Q2_adap.F quadpack.o
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS) quadpack.o

crs_min_single: crs_min_single.F
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS)

crs_newff_adap: crs_newff_adap.F quadpack.o
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS) quadpack.o

crs_newff_W_adap: crs_newff_W_adap.F quadpack.o
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS) quadpack.o

crs_newff_Q2_adap: crs_newff_Q2_adap.F quadpack.o
	$(FC) $(FCOPTIONS) -o $@ $< $(INCDIRS) -Xlinker -R`pwd` $(LIBS) quadpack.o

%.o:	%.F
	$(FC) $(FCOPTIONS) -c $< $(NEUTINC)

clean:
	$(RM) -v *.o *~ ${TARGET}
