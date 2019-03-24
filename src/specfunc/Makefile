## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

LIBBASE	= specfunc
LIBNAME = ${LIBBASE}_$(SPECFUNCVER)

## For library
OBJS = fourVector.o calcLH.o selectSfValues.o buildSf.o maxDiff.o \
	   readXsecData.o sfevent.o selectEffSFValues.o effsfevent.o  \
	   buildEffSf.o effSFMaxDiff.o readEffSFXsecData.o

## executable objectss
MAKETABLES      = makeTables.o combineTables.o
MAKEEFFSFTABLES = makeEffSFTables.o combineEffSFTables.o

OBJD = tauvect_dummy.o

library: lib${LIBNAME}.a

exec: makeTables makeEffSFTables

makeTables: $(MAKETABLES) $(OBJD)
	$(FC) $(FCOPTIONS) $(FCOPTFLAGS) -o makeTables.exe $(MAKETABLES) $(OBJD) \
			$(NEUTLIBS) $(ROOTLIBS) $(CERNLIBS) 

makeEffSFTables: $(MAKEEFFSFTABLES) $(OBJD)
	$(FC) $(FLAGS) -o makeEffSFTables.exe $(MAKEEFFSFTABLES) $(OBJD) \
			$(NEUTLIBS) $(ROOTLIBS) $(CERNLIBS)

lib${LIBNAME}.a: ${OBJS}
	${AR}     rv $@ ${OBJS}
	${RANLIB} $@

install.include: ${INCS}
	for i in $? ;do ${FH2H} $$i > `echo $$i | sed "s/\.h$$/C.h/"`; done
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS} ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS:.h=C.h} ${NEUT_ROOT}/include

install.library: lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib

