FC        = gfortran
AR        = /bin/ar
RANLIB    = /bin/ranlib
LIBDIR    = ${NEUT_ROOT}/lib


FCOPTIONS = -fno-backslash -fno-automatic -fno-second-underscore \
            -ffixed-line-length-132 -fdollar-ok -fPIC -fbounds-check \
			-Wall -W -ffpe-trap=invalid,zero,overflow \
			-DstrictF77 -DgFortran -Df2cFortran 
FCOPTFLAGS= -g -fPIC

.F.o:
	$(FC) $(DEFS) -c $(FCOPTFLAGS) $(FCOPTIONS) $(FCFLAGS) \
		   $(NEUTINCDIR) $(ROOTINCDIR) $(CERNINCDIR) -o $@ $<

.f.o:
	$(FC) $(DEFS) -c $(FCOPTFLAGS) $(FCOPTIONS) $(FCFLAGS) \
		   $(NEUTINCDIR) $(ROOTINCDIR) $(CERNINCDIR) -o $@ $<

#
#  Objects
#

KZO =	kzaddn.o kzbcr0.o kzbcr1.o kzbcr2.o kzbcr3.o kzbcre.o \
		kzbdel.o kzbloc.o kzbtyp.o kzchg0.o kzcldt.o kzcmpn.o \
		kzcrb0.o kzcrb1.o kzcrb2.o kzcrb3.o kzcrd0.o kzcrd1.o \
		kzcrd2.o kzcrd3.o kzcreh.o kzcrem.o	kzcrmh.o kzddiv.o \
		kzdevt.o kzdsev.o kzdsto.o kzdtab.o kzeclr.o kzecur.o \
		kzefst.o kzejmp.o kzelst.o kzemrg.o kzemrk.o kzend.o \
		kzenew.o kzenum.o kzenxt.o kzeprv.o kzexb0.o kzexb1.o \
		kzexb2.o kzexdi.o kzexdr.o kzfile.o kzform.o kzget.o \
		kzget0.o kzget1.o kzget2.o kzgeth.o kzinit.o kziofm.o \
		kzleng.o kzlimit.o kzlseg.o kzldat.o kzlogl.o kzmcur.o \
		kzmfst.o kzmjmp.o kzmkdt.o kzmlst.o kzmnum.o kzmnxt.o \
		kzmprv.o kzname.o kznseg.o kznsg0.o kznssg.o kzoptn.o \
		kzpdir.o kzput0.o kzput1.o kzput2.o kzputh.o kzread.o \
		kzrep0.o kzrep1.o kzrep2.o kzwpos.o kzwps2.o kzwrit.o \
		kzctoi.o kzitoc.o kzreda.o msgout.o

OBJS = $(KZO)

LIBNAME = zbs

#
#  Rules for building library 
#

.PHONY:  lib$(LIBNAME).a $(LIBDIR)lib$(LIBNAME).a

lib$(LIBNAME).a : $(OBJS)
	$(RM) $@
	$(AR) rv $@ $(OBJS) 
	$(RANLIB) $@

$(LIBDIR)lib$(LIBNAME).a : lib$(LIBNAME).a
	$(RM) $@
	$(INSTALL_LIB) $< $@

#
#  Targets 
#

.PHONY:  clean setup includes install.includes depend lib install.lib exec install.exec

emptyrule:: lib

clean::
	$(RM) *.o *~ *.a *.so core* 

setup::

includes::

install.includes::
	$(INSTALL_INC) zbsdeclr.h $(FORTRAN_INCDIR)
	$(INSTALL_INC) zbsinf.h   $(FORTRAN_INCDIR)
	$(INSTALL_INC) zbssys.h   $(FORTRAN_INCDIR)
	$(INSTALL_INC) zbscf.h    $(INCDIR)

depend::

lib:: lib$(LIBNAME).a

install.lib:: $(LIBDIR)lib$(LIBNAME).a

exec::

install.exec:: 
