## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

FCFLAGS += -DNENUCEVER=$(INUCEFFVER)

LIBBASE	= nuceff
LIBNAME = $(LIBBASE)_$(NUCEFFVER)

OBJS  = efabpath.o efabrho.o    efabsrat.o efcalamp.o efcohsct.o \
		efd13.o    efd15.o      effrmgas.o effrmmom.o effrmmom_l.o efrnmom.o \
		effsquar.o efgammn.o  efgammnt.o efgammp.o    efgammpt.o efgampi.o \
		efhadpro.o \
		efinels.o  efintegral.o efmakmap.o efp11.o    efp11mc.o  efp13.o \
		efp13mc.o  efp31.o      efp31mc.o  efp33.o    efp33mc.o  efparq.o \
		efpfcost.o efpncrsh.o   efpnmltp.o efpqeab.o  efpqeabh.o \
		efqepath.o efranbll.o   efrho.o    efs11.o    efs11mc.o  efs31.o \
		efs31mc.o  efsabrho.o   efscatt.o  efsetpf.o  eftrace.o  \
	    efcoulforce.o eftrcka.o eftrans.o efvertra.o \
		density.o passl.o ainvms.o eftrceta.o \
		gamfct.o typdcy.o rsmass.o  dcy2bd.o sigeta.o dcy3bd.o \
		eftrcomega.o efomegadcy.o efprobom.o efkinom.o efxsctom.o \
		efconv2roxy.o eflibver.o \
		absboost.o    absgen_four.o  absnegenbod.o abspi.o    absrotres2.o \
		absflpsor.o   absgen_three.o absneject.o   abspidpp.o absthreebody.o \
		absfourbody.o absgen_two.o   absnepdk.o    absprob.o  abstwobody.o \
		absdelta.o \
		efgampi0n.o efgampi0nt.o  efgampimp.o efgampiprod.o \
		efgampi0p.o  efgampimpt.o  efgampi0pt.o  efgampipn.o efgampipnt.o \
		efpqeabh_pnsep.o efpncrsh_pnsep.o efhinels_iso.o gdeca2.o \
		efhinels_ps.o efdsg.o efdsgps.o efdsgamp.o efdsgfn.o efhinelscx.o \
		efabpathtun.o efqepathtun.o \
		fsifillpart.o fsifillvert.o evpiprob.o evsetpf.o efclfsi.o \
		efcoulforce.o

INCS = efpion.h efomega.h efpiinth.h fsihist.h efcoulforce.h

all: lib${LIBNAME}.a

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

