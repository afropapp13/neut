## Check default NEUT_ROOT directory
ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif

include ${NEUT_ROOT}/src/neutsmpl/config_neut.gmk

FCFLAGS += -DNECOREVER=${INEUTCOREVER}

LIBBASE	= neutcore
LIBNAME = $(LIBBASE)_$(NEUTCOREVER)

OBJ =	dnela.o fnelscrn.o fnelscrs.o dnelsq2.o fnq2max.o  \
		fnq2min.o fnncq2max.o fnncq2min.o neelsvct.o rnelq2.o rnrot.o \
		rnazi.o neelsspf.o dnels.o fnelspf.o neelsvcp.o fvfrmgev.o \
		neelslvc.o nechadir.o fnnpot.o fnelsint.o fnelspau_select.o \
        nemecvcp.o nemeclpv.o fnmec_select.o fnmec_nieves.o nemechad.o \
        fnelspan_intpol.o neranbll2.o fnels_rpascl.o fnels_erpascl.o\
        fnelspau_sf.o fnelspan_sf.o fnelspau_intpol.o fntotpau.o \
		bbba05.o LH.o fourVector.o fnelspau_effsf.o fnelspan_effsf.o \
		fnncelscrs.o dnncela.o dnncelsq2.o \
		nencelsvnt.o rnncelq2.o dnncelsq2.o \
		nencelsvnp.o nencelslpv.o nencelsspff.o \
		neelsvnp.o neelsvnt.o neclrwrk.o nevect.o nevent.o nemodsel.o \
		fnelspan.o neapivct.o fnapicrs.o rnw.o rnw13.o rnw1320.o rnxy.o \
		dndxdy.o fnf2.o fnxf3.o nedcy2bd.o nekindpa.o necharge.o \
		rnmass.o nemultip.o fnpkno.o fngamfct.o nechkfb.o neffect.o \
		necpnewk.o rsdt33.o nejet.o structg.o qgdisg.o \
		neelslpv.o neelsspff.o dnelsnl.o neomega.o \
		nefrmmom.o nemknpos.o neclpos.o neclrcrs.o \
		fnrworld.o fnr1998.o byabpar.o bykvec.o bykax.o \
		grv94di.o grv98_lo.o pdfset.o structm.o \
		nesettarg.o nrsettarg.o nesetfgparams.o \
        bbba07.o ddifcrs.o evdifcrs.o qedifcrs.o cohma.o necrsmode.o elaxff.o \
		fa3comp.o fazexp.o \
		difcrsne.o ddifmx.o difcrsnerein.o ddifnerein.o nedifvct.o

RSO =	rs3drt.o rsasum.o rsdcrs.o rsclfm.o  rsdist.o rsifl1.o rsifl2.o \
		rsifl3.o rsifl4.o rsifl5.o rsifl6.o rsifl7.o rslzbt.o rsrotvec.o \
		rsspicrs.o rsspivct.o rstcrs.o rstcrsm.o rstcrst.o \
	    rsy00.o rsy20.o rsy21.o rsy22.o \
	    rsclfme.o rsclfmk.o rsdcrse.o rsdcrsk.o rsdiste.o rsdistk.o rsdt33e.o \
	    rsdt33k.o rsifl1e.o rsifl1k.o rsifl2e.o rsifl2k.o rsifl3e.o rsifl3k.o \
	    rslzbte.o rslzbtk.o rssetcrs.o rssetvct.o rsskcrs.o rsskvct.o \
		rstcrset.o rstcrsetm.o rstcrsett.o rstcrsk.o rstcrskm.o rstcrskt.o \
		rspauli.o \
		rsasumg.o  rsclfmg.o  rsdcrsg.o  rsdistg.o  rsdt33g.o \
		rsifl2g.o  rsifl4g.o  rsifl6g.o  rslzbtg.o  rstcrsg.o \
		rssgvct.o  rssgcrs.o  rstcrsgm.o rstcrsgt.o \
		rstcrsetrd.o  rstcrsgamrd.o  rstcrskrd.o  rstcrsrd.o

NEO =	necard.o necohvct.o cohcrsne.o cohcrsnerein.o cohcrsnekart.o \
		dcohnerein.o dcohnekart.o dcohmxnerein.o dcohmxnekart.o \
		cohkartlamda.o cohkarttmin.o cohkartwmax.o cohkartq2min.o \
		fndsigmapidt.o \
		negetr.o necrosin.o necrosto.o \
		nemkgam.o nemodgam.o negamnuc.o \
		nesortdl.o nek0.o neknopar.o fnpkno2.o \
		nudeex.o nudeex_p.o nudeex_n.o \
		negenseed.o nerdseed.o neranbll.o nefillmodel.o nelibver.o \
		nefillver.o lenstr.o dcohneberger.o dcohmxneberger.o necrosel.o \
		cohcrsneberger.o

OBJS = $(OBJ) $(RSO) $(NEO)

INCS = rscons.h mcgenpar.h posinnuc.h neutparams.h neutmodel.h \
	   necard.h nework.h   neutcrs.h  nefillver.h  neutfilepath.h \

all: lib${LIBNAME}.a

lib${LIBNAME}.a: ${OBJS}
	${AR}     rv $@ ${OBJS}
	${RANLIB} $@

install.include: ${INCS}
	for i in $? ;do ${FH2H} $$i > `echo $$i | sed "s/\.h$$/C.h/"`; done
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS} ${NEUT_ROOT}/include
	${INSTALL} -m 0644 ${INCS:.h=C.h} ${NEUT_ROOT}/include
	${INSTALL} -m 0644 neutmodesCPP.h ${NEUT_ROOT}/include

install.library: lib${LIBNAME}.a
	${INSTALL} -d -m 0755 ${NEUT_ROOT}/lib
	${INSTALL} -m 0644 $? ${NEUT_ROOT}/lib




