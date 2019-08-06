include $(top_srcdir)/build.common

AM_FFLAGS = -DNECORE_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DNECORE_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)

LIBNAME = lib@LIBPREFIX@neutcore.@NEUT_VERSION_NUM@.a

if NEUT_LIBCOMPAT
lib_LIBRARIES = $(LIBNAME)
else
noinst_LIBRARIES = $(LIBNAME)
endif

lib@LIBPREFIX@neutcore_@NEUT_VERSION_NUM@_a_SOURCES = bbba05.F fnelspau_sf.F neelsvcp.F pdfset.F rsifl6g.F  fnelspf.F neelsvct.F qgdisg.F rsifl7.F cohcrsneberger.F fnels_rpascl.F neelsvnp.F rnazi.F rslzbte.F cohcrsne.F fnf2.F neelsvnt.F rnelq2.F rslzbt.F cohcrsnekart.F fngamfct.F neffect.F rnmass.F rslzbtg.F cohcrsnerein.F fnmec_nieves.F nefillmodel.F rnncelq2.F rslzbtk.F cohkartlamda.F fnmec_select.F nefillver.F rnrot.F rspauli.F cohkartq2min.F fnncelscrs.F nefrmmom.F rnw1320.F rsrotvec.F cohkarttmin.F fnncq2max.F negamnuc.F rnw13.F rssetcrs.F cohkartwmax.F fnncq2min.F negenseed.F rnw.F rssetvct.F cohma.F fnnpot.F negetr.F rnxy.F rssgcrs.F dcohmxneberger.F fnpkno2.F nejet.F rs3drt.F rssgvct.F dcohmxnekart.F fnpkno.F nek0.F rsasum.F rsskcrs.F dcohmxnerein.F fnq2max.F nekindpa.F rsasumg.F rsskvct.F dcohneberger.F fnq2min.F neknopar.F rsclfme.F rsspicrs.F dcohnekart.F fnrworld.F nelibver.F rsclfm.F rsspivct.F dcohnerein.F fntotpau.F nemechad.F rsclfmg.F rs.spivcttest.F ddifcrs.F fnxf3.F nemeclpv.F rsclfmk.F ddifmx.F fourVector.F nemecvcp.F rsdcrse.F rstcrset.F difcrsne.F fvfrmgev.F nemkgam.F rsdcrs.F rstcrsetm.F difcrsnerein.F grv94di.F nemknpos.F rsdcrsg.F rstcrsetrd.F dndxdy.F grv98_lo.F nemodgam.F rsdcrsk.F rstcrsett.F dnela.F lenstr.F nemodsel.F rsdiste.F rstcrs.F dnels.F nemultip.F rsdist.F rstcrsgamrd.F dnelsnl.F neapivct.F nencelslpv.F rsdistg.F rstcrsg.F dnelsq2.F necard.F nencelsspff.F rsdistk.F rstcrsgm.F dnncela.F nechadir.F nencelsvcp.F rsdt33e.F rstcrsgt.F dnncelsq2.F necharge.F rsdt33.F rstcrsk.F elaxff.F nechkfb.F nencelsvnp.F rsdt33g.F rstcrskm.F evdifcrs.F neclpos.F nencelsvnt.F rsdt33k.F rstcrskrd.F fnapicrs.F neclrcrs.F neomega.F rsifl1e.F rstcrskt.F fndsigmapidt.F neclrwrk.F neranbll2.F rsifl1.F rstcrsm.F fnelscrn.F necohvct.F neranbll.F rsifl1k.F rstcrsrd.F fnelscrs.F necpnewk.F nerdseed.F rsifl2e.F rstcrst.F fnels_erpascl.F necrosel.F nesetfgparams.F rsifl2.F fnelsint.F necrosin.F nesettarg.F rsifl2g.F rsy00.F fnelspan_effsf.F necrosto.F nesortdl.F rsifl2k.F rsy20.F fnelspan.F necrsmode.F nevect.F rsifl3e.F rsy21.F fnelspan_intpol.F nedcy2bd.F nevent.F rsifl3.F rsy22.F fnelspan_sf.F nedifvct.F nrsettarg.F rsifl3k.F structm.F fnelspau_effsf.F neelslpv.F nudeex.F rsifl4.F  neelslvc.F nudeex_n.F rsifl4g.F fnelspau_intpol.F neelsspf.F nudeex_org.F rsifl5.F fnelspau_select.F neelsspff.F nudeex_p.F rsifl6.F structg.f bbba07.cc fa3comp.cc fazexp.cc ddifnerein.c  qedifcrs.c

# Unused sources that live in the directory
UNUSED_srcs = calcqecrs.F rs.tcrsdump.F nencelsvct.F fnelspau.F rs.wdist.F

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
