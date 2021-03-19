include $(top_srcdir)/build.common

AM_FFLAGS = -DNECORE_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)
AM_CXXFLAGS = -DNECORE_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS)
# automake determines that a few of the c files are pure C and thus wont use the CXX flags.
AM_CFLAGS = -DNECORE_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CFLAGS)

LIBNAME = lib@LIBPREFIX@neutcore.@NEUT_VERSION_NUM@.a

if NEUT_LIBCOMPAT
lib_LIBRARIES = $(LIBNAME)
else
noinst_LIBRARIES = $(LIBNAME) libNEUTCERNOverride.a
endif

lib@LIBPREFIX@neutcore_@NEUT_VERSION_NUM@_a_SOURCES = bbba05.F fnelspau_sf.F rsifl6g.F  fnelspf.F neelsvct.F qgdisg.F rsifl7.F cohcrsneberger.F fnels_rpascl.F neelsvnp.F rnazi.F rslzbte.F cohcrsne.F fnf2.F neelsvnt.F rnelq2.F rslzbt.F cohcrsnekart.F fngamfct.F neffect.F rnmass.F rslzbtg.F cohcrsnerein.F nefillmodel.F rnncelq2.F rslzbtk.F cohkartlamda.F nefillver.F rnrot.F rspauli.F cohkartq2min.F fnncelscrs.F nefrmmom.F rnw1320.F rsrotvec.F cohkarttmin.F fnncq2max.F negamnuc.F rnw13.F rssetcrs.F cohkartwmax.F fnncq2min.F negenseed.F rnw.F rssetvct.F cohma.F fnnpot.F negetr.F rnxy.F rssgcrs.F dcohmxneberger.F fnpkno2.F nejet.F rs3drt.F rssgvct.F dcohmxnekart.F fnpkno.F nek0.F rsasum.F rsskcrs.F dcohmxnerein.F fnq2max.F nekindpa.F rsasumg.F rsskvct.F dcohneberger.F fnq2min.F neknopar.F rsclfme.F rsspicrs.F dcohnekart.F fnrworld.F nelibver.F rsclfm.F rsspivct.F dcohnerein.F fntotpau.F rsclfmg.F fnxf3.F rsclfmk.F ddifmx.F rsdcrse.F rstcrset.F difcrsne.F fvfrmgev.F nemkgam.F rsdcrs.F rstcrsetm.F difcrsnerein.F rsdcrsg.F rstcrsetrd.F dndxdy.F nemodgam.F rsdcrsk.F rstcrsett.F dnela.F lenstr.F nemodsel.F rsdiste.F rstcrs.F nemultip.F rsdist.F rstcrsgamrd.F dnelsnl.F neapivct.F nencelslpv.F rsdistg.F rstcrsg.F dnelsq2.F necard.F nencelsspff.F rsdistk.F rstcrsgm.F dnncela.F nechadir.F nencelsvcp.F rsdt33e.F rstcrsgt.F dnncelsq2.F necharge.F rsdt33.F rstcrsk.F elaxff.F nechkfb.F nencelsvnp.F rsdt33g.F rstcrskm.F evdifcrs.F neclpos.F nencelsvnt.F rsdt33k.F rstcrskrd.F fnapicrs.F neclrcrs.F neomega.F rsifl1e.F rstcrskt.F fndsigmapidt.F neclrwrk.F rsifl1.F rstcrsm.F fnelscrn.F necohvct.F neranbll.F rsifl1k.F rstcrsrd.F fnelscrs.F necpnewk.F nerdseed.F rsifl2e.F rstcrst.F fnels_erpascl.F necrosel.F nesetfgparams.F rsifl2.F fnelsint.F necrosin.F nesettarg.F rsifl2g.F fnelspan_effsf.F necrosto.F nesortdl.F rsifl2k.F fnelspan.F necrsmode.F nevect.F rsifl3e.F fnelspan_intpol.F nedcy2bd.F nevent.F rsifl3.F fnelspan_sf.F nedifvct.F nrsettarg.F rsifl3k.F fnelspau_effsf.F neelslpv.F nudeex.F rsifl4.F  neelslvc.F nudeex_n.F rsifl4g.F fnelspau_intpol.F neelsspf.F nudeex_org.F rsifl5.F fnelspau_select.F neelsspff.F nudeex_p.F rsifl6.F structg.f bbba07.cc fa3comp.cc fazexp.cc ddifnerein.c  qedifcrs.c resolvecardname.F mkadler.F mkdcrs.F mkdt33.F res_q2_limits.F rsy.F mk_imode1.cc mk_imode2.cc mk_imode3.cc mk_imode4.cc mk_imode5.cc mk_imode6.cc mk_imode7.cc mk_imode11.cc mk_imode12.cc mk_imode13.cc mk_imode14.cc mk_imode15.cc mk_imode16.cc mk_imode17.cc necard_toml.cc cardistoml.cc bykvec.F fnr1998.F byabpar.F bykax.F

# Only build these sources into neutcore for compatibility, otherwise they are built in the mec directory.
if NEUT_LIBCOMPAT

# Copy these to the build dir so that the object files don't clutter the src tree.
nemechad.F: $(top_srcdir)/mec/nemechad.F
	cp $^ $@

nemeclpv.F: $(top_srcdir)/mec/nemeclpv.F
	cp $^ $@

nemecvcp.F: $(top_srcdir)/mec/nemecvcp.F
	cp $^ $@

neranbll2.F: $(top_srcdir)/mec/neranbll2.F
	cp $^ $@

fnmec_nieves.F: $(top_srcdir)/mec/fnmec_nieves.F
	cp $^ $@

fnmec_select.F: $(top_srcdir)/mec/fnmec_select.F
	cp $^ $@

lib@LIBPREFIX@neutcore_@NEUT_VERSION_NUM@_a_SOURCES += nemechad.F nemeclpv.F nemecvcp.F neranbll2.F fnmec_nieves.F fnmec_select.F

lib@LIBPREFIX@neutcore_@NEUT_VERSION_NUM@_a_SOURCES += pdfset.F structm.F grv94di.F grv98_lo.F

else

#CERN Override
libNEUTCERNOverride_a_SOURCES = pdfset.F structm.F grv94di.F grv98_lo.F

endif

if NEUT_USE_ZBS

lib@LIBPREFIX@neutcore_@NEUT_VERSION_NUM@_a_SOURCES += nemknpos.F

endif

# Unused sources that live in the directory
UNUSED_srcs = calcqecrs.F rs.tcrsdump.F nencelsvct.F fnelspau.F rs.wdist.F rs.spivcttest.F ddifcrs.F neelsvcp.F dnels.F

if NEUT_LIBCOMPAT

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
