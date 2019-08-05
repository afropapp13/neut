include $(top_srcdir)/build.common

AM_FFLAGS = -DNENUCE_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_FFLAGS)

LIBNAME = lib@LIBPREFIX@nuceff.@NEUT_VERSION_MAJ@.@NEUT_VERSION_MIN@.@NEUT_VERSION_REV@.a

lib_LIBRARIES = $(LIBNAME)

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

lib@LIBPREFIX@nuceff_@NEUT_VERSION_MAJ@_@NEUT_VERSION_MIN@_@NEUT_VERSION_REV@_a_SOURCES = absboost.F absdelta.F absflpsor.F absfourbody.F absgen_four.F absgen_three.F absgen_two.F absnegenbod.F absneject.F absnepdk.F abspidpp.F abspi.F absprob.F absrotres2.F absthreebody.F abstwobody.F ainvms.F dcy2bd.F dcy3bd.F density.F efabpath.F efabpathtun.F efabrho.F efabsrat.F efcalamp.F efclfsi.F efcohsct.F efconv2roxy.F efcoulforce.F efd13.F efd15.F efdsgamp.F efdsg.F efdsgfn.F efdsgps.F effrmgas.F effrmmom.F effrmmom_l.F effsquar.F efgammn.F efgammnt.F efgammp.F efgammpt.F efgampi0n.F efgampi0nt.F efgampi0p.F efgampi0pt.F efgampi.F efgampimp.F efgampimpt.F efgampipn.F efgampipnt.F efgampiprod.F efhadpro.F efhinelscx.F efhinels_iso.F efhinels_ps.F efinels.F efintegral.F efkinom.F eflibver.F efmakmap.F efomegadcy.F efp11.F efp11mc.F efp13.F efp13mc.F efp31.F efp31mc.F efp33.F efp33mc.F efparq.F efpfcost.F efpncrsh.F efpncrsh_pnsep.F efpnmltp.F efpqeab.F efpqeabh.F efpqeabh_pnsep.F efprobom.F efqepath.F efqepathtun.F efranbll.F efrho.F efrnmom.F efs11.F efs11mc.F efs31.F efs31mc.F efsabrho.F efscatt.F efsetpf.F eftrace.F eftrans.F eftrceta.F eftrcka.F eftrcomega.F efvertra.F efxsctom.F evpiprob.F evsetpf.F fsifillpart.F fsifillvert.F gamfct.F gdeca2.F passl.F rsmass.F sigeta.F typdcy.F


# Unused sources that live in the directory
UNUSED_srcs =
