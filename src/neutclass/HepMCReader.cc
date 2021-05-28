#include "HepMCReader.h"

#include "necardC.h"
#include "neutcrsC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"
#include "neworkC.h"
#include "vcvrtxC.h"
#include "vcworkC.h"

#include "fsihistC.h"
#include "nucleonfsihistC.h"

#include "posinnucC.h"

#include "NuHepMC/Print.hxx"
#include "NuHepMC/WriterTools"

#include "HepMC3/Print.h"

template <typename T> inline std::string sstostr(T const &t) {
  std::stringstream ss("");
  ss << t;
  return ss.str();
}

void ReadGenRunInfo(std::shared_ptr<HepMC3::GenRunInfo const> gri) {

  neutcard_common &neutcard = neutcard_;

  NuHepMC::FromAttribute(gri, "neutcard.nefrmflg", neutcard.nefrmflg);
  NuHepMC::FromAttribute(gri, "neutcard.nepauflg", neutcard.nepauflg);
  NuHepMC::FromAttribute(gri, "neutcard.nenefo16", neutcard.nenefo16);
  NuHepMC::FromAttribute(gri, "neutcard.nenefmodl", neutcard.nenefmodl);
  NuHepMC::FromAttribute(gri, "neutcard.nenefmodh", neutcard.nenefmodh);
  NuHepMC::FromAttribute(gri, "neutcard.nenefkinh", neutcard.nenefkinh);
  NuHepMC::FromAttribute(gri, "neutcard.nemodflg", neutcard.nemodflg);
  NuHepMC::FromAttribute(gri, "neutcard.neselmod", neutcard.neselmod);

  std::vector<std::remove_reference<decltype(*neutcard.crsneut)>::type> crsneut;
  NuHepMC::FromAttribute(gri, "neutcard.crsneut", crsneut);
  std::copy_n(crsneut.begin(), 30, neutcard.crsneut);

  std::vector<std::remove_reference<decltype(*neutcard.crsneutb)>::type>
      crsneutb;
  NuHepMC::FromAttribute(gri, "neutcard.crsneutb", crsneutb);
  std::copy_n(crsneutb.begin(), 30, neutcard.crsneutb);

  NuHepMC::FromAttribute(gri, "neutcard.itauflgcore", neutcard.itauflgcore);
  NuHepMC::FromAttribute(gri, "neutcard.nusim", neutcard.nusim);
  NuHepMC::FromAttribute(gri, "neutcard.quiet", neutcard.quiet);

  nuceffver_common &nuceffver = nuceffver_;

  NuHepMC::FromAttribute(gri, "nuceffver.nefkinver", nuceffver.nefkinver);

  neutdis_common &neutdis = neutdis_;

  NuHepMC::FromAttribute(gri, "neutdis.nepdf", neutdis.nepdf);
  NuHepMC::FromAttribute(gri, "neutdis.nebodek", neutdis.nebodek);
  NuHepMC::FromAttribute(gri, "neutdis.nemult", neutdis.nemult);

  neut1pi_common &neut1pi = neut1pi_;

  NuHepMC::FromAttribute(gri, "neut1pi.xmanffres", neut1pi.xmanffres);
  NuHepMC::FromAttribute(gri, "neut1pi.xmvnffres", neut1pi.xmvnffres);
  NuHepMC::FromAttribute(gri, "neut1pi.xmarsres", neut1pi.xmarsres);
  NuHepMC::FromAttribute(gri, "neut1pi.xmvrsres", neut1pi.xmvrsres);
  NuHepMC::FromAttribute(gri, "neut1pi.neiff", neut1pi.neiff);
  NuHepMC::FromAttribute(gri, "neut1pi.nenrtype", neut1pi.nenrtype);
  NuHepMC::FromAttribute(gri, "neut1pi.rneca5i", neut1pi.rneca5i);
  NuHepMC::FromAttribute(gri, "neut1pi.rnebgscl", neut1pi.rnebgscl);

  neutdif_common &neutdif = neutdif_;

  NuHepMC::FromAttribute(gri, "neutdif.nedifpi", neutdif.nedifpi);

  neutcoh_common &neutcoh = neutcoh_;

  NuHepMC::FromAttribute(gri, "neutcoh.necohepi", neutcoh.necohepi);

  neutpiabs_common &neutpiabs = neutpiabs_;

  NuHepMC::FromAttribute(gri, "neutpiabs.neabspiemit", neutpiabs.neabspiemit);

  neutpiless_common &neutpiless = neutpiless_;

  NuHepMC::FromAttribute(gri, "neutpiless.ipilessdcy", neutpiless.ipilessdcy);
  NuHepMC::FromAttribute(gri, "neutpiless.rpilessdcy", neutpiless.rpilessdcy);

  neutradcorr_common &neutradcorr = neutradcorr_;

  NuHepMC::FromAttribute(gri, "neutradcorr.iradcorr", neutradcorr.iradcorr);

  nemdls_common &nemdls = nemdls_;

  NuHepMC::FromAttribute(gri, "nemdls.mdlqe", nemdls.mdlqe);
  NuHepMC::FromAttribute(gri, "nemdls.mdlspi", nemdls.mdlspi);
  NuHepMC::FromAttribute(gri, "nemdls.mdldis", nemdls.mdldis);
  NuHepMC::FromAttribute(gri, "nemdls.mdlcoh", nemdls.mdlcoh);
  NuHepMC::FromAttribute(gri, "nemdls.mdldif", nemdls.mdldif);
  NuHepMC::FromAttribute(gri, "nemdls.mdlqeaf", nemdls.mdlqeaf);
  NuHepMC::FromAttribute(gri, "nemdls.xmaqe", nemdls.xmaqe);
  NuHepMC::FromAttribute(gri, "nemdls.xmaspi", nemdls.xmaspi);
  NuHepMC::FromAttribute(gri, "nemdls.xmvqe", nemdls.xmvqe);
  NuHepMC::FromAttribute(gri, "nemdls.xmvspi", nemdls.xmvspi);
  NuHepMC::FromAttribute(gri, "nemdls.kapp", nemdls.kapp);
  NuHepMC::FromAttribute(gri, "nemdls.xmacoh", nemdls.xmacoh);
  NuHepMC::FromAttribute(gri, "nemdls.rad0nu", nemdls.rad0nu);
  NuHepMC::FromAttribute(gri, "nemdls.fa1coh", nemdls.fa1coh);
  NuHepMC::FromAttribute(gri, "nemdls.fb1coh", nemdls.fb1coh);
  NuHepMC::FromAttribute(gri, "nemdls.iffspi", nemdls.iffspi);
  NuHepMC::FromAttribute(gri, "nemdls.nrtypespi", nemdls.nrtypespi);
  NuHepMC::FromAttribute(gri, "nemdls.rca5ispi", nemdls.rca5ispi);
  NuHepMC::FromAttribute(gri, "nemdls.rbgsclspi", nemdls.rbgsclspi);
  NuHepMC::FromAttribute(gri, "nemdls.xmares", nemdls.xmares);
  NuHepMC::FromAttribute(gri, "nemdls.xmvres", nemdls.xmvres);
  NuHepMC::FromAttribute(gri, "nemdls.sccfv", nemdls.sccfv);
  NuHepMC::FromAttribute(gri, "nemdls.sccfa", nemdls.sccfa);
  NuHepMC::FromAttribute(gri, "nemdls.fpqe", nemdls.fpqe);
  NuHepMC::FromAttribute(gri, "nemdls.pfsf", nemdls.pfsf);
  NuHepMC::FromAttribute(gri, "nemdls.xmadif", nemdls.xmadif);
  NuHepMC::FromAttribute(gri, "nemdls.nucvoldif", nemdls.nucvoldif);
  NuHepMC::FromAttribute(gri, "nemdls.axffalpha", nemdls.axffalpha);
  NuHepMC::FromAttribute(gri, "nemdls.axffgamma", nemdls.axffgamma);
  NuHepMC::FromAttribute(gri, "nemdls.axfftheta", nemdls.axfftheta);
  NuHepMC::FromAttribute(gri, "nemdls.axffbeta", nemdls.axffbeta);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpq4", nemdls.axzexpq4);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpnt", nemdls.axzexpnt);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpt0", nemdls.axzexpt0);
  NuHepMC::FromAttribute(gri, "nemdls.axzexptc", nemdls.axzexptc);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa0", nemdls.axzexpa0);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa1", nemdls.axzexpa1);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa2", nemdls.axzexpa2);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa3", nemdls.axzexpa3);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa4", nemdls.axzexpa4);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa5", nemdls.axzexpa5);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa6", nemdls.axzexpa6);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa7", nemdls.axzexpa7);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa8", nemdls.axzexpa8);
  NuHepMC::FromAttribute(gri, "nemdls.axzexpa9", nemdls.axzexpa9);
  NuHepMC::FromAttribute(gri, "nemdls.mdl2p2h", nemdls.mdl2p2h);
  NuHepMC::FromAttribute(gri, "nemdls.xmancel", nemdls.xmancel);

  nenupr_common &nenupr = nenupr_;

  NuHepMC::FromAttribute(gri, "nenupr.iformlen", nenupr.iformlen);
  NuHepMC::FromAttribute(gri, "nenupr.fzmu2", nenupr.fzmu2);
  NuHepMC::FromAttribute(gri, "nenupr.sfebshift", nenupr.sfebshift);
  NuHepMC::FromAttribute(gri, "nenupr.sfebnegbeh", nenupr.sfebnegbeh);

  neffpr_common &neffpr = neffpr_;

  NuHepMC::FromAttribute(gri, "neffpr.fefqe", neffpr.fefqe);
  NuHepMC::FromAttribute(gri, "neffpr.fefqeh", neffpr.fefqeh);
  NuHepMC::FromAttribute(gri, "neffpr.fefinel", neffpr.fefinel);
  NuHepMC::FromAttribute(gri, "neffpr.fefabs", neffpr.fefabs);
  NuHepMC::FromAttribute(gri, "neffpr.fefcoh", neffpr.fefcoh);
  NuHepMC::FromAttribute(gri, "neffpr.fefqehf", neffpr.fefqehf);
  NuHepMC::FromAttribute(gri, "neffpr.fefcohf", neffpr.fefcohf);
  NuHepMC::FromAttribute(gri, "neffpr.fefcx", neffpr.fefcx);
  NuHepMC::FromAttribute(gri, "neffpr.fefcxhf", neffpr.fefcxhf);
  NuHepMC::FromAttribute(gri, "neffpr.fefcxh", neffpr.fefcxh);
  NuHepMC::FromAttribute(gri, "neffpr.fefcoul", neffpr.fefcoul);
  NuHepMC::FromAttribute(gri, "neffpr.fefall", neffpr.fefall);
}

// let the compiler get the types correct
template <typename T, size_t N>
void FromVectorAttribute(std::shared_ptr<HepMC3::GenEvent> evt,
                         std::string const &name, T (&t)[N], size_t n) {
  std::vector<T> vect;
  NuHepMC::FromAttribute(evt, name, vect);
  std::copy_n(vect.begin(), n, t);
}

template <typename T, size_t N, size_t stride>
void FromVectorAttribute(std::shared_ptr<HepMC3::GenEvent> evt,
                         std::string const &name, T (&t)[N][stride], size_t n) {
  std::vector<T> vect;
  NuHepMC::FromAttribute(evt, name, vect);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < stride; ++j) {
      t[i][j] = vect[(i * stride) + j];
    }
  }
}

template <typename T> std::string VectorString(std::vector<T> const &v) {
  std::stringstream ss("");
  ss << "[ ";
  for (auto const &e : v) {
    ss << v << ", ";
  }

  std::string str = ss.str();
  return str.substr(0, str.length() - 2) + " ]";
}

void ReadNuHepMCEvent(std::shared_ptr<HepMC3::GenEvent> evt) {

  nework_.modene = NuHepMC::genevent::GetHardScatterMode(*evt);

  NuHepMC::FromAttribute(evt, "Totcrs", neutcrscom_.totcrsne);
  NuHepMC::FromAttribute(evt, "CrsEnergy", neutcrscom_.crsenergy);

  FromVectorAttribute(evt, "DifCrsNE", neutcrscom_.difcrsne, 8);

  NuHepMC::FromAttribute(evt, "Crsx", neutcrscom_.crsx);
  NuHepMC::FromAttribute(evt, "Crsy", neutcrscom_.crsy);
  NuHepMC::FromAttribute(evt, "Crsz", neutcrscom_.crsz);
  NuHepMC::FromAttribute(evt, "Crsphi", neutcrscom_.crsphi);
  NuHepMC::FromAttribute(evt, "Crsq2", neutcrscom_.crsq2);

  NuHepMC::FromAttribute(evt, "numbndn", neuttarget_.numbndn);
  NuHepMC::FromAttribute(evt, "numbndp", neuttarget_.numbndp);
  NuHepMC::FromAttribute(evt, "numfrep", neuttarget_.numfrep);
  NuHepMC::FromAttribute(evt, "numatom", neuttarget_.numatom);

  NuHepMC::FromAttribute(evt, "pfsurf", nenupr_.pfsurf);
  NuHepMC::FromAttribute(evt, "pfmax", nenupr_.pfmax);
  NuHepMC::FromAttribute(evt, "vnuini", nenupr_.vnuini);
  NuHepMC::FromAttribute(evt, "vnufin", nenupr_.vnufin);

  NuHepMC::FromAttribute(evt, "ibound", posinnuc_.ibound);

  NuHepMC::FromAttribute(evt, "nvtxvc", vcvrtx_.nvtxvc);
  FromVectorAttribute(evt, "pvtxvc", vcvrtx_.pvtxvc, vcvrtx_.nvtxvc);
  FromVectorAttribute(evt, "iflvvc", vcvrtx_.iflvvc, vcvrtx_.nvtxvc);
  FromVectorAttribute(evt, "iparvc", vcvrtx_.iparvc, vcvrtx_.nvtxvc);
  FromVectorAttribute(evt, "timvvc", vcvrtx_.timvvc, vcvrtx_.nvtxvc);
  NuHepMC::FromAttribute(evt, "nvc", vcwork_.nvc);
  FromVectorAttribute(evt, "posvc", vcwork_.posvc, 3);
  FromVectorAttribute(evt, "ipvc", vcwork_.ipvc, vcwork_.nvc);
  FromVectorAttribute(evt, "amasvc", vcwork_.amasvc, vcwork_.nvc);
  FromVectorAttribute(evt, "pvc", vcwork_.pvc, vcwork_.nvc);
  FromVectorAttribute(evt, "iorgvc", vcwork_.iorgvc, vcwork_.nvc);
  FromVectorAttribute(evt, "iflgvc", vcwork_.iflgvc, vcwork_.nvc);
  FromVectorAttribute(evt, "icrnvc", vcwork_.icrnvc, vcwork_.nvc);
  FromVectorAttribute(evt, "timvc", vcwork_.timvc, vcwork_.nvc);
  FromVectorAttribute(evt, "posivc", vcwork_.posivc, vcwork_.nvc);
  FromVectorAttribute(evt, "ivtivc", vcwork_.ivtivc, vcwork_.nvc);
  FromVectorAttribute(evt, "posfvc", vcwork_.posfvc, vcwork_.nvc);
  FromVectorAttribute(evt, "ivtfvc", vcwork_.ivtfvc, vcwork_.nvc);
  NuHepMC::FromAttribute(evt, "numne", nework_.numne);
  NuHepMC::FromAttribute(evt, "modene", nework_.modene);

  FromVectorAttribute(evt, "ipne", nework_.ipne, nework_.numne);
  FromVectorAttribute(evt, "pne", nework_.pne, nework_.numne);
  FromVectorAttribute(evt, "iorgne", nework_.iorgne, nework_.numne);
  FromVectorAttribute(evt, "iflgne", nework_.iflgne, nework_.numne);
  FromVectorAttribute(evt, "icrnne", nework_.icrnne, nework_.numne);
  NuHepMC::FromAttribute(evt, "nvert", fsihist_.nvert);

  FromVectorAttribute(evt, "posvert", fsihist_.posvert, fsihist_.nvert);
  FromVectorAttribute(evt, "iflgvert", fsihist_.iflgvert, fsihist_.nvert);

  NuHepMC::FromAttribute(evt, "nvcvert", fsihist_.nvcvert);

  FromVectorAttribute(evt, "dirvert", fsihist_.dirvert, fsihist_.nvcvert);
  FromVectorAttribute(evt, "abspvert", fsihist_.abspvert, fsihist_.nvcvert);
  FromVectorAttribute(evt, "abstpvert", fsihist_.abstpvert, fsihist_.nvcvert);
  FromVectorAttribute(evt, "ipvert", fsihist_.ipvert, fsihist_.nvcvert);
  FromVectorAttribute(evt, "iverti", fsihist_.iverti, fsihist_.nvcvert);
  FromVectorAttribute(evt, "ivertf", fsihist_.ivertf, fsihist_.nvcvert);
  NuHepMC::FromAttribute(evt, "nfnvert", nucleonfsihist_.nfnvert);
  FromVectorAttribute(evt, "nfiflag", nucleonfsihist_.nfiflag,
                      nucleonfsihist_.nfnvert);

  FromVectorAttribute(evt, "nfx", nucleonfsihist_.nfx, nucleonfsihist_.nfnvert);
  FromVectorAttribute(evt, "nfy", nucleonfsihist_.nfy, nucleonfsihist_.nfnvert);
  FromVectorAttribute(evt, "nfz", nucleonfsihist_.nfz, nucleonfsihist_.nfnvert);

  FromVectorAttribute(evt, "nfpz", nucleonfsihist_.nfpz,
                      nucleonfsihist_.nfnvert);
  FromVectorAttribute(evt, "nfpy", nucleonfsihist_.nfpy,
                      nucleonfsihist_.nfnvert);
  FromVectorAttribute(evt, "nfpz", nucleonfsihist_.nfpz,
                      nucleonfsihist_.nfnvert);
  FromVectorAttribute(evt, "nfe", nucleonfsihist_.nfe, nucleonfsihist_.nfnvert);
  FromVectorAttribute(evt, "nffirststep", nucleonfsihist_.nffirststep,
                      nucleonfsihist_.nfnvert);

  NuHepMC::FromAttribute(evt, "nfnstep", nucleonfsihist_.nfnstep);

  FromVectorAttribute(evt, "nfecms2", nucleonfsihist_.nfecms2,
                      nucleonfsihist_.nfnstep);
  FromVectorAttribute(evt, "nfptot", nucleonfsihist_.nfptot,
                      nucleonfsihist_.nfnstep);
}
