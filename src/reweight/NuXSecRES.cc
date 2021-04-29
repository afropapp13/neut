#include "NuXSecRES.h"

#include "NModeDefn.h"

#include "CommonBlockIFace.h"

#include "neworkC.h"

#include <iostream>

//#define _N_REWEIGHT_RES_DEBUG_

namespace neut {
namespace rew {

NuXSecRESEngine::NuXSecRESEngine() {
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  // From nefillmodel.F
  if (cbfa.fneut1pi_gen.neiff != 0) {
    MaRES = RegisterDial("MaRES", cbfa.fneut1pi_gen.xmanffres, 0.15, 0.15);
    MvRES = RegisterDial("MvRES", cbfa.fneut1pi_gen.xmvnffres);
  } else {
    MaRES = RegisterDial("MaRES", cbfa.fneut1pi_gen.xmarsres);
    MvRES = RegisterDial("MvRES", cbfa.fneut1pi_gen.xmvrsres);
  }

  CA5RES = RegisterDial("CA5RES", cbfa.fneut1pi_gen.rneca5i, 0.15, 0.15);
  BgSclRES = RegisterDial("BgSclRES", cbfa.fneut1pi_gen.rnebgscl, 0.15, 0.15);
  MDLSPiEj = RegisterDial("MDLSPiEj", cbfa.fnemdls_gen.mdlspiej);

  UseSeparateBgSclLMCPiBar = RegisterDial("UseSeparateBgSclLMCPiBar", 0);
  BgSclLMCPiBarRES =
      RegisterDial("BgSclLMCPiBarRES", cbfa.fneut1pi_gen.rnebgscl, 1.3, 1.3);
}

void NuXSecRESEngine::Reconfigure(void) { ApplyTweaks(); }

namespace {
bool EventIsLMCPiBar() {
  if (nework_.ipne[0] > 0) { // Is anti-nu
    return false;
  }

  size_t NParts = nework_.numne;

  size_t NFSCPions = 0;
  // Use primary pion (i.e. the first one you see)
  double FirstCPiMom = 0;

  // std::cout << "Mode = " << nework_.modene << ", NParts = " << NParts
  // << std::endl;

  for (size_t p_it = 0; p_it < NParts; ++p_it) {
    int PID = nework_.ipne[p_it];
    if (abs(PID) != 211) {
      continue;
    }
    NFSCPions++;

    if (NFSCPions == 1) {
      // vcwork is in MeV
      FirstCPiMom = sqrt(nework_.pne[p_it][0] * nework_.pne[p_it][0] +
                         nework_.pne[p_it][1] * nework_.pne[p_it][1] +
                         nework_.pne[p_it][2] * nework_.pne[p_it][2]) *
                    1E-3;
    }
  }
  return (NFSCPions > 0) && (FirstCPiMom < 0.2);
}
} // namespace

double NuXSecRESEngine::CalcWeight() {
  if (!NModeDefn::isRES(nework_.modene)) {
    return 1;
  }

  if (!AnyTweaked()) {
    return 1;
  }

  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();
  cbfa.ResetGenValues();

  double old_xsec = NEUTGetXSec();
  double piej_old = NModeDefn::is1PI(nework_.modene) ? NEUTGetPiEj() : 1;

  if (old_xsec == 0) {
    std::cout << "NuXSecRESEngine::CalcWeight() Warning: old_xsec==0, setting "
                 "weight to 1"
              << std::endl;
    return 1;
  }

  if (piej_old == 0) { 
    std::cout << "NuXSecRESEngine::CalcWeight() Warning: piej_old==0, setting "
                 "piej_old to 1" << std::endl;
    piej_old = 1;
  }

  if (cbfa.fneut1pi_gen.neiff != 0) {
    neut1pi_.xmanffres = fDials[MaRES].ToValue;
    neut1pi_.xmvnffres = fDials[MvRES].ToValue;
  } else {
    neut1pi_.xmarsres = fDials[MaRES].ToValue;
    neut1pi_.xmvrsres = fDials[MvRES].ToValue;
  }

  neut1pi_.rneca5i = fDials[CA5RES].ToValue;
  nemdls_.mdlspiej = fDials[MDLSPiEj].ToValue;

  if (fDials[UseSeparateBgSclLMCPiBar].ToValue && EventIsLMCPiBar()) {
    neut1pi_.rnebgscl = fDials[BgSclLMCPiBarRES].ToValue;
  } else {
    neut1pi_.rnebgscl = fDials[BgSclRES].ToValue;
  }

  NEUTSetParams();

  double new_xsec = NEUTGetXSec();
  double piej_new = NModeDefn::is1PI(nework_.modene) ? NEUTGetPiEj() : 1;

  if (piej_new == 0) { 
    std::cout << "NuXSecRESEngine::CalcWeight() Warning: piej_new==0" << std::endl;
  }

#ifdef _N_REWEIGHT_RES_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << endl;
  cout << "differential cross section (new) = " << new_xsec << endl;
#endif

  return CheckReturnWeight(new_xsec / old_xsec) *
         CheckReturnWeight(piej_new / piej_old);
}

} // namespace rew
} // namespace neut
