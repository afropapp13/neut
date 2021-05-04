#include "CascadePion.h"

#include "CommonBlockIFace.h"

#include "vcworkC.h"

#include "posinnucC.h"

#include <cmath>
#include <iostream>

//#define _N_REWEIGHT_CASC_DEBUG_

namespace neut {
namespace rew {

CascadePionEngine::CascadePionEngine() {
  PionFSI_AbsProb =
      RegisterDial("PionFSI_AbsProb", neffpr_.fefabs, 0.432, 0.432);
  PionFSI_QELowMomProb =
      RegisterDial("PionFSI_QELowMomProb", neffpr_.fefqe, 0.313, 0.313);
  PionFSI_CExLowMomProb =
      RegisterDial("PionFSI_CExLowMomProb", neffpr_.fefcx, 0.305, 0.305);
  PionFSI_QEHighMomProb =
      RegisterDial("PionFSI_QEHighMomProb", neffpr_.fefqeh, 0.859, 0.859);
  PionFSI_CExHighMomProb =
      RegisterDial("PionFSI_CExHighMomProb", neffpr_.fefcxh);
  PionFSI_InelProb =
      RegisterDial("PionFSI_InelProb", neffpr_.fefinel, 1.101, 1.101);
}

void CascadePionEngine::Reconfigure(void) {
  ApplyTweaks();
  for (auto &dial : fDials) {
    if (dial.second.ToValue < 0) {
      std::cout << "[ERROR]: Dial: " << NSyst::DialAsString(dial.first)
                << " set unphysical value: " << dial.second.ToValue
                << " (CV: " << dial.second.FromValue
                << ", Uncert: " << dial.second.GetOneSigma(dial.second.Tweak)
                << ", Tweak: " << dial.second.Tweak << ")" << std::endl;
      abort();
    }
  }
}

double CascadePionEngine::CalcWeight() {

  // Not bound
  if (!posinnuc_.ibound) {
    return 1;
  }

  if (!AnyTweaked()) {
    return 1;
  }

  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();
  cbfa.ResetGenValues();

  double old_xsec = NEUTGetPiCascProb();

  if (old_xsec <= 0) {
    std::cout
        << "[WARN]: CascadePionEngine() evpiprob old_xsec <= 0, returning "
           "weight = 1"
        << std::endl;
    return 1;
  }

  neffpr_.fefabs = fDials[PionFSI_AbsProb].ToValue;

  neffpr_.fefqe = fDials[PionFSI_QELowMomProb].ToValue;
  neffpr_.fefqeh = fDials[PionFSI_QEHighMomProb].ToValue;

  neffpr_.fefcx = fDials[PionFSI_CExLowMomProb].ToValue;
  neffpr_.fefcxh = fDials[PionFSI_CExHighMomProb].ToValue;

  neffpr_.fefinel = fDials[PionFSI_InelProb].ToValue;

  NEUTSetParams();

  double new_xsec = NEUTGetPiCascProb();

#ifdef _N_REWEIGHT_CASC_DEBUG_
  std::cout << "pion cascade probability (old) = " << old_xsec << std::endl;
  std::cout << "pion cascade probability (new) = " << new_xsec << std::endl;
#endif

  return CheckReturnWeight(new_xsec / old_xsec);
}

} // namespace rew
} // namespace neut
