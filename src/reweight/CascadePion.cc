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
  FrAbs_pi = RegisterDial("FrAbs_pi", neffpr_.fefabs, 0.432, 0.432);
  FrInelLow_pi = RegisterDial("FrInelLow_pi", neffpr_.fefqe, 0.313, 0.313);
  FrCExLow_pi = RegisterDial("FrCExLow_pi", neffpr_.fefcx, 0.305, 0.305);
  FrInelHigh_pi = RegisterDial("FrInelHigh_pi", neffpr_.fefqeh, 0.859, 0.859);
  FrCExHigh_pi = RegisterDial("FrCExHigh_pi", neffpr_.fefcxh);
  FrPiProd_pi = RegisterDial("FrPiProd_pi", neffpr_.fefinel, 1.101, 1.101);
}

void CascadePionEngine::Reconfigure(void) { ApplyTweaks(); }

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

  neffpr_.fefabs = fDials[FrAbs_pi].ToValue;

  neffpr_.fefqe = fDials[FrInelLow_pi].ToValue;
  neffpr_.fefqeh = fDials[FrInelHigh_pi].ToValue;

  neffpr_.fefcx = fDials[FrCExLow_pi].ToValue;
  neffpr_.fefcxh = fDials[FrCExHigh_pi].ToValue;

  neffpr_.fefinel = fDials[FrPiProd_pi].ToValue;

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
