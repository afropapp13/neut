#include "NuXSecCCQE.h"

#include "NModeDefn.h"

#include "CommonBlockIFace.h"

#include "neutcrsC.h"
#include "neworkC.h"

#include <iostream>

namespace neut {
namespace rew {

NuXSecCCQEEngine::NuXSecCCQEEngine() { // Get the parameter values at generation
                                       // time and store them as the 'default'
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  MaCCQE = RegisterDial("MaCCQE", cbfa.fnemdls_gen.xmaqe, 0.196,0.196);
  AxlFFCCQE = RegisterDial("AxlFFCCQE", cbfa.fnemdls_gen.mdlqeaf);
  DocumentDial(AxlFFCCQE, "Acceptable values: 1. Dipole, 2. BBBA07, 3. 2 Comp, 4. 3 Comp., 5. Z Exp.");
  SCCVecQE = RegisterDial("SCCVecQE", cbfa.fnemdls_gen.sccfv);
  SCCAxlQE = RegisterDial("SCCAxlQE", cbfa.fnemdls_gen.sccfa);
  PsFF = RegisterDial("PsFF", cbfa.fnemdls_gen.fpqe);

  FAxlCCQEAlpha = RegisterDial("FAxlCCQEAlpha", cbfa.fnemdls_gen.axffalpha,0.123,0.123);
  FAxlCCQEGamma = RegisterDial("FAxlCCQEGamma", cbfa.fnemdls_gen.axffgamma,0.121,0.121);
  FAxlCCQEBeta = RegisterDial("FAxlCCQEBeta", cbfa.fnemdls_gen.axffbeta),0.178,0.178;
  FAxlCCQETheta = RegisterDial("FAxlCCQETheta", cbfa.fnemdls_gen.axfftheta,0.031,0.031);

  // Not current used for reweighting, but could be
  // FAZExp_NTerms = RegisterDial("FAZExp_NTerms", cbfa.fnemdls_gen.axzexpnt);
  // FAZExp_TCut = RegisterDial("FAZExp_TCut", cbfa.fnemdls_gen.axzexptc);
  // FAZExp_T0 = RegisterDial("FAZExp_T0", cbfa.fnemdls_gen.axzexpt0);
  // FAZExp_Q4Cut = RegisterDial("FAZExp_Q4Cut", cbfa.fnemdls_gen.axzexpq4);

  FAZExp_A1 = RegisterDial("FAZExp_A1", cbfa.fnemdls_gen.axzexpa1,0.186,0.186);
  FAZExp_A2 = RegisterDial("FAZExp_A2", cbfa.fnemdls_gen.axzexpa2,1.559,1.559);
  FAZExp_A3 = RegisterDial("FAZExp_A3", cbfa.fnemdls_gen.axzexpa3,3.836,3.836);
  FAZExp_A4 = RegisterDial("FAZExp_A4", cbfa.fnemdls_gen.axzexpa4,4.08,4.08);
}

void NuXSecCCQEEngine::Reconfigure() {
  ApplyTweaks();

  if (fDials[AxlFFCCQE].IsTweaked() &&
      ((fDials[AxlFFCCQE].ToValue > 5) || (fDials[AxlFFCCQE].ToValue < 1))) {
    std::cout << "ERROR: AxlFFCCQE can only go between 1 and 5\n\t1. "
                 "Dipole\n\t2. BBBA07\n\t3. 2 Comp\n\t4. 3 Comp.\n\t5. Z Exp."
              << std::endl;
    abort();
  }
}

double NuXSecCCQEEngine::CalcWeight() {
  if (!NModeDefn::isCCQE(nework_.modene)) {
    return 1;
  }
  if (!AnyTweaked()) {
    return 1;
  }

  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();
  cbfa.ResetGenValues();

  double old_xsec = NEUTGetXSec();

  if (old_xsec == 0) {
    std::cout << "NuXSecCCQEEngine::CalcWeight() Warning: old_xsec==0, setting "
                 "weight to 1"
              << std::endl;
    return 1;
  }

  nemdls_.mdlqeaf = fDials[AxlFFCCQE].ToValue;

  nemdls_.xmaqe = fDials[MaCCQE].ToValue;

  nemdls_.sccfv = fDials[SCCVecQE].ToValue;
  nemdls_.sccfa = fDials[SCCAxlQE].ToValue;
  nemdls_.fpqe = fDials[PsFF].ToValue;

  nemdls_.axffalpha = fDials[FAxlCCQEAlpha].ToValue;
  nemdls_.axffgamma = fDials[FAxlCCQEGamma].ToValue;
  nemdls_.axfftheta = fDials[FAxlCCQETheta].ToValue;
  nemdls_.axffbeta = fDials[FAxlCCQEBeta].ToValue;

  nemdls_.axzexpa1 = fDials[FAZExp_A1].ToValue;
  nemdls_.axzexpa2 = fDials[FAZExp_A2].ToValue;
  nemdls_.axzexpa3 = fDials[FAZExp_A3].ToValue;
  nemdls_.axzexpa4 = fDials[FAZExp_A4].ToValue;

  NEUTSetParams();

  double new_xsec = NEUTGetXSec();

#ifdef _N_REWEIGHT_CCQE_DEBUG_
  std::cout << "differential cross section (old) = " << old_xsec << std::endl;
  std::cout << "differential cross section (new) = " << new_xsec << std::endl;
#endif

  return CheckReturnWeight(new_xsec / old_xsec);
}

} // namespace rew
} // namespace neut
