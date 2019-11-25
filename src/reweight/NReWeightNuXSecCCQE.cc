//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
          University of Toronto

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ May 17, 2010 - CA
   Code extracted from NReWeightNuXSec and redeveloped in preparation for
   the Summer 2010 T2K analyses.
 @ Oct 22, 2010 - CA
   Make static consts kModeMa and kModeNormAndMaShape public to aid
   external configuration.
 @ Nov 25, 2010 - CA
   Allow asymmetric errors
 @ Apr 6, 2011 - PD
   Implement for NEUT
 @ Nov 2019 - LP
   Refactor NEUT reweight
*/
//____________________________________________________________________________

#include "NReWeightNuXSecCCQE.h"
#include "CommonBlockIFace.h"
#include "NModeDefn.h"
#include "NSystUncertainty.h"

#include "neworkC.h"

#include <iostream>

//#define _N_REWEIGHT_CCQE_DEBUG_

namespace neut {
namespace rew {

NReWeightNuXSecCCQE::NReWeightNuXSecCCQE() { this->Init(); }

NReWeightNuXSecCCQE::~NReWeightNuXSecCCQE() {}

void NReWeightNuXSecCCQE::Init() {

  // Get the parameter values at generation time and store them as the 'default'
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  NSYST_SETDEF(kXSecTwkDial_MaCCQE, cbfa.fnemdls_gen.xmaqe);
  NSYST_SETDEF(kXSecTwkDial_AxlFFCCQE, cbfa.fnemdls_gen.mdlqeaf);
  NSYST_SETDEF(kXSecTwkDial_SCCVecQE, cbfa.fnemdls_gen.sccfv);
  NSYST_SETDEF(kXSecTwkDial_SCCAxlQE, cbfa.fnemdls_gen.sccfa);
  NSYST_SETDEF(kXSecTwkDial_PsFF, cbfa.fnemdls_gen.fpqe);
  NSYST_SETDEF(kXSecTwkDial_FAxlCCQEAlpha, cbfa.fnemdls_gen.axffalpha);
  NSYST_SETDEF(kXSecTwkDial_FAxlCCQEGamma, cbfa.fnemdls_gen.axffgamma);
  NSYST_SETDEF(kXSecTwkDial_FAxlCCQEBeta, cbfa.fnemdls_gen.axffbeta);
  NSYST_SETDEF(kXSecTwkDial_FAxlCCQETheta, cbfa.fnemdls_gen.axfftheta);
  // NSYST_SETDEF(kXSecTwkDial_FAZExp_NTerms, cbfa.fnemdls_gen.axzexpnt);
  // NSYST_SETDEF(kXSecTwkDial_FAZExp_TCut, cbfa.fnemdls_gen.axzexptc);
  // NSYST_SETDEF(kXSecTwkDial_FAZExp_T0, cbfa.fnemdls_gen.axzexpt0);
  // NSYST_SETDEF(kXSecTwkDial_FAZExp_Q4Cut, cbfa.fnemdls_gen.axzexpq4);
  NSYST_SETDEF(kXSecTwkDial_FAZExp_A1, cbfa.fnemdls_gen.axzexpa1);
  NSYST_SETDEF(kXSecTwkDial_FAZExp_A2, cbfa.fnemdls_gen.axzexpa2);
  NSYST_SETDEF(kXSecTwkDial_FAZExp_A3, cbfa.fnemdls_gen.axzexpa3);
  NSYST_SETDEF(kXSecTwkDial_FAZExp_A4, cbfa.fnemdls_gen.axzexpa4);
}

bool NReWeightNuXSecCCQE::IsHandled(NSyst_t syst) {
  NSYST_USESDIAL(syst, kXSecTwkDial_MaCCQE);
  NSYST_USESDIAL(syst, kXSecTwkDial_AxlFFCCQE);
  NSYST_USESDIAL(syst, kXSecTwkDial_SCCVecQE);
  NSYST_USESDIAL(syst, kXSecTwkDial_SCCAxlQE);
  NSYST_USESDIAL(syst, kXSecTwkDial_PsFF);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAxlCCQEAlpha);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAxlCCQEGamma);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAxlCCQEBeta);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAxlCCQETheta);
  // Cannot vary these by reweighting
  // NSYST_USESDIAL(syst,kXSecTwkDial_FAZExp_NTerms);
  // NSYST_USESDIAL(syst,kXSecTwkDial_FAZExp_TCut);
  // NSYST_USESDIAL(syst,kXSecTwkDial_FAZExp_T0);
  // NSYST_USESDIAL(syst,kXSecTwkDial_FAZExp_Q4Cut);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAZExp_A1);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAZExp_A2);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAZExp_A3);
  NSYST_USESDIAL(syst, kXSecTwkDial_FAZExp_A4);

  return false;
}

void NReWeightNuXSecCCQE::SetSystematic(NSyst_t syst, double twk_dial) {
  NSYST_UPDATEVALUE(kXSecTwkDial_MaCCQE, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_AxlFFCCQE, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_SCCVecQE, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_SCCAxlQE, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_PsFF, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAxlCCQEAlpha, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAxlCCQEGamma, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAxlCCQEBeta, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAxlCCQETheta, syst, twk_dial);
  // NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_NTerms, syst, twk_dial);
  // NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_TCut, syst, twk_dial);
  // NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_T0, syst, twk_dial);
  // NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_Q4Cut, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_A1, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_A2, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_A3, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_FAZExp_A4, syst, twk_dial);
}

void NReWeightNuXSecCCQE::Reset() {
  NSYST_SETTODEF(kXSecTwkDial_MaCCQE);
  NSYST_SETTODEF(kXSecTwkDial_AxlFFCCQE);
  NSYST_SETTODEF(kXSecTwkDial_SCCVecQE);
  NSYST_SETTODEF(kXSecTwkDial_SCCAxlQE);
  NSYST_SETTODEF(kXSecTwkDial_PsFF);
  NSYST_SETTODEF(kXSecTwkDial_FAxlCCQEAlpha);
  NSYST_SETTODEF(kXSecTwkDial_FAxlCCQEGamma);
  NSYST_SETTODEF(kXSecTwkDial_FAxlCCQEBeta);
  NSYST_SETTODEF(kXSecTwkDial_FAxlCCQETheta);
  // NSYST_SETTODEF(kXSecTwkDial_FAZExp_NTerms);
  // NSYST_SETTODEF(kXSecTwkDial_FAZExp_TCut);
  // NSYST_SETTODEF(kXSecTwkDial_FAZExp_T0);
  // NSYST_SETTODEF(kXSecTwkDial_FAZExp_Q4Cut);
  NSYST_SETTODEF(kXSecTwkDial_FAZExp_A1);
  NSYST_SETTODEF(kXSecTwkDial_FAZExp_A2);
  NSYST_SETTODEF(kXSecTwkDial_FAZExp_A3);
  NSYST_SETTODEF(kXSecTwkDial_FAZExp_A4);

  Reconfigure();
}

void NReWeightNuXSecCCQE::Reconfigure() {
  NSystUncertainty *fracerr = NSystUncertainty::Instance();

  NSYST_RECONFCURRVALUE(kXSecTwkDial_MaCCQE, fracerr);
  NSYST_RECONFCURRVALUE_NOUNCERT(kXSecTwkDial_AxlFFCCQE);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_SCCVecQE, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_SCCAxlQE, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_PsFF, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAxlCCQEAlpha, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAxlCCQEGamma, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAxlCCQEBeta, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAxlCCQETheta, fracerr);
  // NSYST_RECONFCURRVALUE_NOUNCERT(kXSecTwkDial_FAZExp_NTerms);
  // NSYST_RECONFCURRVALUE(kXSecTwkDial_FAZExp_TCut, fracerr);
  // NSYST_RECONFCURRVALUE(kXSecTwkDial_FAZExp_T0, fracerr);
  // NSYST_RECONFCURRVALUE(kXSecTwkDial_FAZExp_Q4Cut, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAZExp_A1, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAZExp_A2, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAZExp_A3, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_FAZExp_A4, fracerr);

  if (NSYST_ISTWKD(kXSecTwkDial_AxlFFCCQE) &&
      ((NSYST_CURRVAR(kXSecTwkDial_AxlFFCCQE) > 5) ||
       (NSYST_CURRVAR(kXSecTwkDial_AxlFFCCQE) < 1))) {
    std::cout << "ERROR: AxlFFCCQE can only go between 1 and 5\n\t1. "
                 "Dipole\n\t2. BBBA07\n\t3. 2 Comp\n\t4. 3 Comp.\n\t5. Z Exp."
              << std::endl;
    throw;
  }
}

double NReWeightNuXSecCCQE::CalcChisq() {
  return std::pow(NSYST_TWKVAR(kXSecTwkDial_MaCCQE), 2);
}

double NReWeightNuXSecCCQE::CalcWeight() {
  if (!NModeDefn::isCCQE(nework_.modene)) {
    return 1;
  }
  bool tweaked = NSYST_ISTWKD(kXSecTwkDial_MaCCQE) ||
                 NSYST_ISTWKD(kXSecTwkDial_AxlFFCCQE) ||
                 NSYST_ISTWKD(kXSecTwkDial_SCCVecQE) ||
                 NSYST_ISTWKD(kXSecTwkDial_SCCAxlQE) ||
                 NSYST_ISTWKD(kXSecTwkDial_PsFF) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAxlCCQEAlpha) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAxlCCQEGamma) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAxlCCQEBeta) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAxlCCQETheta) ||
                 //  NSYST_ISTWKD(kXSecTwkDial_FAZExp_NTerms) ||
                 //  NSYST_ISTWKD(kXSecTwkDial_FAZExp_TCut) ||
                 //  NSYST_ISTWKD(kXSecTwkDial_FAZExp_T0) ||
                 //  NSYST_ISTWKD(kXSecTwkDial_FAZExp_Q4Cut) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAZExp_A1) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAZExp_A2) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAZExp_A3) ||
                 NSYST_ISTWKD(kXSecTwkDial_FAZExp_A4);

  if (!tweaked) {
    return 1;
  }

  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();
  cbfa.ResetGenValues();

  double old_xsec = NEUTGetXSec();

  if (old_xsec == 0) {
    std::cout
        << "NReWeightNuXSecCCQE::CalcWeight() Warning: old_xsec==0, setting "
           "weight to 1"
        << std::endl;
    return 1;
  }

  NSYST_ASSIGNIFTWKD(nemdls_.mdlqeaf, kXSecTwkDial_AxlFFCCQE);

  NSYST_ASSIGNIFTWKD(nemdls_.xmaqe, kXSecTwkDial_MaCCQE);

  NSYST_ASSIGNIFTWKD(nemdls_.sccfv, kXSecTwkDial_SCCVecQE);
  NSYST_ASSIGNIFTWKD(nemdls_.sccfa, kXSecTwkDial_SCCAxlQE);
  NSYST_ASSIGNIFTWKD(nemdls_.fpqe, kXSecTwkDial_PsFF);

  NSYST_ASSIGNIFTWKD(nemdls_.axffalpha, kXSecTwkDial_FAxlCCQEAlpha);
  NSYST_ASSIGNIFTWKD(nemdls_.axffgamma, kXSecTwkDial_FAxlCCQEGamma);

  NSYST_ASSIGNIFTWKD(nemdls_.axffalpha, kXSecTwkDial_FAxlCCQEAlpha);
  NSYST_ASSIGNIFTWKD(nemdls_.axffgamma, kXSecTwkDial_FAxlCCQEGamma);
  NSYST_ASSIGNIFTWKD(nemdls_.axfftheta, kXSecTwkDial_FAxlCCQETheta);
  NSYST_ASSIGNIFTWKD(nemdls_.axffbeta, kXSecTwkDial_FAxlCCQEBeta);

  NSYST_ASSIGNIFTWKD(nemdls_.axzexpa1, kXSecTwkDial_FAZExp_A1);
  NSYST_ASSIGNIFTWKD(nemdls_.axzexpa2, kXSecTwkDial_FAZExp_A2);
  NSYST_ASSIGNIFTWKD(nemdls_.axzexpa3, kXSecTwkDial_FAZExp_A3);
  NSYST_ASSIGNIFTWKD(nemdls_.axzexpa4, kXSecTwkDial_FAZExp_A4);

  NEUTSetParams();

  double new_xsec = NEUTGetXSec();

#ifdef _N_REWEIGHT_CCQE_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << endl;
  cout << "differential cross section (new) = " << new_xsec << endl;
#endif
  NREWCHECKRETURN(new_xsec / old_xsec);
}

} // namespace rew
} // namespace neut
