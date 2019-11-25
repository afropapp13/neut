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
 @ Oct 20, 2010 - CA
   Made static consts `kModeMaMv' and `kModeNormAndMaMvShape' public to
   aid external configuration.
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "CommonBlockIFace.h"
#include "NModeDefn.h"
#include "NReWeightNuXSecRES.h"
#include "NSystUncertainty.h"

#include "neworkC.h"

//#define _N_REWEIGHT_RES_DEBUG_

namespace neut {
namespace rew {

NReWeightNuXSecRES::NReWeightNuXSecRES() { this->Init(); }

NReWeightNuXSecRES::~NReWeightNuXSecRES() {}

void NReWeightNuXSecRES::Init(void) {
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  // From nefillmodel.F
  if (cbfa.fneut1pi_gen.neiff != 0) {
    NSYST_SETDEF(kXSecTwkDial_MaRES, cbfa.fneut1pi_gen.xmanffres);
    NSYST_SETDEF(kXSecTwkDial_MvRES, cbfa.fneut1pi_gen.xmvnffres);
  } else {
    NSYST_SETDEF(kXSecTwkDial_MaRES, cbfa.fneut1pi_gen.xmarsres);
    NSYST_SETDEF(kXSecTwkDial_MvRES, cbfa.fneut1pi_gen.xmvrsres);
  }

  NSYST_SETDEF(kXSecTwkDial_CA5RES, cbfa.fneut1pi_gen.rneca5i);
  NSYST_SETDEF(kXSecTwkDial_BgSclRES, cbfa.fneut1pi_gen.rnebgscl);
  NSYST_SETDEF(kXSecTwkDial_UseSeparateBgSclLMCPiBar, 0);
  NSYST_SETDEF(kXSecTwkDial_BgSclLMCPiBarRES, cbfa.fneut1pi_gen.rnebgscl);
}

bool NReWeightNuXSecRES::IsHandled(NSyst_t syst) {

  NSYST_USESDIAL(syst, kXSecTwkDial_MaRES);
  NSYST_USESDIAL(syst, kXSecTwkDial_MvRES);
  NSYST_USESDIAL(syst, kXSecTwkDial_CA5RES);
  NSYST_USESDIAL(syst, kXSecTwkDial_BgSclRES);
  NSYST_USESDIAL(syst, kXSecTwkDial_UseSeparateBgSclLMCPiBar);
  NSYST_USESDIAL(syst, kXSecTwkDial_BgSclLMCPiBarRES);

  return false;
}

void NReWeightNuXSecRES::SetSystematic(NSyst_t syst, double twk_dial) {

  NSYST_UPDATEVALUE(kXSecTwkDial_MaRES, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_MvRES, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_CA5RES, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_BgSclRES, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_UseSeparateBgSclLMCPiBar, syst, twk_dial);
  NSYST_UPDATEVALUE(kXSecTwkDial_BgSclLMCPiBarRES, syst, twk_dial);
}

void NReWeightNuXSecRES::Reset(void) {
  NSYST_SETTODEF(kXSecTwkDial_MaRES);
  NSYST_SETTODEF(kXSecTwkDial_MvRES);
  NSYST_SETTODEF(kXSecTwkDial_CA5RES);
  NSYST_SETTODEF(kXSecTwkDial_BgSclRES);
  NSYST_SETTODEF(kXSecTwkDial_UseSeparateBgSclLMCPiBar);
  NSYST_SETTODEF(kXSecTwkDial_BgSclLMCPiBarRES);
}

void NReWeightNuXSecRES::Reconfigure(void) {
  NSystUncertainty *fracerr = NSystUncertainty::Instance();

  NSYST_RECONFCURRVALUE(kXSecTwkDial_MaRES, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_MvRES, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_CA5RES, fracerr);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_BgSclRES, fracerr);
  NSYST_RECONFCURRVALUE_NOUNCERT(kXSecTwkDial_UseSeparateBgSclLMCPiBar);
  NSYST_RECONFCURRVALUE(kXSecTwkDial_BgSclLMCPiBarRES, fracerr);
}

double NReWeightNuXSecRES::CalcChisq() {
  double chisq = 0.;
  chisq += std::pow(NSYST_TWKVAR(kXSecTwkDial_MaRES), 2);
  chisq += std::pow(NSYST_TWKVAR(kXSecTwkDial_MvRES), 2);
  return chisq;
}

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

double NReWeightNuXSecRES::CalcWeight() {
  if (!NModeDefn::isRES(nework_.modene)) {
    return 1;
  }

  bool tweaked = NSYST_ISTWKD(kXSecTwkDial_MaRES) ||
                 NSYST_ISTWKD(kXSecTwkDial_MvRES) ||
                 NSYST_ISTWKD(kXSecTwkDial_CA5RES) ||
                 NSYST_ISTWKD(kXSecTwkDial_BgSclRES) ||
                 NSYST_ISTWKD(kXSecTwkDial_UseSeparateBgSclLMCPiBar) ||
                 NSYST_ISTWKD(kXSecTwkDial_BgSclLMCPiBarRES);

  if (!tweaked) {
    return 1.0;
  }

  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();
  cbfa.ResetGenValues();

  double old_xsec = NEUTGetXSec();

  if (old_xsec == 0) {
    std::cout
        << "NReWeightNuXSecRes::CalcWeight() Warning: old_xsec==0, setting "
           "weight to 1"
        << std::endl;
    return 1;
  }

  if (cbfa.fneut1pi_gen.neiff != 0) {
    neut1pi_.xmanffres = NSYST_CURRVAR(kXSecTwkDial_MaRES);
    neut1pi_.xmvnffres = NSYST_CURRVAR(kXSecTwkDial_MvRES);
  } else {
    neut1pi_.xmarsres = NSYST_CURRVAR(kXSecTwkDial_MaRES);
    neut1pi_.xmvrsres = NSYST_CURRVAR(kXSecTwkDial_MvRES);
  }

  neut1pi_.rneca5i = NSYST_CURRVAR(kXSecTwkDial_CA5RES);
  if (NSYST_CURRVAR(kXSecTwkDial_UseSeparateBgSclLMCPiBar) &&
      EventIsLMCPiBar()) {
    neut1pi_.rnebgscl = NSYST_CURRVAR(kXSecTwkDial_BgSclLMCPiBarRES);
  } else {
    neut1pi_.rnebgscl = NSYST_CURRVAR(kXSecTwkDial_BgSclRES);
  }

  NEUTSetParams();

  double new_xsec = NEUTGetXSec();

#ifdef _N_REWEIGHT_RES_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << endl;
  cout << "differential cross section (new) = " << new_xsec << endl;
#endif

  NREWCHECKRETURN(new_xsec / old_xsec);
}

} // namespace rew
} // namespace neut
