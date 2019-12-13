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
 @ Sep 10, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Apr 6, 2011 - PD
   Implemented for NEUT

*/
//____________________________________________________________________________

#include "NReWeightCasc.h"
#include "CommonBlockIFace.h"
#include "NSystUncertainty.h"

#include "vcworkC.h"

#include "posinnucC.h"

#include <cmath>
#include <iostream>

//#define _N_REWEIGHT_CASC_DEBUG_

namespace neut {
namespace rew {

NReWeightCasc::NReWeightCasc() { this->Init(); }
NReWeightCasc::~NReWeightCasc() {}

void NReWeightCasc::Init(void) {
  NSYST_SETDEF(kCascTwkDial_FrAbs_pi, neffpr_.fefabs);
  NSYST_SETDEF(kCascTwkDial_FrInelLow_pi, neffpr_.fefqe);
  NSYST_SETDEF(kCascTwkDial_FrCExLow_pi, neffpr_.fefcx);
  NSYST_SETDEF(kCascTwkDial_FrInelHigh_pi, neffpr_.fefinel);
  NSYST_SETDEF(kCascTwkDial_FrCExHigh_pi, neffpr_.fefcxh);
  NSYST_SETDEF(kCascTwkDial_FrPiProd_pi, neffpr_.fefqeh);
}

bool NReWeightCasc::IsHandled(NSyst_t syst) {
  NSYST_USESDIAL(syst, kCascTwkDial_FrAbs_pi);
  NSYST_USESDIAL(syst, kCascTwkDial_FrInelLow_pi);
  NSYST_USESDIAL(syst, kCascTwkDial_FrCExLow_pi);
  NSYST_USESDIAL(syst, kCascTwkDial_FrInelHigh_pi);
  NSYST_USESDIAL(syst, kCascTwkDial_FrCExHigh_pi);
  NSYST_USESDIAL(syst, kCascTwkDial_FrPiProd_pi);

  return false;
}

void NReWeightCasc::SetSystematic(NSyst_t syst, double twk_dial) {
  NSYST_UPDATEVALUE(kCascTwkDial_FrAbs_pi, syst, twk_dial);
  NSYST_UPDATEVALUE(kCascTwkDial_FrInelLow_pi, syst, twk_dial);
  NSYST_UPDATEVALUE(kCascTwkDial_FrCExLow_pi, syst, twk_dial);
  NSYST_UPDATEVALUE(kCascTwkDial_FrInelHigh_pi, syst, twk_dial);
  NSYST_UPDATEVALUE(kCascTwkDial_FrCExHigh_pi, syst, twk_dial);
  NSYST_UPDATEVALUE(kCascTwkDial_FrPiProd_pi, syst, twk_dial);
}

void NReWeightCasc::Reset(void) {
  NSYST_SETTODEF(kCascTwkDial_FrAbs_pi);
  NSYST_SETTODEF(kCascTwkDial_FrInelLow_pi);
  NSYST_SETTODEF(kCascTwkDial_FrCExLow_pi);
  NSYST_SETTODEF(kCascTwkDial_FrInelHigh_pi);
  NSYST_SETTODEF(kCascTwkDial_FrCExHigh_pi);
  NSYST_SETTODEF(kCascTwkDial_FrPiProd_pi);

  Reconfigure();
}

void NReWeightCasc::Reconfigure(void) {

  // For 1-sigma change
  NSystUncertainty *fracerr = NSystUncertainty::Instance();

  NSYST_RECONFCURRVALUE(kCascTwkDial_FrAbs_pi, fracerr);
  NSYST_RECONFCURRVALUE(kCascTwkDial_FrInelLow_pi, fracerr);
  NSYST_RECONFCURRVALUE(kCascTwkDial_FrCExLow_pi, fracerr);
  NSYST_RECONFCURRVALUE(kCascTwkDial_FrInelHigh_pi, fracerr);
  NSYST_RECONFCURRVALUE(kCascTwkDial_FrCExHigh_pi, fracerr);
  NSYST_RECONFCURRVALUE(kCascTwkDial_FrPiProd_pi, fracerr);
}

double NReWeightCasc::CalcWeight() {

  bool tweaked = NSYST_ISTWKD(kCascTwkDial_FrAbs_pi) ||
                 NSYST_ISTWKD(kCascTwkDial_FrInelLow_pi) ||
                 NSYST_ISTWKD(kCascTwkDial_FrCExLow_pi) ||
                 NSYST_ISTWKD(kCascTwkDial_FrInelHigh_pi) ||
                 NSYST_ISTWKD(kCascTwkDial_FrCExHigh_pi) ||
                 NSYST_ISTWKD(kCascTwkDial_FrPiProd_pi);

  if (!tweaked) {
    return 1.0;
  }

  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();
  cbfa.ResetGenValues();

  // Not bound
  if (!posinnuc_.ibound) {
    return 1;
  }

  double old_xsec = NEUTGetPiCascProb();

  if (old_xsec <= 0) {
    std::cout << "[WARN]: NReWeightCasc() evpiprob old_xsec <= 0, returning "
                 "weight = 1"
              << std::endl;
    return 1;
  }
  // else if (old_xsec != 1) {
  //   std::cout << "[INFO]: evpiprob old_xsec non zero! " << old_xsec
  //             << std::endl;
  // }

  neffpr_.fefabs = NSYST_CURRVAR(kCascTwkDial_FrAbs_pi);

  neffpr_.fefqe = NSYST_CURRVAR(kCascTwkDial_FrInelLow_pi);
  neffpr_.fefqeh = NSYST_CURRVAR(kCascTwkDial_FrInelHigh_pi);

  neffpr_.fefcx = NSYST_CURRVAR(kCascTwkDial_FrCExLow_pi);
  neffpr_.fefcxh = NSYST_CURRVAR(kCascTwkDial_FrCExHigh_pi);

  neffpr_.fefinel = NSYST_CURRVAR(kCascTwkDial_FrPiProd_pi);

  NEUTSetParams();

  double new_xsec = NEUTGetPiCascProb();

  // if (old_xsec != 1) {
  //   std::cout << "[INFO]: evpiprob new_xsec = " << new_xsec << ", weight = "
  //   << (new_xsec/old_xsec) << std::endl;
  // }
#ifdef _N_REWEIGHT_CASC_DEBUG_
  cout << "pion cascade probability (old) = " << old_xsec << endl;
  cout << "pion cascade probability (new) = " << new_xsec << endl;
#endif
  NREWCHECKRETURN(new_xsec / old_xsec);
}

double NReWeightCasc::CalcChisq(void) {
  double chisq = 0.;

  // For 1-sigma change
  chisq += std::pow(NSYST_TWKVAR(kCascTwkDial_FrAbs_pi), 2.);
  chisq += std::pow(NSYST_TWKVAR(kCascTwkDial_FrInelLow_pi), 2.);
  chisq += std::pow(NSYST_TWKVAR(kCascTwkDial_FrCExLow_pi), 2.);
  chisq += std::pow(NSYST_TWKVAR(kCascTwkDial_FrInelHigh_pi), 2.);
  chisq += std::pow(NSYST_TWKVAR(kCascTwkDial_FrCExHigh_pi), 2.);
  chisq += std::pow(NSYST_TWKVAR(kCascTwkDial_FrPiProd_pi), 2.);

  return chisq;
}

std::vector<double> NReWeightCasc::GetCurrParVals(void) {
  std::vector<double> parVals;
  parVals.push_back(NSYST_TWKVAR(kCascTwkDial_FrAbs_pi));
  parVals.push_back(NSYST_TWKVAR(kCascTwkDial_FrInelLow_pi));
  parVals.push_back(NSYST_TWKVAR(kCascTwkDial_FrCExLow_pi));
  parVals.push_back(NSYST_TWKVAR(kCascTwkDial_FrInelHigh_pi));
  parVals.push_back(NSYST_TWKVAR(kCascTwkDial_FrCExHigh_pi));
  parVals.push_back(NSYST_TWKVAR(kCascTwkDial_FrPiProd_pi));

  return parVals;
}

} // namespace rew
} // namespace neut
