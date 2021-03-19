//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
          University of Toronto

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Apr 27, 2010 - CA
   Included new parameters in preparation for the Summer 2010 T2K analyses.
   Added option to override the default 1\sigma errors.
 @ Nov 25, 2010 - CA
   Allow for asymmetric 1 sigma fractional errors.
 @ Apr. 6, 2011 - PD
   Implement for NEUT
*/
//____________________________________________________________________________

// The order here allows us to leave the xmacro defined, do not move the NSyst
// include after the NSystUncertainty.h include

#define NEUTREWEIGHT_LEAVE_DIALS_DEFINED
#include "NSyst.h"

#include "NSystUncertainty.h"

#include <iostream>

namespace neut {
namespace rew {

NSystUncertainty *NSystUncertainty::fInstance = 0;
NSystUncertainty::NSystUncertainty() {}
NSystUncertainty::~NSystUncertainty() { fInstance = 0; }
NSystUncertainty *NSystUncertainty::Instance() {
  if (fInstance == 0) {
    fInstance = new NSystUncertainty();
    fInstance->SetDefaults();
  }
  return fInstance;
}

double NSystUncertainty::OneSigmaErr(NSyst_t s, int sign) const {
  if (sign > 0) {
    std::map<NSyst_t, double>::const_iterator it = fOneSigPlusErrMap.find(s);
    if (it != fOneSigPlusErrMap.end()) {
      return it->second;
    }
    return 0;
  } else if (sign < 0) {
    std::map<NSyst_t, double>::const_iterator it = fOneSigMnusErrMap.find(s);
    if (it != fOneSigMnusErrMap.end()) {
      return it->second;
    }
    return 0;
  } else {
    // Handle default argument (sign=0)
    // Case added for compatibility purposes since most existing weight
    // calcutators call NSystUncertainty::OneSigmaErr(NSyst_t) and the error
    // on most NSyst_t params is symmetric.
    double err = 0.5 * (this->OneSigmaErr(s, +1) + this->OneSigmaErr(s, -1));
    return err;
  }
}

void NSystUncertainty::SetUncertainty(NSyst_t s, double plus_err,
                                      double minus_err) {
  fOneSigPlusErrMap[s] = plus_err;
  fOneSigMnusErrMap[s] = minus_err;
}

void NSystUncertainty::GetUncertainty(NSyst_t s, double &plus_err,
                                      double &minus_err) {
  plus_err = fOneSigPlusErrMap[s];
  minus_err = fOneSigMnusErrMap[s];
}

void NSystUncertainty::SetDefaults(void) {

#define X(A, B, C, D) this->SetUncertainty(A, C, D);
  NEUTREWEIGHT_DIAL_LIST
#undef X
}

} // namespace rew
} // namespace neut
