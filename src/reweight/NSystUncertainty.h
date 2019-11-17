//____________________________________________________________________________
/*!

\class    neut::rew::NSystUncertainty

\brief

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
          University of Toronto

\created  Sep 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_SYST_UNCERTAINTY_H_
#define _N_SYST_UNCERTAINTY_H_

#include "NSyst.h"

#include <map>

namespace neut {
namespace rew {

class NSystUncertainty {

public:
  static NSystUncertainty *Instance(void);

  double OneSigmaErr(NSyst_t syst, int sign = 0) const;
  void SetUncertainty(NSyst_t syst, double plus_err, double minus_err);
  void GetUncertainty(NSyst_t syst, double &plus_err, double &minus_err);

private:
  void SetDefaults(void);

  std::map<NSyst_t, double> fOneSigPlusErrMap; // + err
  std::map<NSyst_t, double> fOneSigMnusErrMap; // - err

  NSystUncertainty();
  NSystUncertainty(const NSystUncertainty &err);
  ~NSystUncertainty();

  static NSystUncertainty *fInstance;
};

} // namespace rew
} // namespace neut

#endif
