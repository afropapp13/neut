//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightI

\brief    NEUT event reweighting engine ABC

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_ABC_H_
#define _N_REWEIGHT_ABC_H_

#include "NSyst.h"

namespace neut {

namespace rew {

class NReWeightI {
public:
  virtual ~NReWeightI() {}

  /// does the current weight calculator handle the input nuisance param?
  virtual bool IsHandled(NSyst_t syst) = 0;

  /// update the value for the specified nuisance param
  virtual void SetSystematic(NSyst_t syst, double val) = 0;

  ///  set all nuisance parameters to default values
  virtual void Reset(void) = 0;

  /// propagate updated nuisance parameter values to actual MC, etc
  virtual void Reconfigure(void) = 0;

  /// calculate a weight for the input event using the current nuisance param
  /// values
  virtual double CalcWeight() = 0;

  /// calculate penalty factors
  virtual double CalcChisq(void) = 0;

  // Get the variation in sigma units for a given absolute dial value
  virtual double GetTwkForAbs(NSyst_t syst, double val) { return 0; }

  virtual std::string ToString() {
    return "No ToString implemented for this NReWeightI";
  };
};

} // namespace rew
} // namespace neut

#define NREWCHECKRETURN(a)                                                     \
  if ((a) && !std::isnormal(a)) {                                              \
    std::cout << "[WARN]: Abnormal weight being returned from " << __FILE__    \
              << ":" << __LINE__ << std::endl;                                 \
  }                                                                            \
  return (a)

#endif
