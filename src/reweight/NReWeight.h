//____________________________________________________________________________
/*!

\class    neut::rew::NReWeight

\brief    Interface to the NEUT event reweighting engines

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
          University of Toronto

\created  Apr 5, 2011

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_H_
#define _N_REWEIGHT_H_

#include <map>
#include <string>

#include "NReWeightI.h"
#include "NSystSet.h"

namespace neut {
namespace rew {

class NReWeight {
public:
  NReWeight();
  ~NReWeight();

  /// Add concrete weight calculator, transfers ownership
  void AdoptWghtCalc(std::string name, NReWeightI *wcalc);

  /// access a weight calculator by name
  NReWeightI *WghtCalc(std::string name);
  /// set of enabled systematic params & values
  NSystSet &Systematics(void);
  /// reconfigure weight calculators with new params
  void Reconfigure(void);
  /// calculate weight for input event
  double CalcWeight();
  /// calculate penalty chisq for current values of tweaking dials
  double CalcChisq(void);
  void Print(void);

private:
  void CleanUp(void);

  /// set of enabled nuisance parameters
  NSystSet fSystSet;
  /// concrete weight calculators
  std::map<std::string, NReWeightI *> fWghtCalc;
};

} // namespace rew
} // namespace neut

#endif
