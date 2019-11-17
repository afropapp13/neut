//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightCasc

\brief    Reweighting NEUT FSI Cascade Model

          See T2K-TN-033

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
          University of Toronto

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_CASC_H_
#define _N_REWEIGHT_CASC_H_

#include <vector>

#include "NReWeightI.h"

namespace neut {
namespace rew {

class NReWeightCasc : public NReWeightI {
public:
  NReWeightCasc();
  ~NReWeightCasc();

  // implement the NReWeightI interface
  bool IsHandled(NSyst_t syst);
  void SetSystematic(NSyst_t syst, double val);
  void Reset(void);
  void Reconfigure(void);
  double CalcWeight();
  double CalcChisq(void);

  std::vector<double> GetCurrParVals(void);

private:
  void Init(void);

  NSYST_DECLAREDIALVARIABLES(kCascTwkDial_FrAbs_pi);
  NSYST_DECLAREDIALVARIABLES(kCascTwkDial_FrInelLow_pi);
  NSYST_DECLAREDIALVARIABLES(kCascTwkDial_FrCExLow_pi);
  NSYST_DECLAREDIALVARIABLES(kCascTwkDial_FrInelHigh_pi);
  NSYST_DECLAREDIALVARIABLES(kCascTwkDial_FrCExHigh_pi);
  NSYST_DECLAREDIALVARIABLES(kCascTwkDial_FrPiProd_pi);
};

} // namespace rew
} // namespace neut

#endif
