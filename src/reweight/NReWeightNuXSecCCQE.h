//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightNuXSecCCQE

\brief    Reweighting CCQE NEUT neutrino cross sections

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
          University of Toronto

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_NU_XSEC_CCQE_H_
#define _N_REWEIGHT_NU_XSEC_CCQE_H_

#include <iostream>

#include "NReWeightI.h"

namespace neut {
namespace rew {

class NReWeightNuXSecCCQE : public NReWeightI {
public:
  NReWeightNuXSecCCQE();
  ~NReWeightNuXSecCCQE();

  // implement the NReWeightI interface
  bool IsHandled(NSyst_t syst);
  void SetSystematic(NSyst_t syst, double val);
  void Reset();
  void Reconfigure();
  double CalcWeight();
  double CalcChisq();

private:
  void Init();
  double CalcWeightMa();

  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MaCCQE);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_AxlFFCCQE);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_VecFFCCQE);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_SCCVecQE);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_SCCAxlQE);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_PsFF);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_AxlDipToAlt);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAxlCCQEAlpha);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAxlCCQEGamma);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAxlCCQEBeta);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAxlCCQETheta);
  // NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_NTerms);
  // NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_TCut);
  // NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_T0);
  // NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_Q4Cut);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_A1);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_A2);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_A3);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FAZExp_A4);

};

} // namespace rew
} // namespace neut

#endif
