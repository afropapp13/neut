//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightNuXSecRES

\brief    Reweight NEUT CC resonance neutrino-production

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

#ifndef _N_REWEIGHT_NU_XSEC_RES_H_
#define _N_REWEIGHT_NU_XSEC_RES_H_

#include "NReWeightI.h"
#include "NSyst.h"

namespace neut {
namespace rew {

class NReWeightNuXSecRES : public NReWeightI {
public:
  NReWeightNuXSecRES();
  ~NReWeightNuXSecRES();

  // implement the NReWeightI interface
  bool IsHandled(NSyst_t syst);
  void SetSystematic(NSyst_t syst, double val);
  void Reset(void);
  void Reconfigure(void);
  double CalcWeight();
  double CalcChisq(void);
  void SetUseAngular(bool val) { 
    std::cout << "Using angular information in reweight: " << val << std::endl;
    UseAngular = val; 
  };

private:
  void Init(void);
  bool UseAngular;

  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MaRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MvRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_FFRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_TypeRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_CA5RES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_BgSclRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MaNFFRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MvNFFRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MaRSRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MvRSRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_UseSeparateBgSclLMCPiBar);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_BgSclLMCPiBarRES);
  NSYST_DECLAREDIALVARIABLES(kXSecTwkDial_MDLSPiEj);
};

} // namespace rew
} // namespace neut

#endif
