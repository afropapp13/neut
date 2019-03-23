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
 @ May 25, 2010 - CA
   First included in v2.7.1.
 @ Apr 13, 2011 - PD
   Implemented for NEUT
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NReWeightNuXSecNC.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

//_______________________________________________________________________________________
NReWeightNuXSecNC::NReWeightNuXSecNC() 
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuXSecNC::~NReWeightNuXSecNC()
{

}

//_______________________________________________________________________________________
void NReWeightNuXSecNC::Init(void)
{  
  fortFns = NFortFns::Instance();

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;

  this->RewNue       (true);
  this->RewNuebar    (true);
  this->RewNumu      (true);
  this->RewNumubar   (true);

  this->RewQE  (true );
  this->RewRES (false); // assume NReWeightNuXSecNCRES is going to be your 1st choice
  this->RewDIS (false); // assume NReWeightNuXSecDIS   is going to be your 1st choice

}
//_______________________________________________________________________________________
bool NReWeightNuXSecNC::IsHandled(NSyst_t syst)
{
   switch(syst) {
     case ( kXSecTwkDial_NC ) : 
       return true;
       break;
     default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
void NReWeightNuXSecNC::SetSystematic(NSyst_t syst, double twk_dial)
{
   switch(syst) {
     case ( kXSecTwkDial_NC ) : 
       fNormTwkDial = twk_dial;
       break;
     default:
       return;
       break;
   }
}
//_______________________________________________________________________________________
void NReWeightNuXSecNC::Reset(void)
{
  fNormTwkDial = 0.;
}
//_______________________________________________________________________________________
void NReWeightNuXSecNC::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  int    sign_normtwk = utils::rew::Sign(fNormTwkDial);
  double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NC,    sign_normtwk);
  fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
  
  fNormCurr = TMath::Max(0., fNormCurr);
}
//_______________________________________________________________________________________
double NReWeightNuXSecNC::CalcWeight() 
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  bool is_nc  = modeDefn.isNC(nework_.modene);
  if(!is_nc) return 1.;

  bool is_ncqe = modeDefn.isNCEL(nework_.modene);
  if(is_ncqe && !fRewQE) return 1.;

  bool is_ncres = modeDefn.isNCRES(nework_.modene);
  if(is_ncres && !fRewRES) return 1.;

  bool is_dis = ( modeDefn.isDIS(nework_.modene) || modeDefn.isMPI(nework_.modene) );
  if(is_dis && !fRewDIS) return 1.;

  int nupdg = nework_.ipne[0];
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  return fNormCurr;
}
//_______________________________________________________________________________________
double NReWeightNuXSecNC::CalcChisq(void)
{
  double chisq = TMath::Power(fNormTwkDial,  2.);
  return chisq;
}
//_______________________________________________________________________________________


