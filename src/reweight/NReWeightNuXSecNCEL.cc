//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 25, 2010 - CA
   First included in v2.7.1. Added option to reweight NCEL axial mass M_{A} 
   and the strange axial form factor eta.
 @ Apr 13, 2011 - PD
   Implemented for NEUT
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NReWeightNuXSecNCEL.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_NCEL_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

//_______________________________________________________________________________________
NReWeightNuXSecNCEL::NReWeightNuXSecNCEL() 
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuXSecNCEL::~NReWeightNuXSecNCEL()
{
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Init(void)
{
  fortFns = NFortFns::Instance();

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fMaTwkDial   = 0.; 
  fMaDef       = fortFns->XMAQEdef;
  fMaCurr      = fMaDef;
  //fEtaTwkDial  = 0.; 
  //fEtaDef      = fortFns->XETAdef;
  //fEtaCurr     = fEtaDef;

}
//_______________________________________________________________________________________
bool NReWeightNuXSecNCEL::IsHandled(NSyst_t syst)
{
   bool handle;
   switch(syst) {
     case ( kXSecTwkDial_MaNCEL  ) : handle = true; break;
       //case ( kXSecTwkDial_EtaNCEL ) : handle = true; break;
     default:
          handle = false;
          break;
   }
   return handle;
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::SetSystematic(NSyst_t syst, double twk_dial)
{
  switch(syst) {
    case ( kXSecTwkDial_MaNCEL ) :
      fMaTwkDial = twk_dial;
      break;
    //case ( kXSecTwkDial_EtaNCEL ) :
    //  fEtaTwkDial = twk_dial;
    //  break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Reset(void)
{
  fMaTwkDial   = 0.; 
  fMaCurr      = fortFns->XMAQEdef;
  //fEtaTwkDial  = 0.; 
  //fEtaCurr     = fEtaDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  int sign_matwk  = utils::rew::Sign(fMaTwkDial );
  //int sign_etatwk = utils::rew::Sign(fEtaTwkDial);
  
  double fracerr_ma  = fracerr->OneSigmaErr(kXSecTwkDial_MaNCEL,  sign_matwk );
  //double fracerr_eta = fracerr->OneSigmaErr(kXSecTwkDial_EtaNCEL, sign_etatwk);
  
  fMaCurr  = fMaDef  * (1. + fMaTwkDial  * fracerr_ma);
  //fEtaCurr = fEtaDef * (1. + fEtaTwkDial * fracerr_eta);

  fMaCurr  = TMath::Max(0., fMaCurr  );
  //fEtaCurr = TMath::Max(0., fEtaCurr );

}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeight() 
{

  bool is_ncqe = modeDefn.isNCEL(nework_.modene);
  if(!is_ncqe) return 1.;

  int nupdg = nework_.ipne[0];
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  bool tweaked  = (TMath::Abs(fMaTwkDial ) > controls::kASmallNum);
  // || (TMath::Abs(fEtaTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;


  fortFns->SetDefaults();
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif
  
  float old_xsec     = fortFns->evdifcrs();
  //float old_xsec   = event.DiffXSec();
  //float old_weight = event.Weight();

  if (old_xsec==0) {
    cout << "NReWeightNuXSecNCEL::CalcWeight() Warning: old_xsec==0, setting weight to 1" << endl;
    return 1;
  }

  nemdls_.xmaqe = fMaCurr;
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  fortFns->print_allparams();
#endif
  float new_xsec   = fortFns->evdifcrs();
  float new_weight = (new_xsec/old_xsec);
  //float new_weight = old_weight * (new_xsec/old_xsec);

  if (isinf(new_weight) || isnan(new_weight)) {
    cout << "NReWeightNuXSecNCEL::CalcWeightMa() Warning: new_weight is infinite, setting to 1" << endl;
    new_weight = 1;
  }
  
#ifdef _N_REWEIGHT_NCEL_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << endl;
  cout << "differential cross section (new) = " << new_xsec << endl;
  //cout << "event generation weight = " << old_weight << endl;
  cout << "new weight = " << new_weight << endl;
#endif 

  return new_weight;
}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcChisq()
{
  double chisq = 0;
  chisq += TMath::Power(fMaTwkDial,   2.);
  //chisq += TMath::Power(fEtaTwkDial,  2.);
  return chisq;
}
//_______________________________________________________________________________________
