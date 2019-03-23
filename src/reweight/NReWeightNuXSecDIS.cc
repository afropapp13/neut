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
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code. 
   First included in v2.5.1.
 @ May 17, 2010 - CA
   Code extracted from NReWeightNuXSec and redeveloped in preparation for 
   the Summer 2010 T2K analyses.
 @ Oct 20, 2010 - CA
   Make static consts kModeABCV12u and kModeABCV12uShape public so as to
   aid external configuration.
 @ Apr 13, 2011 - PD
   Implemented for NEUT
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NReWeightNuXSecDIS.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_DIS_DEBUG_

using namespace neut;
using namespace neut::rew;
 
using std::cout;

//_______________________________________________________________________________________
NReWeightNuXSecDIS::NReWeightNuXSecDIS() 
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuXSecDIS::~NReWeightNuXSecDIS()
{
}
//_______________________________________________________________________________________
void NReWeightNuXSecDIS::Init(void)
{
  fortFns = NFortFns::Instance();

  this->SetMode(0);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;

  fBYOnOffTwkDial = 1;  // On by default
  fBYOnOffDef     = fortFns->NEBODEKdef;
  fBYOnOffCurr    = fBYOnOffDef;
}
//_______________________________________________________________________________________
bool NReWeightNuXSecDIS::IsHandled(NSyst_t syst)
{
   bool handle = false;

   switch(syst) {

   case ( kXSecTwkDial_NormDIS ) : 
     handle = true;
     break;

   case ( kXSecTwkDial_BYOnOffDIS ) : 
     handle = true;
     break;

   default:
     handle = false;
     break;
   }

   return handle;
}
//_______________________________________________________________________________________
void NReWeightNuXSecDIS::SetSystematic(NSyst_t syst, double twk_dial)
{
   if(! this->IsHandled(syst)) return;

   switch(syst) {

   case ( kXSecTwkDial_NormDIS ) :
     fNormTwkDial = twk_dial;
     break;

   case ( kXSecTwkDial_BYOnOffDIS ) :
     fBYOnOffTwkDial = twk_dial;
     break;

   default:
     return;
     break;
   }
}
//_______________________________________________________________________________________
void NReWeightNuXSecDIS::Reset(void)
{
  fNormTwkDial = 0.;        

  fBYOnOffTwkDial = 1;
  fBYOnOffCurr = fBYOnOffDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuXSecDIS::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormDIS);
  fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
  fNormCurr = TMath::Max(0., fNormCurr);

  fBYOnOffCurr = fBYOnOffTwkDial;

}
//_______________________________________________________________________________________
double NReWeightNuXSecDIS::CalcWeight() 
{
  bool is_dis = modeDefn.isDIS(nework_.modene);
  bool is_mpi = modeDefn.isMPI(nework_.modene);

  if(fMode==0 && !is_dis && !is_mpi) return 1.;
  if(fMode==1 && !is_dis) return 1.;
  if(fMode==2 && !is_mpi) return 1.;

  bool is_cc  = modeDefn.isCC(nework_.modene);
  bool is_nc  = modeDefn.isNC(nework_.modene);
  if(is_cc && !fRewCC) return 1.;
  if(is_nc && !fRewNC) return 1.;

  int nupdg = nework_.ipne[0];
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  double wght = this->CalcWeightNorm()*this->CalcWeightBYOnOff();
  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuXSecDIS::CalcWeightNorm() 
{
  bool tweaked = 
    (TMath::Abs(fNormTwkDial)   > controls::kASmallNum) ;
  if(!tweaked) return 1.0;

  double wght = fNormCurr;
  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuXSecDIS::CalcWeightBYOnOff() 
{
  bool tweaked = 
    (TMath::Abs(fBYOnOffTwkDial-1)   > controls::kASmallNum); // Nominal is 1
  if(!tweaked) return 1.0;

  fortFns->SetDefaults();
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_DIS_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif

  float old_xsec     = fortFns->evdifcrs();
  if (old_xsec==0) {
    return 1.;
  }
  
  neutdis_.nebodek = 0; //fBYOnOffCurr;

  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_DIS_DEBUG_
  fortFns->print_allparams();
#endif

  float new_xsec   = fortFns->evdifcrs();

#ifdef _N_REWEIGHT_DIS_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << '\n';
  cout << "differential cross section (new) = " << new_xsec << '\n';
#endif 

  float new_weight = (new_xsec/old_xsec);

  if (isinf(new_weight) || isnan(new_weight)) {
#ifdef _N_REWEIGHT_DIS_DEBUG_
    cout << "NReWeightNuXSecDIS::CalcWeightBYOnOff() Warning: new_weight is infinite, setting to 1" << '\n';
#endif
    return 1;
  }

  // Linear interpolation
  double wght = TMath::Max(0., (1-new_weight)*fBYOnOffCurr + new_weight);
    
  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuXSecDIS::CalcChisq(void)
{
  double chisq = 0.; 
      chisq += TMath::Power(fNormTwkDial,  2.);
  return chisq;
}
//_______________________________________________________________________________________


