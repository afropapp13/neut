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
 @ Apr 13, 2011 - PD
   Implemented for NEUT

*/
//____________________________________________________________________________

//#include <NeutRootHandlers.h>

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NReWeightNuXSecCOH.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_COH_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;
//_______________________________________________________________________________________
NReWeightNuXSecCOH::NReWeightNuXSecCOH() 
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuXSecCOH::~NReWeightNuXSecCOH()
{

}
//_______________________________________________________________________________________
void NReWeightNuXSecCOH::Init(void)
{
  fortFns = NFortFns::Instance();

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  // dial =1, reweight the default model (RS) with BS model
  // dial =0, do not change the model
  fNECOHEPIDef = 0;
  fNECOHEPICurr      = fNECOHEPIDef;
  fNECOHEPITwkDial   = 0; 

  fMaTwkDial   = 0; 
  fMaDef       = 0.; 
  fMaCurr      = fMaDef;
  fR0TwkDial   = 0.; 
  fR0Def       = 0.; 
  fR0Curr      = fR0Def;
  fA1TwkDial   = 0.; 
  fA1Def       = 0.; 
  fA1Curr      = fA1Def;
  fb1TwkDial   = 0.; 
  fb1Def       = 0.; 
  fb1Curr      = fb1Def;
}
//_______________________________________________________________________________________
bool NReWeightNuXSecCOH::IsHandled(NSyst_t syst)
{
   bool handle;

   switch(syst) {
     case ( kXSecTwkDial_NECOHEPI ) : 
     case ( kXSecTwkDial_MaCOHpi ) : 
     case ( kXSecTwkDial_R0COHpi ) : 
     case ( kXSecTwkDial_fA1COHpi ) :
     case ( kXSecTwkDial_fb1COHpi ) :
       handle = true;  
       break;
     default:
       handle = false;
       break;
   }
   return handle;
}
//_______________________________________________________________________________________
void NReWeightNuXSecCOH::SetSystematic(NSyst_t syst, double twk_dial)
{
   switch(syst) {
       case ( kXSecTwkDial_NECOHEPI ) : 
       fNECOHEPITwkDial = twk_dial;
       break;
       case ( kXSecTwkDial_MaCOHpi ) : 
       fMaTwkDial = twk_dial;
       break;
       case ( kXSecTwkDial_R0COHpi ) : 
       fR0TwkDial = twk_dial;
       break;
       case ( kXSecTwkDial_fA1COHpi ) : 
       fA1TwkDial = twk_dial;
       break;
       case ( kXSecTwkDial_fb1COHpi ) : 
       fb1TwkDial = twk_dial;
       break;
     default:
       break;
   }
}
//_______________________________________________________________________________________
void NReWeightNuXSecCOH::Reset(void)
{
  fNECOHEPITwkDial = 0.;
  fNECOHEPICurr    = fNECOHEPIDef;
  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;
  fR0TwkDial   = 0.; 
  fR0Curr      = fR0Def;
  fA1TwkDial   = 0.; 
  fA1Curr      = fA1Def;
  fb1TwkDial   = 0.; 
  fb1Curr      = fb1Def;
  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuXSecCOH::Reconfigure(void)
{  
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

//  cout<< "In NReweightNuXSecCOH  fNECOHEPICurr = "<<fNECOHEPICurr<<" fNECOHEPIDef = "<<fNECOHEPIDef <<" fNECOHEPITwkDial = "<< fNECOHEPITwkDial<<endl;
 // if (fNECOHEPICurr != fNECOHEPIDef) fortFns->SetMCDefaultVal(kXSecTwkDial_NECOHEPI, fNECOHEPICurr);
  
  double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCOHpi);
  double fracerr_r0 = fracerr->OneSigmaErr(kXSecTwkDial_R0COHpi);
  double fracerr_a1 = fracerr->OneSigmaErr(kXSecTwkDial_fA1COHpi);
  double fracerr_b1 = fracerr->OneSigmaErr(kXSecTwkDial_fb1COHpi);
  /*
  fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
  fR0Curr = fR0Def * (1. + fR0TwkDial * fracerr_r0);
  fA1Curr = fA1Def * (fA1TwkDial * fracerr_a1);
  fb1Curr = fb1Def * (fb1TwkDial * fracerr_b1);
  */

  fMaCurr = (1. + fMaTwkDial * fracerr_ma);
  fR0Curr = (1. + fR0TwkDial * fracerr_r0);
  fA1Curr = (fA1TwkDial * fracerr_a1);
  fb1Curr = (fb1TwkDial * fracerr_b1);

  double fracerr_necohepi = fracerr->OneSigmaErr(kXSecTwkDial_NECOHEPI);
  //fNECOHEPICurr = int(fNECOHEPIDef + fNECOHEPITwkDial * fracerr_necohepi) ; 
  fNECOHEPICurr = int(fNECOHEPITwkDial * fracerr_necohepi) ; 

  /*
  fNECOHEPICurr = TMath::Max(0.,fNECOHEPICurr );
  fMaCurr = TMath::Max(0., fMaCurr  );
  fR0Curr = TMath::Max(0., fR0Curr  );
  fA1Curr = TMath::Max(0., fA1Curr  );
  fb1Curr = TMath::Max(0., fb1Curr  );
  */
  }
//_______________________________________________________________________________________
double NReWeightNuXSecCOH::CalcWeight() 
{
  bool is_coh = modeDefn.isCOH(nework_.modene);
  if(!is_coh) return 1.;

  bool tweaked =
    (TMath::Abs(fNECOHEPITwkDial) > controls::kASmallNum) || 
    (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
    (TMath::Abs(fR0TwkDial) > controls::kASmallNum) ||
    (TMath::Abs(fA1TwkDial) > controls::kASmallNum) ||
    (TMath::Abs(fb1TwkDial) > controls::kASmallNum); 
  if(!tweaked) return 1.0;

  bool is_cc  = modeDefn.isCC(nework_.modene);
  bool is_nc  = modeDefn.isNC(nework_.modene);
  if(is_cc && !fRewCC) return 1.;
  if(is_nc && !fRewNC) return 1.;

  int nupdg = nework_.ipne[0];
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;
  
  fortFns->SetDefaults();
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_COH_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif

  float old_xsec     = fortFns->evdifcrs();
  //float old_xsec   = event.DiffXSec();
  //float old_weight = event.Weight();
  if (old_xsec==0) {
    cout << "NReWeightNuXSecCOH::CalcWeight() Warning: old_xsec==0, setting weight to 1" << endl;
    return 1;
  }

  fNECOHEPIDef = neutcoh_.necohepi;
  fMaDef = nemdls_.xmacoh;
  fR0Def = nemdls_.rad0nu;
  fA1Def = nemdls_.fa1coh;
  fb1Def = nemdls_.fb1coh;

  neutcoh_.necohepi = int(fNECOHEPICurr + fNECOHEPIDef);
  nemdls_.xmacoh = fMaDef * fMaCurr;
  nemdls_.rad0nu = fR0Def * fR0Curr;
  nemdls_.fa1coh = fA1Def * fA1Curr;
  nemdls_.fb1coh = fb1Def * fb1Curr;

  fNECOHEPICurr = TMath::Max(0.,fNECOHEPICurr );
  fMaCurr = TMath::Max(0., fMaCurr  );
  fR0Curr = TMath::Max(0., fR0Curr  );
  fA1Curr = TMath::Max(0., fA1Curr  );
  fb1Curr = TMath::Max(0., fb1Curr  );//can argue this one isn't needed

  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_COH_DEBUG_
  fortFns->print_allparams();
#endif

  float new_xsec   = fortFns->evdifcrs();
  float new_weight = (new_xsec/old_xsec);
  //float new_weight = old_weight * (new_xsec/old_xsec);

  if (isinf(new_weight) || isnan(new_weight)) {
    cout << "NReWeightNuXSecCOH::CalcWeightMa() Warning: new_weight is infinite, setting to 1" << endl;
    new_weight = 1;
  }

#ifdef _N_REWEIGHT_COH_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << endl;
  cout << "differential cross section (new) = " << new_xsec << endl;
  //cout << "event generation weight = " << old_weight << endl;
  cout << "new weight = " << new_weight << endl;
#endif

  return new_weight;
}
//_______________________________________________________________________________________
double NReWeightNuXSecCOH::CalcChisq()
{
  double chisq = 
    TMath::Power(fMaTwkDial, 2.) +
    TMath::Power(fR0TwkDial, 2.) +
    TMath::Power(fA1TwkDial, 2.) +
    TMath::Power(fb1TwkDial, 2.) ;
  return chisq;
}
//_______________________________________________________________________________________


