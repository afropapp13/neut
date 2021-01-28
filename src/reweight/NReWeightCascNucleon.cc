//____________________________________________________________________________
/*
Based on code for Pion reweighting, to perform nucleon reweighting in place of INuke

 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
    (pion version of the code)

    Toby Nonnenmacher, Nucleon reweighting code

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Apr 6, 2011 - PD
   Implemented for NEUT



*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "NReWeightCascNucleon.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"


//#define _N_REWEIGHT_CASC_DEBUG_
//#define _N_REWEIGHT_CASC_SUPER_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

//_______________________________________________________________________________________
NReWeightCascNucleon::NReWeightCascNucleon()
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightCascNucleon::~NReWeightCascNucleon()
{
}
//_______________________________________________________________________________________
void NReWeightCascNucleon::Init(void)
{
  fortFns = NFortFns::Instance();

  /*  fInelLowDef  = 0.;
  fInelHighDef = 0.;
  fPiProdDef   = 0.;
  fAbsDef      = 0.;
  fCExLowDef   = 0.;
  fCExHighDef  = 0.;
  fAllDef      = 0.;
  
  fInelLowCurr  = fInelLowDef  ;   
  fInelHighCurr = fInelHighDef ;  
  fPiProdCurr   = fPiProdDef   ;      
  fAbsCurr      = fAbsDef      ;   
  fCExLowCurr   = fCExLowDef   ;  
  fCExHighCurr  = fCExHighDef  ; 
  fAllCurr      = fAllDef  ;

  // For 1-sigma change
  fInelLowTwkDial  = 0;   
  fInelHighTwkDial = 0;  
  fPiProdTwkDial   = 0;      
  fAbsTwkDial      = 0;   
  fCExLowTwkDial   = 0;  
  fCExHighTwkDial  = 0; 
  fAllTwkDial      = 0; 
  */
  fTotalProbDef = 0;
  fElasticProbDef = 0;
  fSinglePiProbDef = 0;
  fDoublePiProbDef = 0;

  fTotalProbCurr = 0;
  fElasticProbCurr = 0;
  fSinglePiProbCurr = 0;
  fDoublePiProbCurr = 0;

  fTotalProbTwkDial = 0;
  fElasticProbTwkDial = 0;
  fSinglePiProbTwkDial = 0;
  fDoublePiProbTwkDial = 0;



  // For absolute change
  //fInelLowTwkDial  = fInelLowDef ;   
  //fInelHighTwkDial = fInelHighDef; 
  //fPiProdTwkDial   = fPiProdDef  ; 
  //fAbsTwkDial      = fAbsDef     ; 
  //fCExLowTwkDial   = fCExLowDef  ; 
  //fCExHighTwkDial  = fCExHighDef ; 
}
//_______________________________________________________________________________________
bool NReWeightCascNucleon::IsHandled(NSyst_t syst)
{
   bool handle;



   switch(syst) {

     case ( kCascTwkDial_TotalProb    ) :
     case ( kCascTwkDial_ElasticProb   ) :
     case ( kCascTwkDial_SinglePiProb       ) :
     case ( kCascTwkDial_DoublePiProb    ) :

       //            case ( kCascTwkDial_FrCExHigh_pi   ) :
       //case ( kCascTwkDial_FrInelHigh_pi  ) :
       //case ( kCascTwkDial_All_pi         ) :
     //case ( kINukeTwkDial_MFP_N       ) :
     //case ( kINukeTwkDial_FrCEx_N     ) :
     //case ( kINukeTwkDial_FrElas_N    ) :
     //case ( kINukeTwkDial_FrInel_N    ) :
     //case ( kINukeTwkDial_FrAbs_N     ) :
     //case ( kINukeTwkDial_FrPiProd_N  ) :
          handle = true;
          break;

     default:
          handle = false;
    break;
   }

   return handle;
}
//_______________________________________________________________________________________

void NReWeightCascNucleon::SetSystematic(NSyst_t syst, double twk_dial)
{
  //4 dials for nucleon reweighting: total/elastic/single pi/double pi
  std::cout << "start of SetSystematic" <<std::endl;
  switch(syst) {
  case ( kCascTwkDial_TotalProb   ) :
    std::cout << "case kCascTwkDial_TotalProb" << std::endl;
    fTotalProbTwkDial   =  twk_dial ;  
    std::cout << "twk_dial for total prob =" << twk_dial;
    break;

  case ( kCascTwkDial_ElasticProb   ) :
    std::cout << "case kCascTwkDial_ElasticProb" << std::endl;
    fElasticProbTwkDial  =  twk_dial ;   
    break;

  case ( kCascTwkDial_SinglePiProb      ) :
    std::cout << "case kCascTwkDial_SinglePiProb" << std::endl;
    fSinglePiProbTwkDial      =  twk_dial ;   
    break;

  case ( kCascTwkDial_DoublePiProb    ) :
    std::cout << "case kCascTwkDial_DoublePiProb" << std::endl;
    fDoublePiProbTwkDial   =  twk_dial;      
    break;

    /*
  case ( kCascTwkDial_DoublePiProb   ) :
    fCExHighTwkDial  =  twk_dial ; 
    break;

  case ( kCascTwkDial_FrInelHigh_pi  ) :
    fInelHighTwkDial =  twk_dial  ;  
    break;

  case ( kCascTwkDial_All_pi  ) :
    fAllTwkDial =  twk_dial  ;  
    break;
    */
  default:
    break;
  }
}
//_______________________________________________________________________________________
void NReWeightCascNucleon::Reset(void)
{




  fTotalProbCurr = fTotalProbDef;
  fElasticProbCurr = fElasticProbDef;
  fSinglePiProbCurr = fSinglePiProbDef;
  fDoublePiProbCurr = fDoublePiProbDef;
  
  
  //For 1 sigma changes
  fTotalProbTwkDial = 0;
  fElasticProbTwkDial = 0;
  fSinglePiProbTwkDial = 0;
  fDoublePiProbTwkDial = 0;


  //For absolute changes
  //fTotalProbTwkDial = fTotalProbDef;
  //fElasticProbTwkDial = fElasticProbDef;
  //fSinglePiProbTwkDial = fSinglePiProbDef;
  //  fDoublePiProbTwkDial = fDoublePiProbDef;




  //fInelLowCurr  = fInelLowDef ;   
  //fInelHighCurr = fInelHighDef;  
  //fPiProdCurr   = fPiProdDef  ;      
  //fAbsCurr      = fAbsDef     ;   
  //fCExLowCurr   = fCExLowDef  ;  
  //fCExHighCurr  = fCExHighDef ; 
  //fAllCurr      = fAllDef ; 

  // For 1-sigma change
  //fInelLowTwkDial  = 0;   
  //fInelHighTwkDial = 0;  
  //fPiProdTwkDial   = 0;      
  //fAbsTwkDial      = 0;   
  //fCExLowTwkDial   = 0;  
  //fCExHighTwkDial  = 0;
  //fAllTwkDial      = 0;

  // For absolute change
  //fInelLowTwkDial  = fInelLowDef ;   
  //fInelHighTwkDial = fInelHighDef; 
  //fPiProdTwkDial   = fPiProdDef  ; 
  //fAbsTwkDial      = fAbsDef     ; 
  //fCExLowTwkDial   = fCExLowDef  ; 
  //fCExHighTwkDial  = fCExHighDef ; 


 
  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightCascNucleon::Reconfigure(void)
{
  fTotalProbDef = 1;
  fElasticProbDef = 1;
  fSinglePiProbDef=1;
  fDoublePiProbDef=1;
  std::cout << "ReconfTotalProbDef " << fTotalProbDef << std::endl;

  // For 1-sigma change
  NSystUncertainty * fracerr = NSystUncertainty::Instance();
  //  std::cout << "fracerr = " << fracerr << std::endl;

  int    sign_TotProb = utils::rew::Sign(fTotalProbTwkDial);
  std::cout << "ReconfkCascTwkDial_TotalProb " << kCascTwkDial_TotalProb << std::endl;
  std::cout << "Reconsign_TotProb " << sign_TotProb << std::endl;

  double fracerr_TotProb = fracerr->OneSigmaErr(kCascTwkDial_TotalProb, sign_TotProb);
  std::cout << "ReconfTotalProbTwkDial " << fTotalProbTwkDial << std::endl;
  std::cout << "ReconfracerrTotProb " << fracerr_TotProb << std::endl;
  
  //    fracerr_TotProb = 0.0;
  
  
  fTotalProbCurr = fTotalProbDef * (1. + fTotalProbTwkDial * fracerr_TotProb);
  fTotalProbCurr = TMath::Max(0., fTotalProbCurr  );
  std::cout << "fTotalProbCurr = " << fTotalProbCurr << std::endl;
  std::cout << "fTotalProbTwkDial = " << fTotalProbTwkDial << std::endl;
  std::cout << "fTotalProbDef = " << fTotalProbDef << std::endl;

  int    sign_ElasticProb = utils::rew::Sign(fElasticProbTwkDial);
  double fracerr_ElasticProb = fracerr->OneSigmaErr(kCascTwkDial_ElasticProb, sign_ElasticProb);
  //fracerr_ElasticProb = 0.3;

  fElasticProbCurr = fElasticProbDef * (1. + fElasticProbTwkDial * fracerr_ElasticProb);
  fElasticProbCurr = TMath::Max(0., fElasticProbCurr  );


  int    sign_SinglePiProb = utils::rew::Sign(fSinglePiProbTwkDial);
  double fracerr_SinglePiProb = fracerr->OneSigmaErr(kCascTwkDial_SinglePiProb, sign_SinglePiProb);
  fSinglePiProbCurr = fSinglePiProbDef * (1. + fSinglePiProbTwkDial * fracerr_SinglePiProb);
  fSinglePiProbCurr = TMath::Max(0., fSinglePiProbCurr  );


  int    sign_DoublePiProb = utils::rew::Sign(fDoublePiProbTwkDial);
  double fracerr_DoublePiProb = fracerr->OneSigmaErr(kCascTwkDial_DoublePiProb, sign_DoublePiProb);
  fDoublePiProbCurr = fDoublePiProbDef * (1. + fDoublePiProbTwkDial * fracerr_DoublePiProb);
  fDoublePiProbCurr = TMath::Max(0., fDoublePiProbCurr  );

  /*
  int    sign_InelHightwk = utils::rew::Sign(fInelHighTwkDial);
  double fracerr_InelHigh = fracerr->OneSigmaErr(kCascTwkDial_FrInelHigh_pi, sign_InelHightwk);
  fInelHighCurr = fInelHighDef * (1. + fInelHighTwkDial * fracerr_InelHigh);
  fInelHighCurr   = TMath::Max(0., fInelHighCurr  );
  
  int    sign_PiProdtwk = utils::rew::Sign(fPiProdTwkDial);
  double fracerr_PiProd = fracerr->OneSigmaErr(kCascTwkDial_FrPiProd_pi, sign_PiProdtwk);
  fPiProdDef * (1. + fPiProdTwkDial * fracerr_PiProd);
  fPiProdCurr   = TMath::Max(0., fPiProdCurr  );
  
  int    sign_Abstwk = utils::rew::Sign(fAbsTwkDial);
  double fracerr_Abs = fracerr->OneSigmaErr(kCascTwkDial_FrAbs_pi, sign_Abstwk);
  fAbsCurr = fAbsDef * (1. + fAbsTwkDial * fracerr_Abs);
  fAbsCurr   = TMath::Max(0., fAbsCurr  );
  
  int    sign_CExLowtwk = utils::rew::Sign(fCExLowTwkDial);
  double fracerr_CExLow = fracerr->OneSigmaErr(kCascTwkDial_FrCExLow_pi, sign_CExLowtwk);
  fCExLowCurr = fCExLowDef * (1. + fCExLowTwkDial * fracerr_CExLow);
  fCExLowCurr   = TMath::Max(0., fCExLowCurr  );  

  int    sign_CExHightwk = utils::rew::Sign(fCExHighTwkDial);
  double fracerr_CExHigh = fracerr->OneSigmaErr(kCascTwkDial_FrCExHigh_pi, sign_CExHightwk);
  fCExHighDef * (1. + fCExHighTwkDial * fracerr_CExHigh);
  fCExHighCurr   = TMath::Max(0., fCExHighCurr  );

  int    sign_Alltwk = utils::rew::Sign(fAllTwkDial);
  double fracerr_All = fracerr->OneSigmaErr(kCascTwkDial_All_pi, sign_Alltwk);
  fAllCurr = fAllDef * (1. + fAllTwkDial * fracerr_All);
  fAllCurr   = TMath::Max(0., fAllCurr  );
  */

  // For absolute change
  //fInelLowCurr =  fInelLowTwkDial ;
  //
  //fInelHighCurr = fInelHighTwkDial;
  //
  //fPiProdCurr =  fPiProdTwkDial   ;
  //
  //fAbsCurr = fAbsTwkDial          ;
  //
  //fCExLowCurr =  fCExLowTwkDial   ;
  //
  //fCExHighCurr = fCExHighTwkDial;  

}
//_______________________________________________________________________________________
double NReWeightCascNucleon::CalcWeight() 
{ 

  bool tweaked = ( (TMath::Abs(fTotalProbTwkDial) > controls::kASmallNum) ||
                   (TMath::Abs(fElasticProbTwkDial) > controls::kASmallNum) ||
                   (TMath::Abs(fSinglePiProbTwkDial) > controls::kASmallNum) ||
                   (TMath::Abs(fDoublePiProbTwkDial) > controls::kASmallNum) );

  #ifdef _N_REWEIGHT_CASC_DEBUG_
    cout << "Starting NReWeightCascNucleon::CalcWeight, nrint_.pcascprob = " << nrint_.pcascprob << endl;
    cout << "fsihist_.fsiprob is: " << fsihist_.fsiprob << endl;
  #endif


  // For absolute change
  //bool tweaked = ( fInelLowCurr != fInelLowDef ||
  //       fInelHighCurr != fInelHighDef ||
  //       fAbsCurr != fAbsDef ||
  //       fPiProdCurr != fPiProdDef ||
  //       fCExLowCurr != fCExLowDef ||
  //       fCExHighCurr != fCExHighDef );  
 if(!tweaked) {
    std::cout << "no tweaked" << std::endl;
    return 1.0;
  }

  fortFns->SetDefaults();
  fortFns->Reconfigure();


  #ifdef _N_REWEIGHT_CASC_SUPER_DEBUG_
    fortFns->print_allevent();
    fortFns->print_allparams();
  #endif 

  //weight will be new_xsec/old_xsec
  float old_xsec = nrint_.pcascprob;

  // Toby's Macro and NUISANCE read different values for fsiprob ... NUISANCE often triggers
  // this line. For the moment I'm going to comment it out, assuming that Toby used his macro
  // to show sensible regen vs reweight results.  
  // if (fsihist_.fsiprob == 1) return 1;

  #ifdef _N_REWEIGHT_CASC_DEBUG_
    cout << "Filling old_xsec with nrint_.pcascprob = " << nrint_.pcascprob << endl;
    cout << "fsihist_.fsiprob is: " << fsihist_.fsiprob << endl;
    if (fabs( fsihist_.fsiprob - old_xsec ) > 0.00005)
      cout << "NReWeightCascNucleon() Error: Previously calculated FSIPROB inconsistent" << endl;
  #endif

  if (old_xsec<=0) {
    cout << "NReWeightCascNucleon() Warning: evpiprob old_xsec <= 0, returning weight = 1" << endl;
    //return 1;
  }

  
  nucres_.xnucfact = 1;
  nucres_.xnucelafact =1;
  nucres_.xnucspifact=1;
  nucres_.xnucdpifact=1;

  fTotalProbDef = nucres_.xnucfact;
  fElasticProbDef = nucres_.xnucelafact;
  fSinglePiProbDef = nucres_.xnucspifact;
  fDoublePiProbDef = nucres_.xnucdpifact;

  nucres_.xnucfact = TMath::Max(float(0.),float(fTotalProbCurr*fTotalProbCurr));
  nucres_.xnucelafact = TMath::Max(float(0.),float(fElasticProbCurr*fElasticProbCurr));
  nucres_.xnucspifact = TMath::Max(float(0.),float(fSinglePiProbCurr*fSinglePiProbCurr));
  nucres_.xnucdpifact = TMath::Max(float(0.),float(fDoublePiProbCurr*fDoublePiProbCurr));
  
  //set up to perform reweighting
  fortFns->Reconfigure();
  //reinitialise pcascprob to 1 (start of cascade)
  nrint_.pcascprob = 1;

  #ifdef _N_REWEIGHT_CASC_DEBUG_
    cout << "Reconfigured to calc new_xsec" << endl;
    cout << "Calling nrprton" << endl;
  #endif
  
  //rerun nrprton.F to get pcascprob for reweight
  fortFns->nrprton(); 

  #ifdef _N_REWEIGHT_CASC_SUPER_DEBUG_
    fortFns->print_allparams();
  #endif

  //set new_xsec
  float new_xsec  = nrint_.pcascprob;

  //calculate weight
  float new_weight = 1.0; 
  if (old_xsec<=0) new_weight = 1.0;  
  else new_weight = (new_xsec/old_xsec);
  

  //if there was no cascade, just set weight to 1:
  if(old_xsec == 1) new_weight =1;
  
  #ifdef _N_REWEIGHT_CASC_DEBUG_
    cout << "Filling new_xsec with nrint_.pcascprob = " << nrint_.pcascprob << endl;
    cout << "fsihist_.fsiprob is: " << fsihist_.fsiprob << endl;
  #endif

  if (isinf(new_weight) || isnan(new_weight) || new_weight<0) {
    cout << "NReWeightCascNucleon::CalcWeight() Warning: new_weight is infinite or negative (" << new_weight << "), setting to -9999" << endl;
    new_weight = -9999;
  }

  return new_weight;
}
//_______________________________________________________________________________________
double NReWeightCascNucleon::CalcChisq(void)
{
  double chisq = 0.;


  chisq += TMath::Power(fTotalProbTwkDial, 2.);
  chisq += TMath::Power(fElasticProbTwkDial, 2.);
  chisq += TMath::Power(fSinglePiProbTwkDial, 2.);
  chisq += TMath::Power(fDoublePiProbTwkDial, 2.);

  return chisq;
}
//_______________________________________________________________________________________
vector<double> NReWeightCascNucleon::GetCurrParVals(void)
{
  vector<double> parVals;


  parVals.push_back(fTotalProbCurr);
  parVals.push_back(fElasticProbCurr);
  parVals.push_back(fSinglePiProbCurr);
  parVals.push_back(fDoublePiProbCurr); 
  return parVals;  
}
//_______________________________________________________________________________________
