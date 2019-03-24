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
#include "NReWeightCasc.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_CASC_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

//_______________________________________________________________________________________
NReWeightCasc::NReWeightCasc()
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightCasc::~NReWeightCasc()
{
}
//_______________________________________________________________________________________
void NReWeightCasc::Init(void)
{
  fortFns = NFortFns::Instance();

  fInelLowDef  = fortFns->FEFQEdef   ;   
  fInelHighDef = fortFns->FEFQEHdef  ;  
  fPiProdDef   = fortFns->FEFINELdef ;      
  fAbsDef      = fortFns->FEFABSdef  ;   
  fCExLowDef   = fortFns->FEFCXdef   ;  
  fCExHighDef  = fortFns->FEFCXHdef  ; 
  fAllDef  = fortFns->FEFALLdef  ; 

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


  // For absolute change
  //fInelLowTwkDial  = fInelLowDef ;   
  //fInelHighTwkDial = fInelHighDef; 
  //fPiProdTwkDial   = fPiProdDef  ; 
  //fAbsTwkDial      = fAbsDef     ; 
  //fCExLowTwkDial   = fCExLowDef  ; 
  //fCExHighTwkDial  = fCExHighDef ; 

}
//_______________________________________________________________________________________
bool NReWeightCasc::IsHandled(NSyst_t syst)
{
   bool handle;

   switch(syst) {
     case ( kCascTwkDial_FrCExLow_pi    ) :
     case ( kCascTwkDial_FrInelLow_pi   ) :
     case ( kCascTwkDial_FrAbs_pi       ) :
     case ( kCascTwkDial_FrPiProd_pi    ) :
     case ( kCascTwkDial_FrCExHigh_pi   ) :
     case ( kCascTwkDial_FrInelHigh_pi  ) :
     case ( kCascTwkDial_All_pi         ) :
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
void NReWeightCasc::SetSystematic(NSyst_t syst, double twk_dial)
{
  switch(syst) {
  case ( kCascTwkDial_FrCExLow_pi    ) :
    fCExLowTwkDial   =  twk_dial ;  
    break;

  case ( kCascTwkDial_FrInelLow_pi   ) :
    fInelLowTwkDial  =  twk_dial ;   
    break;

  case ( kCascTwkDial_FrAbs_pi       ) :
    fAbsTwkDial      =  twk_dial ;   
    break;

  case ( kCascTwkDial_FrPiProd_pi    ) :
    fPiProdTwkDial   =  twk_dial ;      
    break;

  case ( kCascTwkDial_FrCExHigh_pi   ) :
    fCExHighTwkDial  =  twk_dial ; 
    break;

  case ( kCascTwkDial_FrInelHigh_pi  ) :
    fInelHighTwkDial =  twk_dial  ;  
    break;

  case ( kCascTwkDial_All_pi  ) :
    fAllTwkDial =  twk_dial  ;  
    break;

  default:
    break;
  }
}
//_______________________________________________________________________________________
void NReWeightCasc::Reset(void)
{
  fInelLowCurr  = fInelLowDef ;   
  fInelHighCurr = fInelHighDef;  
  fPiProdCurr   = fPiProdDef  ;      
  fAbsCurr      = fAbsDef     ;   
  fCExLowCurr   = fCExLowDef  ;  
  fCExHighCurr  = fCExHighDef ; 
  fAllCurr      = fAllDef ; 

  // For 1-sigma change
  fInelLowTwkDial  = 0;   
  fInelHighTwkDial = 0;  
  fPiProdTwkDial   = 0;      
  fAbsTwkDial      = 0;   
  fCExLowTwkDial   = 0;  
  fCExHighTwkDial  = 0;
  fAllTwkDial      = 0;

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
void NReWeightCasc::Reconfigure(void)
{

  // For 1-sigma change
  NSystUncertainty * fracerr = NSystUncertainty::Instance();
  
  int    sign_InelLowtwk = utils::rew::Sign(fInelLowTwkDial);
  double fracerr_InelLow = fracerr->OneSigmaErr(kCascTwkDial_FrInelLow_pi, sign_InelLowtwk);
  fInelLowCurr = fInelLowDef * (1. + fInelLowTwkDial * fracerr_InelLow);
  fInelLowCurr   = TMath::Max(0., fInelLowCurr  );
  
  int    sign_InelHightwk = utils::rew::Sign(fInelHighTwkDial);
  double fracerr_InelHigh = fracerr->OneSigmaErr(kCascTwkDial_FrInelHigh_pi, sign_InelHightwk);
  fInelHighCurr = fInelHighDef * (1. + fInelHighTwkDial * fracerr_InelHigh);
  fInelHighCurr   = TMath::Max(0., fInelHighCurr  );
  
  int    sign_PiProdtwk = utils::rew::Sign(fPiProdTwkDial);
  double fracerr_PiProd = fracerr->OneSigmaErr(kCascTwkDial_FrPiProd_pi, sign_PiProdtwk);
  fPiProdCurr = fPiProdDef * (1. + fPiProdTwkDial * fracerr_PiProd);
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
  fCExHighCurr = fCExHighDef * (1. + fCExHighTwkDial * fracerr_CExHigh);
  fCExHighCurr   = TMath::Max(0., fCExHighCurr  );

  int    sign_Alltwk = utils::rew::Sign(fAllTwkDial);
  double fracerr_All = fracerr->OneSigmaErr(kCascTwkDial_All_pi, sign_Alltwk);
  fAllCurr = fAllDef * (1. + fAllTwkDial * fracerr_All);
  fAllCurr   = TMath::Max(0., fAllCurr  );

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
double NReWeightCasc::CalcWeight() 
{ 
  // For 1-sigma change
  bool tweaked = ( (TMath::Abs(fInelLowTwkDial) > controls::kASmallNum) ||
  		   (TMath::Abs(fInelHighTwkDial) > controls::kASmallNum) ||
  		   (TMath::Abs(fAbsTwkDial) > controls::kASmallNum) ||
  		   (TMath::Abs(fPiProdTwkDial) > controls::kASmallNum) ||
  		   (TMath::Abs(fCExLowTwkDial) > controls::kASmallNum) ||
  		   (TMath::Abs(fCExHighTwkDial) > controls::kASmallNum) ||
  		   (TMath::Abs(fAllTwkDial) > controls::kASmallNum) );


  //cout << fInelLowCurr << " " <<  fInelLowDef    << endl;
  //cout << fInelHighCurr << " " <<  fInelHighDef  << endl;
  //cout << fAbsCurr << " " <<  fAbsDef  		 << endl;
  //cout << fPiProdCurr << " " <<  fPiProdDef 	 << endl;  
  //cout << fCExLowCurr << " " <<  fCExLowDef 	 << endl;  
  //cout << fCExHighCurr << " " <<  fCExHighDef    << endl;

  // For absolute change
  //bool tweaked = ( fInelLowCurr != fInelLowDef ||
  //		   fInelHighCurr != fInelHighDef ||
  //		   fAbsCurr != fAbsDef ||
  //		   fPiProdCurr != fPiProdDef ||
  //		   fCExLowCurr != fCExLowDef ||
  //		   fCExHighCurr != fCExHighDef );  
  if(!tweaked) return 1.0;


  fortFns->SetDefaults();
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_CASC_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif 

  // Not bound
  if (posinnuc_.ibound == 0) return 1;
  
  float old_xsec;

  // Assume FSIPROB <= 0 if not filled properly
  // Currently in ND280 MCP4 only (set in ${T2KREWEIGHT}/src/T2KNeutUtils.cxx)
  if (fsihist_.fsiprob <= 0) {  
    old_xsec = fortFns->evpiprob();
  
  // evpiprob not passing properly with g77 64-bit so grab from common block
#if defined(f2cFortran)&&!defined(gFortran)
    old_xsec = fsihist_.fsiprob;
#endif
  }
  
  // Following should work for SK >=11b (NEUT >= 5.1.2) and neutroot/piscat

  // FSIPROB was not pre-calculated properly during event generation
  else if (fsihist_.fsiprob == 1) 
    return 1;
  
  // Good FSIPROB pre-calculation
  else
    old_xsec = fsihist_.fsiprob;

  if (old_xsec<=0) {
    cout << "NReWeightCasc() Warning: evpiprob old_xsec <= 0, returning weight = 1" << endl;
    return 1;
  }


#ifdef _N_REWEIGHT_CASC_DEBUG_
  cout << "FSI Probability (old) = " << old_xsec << endl;
  fortFns->evpiprob();
  cout << "FSI Probability (old, evpiprob) = " << fsihist_.fsiprob << endl;

  if (fabs( fsihist_.fsiprob - old_xsec ) > 0.00005) {
    cout << "NReWeightCasc() Error: Previously calculated FSIPROB inconsistent" << endl;
    //int cin_tmp;
    //std::cin >> cin_tmp;
  }
#endif

  neffpr_.fefqe   = (float)fInelLowCurr;
  neffpr_.fefqeh  = (float)fInelHighCurr;
  neffpr_.fefinel = (float)fPiProdCurr;
  neffpr_.fefabs  = (float)fAbsCurr;
  neffpr_.fefcx   = (float)fCExLowCurr;
  neffpr_.fefcxh  = (float)fCExHighCurr;  
  neffpr_.fefall  = (float)fAllCurr;  

  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_CASC_DEBUG_
  fortFns->print_allparams();
#endif
  
  //cout << "NReWeightCasc(): evpiprob = " <<  fortFns->evpiprob() << ", FSIPROB= " << fsihist_.fsiprob << endl;
  float new_xsec   = fortFns->evpiprob();

  // evpiprob not passing properly with g77 64-bit so grab from common block
#if defined(f2cFortran)&&!defined(gFortran)
  //if (new_xsec) {
  //cout << "NReWeightCasc() Error: new_xsec=" << new_xsec << " when 0 expected for g77 compiled Fortran code on 64-bit machine" << endl;
  //exit (-1);
  //}
  new_xsec = fsihist_.fsiprob;
#endif
  float new_weight = (new_xsec/old_xsec);

#ifdef _N_REWEIGHT_CASC_DEBUG_
  cout << "FSI Probability (new) = " << new_xsec << endl;
#endif

  if (isinf(new_weight) || isnan(new_weight) || new_weight<0) {
    cout << "NReWeightCasc::CalcWeight() Warning: new_weight is infinite or negative (" << new_weight << "), setting to 1" << endl;
    new_weight = 1;
  }

#ifdef _N_REWEIGHT_CASC_DEBUG_  
  cout << "new weight = " << new_weight << endl;
#endif

  return new_weight;
}
//_______________________________________________________________________________________
double NReWeightCasc::CalcChisq(void)
{
  double chisq = 0.;

  // For 1-sigma change
  chisq += TMath::Power(fInelLowTwkDial, 2.);
  chisq += TMath::Power(fInelHighTwkDial, 2.);
  chisq += TMath::Power(fAbsTwkDial, 2.);
  chisq += TMath::Power(fPiProdTwkDial, 2.);
  chisq += TMath::Power(fCExLowTwkDial, 2.);
  chisq += TMath::Power(fCExHighTwkDial, 2.);
  chisq += TMath::Power(fAllTwkDial, 2.);

  return chisq;
}
//_______________________________________________________________________________________
vector<double> NReWeightCasc::GetCurrParVals(void)
{
  vector<double> parVals;
  // This must be in the same order as ${T2KREWEIGHT}/src/T2KSyst.h
  parVals.push_back(fAbsCurr);   
  parVals.push_back(fInelLowCurr); 
  parVals.push_back(fCExLowCurr);  
  parVals.push_back(fInelHighCurr);
  parVals.push_back(fCExHighCurr);   
  parVals.push_back(fPiProdCurr);	   
  parVals.push_back(fAllCurr);
  return parVals;  
}
//_______________________________________________________________________________________
