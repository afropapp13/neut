//____________________________________________________________________________
/*
Based on code for Pion reweighting, to perform nucleon reweighting in place of INuke

 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

	  Toby Nucleon

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

  fInelLowDef  = 0.;
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

  std::cout <<"fTotalProbTwkDial = "<< TMath::Abs(fTotalProbTwkDial) << " fElasticProbTwkDial = " << TMath::Abs(fElasticProbTwkDial)<<" fSinglePiProbTwkDial= " << TMath::Abs(fSinglePiProbTwkDial)<<" fDoublePiProbTwkDial= " << TMath::Abs(fDoublePiProbTwkDial) << std::endl;
  // For 1-sigma change  
bool tweaked = ( (TMath::Abs(fTotalProbTwkDial) > controls::kASmallNum) ||
		 (TMath::Abs(fElasticProbTwkDial) > controls::kASmallNum) ||
		 (TMath::Abs(fSinglePiProbTwkDial) > controls::kASmallNum) ||
		 (TMath::Abs(fDoublePiProbTwkDial) > controls::kASmallNum) );
// std::cout << "nfnstep calcweight = " << nucleonfsihist_.nfnstep << std::endl;


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
 if(!tweaked) 
   {
     std::cout << "no tweaked" << std::endl;
     return 1.0;
   }

 fortFns->SetDefaults();


 fortFns->Reconfigure();

#ifdef _N_REWEIGHT_CASC_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif 

  //maybe comment out
  // Not bound
  //  if (posinnuc_.ibound == 0)
  //  {
  //    std::cout << "not bound" << std::endl;
  //   return 1;
  // }

  float old_xsec;

  //cout << "NReWeightCascNucleon FsiProb = " << fsihist_.fsiprob << endl;

  // Assume FSIPROB <= 0 if not filled properly
  // Currently in ND280 MCP4 only (set in ${T2KREWEIGHT}/src/T2KNeutUtils.cxx)


    std::cout << "fsihist_.fsiprob = " << fsihist_.fsiprob << std::endl;
  if (fsihist_.fsiprob <= 0) {  

    //    fortFns->nrprton();
    //    fortFns->nrnuc();
    //    fortFns->evpiprob_();
    std::cout << "nrint_.pcascprob = " << nrint_.pcascprob << std::endl;
    
    old_xsec = nrint_.pcascprob;
  // evpiprob not passing properly with g77 64-bit so grab from common block
#if defined(f2cFortran)&&!defined(gFortran)
    old_xsec = nrint_.pcascprob;
#endif
  }
  
  // Following should work for SK >=11b (NEUT >= 5.1.2) and neutroot/piscat

  // FSIPROB was not pre-calculated properly during event generation
  else if (fsihist_.fsiprob == 1) 
    {
      //    std::cout << "fsihist_.fsiprob = 1" << std::endl;
    return 1;
    }
  // Good FSIPROB pre-calculation
  else
    //std::cout << "nrint_.pcascprob = " << nrint_.pcascprob << std::endl;
    //std::cout << "nrint_.pnuccounter = " << nrint_.pnuccounter << std::endl;
    //std::cout << "nrint_.survivalcounter = " << nrint_.survivalcounter << std::endl;
    old_xsec = nrint_.pcascprob;

  if (old_xsec<=0) {
    cout << "NReWeightCascNucleon() Warning: evpiprob old_xsec <= 0, returning weight = 1" << endl;
    return 1;
  }


#ifdef _N_REWEIGHT_CASC_DEBUG_
  cout << "FSI Probability (old) = " << old_xsec << endl;
  //  int dummy = 8;
  //fortFns->nrnuc(&dummy);
  cout << "FSI Probability (old, evpiprob) = " << fsihist_.fsiprob << endl;

  if (fabs( fsihist_.fsiprob - old_xsec ) > 0.00005) {
    cout << "NReWeightCascNucleon() Error: Previously calculated FSIPROB inconsistent" << endl;
    //int cin_tmp;
    //std::cin >> cin_tmp;
  }
#endif

  fInelLowDef = neffpr_.fefqe;
  fInelHighDef= neffpr_.fefqeh;
  fPiProdDef  = neffpr_.fefinel;
  fAbsDef     = neffpr_.fefabs;
  fCExLowDef  = neffpr_.fefcx;
  fCExHighDef = neffpr_.fefcxh;
  fAllDef     = neffpr_.fefall;

   
  //  std::cout << "xnucfact in NReWeightCascNucleon = " << nucres_.xnucfact << std::endl;
  //std::cout << "xnucelafact in NReWeightCascNucleon = " << nucres_.xnucelafact << std::endl;
  //std::cout << "xnucspifact in NReWeightCascNucleon = " << nucres_.xnucspifact << std::endl;
  //std::cout << "xnucdpifact in NReWeightCascNucleon = " << nucres_.xnucdpifact << std::endl;

  nucres_.xnucfact = 1;
  nucres_.xnucelafact =1;
  nucres_.xnucspifact=1;
  nucres_.xnucdpifact=1;

    fTotalProbDef = nucres_.xnucfact;
  fElasticProbDef = nucres_.xnucelafact;
  fSinglePiProbDef = nucres_.xnucspifact;
  fDoublePiProbDef = nucres_.xnucdpifact;

  // std::cout << "fTotalProbDef " << fTotalProbDef << std::endl;

  
   std::cout << "2xnucfact in NReWeightCascNucleon = " << nucres_.xnucfact << std::endl;
  //std::cout << "2xnucelafact in NReWeightCascNucleon = " << nucres_.xnucelafact << std::endl;
  //std::cout << "2xnucspifact in NReWeightCascNucleon = " << nucres_.xnucspifact << std::endl;
  //std::cout << "2xnucdpifact in NReWeightCascNucleon = " << nucres_.xnucdpifact << std::endl; 
  neffpr_.fefqe   = TMath::Max(float(0.),float(fInelLowCurr*fInelLowCurr));
  neffpr_.fefqeh  = TMath::Max(float(0.),float(fInelHighCurr*fInelHighCurr));
  neffpr_.fefinel = TMath::Max(float(0.),float(fPiProdCurr*fPiProdCurr));
  neffpr_.fefabs  = TMath::Max(float(0.),float(fAbsCurr*fAbsCurr));
  neffpr_.fefcx   = TMath::Max(float(0.),float(fCExLowCurr*fCExLowCurr));
  neffpr_.fefcxh  = TMath::Max(float(0.),float(fCExHighCurr*fCExHighCurr));  
  //  neffpr_.feffall = TMath::Max(float(0.),float(fAllCurr*fAllCurr));  
  
  nucres_.xnucfact = TMath::Max(float(0.),float(fTotalProbCurr*fTotalProbCurr));
  nucres_.xnucelafact = TMath::Max(float(0.),float(fElasticProbCurr*fElasticProbCurr));
  nucres_.xnucspifact = TMath::Max(float(0.),float(fSinglePiProbCurr*fSinglePiProbCurr));
  nucres_.xnucdpifact = TMath::Max(float(0.),float(fDoublePiProbCurr*fDoublePiProbCurr));

  std::cout << "3xnucfact in NReWeightCascNucleon = " << nucres_.xnucfact << std::endl;
  //std::cout << "3xnucelafact in NReWeightCascNucleon = " << nucres_.xnucelafact << std::endl;
  //std::cout << "3xnucspifact in NReWeightCascNucleon = " << nucres_.xnucspifact << std::endl;
  //std::cout << "3xnucdpifact in NReWeightCascNucleon = " << nucres_.xnucdpifact << std::endl;
  std::cout << "fTotalProbDef " << fTotalProbDef << std::endl;

  fortFns->Reconfigure();
  //  nucleonfsihist_.nfreweightnucleonflag = 1;
 
  //nrint_.ptotnuc = 5;
  //nucleonfsihist_.nfnvert = 3876;

  //fortFns->evdifcrs();
  //fortFns->evpiprob();
  nrint_.pcascprob = 1;
  std::cout << "nrint_.pcascprob after reconfigure = " << nrint_.pcascprob << std::endl;
  
  fortFns->nrprton(); //runs nfortfns code here to run nrprton and set value for nfnvert
  //fortFns->nrnuc();
#ifdef _N_REWEIGHT_CASC_DEBUG_
  fortFns->print_allparams();
#endif
  std::cout << "nrint_.pcascprob after nrprton = " << nrint_.pcascprob << std::endl;
  //cout << "NReWeightCascNucleon(): evpiprob = " <<  fortFns->evpiprob() << ", FSIPROB= " << fsihist_.fsiprob << endl;
  float new_xsec  = nrint_.pcascprob;
  
  // evpiprob not passing properly with g77 64-bit so grab from common block
#if defined(f2cFortran)&&!defined(gFortran)
  //if (new_xsec) {
  //cout << "NReWeightCascNucleon() Error: new_xsec=" << new_xsec << " when 0 expected for g77 compiled Fortran code on 64-bit machine" << endl;
  //exit (-1);
  //}
  new_xsec = nrint_.pcascprob;
#endif
  float new_weight = (new_xsec/old_xsec);
  if(old_xsec == 1)
    {


      new_weight =1;
    }
  if(new_xsec!=1 ||old_xsec!=1)
    {

      std::cout << "for this mode, there is a weight not equal to 1 " << std::endl;
      std::cout << "new_xsec = " << new_xsec << "old_xsec = " << old_xsec << std::endl;
    }
  
#ifdef _N_REWEIGHT_CASC_DEBUG_
  cout << "FSI Probability (new) = " << new_xsec << endl;
#endif

  if (isinf(new_weight) || isnan(new_weight) || new_weight<0) {
    cout << "NReWeightCascNucleon::CalcWeight() Warning: new_weight is infinite or negative (" << new_weight << "), setting to 1" << endl;
    //    new_weight = 1;
    new_weight = -9999;
  }


  cout << "new weight = " << new_weight << endl;
  cout << "Neut Mode = " << 

  std::cout << "old prob = " << old_xsec << " new prob = " << new_xsec << std::endl;
  return new_weight;
}
//_______________________________________________________________________________________
double NReWeightCascNucleon::CalcChisq(void)
{
  double chisq = 0.;

  // For 1-sigma change
  //  chisq += TMath::Power(fInelLowTwkDial, 2.);
  //chisq += TMath::Power(fInelHighTwkDial, 2.);
  //chisq += TMath::Power(fAbsTwkDial, 2.);
  //chisq += TMath::Power(fPiProdTwkDial, 2.);
  //chisq += TMath::Power(fCExLowTwkDial, 2.);
  //chisq += TMath::Power(fCExHighTwkDial, 2.);
  //chisq += TMath::Power(fAllTwkDial, 2.);
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
  // This must be in the same order as ${T2KREWEIGHT}/src/T2KSyst.h
  //parVals.push_back(fAbsCurr);   
  //parVals.push_back(fInelLowCurr); 
  //parVals.push_back(fCExLowCurr);  
  //parVals.push_back(fInelHighCurr);
  //parVals.push_back(fCExHighCurr);   
  //parVals.push_back(fPiProdCurr);	   
  //parVals.push_back(fAllCurr);
  parVals.push_back(fTotalProbCurr);
  parVals.push_back(fElasticProbCurr);
  parVals.push_back(fSinglePiProbCurr);
  parVals.push_back(fDoublePiProbCurr); 
  return parVals;  
}
//_______________________________________________________________________________________
