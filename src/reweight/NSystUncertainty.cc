//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Apr 27, 2010 - CA
   Included new parameters in preparation for the Summer 2010 T2K analyses.
   Added option to override the default 1\sigma errors.
 @ Nov 25, 2010 - CA
   Allow for asymmetric 1 sigma fractional errors.
 @ Apr. 6, 2011 - PD
   Implement for NEUT
*/
//____________________________________________________________________________

#include <iostream>

#include "NSystUncertainty.h"

using namespace neut;
using namespace neut::rew;

NSystUncertainty * NSystUncertainty::fInstance = 0;
//____________________________________________________________________________
NSystUncertainty::NSystUncertainty()
{
//  fInstance = 0;
}
//____________________________________________________________________________
NSystUncertainty::~NSystUncertainty()
{
  fInstance = 0;
}
//____________________________________________________________________________
NSystUncertainty * NSystUncertainty::Instance()
{
  if(fInstance == 0) {
    std::cout << std::endl << "NSystUncertainty late initialization" << std::endl;
    static NSystUncertainty::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new NSystUncertainty;
    fInstance->SetDefaults();
  }
  return fInstance;
}
//____________________________________________________________________________
double NSystUncertainty::OneSigmaErr(NSyst_t s, int sign) const
{
  if(sign > 0) {
    map<NSyst_t,double>::const_iterator it = fOneSigPlusErrMap.find(s);
    if(it != fOneSigPlusErrMap.end()) return it->second;
    return 0;
  } 
  else 
  if(sign < 0) {
    map<NSyst_t,double>::const_iterator it = fOneSigMnusErrMap.find(s);
    if(it != fOneSigMnusErrMap.end()) return it->second;
    return 0;
  } 
  else {
    // Handle default argument (sign=0)
    // Case added for compatibility purposes since most existing weight 
    // calcutators call NSystUncertainty::OneSigmaErr(NSyst_t) and the error 
    // on most NSyst_t params is symmetric.
    double err = 0.5 * (
        this->OneSigmaErr(s, +1) + this->OneSigmaErr(s, -1));
    return err;
  }
}
//____________________________________________________________________________
void NSystUncertainty::SetUncertainty(
   NSyst_t s, double plus_err, double minus_err)
{
  //fOneSigPlusErrMap.insert( map<NSyst_t,double>::value_type(s, plus_err ) );
  //fOneSigMnusErrMap.insert( map<NSyst_t,double>::value_type(s, minus_err) );

  fOneSigPlusErrMap[s] = plus_err;
  fOneSigMnusErrMap[s] =  minus_err;

}
//____________________________________________________________________________
void NSystUncertainty::GetUncertainty(
   NSyst_t s, double& plus_err, double& minus_err)
{
 
  plus_err = fOneSigPlusErrMap[s];
  minus_err = fOneSigMnusErrMap[s];

}//____________________________________________________________________________
void NSystUncertainty::SetDefaults(void)
{
  this->SetUncertainty( kXSecTwkDial_NormNCEL,       0.15, 0.15);
  this->SetUncertainty( kXSecTwkDial_MaNCELshape,    0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_1overMaNCEL2,    0.2, 0.2);
  this->SetUncertainty( kXSecTwkDial_MaNCEL,         0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_AxlFFNCEL,    0, 1);
  this->SetUncertainty( kXSecTwkDial_VecFFNCEL,    0, 1);
  //this->SetUncertainty( kXSecTwkDial_EtaNCEL,        0.30, 0.30);
  this->SetUncertainty( kXSecTwkDial_NormCCQE,       0.15, 0.15);
  this->SetUncertainty( kXSecTwkDial_MaCCQEshape,    0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_1overMaCCQE2,    0.2, 0.2);
  this->SetUncertainty( kXSecTwkDial_MaCCQE,         0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_AxlFFCCQE,    0, 1);
  this->SetUncertainty( kXSecTwkDial_VecFFCCQE,    0, 1);
  this->SetUncertainty( kXSecTwkDial_VecFFCCQE_out, 0, 1);
  this->SetUncertainty( kXSecTwkDial_SCCVecQE,    1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_SCCAxlQE,    1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_PsFF,    0.03, 0.03);

  this->SetUncertainty( kXSecTwkDial_NormRES,      0.20, 0.20);
  this->SetUncertainty( kXSecTwkDial_MaRESshape,   0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_MaRES,        0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_MvRES,        0.119047619, 0.119047619);

  this->SetUncertainty( kXSecTwkDial_FFRES,       1,1);
  this->SetUncertainty( kXSecTwkDial_TypeRES,     1,1);
  this->SetUncertainty( kXSecTwkDial_CA5RES,      0.247524752,0.247524752);
  this->SetUncertainty( kXSecTwkDial_BgSclRES,    0.153846154,0.153846154);
  this->SetUncertainty( kXSecTwkDial_MaNFFRES,    0.157894737,0.157894737);
  this->SetUncertainty( kXSecTwkDial_MvNFFRES,    0.119047619, 0.119047619);
  this->SetUncertainty( kXSecTwkDial_MaRSRES,     0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_MvRSRES,     0.119047619, 0.119047619);

  this->SetUncertainty( kXSecTwkDial_NormCCRES,      0.20, 0.20);
  this->SetUncertainty( kXSecTwkDial_MaCCRESshape,   0.165289256, 0.165289256);
  //this->SetUncertainty( kXSecTwkDial_MvCCRESshape,   0.05, 0.05);
  this->SetUncertainty( kXSecTwkDial_MaCCRES,        0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_MvCCRES,        0.119047619, 0.119047619);

  this->SetUncertainty( kXSecTwkDial_FFCCRES,       1,1);
  this->SetUncertainty( kXSecTwkDial_TypeCCRES,     1,1);
  this->SetUncertainty( kXSecTwkDial_CA5CCRES,      0.247524752,0.247524752);
  this->SetUncertainty( kXSecTwkDial_BgSclCCRES,    0.153846154,0.153846154);
  this->SetUncertainty( kXSecTwkDial_MaNFFCCRES,    0.157894737,0.157894737);
  this->SetUncertainty( kXSecTwkDial_MvNFFCCRES,    0.119047619, 0.119047619);
  this->SetUncertainty( kXSecTwkDial_MaRSCCRES,     0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_MvRSCCRES,     0.119047619, 0.119047619);

  this->SetUncertainty( kXSecTwkDial_NormNCRES,      0.20, 0.20);
  this->SetUncertainty( kXSecTwkDial_MaNCRESshape,   0.165289256, 0.165289256);
  //this->SetUncertainty( kXSecTwkDial_MvNCRESshape,   0.05, 0.05);
  this->SetUncertainty( kXSecTwkDial_MaNCRES,        0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_MvNCRES,        0.119047619, 0.119047619);

  this->SetUncertainty( kXSecTwkDial_FFNCRES,       1,1);
  this->SetUncertainty( kXSecTwkDial_TypeNCRES,     1,1);
  this->SetUncertainty( kXSecTwkDial_CA5NCRES,      0.247524752,0.247524752);
  this->SetUncertainty( kXSecTwkDial_BgSclNCRES,    0.153846154,0.153846154);
  this->SetUncertainty( kXSecTwkDial_MaNFFNCRES,    0.157894737,0.157894737);
  this->SetUncertainty( kXSecTwkDial_MvNFFNCRES,    0.119047619, 0.119047619);
  this->SetUncertainty( kXSecTwkDial_MaRSNCRES,     0.165289256, 0.165289256);
  this->SetUncertainty( kXSecTwkDial_MvRSNCRES,     0.119047619, 0.119047619);

  this->SetUncertainty( kXSecTwkDial_NECOHEPI,          0,    1);
  this->SetUncertainty( kXSecTwkDial_MaCOHpi,        0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_R0COHpi,        0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_fA1COHpi,       0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_fb1COHpi,       0.10, 0.10);

  //this->SetUncertainty( kXSecTwkDial_RvpCC1pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvpCC2pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvpNC1pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvpNC2pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvnCC1pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvnCC2pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvnNC1pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvnNC2pi,       0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarpCC1pi,    0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarpCC2pi,    0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarpNC1pi,    0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarpNC2pi,    0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarnCC1pi,    0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarnCC2pi,    0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarnNC1pi,    0.50, 0.50);
  //this->SetUncertainty( kXSecTwkDial_RvbarnNC2pi,    0.50, 0.50);

  // From Debdatta's thesis: 
  //   Aht  = 0.538 +/- 0.134	
  //   Bht  = 0.305 +/- 0.076
  //   CV1u = 0.291 +/- 0.087
  //   CV2u = 0.189 +/- 0.076

  //this->SetUncertainty( kXSecTwkDial_AhtBY,          0.25, 0.25);
  //this->SetUncertainty( kXSecTwkDial_BhtBY,          0.25, 0.25);
  //this->SetUncertainty( kXSecTwkDial_CV1uBY,         0.30, 0.30);
  //this->SetUncertainty( kXSecTwkDial_CV2uBY,         0.40, 0.40);

  //this->SetUncertainty( kXSecTwkDial_AhtBYshape,     0.25, 0.25);
  //this->SetUncertainty( kXSecTwkDial_BhtBYshape,     0.25, 0.25);
  //this->SetUncertainty( kXSecTwkDial_CV1uBYshape,    0.30, 0.30);
  //this->SetUncertainty( kXSecTwkDial_CV2uBYshape,    0.40, 0.40);
  this->SetUncertainty( kXSecTwkDial_NormDIS,          0.30, 0.30);
  this->SetUncertainty( kXSecTwkDial_BYOnOffDIS,          0, 1);
  //this->SetUncertainty( kXSecTwkDial_RnubarnuCC,       0.05, 0.05);
  //this->SetUncertainty( kXSecTwkDial_DISNuclMod,     1.00, 1.00);
  this->SetUncertainty( kXSecTwkDial_NC,               0.30, 0.30);

  this->SetUncertainty( kXSecTwkDial_FAxlCCQEAlpha,  1.0, 1.0 );
  this->SetUncertainty( kXSecTwkDial_FAxlCCQEGamma,  1.0, 1.0 );
  this->SetUncertainty( kXSecTwkDial_FAxlCCQEBeta,   1.0, 1.0 );
  this->SetUncertainty( kXSecTwkDial_FAxlCCQETheta,  1.0, 1.0 );

  //this->SetUncertainty( kHadrAGKYTwkDial_xF1pi,      0.20, 0.20);
  //this->SetUncertainty( kHadrAGKYTwkDial_pT1pi,      0.03, 0.03);
  //this->SetUncertainty( kHadrNuclTwkDial_FormZone,   0.50, 0.50);

  // From NEUT Cascade pi+A and N+A mode comparisons with hadron scattering data:
  //
  this->SetUncertainty( kCascTwkDial_FrAbs_pi,         0.50, 0.50);
  this->SetUncertainty( kCascTwkDial_FrInelLow_pi,     0.50, 0.50);
  this->SetUncertainty( kCascTwkDial_FrCExLow_pi,      0.50, 0.50);
  this->SetUncertainty( kCascTwkDial_FrInelHigh_pi,    0.30, 0.30);
  this->SetUncertainty( kCascTwkDial_FrCExHigh_pi,     0.30, 0.30);
  this->SetUncertainty( kCascTwkDial_FrPiProd_pi,      0.50, 0.50);
  this->SetUncertainty( kCascTwkDial_All_pi,           0.50, 0.50);
  //this->SetUncertainty( kINukeTwkDial_MFP_N,         0.20, 0.20);
  //this->SetUncertainty( kINukeTwkDial_FrCEx_N,       0.50, 0.50);
  //this->SetUncertainty( kINukeTwkDial_FrElas_N,      0.30, 0.30);
  //this->SetUncertainty( kINukeTwkDial_FrInel_N,      0.40, 0.40);
  //this->SetUncertainty( kINukeTwkDial_FrAbs_N,       0.20, 0.20);
  //this->SetUncertainty( kINukeTwkDial_FrPiProd_N,    0.20, 0.20);
  //
  this->SetUncertainty( kSystNucl_CCQEPauliSupViaKF, 0.01, 0.01);
  this->SetUncertainty( kSystNucl_CCQEFermiSurfMom,  0.15, 0.15);
  this->SetUncertainty( kSystNucl_CCQEBindingEnergy, 0.50, 0.50);

  this->SetUncertainty( kXSecTwkDial_FAZExp_TCut      , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_T0        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A0        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A1        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A2        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A3        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A4        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A5        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A6        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A7        , 1.0, 1.0);
  this->SetUncertainty( kXSecTwkDial_FAZExp_A8        , 1.0, 1.0); 
  this->SetUncertainty( kXSecTwkDial_FAZExp_A9        , 1.0, 1.0); 
  
  //this->SetUncertainty( kRDcyTwkDial_BR1gamma,       0.50, 0.50);
  //this->SetUncertainty( kRDcyTwkDial_BR1eta,         0.50, 0.50);
  //this->SetUncertainty( kRDcyTwkDial_Theta_Delta2Npi,         0.50, 0.50);

  this->SetUncertainty( kSystNucl_PilessDcyRES, 0.20, 0.20);

  
}
//____________________________________________________________________________
