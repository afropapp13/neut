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
 @ Oct 22, 2010 - CA
   Make static consts kModeMa and kModeNormAndMaShape public to aid
   external configuration.
 @ Nov 25, 2010 - CA
   Allow asymmetric errors
 @ Apr 6, 2011 - PD
   Implement for NEUT
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NReWeightNuXSecCCQE.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_CCQE_DEBUG_

using namespace neut;
using namespace rew;

using std::cout;
using std::endl;

const int NReWeightNuXSecCCQE::kModeMa;
const int NReWeightNuXSecCCQE::kModeNormAndMaShape;

//_______________________________________________________________________________________
NReWeightNuXSecCCQE::NReWeightNuXSecCCQE()
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuXSecCCQE::~NReWeightNuXSecCCQE()
{

}
//_______________________________________________________________________________________
void NReWeightNuXSecCCQE::Init(void)
{
  fortFns = NFortFns::Instance();

  neutTotCrs = NTotCrs::Instance();

  this->SetMode(kModeMa);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;

  fMaTwkDial   = 0.;
  fMaDef       = 0.;
  fMaCurr      = fMaDef;

  fAxlFFTwkDial   = 0.;
  fAxlFFDef       = 0.;
  fAxlFFCurr      = fAxlFFDef;

  fVecFFDef       = 0.;
  fVecFFCurr      = fVecFFDef;
  fVecFFTwkDial   = 0.;

  fVecFFOutDef     = 0.;
  fVecFFOutCurr    = fVecFFOutDef;
  fVecFFOutTwkDial = 0.;

  fPfTwkDial   = 0.;
  fPfDef       = 0;  // Set in nesetfgparms.F
  fPfCurr      = fPfDef;

  fEbTwkDial   = 0.;
  fEbDef       = 0;  // Set in nesetfgparms.F
  fEbCurr      = fEbDef;

  fKapTwkDial   = 0.;
  fKapDef       = 0.;
  fKapCurr      = fKapDef;

  fSCCVecTwkDial  = 0;
  fSCCVecDef      = fortFns->SCCFVdef;
  fSCCVecCurr     = fSCCVecDef;

  fSCCAxlTwkDial  = 0;
  fSCCAxlDef      = fortFns->SCCFAdef;
  fSCCAxlCurr     = fSCCAxlDef;

  fPsFFTwkDial  = 0;
  fPsFFDef      = fortFns->FPQEdef;
  fPsFFCurr     = fPsFFDef;

  // Alternative Form Factors
  fAltMDLQEAF = 1; // Default - Dipole

  // 2/3 Comp
  fTwk_3CompAlpha  = 0.0;
  fDef_3CompAlpha  = 0.95;
  fCurr_3CompAlpha = fDef_3CompAlpha;
  fErr_3CompAlpha  = 0.0;

  fTwk_3CompGamma  = 0.0;
  fDef_3CompGamma  = 0.515;
  fCurr_3CompGamma = fDef_3CompGamma;
  fErr_3CompGamma  = 0.0;

  fTwk_3CompBeta  = 0.0;
  fDef_3CompBeta  = 2.0;
  fCurr_3CompBeta = fDef_3CompBeta;
  fErr_3CompBeta  = 0.0;

  fTwk_3CompTheta  = 0.0;
  fDef_3CompTheta  = -0.15;
  fCurr_3CompTheta = fDef_3CompTheta;
  fErr_3CompTheta  = 0.0;

  // Z expansion
  kMaxZExpA = 10;

  fDef_ZExp_NTerms  = 4;
  fCurr_ZExp_NTerms = fDef_ZExp_NTerms;
  fDef_ZExp_Q4Cut   = 1;
  fCurr_ZExp_Q4Cut  = fDef_ZExp_Q4Cut;

  fTwk_ZExp_TCut  = 0.0;
  fDef_ZExp_TCut  = 0.1764;
  fCurr_ZExp_TCut = fDef_ZExp_TCut;
  fErr_ZExp_TCut  = 0.0;

  fTwk_ZExp_T0  = 0.0;
  fDef_ZExp_T0  = -0.28;
  fCurr_ZExp_T0 = fDef_ZExp_T0;
  fErr_ZExp_T0  = 0.0;

  for (int i = 0; i < kMaxZExpA; i++){
    fTwk_ZExp_ATerms[i]  = 0.0;
    fDef_ZExp_ATerms[i]  = 1.0;
    fCurr_ZExp_ATerms[i] = 1.0;
    fErr_ZExp_ATerms[i]  = 0.0;
  }

  // Set 4 Starting  ZExp Terms
  fDef_ZExp_ATerms[1] =  2.3;
  fDef_ZExp_ATerms[2] = -0.6;
  fDef_ZExp_ATerms[3] = -3.8;
  fDef_ZExp_ATerms[4] =  2.3;
}
//_______________________________________________________________________________________
bool NReWeightNuXSecCCQE::IsHandled(NSyst_t syst)
{
  bool handle;

  switch(syst) {

    //case ( kXSecTwkDial_NormCCQE    ) :
    //case ( kXSecTwkDial_MaCCQEshape ) :
    //	 if(fMode==kModeNormAndMaShape) {
    //	   handle = true;
    //	 } else {
    //	   handle = false;
    //	 }
    //  break;

  case ( kXSecTwkDial_MaCCQEshape ) :
  case ( kXSecTwkDial_1overMaCCQE2 ) :
  case ( kXSecTwkDial_NormCCQE    ) :  // Placed here temporarily, in GENIE it's only used with Shape as above
  case ( kXSecTwkDial_MaCCQE ) :
    //if(fMode==kModeMa) {
      handle = true;
  //} else {
  //  handle = false;
  //}
  break;

  case ( kXSecTwkDial_AxlFFCCQE ) :
  case ( kXSecTwkDial_VecFFCCQE ) :
  case ( kXSecTwkDial_VecFFCCQE_out ) :
  case ( kSystNucl_CCQEPauliSupViaKF ) :
  case ( kSystNucl_CCQEFermiSurfMom  ) :
  case ( kSystNucl_CCQEBindingEnergy ) :
    handle = true;
  break;

  case ( kXSecTwkDial_SCCVecQE ) :
  case ( kXSecTwkDial_SCCAxlQE ) :
  case ( kXSecTwkDial_PsFF ) :
    handle = true;

  case ( kXSecTwkDial_AxlDipToAlt ) :
    handle = true;
  break;

  // 2/3 Comp
  case ( kXSecTwkDial_FAxlCCQEAlpha ) :
  case ( kXSecTwkDial_FAxlCCQEGamma ) :
  case ( kXSecTwkDial_FAxlCCQEBeta  ) :
  case ( kXSecTwkDial_FAxlCCQETheta ) :
    handle = true;
    break;

    // ZExp
  case (kXSecTwkDial_FAZExp_NTerms) :
  case (kXSecTwkDial_FAZExp_TCut) :
  case (kXSecTwkDial_FAZExp_T0 ) :
  case (kXSecTwkDial_FAZExp_Q4Cut ) :
  case (kXSecTwkDial_FAZExp_A0 ) :
  case (kXSecTwkDial_FAZExp_A1 ) :
  case (kXSecTwkDial_FAZExp_A2 ) :
  case (kXSecTwkDial_FAZExp_A3 ) :
  case (kXSecTwkDial_FAZExp_A4 ) :
  case (kXSecTwkDial_FAZExp_A5 ) :
  case (kXSecTwkDial_FAZExp_A6 ) :
  case (kXSecTwkDial_FAZExp_A7 ) :
  case (kXSecTwkDial_FAZExp_A8 ) :
  case (kXSecTwkDial_FAZExp_A9 ) :
    handle = true;
  break;

  default:
    handle = false;
    break;
  }

  return handle;
}
//_______________________________________________________________________________________
void NReWeightNuXSecCCQE::SetSystematic(NSyst_t syst, double twk_dial)
{
  switch(syst) {
  case ( kXSecTwkDial_NormCCQE ) :
    fNormTwkDial = twk_dial;
    break;

    // Following kludge necessary when using MaCCQEshape and MaCCQE in same
    // instantiation (i.e. for weight tree generation)
  case ( kXSecTwkDial_MaCCQEshape ) :
    if (fMode==kModeNormAndMaShape)
      fMaTwkDial = twk_dial;
  case ( kXSecTwkDial_1overMaCCQE2 ) :
    if (fMode==kMode1overMa2)
      fMaTwkDial = twk_dial;
  case ( kXSecTwkDial_MaCCQE ) :
    if (fMode==kModeMa)
      fMaTwkDial = twk_dial;
    break;

   case ( kXSecTwkDial_AxlFFCCQE ) :
    fAxlFFTwkDial = twk_dial;
    break;
   case ( kXSecTwkDial_VecFFCCQE ) :
    fVecFFTwkDial = twk_dial;
    break;
  case ( kXSecTwkDial_VecFFCCQE_out ) :
    fVecFFOutTwkDial = twk_dial;
    break;
  case ( kSystNucl_CCQEPauliSupViaKF ) :
    fKapTwkDial = twk_dial;
    break;
  case ( kSystNucl_CCQEFermiSurfMom ) :
    fPfTwkDial = twk_dial;
    break;
  case ( kSystNucl_CCQEBindingEnergy ) :
    fEbTwkDial = twk_dial;
    break;

  case ( kXSecTwkDial_SCCVecQE ) :
    fSCCVecTwkDial = twk_dial;
    break;
  case ( kXSecTwkDial_SCCAxlQE ) :
    fSCCAxlTwkDial = twk_dial;
    break;
  case ( kXSecTwkDial_PsFF ) :
    fPsFFTwkDial = twk_dial;
    break;

  case ( kXSecTwkDial_AxlDipToAlt ) :
    fAltMDLQEAF = twk_dial;
    break;

  case ( kXSecTwkDial_FAxlCCQEAlpha ) :
    fTwk_3CompAlpha = twk_dial;
    break;
  case ( kXSecTwkDial_FAxlCCQEGamma ) :
    fTwk_3CompGamma = twk_dial;
    break;
  case ( kXSecTwkDial_FAxlCCQETheta ) :
    fTwk_3CompTheta = twk_dial;
    break;
  case ( kXSecTwkDial_FAxlCCQEBeta ) :
    fTwk_3CompBeta = twk_dial;
    break;


  // ZExp Info Fill IN
  case (kXSecTwkDial_FAZExp_NTerms) :
    fCurr_ZExp_NTerms = twk_dial;
    break;
  case (kXSecTwkDial_FAZExp_Q4Cut ) :
    fCurr_ZExp_Q4Cut = twk_dial;
    break;
  case (kXSecTwkDial_FAZExp_TCut) :
    fTwk_ZExp_TCut = twk_dial;
    break;
  case (kXSecTwkDial_FAZExp_T0 ) :
    fTwk_ZExp_T0 = twk_dial;
    break;

  default:
    break;
  }

  // Fill 10 ZExp ATerms in loop
  int zexpenum = (int) kXSecTwkDial_FAZExp_A0;
  for (int i = 0; i < kMaxZExpA; i++){
    if ( syst == zexpenum + i ){
      fTwk_ZExp_ATerms[i] = twk_dial;
      break;
    }
  }

}
//_______________________________________________________________________________________
void NReWeightNuXSecCCQE::Reset(void)
{
  fNormTwkDial = 0.;
  fNormCurr    = fNormDef;

  fMaTwkDial   = 0.;
  fMaCurr      = fMaDef;

  fAxlFFTwkDial   = 0.;
  fAxlFFCurr      = fAxlFFDef;

  fVecFFTwkDial   = 0.;
  fVecFFCurr      = fVecFFDef;

  fVecFFOutTwkDial = 0.;
  fVecFFOutCurr    = fVecFFOutDef;

  fPfTwkDial   = 0.;
  fPfCurr      = fPfDef;

  fEbTwkDial   = 0.;
  fEbCurr      = fEbDef;

  fKapTwkDial   = 0.;
  fKapCurr      = fKapDef;

  fSCCVecTwkDial = 0.;
  fSCCVecCurr = fSCCVecDef;

  fSCCAxlTwkDial = 0.;
  fSCCAxlCurr = fSCCAxlDef;

  fPsFFTwkDial = 0.;
  fPsFFCurr = fPsFFDef;

  // Alt FF
  fAltMDLQEAF = 1;

  fTwk_3CompAlpha  = 0.0;
  fCurr_3CompAlpha = fDef_3CompAlpha;

  fTwk_3CompGamma  = 0.0;
  fCurr_3CompGamma = fDef_3CompGamma;

  fTwk_3CompBeta  = 0.0;
  fCurr_3CompBeta = fDef_3CompBeta;

  fTwk_3CompTheta  = 0.0;
  fCurr_3CompTheta = fDef_3CompTheta;

  // ZEXP
  // Z expansion
  fCurr_ZExp_NTerms = fDef_ZExp_NTerms;
  fCurr_ZExp_Q4Cut  = fDef_ZExp_Q4Cut;

  fTwk_ZExp_TCut  = 0.0;
  fCurr_ZExp_TCut = fDef_ZExp_TCut;

  fTwk_ZExp_T0  = 0.0;
  fCurr_ZExp_T0 = fDef_ZExp_T0;

  for (int i = 0; i < kMaxZExpA; i++){
    fTwk_ZExp_ATerms[i]  = 0.0;
    fCurr_ZExp_ATerms[i] = fDef_ZExp_ATerms[i];
  }

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuXSecCCQE::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  if(fMode==kModeMa) {
    int    sign_matwk = utils::rew::Sign(fMaTwkDial);
    double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQE, sign_matwk);
    //fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
    fMaCurr = (1. + fMaTwkDial * fracerr_ma);
  }
  else if(fMode==kModeNormAndMaShape) {
    int    sign_mashtwk = utils::rew::Sign(fMaTwkDial  );
    double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaCCQEshape, sign_mashtwk);
    fMaCurr   = (1. + fMaTwkDial   * fracerr_mash);
  }
  else if(fMode==kMode1overMa2) {
    int    sign_mashtwk = utils::rew::Sign(fMaTwkDial  );
    double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_1overMaCCQE2, sign_mashtwk);
    fMaCurr   = 1. / sqrt(1. + fMaTwkDial   * fracerr_mash);
  }

  //else
  //if(fMode==kModeNormAndMaShape) {
  int    sign_normtwk = utils::rew::Sign(fNormTwkDial);
  double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormCCQE,    sign_normtwk);
  fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
  //}

  fNormCurr = TMath::Max(0., fNormCurr);
  fMaCurr   = TMath::Max(0., fMaCurr  );

  int    sign_affshtwk = utils::rew::Sign(fAxlFFTwkDial  );
  double fracerr_affsh = fracerr->OneSigmaErr(kXSecTwkDial_AxlFFCCQE, sign_affshtwk);
  fAxlFFCurr   = (1. + fAxlFFTwkDial   * fracerr_affsh);

  int    sign_vffshtwk = utils::rew::Sign(fVecFFTwkDial  );
  double fracerr_vffsh = fracerr->OneSigmaErr(kXSecTwkDial_VecFFCCQE, sign_vffshtwk);
  fVecFFCurr   = (1. + fVecFFTwkDial   * fracerr_vffsh);

  int    sign_vffoshtwk = utils::rew::Sign(fVecFFOutTwkDial  );
  double fracerr_vffosh = fracerr->OneSigmaErr(kXSecTwkDial_VecFFCCQE_out, sign_vffoshtwk);
  fVecFFOutCurr   = (1. + fVecFFOutTwkDial   * fracerr_vffosh);

  /*
  if (fVecFFTwkDial) fVecFFCurr = fVecFFTwkDial;
  if (fVecFFCurr != fVecFFDef) fortFns->SetMCDefaultVal(kXSecTwkDial_VecFFCCQE, fVecFFCurr);

  if (fVecFFOutTwkDial) fVecFFOutCurr = fVecFFOutTwkDial;
  else fVecFFOutCurr = fVecFFCurr;
  */

  int    sign_kaptwk = utils::rew::Sign(fKapTwkDial);
  double fracerr_kap = fracerr->OneSigmaErr(kSystNucl_CCQEPauliSupViaKF,    sign_kaptwk);
  fKapCurr = (1. + fKapTwkDial * fracerr_kap);
  fKapCurr = TMath::Max(0., fKapCurr);

  int    sign_pftwk = utils::rew::Sign(fPfTwkDial);
  double fracerr_pf = fracerr->OneSigmaErr(kSystNucl_CCQEFermiSurfMom,    sign_pftwk);
  // Save fractional change only since default depends on target nucleus of given event
  // Multiplied by default later in CalcWeightMa()
  fPfCurr = 1. + fPfTwkDial * fracerr_pf;
  fPfCurr = TMath::Max(0., fPfCurr);

  int    sign_ebtwk = utils::rew::Sign(fEbTwkDial);
  double fracerr_eb = fracerr->OneSigmaErr(kSystNucl_CCQEBindingEnergy,    sign_ebtwk);
  // Save fractional change only since default depends on target nucleus of given event
  // Multiplied by default later in CalcWeightMa()
  fEbCurr = 1. + fEbTwkDial * fracerr_eb;
  fEbCurr = TMath::Max(0., fEbCurr);

  int sign_sccvec = utils::rew::Sign(fSCCVecTwkDial);
  double fracerr_sccvec = fracerr->OneSigmaErr(kXSecTwkDial_SCCVecQE, sign_sccvec);
  fSCCVecCurr = fSCCVecDef + fSCCVecTwkDial * fracerr_sccvec;

  int sign_sccaxl = utils::rew::Sign(fSCCAxlTwkDial);
  double fracerr_sccaxl = fracerr->OneSigmaErr(kXSecTwkDial_SCCAxlQE, sign_sccaxl);
  fSCCAxlCurr = fSCCAxlDef + fSCCAxlTwkDial * fracerr_sccaxl;

  int sign_psff = utils::rew::Sign(fPsFFTwkDial);
  double fracerr_psff = fracerr->OneSigmaErr(kXSecTwkDial_PsFF, sign_psff);
  fPsFFCurr = fPsFFDef * (1. + fPsFFTwkDial * fracerr_psff);

  // Make sure MDLQEAF is set correctly
  if (fAltMDLQEAF == 0) fAltMDLQEAF = 1;
  if (fAltMDLQEAF > 5) {
    cout << "ERROR: MDLQEAF can only go between 1 and 5" << endl;
    cout << "1. Dipole, 2. BBBA07, 3. 2 Comp, 4. 3 Comp. 5. Z Exp." << endl;
    throw;
  }

  // Set 2/3 Comp Values
  int sign_3compalpha = utils::rew::Sign(fTwk_3CompAlpha);
  fErr_3CompAlpha = fracerr->OneSigmaErr(kXSecTwkDial_FAxlCCQEAlpha, sign_3compalpha);
  fCurr_3CompAlpha = fDef_3CompAlpha * ( 1. + fTwk_3CompAlpha * fErr_3CompAlpha );

  int sign_3compgamma = utils::rew::Sign(fTwk_3CompGamma);
  fErr_3CompGamma = fracerr->OneSigmaErr(kXSecTwkDial_FAxlCCQEGamma, sign_3compgamma);
  fCurr_3CompGamma = fDef_3CompGamma * ( 1. + fTwk_3CompGamma * fErr_3CompGamma );

  int sign_3compbeta = utils::rew::Sign(fTwk_3CompBeta);
  fErr_3CompBeta = fracerr->OneSigmaErr(kXSecTwkDial_FAxlCCQEBeta, sign_3compbeta);
  fCurr_3CompBeta = fDef_3CompBeta * ( 1. + fTwk_3CompBeta * fErr_3CompBeta );

  int sign_3comptheta = utils::rew::Sign(fTwk_3CompTheta);
  fErr_3CompTheta = fracerr->OneSigmaErr(kXSecTwkDial_FAxlCCQETheta, sign_3comptheta);
  fCurr_3CompTheta = fDef_3CompTheta * ( 1. + fTwk_3CompTheta * fErr_3CompTheta );

  // Set Z Expansion Values
  // TCut
  fErr_ZExp_TCut  = fracerr->OneSigmaErr( kXSecTwkDial_FAZExp_TCut,
					  utils::rew::Sign(fTwk_ZExp_TCut) );
  fCurr_ZExp_TCut = fDef_ZExp_TCut * ( 1.0 + fTwk_ZExp_TCut * fErr_ZExp_TCut );

  // T0
  fErr_ZExp_T0  = fracerr->OneSigmaErr( kXSecTwkDial_FAZExp_T0,
					  utils::rew::Sign(fTwk_ZExp_T0) );
  fCurr_ZExp_T0 = fDef_ZExp_T0 * ( 1.0 + fTwk_ZExp_T0 * fErr_ZExp_T0 );

  // A Terms
  int zexpenum = (int) kXSecTwkDial_FAZExp_A0;
  for (int i = 0; i < kMaxZExpA; i++){
    fErr_ZExp_ATerms[i] = fracerr->OneSigmaErr( NSyst_t(zexpenum + i),
						utils::rew::Sign(fTwk_ZExp_ATerms[i]) );
    fCurr_ZExp_ATerms[i]
      = fDef_ZExp_ATerms[i] * ( 1.0 + fTwk_ZExp_ATerms[i] * fErr_ZExp_ATerms[i] );
  }


}
//_______________________________________________________________________________________
double NReWeightNuXSecCCQE::CalcWeight()
{
  bool is_ccqe = modeDefn.isCCQE(nework_.modene);
  if(!is_ccqe) return 1.;

  int nupdg = nework_.ipne[0];
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  //if(fMode==kModeMa) {
  double wght = this->CalcWeightMa() * this->CalcWeightNorm();
  return wght;
     //}
  //else
  //if(fMode==kModeNormAndMaShape) {
  //   double wght =
  //       this->CalcWeightNorm   () *
  //       this->CalcWeightMaShape();
  //   return wght;
  //}

  return 1.;
}
//_______________________________________________________________________________________
double NReWeightNuXSecCCQE::CalcChisq()
{
  double chisq = 0.;
  //if(fMode==kModeMa) {
     chisq += TMath::Power(fMaTwkDial, 2.);
     //}
  //else
  //if(fMode==kModeNormAndMaShape) {
     chisq += TMath::Power(fNormTwkDial, 2.);
  //   chisq += TMath::Power(fMaTwkDial,   2.);
  //}
  return chisq;
}
//_______________________________________________________________________________________
double NReWeightNuXSecCCQE::CalcWeightNorm()
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;

  return wght;

}
//_______________________________________________________________________________________
double NReWeightNuXSecCCQE::CalcWeightMa()
{
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fKapTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fPfTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fEbTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fAxlFFTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fVecFFTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fVecFFOutTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fSCCVecTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fSCCAxlTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fPsFFTwkDial) > controls::kASmallNum ||
		  (fAltMDLQEAF > 1)
		  );

  if(!tweaked) return 1.0;


  fortFns->SetDefaults();
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_CCQE_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif

  // Double Check Dipole -> Alternative is set.
  if (fAltMDLQEAF > 1) nemdls_.mdlqeaf = 1;

  float old_xsec     = fortFns->evdifcrs();
  //float old_xsec   = event.DiffXSec();
  //float old_weight = event.Weight();

  // Correct for the units change between RFG and SF
  if (nemdls_.mdlqe == 402) old_xsec *= 1E-23;

  if (old_xsec==0) {
    cout << "NReWeightNuXSecCCQE::CalcWeightMa() Warning: old_xsec==0, setting weight to 1" << endl;
    return 1;
  }

  //nemdls_.xmaqe = fMaCurr;//RT: comment out since now common block
  //nemdls_.kapp  = fKapCurr;//should be passing this info through
  fMaDef = nemdls_.xmaqe;
  nemdls_.xmaqe = fMaDef * fMaCurr;

  fKapDef = nemdls_.kapp;
  nemdls_.kapp = fKapDef * fKapCurr;

  // FUNCTION SHOULD
  // - Take dipole cross-section to start
  // - Then calculate weights for other models.

// MDLQE    : CC Quasi-elastic / NC elastic model
//          : x1 : Smith-Moniz for CC
//          : x2 : Smith-Moniz for CC with BBBA05
//          : 0x : Scaling to CCQE     ( same as 5.0.x )
//          : 1x : Scaling to Spectrum func. with Dipole (prior to v5.1.2)
//          : 2x : Scaling to Spectrum func. with BBBA05 (default from v5.1.2)
/*
  if (fAxlFFCurr) nemdls_.mdlqeaf = fAxlFFCurr;
*/
  fVecFFDef = nemdls_.mdlqe;
  //fVecFFCurr = fVecFFDef;
  fAxlFFDef = nemdls_.mdlqeaf;
  //fAxlFFCurr = TMath::Max(0.,fAxlFFTwkDial);

  // If this dial has been tweaked, if changes the default in NFortFns.cc
  nemdls_.mdlqeaf = fAxlFFDef * fAxlFFCurr;
  nemdls_.mdlqe = fVecFFDef * fVecFFCurr;

  //if (fVecFFOutCurr != fVecFFCurr) nemdls_.mdlqe = fVecFFOutCurr;

  if (fSCCVecCurr) nemdls_.sccfv = fSCCVecCurr;
  if (fSCCAxlCurr) nemdls_.sccfa = fSCCAxlCurr;
  if (fPsFFCurr) nemdls_.fpqe = fPsFFCurr;

  // Dipole to alternative RW Functions
  if (fAltMDLQEAF > 1) {
    //    std::cout << "Changing MDLQEAF " << fAltMDLQEAF << std::endl;
    nemdls_.mdlqeaf = fAltMDLQEAF;
  }

  // 2 Component
  if (fAltMDLQEAF == 3){
    nemdls_.axffalpha = fCurr_3CompAlpha;
    nemdls_.axffgamma = fCurr_3CompGamma;
  }

  // 3 Component
  if (fAltMDLQEAF == 4){
    nemdls_.axffalpha = fCurr_3CompAlpha;
    nemdls_.axffgamma = fCurr_3CompGamma;
    nemdls_.axfftheta = fCurr_3CompTheta;
    nemdls_.axffbeta  = fCurr_3CompBeta;
  }

  // Z Expansion
  // Calling zexpconfig handled in NFortFns
  if (fAltMDLQEAF == 5){

    nemdls_.axzexpnt = fCurr_ZExp_NTerms;
    nemdls_.axzexpq4 = fCurr_ZExp_Q4Cut;

    nemdls_.axzexptc = fCurr_ZExp_TCut;
    nemdls_.axzexpt0 = fCurr_ZExp_T0;

    nemdls_.axzexpa0 = fCurr_ZExp_ATerms[0];
    nemdls_.axzexpa1 = fCurr_ZExp_ATerms[1];
    nemdls_.axzexpa2 = fCurr_ZExp_ATerms[2];
    nemdls_.axzexpa3 = fCurr_ZExp_ATerms[3];
    nemdls_.axzexpa4 = fCurr_ZExp_ATerms[4];
    nemdls_.axzexpa5 = fCurr_ZExp_ATerms[5];
    nemdls_.axzexpa6 = fCurr_ZExp_ATerms[6];
    nemdls_.axzexpa7 = fCurr_ZExp_ATerms[7];
    nemdls_.axzexpa8 = fCurr_ZExp_ATerms[8];
    nemdls_.axzexpa9 = fCurr_ZExp_ATerms[9];

  }

  fortFns->Reconfigure();

  // Grab FG parameters after they've been set via nesetgfparams depending on target nucleus
  fPfDef = nenupr_.pfsurf;  // Unit conversion is taken care of in qedifcrs.c
  fEbDef = nenupr_.vnuini;  // Unit conversion is taken care of in qedifcrs.c
  nenupr_.pfsurf = fPfDef * fPfCurr;
  nenupr_.vnuini = fEbDef * fEbCurr;

#ifdef _N_REWEIGHT_CCQE_DEBUG_
  fortFns->print_allparams();
#endif
  float new_xsec   = fortFns->evdifcrs();

  // Correct for the units change between RFG and SF
  if (nemdls_.mdlqe == 402) new_xsec *= 1E-23;

  float new_weight = (new_xsec/old_xsec);
  //float new_weight = old_weight * (new_xsec/old_xsec);

#ifdef _N_REWEIGHT_CCQE_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << endl;
  cout << "differential cross section (new) = " << new_xsec << endl;
#endif

  // Normalize out change in total cross section from MA variations for MA_shape variation
  if(fMode==kModeNormAndMaShape) {

    // Get neutrino type
    int inu=-1;
    if (nework_.ipne[0] == 14) inu = neutTotCrs->numu;
    else if (nework_.ipne[0] == -14) inu = neutTotCrs->numub;
    else if (nework_.ipne[0] == 12) inu = neutTotCrs->nue;
    else if (nework_.ipne[0] == -12) inu = neutTotCrs->nueb;
    else {
      cout << "NReWeightNuXSecCCQE::CalcWeightMa() Error: Cannot MA-shape reweight this neutrino type = "
	   << nework_.ipne[0] << endl;
      exit (-1);
    }

    // Calculate neutrino energy
    float Enu = 0;  // GeV
    for (int i=0; i<3; i++)
      Enu += nework_.pne[0][i]*nework_.pne[0][i];
    Enu = sqrt(Enu);

    // Warning: Does not take into account variations of PFSurf
    float old_tot_xsec = neutTotCrs->ccqe_crs[inu]->Interpolate(Enu, fPfDef, fMaDef);
    float new_tot_xsec = neutTotCrs->ccqe_crs[inu]->Interpolate(Enu, fPfDef, fMaCurr);

    if (new_tot_xsec==0) {
      cout << "NReWeightNuXSecCCQE::CalcWeightMa() Warning: new_tot_xsec==0, setting weight to 1" << endl;
      return 1;
    }

#ifdef _N_REWEIGHT_CCQE_DEBUG_
    cout << "total cross section (old) = " << old_tot_xsec << endl;
    cout << "total cross section (new) = " << new_tot_xsec << endl;
#endif

    new_weight *= old_tot_xsec / new_tot_xsec ;
  }

  if (isinf(new_weight) || isnan(new_weight)) {
    cout << "NReWeightNuXSecCCQE::CalcWeightMa() Warning: new_weight is infinite, setting to 1" << endl;
    return 1;
  }

#ifdef _N_REWEIGHT_CCQE_DEBUG_
  cout << "new weight = " << new_weight << endl;
#endif

  return new_weight;

}
//_______________________________________________________________________________________
double NReWeightNuXSecCCQE::CalcWeightMaShape()
{
  return 1.0;

  /* //GENIE Implementation below
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();
  interaction->SetBit(kIAssumeFreeNucleon);

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

//double old_integrated_xsec = event.XSec();
  double old_integrated_xsec = fXSecModelDef -> Integral(interaction);
  double new_integrated_xsec = fXSecModel    -> Integral(interaction);
  assert(new_integrated_xsec > 0);
  new_weight *= (old_integrated_xsec/new_integrated_xsec);

//LOG("ReW", pDEBUG) << "integrated cross section (old) = " << old_integrated_xsec;
//LOG("ReW", pDEBUG) << "integrated cross section (new) = " << new_integrated_xsec;
//LOG("ReW", pDEBUG) << "new weight (normalized to const integral) = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();
  interaction->ResetBit(kIAssumeFreeNucleon);

#ifdef _N_REWEIGHT_CCQE_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(E,Q2,new_weight);
#endif

  return new_weight;
  */
}
//_______________________________________________________________________________________
