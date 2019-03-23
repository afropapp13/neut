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
#include "NReWeightNuXSecNCEL.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_NCEL_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;

const int NReWeightNuXSecNCEL::kModeMa;
const int NReWeightNuXSecNCEL::kModeNormAndMaShape;

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
  fMaDef       = fortFns->XMAQEdef;
  fMaCurr      = fMaDef;

  fAxlFFTwkDial   = 0.; 
  fAxlFFDef       = fortFns->MDLQEAFdef;
  fAxlFFCurr      = fAxlFFDef;

  fVecFFTwkDial   = 0.; 
  fVecFFDef       = fortFns->MDLQEdef;
  fVecFFCurr      = fVecFFDef;

  fPfTwkDial   = 0.; 
  fPfDef       = 0;  // Set in nesetfgparms.F  
  fPfCurr      = fPfDef;

  fEbTwkDial   = 0.; 
  fEbDef       = 0;  // Set in nesetfgparms.F    
  fEbCurr      = fEbDef;

  fKapTwkDial   = 0.; 
  fKapDef       = fortFns->KAPPdef;
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
}
//_______________________________________________________________________________________
bool NReWeightNuXSecNCEL::IsHandled(NSyst_t syst)
{
  bool handle;

  switch(syst) {

    //case ( kXSecTwkDial_NormNCEL    ) :
    //case ( kXSecTwkDial_MaNCELshape ) :
    //	 if(fMode==kModeNormAndMaShape) { 
    //	   handle = true;  
    //	 } else { 
    //	   handle = false; 
    //	 }
    //  break;

  case ( kXSecTwkDial_MaNCELshape ) :  
  case ( kXSecTwkDial_1overMaNCEL2 ) :  
  case ( kXSecTwkDial_NormNCEL    ) :  // Placed here temporarily, in GENIE it's only used with Shape as above
  case ( kXSecTwkDial_MaNCEL ) :
    //if(fMode==kModeMa) { 
      handle = true;  
  //} else { 
  //  handle = false; 
  //}
  break;

  case ( kXSecTwkDial_AxlFFNCEL ) :
  case ( kXSecTwkDial_VecFFNCEL ) :
    //case ( kSystNucl_CCQEPauliSupViaKF ) :  
    //case ( kSystNucl_CCQEFermiSurfMom  ) :  
    //case ( kSystNucl_CCQEBindingEnergy ) :  
    handle = true;  
  break;

  case ( kXSecTwkDial_SCCVecQE ) :
  case ( kXSecTwkDial_SCCAxlQE ) :
  case ( kXSecTwkDial_PsFF ) :
    handle = true;  
  break;

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
  case ( kXSecTwkDial_NormNCEL ) :
    fNormTwkDial = twk_dial;
    break;

    // Following kludge necessary when using MaNCELshape and MaNCEL in same
    // instantiation (i.e. for weight tree generation)
  case ( kXSecTwkDial_MaNCELshape ) :
    if (fMode==kModeNormAndMaShape)
      fMaTwkDial = twk_dial;
  case ( kXSecTwkDial_1overMaNCEL2 ) :
    if (fMode==kMode1overMa2)
      fMaTwkDial = twk_dial;
  case ( kXSecTwkDial_MaNCEL ) :
    if (fMode==kModeMa)
      fMaTwkDial = twk_dial;
    break;

   case ( kXSecTwkDial_AxlFFNCEL ) :
    fAxlFFTwkDial = twk_dial;
    break;
   case ( kXSecTwkDial_VecFFNCEL ) :
    fVecFFTwkDial = twk_dial;
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
  default:
    break;
  }
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Reset(void)
{
  fNormTwkDial = 0.;
  fNormCurr    = fNormDef;

  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;

  fAxlFFTwkDial   = 0.; 
  fAxlFFCurr      = fAxlFFDef;

  fVecFFTwkDial   = 0.; 
  fVecFFCurr      = fVecFFDef;

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

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  if(fMode==kModeMa) {   
    int    sign_matwk = utils::rew::Sign(fMaTwkDial);
    double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaNCEL, sign_matwk);
    fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
  }
  else if(fMode==kModeNormAndMaShape) {
    int    sign_mashtwk = utils::rew::Sign(fMaTwkDial  );
    double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaNCELshape, sign_mashtwk);
    fMaCurr   = fMaDef * (1. + fMaTwkDial   * fracerr_mash);
  }
  else if(fMode==kMode1overMa2) {
    int    sign_mashtwk = utils::rew::Sign(fMaTwkDial  );
    double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_1overMaNCEL2, sign_mashtwk);
    fMaCurr   = fMaDef / sqrt(1. + fMaTwkDial   * fracerr_mash);
  }

  //else
  //if(fMode==kModeNormAndMaShape) { 
  int    sign_normtwk = utils::rew::Sign(fNormTwkDial);
  double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormNCEL,    sign_normtwk);
  fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);     
  //}

  fNormCurr = TMath::Max(0., fNormCurr);
  fMaCurr   = TMath::Max(0., fMaCurr  );

  fAxlFFCurr = TMath::Max(0,fAxlFFTwkDial);
  fVecFFCurr = TMath::Max(0,fVecFFTwkDial);

  //int    sign_kaptwk = utils::rew::Sign(fKapTwkDial);
  //double fracerr_kap = fracerr->OneSigmaErr(kSystNucl_CCQEPauliSupViaKF,    sign_kaptwk);
  //fKapCurr = fKapDef * (1. + fKapTwkDial * fracerr_kap);
  //fKapCurr = TMath::Max(0., fKapCurr);

  //int    sign_pftwk = utils::rew::Sign(fPfTwkDial);
  //double fracerr_pf = fracerr->OneSigmaErr(kSystNucl_CCQEFermiSurfMom,    sign_pftwk);
  // Save fractional change only since default depends on target nucleus of given event
  // Multiplied by default later in CalcWeightMa()
  //fPfCurr = 1. + fPfTwkDial * fracerr_pf;
  //fPfCurr = TMath::Max(0., fPfCurr);

  //int    sign_ebtwk = utils::rew::Sign(fEbTwkDial);
  //double fracerr_eb = fracerr->OneSigmaErr(kSystNucl_CCQEBindingEnergy,    sign_ebtwk);
  // Save fractional change only since default depends on target nucleus of given event
  // Multiplied by default later in CalcWeightMa()
  //fEbCurr = 1. + fEbTwkDial * fracerr_eb;
  //fEbCurr = TMath::Max(0., fEbCurr);

  int sign_sccvec = utils::rew::Sign(fSCCVecTwkDial);
  double fracerr_sccvec = fracerr->OneSigmaErr(kXSecTwkDial_SCCVecQE, sign_sccvec);
  fSCCVecCurr = fSCCVecDef + fSCCVecTwkDial * fracerr_sccvec;//RT: not sure if this is correct just yet....

  int sign_sccaxl = utils::rew::Sign(fSCCAxlTwkDial);
  double fracerr_sccaxl = fracerr->OneSigmaErr(kXSecTwkDial_SCCAxlQE, sign_sccaxl);
  fSCCAxlCurr = fSCCAxlDef + fSCCAxlTwkDial * fracerr_sccaxl;//RT: not sure if this is correct just yet....

  int sign_psff = utils::rew::Sign(fPsFFTwkDial);
  double fracerr_psff = fracerr->OneSigmaErr(kXSecTwkDial_PsFF, sign_psff);
  fPsFFCurr = fPsFFDef * (1. + fPsFFTwkDial * fracerr_psff);//RT: not sure if this is correct just yet....

}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeight() 
{
  bool is_ncel = modeDefn.isNCEL(nework_.modene);
  if(!is_ncel) return 1.;

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
double NReWeightNuXSecNCEL::CalcChisq()
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
double NReWeightNuXSecNCEL::CalcWeightNorm() 
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;

  return wght;

}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeightMa() 
{
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fKapTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fPfTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fEbTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fAxlFFTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fVecFFTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fSCCVecTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fSCCAxlTwkDial) > controls::kASmallNum ||
		  TMath::Abs(fPsFFTwkDial) > controls::kASmallNum
		  );
  
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
#ifdef _N_REWEIGHT_NCEL_DEBUG_
    cout << "NReWeightNuXSecNCEL::CalcWeightMa() Warning: old_xsec==0, setting weight to 1" << '\n';
#endif
    return 1;
  }
  
  nemdls_.xmaqe = fMaCurr;
  nemdls_.kapp  = fKapCurr;

// MDLQE    : CC Quasi-elastic / NC elastic model
//          : x1 : Smith-Moniz for CC
//          : x2 : Smith-Moniz for CC with BBBA05
//          : 0x : Scaling to NCEL     ( same as 5.0.x )
//          : 1x : Scaling to Spectrum func. with Dipole (prior to v5.1.2)
//          : 2x : Scaling to Spectrum func. with BBBA05 (default from v5.1.2)
  if (fAxlFFCurr) nemdls_.mdlqeaf = fAxlFFCurr;
  if (fVecFFCurr) nemdls_.mdlqe = fVecFFCurr;

  if (fSCCVecCurr) nemdls_.sccfv = fSCCVecCurr;
  if (fSCCAxlCurr) nemdls_.sccfa = fSCCAxlCurr;
  if (fPsFFCurr) nemdls_.fpqe = fPsFFCurr;

  fortFns->Reconfigure();

  // Grab FG parameters after they've been set via nesetgfparams depending on target nucleus
  fPfDef = nenupr_.pfsurf;  // Unit conversion is taken care of in qedifcrs.c
  fEbDef = nenupr_.vnuini;  // Unit conversion is taken care of in qedifcrs.c
  nenupr_.pfsurf = fPfDef * fPfCurr;
  nenupr_.vnuini = fEbDef * fEbCurr;

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  fortFns->print_allparams();
#endif
  float new_xsec   = fortFns->evdifcrs();
  float new_weight = (new_xsec/old_xsec);
  //float new_weight = old_weight * (new_xsec/old_xsec);


#ifdef _N_REWEIGHT_NCEL_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << '\n';
  cout << "differential cross section (new) = " << new_xsec << '\n';
#endif 

  // Normalize out change in total cross section from MA variations for MA_shape variation
  if(fMode==kModeNormAndMaShape) {
    
    // Get neutrino type
    int inu=-1;

    if (nework_.ipne[0] > 0 ) inu = neutTotCrs->nue;
    else if (nework_.ipne[0] < 0) inu = neutTotCrs->nueb;
    else {
#ifdef _N_REWEIGHT_NCEL_DEBUG_
      cout << "NReWeightNuXSecNCEL::CalcWeightMa() Error: Cannot MA-shape reweight this neutrino type = " 
	       << nework_.ipne[0] << '\n';
#endif
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
#ifdef _N_REWEIGHT_NCEL_DEBUG_
      cout << "NReWeightNuXSecNCEL::CalcWeightMa() Warning: new_tot_xsec==0, setting weight to 1" << '\n';
#endif
      return 1;
    }

#ifdef _N_REWEIGHT_NCEL_DEBUG_
    cout << "total cross section (old) = " << old_tot_xsec << '\n';
    cout << "total cross section (new) = " << new_tot_xsec << '\n';
#endif 

    new_weight *= old_tot_xsec / new_tot_xsec ;
  }

  if (isinf(new_weight) || isnan(new_weight)) {
#ifdef _N_REWEIGHT_NCEL_DEBUG_
    cout << "NReWeightNuXSecNCEL::CalcWeightMa() Warning: new_weight is infinite, setting to 1" << '\n';
#endif
    return 1;
  }
  
#ifdef _N_REWEIGHT_NCEL_DEBUG_
  cout << "new weight = " << new_weight << '\n';
#endif 

  return new_weight;
  
}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeightMaShape() 
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

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(E,Q2,new_weight);
#endif

  return new_weight;
  */
}
//_______________________________________________________________________________________
