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

#include "NReWeightControls.h"
#include "NReWeightNuXSecNCEL.h"
#include "NReWeightUtils.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "PDGCodes.h"

//#define _N_REWEIGHT_NCEL_DEBUG_

using namespace neut;
using namespace rew;

using std::cout;
using std::endl;

const int NReWeightNuXSecNCEL::kModeMa;
const int NReWeightNuXSecNCEL::kModeNormAndMaShape;

//_______________________________________________________________________________________
NReWeightNuXSecNCEL::NReWeightNuXSecNCEL() { this->Init(); }
//_______________________________________________________________________________________
NReWeightNuXSecNCEL::~NReWeightNuXSecNCEL() {}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Init(void) {
  fortFns = NFortFns::Instance();

  neutTotCrs = NTotCrs::Instance();

  this->SetMode(kModeMa);

  this->RewNue(true);
  this->RewNuebar(true);
  this->RewNumu(true);
  this->RewNumubar(true);

  fNormTwkDial = 0.;
  fNormDef = 1.;
  fNormCurr = fNormDef;

  fMaTwkDial = 0.;
  fMaDef = 0.;
  fMaCurr = fMaDef;

  fAxlFFTwkDial = 0.;
  fAxlFFDef = 0.;
  fAxlFFCurr = fAxlFFDef;

  fVecFFDef = 0.;
  fVecFFCurr = fVecFFDef;
  fVecFFTwkDial = 0.;
}
//_______________________________________________________________________________________
bool NReWeightNuXSecNCEL::IsHandled(NSyst_t syst) {
  bool handle;

  switch (syst) {

    // case ( kXSecTwkDial_NormNCEL    ) :
    // case ( kXSecTwkDial_MaNCELshape ) :
    //	 if(fMode==kModeNormAndMaShape) {
    //	   handle = true;
    //	 } else {
    //	   handle = false;
    //	 }
    //  break;

  case (kXSecTwkDial_MaNCELshape):
  case (kXSecTwkDial_1overMaNCEL2):
  case (kXSecTwkDial_NormNCEL): // Placed here temporarily, in GENIE it's only
                                // used with Shape as above
  case (kXSecTwkDial_MaNCEL):
    // if(fMode==kModeMa) {
    handle = true;
    //} else {
    //  handle = false;
    //}
    break;

  case (kXSecTwkDial_AxlFFNCEL):
  case (kXSecTwkDial_VecFFNCEL):
    handle = true;
    break;

  default:
    handle = false;
    break;
  }

  return handle;
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::SetSystematic(NSyst_t syst, double twk_dial) {
  switch (syst) {
  case (kXSecTwkDial_NormNCEL):
    fNormTwkDial = twk_dial;
    break;

    // Following kludge necessary when using MaNCELshape and MaNCEL in same
    // instantiation (i.e. for weight tree generation)
  case (kXSecTwkDial_MaNCELshape):
    if (fMode == kModeNormAndMaShape)
      fMaTwkDial = twk_dial;
  case (kXSecTwkDial_1overMaNCEL2):
    if (fMode == kMode1overMa2)
      fMaTwkDial = twk_dial;
  case (kXSecTwkDial_MaNCEL):
    if (fMode == kModeMa)
      fMaTwkDial = twk_dial;
    break;

  case (kXSecTwkDial_AxlFFNCEL):
    fAxlFFTwkDial = twk_dial;
    break;
  case (kXSecTwkDial_VecFFNCEL):
    fVecFFTwkDial = twk_dial;
    break;

  default:
    break;
  }
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Reset(void) {
  fNormTwkDial = 0.;
  fNormCurr = fNormDef;

  fMaTwkDial = 0.;
  fMaCurr = fMaDef;

  fAxlFFTwkDial = 0.;
  fAxlFFCurr = fAxlFFDef;

  fVecFFTwkDial = 0.;
  fVecFFCurr = fVecFFDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuXSecNCEL::Reconfigure(void) {
  NSystUncertainty *fracerr = NSystUncertainty::Instance();

  if (fMode == kModeMa) {
    int sign_matwk = utils::rew::Sign(fMaTwkDial);
    double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaNCEL, sign_matwk);
    // fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
    fMaCurr = (1. + fMaTwkDial * fracerr_ma);
  } else if (fMode == kModeNormAndMaShape) {
    int sign_mashtwk = utils::rew::Sign(fMaTwkDial);
    double fracerr_mash =
        fracerr->OneSigmaErr(kXSecTwkDial_MaNCELshape, sign_mashtwk);
    fMaCurr = (1. + fMaTwkDial * fracerr_mash);
  } else if (fMode == kMode1overMa2) {
    int sign_mashtwk = utils::rew::Sign(fMaTwkDial);
    double fracerr_mash =
        fracerr->OneSigmaErr(kXSecTwkDial_1overMaNCEL2, sign_mashtwk);
    fMaCurr = 1. / sqrt(1. + fMaTwkDial * fracerr_mash);
  }

  // else
  // if(fMode==kModeNormAndMaShape) {
  int sign_normtwk = utils::rew::Sign(fNormTwkDial);
  double fracerr_norm =
      fracerr->OneSigmaErr(kXSecTwkDial_NormNCEL, sign_normtwk);
  fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
  //}

  fNormCurr = TMath::Max(0., fNormCurr);
  fMaCurr = TMath::Max(0., fMaCurr);

  int sign_affshtwk = utils::rew::Sign(fAxlFFTwkDial);
  double fracerr_affsh =
      fracerr->OneSigmaErr(kXSecTwkDial_AxlFFNCEL, sign_affshtwk);
  fAxlFFCurr = (1. + fAxlFFTwkDial * fracerr_affsh);

  int sign_vffshtwk = utils::rew::Sign(fVecFFTwkDial);
  double fracerr_vffsh =
      fracerr->OneSigmaErr(kXSecTwkDial_VecFFNCEL, sign_vffshtwk);
  fVecFFCurr = (1. + fVecFFTwkDial * fracerr_vffsh);

  /*
  if (fVecFFTwkDial) fVecFFCurr = fVecFFTwkDial;
  if (fVecFFCurr != fVecFFDef) fortFns->SetMCDefaultVal(kXSecTwkDial_VecFFNCEL,
  fVecFFCurr);

  if (fVecFFOutTwkDial) fVecFFOutCurr = fVecFFOutTwkDial;
  else fVecFFOutCurr = fVecFFCurr;
  */
}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeight() {
  bool is_ccqe = modeDefn.isNCEL(nework_.modene);
  if (!is_ccqe)
    return 1.;

  int nupdg = nework_.ipne[0];
  if (nupdg == kPdgNuMu && !fRewNumu)
    return 1.;
  if (nupdg == kPdgAntiNuMu && !fRewNumubar)
    return 1.;
  if (nupdg == kPdgNuE && !fRewNue)
    return 1.;
  if (nupdg == kPdgAntiNuE && !fRewNuebar)
    return 1.;

  // if(fMode==kModeMa) {
  double wght = this->CalcWeightMa() * this->CalcWeightNorm();
  return wght;
  //}
  // else
  // if(fMode==kModeNormAndMaShape) {
  //   double wght =
  //       this->CalcWeightNorm   () *
  //       this->CalcWeightMaShape();
  //   return wght;
  //}

  return 1.;
}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcChisq() {
  double chisq = 0.;
  // if(fMode==kModeMa) {
  chisq += TMath::Power(fMaTwkDial, 2.);
  //}
  // else
  // if(fMode==kModeNormAndMaShape) {
  chisq += TMath::Power(fNormTwkDial, 2.);
  //   chisq += TMath::Power(fMaTwkDial,   2.);
  //}
  return chisq;
}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeightNorm() {
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if (!tweaked)
    return 1.0;

  double wght = fNormCurr;

  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeightMa() {
  bool tweaked = (TMath::Abs(fMaTwkDial) > controls::kASmallNum ||
                  TMath::Abs(fAxlFFTwkDial) > controls::kASmallNum ||
                  TMath::Abs(fVecFFTwkDial) > controls::kASmallNum);

  if (!tweaked)
    return 1.0;

  fortFns->SetDefaults();
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif
  float old_xsec = fortFns->evdifcrs();
  // float old_xsec   = event.DiffXSec();
  // float old_weight = event.Weight();

  // Correct for the units change between RFG and SF
  if (nemdls_.mdlqe == 402)
    old_xsec *= 1E-23;

  if (old_xsec == 0) {
    // cout << "NReWeightNuXSecNCEL::CalcWeightMa() Warning: old_xsec==0,
    // setting weight to 1" << endl;
    return 1;
  }

  // nemdls_.xmaqe = fMaCurr;//RT: comment out since now common block
  // nemdls_.kapp  = fKapCurr;//should be passing this info through
  fMaDef = nemdls_.xmaqe;
  nemdls_.xmaqe = fMaDef * fMaCurr;

  // MDLQE    : CC Quasi-elastic / NC elastic model
  //          : x1 : Smith-Moniz for CC
  //          : x2 : Smith-Moniz for CC with BBBA05
  //          : 0x : Scaling to NCEL     ( same as 5.0.x )
  //          : 1x : Scaling to Spectrum func. with Dipole (prior to v5.1.2)
  //          : 2x : Scaling to Spectrum func. with BBBA05 (default from v5.1.2)
  /*
    if (fAxlFFCurr) nemdls_.mdlqeaf = fAxlFFCurr;
  */
  fVecFFDef = nemdls_.mdlqe;
  // fVecFFCurr = fVecFFDef;
  fAxlFFDef = nemdls_.mdlqeaf;
  // fAxlFFCurr = TMath::Max(0.,fAxlFFTwkDial);

  // If this dial has been tweaked, if changes the default in NFortFns.cc
  nemdls_.mdlqeaf = fAxlFFDef * fAxlFFCurr;
  nemdls_.mdlqe = fVecFFDef * fVecFFCurr;

  // if (fVecFFOutCurr != fVecFFCurr) nemdls_.mdlqe = fVecFFOutCurr;

  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  fortFns->print_allparams();
#endif
  float new_xsec = fortFns->evdifcrs();

  // Correct for the units change between RFG and SF
  if (nemdls_.mdlqe == 402)
    new_xsec *= 1E-23;

  float new_weight = (new_xsec / old_xsec);
  // float new_weight = old_weight * (new_xsec/old_xsec);

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  cout << "differential cross section (old) = " << old_xsec << endl;
  cout << "differential cross section (new) = " << new_xsec << endl;
#endif

  // Normalize out change in total cross section from MA variations for MA_shape
  // variation
  if (fMode == kModeNormAndMaShape) {

    // Get neutrino type
    int inu = -1;
    if (nework_.ipne[0] == 14)
      inu = neutTotCrs->numu;
    else if (nework_.ipne[0] == -14)
      inu = neutTotCrs->numub;
    else if (nework_.ipne[0] == 12)
      inu = neutTotCrs->nue;
    else if (nework_.ipne[0] == -12)
      inu = neutTotCrs->nueb;
    else {
      cout << "NReWeightNuXSecNCEL::CalcWeightMa() Error: Cannot MA-shape "
              "reweight this neutrino type = "
           << nework_.ipne[0] << endl;
      exit(-1);
    }

    // Calculate neutrino energy
    float Enu = 0; // GeV
    for (int i = 0; i < 3; i++)
      Enu += nework_.pne[0][i] * nework_.pne[0][i];
    Enu = sqrt(Enu);

    // Warning: Does not take into account variations of PFSurf
    float old_tot_xsec =
        0.; // neutTotCrs->ccqe_crs[inu]->Interpolate(Enu, fPfDef, fMaDef);
    float new_tot_xsec =
        0.; // neutTotCrs->ccqe_crs[inu]->Interpolate(Enu, fPfDef, fMaCurr);

    if (new_tot_xsec == 0) {
      cout << "NReWeightNuXSecNCEL::CalcWeightMa() Warning: new_tot_xsec==0, "
              "setting weight to 1"
           << endl;
      return 1;
    }

#ifdef _N_REWEIGHT_NCEL_DEBUG_
    cout << "total cross section (old) = " << old_tot_xsec << endl;
    cout << "total cross section (new) = " << new_tot_xsec << endl;
#endif

    new_weight *= old_tot_xsec / new_tot_xsec;
  }

  if (isinf(new_weight) || isnan(new_weight)) {
    cout << "NReWeightNuXSecNCEL::CalcWeightMa() Warning: new_weight is "
            "infinite, setting to 1"
         << endl;
    return 1;
  }

#ifdef _N_REWEIGHT_NCEL_DEBUG_
  cout << "new weight = " << new_weight << endl;
#endif

  return new_weight;
}
//_______________________________________________________________________________________
double NReWeightNuXSecNCEL::CalcWeightMaShape() {
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

//LOG("ReW", pDEBUG) << "integrated cross section (old) = " <<
old_integrated_xsec;
//LOG("ReW", pDEBUG) << "integrated cross section (new) = " <<
new_integrated_xsec;
//LOG("ReW", pDEBUG) << "new weight (normalized to const integral) = " <<
new_weight;

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
