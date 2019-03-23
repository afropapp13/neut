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
   Made static consts `kModeMaMv' and `kModeNormAndMaMvShape' public to
   aid external configuration.
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NReWeightNuXSecCCRES.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_CCRES_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

const int NReWeightNuXSecCCRES::kModeMaMv;
const int NReWeightNuXSecCCRES::kModeNormAndMaShape;

//_______________________________________________________________________________________
NReWeightNuXSecCCRES::NReWeightNuXSecCCRES() 
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuXSecCCRES::~NReWeightNuXSecCCRES()
{
}
//_______________________________________________________________________________________
void NReWeightNuXSecCCRES::Init(void)
{
  fortFns = NFortFns::Instance();

  neutTotCrs = NTotCrs::Instance();

  this->SetMode(kModeMaMv);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;
  /*
  fMaTwkDial   = 0.; 
  fMaDef       = fortFns->XMACCRESdef;
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.; 
  fMvDef       = fortFns->XMVCCRESdef;
  fMvCurr      = fMvDef;

  fMaTwkDial   = 0.; 
  fMaDef       = 0;
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.; 
  fMvDef       = 0;
  fMvCurr      = fMvDef;
  */
  fMaNFFTwkDial   = 0.; 
  fMaNFFDef       = fortFns->XMANFFRESdef;
  fMaNFFCurr      = fMaNFFDef;
  fMvNFFTwkDial   = 0.; 
  fMvNFFDef       = fortFns->XMVNFFRESdef;
  fMvNFFCurr      = fMvNFFDef;
  fMaRSTwkDial   = 0.; 
  fMaRSDef       = fortFns->XMARSRESdef;
  fMaRSCurr      = fMaRSDef;
  fMvRSTwkDial   = 0.; 
  fMvRSDef       = fortFns->XMVRSRESdef;
  fMvRSCurr      = fMvRSDef;
  fCA5TwkDial   = 0.; 
  fCA5Def       = fortFns->RNECA5Idef;
  fCA5Curr      = fCA5Def;
  fBgSclTwkDial   = 0.; 
  fBgSclDef       = fortFns->RNEBGSCLdef;
  fBgSclCurr      = fBgSclDef;
  fIFFTwkDial   = 0; 
  fIFFDef       = fortFns->NEIFFdef;
  fIFFCurr      = fIFFDef;
  fNRTypeTwkDial   = 0; 
  fNRTypeDef       = fortFns->NENRTYPEdef;
  fNRTypeCurr      = fNRTypeDef;


}
//_______________________________________________________________________________________
bool NReWeightNuXSecCCRES::IsHandled(NSyst_t syst)
{
   bool handle;
  
   switch(syst) {

     //case ( kXSecTwkDial_NormCCRES    ) :
     //  //case ( kXSecTwkDial_MaCCRESshape ) :
     //  //case ( kXSecTwkDial_MvCCRESshape ) :
     //  if(fMode==kModeNormAndMaShape) { 
     //     handle = true;  
     //  } else { 
     //     handle = false; 
     //  }
     //  break;

     case ( kXSecTwkDial_MaCCRESshape ) :
       //case ( kXSecTwkDial_MaCCRES ) :
       //case ( kXSecTwkDial_MvCCRES ) :
     case ( kXSecTwkDial_NormCCRES    ) : // Placed here temporarily, in GENIE it's only used with Shape as above
       //if(fMode==kModeMaMv) { 
          handle = true;  
       //} else { 
       //   handle = false; 
       //}
       break;

     //case ( kSystNucl_PilessDcyCCRES ) :
     //  handle = true;
     //  break;

     case ( kXSecTwkDial_MaNFFCCRES ) :
     case ( kXSecTwkDial_MvNFFCCRES ) :
     case ( kXSecTwkDial_MaRSCCRES ) :
     case ( kXSecTwkDial_MvRSCCRES ) :
     case ( kXSecTwkDial_CA5CCRES ) :
     case ( kXSecTwkDial_BgSclCCRES ) :
     case ( kXSecTwkDial_FFCCRES ) :
     case ( kXSecTwkDial_TypeCCRES ) :
       handle = true;
     break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
void NReWeightNuXSecCCRES::SetSystematic(NSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_NormCCRES ) :
      fNormTwkDial = twk_dial;
      break;


    // Following kludge necessary when using MaCCQEshape and MaCCQE in same
      /*
    // instantiation (i.e. for weight tree generation      
    case ( kXSecTwkDial_MaCCRESshape ) :
      if (fMode==kModeNormAndMaShape)
	fMaTwkDial = twk_dial;
      break;

     case ( kXSecTwkDial_MaCCRES ) :
       if (fMode==kModeMaMv)
	 fMaTwkDial = twk_dial;
       break;

      //case ( kXSecTwkDial_MvCCRESshape ) :
     case ( kXSecTwkDial_MvCCRES ) :
       fMvTwkDial = twk_dial;
       break;
      */
     case ( kXSecTwkDial_MaNFFCCRES ) :
       fMaNFFTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_MvNFFCCRES ) :
       fMvNFFTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_MaRSCCRES ) :
       fMaRSTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_MvRSCCRES ) :
       fMvRSTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_CA5CCRES ) :
       fCA5TwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_BgSclCCRES ) :
       fBgSclTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_FFCCRES ) :
       fIFFTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_TypeCCRES ) :
       fNRTypeTwkDial = twk_dial;
       break;

    default:
      break;
  }
}
//_______________________________________________________________________________________
void NReWeightNuXSecCCRES::Reset(void)
{
  fNormTwkDial = 0.;
  fNormCurr    = fNormDef;
  /*
  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.; 
  fMvCurr      = fMvDef;
  */
  fMaNFFTwkDial   = 0.; 
  fMaNFFCurr      = fMaNFFDef;
  fMvNFFTwkDial   = 0.; 
  fMvNFFCurr      = fMvNFFDef;
  fMaRSTwkDial   = 0.; 
  fMaRSCurr      = fMaRSDef;
  fMvRSTwkDial   = 0.; 
  fMvRSCurr      = fMvRSDef;
  fCA5TwkDial   = 0.; 
  fCA5Curr      = fCA5Def;
  fBgSclTwkDial   = 0.; 
  fBgSclCurr      = fBgSclDef;
  fIFFTwkDial   = 0; 
  fIFFCurr      = fIFFDef;
  fNRTypeTwkDial   = 0; 
  fNRTypeCurr      = fNRTypeDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuXSecCCRES::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  if(fMode==kModeMaMv) {   
     double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCCRES);
     double fracerr_mv = fracerr->OneSigmaErr(kXSecTwkDial_MvCCRES);
     fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
     fMvCurr = fMvDef * (1. + fMvTwkDial * fracerr_mv);
  }

  else if(fMode==kModeNormAndMaShape) { 
    double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaCCRESshape);
    fMaCurr   = fMaDef   * (1. + fMaTwkDial   * fracerr_mash);    
    //double fracerr_mvsh = fracerr->OneSigmaErr(kXSecTwkDial_MvCCRESshape);
    //fMvCurr   = fMvDef   * (1. + fMvTwkDial   * fracerr_mvsh);
  }

     double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormCCRES);
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);     


     double fracerr_nffma = fracerr->OneSigmaErr(kXSecTwkDial_MaNFFCCRES);
     double fracerr_nffmv = fracerr->OneSigmaErr(kXSecTwkDial_MvNFFCCRES);
     fMaNFFCurr = fMaNFFDef * (1. + fMaNFFTwkDial * fracerr_nffma);
     fMvNFFCurr = fMvNFFDef * (1. + fMvNFFTwkDial * fracerr_nffmv);

     double fracerr_rsma = fracerr->OneSigmaErr(kXSecTwkDial_MaRSCCRES);
     double fracerr_rsmv = fracerr->OneSigmaErr(kXSecTwkDial_MvRSCCRES);
     fMaRSCurr = fMaRSDef * (1. + fMaRSTwkDial * fracerr_rsma);
     fMvRSCurr = fMvRSDef * (1. + fMvRSTwkDial * fracerr_rsmv);

     double fracerr_ca5 = fracerr->OneSigmaErr(kXSecTwkDial_CA5CCRES);
     double fracerr_bgscl = fracerr->OneSigmaErr(kXSecTwkDial_BgSclCCRES);
     fCA5Curr = fCA5Def * (1. + fCA5TwkDial * fracerr_ca5);
     fBgSclCurr = fBgSclDef * (1. + fBgSclTwkDial * fracerr_bgscl);

     double fracerr_iff = fracerr->OneSigmaErr(kXSecTwkDial_FFCCRES);
     fIFFCurr = int(fIFFDef + fracerr_iff * fIFFTwkDial);

     double fracerr_type = fracerr->OneSigmaErr(kXSecTwkDial_TypeCCRES); 
     fNRTypeCurr = int(fNRTypeDef + fracerr_type * fNRTypeTwkDial);
     //fIFFCurr   = TMath::Max(0, fIFFTwkDial  );
     //fNRTypeCurr   = TMath::Max(0, fNRTypeTwkDial  );
     fIFFCurr = TMath::Max(0, fIFFCurr);
  fNormCurr = TMath::Max(0., fNormCurr);
  fMaCurr   = TMath::Max(0., fMaCurr  );
  fMvCurr   = TMath::Max(0., fMvCurr  );
  fMaNFFCurr   = TMath::Max(0., fMaNFFCurr  );
  fMvNFFCurr   = TMath::Max(0., fMvNFFCurr  );
  fMaRSCurr   = TMath::Max(0., fMaRSCurr  );
  fMvRSCurr   = TMath::Max(0., fMvRSCurr  );
  fCA5Curr   = TMath::Max(0., fCA5Curr  );
  fBgSclCurr   = TMath::Max(0., fBgSclCurr  );

}
//_______________________________________________________________________________________
double NReWeightNuXSecCCRES::CalcWeight() 
{
  bool is_res = modeDefn.isCCRES(nework_.modene);
  
  if(!is_res) return 1.;

  double wght = 1;

  int nupdg = nework_.ipne[0];
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  //if(fMode==kModeMaMv) {
    wght *= this->CalcWeightMaMv()*this->CalcWeightNorm();
  //} 
  //else 
  //if(fMode==kModeNormAndMaShape) {
  //   wght *= 
  //       this->CalcWeightNorm      () *
  //       this->CalcWeightMaMvShape ();
  //}

  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuXSecCCRES::CalcChisq()
{
  double chisq = 0.;
  //if(fMode==kModeMaMv) {   
     chisq += TMath::Power(fMaTwkDial, 2.);
     chisq += TMath::Power(fMvTwkDial, 2.);
  //}
  //else
  //if(fMode==kModeNormAndMaShape) { 
     chisq += TMath::Power(fNormTwkDial, 2.);
  //   chisq += TMath::Power(fMaTwkDial,   2.);
  //   chisq += TMath::Power(fMvTwkDial,   2.);
  //}  

  return chisq;
}
//_______________________________________________________________________________________
double NReWeightNuXSecCCRES::CalcWeightNorm() 
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;
  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuXSecCCRES::CalcWeightMaMv() 
{
  bool tweaked = 
     (TMath::Abs(fMaNFFTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvNFFTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMaRSTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvRSTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fCA5TwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fBgSclTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fIFFTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fNRTypeTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  
  fortFns->SetDefaults();
  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_CCRES_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif
  
  float old_xsec     = fortFns->evdifcrs();
  //float old_xsec   = event.DiffXSec();
  //float old_weight = event.Weight();

  if (old_xsec==0) {
    //cout << "NReWeightNuXSecCCRES::CalcWeightMa() Warning: old_xsec==0, setting weight to 1" << endl;
    return 1;
  }

  //nemdls_.xmaspi = fMaCurr;
  //nemdls_.xmvspi = fMvCurr;
  neut1pi_.xmanffres = fMaNFFCurr;
  neut1pi_.xmvnffres = fMvNFFCurr;
  neut1pi_.xmarsres = fMaRSCurr;
  neut1pi_.xmvrsres = fMvRSCurr;
  neut1pi_.neiff    = fIFFCurr;
  neut1pi_.nenrtype = fNRTypeCurr;
  neut1pi_.rneca5i  = fCA5Curr;
  neut1pi_.rnebgscl = fBgSclCurr;

  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_CCRES_DEBUG_
  fortFns->print_allparams();
#endif

  float new_xsec   = fortFns->evdifcrs();
  float new_weight = (new_xsec/old_xsec);
  //float new_weight = old_weight * (new_xsec/old_xsec);

#ifdef _N_REWEIGHT_CCRES_DEBUG_
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
      cout << "NReWeightNuXSecCCRES::CalcWeightMa() Error: Cannot MA-shape reweight this neutrino type = " 
	   << nework_.ipne[0] << endl;
      exit (-1);
    }
    
    // Calculate neutrino energy
    float Enu = 0;  // GeV
    for (int i=0; i<3; i++)
      Enu += nework_.pne[0][i]*nework_.pne[0][i];
    Enu = sqrt(Enu);

    // Get reaction mode
    int imode=-1;
    if      (abs(nework_.modene)==11) imode = neutTotCrs->neutmode11;
    else if (abs(nework_.modene)==12) imode = neutTotCrs->neutmode12;
    else if (abs(nework_.modene)==13) imode = neutTotCrs->neutmode13;
    else if (abs(nework_.modene)==31) imode = neutTotCrs->neutmode31;
    else if (abs(nework_.modene)==32) imode = neutTotCrs->neutmode32;
    else if (abs(nework_.modene)==33) imode = neutTotCrs->neutmode33;
    else if (abs(nework_.modene)==34) imode = neutTotCrs->neutmode34;
    else {
      //  cout << "NReWeightNuXSecCCRES::CalcWeightMa() Warning: Cannot reweight MaShape for mode = " 
      //       << nework_.modene << ", setting weight = 1" << endl;
      return 1;
    }
    
    float old_tot_xsec = neutTotCrs->resspi_crs[inu][imode]->Interpolate(Enu, fMaDef);
    float new_tot_xsec = neutTotCrs->resspi_crs[inu][imode]->Interpolate(Enu, fMaCurr);

    if (new_tot_xsec==0) {
      cout << "NReWeightNuXSecCCRES::CalcWeightMa() Warning: new_tot_xsec==0, setting weight to 1" << endl;
      return 1;
    }

#ifdef _N_REWEIGHT_CCRES_DEBUG_
    cout << "total cross section (old) = " << old_tot_xsec << endl;
    cout << "total cross section (new) = " << new_tot_xsec << endl;
#endif 

    new_weight *= old_tot_xsec / new_tot_xsec ;
  }

  if (isinf(new_weight) || isnan(new_weight)) {
    cout << "NReWeightNuXSecCCRES::CalcWeightMa() Warning: new_weight is infinite, setting to 1" << endl;
    new_weight = 1;
  }

#ifdef _N_REWEIGHT_CCRES_DEBUG_
  cout << "new weight = " << new_weight << endl;
#endif

  return new_weight;
}

//_______________________________________________________________________________________
double NReWeightNuXSecCCRES::CalcWeightMaMvShape() 
{
  return 1.0;

  /* //GENIE Implementation below
  bool tweaked = 
     (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSWQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

//double old_integrated_xsec = event.XSec();
  double old_integrated_xsec = fXSecModelDef -> Integral(interaction);
  double twk_integrated_xsec = fXSecModel    -> Integral(interaction);
  assert(twk_integrated_xsec > 0);
  new_weight *= (old_integrated_xsec/twk_integrated_xsec);

//LOG("ReW", pDEBUG) << "integrated cross section (old) = " << old_integrated_xsec;
//LOG("ReW", pDEBUG) << "integrated cross section (twk) = " << twk_integrated_xsec;
//LOG("ReW", pDEBUG) << "new weight (normalized to const integral) = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();

#ifdef _N_REWEIGHT_CCRES_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  double W  = interaction->Kine().W(true);
  fTestNtp->Fill(E,Q2,W,new_weight);
#endif

  return new_weight;
  */
}
