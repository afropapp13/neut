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
#include "NReWeightNuXSecRES.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_RES_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

const int NReWeightNuXSecRES::kModeMaMv;
const int NReWeightNuXSecRES::kModeNormAndMaShape;

//_______________________________________________________________________________________
NReWeightNuXSecRES::NReWeightNuXSecRES() 
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuXSecRES::~NReWeightNuXSecRES()
{
}
//_______________________________________________________________________________________
void NReWeightNuXSecRES::Init(void)
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
  fMaDef       = fortFns->XMARESdef;
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.; 
  fMvDef       = fortFns->XMVRESdef;
  fMvCurr      = fMvDef;
  */
  fMaNFFTwkDial   = 0.; 
  fMaNFFDef       = 0.; 
  fMaNFFCurr      = fMaNFFDef;
  fMvNFFTwkDial   = 0.; 
  fMvNFFDef       = 0.; 
  fMvNFFCurr      = fMvNFFDef;
  fMaRSTwkDial   = 0.; 
  fMaRSDef       = 0.; 
  fMaRSCurr      = fMaRSDef;
  fMvRSTwkDial   = 0.; 
  fMvRSDef       = 0.; 
  fMvRSCurr      = fMvRSDef;
  fCA5TwkDial   = 0.; 
  fCA5Def       = 0.; 
  fCA5Curr      = fCA5Def;
  fBgSclTwkDial   = 0.; 
  fBgSclDef       = 0.; 
  fBgSclCurr      = fBgSclDef;
  fIFFTwkDial   = 0; 
  fIFFDef       = 0.; 
  fIFFCurr      = fIFFDef;
  fNRTypeTwkDial   = 0; 
  fNRTypeDef       = 0.; 
  fNRTypeCurr      = fNRTypeDef;


}
//_______________________________________________________________________________________
bool NReWeightNuXSecRES::IsHandled(NSyst_t syst)
{
   bool handle;
  
   switch(syst) {

     //case ( kXSecTwkDial_NormRES    ) :
     //  //case ( kXSecTwkDial_MaRESshape ) :
     //  //case ( kXSecTwkDial_MvRESshape ) :
     //  if(fMode==kModeNormAndMaShape) { 
     //     handle = true;  
     //  } else { 
     //     handle = false; 
     //  }
     //  break;

     case ( kXSecTwkDial_MaRESshape ) :
       //case ( kXSecTwkDial_MaRES ) :
       //case ( kXSecTwkDial_MvRES ) :
     case ( kXSecTwkDial_NormRES    ) : // Placed here temporarily, in GENIE it's only used with Shape as above
       //if(fMode==kModeMaMv) { 
          handle = true;  
       //} else { 
       //   handle = false; 
       //}
       break;

     //case ( kSystNucl_PilessDcyRES ) :
     //  handle = true;
     //  break;

     case ( kXSecTwkDial_MaNFFRES ) :
     case ( kXSecTwkDial_MvNFFRES ) :
     case ( kXSecTwkDial_MaRSRES ) :
     case ( kXSecTwkDial_MvRSRES ) :
     case ( kXSecTwkDial_CA5RES ) :
     case ( kXSecTwkDial_BgSclRES ) :
     case ( kXSecTwkDial_FFRES ) :
     case ( kXSecTwkDial_TypeRES ) :
       handle = true;
     break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
void NReWeightNuXSecRES::SetSystematic(NSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_NormRES ) :
      fNormTwkDial = twk_dial;
      break;


    // Following kludge necessary when using MaCCQEshape and MaCCQE in same
      /*
    // instantiation (i.e. for weight tree generation      
    case ( kXSecTwkDial_MaRESshape ) :
      if (fMode==kModeNormAndMaShape)
	fMaTwkDial = twk_dial;
      break;

     case ( kXSecTwkDial_MaRES ) :
       if (fMode==kModeMaMv)
	 fMaTwkDial = twk_dial;
       break;

      //case ( kXSecTwkDial_MvRESshape ) :
     case ( kXSecTwkDial_MvRES ) :
       fMvTwkDial = twk_dial;
       break;
      */
     case ( kXSecTwkDial_MaNFFRES ) :
       fMaNFFTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_MvNFFRES ) :
       fMvNFFTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_MaRSRES ) :
       fMaRSTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_MvRSRES ) :
       fMvRSTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_CA5RES ) :
       fCA5TwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_BgSclRES ) :
       fBgSclTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_FFRES ) :
       fIFFTwkDial = twk_dial;
       break;

     case ( kXSecTwkDial_TypeRES ) :
       fNRTypeTwkDial = twk_dial;
       break;

    default:
      break;
  }
}
//_______________________________________________________________________________________
void NReWeightNuXSecRES::Reset(void)
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
void NReWeightNuXSecRES::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  if(fMode==kModeMaMv) {   
     double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaRES);
     double fracerr_mv = fracerr->OneSigmaErr(kXSecTwkDial_MvRES);
     fMaCurr = (1. + fMaTwkDial * fracerr_ma);
     fMvCurr = (1. + fMvTwkDial * fracerr_mv);
  }

  else if(fMode==kModeNormAndMaShape) { 
    double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaRESshape);
    fMaCurr   = fMaDef   * (1. + fMaTwkDial   * fracerr_mash);    
    //double fracerr_mvsh = fracerr->OneSigmaErr(kXSecTwkDial_MvRESshape);
    //fMvCurr   = fMvDef   * (1. + fMvTwkDial   * fracerr_mvsh);
  }

     double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormRES);
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);     


     double fracerr_nffma = fracerr->OneSigmaErr(kXSecTwkDial_MaNFFRES);
     double fracerr_nffmv = fracerr->OneSigmaErr(kXSecTwkDial_MvNFFRES);
     fMaNFFCurr = (1. + fMaNFFTwkDial * fracerr_nffma);
     fMvNFFCurr = (1. + fMvNFFTwkDial * fracerr_nffmv);

     double fracerr_rsma = fracerr->OneSigmaErr(kXSecTwkDial_MaRSRES);
     double fracerr_rsmv = fracerr->OneSigmaErr(kXSecTwkDial_MvRSRES);
     fMaRSCurr = (1. + fMaRSTwkDial * fracerr_rsma);
     fMvRSCurr = (1. + fMvRSTwkDial * fracerr_rsmv);

     double fracerr_ca5 = fracerr->OneSigmaErr(kXSecTwkDial_CA5RES);
     double fracerr_bgscl = fracerr->OneSigmaErr(kXSecTwkDial_BgSclRES);
     fCA5Curr = (1. + fCA5TwkDial * fracerr_ca5);
     fBgSclCurr = (1. + fBgSclTwkDial * fracerr_bgscl);

     double fracerr_iff = fracerr->OneSigmaErr(kXSecTwkDial_FFRES);
     fIFFCurr = int(fracerr_iff * fIFFTwkDial);

     double fracerr_type = fracerr->OneSigmaErr(kXSecTwkDial_TypeRES); 
     fNRTypeCurr = int(fracerr_type * fNRTypeTwkDial);
     //fIFFCurr   = TMath::Max(0, fIFFTwkDial  );
     //fNRTypeCurr   = TMath::Max(0, fNRTypeTwkDial  );
     /*
     fIFFCurr = TMath::Max(-1, fIFFCurr);
  fNormCurr = TMath::Max(-1., fNormCurr);
  fMaCurr   = TMath::Max(-1., fMaCurr  );
  fMvCurr   = TMath::Max(-1., fMvCurr  );
  fMaNFFCurr   = TMath::Max(-1., fMaNFFCurr  );
  fMvNFFCurr   = TMath::Max(-1., fMvNFFCurr  );
  fMaRSCurr   = TMath::Max(-1., fMaRSCurr  );
  fMvRSCurr   = TMath::Max(-1., fMvRSCurr  );
  fCA5Curr   = TMath::Max(-1., fCA5Curr  );
  fBgSclCurr   = TMath::Max(-1., fBgSclCurr  );
     */
}
//_______________________________________________________________________________________
double NReWeightNuXSecRES::CalcWeight() 
{
  bool is_res = modeDefn.isRES(nework_.modene);
  
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
double NReWeightNuXSecRES::CalcChisq()
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
double NReWeightNuXSecRES::CalcWeightNorm() 
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;
  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuXSecRES::CalcWeightMaMv() 
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

#ifdef _N_REWEIGHT_RES_DEBUG_
  fortFns->print_allevent();
  fortFns->print_allparams();
#endif
  
  float old_xsec     = fortFns->evdifcrs();
  //float old_xsec   = event.DiffXSec();
  //float old_weight = event.Weight();

  if (old_xsec==0) {
    //cout << "NReWeightNuXSecRES::CalcWeightMa() Warning: old_xsec==0, setting weight to 1" << endl;
    return 1;
  }

  //nemdls_.xmaspi = fMaCurr;
  //nemdls_.xmvspi = fMvCurr;
  fMaNFFDef = neut1pi_.xmanffres;
  fMvNFFDef = neut1pi_.xmvnffres;
  fMaRSDef = neut1pi_.xmarsres;
  fMvRSDef = neut1pi_.xmvrsres;
  fIFFDef = neut1pi_.neiff;
  fNRTypeDef = neut1pi_.nenrtype;
  fCA5Def = neut1pi_.rneca5i;
  fBgSclDef = neut1pi_.rnebgscl;

  neut1pi_.xmanffres = TMath::Max(0.,fMaNFFCurr * fMaNFFDef);
  neut1pi_.xmvnffres = TMath::Max(0.,fMvNFFCurr * fMvNFFDef);
  neut1pi_.xmarsres = TMath::Max(0.,fMaRSCurr * fMaRSDef);
  neut1pi_.xmvrsres = TMath::Max(0.,fMvRSCurr * fMvRSDef);
  neut1pi_.neiff    = TMath::Max(0,fIFFCurr + fIFFDef);
  neut1pi_.nenrtype = TMath::Max(0,fNRTypeCurr + fNRTypeDef);
  neut1pi_.rneca5i  = TMath::Max(0.,fCA5Curr * fCA5Def);
  neut1pi_.rnebgscl = TMath::Max(0.,fBgSclCurr * fBgSclDef);

  fortFns->Reconfigure();

#ifdef _N_REWEIGHT_RES_DEBUG_
  fortFns->print_allparams();
#endif

  float new_xsec   = fortFns->evdifcrs();
  float new_weight = (new_xsec/old_xsec);
  //float new_weight = old_weight * (new_xsec/old_xsec);

#ifdef _N_REWEIGHT_RES_DEBUG_
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
      cout << "NReWeightNuXSecRES::CalcWeightMa() Error: Cannot MA-shape reweight this neutrino type = " 
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
      //  cout << "NReWeightNuXSecRES::CalcWeightMa() Warning: Cannot reweight MaShape for mode = " 
      //       << nework_.modene << ", setting weight = 1" << endl;
      return 1;
    }
    
    float old_tot_xsec = neutTotCrs->resspi_crs[inu][imode]->Interpolate(Enu, fMaDef);
    float new_tot_xsec = neutTotCrs->resspi_crs[inu][imode]->Interpolate(Enu, fMaCurr);

    if (new_tot_xsec==0) {
      cout << "NReWeightNuXSecRES::CalcWeightMa() Warning: new_tot_xsec==0, setting weight to 1" << endl;
      return 1;
    }

#ifdef _N_REWEIGHT_RES_DEBUG_
    cout << "total cross section (old) = " << old_tot_xsec << endl;
    cout << "total cross section (new) = " << new_tot_xsec << endl;
#endif 

    new_weight *= old_tot_xsec / new_tot_xsec ;
  }

  if (isinf(new_weight) || isnan(new_weight)) {
    cout << "NReWeightNuXSecRES::CalcWeightMa() Warning: new_weight is infinite, setting to 1" << endl;
    new_weight = 1;
  }

#ifdef _N_REWEIGHT_RES_DEBUG_
  cout << "new weight = " << new_weight << endl;
#endif

  return new_weight;
}

//_______________________________________________________________________________________
double NReWeightNuXSecRES::CalcWeightMaMvShape() 
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

#ifdef _N_REWEIGHT_RES_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  double W  = interaction->Kine().W(true);
  fTestNtp->Fill(E,Q2,W,new_weight);
#endif

  return new_weight;
  */
}
