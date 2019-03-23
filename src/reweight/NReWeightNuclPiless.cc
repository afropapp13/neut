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
 @ Aug 22, 2011 - PD
   Implemented for NEUT with pionless Delta decay
*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NReWeightNuclPiless.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"

//#define _N_REWEIGHT_NUCLPILESS_DEBUG_

using namespace neut;
using namespace neut::rew;

using std::cout;

//_______________________________________________________________________________________
NReWeightNuclPiless::NReWeightNuclPiless() 
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightNuclPiless::~NReWeightNuclPiless()
{
}
//_______________________________________________________________________________________
void NReWeightNuclPiless::Init(void)
{
  fortFns = NFortFns::Instance();

  this->SetModeEnu(0);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fPilessDcyTwkDial = 0.;   
  fPilessDcyDef     = 0.2; // 20% NEUT default
  fPilessDcyCurr    = fPilessDcyDef;    


}
//_______________________________________________________________________________________
bool NReWeightNuclPiless::IsHandled(NSyst_t syst)
{
   bool handle;
  
   switch(syst) {

     case ( kSystNucl_PilessDcyRES ) :
       handle = true;
       break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
void NReWeightNuclPiless::SetSystematic(NSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;
  switch(syst) {
    case ( kSystNucl_PilessDcyRES ) :
      fPilessDcyTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void NReWeightNuclPiless::Reset(void)
{

  fPilessDcyTwkDial = 0.;   		
  fPilessDcyCurr    = fPilessDcyDef;    

  this->Reconfigure();
}
//_______________________________________________________________________________________
void NReWeightNuclPiless::Reconfigure(void)
{
  NSystUncertainty * fracerr = NSystUncertainty::Instance();

  double fracerr_pilessdcy = fracerr->OneSigmaErr(kSystNucl_PilessDcyRES);
  fPilessDcyCurr = fPilessDcyDef + fPilessDcyTwkDial * fracerr_pilessdcy;  // Adding fractions
  fPilessDcyCurr = TMath::Max(0., fPilessDcyCurr  );

}
//_______________________________________________________________________________________
double NReWeightNuclPiless::CalcWeight() 
{
  bool is_res = modeDefn.isRES(nework_.modene);
  
  if(!is_res) return 1.;

  double wght = 1;

  int nupdg = nework_.ipne[0];
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  // Pionless Delta decay reweighting, but only pion modes
  wght = CalcWeightPilessDcy();

  return wght;
}
//_______________________________________________________________________________________
double NReWeightNuclPiless::CalcChisq()
{
  double chisq = 0.;

  chisq += TMath::Power(fPilessDcyTwkDial, 2.);

  return chisq;
}

//_______________________________________________________________________________________
double NReWeightNuclPiless::CalcWeightPilessDcy() 
{
    bool tweaked = 
     (TMath::Abs(fPilessDcyTwkDial) > controls::kASmallNum) ;
  if(!tweaked) return 1.0;
  
#ifdef _N_REWEIGHT_NUCLPILESS_DEBUG_
  cout << "CalcWeightPilessDcy(): mode = " << nework_.modene << ", ibound = " << posinnuc_.ibound << '\n';
#endif

  // Only 1-pi production events
  if(!modeDefn.is1PI(nework_.modene)) return 1;

  // Do not weight interactions not on bound nucleon
  if (posinnuc_.ibound == 0) return 1;
  
  double new_weight = 1;
  
  // Calculate neutrino energy
  float Enu = 0;  // GeV
  for (int i=0; i<3; i++)
    Enu += nework_.pne[0][i]*nework_.pne[0][i];
  Enu = sqrt(Enu);

  // Find Pionless Delta decay events
  if (nework_.ipne[3]==kPdgP33m1232_DeltaPP || nework_.ipne[3]==kPdgP33m1232_DeltaP
	|| nework_.ipne[3]==kPdgP33m1232_Delta0 || nework_.ipne[3]==kPdgP33m1232_DeltaM) {
    
    if (fabs(fModeEnu)<0.001 || Enu<fabs(fModeEnu)) 
      new_weight = fPilessDcyCurr/fPilessDcyDef;
    
  } else if (fModeEnu>=0) {
    
    if (fabs(fModeEnu)<0.001 || Enu<fabs(fModeEnu)) 
      new_weight = (1-fPilessDcyCurr)/(1-fPilessDcyDef);
    
  }
  
#ifdef _N_REWEIGHT_NUCLPILESS_DEBUG_
  cout << "NReWeightNuclPiless::CalcWeightPilessDcy(): fPilessDcyCurr = " << fPilessDcyCurr << ", fPilessDcyDef = " << fPilessDcyDef << '\n';
  cout << "NReWeightNuclPiless::CalcWeightPilessDcy(): new_weight = " << new_weight << '\n';
#endif


  return new_weight;
}
//_______________________________________________________________________________________
