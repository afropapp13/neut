#include "CascadeNucleon.h"

#include "CommonBlockIFace.h"

#include "vcworkC.h"

#include "posinnucC.h"

#include "nrintC.h"
#include "nucleonfsihistC.h"

#include <cmath>
#include <iostream>

//#define _N_REWEIGHT_CASC_DEBUG_


extern "C" {
  void nrprton_(float*,float[3],float[4],int*,float[4],int*,float[3],int*,int*,int*);
}

namespace neut {
namespace rew {

CascadeNucleonEngine::CascadeNucleonEngine() {
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  NucleonFSI_Total   = RegisterDial("NucleonFSI_Total",
                                      nucres_.xnucfact, 0.5, 0.5);
  NucleonFSI_Elastic = RegisterDial("NucleonFSI_Elastic",
                                      nucres_.xnucelafact, 0.5, 0.5);
  NucleonFSI_Single  = RegisterDial("NucleonFSI_Single",
                                      nucres_.xnucspifact, 0.5, 0.5);
  NucleonFSI_Double  = RegisterDial("NucleonFSI_Double",
                                      nucres_.xnucdpifact, 0.5, 0.5);

  #ifdef _N_REWEIGHT_CASC_DEBUG_
    std::cout << "[DEBUG] Initiating CascadeNucleonEngine" << std::endl;
    std::cout << "[DEBUG] Initial  xnucfact    = " << nucres_.xnucfact    << std::endl;
    std::cout << "[DEBUG] Initial  xnucelafact = " << nucres_.xnucelafact << std::endl;
    std::cout << "[DEBUG] Initial  xnucspifact = " << nucres_.xnucspifact << std::endl;
    std::cout << "[DEBUG] Initial  xnucdpifact = " << nucres_.xnucdpifact << std::endl;
  #endif

}

void CascadeNucleonEngine::Reconfigure(void) {
  ApplyTweaks();
  for (auto &dial : fDials) {
    if (dial.second.ToValue < 0) {
      std::cout << "[ERROR]: Dial: " << NSyst::DialAsString(dial.first)
                << " set unphysical value: " << dial.second.ToValue
                << " (CV: " << dial.second.FromValue
                << ", Uncert: " << dial.second.GetOneSigma(dial.second.Tweak)
                << ", Tweak: " << dial.second.Tweak << ")" << std::endl;
      abort();
    }
  }
}

double CascadeNucleonEngine::CalcWeight() {

  // Not bound
  if (!posinnuc_.ibound) {
    return 1;
  }

  if (!AnyTweaked()) {
    return 1;
  }

  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();
  cbfa.ResetGenValues();


  #ifdef _N_REWEIGHT_CASC_DEBUG_
    std::cout << "[DEBUG] Starting CascadeNucleonEngine::CalcWeight()" << std::endl;
    std::cout << "[DEBUG] Pre-weight casc prob   = " << nrint_.pcascprob    << std::endl;
    std::cout << "[DEBUG] Pre-weight xnucfact    = " << nucres_.xnucfact    << std::endl;
    std::cout << "[DEBUG] Pre-weight xnucelafact = " << nucres_.xnucelafact << std::endl;
    std::cout << "[DEBUG] Pre-weight xnucspifact = " << nucres_.xnucspifact << std::endl;
    std::cout << "[DEBUG] Pre-weight xnucdpifact = " << nucres_.xnucdpifact << std::endl;
  #endif

  // Rerun the cascade
  //reinitialise pcascprob to 1 (start of cascade)
  nrint_.pcascprob = 1.0;
  //nrprton just needs some arguments, none of these numbers actually used in the reweight
  float anuc = -999; float stpi[3] = {-999,-999,-999}; float pi[4] = {-999,-999,-999,-999};
  int idmc =-999; float po[4][20]; int ido; int no; float stpo[3][20]; int imode; int icont;
  //tell the code it's reweight o'clock
  nucleonfsihist_.nreweightflag = 1;
  //run nrprton.F for reweighting
  nrprton_(&anuc, stpi, pi, &idmc, *po, &ido, *stpo, &no, &imode, &icont);

  double old_xsec = nrint_.pcascprob;

  if (old_xsec <= 0) {
    std::cout
        << "[WARN]: CascadeNucleonEngine() pcascprob old_xsec <= 0, returning "
           "weight = 1. Updated pcascprob = " << old_xsec << "  " 
        << std::endl;
    return 1;
  }

  nucres_.xnucfact    = TMath::Max(float(0.), float(fDials[NucleonFSI_Total].ToValue));
  nucres_.xnucelafact = TMath::Max(float(0.), float(fDials[NucleonFSI_Elastic].ToValue));
  nucres_.xnucspifact = TMath::Max(float(0.), float(fDials[NucleonFSI_Single].ToValue));
  nucres_.xnucdpifact = TMath::Max(float(0.), float(fDials[NucleonFSI_Double].ToValue));

  NEUTSetParams();

  #ifdef _N_REWEIGHT_CASC_DEBUG_
    std::cout << "[DEBUG] Set reweight parameters to: " << std::endl;
    std::cout << "[DEBUG] xnucfact    = " << nucres_.xnucfact    << std::endl;
    std::cout << "[DEBUG] xnucelafact = " << nucres_.xnucelafact << std::endl;
    std::cout << "[DEBUG] xnucspifact = " << nucres_.xnucspifact << std::endl;
    std::cout << "[DEBUG] xnucdpifact = " << nucres_.xnucdpifact << std::endl;
  #endif

  // Rerun the cascade
  //reinitialise pcascprob to 1 (start of cascade)
  nrint_.pcascprob = 1.0;
  //nrprton just needs some arguments, none of these numbers actually used in the reweight
  //float anuc = -999; float stpi[3] = {-999,-999,-999}; float pi[4] = {-999,-999,-999,-999};
  //int idmc =-999; float po[4][20]; int ido; int no; float stpo[3][20]; int imode; int icont;
  //tell the code it's reweight o'clock
  nucleonfsihist_.nreweightflag = 1;
  //run nrprton.F for reweighting
  nrprton_(&anuc, stpi, pi, &idmc, *po, &ido, *stpo, &no, &imode, &icont);

  double new_xsec = nrint_.pcascprob;

#ifdef _N_REWEIGHT_CASC_DEBUG_
  std::cout << "[DEBUG] Nucleon cascade probability (old) = " << old_xsec << std::endl;
  std::cout << "[DEBUG] Nucleon cascade probability (new) = " << new_xsec << std::endl;
#endif

  return CheckReturnWeight(new_xsec / old_xsec);
}

} // namespace rew
} // namespace neut
