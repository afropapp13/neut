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

#include "NReWeightBindingEnergy.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>

#include <TMath.h>
#include <TTree.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TRandom3.h>

#include <specfuncC.h>

#include "Controls.h"
#include "PDGCodes.h"
#include "NSystSet.h"
#include "NSystUncertainty.h"
#include "NReWeightUtils.h"
#include "neutread.h"

extern "C" {
  double ddifcrs_(float *enu, float *co, float *elep, int *ipar);
  void mcmass_(int *pdgcode, float *mass);
  void buildsf_();
  void selectsfvalues_(double *p);
  void pickdirection_(double *fourVector, double *mass);
  void boost_(double *vec, double *v);
  double dotproduct_(double *a, double* b);
  double jacobian_(double *kPrime, double *kPrimeCom, double *pPrime, double *comVel);
  double lh_(double *k, double *kPrime, double *p, double *qTilde,
             int *nuPdg, bool *nc, int *nucPdg,  int *numAtom);
}

namespace neut { namespace rew {

namespace binding_energy {
TargetMaterial::TargetMaterial(int atomic_number_, double default_value_, double uncertainty_)
  : atomic_number(atomic_number_)
  , default_value(default_value_)
  , fractional_uncertainty(uncertainty_/default_value_)
{ }

CarbonTarget::CarbonTarget(double default_value_, double uncertainty_)
  : TargetMaterial(12, default_value_, uncertainty_)
{ }

OxygenTarget::OxygenTarget(double default_value_, double uncertainty_)
  : TargetMaterial(16, default_value_, uncertainty_)
{ }

AluminumTarget::AluminumTarget(double default_value_, double uncertainty_)
  : TargetMaterial(27, default_value_, uncertainty_)
{ }

IronTarget::IronTarget(double default_value_, double uncertainty_)
  : TargetMaterial(56, default_value_, uncertainty_)
{ }

CopperTarget::CopperTarget(double default_value_, double uncertainty_)
  : TargetMaterial(63, default_value_, uncertainty_)
{ }

ZincTarget::ZincTarget(double default_value_, double uncertainty_)
  : TargetMaterial(64, default_value_, uncertainty_)
{ }

LeadTarget::LeadTarget(double default_value_, double uncertainty_)
  : TargetMaterial(208, default_value_, uncertainty_)
{ }
/*namespace binding_energy*/}


NReWeightBindingEnergy::NReWeightBindingEnergy(const binding_energy::TargetMaterial& target,
                                               int desired_mdlqe, double weight_cap, bool debug, bool trace)
  : NReWeightI()
  , fRewNue(true)
  , fRewNuebar(true)
  , fRewNumu(true)
  , fRewNumubar(true)
  , fAbsolute(false)
  , fDebug(debug)
  , fTrace(trace)
  , fDesiredMDLQE(desired_mdlqe)
  , fWeightCap(weight_cap)
  , fEbTwkDial(0)
  , fEbCurr(float(target.default_value))
  , fEbDef(float(target.default_value))
  , fEbUncertainty(float(target.fractional_uncertainty))
  , fNumAtom(target.atomic_number)
  , fSystematic(kNullSystematic)
{
  this->Init();
}
//_______________________________________________________________________________________
NReWeightBindingEnergy::~NReWeightBindingEnergy()
{ }
//_______________________________________________________________________________________
void NReWeightBindingEnergy::Init()
{
  fEbDef *= -1.e-3; //< Convert to -MeV
  fSystematic = GetSystematic();
  NSystUncertainty::Instance()->SetUncertainty(fSystematic, fEbUncertainty, fEbUncertainty);

  if (!IsQEModelValid()) {
    if (fDebug) {
      std::cout << "fDesiredMDLQE = " << fDesiredMDLQE << ", is invalid. Defaulting to 2 (RFG)\n";
    }
    fDesiredMDLQE = 2;
  }

  if (const char* env_neutroot = std::getenv("NEUT_ROOT")) {
      // Yes, only one '='.
      // That assigns env_neutroot if the environment variable exists,
      // and does nothing if it doesn't exist.
      // Furthermore, env_neutroot goes out of scope with the if-block.
    std::string crstblpath = env_neutroot;
    crstblpath += "/src/neutsmpl/";
    std::strncpy(neutfilepath_.crstblpath, crstblpath.c_str(), crstblpath.size());
    std::fill(&neutfilepath_.crstblpath[crstblpath.size()], &neutfilepath_.crstblpath[1024], ' ');
  }
}
//_______________________________________________________________________________________
NSyst_t NReWeightBindingEnergy::GetSystematic() const
{
  switch (fNumAtom) {
    case  12: return kSystNucl_CCQEBindingEnergy_C12;
    case  16: return kSystNucl_CCQEBindingEnergy_O16;
    case  27: return kSystNucl_CCQEBindingEnergy_Al27;
    case  56: return kSystNucl_CCQEBindingEnergy_Fe56;
    case  63: return kSystNucl_CCQEBindingEnergy_Cu63;
    case  64: return kSystNucl_CCQEBindingEnergy_Zn64;
    case 208: return kSystNucl_CCQEBindingEnergy_Pb208;
    default : return kNullSystematic;
  }
}
//_______________________________________________________________________________________
bool NReWeightBindingEnergy::IsQEModelValid() const {
  const int      ones = (fDesiredMDLQE %    10) /    1;
  const int      tens = (fDesiredMDLQE %   100) /   10;
  const int  hundreds = (fDesiredMDLQE %  1000) /  100;
  const int thousands = (fDesiredMDLQE % 10000) / 1000;

  bool is_valid = false;
  if (  (     ones == 1 ||      ones == 2 || ones == 3)
      // Form factor: 1 := dipole,      2 := bbba05, 3 := bbba07

      &&(     tens == 0 ||      tens == 1)
      // MEC:         0 := none,        1 := FGM *= sqrt(1+TRCORA*Q2*exp(-1*Q2/TRCORA));

      &&( hundreds == 0 ) //||  hundreds == 4 || hundreds == 6 || hundreds == 7)
      // TODO: permit non-zero values of SF tuning once we get that part of the dial working
      // Various spectral functions: 0 := off

      &&(thousands == 0 || thousands == 1)
      // RPA:         0 := off          1 := Non-relativistic
     )
  {
    is_valid = true;
  }
  return is_valid;
}

//_______________________________________________________________________________________
bool NReWeightBindingEnergy::IsHandled(NSyst_t syst)
{
  bool is_handled = (syst == fSystematic);
  return is_handled;
}
//_______________________________________________________________________________________
void NReWeightBindingEnergy::SetSystematic(NSyst_t syst, double twk_value)
{
  if (IsHandled(syst)) {
    fEbTwkDial = twk_value;
  } else if (fTrace) {
      std::cout << "["<<__FILE__<<":L"<<__LINE__<<"] "
         << "Systematic '" << NSyst::AsString(syst) << "' is not handled by this dial.\n";
  }
}
//_______________________________________________________________________________________
void NReWeightBindingEnergy::Reset()
{
  fEbTwkDial = 0.;
  NSystUncertainty::Instance()->SetUncertainty(fSystematic, fEbUncertainty, fEbUncertainty);
  this->Reconfigure();
}

//_______________________________________________________________________________________
void NReWeightBindingEnergy::Reconfigure()
{
  const int    sign_ebtwk = utils::rew::Sign(fEbTwkDial);
  const double fracerr_eb = NSystUncertainty::Instance()->OneSigmaErr(fSystematic, sign_ebtwk);
  fEbCurr = 1. + fEbTwkDial * fracerr_eb;
}
//_______________________________________________________________________________________
double NReWeightBindingEnergy::CalcWeight(TTree* tree, long long event) {
  readtree(tree, event);
  if (fTrace) { PrintCommonBlocks(); }
  double weight = CalcWeight();
  return weight;
}
//_______________________________________________________________________________________
double NReWeightBindingEnergy::CalcWeight()
{
  double weight = 1.;
  const bool tweaked = (std::fabs(fEbTwkDial) > 1.e-9);
  if (!tweaked) { return weight; }

  const bool is_ccqe = modeDefn.isCCQE(nework_.modene);
  const bool use_spectral_functions = (((fDesiredMDLQE % 1000) / 100) == 4);
  const bool correct_neutrino = GoodNeutrino();
  const bool nucleon_is_bound = (posinnuc_.ibound != 0);
  const bool correct_target = nucleon_is_bound && (fNumAtom == neuttarget_.numatom);
  if (is_ccqe && correct_target && correct_neutrino && !use_spectral_functions) {
    // The desired mode cannot be one of the spectral functions,
    //  they do not have any dependance on Eb, so it makes no sense to reweight them.
    if (EventIsOkayForDesiredMDLQE()) {

      //{
        // TODO: remove this block once SF->RFG works
        const int actual_mdlqe = nemdls_.mdlqe;
        const int      ones = (actual_mdlqe %    10) /    1;
        const int      tens = (actual_mdlqe %   100) /   10;
        const int  hundreds = (actual_mdlqe %  1000) /  100; // Spectral functions
        const int thousands = (actual_mdlqe % 10000) / 1000;

        // Intentionally omit the hundreds, because we can't currently reweight from SF to RFG
        nemdls_.mdlqe = thousands*1000 + tens*10 + ones*1;
      //}
      weight = CalcWeightI();
      //{
        // TODO: remove this block once SF->RFG works
        nemdls_.mdlqe = actual_mdlqe;
      //}
    } else {
      // If the current event is "CCQE" but not allowed for the desired mode, then
      //  the cross section for the desired mode is defined to be zero.
      weight = 0;
    }
  } else if (fTrace) {
    std::cerr << std::boolalpha << "This event is bad, specifically:\n"
              << "\t                is_ccqe = " << is_ccqe << ",\n"
              << "\t!use_spectral_functions = " << !use_spectral_functions << ",\n"
              << "\t       correct_neutrino = " << correct_neutrino << ",\n"
              << "\t         correct_target = " << correct_target << "(fNumAtom="<<fNumAtom<<" == neuttarget_.numatom="<<neuttarget_.numatom<< ")\n"
              << "All of the above must be true to actually call CalcWeightI()\n";
  }

  if (weight > fWeightCap) { weight = fWeightCap+1.; }
  return weight;
}
//_______________________________________________________________________________________
bool NReWeightBindingEnergy::EventIsOkayForDesiredMDLQE() const
{
  bool good = false;
  const bool use_spectral_functions = (((fDesiredMDLQE % 1000) / 100) == 4);
  if (use_spectral_functions) {
    // Spectral functions can have correlated nucleons in the initial state.
    good = true;

  } else {
    const size_t outgoing_lepton_index = GetIndexOfOutgoingLepton();
    good = (outgoing_lepton_index == 2);
  }
  return good;
}
//_______________________________________________________________________________________
size_t NReWeightBindingEnergy::GetIndexOfOutgoingLepton() const {
  // We will verify that the given particle is not in the initial state
  const int status_for_initial_particle = -1;

  // Skip over the initial neutrino and the first nucleon
  bool still_looking = false;
  size_t index = 2;
  //for (; index != nework_.numne; ++index) {
  while ((still_looking) && (index != nework_.numne)) {
    const int& Status = nework_.iflgne[index];
    const bool good_status = (Status != status_for_initial_particle);

    const int& PID = nework_.ipne[index];
    const bool good_pid = ((PID == kPdgElectron)
                        || (PID == kPdgPositron)
                        || (PID == kPdgMuon)
                        || (PID == kPdgAntiMuon)
                        || (PID == kPdgTau)
                        || (PID == kPdgAntiTau)
                          );

    const bool good = good_status && good_pid;
    if (good && fTrace) {
      std::cout << '\n'
                << "nework_.iflgne[" << index << "] = " << Status << " ( != " << status_for_initial_particle << " )\n"
                << "nework_.ipne[" << index << "]   = " << PID << '\n'
                << std::boolalpha << "EventIsOkayForDesiredMDLQE => " << good << '\n';
    }

    // update the loop conditions (don't update index if we are good)
    still_looking = !good;
    if (still_looking) { ++index; }
  }
  return index;
}
//_______________________________________________________________________________________
bool NReWeightBindingEnergy::GoodNeutrino() const
{
  const int& nupdg = nework_.ipne[0];
  bool good = ((nupdg==kPdgNuMu     && fRewNumu   )
            || (nupdg==kPdgAntiNuMu && fRewNumubar)
            || (nupdg==kPdgNuE      && fRewNue    )
            || (nupdg==kPdgAntiNuE  && fRewNuebar ));
  return good;
}
//_______________________________________________________________________________________
double NReWeightBindingEnergy::CalcChisq()
{
  return 0.;
}
//_______________________________________________________________________________________
double NReWeightBindingEnergy::CalcWeightI()
{
  double new_weight = 1;
  double old_xsec = 1;

  bool use_spectral_functions = (((nemdls_.mdlqe % 1000) / 100) == 4);
  if (use_spectral_functions) {
    old_xsec = CalcXsecSF();
  } else {
    //
    // New SF to RFG reweighting tables include a factor of 10^-9 added to the RFG bins
    // to provide information about the underlying SF cross section in regions where
    // the RFG model is kinematically forbidden.
    //
    // Cancel that factor of 10^-9 here.
    old_xsec = CalcXsecRFG() + 1.e-9;
  }

  //
  // Update the QE mode
  //
  const float original_vnuini = nenupr_.vnuini;
  const int original_mdlqe = nemdls_.mdlqe;
  nemdls_.mdlqe = fDesiredMDLQE;

  // Update the value of the nuclear binding energy, in MeV
  nenupr_.vnuini = fEbDef * fEbCurr;
  if (fDebug) {
    const int    sign_ebtwk = utils::rew::Sign(fEbTwkDial);
    const double fracerr_eb = NSystUncertainty::Instance()->OneSigmaErr(fSystematic, sign_ebtwk);
    std::cout << std::fixed << std::setprecision(3)
              << "Eb (new) = " << nenupr_.vnuini
              << " = " << fEbDef << " *(1 + ( " << fEbTwkDial << " * " << fracerr_eb << "))\n";
  }

  double new_xsec = CalcXsecRFG();

  if (fDebug) {
    std::cout << std::fixed << std::setprecision(3)
              << "NReWeightBindingEnergy::CalcWeightI(): differential cross section (old) = " << old_xsec << '\n'
              << "NReWeightBindingEnergy::CalcWeightI(): differential cross section (new) = " << new_xsec << '\n';
  }
  new_weight = (new_xsec/old_xsec);

  if (isinf(new_weight) || isnan(new_weight)) {
    if (fDebug) {
      std::cout << "NReWeightBindingEnergy::calc_weight() Warning: new_weight is infinite, setting to fWeightCap+1\n";
    }
    new_weight = fWeightCap+1.;
  }
  if (new_weight < 0) {
    if (fTrace) {
      std::cout << "NReWeightBindingEnergy::calc_weight() Warning: new_weight is negative, setting to 0\n"
                << std::scientific << std::setprecision(2)
                << "\t" << new_weight << " = "
                << "[differential cross section (new) = " << new_xsec << "] / "
                << "[differential cross section (old) = " << old_xsec << "]\n";
    }
    new_weight = 0;
  }

  //
  // Set QE mode back to original for the next event
  //
  nenupr_.vnuini = original_vnuini;
  nemdls_.mdlqe = original_mdlqe;


  if (fDebug) {
    std::cout << std::fixed << std::setprecision(4)
              << "NReWeightBindingEnergy::CalcWeightI(): new weight = " << new_weight << "\n\n\n";
  }

  return new_weight;
}
//_______________________________________________________________________________________
double NReWeightBindingEnergy::CalcXsecRFG() {
  // cache the InputLeptonIndex and OutputLeptonIndex
  const size_t ili = 0;
  const size_t oli = GetIndexOfOutgoingLepton();
  float p_nu = std::sqrt(nework_.pne[ili][0]*nework_.pne[ili][0] +
                         nework_.pne[ili][1]*nework_.pne[ili][1] +
                         nework_.pne[ili][2]*nework_.pne[ili][2]);

  float p_mu = std::sqrt(nework_.pne[oli][0]*nework_.pne[oli][0] +
                         nework_.pne[oli][1]*nework_.pne[oli][1] +
                         nework_.pne[oli][2]*nework_.pne[oli][2]);

  const float costheta = (nework_.pne[ili][0]*nework_.pne[oli][0] +
                          nework_.pne[ili][1]*nework_.pne[oli][1] +
                          nework_.pne[ili][2]*nework_.pne[oli][2]) / (p_nu * p_mu);

  const float MeV_per_GeV = 1.e3;
  p_nu *= MeV_per_GeV;
  p_mu *= MeV_per_GeV;

  const float m_nu = 0.;
  const float e_nu = std::sqrt(p_nu*p_nu + m_nu*m_nu);

  float m_mu = 0;
  const int& mu_pdg = nework_.ipne[oli];
  mcmass_(&mu_pdg, &m_mu);
  const float e_mu = std::sqrt(p_mu*p_mu + m_mu*m_mu);

  const int& nu_pdg = nework_.ipne[ili];

  // Convert from [fm^2 / nucleon / MeV] to [fb / nucleon / GeV/c / sr]
  const double fb_per_fm2 = 1.e7;
  const double dcosth_dOmega = 0.5 / std::acos(-1); //< = 1/(2\pi)
  const double dE_dp = (p_mu/MeV_per_GeV) / e_mu;
  const double magic_power_of_ten = 1.e11; //< Needed to get units to match my tables.
  const double xsec = ddifcrs_(&e_nu, &costheta, &e_mu, &nu_pdg) * fb_per_fm2 * dcosth_dOmega * dE_dp * magic_power_of_ten;
  //                          /______/__________/______/
  //                         /
  //                  Fortran wants pointers, so pass address-of values

  if (fDebug) {
    std::cout << std::scientific << std::setprecision(3)
              << "xsec [fb / nucleon / GeV/c / sr] = " << xsec
              << " [= ddifcrs_(enu=" << e_nu
              << ", costh=" << costheta
              << ", emu=" << e_mu
              << ", nupdg=" << nu_pdg
              << ")]\n";
  }

  return xsec;
}

//_______________________________________________________________________________________
double NReWeightBindingEnergy::CalcXsecSF()
{
  const double MeV_per_GeV = 1.e3;
  double k[4] = {0,0,0,0};
  k[1] = nework_.pne[0][0] * MeV_per_GeV;
  k[2] = nework_.pne[0][1] * MeV_per_GeV;
  k[3] = nework_.pne[0][2] * MeV_per_GeV;
  k[0] = std::sqrt((0.) // neutrino is (effectively) massless
                   +(k[1]*k[1]) +(k[2]*k[2]) +(k[3]*k[3]));

  int in_nucleon_pdg  = nework_.ipne[1];
  float dummy = 0;
  mcmass_(&in_nucleon_pdg, &dummy);
  const double in_nucleon_mass = static_cast<double>(dummy);
  double p[4] = {0,0,0,0};
  p[1] = nework_.pne[1][0] * MeV_per_GeV;
  p[2] = nework_.pne[1][1] * MeV_per_GeV;
  p[3] = nework_.pne[1][2] * MeV_per_GeV;
  p[0] = std::sqrt((in_nucleon_mass*in_nucleon_mass)
                   +(p[1]*p[1]) +(p[2]*p[2]) +(p[3]*p[3]));

  const size_t oli = GetIndexOfOutgoingLepton();
  const int& lepton_pdg = nework_.ipne[oli];
  mcmass_(&lepton_pdg, &dummy);
  const double lepton_mass = static_cast<double>(dummy);
  double kPrime[4] = {0,0,0,0};
  kPrime[1] = nework_.pne[oli][0] * MeV_per_GeV;
  kPrime[2] = nework_.pne[oli][1] * MeV_per_GeV;
  kPrime[3] = nework_.pne[oli][2] * MeV_per_GeV;
  kPrime[0] = std::sqrt((lepton_mass*lepton_mass)
                        +(kPrime[1]*kPrime[1]) +(kPrime[2]*kPrime[2]) +(kPrime[3]*kPrime[3]));

  const int& out_nucleon_pdg = nework_.ipne[oli+1];
  mcmass_(&out_nucleon_pdg, &dummy);
  const double out_nucleon_mass = dummy;
  double pPrime[4] = {0,0,0,0};
  pPrime[1] = nework_.pne[oli+1][0] * MeV_per_GeV;
  pPrime[2] = nework_.pne[oli+1][1] * MeV_per_GeV;
  pPrime[3] = nework_.pne[oli+1][2] * MeV_per_GeV;
  pPrime[0] = std::sqrt((out_nucleon_mass*out_nucleon_mass)
                        +(pPrime[1]*pPrime[1]) +(pPrime[2]*pPrime[2]) +(pPrime[3]*pPrime[3]));

  double xsec = CalcXsecSF(k, p, in_nucleon_mass, kPrime, pPrime, out_nucleon_mass);
  return xsec;
}

double NReWeightBindingEnergy::CalcXsecSF(const double k[], const double p[], const double in_nucleon_mass,
                                          const double kPrime[], const double pPrime[], const double out_nucleon_mass)
{
  // Equation 1 from TN184
  double ETilde = out_nucleon_mass + k[0] - pPrime[0];
  double com_momentum[4] = {(k[0] + in_nucleon_mass - ETilde),
                            (k[1] + p[1]),
                            (k[2] + p[2]),
                            (k[3] + p[3])};

  double comVel[3] = {(com_momentum[1] / com_momentum[0]),
                      (com_momentum[2] / com_momentum[0]),
                      (com_momentum[3] / com_momentum[0])};

  double kPrimeCom[4] = {(kPrime[0]),
                         (kPrime[1]),
                         (kPrime[2]),
                         (kPrime[3])};
  boost_(kPrimeCom, comVel);

  // Reduced 4-momentum transfer
  double qTilde[4] = {(pPrime[0] - p[0]),
                      (k[1] - kPrime[1]),
                      (k[2] - kPrime[2]),
                      (k[3] - kPrime[3])};


  const double Gf = 1.16639E-11;
  const double cosThetaC = 0.97418;
  const double Pi = 3.14159265358;

  double initialTerm = Gf*Gf * cosThetaC*cosThetaC
                      / ( 8*Pi*Pi*k[0]*kPrime[0]*p[0]*pPrime[0] );

  double square_momentum_com = kPrimeCom[1]*kPrimeCom[1]
                             + kPrimeCom[2]*kPrimeCom[2]
                             + kPrimeCom[3]*kPrimeCom[3];
  double volumeTerm = 1
                      * 4*Pi*square_momentum_com
                      * jacobian_(kPrime, kPrimeCom, pPrime, comVel);

  bool nc = false;
  double crossSection = initialTerm*volumeTerm;
  crossSection *= lh_(k, kPrime, p, qTilde,
                      &nework_.ipne[0], &nc,
                      &nework_.ipne[1], &neuttarget_.numatom);

  // ----- turn cross-sections into per-nucleon cross-sections in cm^2
  const double hbarc_sq = 3.89e16; // MeV^2 (10^-38 cm^2)
  crossSection *= hbarc_sq;

  return crossSection;
}

//TH3F
TH2F
NReWeightBindingEnergy::GenerateLookupTable(TTree* tree)
{
  TH1::SetDefaultSumw2(true);

  readtree(tree, 0);
  const size_t oli = GetIndexOfOutgoingLepton();
  const int& lepton_pdg = nework_.ipne[oli];
  float lepton_mass = 0;
  mcmass_(&lepton_pdg, &lepton_mass);

  //TODO: use binning from flux histogram
  const int nbins_e = 130;
  const double emin = 2.06; // ~= std::log10(1.1e0 * lepton_mass);
  const double emax = 3.36; // ~= std::log10(2.2e1 * lepton_mass);
  const double de = (emax - emin) / nbins_e;

  const int nbins_p = 130;
  const double pmin = 2.06; // ~= std::log10(1.1e0 * lepton_mass);
  const double pmax = 3.36; // ~= std::log10(2.2e1 * lepton_mass);
  const double dp = (pmax - pmin) / nbins_p;

  const int nbins_t = 90;
  const double tmin = 0.;
  const double tmax = std::acos(-1); // pi
  const double dt = (tmax - tmin) / nbins_t;

  std::string target="Unknown";
  const int& atomic_number = neuttarget_.numatom;
  if      (atomic_number == 12) { target="C12"; }
  else if (atomic_number == 16) { target="O16"; }
  else if (atomic_number == 56) { target="Fe56"; }

  std::string neutrino="Unknown";
  const int& nu_pdg = nework_.ipne[0];
  if      (nu_pdg == kPdgNuE     ) { neutrino="nue";     }
  else if (nu_pdg == kPdgAntiNuE ) { neutrino="nuebar";  }
  else if (nu_pdg == kPdgNuMu    ) { neutrino="numu";    }
  else if (nu_pdg == kPdgAntiNuMu) { neutrino="numubar"; }

  nemdls_.mdlqe = 402;
  //TH3F
  TH2F
    table(Form("%s_%s_%04d", target.c_str(), neutrino.c_str(), nemdls_.mdlqe),
             //";log10(enu [MeV])"
             ";log10(pmu [MeV/c])"
             ";theta [rad]",
              //nbins_e, emin, emax,
              nbins_p, pmin, pmax,
              nbins_t, tmin, tmax);

  TH2F
    raw(Form("raw_%s_%s_%04d", target.c_str(), neutrino.c_str(), nemdls_.mdlqe),
             //";log10(enu [MeV])"
             ";log10(pmu [MeV/c])"
             ";theta [rad]",
              //nbins_e, emin, emax,
              nbins_p, pmin, pmax,
              nbins_t, tmin, tmax);


  const long long entries = tree->GetEntriesFast();
  for (long long entry=0; entry!=entries; ++entry) {
      readtree(tree, entry);

      // Neutrino energy [MeV] == Neutrino momentum [MeV/c]
      double p_nu = std::sqrt(nework_.pne[0][0]*nework_.pne[0][0] +
                              nework_.pne[0][1]*nework_.pne[0][1] +
                              nework_.pne[0][2]*nework_.pne[0][2]);

      // Muon momentum [MeV/c]
      double p_mu = std::sqrt(nework_.pne[oli][0]*nework_.pne[oli][0] +
                              nework_.pne[oli][1]*nework_.pne[oli][1] +
                              nework_.pne[oli][2]*nework_.pne[oli][2]);

      // Opening angle between muon and neutrino
      double costheta = (nework_.pne[0][0]*nework_.pne[oli][0] +
                         nework_.pne[0][1]*nework_.pne[oli][1] +
                         nework_.pne[0][2]*nework_.pne[oli][2]) / (p_nu * p_mu);
      double theta = std::acos(costheta);

      p_mu *= 1.e3; // GeV/c -> MeV/c
      p_nu *= 1.e3; // GeV/c -> MeV/c

      double xsec = CalcXsecSF();
      table.Fill(//std::log10(p_nu),
                 std::log10(p_mu),
                 theta,
                 xsec);
      raw.Fill(//std::log10(p_nu),
                 std::log10(p_mu),
                 theta,
                 1);

      if (entry % 1000 == 0) {
        std::cout << std::scientific << std::setprecision(2)
                  << "table->Fill(p_nu= " << p_nu << ", p_mu=" << p_mu
                  << ", theta=" << theta << ", xsec=" << xsec << ");\n";
      }
  }

  table.Divide(&raw);
  return table;
}

TH3F NReWeightBindingEnergy::GenerateLookupTable() {
  TH1::SetDefaultSumw2(true);

  float dummy = 0;
  const int in_nucleon_pdg = 2112; // N
  mcmass_(&in_nucleon_pdg, &dummy);
  double in_nucleon_mass = static_cast<double>(dummy);

  const int out_nucleon_pdg = 2112; // P
  mcmass_(&out_nucleon_pdg, &dummy);
  double out_nucleon_mass = static_cast<double>(dummy);

  const int neutrino_pdg = 14; // numu
  const int lepton_pdg = 13; // mu^-
  mcmass_(&lepton_pdg, &dummy);
  double lepton_mass = static_cast<double>(dummy);

  nework_.numne = 4;
  nework_.ipne[0] = neutrino_pdg;
  nework_.ipne[1] = in_nucleon_pdg;
  nework_.ipne[2] = lepton_pdg;
  nework_.ipne[3] = out_nucleon_pdg;

  //TODO: use binning from flux histogram
  const int nbins_e = 100;
  const double emin = std::log10(500);
  const double emax = 3.;
  //const double emin = 2.06; // ~= std::log10(1.1e0 * lepton_mass);
  //const double emax = 3.36; // ~= std::log10(2.2e1 * lepton_mass);
  const double de = (emax - emin) / nbins_e;

  const int nbins_p = 130;
  const double pmin = 2.06; // ~= std::log10(1.1e0 * lepton_mass);
  const double pmax = 3.36; // ~= std::log10(2.2e1 * lepton_mass);
  const double dp = (pmax - pmin) / nbins_p;

  const int nbins_t = 90;
  const double tmin = 0.;
  const double tmax = std::acos(-1); // pi
  const double dt = (tmax - tmin) / nbins_t;

  std::string target="Unknown";
  neuttarget_.numatom = fNumAtom;
  if      (fNumAtom == 12) {
    target="C12";
    neuttarget_.numbndn = 6;
    neuttarget_.numbndp = 6;
    neuttarget_.numfrep = 0;
  }
  else if (fNumAtom == 16) {
    target="O16";
    neuttarget_.numbndn = 8;
    neuttarget_.numbndp = 8;
    neuttarget_.numfrep = 0;
  }
  else if (fNumAtom == 56) {
    target="Fe56";
    neuttarget_.numbndn = 30;
    neuttarget_.numbndp = 26;
    neuttarget_.numfrep = 0;
  }

  std::string neutrino="Unknown";
  const int& nu_pdg = 14;
  if      (nu_pdg == kPdgNuE     ) { neutrino="nue";     }
  else if (nu_pdg == kPdgAntiNuE ) { neutrino="nuebar";  }
  else if (nu_pdg == kPdgNuMu    ) { neutrino="numu";    }
  else if (nu_pdg == kPdgAntiNuMu) { neutrino="numubar"; }

  nework_.modene = neutcard_.nemodflg;

  nemdls_.mdlqe   = 402;
  nemdls_.mdlqeaf = 1;
  nemdls_.xmaqe   = 1.21;
  nemdls_.xmvqe   = 0.84;
  nemdls_.sccfv   = 0.;
  nemdls_.sccfa   = 0.;
#if (NECOREVER != 5321) && (NECOREVER != 5322) && (NECOREVER >= 534)
  nemdls_.fpqe    = 1.;
#endif
  nemdls_.pfsf    = -1.;
  nemdls_.kapp    = 1.0;

  neutcard_.nefrmflg  = 0;
  neutcard_.nepauflg  = 0;
  neutcard_.nenefo16  = 0;
  neutcard_.nenefmodh = 1;
  neutcard_.nenefmodl = 1;
  neutcard_.nenefkinh = 1;
  neutcard_.nemodflg  = 1; // != 0
  neutcard_.quiet     = 2;

  nenupr_.pfsurf = 0.225;
  nenupr_.pfmax = 0.225;
  nenupr_.vnuini = fEbDef;
  nenupr_.vnufin = 0.000;
  nenupr_.iformlen = 1;
  nenupr_.fzmu2 = 0.08;

  neffpr_.fefqe   = 1.;
  neffpr_.fefqeh  = 1.8;
  neffpr_.fefinel = 1.;
  neffpr_.fefabs  = 1.1;
  neffpr_.fefcoh  = 1.;
  neffpr_.fefqehf = 1.;
  neffpr_.fefcohf = 0.;
  neffpr_.fefcx   = 1.;
  neffpr_.fefcxh  = 1.8;
  neffpr_.fefcxhf = 0.;
  neffpr_.fefcoul = 0.;
#if (NECOREVER != 5321) && (NECOREVER != 5322) && (NECOREVER >= 535)
  neffpr_.fefall  = 1.;
#endif

  neutpiless_.ipilessdcy = 0;
  neutpiless_.rpilessdcy = 0.;

  neut::PrintCommonBlocks();
  buildsf_();

  TH3F
    table(Form("%s_%s_%04d", target.c_str(), neutrino.c_str(), nemdls_.mdlqe),
        ";log10(enu [MeV])"
        ";log10(pmu [MeV/c])"
        ";theta [rad]",
        nbins_e, emin, emax,
        nbins_p, pmin, pmax,
        nbins_t, tmin, tmax);

  TH3F
    raw(Form("raw_%s_%s_%04d", target.c_str(), neutrino.c_str(), nemdls_.mdlqe),
        ";log10(enu [MeV])"
        ";log10(pmu [MeV/c])"
        ";theta [rad]",
        nbins_e, emin, emax,
        nbins_p, pmin, pmax,
        nbins_t, tmin, tmax);

  float dir[3] = {0,0,1};
  const int events_per_energy_bin = 1e7;
  TRandom3 rng(0);
  for (int ie=1; ie!=nbins_e; ++ie) {
    const double this_min = table.GetXaxis()->GetBinLowEdge(ie);
    const double this_max = table.GetXaxis()->GetBinUpEdge(ie);
    //const double loge = table.GetXaxis()->GetBinCenterLog(ie);
    if (fTrace) {
      std::cout << "E_{#nu}/MeV \\in [" << std::pow(10,this_min) << ", " << std::pow(10,this_max) << ")\n";
    }

    int successful_events = 0;
    while (successful_events != events_per_energy_bin) {
      bool successful_event_so_far = true;
      const double loge = rng.Uniform(this_min, this_max);
      const double e = std::pow(10, loge);

      // Define the neutrino 4-momentum
      double k[4] = {e,e*dir[0],e*dir[1],e*dir[2]};

      // Define the input nucleon 4-momentum
      double p[4] = {0,0,0,0};
      specfunc_.etilde = 0;
      specfunc_.fermimomentum = 0;
      selectsfvalues_(p);
      p[0] = std::sqrt(in_nucleon_mass*in_nucleon_mass
          +(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]));

      // COM momentum
      double com[4] = {(k[0] + in_nucleon_mass + specfunc_.etilde),
        (k[1] + p[1]),
        (k[2] + p[2]),
        (k[3] + p[3])};

      const double com_dot_com = dotproduct_(com, com);
      const double output_mass_sq = lepton_mass*lepton_mass
        + out_nucleon_mass*out_nucleon_mass;
      if (com_dot_com < output_mass_sq) {
        successful_event_so_far = false;
        //continue;
      }

      // Define output lepton 4-momentum
      double kPrime[4] = {0,0,0,0};
      kPrime[0] = (com_dot_com - out_nucleon_mass*out_nucleon_mass + lepton_mass*lepton_mass)
        / (2.*std::sqrt(com_dot_com));
      pickdirection_(kPrime, &lepton_mass);

      double pPrime[4] = {0, -kPrime[1], -kPrime[2], -kPrime[3]};
      pPrime[0] = std::sqrt(out_nucleon_mass*out_nucleon_mass
          +(pPrime[1]*pPrime[1] + pPrime[2]*pPrime[2] + pPrime[3]*pPrime[3]));

      // boost the output particles out of the COM frame into the lab frame
      double com_vel[4] = {0, -com[1]/com[0], -com[2]/com[0], -com[3]/com[0]};
      double kPrimeCom[4] = {kPrime[0], kPrime[1], kPrime[2], kPrime[3]};
      boost_(kPrime, com_vel);
      boost_(pPrime, com_vel);

      // check if pauli blocking applies
      const double final_nucleon_momentum = std::sqrt(pPrime[1]*pPrime[1]
          +pPrime[2]*pPrime[2]
          +pPrime[3]*pPrime[3]);
      if (successful_event_so_far && (final_nucleon_momentum < specfunc_.fermimomentum)) {
        successful_event_so_far = false;
        //continue;
      }

      double xsec = 0;
      if (successful_event_so_far) {
        xsec = CalcXsecSF(k, p, in_nucleon_mass, kPrime, pPrime, out_nucleon_mass);

        if ((xsec < 0) || std::isnan(xsec) || std::isinf(xsec)) {
          successful_event_so_far = false;
          xsec = 0;
        }
      }

      // Neutrino energy [MeV] == neutrino momentum because m_nu := 0.
      double e_nu = k[0];

      // Muon momentum [MeV/c]
      double p_mu = std::sqrt(kPrime[1]*kPrime[1]
          +kPrime[2]*kPrime[2]
          +kPrime[3]*kPrime[3]);

      // Opening angle between muon and neutrino
      double costheta = (k[1]*kPrime[1]
          +k[2]*kPrime[2]
          +k[3]*kPrime[3]) / (e_nu * p_mu);
      double theta = std::acos(costheta);

      table.Fill(std::log10(e_nu),
          std::log10(p_mu),
          theta,
          xsec);

      raw.Fill(std::log10(e_nu),
          std::log10(p_mu),
          theta,
          1);

      if ( fTrace && ((successful_events % 1000) == 0) ) {
        std::cout << std::fixed << std::setprecision(0)
                  << "table->Fill(e_nu= " << std::setw(5) << e_nu << ", p_mu=" << std::setw(5) << p_mu
                  << ", theta (deg)=" << std::setw(5) << theta*(180./std::acos(-1))
                  << std::scientific << std::setprecision(2) << ", xsec=" << xsec << ");\n";
      }

      if (successful_event_so_far) {
        ++successful_events;
      }
    }
  }

  table.Divide(&raw);
  return table;
}

/*namespace rew*/ } /*namespace neut*/ }

