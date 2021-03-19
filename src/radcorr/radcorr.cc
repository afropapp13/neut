// Konosuke Iwamoto
// June 8, 2014

// DESCRIPTION:

// Script to read in muon & neutrino kinematics in text file, attach photon
// under certain probability

// 3/31/2014 17:43 EST - Completed radiative correction script using
// interpolation and integrals of 2D distributions and splines that are already
// made; no need to run over text file multiple times. 6/8/2014 17:43 EST -
// Completed radiative correction script using interpolation for probability
// density of photon as well.
// 7/8/2014  14:58 EST - Turned off FV constraint
// 7/30/2014 15:37 EST - nue rad CCQE added
// 10/29/2014 16:56 JST - Modified for NEUT
// 11/11/2014 16:32 JST - Completed the modification for NEUT
// 2020/11/20 12:09 CT - Luke began rewrite

#include "spline_definitions.h"

#include "necardC.h"
#include "neworkC.h"
#include "vcworkC.h"

#include "posinnucC.h"

#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include <iostream>
#include <sstream>

// #define RADCORRDEBUG

extern "C" {
void radcorr();
void radcorr_();
}

void radcorr_() {
  switch (neutradcorr_.iradcorr) {
  case 1:
    radcorr();
    break;
  case 0:
    break;

  default:
    std::cerr << "Undefined iradcorr ( neut card ) value was set"
              << neutradcorr_.iradcorr << std::endl;
    exit(1);
  }
  return;
}

bool isnu(int pdg) {
  pdg = std::abs(pdg);
  return (pdg == 12) || (pdg == 14) || (pdg == 16);
}

bool isclep(int pdg) {
  pdg = std::abs(pdg);
  return (pdg == 11) || (pdg == 13) || (pdg == 15);
}

struct radcorr_inputs {

  TGraph2D *colinear_gamma_total_prob[2];
  TGraph2D *colinear_Egamma_prob[24][2];

  double egamma_max;
  double egamma_min;

  void Init() {
    egamma_max = 3.5;
    egamma_min = 0.03;
    Read();

    for (int is_ele = 0; is_ele < 2; ++is_ele) {
      colinear_gamma_total_prob[is_ele]->SetMarginBinsContent(-999);
      for (int i = 0; i < 24; ++i) {
        colinear_Egamma_prob[is_ele][i]->SetMarginBinsContent(-999);
      }
    }
  }

  double GetProb(double enu, int pdglep, double thetalep) {
    int is_ele = (pdglep == 11);

    if (enu < egamma_knots[0]) {
      return 0;
    }

#ifdef RADCORRDEBUG
    // std::cout << "[RADCORR]: Prob for enu: " << enu << ", pdglep: " << pdglep
    //           << ", thetalep: " << thetalep << " = "
    //           << colinear_gamma_total_prob[is_ele]->Interpolate(enu,
    //                                                             cos(thetalep))
    //           << std::endl;
#endif
    double prob =
        colinear_gamma_total_prob[is_ele]->Interpolate(enu, cos(thetalep));

    if (prob >= 0) {
      return prob;
    }
    // Here we're in the margins, so lets do a 1D extrapolation

    TGraph FixedEnuSpline;
    double min_ct = -0.95, max_ct = 0.99;
    int nknots = 25;
    for (int i = 0; i < nknots; ++i) {
      double ct = min_ct + double(i) * (max_ct - min_ct) / double(nknots);
      FixedEnuSpline.SetPoint(
          i, ct, colinear_gamma_total_prob[is_ele]->Interpolate(enu, ct));
    }
    prob = FixedEnuSpline.Eval(cos(thetalep));
    // std::cout << "[INFO]: At costhetalep: " << cos(thetalep)
    //           << " was outside bound, used 1D extrapolation to determine
    //           prob: "
    //           << prob << "(prob at boundary: "
    //           << colinear_gamma_total_prob[is_ele]->Interpolate(enu, max_ct)
    //           << ")" << std::endl;
    return prob;
  }

  double GetEGamma(double enu, int pdglep, double thetalep, double EGammaMax) {

    if (EGammaMax <= egamma_min) {
      return 0;
    }

    int is_ele = (pdglep == 11);

    if (enu < egamma_knots[0]) {
      std::cout << "[ERROR]: Attempted to get radiative EGamma distribution "
                   "for neutrino of energy: "
                << enu << " GeV." << std::endl;
      abort();
    }

    TGraph EGamma_dist;
    double acc_rej_ceiling = 0;
    for (int i = 0; i < negamma_knots; ++i) {
      double prob =
          colinear_Egamma_prob[is_ele][i]->Interpolate(enu, cos(thetalep));
      if (prob < 0) { // Need to extrapolate out
        TGraph FixedEnuSpline;
        double min_ct = -0.95, max_ct = 0.99;
        int nknots = 25;
        for (int i = 0; i < nknots; ++i) {
          double ct = min_ct + double(i) * (max_ct - min_ct) / double(nknots);
          FixedEnuSpline.SetPoint(
              i, ct, colinear_Egamma_prob[is_ele][i]->Interpolate(enu, ct));
        }
        prob = FixedEnuSpline.Eval(cos(thetalep));
      }
      EGamma_dist.SetPoint(i, egamma_knots[i], prob);
      acc_rej_ceiling = std::max(acc_rej_ceiling, prob * 1.1);
    }
    if (acc_rej_ceiling <= 0) {
      std::cout
          << "[ERROR]: Failed to build a valid EGamma distribution for Enu: "
          << enu << ", thetalep: " << thetalep
          << ", acc_rej_ceiling: " << acc_rej_ceiling << std::endl;
      abort();
    }

    double this_egamma_max = std::min(egamma_max, EGammaMax);

    size_t attempts = 0;
    while (attempts < 1E5) {
      double egamma_trial =
          egamma_min + gRandom->Uniform(0, 1) * (this_egamma_max - egamma_min);
      double rngthrow = gRandom->Uniform(0, 1) * acc_rej_ceiling;
      if (EGamma_dist.Eval(egamma_trial) > rngthrow) {

#ifdef RADCORRDEBUG
        std::cout << "[RADCORR]: Chose egamma: " << egamma_trial << " after "
                  << attempts
                  << " attemps with prob: " << EGamma_dist.Eval(egamma_trial)
                  << "/" << acc_rej_ceiling << " (throw: " << rngthrow << ")"
                  << std::endl;
#endif

        return egamma_trial;
      }
      attempts++;
    }
    std::cout << "[WARN]: Failed to choose an EGamma for Enu: " << enu
              << ", thetalep: " << thetalep << " after " << attempts
              << " attempts." << std::endl;
    return 0;
  }

  void Read() {
    char *NEUT_CRSPATH = getenv("NEUT_CRSPATH");
    if (!NEUT_CRSPATH) {
      std::cout << "[ERROR]: Cannot read environment variable NEUT_CRSPATH, is "
                   "the NEUT environment set up?"
                << std::endl;
      abort();
    }

    TFile *MCfile = new TFile(
        (std::string(NEUT_CRSPATH) + "/radcorr/radcorr_inputs.root").c_str());

    std::stringstream ss("");
    for (int is_ele = 0; is_ele < 2; ++is_ele) {
      ss.str("");
      ss << (is_ele ? "electron" : "muon") << "_colinear_gamma_total_prob";

      MCfile->GetObject(ss.str().c_str(), colinear_gamma_total_prob[is_ele]);
      if (!colinear_gamma_total_prob[is_ele]) {
        std::cout << "[ERROR]: Failed to read radcorr inputs TGraph2D: "
                  << ss.str() << " from " << std::string(NEUT_CRSPATH)
                  << "/radcorr/radcorr_inputs.root" << std::endl;
        abort();
      }

      for (int i = 0; i < negamma_knots; ++i) {
        ss.str("");
        ss << (is_ele ? "electron" : "muon") << "_colinear_gamma_prob_"
           << egamma_knots[i] << "_GeV";

        MCfile->GetObject(ss.str().c_str(), colinear_Egamma_prob[is_ele][i]);
        if (!colinear_Egamma_prob[is_ele][i]) {
          std::cout << "[ERROR]: Failed to read radcorr inputs plot: "
                    << ss.str() << " from " << std::string(NEUT_CRSPATH)
                    << "/radcorr/radcorr_inputs.root" << std::endl;
          abort();
        }
      }
    }
  }
};

void radcorr() {
  // Only for CCQE
  if (std::abs(nework_.modene) != 1) {
    return;
  }

  static radcorr_inputs inputs;
  static bool first = true;
  if (first) {
    inputs.Init();
    first = false;
  }
  // Store necessary variables to calculate probability and photon energy

  // Find important particles from vcwork block
  int nu_idx = -1, lep_idx = -1;
  for (int i = 0; i < vcwork_.nvc; ++i) {

    // std::cout << "[part] " << i << ", p: [" << vcwork_.pvc[i][0] << ", "
    //           << vcwork_.pvc[i][1] << ", " << vcwork_.pvc[i][2]
    //           << "], pdg: " << vcwork_.ipvc[i]
    //           << " status: " << vcwork_.iflgvc[i] << ", " <<
    //           vcwork_.icrnvc[i]
    //           << std::endl;

    // look for initial state neutrino until you've found it
    if ((vcwork_.iflgvc[i] == -1) && (nu_idx == -1)) {
      if (isnu(vcwork_.ipvc[i])) {
        nu_idx = i;
        // std::cout << "--found nu" << std::endl;
      }
      // look for final state charged lepton until you've found one
    } else if ((vcwork_.iflgvc[i] == 0) && (vcwork_.icrnvc[i] == 1) &&
               (lep_idx == -1)) {
      if (isclep(vcwork_.ipvc[i])) {
        lep_idx = i;
        // std::cout << "--found lep" << std::endl;
      }
    }
  }
  if ((nu_idx == -1) || (lep_idx == -1)) {
    std::cout << "[WARN]: Failed to find: nu " << (nu_idx == -1)
              << ", lep: " << (lep_idx == -1) << std::endl;
    return;
  }

  // Only for numu and nue
  bool isnue = (vcwork_.ipvc[nu_idx] == 12);
  bool isnumu = (vcwork_.ipvc[nu_idx] == 14);

  if (!isnue && !isnumu) {
    return;
  }

  TLorentzVector p4_nu_mev, p4_lep_mev;
  p4_nu_mev.SetXYZM(vcwork_.pvc[nu_idx][0], vcwork_.pvc[nu_idx][1],
                    vcwork_.pvc[nu_idx][2], 0);
  p4_lep_mev.SetXYZM(vcwork_.pvc[lep_idx][0], vcwork_.pvc[lep_idx][1],
                     vcwork_.pvc[lep_idx][2], vcwork_.amasvc[lep_idx]);

  // Do not attach
  if (gRandom->Uniform(0, 1) > inputs.GetProb(p4_nu_mev.E() * 1.0E-3,
                                              vcwork_.ipvc[nu_idx],
                                              p4_lep_mev.Theta())) {
    return;
  }

  double KELep = p4_lep_mev.E() - p4_lep_mev.M();
  double EGamma = inputs.GetEGamma(p4_nu_mev.E() * 1.0E-3, vcwork_.ipvc[nu_idx],
                                   p4_lep_mev.Theta(), KELep * 1E-3) *
                  1.0E3;
  if (!EGamma) {
    return;
  }

#ifdef RADCORRDEBUG
  std::cout << "[RADCORR]: before plep: " << p4_lep_mev.Vect().Mag()
            << ", mlep: " << p4_lep_mev.M() << std::endl;
#endif

#ifdef RADCORRDEBUG
  std::cout << "KELep before: " << KELep << std::endl;
#endif
  KELep -= EGamma;
#ifdef RADCORRDEBUG
  std::cout << "KELep after: " << KELep << std::endl;
#endif
  TVector3 plep = p4_lep_mev.Vect().Unit() *
                  sqrt(pow(KELep + p4_lep_mev.M(), 2) - p4_lep_mev.M2());
  p4_lep_mev.SetXYZM(plep[0], plep[1], plep[2], p4_lep_mev.M());

  TLorentzVector p4_gamma_mev =
      TLorentzVector(p4_lep_mev.Vect().Unit() * EGamma, EGamma);

#ifdef RADCORRDEBUG
  std::cout << "[RADCORR]: after plep: " << p4_lep_mev.Vect().Mag()
            << ", mlep: " << p4_lep_mev.M() << std::endl;
#endif

  vcwork_.ipvc[vcwork_.nvc] = 22;
  vcwork_.amasvc[vcwork_.nvc] = 0.;
  vcwork_.iorgvc[vcwork_.nvc] = 3;
  vcwork_.iflgvc[vcwork_.nvc] = 0;
  vcwork_.icrnvc[vcwork_.nvc] = 1;
  vcwork_.timvc[vcwork_.nvc] = vcwork_.timvc[lep_idx];
  vcwork_.ivtivc[vcwork_.nvc] = vcwork_.ivtivc[lep_idx];
  vcwork_.ivtfvc[vcwork_.nvc] = vcwork_.ivtfvc[lep_idx];

  for (int i = 0; i < 3; i++) {
    vcwork_.pvc[vcwork_.nvc][i] = p4_gamma_mev[i];
    vcwork_.posivc[vcwork_.nvc][i] = vcwork_.posivc[lep_idx][i];
    vcwork_.posfvc[vcwork_.nvc][i] = vcwork_.posfvc[lep_idx][i];
    posinnuc_.posnuc[vcwork_.nvc][i] = posinnuc_.posnuc[lep_idx][i];
    vcwork_.pvc[lep_idx][i] = p4_lep_mev[i];
  }

  vcwork_.nvc += 1;
}