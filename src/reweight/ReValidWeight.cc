#include "CommonBlockIFace.h"
#include "EventSummaryHelper.h"
#include "NReWeightFactory.h"
#include "NReWeightNuXSecCCQE.h"
#include "neutrootTreeSingleton.h"

#include "TFile.h"

#include <iostream>
#include <sstream>

using namespace neut;
using namespace neut::rew;

int main(int argc, char const *argv[]) {

  if (argc != 4) {
    std::cout << "[ERROR]: Expected to be passed: <NEUTCard.card> "
                 "<NEUTVector.root> <outfile.root>"
              << std::endl;
    return 1;
  }

  CommonBlockIFace::Initialize(argv[1]);
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  NReWeight *NRW = DefaultConfiguredNReWeightFactory();

  double MAQE103_twkval = NRW->WghtCalc("NReWeightNuXSecCCQE")
                              ->GetTwkForAbs(kXSecTwkDial_MaCCQE, 1.03);
  std::cout << "twkval for MAQE103 = " << MAQE103_twkval << std::endl;

  // double alpha_nom = NRW->WghtCalc("NReWeightNuXSecCCQE")
  //                        ->GetTwkForAbs(kXSecTwkDial_FAxlCCQEAlpha, 2.19);
  // double gamma_nom = NRW->WghtCalc("NReWeightNuXSecCCQE")
  //                        ->GetTwkForAbs(kXSecTwkDial_FAxlCCQEGamma, 0.93);
  // double beta_nom = NRW->WghtCalc("NReWeightNuXSecCCQE")
  //                       ->GetTwkForAbs(kXSecTwkDial_FAxlCCQEBeta, -7.67);
  // double theta_nom = NRW->WghtCalc("NReWeightNuXSecCCQE")
  //                        ->GetTwkForAbs(kXSecTwkDial_FAxlCCQETheta, -4.67);
  NeutrootTreeSingleton &neutRdr = NeutrootTreeSingleton::Initialize(argv[2]);

  Long64_t Nents = neutRdr.GetEntries();

  Nents = std::min(Nents, Long64_t(100000));

  std::cout << cbfa.ParamsToString() << std::endl;

  std::vector<std::vector<std::map<int, HistBlob> > > AllTheHists;
  std::vector<NSyst_t> Dials;
  std::vector<std::string> DialNames;

  Dials.push_back(kXSecTwkDial_MaCCQE);

  // Dials.push_back(kXSecTwkDial_MaCCQE);
  // Dials.push_back(kXSecTwkDial_MaRES);
  // Dials.push_back(kXSecTwkDial_CA5RES);
  // Dials.push_back(kXSecTwkDial_BgSclRES);
  // Dials.push_back(kXSecTwkDial_BgSclLMCPiBarRES);
  // Dials.push_back(kCascTwkDial_FrAbs_pi);
  // Dials.push_back(kCascTwkDial_FrInelLow_pi);
  // Dials.push_back(kCascTwkDial_FrInelHigh_pi);
  // Dials.push_back(kCascTwkDial_FrCExLow_pi);
  // Dials.push_back(kCascTwkDial_FrCExHigh_pi);
  // Dials.push_back(kCascTwkDial_FrPiProd_pi);

  double fweight = neutRdr.GetXsecWeight();

  size_t NDials = Dials.size();

  AllTheHists.resize(NDials);

  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    DialNames.push_back(NSyst::AsString(Dials[d_it]));

    if (Dials[d_it] == kXSecTwkDial_AxlFFCCQE) {
      AllTheHists[d_it].resize(4);
    } else {
      AllTheHists[d_it].resize(7);
    }
  }

  size_t shout_every = Nents / 100;

  NeutVect *nvct = neutRdr.GetNeutVectAddress();
  for (Long64_t ent = 0; ent < Nents; ++ent) {
    neutRdr.GetEntry(ent);

    if (ent && !(ent % shout_every)) {
      std::cout << "Processed " << ent << "/" << Nents << std::endl;
    }

    Kine const &k = GetKine(nvct);
    if (!k.foundlep) {
      continue;
    }

    CommonBlockIFace::ReadVect(nvct);

    for (size_t d_it = 0; d_it < NDials; ++d_it) {
      int numtweaks = AllTheHists[d_it].size();
      for (int dv = -3; dv < (numtweaks - 3); ++dv) {
        NRW->Systematics().Clear();

        if (Dials[d_it] == kXSecTwkDial_BgSclLMCPiBarRES) {
          NRW->Systematics().Init(kXSecTwkDial_UseSeparateBgSclLMCPiBar, 1);
          NRW->Systematics().Init(kXSecTwkDial_BgSclLMCPiBarRES, dv);

        } else if (Dials[d_it] == kXSecTwkDial_AxlFFCCQE) {
          if (dv == -3) {
            // std::cout << "@@Dipole:" << std::endl;
            NRW->Systematics().Init(kXSecTwkDial_AxlFFCCQE, 1);
            NRW->Systematics().Init(kXSecTwkDial_MaCCQE, MAQE103_twkval);
          } else if (dv == -2) {
            // std::cout << "@@3comp:" << std::endl;
            NRW->Systematics().Init(kXSecTwkDial_AxlFFCCQE, 4);
          } else if (dv == -1) {
            // NRW->Systematics().Init(kXSecTwkDial_AxlFFCCQE, 4);
            // NRW->Systematics().Init(kXSecTwkDial_FAxlCCQEAlpha, alpha_nom);
            // NRW->Systematics().Init(kXSecTwkDial_FAxlCCQEGamma, gamma_nom);
            // NRW->Systematics().Init(kXSecTwkDial_FAxlCCQEBeta, beta_nom);
            // NRW->Systematics().Init(kXSecTwkDial_FAxlCCQETheta, theta_nom);

          } else if (dv == 0) {
            // std::cout << "@@zexp:" << std::endl;
            NRW->Systematics().Init(kXSecTwkDial_AxlFFCCQE, 5);
          }

        } else {
          NRW->Systematics().Init(Dials[d_it], dv);
        }

        NRW->Reconfigure();
        // std::cout << "Enu: " << k.Enu << ", Q2: " << k.Q2 << std::endl;
        double w = std::min(NRW->CalcWeight(), 10.0);
        if (!std::isnormal(w)) {
          w = 0;
        }
        AllTheHists[d_it][dv + 3][k.mode].Fill(k, fweight * w);
        AllTheHists[d_it][dv + 3][0].Fill(k, fweight * w);
      }
    }
  }

  TFile *f = new TFile(argv[3], "RECREATE");
  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    TDirectory *dial_dir = f->mkdir(DialNames[d_it].c_str());
    for (std::map<int, HistBlob>::iterator it = AllTheHists[d_it][0].begin();
         it != AllTheHists[d_it][0].end(); ++it) {
      std::stringstream ss("");
      ss << "mode_" << it->first;
      TDirectory *mode_dir = dial_dir->mkdir(ss.str().c_str());
      int numtweaks = AllTheHists[d_it].size();

      for (int dv = -3; dv < (numtweaks - 3); ++dv) {
        ss.str("");
        if (Dials[d_it] == kXSecTwkDial_AxlFFCCQE) {
          if (dv == -3) {
            ss << "tweak_dipole";
          } else if (dv == -2) {
            ss << "tweak_3comp";
          } else if (dv == -1) {
            ss << "tweak_3comp_altfit";
          } else if (dv == 0) {
            ss << "tweak_zexp";
          }
        } else {
          ss << "tweak_" << (dv < 0 ? "m" : "") << std::abs(dv);
        }
        AllTheHists[d_it][dv + 3][it->first].Write(mode_dir, ss.str());
      }
    }
  }
  f->Write();
  f->Close();
}
