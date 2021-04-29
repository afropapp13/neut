#include "CommonBlockIFace.h"
#include "NReWeightFactory.h"
#include "neutrootTreeSingleton.h"

#include "TFile.h"

#include <iostream>
#include <sstream>

using namespace neut;
using namespace neut::rew;

int main(int argc, char const *argv[]) {

  if (argc < 4) {
    std::cout << "[ERROR]: Expected to be passed: <NEUTCard.card> "
                 "<NEUTVector.root> <outfile.root> [nevents]"
              << std::endl;
    return 1;
  }

  CommonBlockIFace::Initialize(argv[1]);
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  auto NRW = MakeNReWeightInstance();

  NeutrootTreeSingleton &neutRdr = NeutrootTreeSingleton::Initialize(argv[2]);

  Long64_t Nents = neutRdr.GetEntries();

  Nents = 100000;

  if (argc > 4) {
    Nents = std::atoi(argv[4]);
  }

  std::cout << cbfa.ParamsToString() << std::endl;

  std::vector<NSyst_t> Dials;
  std::vector<std::vector<double> > Dial_Responses;
  std::vector<std::vector<double> > Dial_Knots;
  std::vector<std::string> DialNames;

  double PionFSI_QEHighMomProb_vals[] = {0.1824, 0.7296, 1.2768, 1.824,
                                         2.3712, 2.9184, 3.4656};
  double PionFSI_InelProb_vals[] = {
      -0.501, 0, 0.501, 1.002, 1.503, 2.004, 2.505,
  };

  Dials.push_back(NRW->DialFromString("PionFSI_QEHighMomProb"));
  Dials.push_back(NRW->DialFromString("PionFSI_InelProb"));

  double fweight = neutRdr.GetXsecWeight();

  size_t NDials = Dials.size();

  TFile *f = new TFile(argv[3], "RECREATE");
  TTree *t = new TTree("EventSplines", "");
  TTree *k = new TTree("SplineKnots", "");

  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    DialNames.push_back(NRW->DialAsString(Dials[d_it]));
    Dial_Responses.push_back(std::vector<double>(7, 0));
    Dial_Knots.push_back(std::vector<double>(7, 0));
  }

  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    t->Branch((DialNames[d_it] + "_Response").c_str(), &Dial_Responses[d_it]);
    k->Branch((DialNames[d_it] + "_Knots").c_str(), &Dial_Knots[d_it]);
  }

  size_t shout_every = Nents / 100;

  NeutVect *nvct = neutRdr.GetNeutVectAddress();
  for (Long64_t ent = 0; ent < Nents; ++ent) {
    neutRdr.GetEntry(ent);

    if (ent && !(ent % shout_every)) {
      std::cout << "Processed " << ent << "/" << Nents << std::endl;
    }

    CommonBlockIFace::ReadVect(nvct);

    for (size_t d_it = 0; d_it < NDials; ++d_it) {
      for (int dv = -3; dv < 4; ++dv) {
        NRW->Reset();
        // NRW->SetDial_NumberOfSigmas(Dials[d_it], dv);
        if (d_it == 0) {
          NRW->SetDial_To_Value(Dials[d_it], PionFSI_QEHighMomProb_vals[dv + 3]);
        } else {
          NRW->SetDial_To_Value(Dials[d_it], PionFSI_InelProb_vals[dv + 3]);
        }
        NRW->Reconfigure();
        if (!ent) {
          Dial_Knots[d_it][dv + 3] = NRW->GetDial_To_Value(Dials[d_it]);
        }
        Dial_Responses[d_it][dv + 3] = NRW->CalcWeight();
      }
    }
    t->Fill();
  }

  k->Fill();

  f->Write();
  f->Close();
}
