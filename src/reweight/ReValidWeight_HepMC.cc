#include "CommonBlockIFace.h"
#include "EventSummaryHelper.h"
#include "NReWeightFactory.h"
#include "NuHepMC/ReaderTools"
#include "neutrootTreeSingleton.h"

#include "TFile.h"

#include <iostream>
#include <sstream>

using namespace neut;
using namespace neut::rew;

int main(int argc, char const *argv[]) {

  if (argc != 3) {
    std::cout
        << "[ERROR]: Expected to be passed: <NEUTVector.root> <outfile.root>"
        << std::endl;
    return 1;
  }

  auto fNuHepMCReader = std::make_unique<NuHepMC::ReaderRootTree>(argv[1]);
  auto fGenInfoHelper = std::make_unique<NuHepMC::genruninfo::GRIHelper>(
      fNuHepMCReader->run_info());

  double fFileWeight = fGenInfoHelper->GetFluxAverageTotalCrossSection() /
                       double(fNuHepMCReader->get_entries());

  CommonBlockIFace::Initialize(fNuHepMCReader->run_info());
  CommonBlockIFace const &cbfa = CommonBlockIFace::Get();

  NReWeight *NRW = DefaultConfiguredNReWeightFactory();

  Long64_t Nents = fNuHepMCReader->get_entries();

  std::cout << cbfa.ParamsToString() << std::endl;

  std::vector<std::vector<std::map<int, HistBlob> > > AllTheHists;
  std::vector<NSyst_t> Dials;
  std::vector<std::string> DialNames;

  Dials.push_back(kXSecTwkDial_MaCCQE);
  Dials.push_back(kXSecTwkDial_MaRES);
  Dials.push_back(kXSecTwkDial_CA5RES);
  Dials.push_back(kXSecTwkDial_BgSclRES);
  Dials.push_back(kXSecTwkDial_BgSclLMCPiBarRES);
  Dials.push_back(kCascTwkDial_FrAbs_pi);
  Dials.push_back(kCascTwkDial_FrInelLow_pi);
  Dials.push_back(kCascTwkDial_FrInelHigh_pi);
  Dials.push_back(kCascTwkDial_FrCExLow_pi);
  Dials.push_back(kCascTwkDial_FrCExHigh_pi);
  Dials.push_back(kCascTwkDial_FrPiProd_pi);

  size_t NDials = Dials.size();

  AllTheHists.resize(NDials);

  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    DialNames.push_back(NSyst::AsString(Dials[d_it]));

    AllTheHists[d_it].resize(7);
  }

  size_t shout_every = Nents / 100;

  std::shared_ptr<NuHepMC::GenEvent> fGenEvent(new NuHepMC::GenEvent());
  for (Long64_t ent = 0; ent < Nents; ++ent) {
    fNuHepMCReader->read_event(*fGenEvent);

    if (ent && shout_every && !(ent % shout_every)) {
      std::cout << "Processed " << ent << "/" << Nents << std::endl;
    }

    Kine const &k = GetKine(fGenEvent);
    if (!k.foundlep) {
      continue;
    }

    CommonBlockIFace::ReadEvent(fGenEvent);

    for (size_t d_it = 0; d_it < NDials; ++d_it) {
      for (int dv = -3; dv < 4; ++dv) {
        NRW->Systematics().Clear();
        NRW->Systematics().Init(Dials[d_it], dv);
        if (Dials[d_it] == kXSecTwkDial_BgSclLMCPiBarRES) {
          NRW->Systematics().Init(kXSecTwkDial_UseSeparateBgSclLMCPiBar, 1);
        }
        NRW->Reconfigure();
        double w = std::min(NRW->CalcWeight(), 10.0);
        AllTheHists[d_it][dv + 3][k.mode].Fill(k, fFileWeight * w);
        AllTheHists[d_it][dv + 3][0].Fill(k, fFileWeight * w);
      }
    }
  }

  TFile *f = new TFile(argv[2], "RECREATE");
  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    TDirectory *dial_dir = f->mkdir(DialNames[d_it].c_str());
    for (std::map<int, HistBlob>::iterator it = AllTheHists[d_it][0].begin();
         it != AllTheHists[d_it][0].end(); ++it) {
      std::stringstream ss("");
      ss << "mode_" << it->first;
      TDirectory *mode_dir = dial_dir->mkdir(ss.str().c_str());
      for (int dv = -3; dv < 4; ++dv) {
        ss.str("");
        ss << "tweak_" << (dv < 0 ? "m" : "") << std::abs(dv);
        AllTheHists[d_it][dv + 3][it->first].Write(mode_dir, ss.str());
      }
    }
  }
  f->Write();
  f->Close();
}
