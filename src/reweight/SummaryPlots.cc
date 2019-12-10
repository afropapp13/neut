#include "EventSummaryHelper.h"
#include "neutrootTreeSingleton.h"

#include "TFile.h"

#include <iostream>
#include <sstream>

int main(int argc, char const *argv[]) {

  if (argc != 3) {
    std::cout
        << "[ERROR]: Expected to be passed: <NEUTVector.root> <outfile.root>"
        << std::endl;
    return 1;
  }

  NeutrootTreeSingleton &neutRdr = NeutrootTreeSingleton::Initialize(argv[1]);

  Long64_t Nents = neutRdr.GetEntries();

  std::map<int, HistBlob> AllTheHists;

  double fweight = neutRdr.GetXsecWeight();

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

    AllTheHists[k.mode].Fill(k, fweight);
    AllTheHists[0].Fill(k, fweight);
  }

  TFile *f = new TFile(argv[2], "RECREATE");
  for (std::map<int, HistBlob>::iterator it = AllTheHists.begin();
       it != AllTheHists.end(); ++it) {
    std::stringstream ss("");
    ss << "mode_" << it->first;
    TDirectory *mode_dir = f->mkdir(ss.str().c_str());
    AllTheHists[it->first].Write(mode_dir, ss.str());
  }
  f->Write();
  f->Close();
}
