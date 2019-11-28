#include "CommonBlockIFace.h"
#include "NReWeightFactory.h"
#include "neutrootTreeSingleton.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>

using namespace neut;
using namespace neut::rew;

struct Kine {
  Kine()
      : Q2(0), W(0), Enu(0), plep(0), coslep(0), phmfspi(0), coshmfspi(0),
        Tpi(0), npi(0), mode(0), foundlep(false), foundpi(false) {}
  double Q2;
  double W;
  double Enu;
  double plep;
  double coslep;
  double phmfspi;
  double coshmfspi;
  double Tpi;
  double npi;
  int mode;
  bool foundlep;
  bool foundpi;
};

struct HistBlob {
  HistBlob() {
    Enu = new TH1D("Enu", ";E_{#nu} (GeV);Count", 100, 0, 6);
    Enu->SetDirectory(0);
    Q2 = new TH1D("Q2", ";Q^{2} (GeV^{2});Count", 100, 0, 6);
    Q2->SetDirectory(0);
    W = new TH1D("W", ";W (GeV);Count", 100, 0, 4);
    W->SetDirectory(0);
    pthetalep =
        new TH2D("pthetalep", ";#it{p}_{lep.} (GeV);#theta_{lep.};Count", 100,
                 0, 3, 100, -1, 1);
    pthetalep->SetDirectory(0);
    pthetahmfspi =
        new TH2D("pthetahmfspi", ";#it{p}_{#pi} (GeV);#theta_{#pi};Count", 100,
                 0, 1.5, 100, -1, 1);
    pthetahmfspi->SetDirectory(0);
    Tpi = new TH1D("Tpi", ";T_{#pi} (GeV);Count", 100, 0, 6);
    Tpi->SetDirectory(0);
  }
  void Write(TDirectory *d, std::string suff) {
    d->WriteTObject(Enu, (std::string("Enu_") + suff).c_str());
    d->WriteTObject(Q2, (std::string("Q2_") + suff).c_str());
    d->WriteTObject(W, (std::string("W_") + suff).c_str());
    d->WriteTObject(pthetalep, (std::string("pthetalep_") + suff).c_str());
    d->WriteTObject(pthetahmfspi,
                    (std::string("pthetahmfspi_") + suff).c_str());
    d->WriteTObject(Tpi, (std::string("Tpi_") + suff).c_str());
  }
  void Fill(Kine const &k, double w) {
    Enu->Fill(k.Enu, w);
    Q2->Fill(k.Q2, w);
    W->Fill(k.W, w);
    pthetalep->Fill(k.plep, k.coslep, w);
    pthetahmfspi->Fill(k.phmfspi, k.coshmfspi, k.foundpi * w);
    Tpi->Fill(k.Tpi, k.foundpi * w);
  }
  TH1D *Enu;
  TH1D *Q2;
  TH1D *W;
  TH2D *pthetalep;
  TH2D *pthetahmfspi;
  TH1D *Tpi;
};

Kine GetKine(NeutVect *nvect) {
  TLorentzVector ISP4 = nvect->PartInfo(0)->fP;

  NeutPart const &sn1 = *nvect->PartInfo(1);
  if ((!sn1.fIsAlive && (sn1.fStatus == -1)) &&
      ((sn1.fPID == 2212) || (sn1.fPID == 2112))) {
    ISP4 += sn1.fP;
  }
  NeutPart const &sn2 = *nvect->PartInfo(2);
  if ((!sn2.fIsAlive && (sn2.fStatus == -1)) &&
      ((sn2.fPID == 2212) || (sn2.fPID == 2112))) {
    ISP4 += sn2.fP;
  }

  Kine k;
  NeutPart const &probe = *nvect->PartInfo(0);
  k.Enu = probe.fP.E() * 1E-3;
  k.mode = nvect->Mode;

  for (size_t i = 2; i < nvect->Npart(); ++i) {
    NeutPart const &np = *nvect->PartInfo(i);
    // NC events are special
    if (std::abs(nvect->Mode) > 30) {
      if (!k.foundlep && (abs(np.fPID) > 10 && abs(np.fPID) <= 16) &&
          (np.fStatus == 2)) {
        // pass
      } else if (!np.fIsAlive || (np.fStatus != 0)) { // Standard conditions
        continue;
      }
    } else {
      if (!np.fIsAlive || (np.fStatus != 0)) {
        continue;
      }
    }
    // Calculate Q2 from the first outgoing lepton in the stack (makes Q2 for
    // MEC events correct)
    if (!k.foundlep && (abs(np.fPID) > 10 && abs(np.fPID) <= 16)) {
      k.Q2 = -(np.fP - probe.fP).Mag2() * 1E-6;
      k.plep = np.fP.Vect().Mag() * 1E-3;
      k.coslep = np.fP.Vect().CosTheta();
      k.W = (ISP4 - np.fP).Mag() * 1E-3;
      k.foundlep = true;
    } else if ((std::abs(np.fPID) == 211) || (np.fPID == 111)) {
      k.Tpi += (np.fP.E() - np.fP.M()) * 1E-3;
      if (np.fP.Vect().Mag() > k.phmfspi) {
        k.phmfspi = np.fP.Vect().Mag() * 1E-3;
        k.coshmfspi = np.fP.Vect().CosTheta();
      }
      k.foundpi = true;
      k.npi++;
    }
  }

  return k;
}

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

  NeutrootTreeSingleton &neutRdr = NeutrootTreeSingleton::Initialize(argv[2]);

  Long64_t Nents = neutRdr.GetEntries();

  std::vector<std::vector<std::map<int, HistBlob>>> AllTheHists;
  std::vector<NSyst_t> Dials;
  std::vector<std::string> DialNames;

  // Dials.push_back(kXSecTwkDial_MaCCQE);
  Dials.push_back(kXSecTwkDial_MaRES);
  // Dials.push_back(kXSecTwkDial_CA5RES);
  // Dials.push_back(kXSecTwkDial_BgSclRES);

  size_t NDials = Dials.size();

  AllTheHists.resize(NDials);

  for (size_t d_it = 0; d_it < NDials; ++d_it) {
    DialNames.push_back(NSyst::AsString(Dials[d_it]));

    AllTheHists[d_it].resize(7);
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

    // std::cout << "ENu: " << k.Enu << ", W2: " << pow(k.W, 2) << ", Q2: " << k.Q2
    //           << std::endl;

    for (size_t d_it = 0; d_it < NDials; ++d_it) {
      for (int dv = -3; dv < 4; ++dv) {
        NRW->Systematics().Clear();
        NRW->Systematics().Init(Dials[d_it], dv);
        NRW->Reconfigure();
        AllTheHists[d_it][dv + 3][k.mode].Fill(k, NRW->CalcWeight());
        AllTheHists[d_it][dv + 3][0].Fill(k, NRW->CalcWeight());
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
