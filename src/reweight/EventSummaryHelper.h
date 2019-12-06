#ifndef SUMMARYTREEHELPER__SEEN
#define SUMMARYTREEHELPER__SEEN

#include "neutpart.h"
#include "neutvect.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>

struct Kine {
  Kine()
      : Q2(0), W(0), Enu(0), plep(0), coslep(0), phmfspi(0), coshmfspi(0),
        Tpi(0), EMiss(0), npi(0), npip(0), npim(0), npi0(0), nfsiparts(0),
        mode(0), foundlep(false), foundpi(false), ibound(false) {}
  double Q2;
  double W;
  double Enu;
  double plep;
  double coslep;
  double phmfspi;
  double coshmfspi;
  double Tpi;
  double EMiss;
  bool ibound;
  size_t npi;
  size_t npip;
  size_t npim;
  size_t npi0;
  size_t nfsiparts;
  int mode;
  bool foundlep;
  bool foundpi;
};

struct HistReadBlob {
  HistReadBlob() {}
  HistReadBlob(TDirectory *d, std::string suff) {
    d->GetObject((std::string("Enu_") + suff).c_str(), Enu);
    d->GetObject((std::string("Q2_") + suff).c_str(), Q2);
    d->GetObject((std::string("EMiss_") + suff).c_str(), EMiss);
    d->GetObject((std::string("W_") + suff).c_str(), W);
    d->GetObject((std::string("pthetalep_") + suff).c_str(), pthetalep);
    d->GetObject((std::string("pthetahmfspi_") + suff).c_str(), pthetahmfspi);
    d->GetObject((std::string("Tpi_") + suff).c_str(), Tpi);
    d->GetObject((std::string("Topo_") + suff).c_str(), Topo);
    d->GetObject((std::string("Topo_freep_") + suff).c_str(), Topo_freep);
  }
  TH1D *Enu;
  TH1D *Q2;
  TH1D *EMiss;
  TH1D *W;
  TH2D *pthetalep;
  TH2D *pthetahmfspi;
  TH1D *Tpi;
  TH1D *Topo;
  TH1D *Topo_freep;
};

struct HistBlob {
  HistBlob() {
    Enu = new TH1D("Enu", ";E_{#nu} (GeV);d#sigma/d (cm^{2} / GeV)", 40, 0, 6);
    Enu->SetDirectory(0);
    Q2 = new TH1D("Q2", ";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2} / GeV^{2})",
                  40, 0, 2);
    Q2->SetDirectory(0);
    EMiss =
        new TH1D("EMiss", ";E_{miss} (MeV);d#sigma/dE_{miss} (cm^{2} / MeV)",
                 40, -20, 60);
    EMiss->SetDirectory(0);
    W = new TH1D("W", ";W (GeV);d#sigma/dQ^{2} (cm^{2} / GeV^{2})", 25, 0.8, 2);
    W->SetDirectory(0);
    pthetalep = new TH2D("pthetalep",
                         ";#it{p}_{lep.} "
                         "(GeV/c);#theta_{lep.};d#sigma/"
                         "d#it{p}_{lep.}#theta_{lep.} (cm^{2} / GeV/c)",
                         40, 0, 3, 40, -1, 1);
    pthetalep->SetDirectory(0);
    pthetahmfspi = new TH2D("pthetahmfspi",
                            ";#it{p}_{#pi} "
                            "(GeV/c);#theta_{#pi};d#sigma/"
                            "d#it{p}_{#pi}#theta_{#pi}(cm^{2} / GeV/c)",
                            40, 0, 1.5, 40, -1, 1);
    pthetahmfspi->SetDirectory(0);
    Tpi = new TH1D("Tpi", ";T_{#pi} (GeV);d#sigma/dT_{#pi} (cm^{2} / GeV)", 40,
                   0, 2);
    Tpi->SetDirectory(0);
    Topo = new TH1D("Topo", ";;#sigma (cm^{2})", 7, 0, 7);
    Topo->GetXaxis()->SetBinLabel(1, "CC0Pi");
    Topo->GetXaxis()->SetBinLabel(2, "CC0PiNoFSI");
    Topo->GetXaxis()->SetBinLabel(3, "CC1Pi+");
    Topo->GetXaxis()->SetBinLabel(4, "CC1Pi-");
    Topo->GetXaxis()->SetBinLabel(5, "CC1Pi0");
    Topo->GetXaxis()->SetBinLabel(6, "CC1PiNoFSI");
    Topo->GetXaxis()->SetBinLabel(7, "CCNPi");
    Topo->SetDirectory(0);

    Topo_freep = new TH1D("Topo_freep", ";;#sigma (cm^{2})", 5, 0, 5);
    Topo_freep->GetXaxis()->SetBinLabel(1, "CC0Pi");
    Topo_freep->GetXaxis()->SetBinLabel(2, "CC1Pi+");
    Topo_freep->GetXaxis()->SetBinLabel(3, "CC1Pi-");
    Topo_freep->GetXaxis()->SetBinLabel(4, "CC1Pi0");
    Topo_freep->GetXaxis()->SetBinLabel(5, "CCNPi");
    Topo_freep->SetDirectory(0);
  }

  void Write(TDirectory *d, std::string suff) {

    // Enu->Scale(1.0, "WIDTH");
    // Q2->Scale(1.0, "WIDTH");
    // W->Scale(1.0, "WIDTH");
    // pthetalep->Scale(1.0, "WIDTH");
    // pthetahmfspi->Scale(1.0, "WIDTH");
    // Tpi->Scale(1.0, "WIDTH");
    // Topo->Scale(1.0, "WIDTH");
    // Topo_freep->Scale(1.0, "WIDTH");

    d->WriteTObject(Enu, (std::string("Enu_") + suff).c_str());
    d->WriteTObject(Q2, (std::string("Q2_") + suff).c_str());
    d->WriteTObject(EMiss, (std::string("EMiss_") + suff).c_str());
    d->WriteTObject(W, (std::string("W_") + suff).c_str());
    d->WriteTObject(pthetalep, (std::string("pthetalep_") + suff).c_str());
    d->WriteTObject(pthetahmfspi,
                    (std::string("pthetahmfspi_") + suff).c_str());
    d->WriteTObject(Tpi, (std::string("Tpi_") + suff).c_str());
    d->WriteTObject(Topo, (std::string("Topo_") + suff).c_str());
    d->WriteTObject(Topo_freep, (std::string("Topo_freep_") + suff).c_str());
  }
  void Fill(Kine const &k, double w) {
    Enu->Fill(k.Enu, w);
    Q2->Fill(k.Q2, w);
    EMiss->Fill(k.EMiss, w);
    W->Fill(k.W, w);
    pthetalep->Fill(k.plep, k.coslep, w);
    pthetahmfspi->Fill(k.phmfspi, k.coshmfspi, k.foundpi * w);
    Tpi->Fill(k.Tpi, k.foundpi * w);

    if (k.ibound) {
      if (!k.npi) {
        Topo->Fill(0.0, w);
        if (!k.nfsiparts) {
          Topo->Fill(1.0, w);
        }
      } else {
        if (k.npi == 1) {
          if (!k.nfsiparts) {
            Topo->Fill(5.0, w);
          }
          if (k.npip) {
            Topo->Fill(2.0, w);
          } else if (k.npim) {
            Topo->Fill(3.0, w);
          } else if (k.npi0) {
            Topo->Fill(4.0, w);
          }
        } else {
          Topo->Fill(6.0, w);
        }
      }
    } else {
      if (!k.npi) {
        Topo_freep->Fill(0.0, w);
      } else {
        if (k.npi == 1) {
          if (k.npip) {
            Topo_freep->Fill(1.0, w);
          } else if (k.npim) {
            Topo_freep->Fill(2.0, w);
          } else if (k.npi0) {
            Topo_freep->Fill(3.0, w);
          }
        } else {
          Topo_freep->Fill(4.0, w);
        }
      }
    }
  }
  TH1D *Enu;
  TH1D *Q2;
  TH1D *EMiss;
  TH1D *W;
  TH2D *pthetalep;
  TH2D *pthetahmfspi;
  TH1D *Tpi;
  TH1D *Topo;
  TH1D *Topo_freep;
};

inline Kine GetKine(NeutVect *nvect) {

  TLorentzVector ISP4 = nvect->PartInfo(0)->fP;

  int NNuc = 0;
  NeutPart const &sn1 = *nvect->PartInfo(1);
  if ((!sn1.fIsAlive && (sn1.fStatus == -1)) &&
      ((sn1.fPID == 2212) || (sn1.fPID == 2112))) {
    ISP4 += sn1.fP;
    NNuc++;
  }
  NeutPart const &sn2 = *nvect->PartInfo(2);
  if ((!sn2.fIsAlive && (sn2.fStatus == -1)) &&
      ((sn2.fPID == 2212) || (sn2.fPID == 2112))) {
    ISP4 += sn2.fP;
    NNuc++;
  }

  Kine k;
  NeutPart const &probe = *nvect->PartInfo(0);
  k.Enu = probe.fP.E() * 1E-3;
  k.mode = nvect->Mode;
  k.ibound = nvect->Ibound;

  double MassRemnant_MeV = (nvect->TargetA - NNuc) * 938;
  double TRemnant_MeV =
      sqrt(pow(MassRemnant_MeV, 2) + ISP4.Vect().Mag2()) - MassRemnant_MeV;
  double EIS = probe.fP.E() + NNuc * 938;
  double EFS = TRemnant_MeV;

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
        if (k.foundlep) {
          k.nfsiparts++;
        }
        continue;
      }
    }
    // Calculate Q2 from the first outgoing lepton in the stack (makes Q2 for
    // MEC events correct)
    EFS += np.fP.E();
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
      k.npip += (np.fPID == 211);
      k.npim += (np.fPID == -211);
      k.npi0 += (np.fPID == 111);
    }
  }

  k.EMiss = EIS - EFS;

  return k;
}

#endif
