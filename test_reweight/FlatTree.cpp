// NEUT ReWeight
#include <iostream>
#include <sstream>

// ROOT 
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

// NEUT class headers
#include "neutvect.h"
#include "neutpart.h"
#include "neutfsipart.h"
#include "neutfsivert.h"
#include "neutrootTreeSingleton.h"

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cerr << "Syntax is: ./tester OUTPUT_NEUT_ROOT_FILE.root" << std::endl;
    return 1;
  }

  std::string InputName = std::string(argv[1]);

  // Make some 1D distributions for each mode
  const int nNEUTmodes = 51;
  TH1D *Weighted[nNEUTmodes];
  TH1D *UnWeighted[nNEUTmodes];
  for (int i = 0; i < nNEUTmodes; ++i) {
    std::stringstream ss;
    ss << "_" << i;
    Weighted[i] = new TH1D((std::string("Weighted")+ss.str()).c_str(), (std::string("Weighted")+ss.str()+";p_{#mu} (GeV);d#sigma/dp_{#mu} (cm^2/nucleon/GeV)").c_str(), 1, 0, 5);
    UnWeighted[i] = new TH1D((std::string("UnWeighted")+ss.str()).c_str(), (std::string("UnWeighted")+ss.str()+";p_{#mu} (GeV);d#sigma/dp_{#mu} (cm^2/nucleon/GeV)").c_str(), 1, 0, 5);
  }

  // Setup the output variables
  double Q2, EnuTrue, W, ppi, Eout, pinit, pn, costhlep, costhpinu, costhpimu, costhpin, weight, ScaleFactor = -999;
  int mode, leppdg, nupdg, pipdg, nucpdg = -999;
  int npart, npi, nnuc, nbad = 0;

  std::stringstream ss;
  ss << "var";
  TFile *OutputFile = new TFile((InputName+"_to_"+ss.str()+".root").c_str(), "RECREATE");
  OutputFile->cd();
  TTree *OutputTree = new TTree("VARS", "VARS");
  // Now set up the branches
  OutputTree->Branch("EnuTrue", &EnuTrue, "EnuTrue/D");
  OutputTree->Branch("Eout", &Eout, "Eout/D");
  OutputTree->Branch("Q2", &Q2, "Q2/D");
  OutputTree->Branch("W", &W, "W/D");
  OutputTree->Branch("ppi", &ppi, "ppi/D");
  OutputTree->Branch("pn", &pn, "pn/D");
  OutputTree->Branch("pinit", &pinit, "pinit/D");
  OutputTree->Branch("costhlep", &costhlep, "costhlep/D");
  OutputTree->Branch("costhpinu", &costhpinu, "costhpinu/D");
  OutputTree->Branch("costhpimu", &costhpimu, "costhpimu/D");
  OutputTree->Branch("costhpin", &costhpin, "costhpin/D");

  OutputTree->Branch("mode", &mode, "mode/I");
  OutputTree->Branch("nupdg", &nupdg, "nupdg/I");
  OutputTree->Branch("leppdg", &leppdg, "leppdg/I");
  OutputTree->Branch("pipdg", &pipdg, "pipdg/I");
  OutputTree->Branch("nucpdg", &nucpdg, "nucpdg/I");

  OutputTree->Branch("npart", &npart, "npart/I");
  OutputTree->Branch("npi", &npi, "npi/I");
  OutputTree->Branch("nnuc", &nnuc, "nnuc/I");
  OutputTree->Branch("nbad", &nbad, "nbad/I");

  OutputTree->Branch("weight", &weight, "weight/D");
  OutputTree->Branch("scalefactor", &ScaleFactor, "ScaleFactor/D");

  // Open the input
  TFile *Infile = new TFile(InputName.c_str(), "OPEN");
  if(!Infile){
    std::cerr << "Cannot open neutroot file!" << std::endl;
    throw;
  }

  // Open a TTree
  TTree * Tree = (TTree*) Infile->Get("neuttree");
  if(!Tree){
    std::cerr << "Cannot find neuttree!" << std::endl;
    throw;
  }

  NeutVect *nvect = NULL;
  Infile->cd();
  Tree->SetBranchAddress("vectorbranch", &nvect);
  int nEvents = Tree->GetEntriesFast();
  // Print width
  int PrintWidth = nEvents/20;

  // Get the scaling factor
  TH1D* Flux = (TH1D*)(Infile->Get("flux_numu")->Clone());
  TH1D* Eventrate = (TH1D*)(Infile->Get("evtrt_numu")->Clone());
  ScaleFactor = ((Eventrate->Integral("width")*1.E-38)/(nEvents))/Flux->Integral("width");

  // And each event
  for (int j = 0; j < nEvents; ++j) {

    //std::cout << "Event " << j << std::endl;

    // Reset variables
    Q2 = -999;
    EnuTrue = -999;
    W = -999;
    ppi = -999;
    Eout = -999;
    pn = -999;
    pinit = -999;
    costhlep = -999;
    costhpinu = -999;
    costhpimu = -999;
    costhpin = -999;
    mode = -999;
    leppdg = -999;
    nupdg = -999;
    pipdg = -999;
    nucpdg = -999;
    weight = -999;

    npart = 0;
    npi = 0;
    nnuc = 0;
    nbad = 0;

    if (j % PrintWidth == 0) {
      std::cout << "On event " << j << "/" << nEvents << " (" << int(double(j)/double(nEvents)*100.0) << "%)" << std::endl;
    }

    Tree->GetEntry(j);
    // Need to fill nework first
    //rw.FillNeutCommons(nvect);
    //rw.PrintNeutAll();
    //weight = rw.CalcWeight();

    // Count how many bad events
    //if (weight == 1.0) nbad++;

    // Now get some distributions
    TLorentzVector Pnu = nvect->PartInfo(0)->fP;
    TLorentzVector Pinit = nvect->PartInfo(1)->fP;
    TLorentzVector Pout = nvect->PartInfo(2)->fP;
    TLorentzVector Ppi;
    TLorentzVector Pn;

    npart = nvect->Npart();
    EnuTrue = Pnu.E()/1.E3;
    Eout = Pout.E()/1.E3;
    costhlep = cos(Pnu.Vect().Angle(Pout.Vect()));
    mode = nvect->Mode;
    nupdg = nvect->PartInfo(0)->fPID;
    leppdg = nvect->PartInfo(2)->fPID;
    Q2 = -1*(Pnu-Pout)*(Pnu-Pout)/1.E6;
    W = sqrt((Pnu+Pinit-Pout)*(Pnu+Pinit-Pout))/1.E3;
    pinit = Pinit.Vect().Mag()/1.E3;

    UnWeighted[mode]->Fill(Eout);
    Weighted[mode]->Fill(Eout, weight);

    // Loop over events
    for (int k = 2; k < nvect->Npart(); ++k) {
      int PID = nvect->PartInfo(k)->fPID;
      if (!(nvect->PartInfo(k))->fIsAlive && (nvect->PartInfo(k))->fStatus != 0) continue;
      npart++;
      // Pion
      if (abs(PID) == 211 || PID == 111) {
        npi++;
        // Check highest energy and set
        if (nvect->PartInfo(k)->fP.E() > Ppi.E()) {
          Ppi = nvect->PartInfo(k)->fP;
          pipdg = nvect->PartInfo(k)->fPID;
        }
      } else if (PID == 2112 || PID == 2212) {
        nnuc++;
        if (nvect->PartInfo(k)->fP.E() > Pn.E()) {
          Pn = nvect->PartInfo(k)->fP;
          nucpdg = nvect->PartInfo(k)->fPID;
        }
      }
    } // Finished scanning for pions and nucleons

    if (npi > 0) {
      ppi = Ppi.Vect().Mag()/1.E3; 
      costhpinu = cos(Ppi.Vect().Angle(Pnu.Vect()));
      costhpimu = cos(Ppi.Vect().Angle(Pout.Vect()));
    }

    if (nnuc > 0) {
      pn = Pn.Vect().Mag()/1.E3;
      costhpin = cos(Ppi.Vect().Angle(Pn.Vect()));
    }

    OutputTree->Fill();
  } // Finished looping over events

  OutputFile->cd();
  OutputTree->Write();
  for (int j = 0; j < nNEUTmodes; ++j) {
    // Set Poisson errors
    for (int i = 0; i < Weighted[j]->GetNbinsX()+1; ++i) {
      double content1 = Weighted[j]->GetBinContent(i+1);
      double error1 = sqrt(content1);
      Weighted[j]->SetBinError(i+1, error1);
      double content2 = UnWeighted[j]->GetBinContent(i+1);
      double error2 = sqrt(content2);
      UnWeighted[j]->SetBinError(i+1, error2);
    }
    Weighted[j]->Scale(ScaleFactor, "width");
    UnWeighted[j]->Scale(ScaleFactor, "width");
    if (Weighted[j]->Integral() > 0) Weighted[j]->Write();
    if (UnWeighted[j]->Integral() > 0) UnWeighted[j]->Write();
  }

  std::cout << "Finished looping events" << std::endl;

  return 0;
}
