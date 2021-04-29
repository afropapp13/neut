// NEUT ReWeight
#include <iostream>
#include <sstream>
#include <iomanip>

#include "NReWeight.h"
#include "NuXSecRES.h"
#include "NReWeightEngineI.h"
#include "NSyst.h"

// Common block variables for NEUT
#include "nefillverC.h"
#include "necardC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"
#include "neworkC.h"
#include "fsihistC.h"
#include "neutcrsC.h"

// NEUT class headers
#include "neutvect.h"
#include "neutpart.h"
#include "neutfsipart.h"
#include "neutfsivert.h"
#include "neutrootTreeSingleton.h"
#include "necardC.h"

// Common block interface
#include "CommonBlockIFace.h"

// ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

// Get Cos theta with Adler angles
double CosThAdler(TLorentzVector Pnu, TLorentzVector Pmu,
                            TLorentzVector Ppi, TLorentzVector Pprot) {
  // Get the "resonance" lorentz vector (pion proton system)
  TLorentzVector Pres = Pprot + Ppi;
  // Boost the particles into the resonance rest frame so we can define the
  // x,y,z axis
  Pnu.Boost(-Pres.BoostVector());
  Pmu.Boost(-Pres.BoostVector());
  Ppi.Boost(-Pres.BoostVector());

  // The z-axis is defined as the axis of three-momentum transfer, \vec{k}
  // Also unit normalise the axis
  TVector3 zAxis =
      (Pnu.Vect() - Pmu.Vect()) * (1.0 / ((Pnu.Vect() - Pmu.Vect()).Mag()));

  // Then the angle between the pion in the "resonance" rest-frame and the
  // z-axis is the theta Adler angle
  double costhAdler = cos(Ppi.Vect().Angle(zAxis));

  return costhAdler;
}

// Get phi with Adler angles, a bit more complicated...
double PhiAdler(TLorentzVector Pnu, TLorentzVector Pmu,
                          TLorentzVector Ppi, TLorentzVector Pprot) {
  // Get the "resonance" lorentz vector (pion proton system)
  TLorentzVector Pres = Pprot + Ppi;
  // Boost the particles into the resonance rest frame so we can define the
  // x,y,z axis
  Pnu.Boost(-Pres.BoostVector());
  Pmu.Boost(-Pres.BoostVector());
  Ppi.Boost(-Pres.BoostVector());

  // The z-axis is defined as the axis of three-momentum transfer, \vec{k}
  // Also unit normalise the axis
  TVector3 zAxis =
      (Pnu.Vect() - Pmu.Vect()) * (1.0 / ((Pnu.Vect() - Pmu.Vect()).Mag()));

  // The y-axis is then defined perpendicular to z and muon momentum in the
  // resonance frame
  TVector3 yAxis = Pnu.Vect().Cross(Pmu.Vect());
  yAxis *= 1.0 / double(yAxis.Mag());

  // And the x-axis is then simply perpendicular to z and x
  TVector3 xAxis = yAxis.Cross(zAxis);
  xAxis *= 1.0 / double(xAxis.Mag());

  // Project the pion on to x and y axes
  double x = Ppi.Vect().Dot(xAxis);
  double y = Ppi.Vect().Dot(yAxis);

  double newphi = atan2(y, x) * (180. / M_PI);
  // Convert negative angles to positive
  if (newphi < 0.0)
    newphi += 360.0;

  return newphi;
}


int main(int argc, char *argv[]) {

  std::cout << std::setprecision(9);

  if (argc != 2) {
    std::cerr << "Syntax is: " << argv[0] << " OUTPUT_NEUT_ROOT_FILE.root" << std::endl;
    return 1;
  }

  std::string InputName = std::string(argv[1]);

  // Get the common block interface
  std::string cardfile = "";
  if (InputName.find("_deltaflat") != std::string::npos) cardfile = "/vols/t2k/users/cvw09/software/NEUT/submit_neut_2020_dev/t2k_rspiej_nuc_nubar_deltaflat.card";
  else if (InputName.find("_all") != std::string::npos) cardfile = "/vols/t2k/users/cvw09/software/NEUT/submit_neut_2020_dev/t2k_rspiej_nuc_nubar_all.card";
  else if (InputName.find("_delta") != std::string::npos) cardfile = "/vols/t2k/users/cvw09/software/NEUT/submit_neut_2020_dev/t2k_rspiej_nuc_nubar_delta.card";
  else if (InputName.find("_iso") != std::string::npos) cardfile = "/vols/t2k/users/cvw09/software/NEUT/submit_neut_2020_dev/t2k_rspiej_nuc_nubar_iso.card";
  //else if (InputName.find("_delta") != std::string::npos) cardfile = "/vols/t2k/users/cvw09/software/NEUT/neut_2020_dev/submit/t2k_rspiej_nuc_nubar_delta.card";
  //else if (InputName.find("_iso") != std::string::npos) cardfile = "/vols/t2k/users/cvw09/software/NEUT/neut_2020_dev/submit/t2k_rspiej_nuc_nubar_iso.card";
  //else if (InputName.find("_all") != std::string::npos) cardfile = "/vols/t2k/users/cvw09/software/NEUT/neut_2020_dev/submit/t2k_rspiej_nuc_nubar_all.card";
  else {
    std::cerr << "Couldn't match filename to cardfile, returning" << std::endl;
    return -1;
  }
  std::cout << "Setting cardfile " << cardfile << " for " << InputName << std::endl;

  neut::CommonBlockIFace::Initialize(cardfile);
  neut::CommonBlockIFace const &cbfa = neut::CommonBlockIFace::Get();

  neut::rew::NReWeight rw;
  rw.AdoptWeightEngine("NuXSecRES", std::unique_ptr<neut::rew::NReWeightEngineI>(new neut::rew::NuXSecRESEngine));

  const int nvars = 4;
  double vals[nvars];
  for (int i = 0; i < nvars; ++i) {
    vals[i] = i;
  }

  // The weights for each variation
  double weights[nvars];

  // Setup the output variables
  double Q2, EnuTrue, W, ppi, Eout, pinit, pn, costhlep, costhpinu, costhpimu, costhpin, ScaleFactor, ppi_prefsi, pnuc_prefsi, costhad, phiad = -999;
  int mode, leppdg, nupdg, pipdg, nucpdg = -999;
  int npart, npi, nnuc = 0;
  bool isbound = false;
  bool pblocked = false;

  std::stringstream ss;
  ss << "var_rew";
  TFile *OutputFile = new TFile((InputName+"_to_MARES_"+ss.str()+".root").c_str(), "RECREATE");
  OutputFile->cd();
  TTree *OutputTree = new TTree("VARS", "VARS");
  // Now set up the branches
  OutputTree->Branch("EnuTrue", &EnuTrue, "EnuTrue/D");
  OutputTree->Branch("Eout", &Eout, "Eout/D");
  OutputTree->Branch("Q2", &Q2, "Q2/D");
  OutputTree->Branch("W", &W, "W/D");
  OutputTree->Branch("ppi", &ppi, "ppi/D");
  OutputTree->Branch("ppi_prefsi", &ppi_prefsi, "ppi_prefsi/D");
  OutputTree->Branch("pnuc_prefsi", &pnuc_prefsi, "pnuc_prefsi/D");
  OutputTree->Branch("pn", &pn, "pn/D");
  OutputTree->Branch("pinit", &pinit, "pinit/D");
  OutputTree->Branch("costhlep", &costhlep, "costhlep/D");
  OutputTree->Branch("costhpinu", &costhpinu, "costhpinu/D");
  OutputTree->Branch("costhpimu", &costhpimu, "costhpimu/D");
  OutputTree->Branch("costhpin", &costhpin, "costhpin/D");

  OutputTree->Branch("costhad", &costhad, "costhad/D");
  OutputTree->Branch("phiad", &phiad, "phiad/D");

  OutputTree->Branch("mode", &mode, "mode/I");
  OutputTree->Branch("nupdg", &nupdg, "nupdg/I");
  OutputTree->Branch("leppdg", &leppdg, "leppdg/I");
  OutputTree->Branch("pipdg", &pipdg, "pipdg/I");
  OutputTree->Branch("nucpdg", &nucpdg, "nucpdg/I");

  OutputTree->Branch("npart", &npart, "npart/I");
  OutputTree->Branch("npi", &npi, "npi/I");
  OutputTree->Branch("nnuc", &nnuc, "nnuc/I");

  OutputTree->Branch("weights", weights, "weights[4]/D");
  OutputTree->Branch("scalefactor", &ScaleFactor, "ScaleFactor/D");

  OutputTree->Branch("isbound", &isbound, "isbound/O");
  OutputTree->Branch("pblocked", &pblocked, "pblocked/O");

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
  //nEvents = 500000;
  // Print width
  int PrintWidth = nEvents/20;

  // Get the scaling factor
  TH1D* Flux = (TH1D*)(Infile->Get("flux_numub")->Clone());
  TH1D* Eventrate = (TH1D*)(Infile->Get("evtrt_numub")->Clone());
  ScaleFactor = ((Eventrate->Integral("width")*1.E-38)/(nEvents))/Flux->Integral("width");

  // And each event
  for (int i = 0; i < nEvents; ++i) {

    //std::cout << "Event " << j << std::endl;

    // Reset variables
    Q2 = -999;
    EnuTrue = -999;
    W = -999;
    ppi = -999;
    ppi_prefsi = -999;
    pnuc_prefsi = -999;
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
    isbound = false;
    pblocked = false;
    costhad = -999;
    phiad = -999;
    for (int j = 0; j < nvars; ++j) weights[j] = -999;

    npart = 0;
    npi = 0;
    nnuc = 0;

    if (i % PrintWidth == 0) {
      std::cout << "On event " << i << "/" << nEvents << " (" << int(double(i)/double(nEvents)*100.0) << "%)" << std::endl;
    }

    Tree->GetEntry(i);
    // Need to fill nework first
    //neut::CommonBlockIFace::ReadVect(nvect);
    cbfa.ReadVect(nvect);

    isbound = nvect->Ibound;

    //std::cout << "***" << std::endl;
    //std::cout << "yfn saved: " << nvect->yfn << std::endl;


    //if (nvect->Mode > 13) continue;

    //std::cout << std::endl;
    //std::cout << "***" << std::endl;
    for (int j = 0; j < nvars; ++j) {
      //std::cout << "*" << std::endl;
      //std::cout << "Variation " << j << std::endl;
      //std::cout << "Setting to " << vals[j] << std::endl;
      rw.SetDial_To_Value(rw.DialFromString("MDLSPiEj"),vals[j]);
      //rw.Systematics().Set(neut::rew::kXSecTwkDial_MaRES, vals[j]);
      //rw.Systematics().Set(neut::rew::kXSecTwkDial_CA5RES, vals[j]);
      // Reconfigure the reweight engine
      rw.Reconfigure(); 
      weights[j] = rw.CalcWeight();
      if (std::isnan(weights[j])) {
        std::cout << "Found nan weight, setting to 1!" << std::endl;
        weights[j] = 1;
      }
      if (weights[j] < 0) {
        std::cout << "Negative weight!" << std::endl;
      }
      //std::cout << "MDLSPIEJ: " << vals[j] << std::endl;
      //double oldx = dynamic_cast<neut::rew::NReWeightNuXSecRES*>(rw.WghtCalc("res"))->GetOld();
      //double newx = dynamic_cast<neut::rew::NReWeightNuXSecRES*>(rw.WghtCalc("res"))->GetNew();
      //std::cout << "yfn old: " << oldx << std::endl;
      //std::cout << "yfn new: " << newx << std::endl;
      //std::cout << "ratio = " << newx/oldx << std::endl;
      //std::cout << "weight = " << weights[j] << std::endl;
    }

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
    pinit = Pinit.Vect().Mag();

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

    if (nvect->PartInfo(3)->fStatus != 5) {
      TLorentzVector prefsi_nuc = nvect->PartInfo(3)->fP;
      pnuc_prefsi = prefsi_nuc.Vect().Mag()/1.E3;
      pblocked = false;
    } else {
      pblocked = true;
    }

    if (nvect->PartInfo(3)->fStatus != 5) {
      TLorentzVector prefsi_pi = nvect->PartInfo(4)->fP;
      ppi_prefsi = prefsi_pi.Vect().Mag()/1.E3;
    }

    if (npi > 0) {
      ppi = Ppi.Vect().Mag()/1.E3; 
      costhpinu = cos(Ppi.Vect().Angle(Pnu.Vect()));
      costhpimu = cos(Ppi.Vect().Angle(Pout.Vect()));
    }

    if (nnuc > 0) {
      pn = Pn.Vect().Mag()/1.E3;
      costhpin = cos(Ppi.Vect().Angle(Pn.Vect()));
    }

    costhad = CosThAdler(nvect->PartInfo(0)->fP, nvect->PartInfo(2)->fP, nvect->PartInfo(4)->fP, nvect->PartInfo(3)->fP);
    phiad = PhiAdler(nvect->PartInfo(0)->fP, nvect->PartInfo(2)->fP, nvect->PartInfo(4)->fP, nvect->PartInfo(3)->fP);

    OutputTree->Fill();
  } // Finished looping over events

  std::cout << "Wrote to " << OutputFile->GetName() << std::endl;
  OutputFile->cd();
  OutputTree->Write();

  std::cout << "Finished looping events" << std::endl;

  return 0;
}
