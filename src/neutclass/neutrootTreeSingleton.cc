#include "neutrootTreeSingleton.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <stdexcept>

NeutrootTreeSingleton *NeutrootTreeSingleton::fInstance = 0;
NeutrootTreeSingleton::NeutrootTreeSingleton() {}
NeutrootTreeSingleton::~NeutrootTreeSingleton() {
  delete tree_neutroot;
  delete nvect;
  fInstance = 0;
}

NeutrootTreeSingleton &
NeutrootTreeSingleton::Initialize(std::string const &filename) {
  if (fInstance == 0) {
    std::cout << "[INFO]: NeutrootTreeSingleton initialization with filename: "
              << filename << std::endl;
    fInstance = new NeutrootTreeSingleton;
    fInstance->LoadTree(filename);
  } else {
    throw std::runtime_error(
        std::string(
            "[ERROR]: NeutrootTreeSingleton already instantiated, cannot "
            "re-instantiate with filename: ") +
        filename);
  }

  return Get();
  ;
}

NeutrootTreeSingleton &NeutrootTreeSingleton::Get() {
  if (!fInstance) {
    throw std::runtime_error(
        "[ERROR]: Attempted to get NeutrootTreeSingleton instance before "
        "NeutrootTreeSingleton::Initialize was called.");
  }

  return *fInstance;
}

void NeutrootTreeSingleton::LoadTree(std::string const &filename) {

  tree_neutroot = new TChain("neuttree", "");
  f_nFiles = tree_neutroot->Add(Form("%s/neuttree", filename.c_str()));

  if (!f_nFiles) {
    std::cout << "[ERROR]: NeutrootTreeSingleton::LoadTree(string): No "
                 "files/neuttree "
                 "found in: "
              << filename << std::endl;
    throw std::runtime_error(
        std::string("[ERROR]: NeutrootTreeSingleton::LoadTree(string): No "
                    "files/neuttree found in: ") +
        filename);
  }

  f_nEvents = tree_neutroot->GetEntries();

  std::cout << "[INFO]: Files added: " << f_nFiles
            << ", with number of events: " << f_nEvents << std::endl;

  tree_neutroot->SetBranchStatus("*", false);
  tree_neutroot->SetBranchStatus("vectorbranch", true);

  br_neutvect = tree_neutroot->GetBranch("vectorbranch");
  nvect = new NeutVect();

  if (br_neutvect) {
    br_neutvect->SetAddress(&nvect);
  } else {
    throw std::runtime_error(
        "[ERROR]: NeutrootTreeSingleton::LoadTree(string) cannot find branch "
        "\"vectorbranch\". Are you using a neutroot generated file?");
  }
}

Int_t NeutrootTreeSingleton::GetEntry(Long64_t entry, Int_t getall) {

  Int_t nbytes = 0;

  if (!tree_neutroot) {
    throw std::runtime_error(
        "[ERROR] NeutrootTreeSingleton::GetEntry(): tree_neutroot is NULL");
  }

  if (entry != tree_neutroot->GetReadEntry()) {
    nbytes = tree_neutroot->GetEntry(entry, getall);
  }

  return nbytes;
}
Long64_t NeutrootTreeSingleton::GetEntries() {
  if (!tree_neutroot) {
    throw std::runtime_error(
        "[ERROR] NeutrootTreeSingleton::GetEntries(): tree_neutroot is NULL");
  }
  return tree_neutroot->GetEntries();
}

NeutVect *NeutrootTreeSingleton::GetNeutVectAddress() {

  if (!nvect) {
    throw std::runtime_error(
        "[ERROR] NeutrootTreeSingleton::GetNeutVectAddress(): Attempting to "
        "return NULL nvect");
  }

  return nvect;
}

double NeutrootTreeSingleton::GetXsecWeight() {
  TFile *cfile = tree_neutroot->GetCurrentFile();
  double weight = 1;
  if (cfile) {
    TH1 *evtrt, *fluxhisto;
    cfile->GetObject("ratehisto", evtrt);
    cfile->GetObject("fluxhisto", fluxhisto);
    if (evtrt && fluxhisto) {
      weight = 1E-38 * (evtrt->Integral() /
                        (fluxhisto->Integral() * double(GetEntries())));
    } else {
      Long64_t cent = tree_neutroot->GetReadEntry();
      NeutVect *nv = GetNeutVectAddress();
      GetEntry(0);
      double e = nv->PartInfo(0)->fP.E();
      bool ismono = true;
      for (int i = 1; i < 1000; ++i) {
        GetEntry(i);
        if(std::abs(e-nv->PartInfo(0)->fP.E())>1E-6){
          ismono = false;
          break;
        }
      }
      if(ismono){
        weight = nv->Totcrs*1E-38/double(GetEntries());
      }
      GetEntry(cent);
    }
  }
  return weight;
}
