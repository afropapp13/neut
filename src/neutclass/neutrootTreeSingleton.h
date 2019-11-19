
#ifndef _N_NEUTROOTTREESINGLETON_H_
#define _N_NEUTROOTTREESINGLETON_H_

#include "TChain.h"
#include "TTree.h"
#include "neutvect.h"

class NeutrootTreeSingleton {

public:
  static NeutrootTreeSingleton &Initialize(std::string const &filename);
  static bool Destroy() { delete fInstance; }
  static NeutrootTreeSingleton &Get();

  NeutVect *GetNeutVectAddress();
  Int_t GetEntry(Long64_t entry = 0, Int_t getall = 0);
  Long64_t GetEntries();

private:
  NeutrootTreeSingleton();
  ~NeutrootTreeSingleton();

  void LoadTree(std::string const &filename);

  TChain *tree_neutroot;
  TBranch *br_neutvect;
  NeutVect *nvect;

  int f_nEvents;
  int f_nFiles;

  static NeutrootTreeSingleton *fInstance;
};

#endif
