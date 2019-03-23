#include <iostream>

#include "neutrootTreeSingleton.h"

using namespace std;
using std::cout;
using std::endl;

NeutrootTreeSingleton * NeutrootTreeSingleton::fInstance = 0;
//____________________________________________________________________________
NeutrootTreeSingleton::NeutrootTreeSingleton()
{
}
//____________________________________________________________________________
NeutrootTreeSingleton::~NeutrootTreeSingleton()
{
  fInstance = 0;
}
//____________________________________________________________________________

NeutrootTreeSingleton * NeutrootTreeSingleton::Instance(string filename)
{
  if(fInstance == 0) {
    cout << "NeutrootTreeSingleton initialization with filename: " << filename.c_str() << endl;
    static NeutrootTreeSingleton::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new NeutrootTreeSingleton;
    fInstance->LoadTree(filename);
  }

  return fInstance;
}
//____________________________________________________________________________

NeutrootTreeSingleton * NeutrootTreeSingleton::Instance(TTree *a_intree)
{
  if(fInstance == 0) {

    cout << "NeutrootTreeSingleton::Instance(TTree) Error: Do not initialize with pre-loaded tree. Pass a filename instead." << endl;
    exit (-1);

    //cout << "NeutrootTreeSingleton late initialization" << endl;
    //static NeutrootTreeSingleton::Cleaner cleaner;
    //cleaner.DummyMethodAndSilentCompiler();
    //fInstance = new NeutrootTreeSingleton;
    //fInstance->LoadTree(a_intree);
  }

  return fInstance;
}
//____________________________________________________________________________
void NeutrootTreeSingleton::LoadTree(string filename) {

  tree_neutroot = new TChain("neuttree","");
  f_nFiles = tree_neutroot->Add(Form("%s/neuttree",filename.c_str()));
      
  if (!f_nFiles) {
    cout << "Error NeutrootTreeSingleton::LoadTree(string): No files/neuttree found in: " << filename.c_str() << endl;
    exit (-1);
  }
      
  f_nEvents = tree_neutroot->GetEntries();

  cout << "Files added: " << f_nFiles << ", with number of events: " << f_nEvents << endl;
    
  br_neutvect = tree_neutroot->GetBranch("vectorbranch");
  nvect = NULL;
    
  if(br_neutvect){
    br_neutvect->SetAddress(&nvect);
      
  } else {
    cout << "Error: NeutrootTreeSingleton::LoadTree(string) cannot find branch \"vectorbranch\". Are you using a neutroot generated file?" << endl;
    exit (-1);
  }

}
//____________________________________________________________________________
void NeutrootTreeSingleton::LoadTree(TTree *a_intree) {

  //tree_neutroot = a_intree;
  //
  //if (!tree_neutroot) {
  //  cout << "Error NeutrootTreeSingleton::LoadTree(TTree*): tree_neutroot is NULL" << endl;
  //  exit (-1);
  //}
  //
  //f_nFiles = tree_neutroot->GetNtrees();
  //    
  //if (!f_nFiles) {
  //  cout << "Error NeutrootTreeSingleton::LoadTree(TTree*): No files/neuttree passed" << endl;
  //  exit (-1);
  //}
  //    
  //f_nEvents = tree_neutroot->GetEntries();
  //
  //cout << "Files added: " << f_nFiles << ", with number of events: " << f_nEvents << endl;
  //  
  //br_neutvect = tree_neutroot->GetBranch("vectorbranch");
  //nvect = NULL;
  //  
  //if(br_neutvect){
  //  br_neutvect->SetAddress(&nvect);
  //    
  //} else {
  //  cout << "Error NeutrootTreeSingleton::LoadTree(TTree*): cannot find branch \"vectorbranch\". Are you using a neutroot generated file?" << endl;
  //  exit (-1);
  //}

}
//____________________________________________________________________________

Int_t NeutrootTreeSingleton::GetEntry(Long64_t entry, Int_t getall) {

  Int_t nbytes = 0;
  
  if (!tree_neutroot) {
    cout << "Error NeutrootTreeSingleton::GetEntry(): tree_neutroot is NULL" << endl;
    exit (-1);
  }
  
  if (entry != tree_neutroot->GetReadEntry()) {
    nbytes = tree_neutroot->GetEntry(entry, getall);
  }
  
  return nbytes;
}
//____________________________________________________________________________

NeutVect* NeutrootTreeSingleton::GetNeutVectAddress() {

  if (!nvect) {
    cout << "Error NeutrootTreeSingleton::GetNeutVectAddress(): Attempting to return NULL nvect" << endl;
    exit (-1);
  }
  
  return nvect;
}

ClassImp(NeutrootTreeSingleton)
