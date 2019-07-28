#ifndef nucltransana_h
#define nucltransana_h

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

#include "pars_nucltrans.h"
#include "../neutgeom/PeriodicTable_const.h"

extern "C"
{
  void nesetfgparams_();
}

class nucltransana {
public :

   TChain          *fChain;   //!pointer to the analyzed TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           vcwork_nvc;
   Float_t         vcwork_posvc[3];
   Int_t           vcwork_ipvc[100];
   Float_t         vcwork_amasvc[100];
   Float_t         vcwork_pvc[100][3];
   Int_t           vcwork_iorgvc[100];
   Int_t           vcwork_iflgvc[100];
   Int_t           vcwork_icrnvc[100];
   Float_t         vcwork_timvc[100];
   Float_t         vcwork_posivc[100][3];
   Int_t           vcwork_ivtivc[100];
   Float_t         vcwork_posfvc[100][3];
   Int_t           vcwork_ivtfvc[100];
   Int_t           nework_modene;
   Int_t           nework_numne;
   Int_t           nework_ipne[100];
   Float_t         nework_pne[100][3];
   Int_t           nework_iorgne[100];
   Int_t           nework_iflgne[100];
   Int_t           nework_icrnne[100];
   Int_t           posinnuc_ibound;
   Float_t         posinnuc_posnuc[100][3];
   Int_t           target_numbndn;
   Int_t           target_numbndp;
   Int_t           target_numfrep;
   Int_t           target_numatom;

   // List of branches
   TBranch        *b_vcwork;   //!
   TBranch        *b_nework;   //!
   TBranch        *b_posinnuc;   //!
   TBranch        *b_target;   //!

   nucltransana(std::string inputfiles,std::string outputfile);
   nucltransana(std::string inputfiles,std::string legNames, int ExpNumEvents);

   virtual ~nucltransana();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);

   void InitFromTree();
   void SetInitNucleon();
   void SetNucleus();

   void InitFromHistos(std::vector<std::string> infiles_parsed, int ExpNumEvents=-1);
 
   virtual void     Loop();
   virtual void     Show(Long64_t entry = -1);
   
   // MC histograms
   void     Init_histo();
   std::vector<TH1F*> piInitSpec[nNucs][kZmax];
   std::vector<TH1F*> piInitIntSpec[nNucs][kZmax][nInts];

   bool *histosLoaded[nNucs][kZmax];

   TFile *outROOT;

   int nConfs;
   std::vector<std::string> legName_parse;

   Float_t abspvc;

 private:
   
   inline void int_xsec();

   int init_nuc;
   
   // For normalizing
   Long64_t nentries;
   
   double targetArea;

   // File names
   std::string theFile,outFile;

};

#endif
