#ifndef nuclscatana_h
#define nuclscatana_h

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

#include "pars_nuclscat.h"
#include "../neutgeom/PeriodicTable_const.h"

extern "C"
{
  void nesetfgparams_();
}

class nuclscatana {
public :

   TChain          *fChain;   //!pointer to the analyzed TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Npvc;
   Int_t           Ipvc[100];   //[Npvc]
   Int_t           Iorgvc[100];   //[Npvc]
   Int_t           Ichvc[100];   //[Npvc]
   Int_t           Iflvc[100];   //[Npvc]
   Float_t         Abspvc[100];   //[Npvc]
   Float_t         Pvc[100][3];   //[Npvc]
   Float_t         Posvc[3];
   Int_t           Nfnvert;
   Int_t           Nfiflag[200];   //[Nfnvert]
   Float_t         X[200];   //[Nfnvert]
   Float_t         Y[200];   //[Nfnvert]
   Float_t         Z[200];   //[Nfnvert]
   Float_t         Px[200];   //[Nfnvert]
   Float_t         Py[200];   //[Nfnvert]
   Float_t         Pz[200];   //[Nfnvert]
   Float_t         E[200];   //[Nfnvert]
   Int_t           Nffirststep[200];   //[Nfnvert]
   Int_t           Nfnstep;
   Float_t         Ecms2[2000];   //[Nfnstep]
   Float_t         Prob[2000];   //[Nfnstep]
   Int_t           Numbndn;
   Int_t           Numbndp;
   Int_t           Numfrep;
   Int_t           Numatom;

   // List of branches
   TBranch        *b_Npvc;   //!
   TBranch        *b_Ipvc;   //!
   TBranch        *b_Iorgvc;   //!
   TBranch        *b_Ichvc;   //!
   TBranch        *b_Iflvc;   //!
   TBranch        *b_Abspvc;   //!
   TBranch        *b_Pvc;   //!
   TBranch        *b_Posvc;   //!
   TBranch        *b_Nfnvert;   //!
   TBranch        *b_Nfiflag;   //!
   TBranch        *b_X;   //!
   TBranch        *b_Y;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_Py;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_Nffirststep;   //!
   TBranch        *b_Nfnstep;   //!
   TBranch        *b_Ecms2;   //!
   TBranch        *b_Prob;   //!
   TBranch        *b_Numbndn;   //!
   TBranch        *b_Numbndp;   //!
   TBranch        *b_Numfrep;   //!
   TBranch        *b_Numatom;   //!

   nuclscatana(std::string inputfiles,std::string outputfile);
   nuclscatana(std::string inputfiles,std::string legNames, int ExpNumEvents);

   virtual ~nuclscatana();
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
