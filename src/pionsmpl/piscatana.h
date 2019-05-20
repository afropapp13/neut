#ifndef piscatana_h
#define piscatana_h

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

#include "pars_piscat.h"
#include "../neutgeom/PeriodicTable_const.h"

extern "C"
{
  void nesetfgparams_();
}

class piscatana {
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
   Int_t           Nvert;
   Float_t         Posvert[150][3];   //[Nvert]
   Int_t           Iflgvert[150];   //[Nvert]
   Int_t           Nvcvert;
   Float_t         Dirvert[900][3];   //[Nvcvert]
   Float_t         Abspvert[900];   //[Nvcvert]
   Float_t         Abstpvert[900];   //[Nvcvert]
   Int_t           Ipvert[900];   //[Nvcvert]
   Int_t           Iverti[900];   //[Nvcvert]
   Int_t           Ivertf[900];   //[Nvcvert]
   Float_t         Fsiprob;
   Int_t           Numbndn;
   Int_t           Numbndp;
   Int_t           Numfrep;
   Int_t           Numatom;

   TBranch        *b_Npvc;   //!
   TBranch        *b_Ipvc;   //!
   TBranch        *b_Iorgvc;   //!
   TBranch        *b_Ichvc;   //!
   TBranch        *b_Iflvc;   //!
   TBranch        *b_Abspvc;   //!
   TBranch        *b_Pvc;   //!
   TBranch        *b_Posvc;   //!
   TBranch        *b_Nvert;   //!
   TBranch        *b_Posvert;   //!
   TBranch        *b_Iflgvert;   //!
   TBranch        *b_Nvcvert;   //!
   TBranch        *b_Dirvert;   //!
   TBranch        *b_Abspvert;   //!
   TBranch        *b_Abstpvert;   //!
   TBranch        *b_Ipvert;   //!
   TBranch        *b_Iverti;   //!
   TBranch        *b_Ivertf;   //!
   TBranch        *b_Fsiprob;   //!
   TBranch        *b_Numbndn;   //!
   TBranch        *b_Numbndp;   //!
   TBranch        *b_Numfrep;   //!
   TBranch        *b_Numatom;   //!

   piscatana(std::string inputfiles,std::string outputfile);
   piscatana(std::string inputfiles,std::string legNames, int ExpNumEvents);

   virtual ~piscatana();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);

   void InitFromTree();
   void InitFromHistos(std::vector<std::string> infiles_parsed, int ExpNumEvents=-1);
 
   virtual void     Loop();
   virtual void     Show(Long64_t entry = -1);
   
   // MC histograms
   void     Init_histo();
   std::vector<TH1F*> piInitSpec[nPions][kZmax];
   std::vector<TH1F*> piInitIntSpec[nPions][kZmax][nInts];
   std::vector<TH1F*> piInitPiprodSpec[nPions][kZmax][nPiprods];
   std::vector<TH1F*> piDiffXsec[nPions][kZmax][3];
   std::vector<TH2F*> outPhaseSpc[nPions][kZmax][3];

   bool *histosLoaded[nPions][kZmax];

   // Data graphs
   void LoadData();
   bool FileExists(std::string filename);
   TGraphErrors *int_graph[nPions][kZmax][nInts];

   TFile *outROOT;

   int nConfs;
   std::vector<std::string> legName_parse;
   bool DrawLegendConfs(std::vector<TH1F*> h_input,Double_t x0=0.5, Double_t y0=0.5, Double_t x1=0.95, Double_t y1=0.95) ;
   void DrawLegendInts(std::vector<TH1F*> *h_input,Double_t x0=0.5, Double_t y0=0.5, Double_t x1=0.95, Double_t y1=0.95, int iNuc=0) ;

   void plot_sep_nuc();
   void plot_sep_pion_nuc_iint();
   void plot_sep_pion_iint();
   void print_canvas(TCanvas *c, std::string outfilename="");

 private:
   
   inline void int_xsec();

   int init_pion;
   
   // For normalizing
   Long64_t nentries;
   
   double targetArea;
   double dtheta;

   void SetInitPion();
   void SetNucleus();
   void fill_necard();
   
   // File names
   std::string theFile,outFile;

};

#endif
