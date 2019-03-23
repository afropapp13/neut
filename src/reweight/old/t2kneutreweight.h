#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

const Int_t kMaxVtx = 16; // Comes from ND280 oaAnalysis

// Number of variable NEUT parameters
// Must change when adding more variable paramters
const Int_t nPars = 14; 
enum parNames {QEMA, QEMV, QEKAPP, SPIMA, SPIMV, COHMA, DISPDF, DISBOD, FSIQE, FSIQEH, FSIINEL, FSIABS, FSICX, FSICXH};

extern "C"
{
  float evdifcrs_();
  float evpiprob_();
  void  nesetfgparams_();
  void  nefillmodel_();
}


class t2kneutreweight {
 private :
   TChain         *fNeutChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent;     //!current Tree number in a TChain

   //_______________________________________________________
   //For standard vector file: SK (floats), old-INGRID(?), ND280-numc (doubles)
   Int_t           Neutmode; 
   TString         fString;    // ND280-numc MCP1
 
   Int_t           Npvc;
   Int_t           Ipvc[100];   //[Npvc]
   Float_t         Pvc[100][3];   //[Npvc]
   Float_t         Abspvc[100];   //[Npvc]

   Float_t         StdHepP4[100][4];  // ND280-numc MCP1
   Int_t           StdHepPdg[100];   // ND280-numc MCP1 [Npvc]

   Int_t           Nvert;
   Float_t         Posvert[100][3];   //[Nvert]
   Int_t           Iflgvert[100];   //[Nvert]
   Int_t           Nvcvert;
   Float_t         Dirvert[300][3];   //[Nvcvert]
   Float_t         Abspvert[300];   //[Nvcvert]
   Float_t         Abstpvert[300];   //[Nvcvert]
   Int_t           Ipvert[300];   //[Nvcvert]
   Int_t           Iverti[300];   //[Nvcvert]
   Int_t           Ivertf[300];   //[Nvcvert]


   //_______________________________________________________
   //For ND280-oaAnalysis file
   Int_t           NVtx;

   Int_t           Vtx_NEneutmode[kMaxVtx];   //[Vtx_]
   Int_t           Vtx_NEnvc[kMaxVtx];   //[Vtx_]
   Int_t           Vtx_NEipvc[kMaxVtx][100];   //[Vtx_]
   Float_t         Vtx_NEpvc[kMaxVtx][100][3];   //[Vtx_]

   Int_t           Vtx_NEnvert[kMaxVtx];   //[Vtx_]
   Float_t         Vtx_NEposvert[kMaxVtx][50][3];   //[Vtx_]
   Int_t           Vtx_NEiflgvert[kMaxVtx][50];   //[Vtx_]
   Int_t           Vtx_NEnvcvert[kMaxVtx];   //[Vtx_]
   Float_t         Vtx_NEdirvert[kMaxVtx][300][3];   //[Vtx_]
   Float_t         Vtx_NEabspvert[kMaxVtx][300];   //[Vtx_]
   Float_t         Vtx_NEabstpvert[kMaxVtx][300];   //[Vtx_]
   Int_t           Vtx_NEipvert[kMaxVtx][300];   //[Vtx_]
   Int_t           Vtx_NEiverti[kMaxVtx][300];   //[Vtx_]
   Int_t           Vtx_NEivertf[kMaxVtx][300];   //[Vtx_]

   // List of branches
   TBranch        *b_NVtx;   //!

   TBranch        *b_Neutmode;   //!
   TBranch        *b_Npvc;   //!
   TBranch        *b_Ipvc;   //!
   TBranch        *b_Pvc;   //!
   TBranch        *b_Abspvc;   //!


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


   // Class variables

   Int_t MCID;  // 0: SK, 1: ND280-numc, 2: ND280-oaAnalysis, 3: Pion-Scattering/Photo-production
   enum eMCID {sk, nd280numc, nd280anal, piscat};

   Int_t    nParSets;
   Double_t **pars;

   TString theFile;
   TString treeKey;

   TFile *mcRwOutFile;

   virtual void     Init();
   
   // ROOT file management
   void openFile();     // Open input file as appendable
   void copyFile();     // Copy input file to new file
   void makeParTree();  // Create paramater tree and new branches 
   
   // Functions to copy ntuple variables to Fortran common blocks
   void fill_nework(Int_t ivtx=0);
   void fill_necard(Int_t iparset=0);
   void fill_fsihist(Int_t ivtx=0);
   void fill_neutparams(Int_t iparset=0);
   void fill_neutmodel(Int_t iparset=0);

   void fill_allparams(Int_t iparset=0);
   void fill_allevent(Int_t ivtx=0);

   // Print 
   int verbose;

   void print_nework();
   void print_fsihist();
   void print_allparams();
   void print_allevent();

 public:
   t2kneutreweight(TString aTheFile, Int_t aMCID, Int_t nParSets, Double_t **apars=NULL, int averbose=0);
   virtual ~t2kneutreweight();
   
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);

   // Functions to call NEUT calculation fortran programs
   void appendWeightsToTree();


};
