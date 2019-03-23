#define piscatana_cxx

//#define VERBOSE

#include <sstream>
#include <sys/stat.h> 

#include <TMath.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>

#include "piscatana.h"

#include "PeriodicTable.h"

#include "necardC.h"
#include "efpionC.h"

// For using old histogram files produced without beam/nucleus type information
// (This information must be hard coded below in function InitFromHistos()
//#define OLD_PISCAT_FILES

// For generating histograms of total cross section only to save time and space
// when doing parameter fits
#define TOT_XSEC_ONLY

using namespace std;

double legTopRight[4] = {0.5,0.65,0.95,0.9};
double legTopLeft[4] = {0.14,0.65,0.39,0.9};
//double legTopLeft[4] = {0.57,0.52,0.93,0.91}; // For large legend in old 1-sig var absorption plot

// Analyse Tree
piscatana::piscatana(string stringInput1, string stringInput2)
{

  theFile = stringInput1; // Input file(s)
  outFile = stringInput2; // Output file
  InitFromTree();
 
}

// Load already analysed histograms
piscatana::piscatana(string stringInput1, string stringInput2, int ExpNumEvents)
{

  // Parse infiles string
   
  // Assumes input string is already space delimited
  std::stringstream infile_ss(stringInput1); 
  std::stringstream legname_ss(stringInput2); 
  
  // Parse and store each file path
  string tmpString;
  vector<string> infile_parse;
  nConfs=0;
  while (infile_ss >> tmpString) {
    infile_parse.push_back(tmpString);

    tmpString = "";
    legname_ss >> tmpString;
    legName_parse.push_back(tmpString);
 
    nConfs++;
  }

  InitFromHistos(infile_parse,ExpNumEvents);

  LoadData();
}

piscatana::~piscatana()
{
   if (!fChain) return;
   
   outROOT->Write();

   cout << endl << "Wrote histograms to: " << outFile.c_str() << endl << endl;

  for (int ipion=0; ipion<nPions; ipion++) {
    for (int iNuc=0; iNuc<kZmax; iNuc++) {
      delete histosLoaded[ipion][iNuc];
    }
  }


   delete fChain->GetCurrentFile();
}

Int_t piscatana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


Long64_t piscatana::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
   }
   return centry;
}



void piscatana::InitFromTree() {

  cout << endl << "Loading piscat tree from file(s): " << theFile.c_str() << endl;
  
  fChain = new TChain("h1","");
  fChain->Add(Form("%s/h1",theFile.c_str()));

  // Set branch addresses and branch pointers
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("Npvc", &Npvc, &b_Npvc);
  fChain->SetBranchAddress("Ipvc", Ipvc, &b_Ipvc);
  fChain->SetBranchAddress("Iorgvc", Iorgvc, &b_Iorgvc);
  fChain->SetBranchAddress("Ichvc", Ichvc, &b_Ichvc);
  fChain->SetBranchAddress("Iflvc", Iflvc, &b_Iflvc);
  fChain->SetBranchAddress("Abspvc", Abspvc, &b_Abspvc);
  fChain->SetBranchAddress("Pvc", Pvc, &b_Pvc);
  fChain->SetBranchAddress("Posvc", Posvc, &b_Posvc);
  fChain->SetBranchAddress("Nvert", &Nvert, &b_Nvert);
  fChain->SetBranchAddress("Posvert", Posvert, &b_Posvert);
  fChain->SetBranchAddress("Iflgvert", Iflgvert, &b_Iflgvert);
  fChain->SetBranchAddress("Nvcvert", &Nvcvert, &b_Nvcvert);
  fChain->SetBranchAddress("Dirvert", Dirvert, &b_Dirvert);
  fChain->SetBranchAddress("Abspvert", Abspvert, &b_Abspvert);
  fChain->SetBranchAddress("Abstpvert", Abstpvert, &b_Abstpvert);
  fChain->SetBranchAddress("Ipvert", Ipvert, &b_Ipvert);
  fChain->SetBranchAddress("Iverti", Iverti, &b_Iverti);
  fChain->SetBranchAddress("Ivertf", Ivertf, &b_Ivertf);
  fChain->SetBranchAddress("Fsiprob", &Fsiprob, &b_Fsiprob);
  fChain->SetBranchAddress("Numbndn", &Numbndn, &b_Numbndn);
  fChain->SetBranchAddress("Numbndp", &Numbndp, &b_Numbndp);
  fChain->SetBranchAddress("Numfrep", &Numfrep, &b_Numfrep);
  fChain->SetBranchAddress("Numatom", &Numatom, &b_Numatom);

  // Set initial pion type and target nucleus
  fChain->GetEntry(0);
  SetInitPion();
  SetNucleus();

  outROOT = new TFile(outFile.c_str(),"RECREATE"); 
  
  Init_histo();
}


void piscatana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
 
void piscatana::SetInitPion() {

  if (Ipvc[4] == 111) {init_pion=0;}
  else if (Ipvc[4] == 211) {init_pion=1;}
  else if (Ipvc[4] == -211) {init_pion=2;}
  else {
    cout << "piscatana::SetInitPion() Error: Unknown pion PID: " << Ipvc[4] << endl;
    exit (-1);
  }

  cout << endl << "Incident pion type: " << pionName[init_pion] << endl;
}

void piscatana::SetNucleus()
{
  neuttarget_.numbndn = Numbndn;
  neuttarget_.numbndp = Numbndp;
  neuttarget_.numfrep = Numfrep;
  neuttarget_.numatom = Numatom;

  // Set nucleus and get parameters for normalization
  nesetfgparams_();  // Check of existence in NEUT contained

  cout << endl << "Target nucleus: " << element_symbol[Numbndp-1] << endl;

  double rMax = 2.5*eftarget_.c;
  targetArea = TMath::Pi()*rMax*rMax*10;  // mb
}

// Must be same naming convention as in Init_histo()
void piscatana::InitFromHistos(vector<string> inputfile_arr, int ExpNumEvents) {
  
  if ( ExpNumEvents <= 0 )
    cout << endl << "piscatana::InitFromHistos() Warning: Input ExpNumEvents <= 0, differential histograms may not be normalized properly" << endl;

  cout << endl << "Loading piscat histos from file(s): " << endl;

  // Initialize flag for histos loaded 
  for (int ipion=0; ipion<nPions; ipion++) {
    for (int iNuc=0; iNuc<kZmax; iNuc++) {
      histosLoaded[ipion][iNuc] = new bool[nConfs];
    }
  }

  // Set Normalization for old files
#ifdef OLD_PISCAT_FILES
  Numbndn=6;
  Numbndp=6;
  Numfrep=0;
  Numatom=12;
  SetNucleus();
#endif

  TFile *infile[nConfs];

  for (int iconf=0; iconf<nConfs; iconf++) {

    infile[iconf] = new TFile(inputfile_arr[iconf].c_str());

    cout << inputfile_arr[iconf].c_str() << endl;

    // Broken or misnamed file
    if (!infile[iconf]->IsOpen()) exit (-1);

    for (int ipion=0; ipion<nPions; ipion++) {
      for (int iNuc=0; iNuc<kZmax; iNuc++) {

#ifdef OLD_PISCAT_FILES
	if (ipion!=1 || iNuc!=5) continue;
#endif
	
	char *charName;
 
	// Initial Spectrum
#ifdef OLD_PISCAT_FILES
	charName = Form("total_init");
#else
	charName = Form("total_init_%s_%s",pionName[ipion],element_symbol[iNuc]);
#endif

	// Check if pion type and nucleus type combination exists
	TH1F *h_temp = (TH1F*)infile[iconf]->Get(charName);
	if (!h_temp) {
	  histosLoaded[ipion][iNuc][iconf] = false;
	  continue;
	}

	// Flag for given pion-nucleon-conf combination loaded
	histosLoaded[ipion][iNuc][iconf] = true;

	piInitSpec[ipion][iNuc].push_back(new TH1F(*(TH1F*)infile[iconf]->Get(charName)));

	// Deduce number of files by using number of triggers
	float nFiles = 1;
	if ( ExpNumEvents > 0 )
	  nFiles = piInitSpec[ipion][iNuc][iconf]->GetEntries() / ExpNumEvents;

	// Interaction separated histograms
	for (int iint=0; iint<nInts; iint++) {
#ifdef OLD_PISCAT_FILES
	  charName = Form("%s",intName[iint]);
#else
	  charName = Form("%s_%s_%s",intName[iint],pionName[ipion],element_symbol[iNuc]);
#endif
	  piInitIntSpec[ipion][iNuc][iint].push_back(new TH1F(*(TH1F*)infile[iconf]->Get(charName)));
#ifndef TOT_XSEC_ONLY // Already normalized assuming TOT_XSEC_ONLY uses only one file
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Divide(piInitSpec[ipion][iNuc][iconf]);
#endif
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Smooth();
#ifdef OLD_PISCAT_FILES
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Scale(targetArea);
#endif
	}

#ifndef TOT_XSEC_ONLY
	// Pi-production 
	for (int ipiprod=0; ipiprod<nPiprods; ipiprod++) {
      
#ifdef OLD_PISCAT_FILES
	  charName = Form("%s",piprodName[ipiprod]);
#else
	  charName = Form("%s_%s_%s",piprodName[ipiprod],pionName[ipion],element_symbol[iNuc]);
#endif
	  piInitPiprodSpec[ipion][iNuc][ipiprod].push_back(new TH1F(*(TH1F*)infile[iconf]->Get(charName)));
	  piInitPiprodSpec[ipion][iNuc][ipiprod][iconf]->Divide(piInitSpec[ipion][iNuc][iconf]);
#ifdef OLD_PISCAT_FILES
	  piInitPiprodSpec[ipion][iNuc][ipiprod][iconf]->Scale(targetArea);;
#endif
	}
        
	for (int pion=0; pion<3; pion++) {
	  // Phase space
#ifdef OLD_PISCAT_FILES
	  charName = Form("out_ps_%s",pionName[pion]);
#else
	  charName = Form("out_ps_%s_%s_%s",pionName[ipion],element_symbol[iNuc],pionName[pion]);
#endif 
	  outPhaseSpc[ipion][iNuc][pion].push_back(new TH2F(*(TH2F*)infile[iconf]->Get(charName)));
	  outPhaseSpc[ipion][iNuc][pion][iconf]->Scale(1/nFiles);
#ifdef OLD_PISCAT_FILES
	  outPhaseSpc[ipion][iNuc][pion][iconf]->Scale(targetArea);
#endif

	  // Differential Cross Section
#ifdef OLD_PISCAT_FILES
	  charName = Form("diff_%s",pionName[pion]);
#else
	  charName = Form("diff_%s_%s_%s",pionName[ipion],element_symbol[iNuc],pionName[pion]);
#endif
	  piDiffXsec[ipion][iNuc][pion].push_back(new TH1F(*(TH1F*)infile[iconf]->Get(charName)));
	  piDiffXsec[ipion][iNuc][pion][iconf]->Scale(1/nFiles);
#ifdef OLD_PISCAT_FILES
	  piDiffXsec[ipion][iNuc][pion][iconf]->Scale(targetArea);
#endif
	}
#endif // TOT_XSEC_ONLY	
      }
    }
  }
  
}


void piscatana::Init_histo()
{
  
  cout << endl << "Initializing piscat MC histograms" << endl;

  // Define Histogram range (must be modified depending on range sample was generated with)
  float dE = 10; // MeV/c
  float s_xMin = 0;
  float s_xMax = 1600; 
  int s_nBins = (s_xMax-s_xMin)/dE;


  // Outgoing Pion Phase Space
  int a_nBins = 90;
  float a_xMin = 0;
  float a_xMax = 180; // degrees
  dtheta = (a_xMax-a_xMin)/a_nBins * TMath::Pi()/180.;

  char *charName,*charTitle,*charLabelX,*charLabelY;
 
  // Initial Spectrum
  charName = Form("total_init_%s_%s",pionName[init_pion],element_symbol[Numbndp-1]);
  charTitle = Form("Total");
  charLabelX = Form("%s Initial Momentum (MeV/c)",pionTitle[init_pion]);
    
  piInitSpec[init_pion][Numbndp-1].push_back(new TH1F(charName,charTitle, s_nBins, s_xMin, s_xMax));
  piInitSpec[init_pion][Numbndp-1][0]->GetXaxis()->SetTitle(charLabelX);
  piInitSpec[init_pion][Numbndp-1][0]->Sumw2();


  // Interaction separated histograms
  for (int iint=0; iint<nInts; iint++) {
    charName = Form("%s_%s_%s",intName[iint],pionName[init_pion],element_symbol[Numbndp-1]);
    charTitle = Form("%s",intTitle[iint]);
    charLabelX = Form("%s Initial Momentum (MeV/c)",pionTitle[init_pion]);
    charLabelY = "#sigma (mb)";
    piInitIntSpec[init_pion][Numbndp-1][iint].push_back(new TH1F(charName,Form("%s;%s;%s",charTitle, charLabelX, charLabelY), s_nBins, s_xMin, s_xMax));
    piInitIntSpec[init_pion][Numbndp-1][iint][0]->Sumw2();
  }

#ifndef TOT_XSEC_ONLY

  // Pi-production 
  for (int ipiprod=0; ipiprod<nPiprods; ipiprod++) {
      
    charName = Form("%s_%s_%s",piprodName[ipiprod],pionName[init_pion],element_symbol[Numbndp-1]);
    charTitle = Form("%s",piprodTitle[ipiprod]);
    charLabelX = Form("%s Initial Momentum (MeV/c)",pionTitle[init_pion]);
    
    piInitPiprodSpec[init_pion][Numbndp-1][ipiprod].push_back(new TH1F(charName,charTitle, s_nBins, s_xMin, s_xMax));
    piInitPiprodSpec[init_pion][Numbndp-1][ipiprod][0]->GetXaxis()->SetTitle(charLabelX);
    piInitPiprodSpec[init_pion][Numbndp-1][ipiprod][0]->Sumw2();
  }
        
  for (int pion=0; pion<3; pion++) {
    // Phase space
    charName = Form("out_ps_%s_%s_%s",pionName[init_pion],element_symbol[Numbndp-1],pionName[pion]);
    charTitle = Form("Outgoing %s Phase Space",pionTitle[pion]);
    charLabelX = Form("%s Momentum (MeV/c)",pionTitle[pion]);
    charLabelY = Form("#theta_{lab,%s}",pionTitle[pion]);
    //char *charLabelZ = Form("#frac{d#sigma}{d#OmegadE} (mb/sr/%dMeV)",(int)((s_xMax-s_xMin)/s_nBins));
    char *charLabelZ = Form("#frac{d#sigma}{d#OmegadE} (mb)");
    
    outPhaseSpc[init_pion][Numbndp-1][pion].push_back(new TH2F(charName,charTitle, s_nBins, s_xMin, s_xMax, a_nBins,a_xMin,a_xMax));
    outPhaseSpc[init_pion][Numbndp-1][pion][0]->GetXaxis()->SetTitle(charLabelX);
    outPhaseSpc[init_pion][Numbndp-1][pion][0]->GetYaxis()->SetTitle(charLabelY);
    outPhaseSpc[init_pion][Numbndp-1][pion][0]->GetZaxis()->SetTitle(charLabelZ);
    outPhaseSpc[init_pion][Numbndp-1][pion][0]->Sumw2();


    // Differential Cross Section
    charName = Form("diff_%s_%s_%s",pionName[init_pion],element_symbol[Numbndp-1],pionName[pion]);
    charTitle = Form("%s #frac{d#sigma}{d#Omega}",pionTitle[pion]);
    charLabelX = Form("#theta_{lab,%s}",pionTitle[pion]);
    charLabelY = Form("#frac{d#sigma}{d#Omega} (mb/sr)");
    
    piDiffXsec[init_pion][Numbndp-1][pion].push_back(new TH1F(charName,charTitle, a_nBins, a_xMin, a_xMax));
    piDiffXsec[init_pion][Numbndp-1][pion][0]->GetXaxis()->SetTitle(charLabelX);
    piDiffXsec[init_pion][Numbndp-1][pion][0]->GetYaxis()->SetTitle(charLabelY);
    piDiffXsec[init_pion][Numbndp-1][pion][0]->Sumw2();
  }
  #endif // TOT_XSEC_ONLY

  
}



// Initial Interaction Spectrum
void piscatana::int_xsec()
{


  Float_t evt_weight;
  evt_weight = 1;
    
  // Survival / Normal
  if (Iflvc[4]==0) {
    piInitIntSpec[init_pion][Numbndp-1][nInts-1][0]->Fill(Abspvc[4],evt_weight);
  }
  else {

    int pion=-1;
    int piID=-1;
    int nPi=0;
    int nPiType[3];
    nPiType[0]=nPiType[1]=nPiType[2]=0;
    
    float cost, px[3], pabs;
    
    // Determine final state
    for (int ipvc=4; ipvc<Npvc; ipvc++) {
      
      if (Ichvc[ipvc]==0) continue;
      
      // Check if it's a pion
      if (Ipvc[ipvc] == 111) {pion=0;}
      else if (Ipvc[ipvc] == 211) {pion=1;}
      else if (Ipvc[ipvc] == -211) {pion=2;}
      else continue;
      
      px[0] = Pvc[ipvc][0];
      px[1] = Pvc[ipvc][1];
      px[2] = Pvc[ipvc][2];
      pabs = Abspvc[ipvc];
      cost = px[0]/pabs;
      
      piID=ipvc;

      nPiType[pion]++;      
      nPi++;
    }
        
    // Absorbed
    if (nPi==0) {
      piInitIntSpec[init_pion][Numbndp-1][1][0]->Fill(Abspvc[4],evt_weight);
    }
    
    // Quasielastic Scattering
    else if (nPi==1 && pion>=0) { // && piID>=0
      
      // Single Charge Exchange (assuming piP/M beam)
      if (nPiType[0]) {
	piInitIntSpec[init_pion][Numbndp-1][3][0]->Fill(Abspvc[4],evt_weight);
	
      }
      // Double Charge Exchange
      else if ( (init_pion==1 && pion==2) || (init_pion==2 && pion==1) )
	piInitIntSpec[init_pion][Numbndp-1][4][0]->Fill(Abspvc[4],evt_weight);
      
      // Elastic Scattering 
      //else if ( Iflvc[Iorgvc[piID]-1] == 9 )
      //	piInitIntSpec[5][0]->Fill(Abspvc[4],evt_weight);
       
      // Quasielastic
      else
	piInitIntSpec[init_pion][Numbndp-1][2][0]->Fill(Abspvc[4],evt_weight);
      
#ifndef TOT_XSEC_ONLY
      // Outgoing Diff. Xsec & Phase space
      float theta = acos(cost)*180/TMath::Pi();
      float sint = sqrt(1-cost*cost);
      if (sint==0) sint=0.000001;
      
      float weight_ang =  1/sint/(2*TMath::Pi()*dtheta);
      
      piDiffXsec[init_pion][Numbndp-1][pion][0]->Fill(theta,weight_ang*evt_weight);
      outPhaseSpc[init_pion][Numbndp-1][pion][0]->Fill(pabs,theta,evt_weight); // /dE);
#endif
    }
    
    // Hadronic Production
    else { // if (nPi>1)
      piInitIntSpec[init_pion][Numbndp-1][6][0]->Fill(Abspvc[4],evt_weight);

#ifndef TOT_XSEC_ONLY
      // 2pi production
      if (nPi==2) {
      
	int ipiprod=0;
      
	piInitPiprodSpec[init_pion][Numbndp-1][ipiprod][0]->Fill(Abspvc[4],evt_weight);
      
	if (nPiType[0]==2) 
	  ipiprod=1; // pi0pi0
	else if (nPiType[1]==1 && nPiType[2]==1) 
	  ipiprod=2; // piPpiM
	else if (nPiType[1]==1 && nPiType[0]==1) 
	  ipiprod=3; // piPpi0
	else if (nPiType[1]==2) 
	  ipiprod=4; // piPpiP
	else if (nPiType[2]==1 && nPiType[0]==1) 
	  ipiprod=5; // piMpi0
	else if (nPiType[2]==2) 
	  ipiprod=6; // piMpiM
      
	piInitPiprodSpec[init_pion][Numbndp-1][ipiprod][0]->Fill(Abspvc[4],evt_weight);    
      
      }
#endif

    }

  }
}

void piscatana::Loop()
{
  if (fChain == 0) return;
  
  // Loop through events
  nentries = fChain->GetEntriesFast();
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //if (jentry%1000==0) cout << "Event: " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    
    // Total Initial Spectrum
    piInitSpec[init_pion][Numbndp-1][0]->Fill(Abspvc[4]);

    // Analysis
    int_xsec();
  }  // End event loop


  // Total and reactive xsec

  piInitIntSpec[init_pion][Numbndp-1][0][0]->Reset();
  piInitIntSpec[init_pion][Numbndp-1][nInts-2][0]->Reset();
  for (int iint=1; iint<nInts-2; iint++) {

    // Total
    piInitIntSpec[init_pion][Numbndp-1][nInts-2][0]->Add(piInitIntSpec[init_pion][Numbndp-1][iint][0]);
    
    // Reactive
    if (iint==5) continue;
    piInitIntSpec[init_pion][Numbndp-1][0][0]->Add(piInitIntSpec[init_pion][Numbndp-1][iint][0]);
  }

    
  // Normalize
  for (int iint=0; iint<nInts; iint++) {
    //piInitIntSpec[iint]->Divide(piInitSpec);
    piInitIntSpec[init_pion][Numbndp-1][iint][0]->Scale(targetArea);
  }

#ifndef TOT_XSEC_ONLY
  for (int ipiprod=0; ipiprod<nPiprods; ipiprod++) {
    //piInitPiprodSpec[ipiprod]->Divide(piInitSpec);
    piInitPiprodSpec[init_pion][Numbndp-1][ipiprod][0]->Scale(targetArea);
  }
  
  for (int ipion=0; ipion<3; ipion++) {
    piDiffXsec[init_pion][Numbndp-1][ipion][0]->Scale(targetArea/(float)nentries);
    outPhaseSpc[init_pion][Numbndp-1][ipion][0]->Scale(targetArea/(float)nentries);
  }
#endif
}



// PLOTTING //

void piscatana::LoadData()
{
  char dirname[500];
  
  if (!getenv("NIWG_DATA")) {
    cout << "piscatana::LoadData() Error: Set \"NIWG_DATA\" environment variable" << endl;
    exit (-1);
  } else {
    sprintf(dirname,"%s",getenv("NIWG_DATA"));
    cout << endl << "Loading Piscat data files from directory: " << dirname << endl;
  }

  // Load data
  //TFile *outdata_file = new TFile(Form("%s/piscat_data.root",extdata_folder),"RECREATE");

  for (int pion=0; pion<nPions; pion++) {
    for (int iNuc=0; iNuc<kZmax; iNuc++) {
      for (int iint=0; iint<nInts; iint++) {

	string element_symbol_lowercase = element_symbol[iNuc];
	element_symbol_lowercase[0] = tolower(element_symbol_lowercase[0]);

	string filename = Form("%s/piscat_xsec/%s_%s_%s.csv",dirname,element_symbol_lowercase.c_str(),intName[iint],pionName[pion]);
	
	int_graph[pion][iNuc][iint] = NULL;
	if (!FileExists(filename)) {
	  continue;
	}

	//cout << filename << endl;
	
	int_graph[pion][iNuc][iint] = new TGraphErrors(filename.c_str(),"%lg %lg %lg");
	
	int_graph[pion][iNuc][iint]->SetName(Form("%s_%s_%s",element_symbol[iNuc],intName[iint],pionName[pion]));
	int_graph[pion][iNuc][iint]->SetTitle(Form("%s %s %s",pionTitle[pion],element_name[iNuc],intTitle[iint]));
	int_graph[pion][iNuc][iint]->GetXaxis()->SetTitle(Form("%s Initial Momentum (MeV/c)",pionTitle[pion]));
	int_graph[pion][iNuc][iint]->GetYaxis()->SetTitle("#sigma (mb)");
	
	int_graph[pion][iNuc][iint]->SetMarkerColor(kGray+2);
	int_graph[pion][iNuc][iint]->SetLineColor(kGray+2);
	
	int_graph[pion][iNuc][iint]->SetMarkerSize(2);
	int_graph[pion][iNuc][iint]->SetMarkerStyle(20);
	int_graph[pion][iNuc][iint]->SetLineWidth(3);

	//if (iint==nInts-2) {
	//  int_graph[pion][iNuc][iint]->SetMarkerColor(tot_color);
	//  int_graph[pion][iNuc][iint]->SetLineColor(tot_color);
	//}

	//int_graph[pion][iNuc][iint]->Write();

      }
    }
  }

}


bool piscatana::FileExists(string strFilename)
{
  struct stat stFileInfo; 
  bool blnReturn; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if(intStat == 0) { 
    blnReturn = true; 
  } else { 
    blnReturn = false; 
  } 
   
  return(blnReturn); 
}



// Plots all interaction channels on a separate canvas for each pion and nucleus type
void piscatana::plot_sep_nuc()
{
  
  for (int ipion=0; ipion<nPions; ipion++) {
    for (int iNuc=0; iNuc<kZmax; iNuc++) {

      if (!histosLoaded[ipion][iNuc][0]) continue;
      
#ifdef OLD_PISCAT_FILES
      string c_name = Form("piscat_hinrgvar_%s",pionName[ipion]);
#else
      string c_name = Form("piscat_oldvstun_%s_%s",pionName[ipion],element_symbol[iNuc]);
#endif
      TCanvas *c = new TCanvas(c_name.c_str(),c_name.c_str(),0,0,1200,900);
      bool isDrawn = false;

      double maxCanvasHeight = 0;

      for (int iint=0; iint<=6; iint++) {
      	for (int iconf=0; iconf<nConfs; iconf++) {
	  if (piInitIntSpec[ipion][iNuc][iint][iconf]->GetMaximum() > maxCanvasHeight) {
	    maxCanvasHeight = piInitIntSpec[ipion][iNuc][iint][iconf]->GetMaximum()*1.05;
	  }
	}
      }

      for (int iint=0; iint<=6; iint++) {

	if (iint==dcx || iint==elas) continue;

	for (int iconf=0; iconf<nConfs; iconf++) {
	  
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Rebin(4);
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Scale(1/4.);
      	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineColor(iint+1);
	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineWidth(2);

	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetMinimum(0);
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetRangeUser(0,maxCanvasHeight);

	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetTitle("");
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetXaxis()->SetTitle(Form("%s Initial Momentum (MeV/c)",pionTitle[ipion]));
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetTitle("#sigma (mb)");

#ifndef OLD_PISCAT_FILES
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetTitleOffset(1.5);
#endif
	  if (iconf==0)
	    piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineStyle(1);
	  else if (iconf==1)
	    piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineStyle(2);

	  if (!isDrawn) {
	    piInitIntSpec[ipion][iNuc][iint][iconf]->Draw("C hist");
	    isDrawn=true;
	  }
	  else
	    piInitIntSpec[ipion][iNuc][iint][iconf]->Draw("C hist same");

      	}
      
	// Data
	if (int_graph[ipion][iNuc][iint]) {
	  
	  Int_t nPoints = int_graph[ipion][iNuc][iint]->GetN();

	  if (nPoints) {

	    for (int i=0; i<nPoints; i++) {
	      Double_t x,y,y_err;
	      int_graph[ipion][iNuc][iint]->GetPoint(i,x,y);
	      y_err = int_graph[ipion][iNuc][iint]->GetErrorY(i);

	      if (y+y_err > maxCanvasHeight) 
		maxCanvasHeight = (y+y_err)*1.02;
	    }

	    int_graph[ipion][iNuc][iint]->SetLineWidth(1);
	    int_graph[ipion][iNuc][iint]->SetLineColor(iint+1);
	    int_graph[ipion][iNuc][iint]->SetMarkerColor(iint+1);
	    int_graph[ipion][iNuc][iint]->SetMarkerStyle(4);

	    int_graph[ipion][iNuc][iint]->Draw("same p");
	  }
	}
      }

      // Reset maximum height change due to data
      for (int iint=0; iint<nInts; iint++) {
	for (int iconf=0; iconf<nConfs; iconf++) {
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetRangeUser(0,maxCanvasHeight);
	}
      }

      DrawLegendInts(piInitIntSpec[ipion][iNuc],legTopRight[0],legTopRight[1],legTopRight[2],legTopRight[3],iNuc);

#ifndef OLD_PISCAT_FILES
      TLatex rxnLabel(500,maxCanvasHeight*.8,Form("%s {}^{%d}_{%d}%s",pionTitle[ipion],a_stable[iNuc],iNuc+1,element_symbol[iNuc]));
      rxnLabel.SetTextSize(0.05);
      rxnLabel.Draw("SAME");
#endif
	
      print_canvas(c);
      
    }
  }


}


// Plots all configurations on a separate canvas for each interaction channel, pion and nucleus type
void piscatana::plot_sep_pion_nuc_iint() 
{

  for (int ipion=0; ipion<nPions; ipion++) {
    for (int iNuc=0; iNuc<kZmax; iNuc++) {

      if (!histosLoaded[ipion][iNuc][0]) continue;
      
      for (int iint=0; iint<=scx; iint++) {
      
      	string c_name = Form("piscat_1sigvars_%s",intName[iint]);
      	TCanvas *c = new TCanvas(c_name.c_str(),c_name.c_str(),0,0,1050,900);
	THStack *h_stack = new THStack();
      
      	for (int iconf=0; iconf<nConfs; iconf++) {
      	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineColor(iconf+1);
      	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineWidth(3);

	  // Rescale for aesthetics
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Rebin(4);
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Scale(1/4.);

      	  h_stack->Add(piInitIntSpec[ipion][iNuc][iint][iconf],"c hist");
      	}
      
      	h_stack->Draw("nostack");
	//h_stack->GetXaxis()->SetRangeUser(0,599);
	h_stack->GetXaxis()->SetTitle(Form("%s Initial Momentum (MeV/c)",pionTitle[ipion]));
	h_stack->GetYaxis()->SetTitle(Form("#sigma_{%s} (mb)",intTitle[iint]));
	h_stack->GetYaxis()->SetTitleOffset(1.2);

	// Data
	if (int_graph[ipion][iNuc][iint]) {
	  if (int_graph[ipion][iNuc][iint]->GetN()) {
	    int_graph[ipion][iNuc][iint]->Draw("same p");
	  }
	}

	if (DrawLegendConfs(piInitIntSpec[ipion][iNuc][iint],legTopLeft[0],legTopLeft[1],legTopLeft[2],legTopLeft[3]))
	  print_canvas(c,"_leg");
	else
	  print_canvas(c);
      }
    }
  }

}

// Plots all nuclei on a separate canvas for each interaction channel and pion type
void piscatana::plot_sep_pion_iint() 
{
  // See neutgeom/PeriodicTable.h
  int nucLineColor[kZmax] = {1};
  nucLineColor[elem_C] = 1;
  nucLineColor[elem_O] = 2;
  nucLineColor[elem_Al] = 4;
  nucLineColor[elem_Fe] = kOrange+7;
  nucLineColor[elem_Nb] = 6;
  nucLineColor[elem_Pb] = kGreen+2;

  // Data point marker styles
  double nucMarkerStyle[kZmax] = {1};
  nucMarkerStyle[elem_C] = 20;
  nucMarkerStyle[elem_O] = 21;
  nucMarkerStyle[elem_Al] = 22;
  nucMarkerStyle[elem_Fe] = 33;
  nucMarkerStyle[elem_Nb] = 23;
  nucMarkerStyle[elem_Pb] = 34;

  // Scale factors for histogram and data to appear nicely on single plot
  double nucScaleFactor[kZmax] = {1};
  nucScaleFactor[elem_C] = 1;
  nucScaleFactor[elem_O] = 1.7;
  nucScaleFactor[elem_Al] = 1.8;
  nucScaleFactor[elem_Fe] = 1.5;
  nucScaleFactor[elem_Nb] = 1.5;
  nucScaleFactor[elem_Pb] = 1.2;

  // x-axis ranges
  double xAxisRange[nPions][kZmax][nInts] = {0};
  for (int iNuc=kZmax-1; iNuc>=0; iNuc--) {

    xAxisRange[pip][iNuc][absb]  = 700;
    xAxisRange[pip][iNuc][dcx]   = 500;
    xAxisRange[pip][iNuc][scatt] = 500;
    xAxisRange[pip][iNuc][scx]   = 700;
    xAxisRange[pim][iNuc][absb]  = 500;
    xAxisRange[pim][iNuc][dcx]   = 500;
    xAxisRange[pim][iNuc][scatt] = 1000;
    xAxisRange[pim][iNuc][scx]   = 400;

  }



  for (int ipion=0; ipion<nPions; ipion++) {
    for (int iint=0; iint<=dcx; iint++) {
    
      string c_name = Form("piscat_%s_%s_ascale",pionName[ipion],intName[iint]);
      TCanvas *c = new TCanvas(c_name.c_str(),c_name.c_str(),0,0,1050,900);
      
      // Set legend position and attributes
      TLegend *leg;
      if ((ipion==pim && iint==scx) || (ipion==pip && iint==scatt))
	leg = new TLegend(0.16,0.5,0.39,0.9);
      else if (ipion==pip && iint==scx) {
	leg = new TLegend(0.16,0.7,0.7,0.9);
	leg->SetNColumns(2);
      }
      else if (iint==dcx)
	leg = new TLegend(0.16,0.65,0.4,0.9);
      else if (iint==rxn) {
	leg = new TLegend(0.35,0.72,0.9,0.9);
	leg->SetNColumns(2);
      }
      else if (ipion==pim && iint==scatt)
	leg = new TLegend(0.6,0.65,0.95,0.9);
      else
	leg = new TLegend(0.6,0.6,0.95,0.9);

      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetLineStyle(0);
      leg->SetLineWidth(0);
      leg->SetLineColor(0);

      double maxCanvasHeight = 0;
      bool isDrawn = false;

      for (int iNuc=kZmax-1; iNuc>=0; iNuc--) {
	
	// Data exists?
	if (!int_graph[ipion][iNuc][iint]) continue;
	if (!int_graph[ipion][iNuc][iint]->GetN()) continue;

	// MC exists?
	if (!histosLoaded[ipion][iNuc][0]) continue;

      	for (int iconf=0; iconf<nConfs; iconf++) {

	  // Rescale for aesthetics
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Rebin(4);
	  piInitIntSpec[ipion][iNuc][iint][iconf]->Scale(1/4.);

	  if (iint!=dcx)
	    piInitIntSpec[ipion][iNuc][iint][iconf]->Scale(nucScaleFactor[iNuc]);

	  // Find maximum height
	  if (piInitIntSpec[ipion][iNuc][iint][iconf]->GetMaximum() > maxCanvasHeight)
	    maxCanvasHeight = piInitIntSpec[ipion][iNuc][iint][iconf]->GetMaximum()*1.03;
	}

		
	for (int iconf=0; iconf<nConfs; iconf++) {
	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineColor(nucLineColor[iNuc]);
	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineStyle(iconf+1);
      	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetLineWidth(3);

	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetTitle(""); 
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetXaxis()->SetTitle(Form("%s Initial Momentum (MeV/c)",pionTitle[ipion])); 
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetTitle(Form("#sigma_{%s} (mb)",intTitle[iint]));	       
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetTitleOffset(1.45);

	  piInitIntSpec[ipion][iNuc][iint][iconf]->SetMinimum(0);
	  piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetRangeUser(0,maxCanvasHeight);

	  if (xAxisRange[ipion][iNuc][iint])
	    piInitIntSpec[ipion][iNuc][iint][iconf]->GetXaxis()->SetRangeUser(0,xAxisRange[ipion][iNuc][iint]);
	  
	  if (!isDrawn)
	    piInitIntSpec[ipion][iNuc][iint][iconf]->Draw("C HIST");
	  else
	    piInitIntSpec[ipion][iNuc][iint][iconf]->Draw("C HIST SAME");
	     
      	}
      
	// Data
	if (int_graph[ipion][iNuc][iint]) {

	  int nPoints = int_graph[ipion][iNuc][iint]->GetN();

	  if (nPoints) {

	    for (int i=0; i<nPoints; i++) {
	      Double_t x,y;
	      int_graph[ipion][iNuc][iint]->GetPoint(i,x,y);

	      if (iint!=dcx) {
		int_graph[ipion][iNuc][iint]->SetPoint(i,x,y*nucScaleFactor[iNuc]);
		int_graph[ipion][iNuc][iint]->SetPointError(i,0,0);
	      }

	      if (y*nucScaleFactor[iNuc] > maxCanvasHeight) 
		maxCanvasHeight = (y*nucScaleFactor[iNuc])*1.02;
	    }
	      
	    // Final canvas height is done here (based on max height of data and MC)
	    for (int iconf=0; iconf<nConfs; iconf++) {
	      piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetRangeUser(0,maxCanvasHeight);
	      
	      if (iint==dcx && ipion==pim)
		piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetRangeUser(0,50);

	      if (iint==scx && ipion==pim)
		piInitIntSpec[ipion][iNuc][iint][iconf]->GetYaxis()->SetRangeUser(0,350);
	    }

	    // Set data graph attributes
	    int_graph[ipion][iNuc][iint]->SetLineColor(nucLineColor[iNuc]);
	    int_graph[ipion][iNuc][iint]->SetMarkerStyle(nucMarkerStyle[iNuc]);
	    int_graph[ipion][iNuc][iint]->SetMarkerColor(nucLineColor[iNuc]);
	    int_graph[ipion][iNuc][iint]->Draw("same p");
	    
	    // Add legend labels
	    if (iint!=dcx && nucScaleFactor[iNuc]!=1)
	      leg->AddEntry(int_graph[ipion][iNuc][iint],Form("^{%d}%s (#times%g)",a_stable[iNuc],element_symbol[iNuc],nucScaleFactor[iNuc]),"lp");
	    else
	      leg->AddEntry(int_graph[ipion][iNuc][iint],Form("^{%d}%s",a_stable[iNuc],element_symbol[iNuc]),"lp");

	  }
	}
	
	isDrawn = true;
	
      }
      
      if (isDrawn) {
	leg->Draw("same");
	print_canvas(c);
	delete leg;
      }
    }
  }

}


void piscatana::print_canvas(TCanvas *c, string outfilename)  {

  string rootname;
  rootname = Form("images/%s%s",c->GetName(),outfilename.c_str());
  
  c->Update();
  //c->Print(Form("%s.pdf",rootname.c_str()));
  //c->Print(Form("%s.eps",rootname.c_str()));
  c->Print(Form("%s.png",rootname.c_str()));
  //cout << "Saved plot: " << Form("%s.png",rootname) << endl;
}



bool piscatana::DrawLegendConfs(vector<TH1F*> h_input,Double_t x0, Double_t y0, Double_t x1, Double_t y1)  {

  if (!legName_parse[0].compare("")) return false;

  TLegend *leg = new TLegend(x0,y0,x1,y1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);
  
  for (int iConf=0; iConf<nConfs; iConf++) {
    leg->AddEntry(h_input[iConf],legName_parse[iConf].c_str(),"l");
  }

  leg->Draw("SAME");
  
  return true;
}




void piscatana::DrawLegendInts(vector<TH1F*> *h_input,Double_t x0, Double_t y0, Double_t x1, Double_t y1, int iNuc)  {

  TLegend *leg = new TLegend(x0,y0,x1,y1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);
  
  int startingIntLabel = rxn;

  if (legName_parse[0].compare("")) {
    for (int iConf=0; iConf<nConfs; iConf++) {
      leg->AddEntry(h_input[startingIntLabel][iConf],Form("%s (%s)",intTitle[startingIntLabel],legName_parse[iConf].c_str()),"l");
    }
    startingIntLabel++;
  }
  
  for (int iint=startingIntLabel; iint<=hadpro; iint++) {
    
    if (iint==dcx || iint==elas) continue;

    // Switch order of abs and QE for light nuclei
    if (iNuc<elem_Nb) {
      if (iint==absb) continue;
      else if (iint==scatt) {
	leg->AddEntry(h_input[scatt][0],Form("%s",intTitle[scatt]),"l");
	leg->AddEntry(h_input[absb][0],Form("%s",intTitle[absb]),"l");
      }
      else  
	leg->AddEntry(h_input[iint][0],Form("%s",intTitle[iint]),"l");
    }
    else  
      leg->AddEntry(h_input[iint][0],Form("%s",intTitle[iint]),"l");
  }

  leg->Draw("SAME");
  
}

