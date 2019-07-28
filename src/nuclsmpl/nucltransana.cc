#define nucltransana_cxx

//#define VERBOSE

#include <sstream>
#include <sys/stat.h> 

#include <TMath.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>

#include "nucltransana.h"

#include "../neutgeom/PeriodicTable.h"


#include "necardC.h"
#include "efpionC.h"

using namespace std;

// Analyse Tree
nucltransana::nucltransana(string stringInput1, string stringInput2)
{

  theFile = stringInput1; // Input file(s)
  outFile = stringInput2; // Output file
  InitFromTree();
 
}

void nucltransana::InitFromTree() {

  for (int inucleon=0; inucleon<nNucs; inucleon++) {
    for (int iNucleus=0; iNucleus<kZmax; iNucleus++) {
	  histosLoaded[inucleon][iNucleus]=0;
    }
  }

  cout << endl << "Loading piscat tree from file(s): " << theFile.c_str() << endl;
  
  fChain = new TChain("neut_tree","");
  fChain->Add(Form("%s/neut_tree",theFile.c_str()));

  // Set branch addresses and branch pointers
  fCurrent = -1;
  fChain->SetMakeClass(1);

   fChain->SetBranchAddress("vcwork", &vcwork_nvc, &b_vcwork);
   fChain->SetBranchAddress("nework", &nework_modene, &b_nework);
   fChain->SetBranchAddress("posinnuc", &posinnuc_ibound, &b_posinnuc);
   fChain->SetBranchAddress("target", &target_numbndn, &b_target);
  
  // Set initial pion type and target nucleus
  fChain->GetEntry(0);
  SetInitNucleon();
  SetNucleus();

  outROOT = new TFile(outFile.c_str(),"RECREATE"); 
  
  Init_histo();
}


// Load already analysed histograms
nucltransana::nucltransana(string stringInput1, string stringInput2, int ExpNumEvents)
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

}

nucltransana::~nucltransana()
{
   if (!fChain) return;
   
   outROOT->Write();

   cout << endl << "Wrote histograms to: " << outFile.c_str() << endl << endl;

  for (int inucleon=0; inucleon<nNucs; inucleon++) {
    for (int iNucleus=0; iNucleus<kZmax; iNucleus++) {
	  if (histosLoaded[inucleon][iNucleus]){
		delete histosLoaded[inucleon][iNucleus];
	  }
    }
  }


   delete fChain->GetCurrentFile();
}

Int_t nucltransana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


Long64_t nucltransana::LoadTree(Long64_t entry)
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


void nucltransana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
 
void nucltransana::SetInitNucleon() {

  if (vcwork_ipvc[3] == 2212) {
	init_nuc=0;
  } else if (vcwork_ipvc[3] == 2112) {
	init_nuc=1;
  } else {
    cout << 
	  "nucltransana::SetInitNucleon() Error: Unknown pion PID: " << 
	  vcwork_ipvc[3] << 
	  endl;
    exit (-1);
  }

  cout << endl << "Incident nucleon type: " << nucName[init_nuc] << endl;
}

void nucltransana::SetNucleus()
{
  neuttarget_.numbndn = target_numbndn;
  neuttarget_.numbndp = target_numbndp;
  neuttarget_.numfrep = target_numfrep;
  neuttarget_.numatom = target_numbndn+target_numbndp;

  // Set nucleus and get parameters for normalization
  nesetfgparams_();  // Check of existence in NEUT contained

  cout << endl << "Target nucleus: " << element_symbol[target_numbndp-1] << endl;

  double rMax = 2*eftarget_.c;
  targetArea = TMath::Pi()*rMax*rMax*10;  // mb
}

// Must be same naming convention as in Init_histo()
void nucltransana::InitFromHistos(vector<string> inputfile_arr, int ExpNumEvents) {
  
  if ( ExpNumEvents <= 0 )
    cout << endl << 
	  "nucltransana::InitFromHistos()" <<
	  "Warning: Input ExpNumEvents <= 0, " <<
	  "differential histograms may not be normalized properly" << endl;

  cout << endl << "Loading nucltrans histos from file(s): " << endl;

  // Initialize flag for histos loaded 
  for (int inucleon=0; inucleon<nNucs; inucleon++) {
    for (int iNucleus=0; iNucleus<kZmax; iNucleus++) {
      histosLoaded[inucleon][iNucleus] = new bool[nConfs];
    }
  }

  // Set Normalization for old files

  TFile *infile[nConfs];

  for (int iconf=0; iconf<nConfs; iconf++) {

    infile[iconf] = new TFile(inputfile_arr[iconf].c_str());

    cout << inputfile_arr[iconf].c_str() << endl;

    // Broken or misnamed file
    if (!infile[iconf]->IsOpen()) exit (-1);

    for (int inucleon=0; inucleon<nNucs; inucleon++) {
      for (int iNucleus=0; iNucleus<kZmax; iNucleus++) {

		char *charName;
 
	// Initial Spectrum
		charName = 
		  Form("total_init_%s_%s",nucName[inucleon],element_symbol[iNucleus]);

	// Check if nucleon type and nucleus type combination exists
		TH1F *h_temp = (TH1F*)infile[iconf]->Get(charName);
		if (!h_temp) {
		  histosLoaded[inucleon][iNucleus][iconf] = false;
		  continue;
		}

		// Flag for given pion-nucleon-conf combination loaded
		histosLoaded[inucleon][iNucleus][iconf] = true;
	
		piInitSpec[inucleon][iNucleus].push_back(new TH1F(*(TH1F*)infile[iconf]->Get(charName)));
	
		// Deduce number of files by using number of triggers
		float nFiles = 1;
		if ( ExpNumEvents > 0 )
		  nFiles = piInitSpec[inucleon][iNucleus][iconf]->GetEntries() / ExpNumEvents;
		
		// Interaction separated histograms
		for (int iint=0; iint<nInts; iint++) {
		  charName = Form("%s_%s_%s",
						  intName[iint],
						  nucName[inucleon],
						  element_symbol[iNucleus]);
	  
		  piInitIntSpec[inucleon][iNucleus][iint].push_back(new TH1F(*(TH1F*)infile[iconf]->Get(charName)));
		  
		  piInitIntSpec[inucleon][iNucleus][iint][iconf]->Smooth();
		  
		}
	  }
	}
  }
}

void nucltransana::Init_histo()
{
  
  cout << endl << "Initializing nucltrans MC histograms" << endl;

  // Define Histogram range (must be modified depending on range sample was generated with)
  float dE = 10; // MeV/c
  float s_xMin = 0;
  float s_xMax = 5000; 
  int s_nBins = (s_xMax-s_xMin)/dE;


  char *charName,*charTitle,*charLabelX,*charLabelY;
 
  // Initial Spectrum
  charName = Form("total_init_%s_%s",
				  nucName[init_nuc],element_symbol[target_numbndp-1]);
  charTitle = Form("Total");
  charLabelX = Form("%s Initial Momentum (MeV/c)",nucTitle[init_nuc]);
    
  piInitSpec[init_nuc][target_numbndp-1].push_back(new TH1F(charName,charTitle, s_nBins, s_xMin, s_xMax));
  piInitSpec[init_nuc][target_numbndp-1][0]->GetXaxis()->SetTitle(charLabelX);
  piInitSpec[init_nuc][target_numbndp-1][0]->Sumw2();


  // Interaction separated histograms
  for (int iint=0; iint<nInts; iint++) {
    charName = Form("%s_%s_%s",
					intName[iint],
					nucName[init_nuc],element_symbol[target_numbndp-1]);
    charTitle = Form("%s",intTitle[iint]);
    charLabelX = Form("%s Initial Momentum (MeV/c)",nucTitle[init_nuc]);
    charLabelY = "#sigma (mb)";
    piInitIntSpec[init_nuc][target_numbndp-1][iint].push_back(new TH1F(charName,Form("%s;%s;%s",charTitle, charLabelX, charLabelY), s_nBins, s_xMin, s_xMax));
    piInitIntSpec[init_nuc][target_numbndp-1][iint][0]->Sumw2();
  }

}

// Initial Interaction Spectrum
void nucltransana::int_xsec()
{

  Float_t evt_weight;
  evt_weight = 1;

  // Survival / Normal
  if (vcwork_iflgvc[3]==0) {
	/* No reaction */
    piInitIntSpec[init_nuc][target_numbndp-1][nInts-1][0]->Fill(abspvc,evt_weight);
  } else {
	/* something happened */
	piInitIntSpec[init_nuc][target_numbndp-1][0][0]->Fill(abspvc,evt_weight);

    // Determine final state
	Int_t inelastic = 0;

    for (int ipvc=4; ipvc<vcwork_nvc; ipvc++) {
      
      if (vcwork_icrnvc[ipvc]==0) continue;
      
      if ((vcwork_ipvc[ipvc] != 2112)&&(vcwork_ipvc[ipvc] != 2212)){
		inelastic = 1;
	  }
    }
        
    if (inelastic==0) {
	  /* no pions produced */
      piInitIntSpec[init_nuc][target_numbndp-1][1][0]->Fill(abspvc,evt_weight);
    }else{
	  piInitIntSpec[init_nuc][target_numbndp-1][2][0]->Fill(abspvc,evt_weight);
	};

  }
}

void nucltransana::Loop()
{
  if (fChain == 0) return;
  
  // Loop through events
  nentries = fChain->GetEntriesFast();
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%1000==0) cout << "Event: " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    
    // Total Initial Spectrum
    piInitSpec[init_nuc][target_numbndp-1][0]->Fill(abspvc);

	abspvc = sqrt(
				  vcwork_pvc[3][0]*vcwork_pvc[3][0] +
				  vcwork_pvc[3][1]*vcwork_pvc[3][1] +
				  vcwork_pvc[3][2]*vcwork_pvc[3][2]);
	Float_t kE_nuc;
	Float_t xmassnuc=998.272;
  
	abspvc = sqrt((abspvc*abspvc)+xmassnuc*xmassnuc)-xmassnuc;
	
    // Analysis
    int_xsec();
  }  // End event loop


  // Normalize
  for (int iint=0; iint<nInts-1; iint++) {
    piInitIntSpec[init_nuc][target_numbndp-1][iint][0]
	  ->Divide(piInitSpec[init_nuc][target_numbndp-1][0]);
    piInitIntSpec[init_nuc][target_numbndp-1][iint][0]->Scale(targetArea);
  }
  /* Last one is transparency : Don't normalize with area*/
  piInitIntSpec[init_nuc][target_numbndp-1][nInts-1][0]
	->Divide(piInitSpec[init_nuc][target_numbndp-1][0]);

}
