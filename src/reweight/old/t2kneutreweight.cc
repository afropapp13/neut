#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "necardC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"
#include "neworkC.h"
#include "fsihistC.h"

#include "t2kneutreweight.h"

using namespace std;



// ROOT file input constructor
t2kneutreweight::t2kneutreweight(TString aTheFile, Int_t aMCID, Int_t aNParSets, Double_t **apars, int averbose)
{
  if (!aNParSets) {
    std::cerr << "t2kneutreweight Error: Initialized with nParSets=0" << std::endl;
    exit (-1);
  } else {
    nParSets = aNParSets;
  }
  
  if (!apars) {
    std::cerr << "t2kneutreweight Error: Initialized with no parameter array" << std::endl;
    exit (-1);
  } else {
    pars = apars;
  }
  
  MCID = aMCID;

  // Developer's warnings
  if (MCID==piscat) {
    std::cout << "t2kneutreweight Warning: piscat nucleus is hardcoded as oxygen" << std::endl;
  }
  
  theFile = aTheFile;
  
  //copyFile();
  openFile();
  makeParTree();

  verbose = averbose;

  Init();
}


t2kneutreweight::~t2kneutreweight()
{
  if (!fNeutChain) return;
  delete fNeutChain->GetCurrentFile();
  
}



void t2kneutreweight::appendWeightsToTree() {

  if (fNeutChain == 0) return;

  // Load input tree
  //TTree *outTree = (TTree*)mcRwOutFile->Get(treeKey);
  TTree *outTree = new TTree("NEUTrw","NuInt and Pion FSI Weights");

  // Weights variables
  const int kParSets = nParSets;
  Float_t NEdifcrs[kParSets];
  Float_t NEfsiprob[kParSets];

  // Assume default parameter set is always the first one
  const int defset=0;


  // Initialize
  for (int iparset=0; iparset<kParSets; iparset++) {
    NEdifcrs[iparset] = -1;
    NEfsiprob[iparset] = -1;
  }

  // New Branches
  TBranch *b_NEdifcrs;
  if (MCID != piscat)
    b_NEdifcrs = outTree->Branch("NEdifcrs", NEdifcrs,Form("NEdifcrs[%d]/F",nParSets));
  
  TBranch *b_NEfsiprob = outTree->Branch("NEfsiprob", NEfsiprob,Form("NEfsiprob[%d]/F",nParSets));
  
  // Loop through events
  Long64_t nentries = fNeutChain->GetEntriesFast();
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fNeutChain->GetEntry(jentry);
    outTree->GetEntry(jentry);

    if (verbose) {
      cout << "------------------------------------" << endl;
      //if (ientry%1000 == 0)
	cout << "Event  " << ientry << endl;
    }
    
    if (MCID != piscat) fill_nework();
    fill_fsihist();
    
    //if (!Nvert) continue;
    

    if (verbose) {
      print_nework();
      print_fsihist();
    }


    // Loop through parameter sets
    for (Int_t iparset=0; iparset<nParSets; iparset++) {
      
      fill_allparams(iparset);
      if (verbose==2) 	{
	cout << "Parameter Set " << iparset;      
	print_allparams();
      }

      // Take care of zeros and nan (but probably should fix in evdifcrs_ !!
      if (iparset!=defset) {

	if (NEdifcrs[defset]<=0) NEdifcrs[iparset]=NEdifcrs[defset];
	else if (MCID!=piscat) NEdifcrs[iparset] = evdifcrs_();

	if (NEfsiprob[defset]<=0) NEfsiprob[iparset]=NEfsiprob[defset];
	else {
	  NEfsiprob[iparset] = evpiprob_();
	}

      } else { // Default set
	if (MCID!=piscat) NEdifcrs[iparset] = evdifcrs_();
	
	NEfsiprob[iparset] = evpiprob_();
      }
      
      if (isinf(NEdifcrs[iparset]) || isnan(NEdifcrs[iparset])) 
	NEdifcrs[iparset] = -1;

      if (isinf(NEfsiprob[iparset]) || isnan(NEfsiprob[iparset])) 
	NEfsiprob[iparset] = -1;

      
      // Return weight instead of absolute number
      if (iparset!=defset) {
	if (NEdifcrs[defset]>0) {
	  NEdifcrs[iparset] = NEdifcrs[iparset]/NEdifcrs[defset];
    
	  // Check (this shouldn't happen)
	  if (isinf(NEdifcrs[iparset]) || isnan(NEdifcrs[iparset])) {
	    cerr << "t2kneutreweight Warning: NEdifcrs is infinite, setting to 1" << endl;
	    NEdifcrs[iparset] = 1;
	  }
	}
	else NEdifcrs[iparset] = 1;  // Default set

	if (NEfsiprob[defset]>0) {
	  NEfsiprob[iparset] = NEfsiprob[iparset]/NEfsiprob[defset];
    
	  // Check (this shouldn't happen)
	  if (isinf(NEfsiprob[iparset]) || isnan(NEfsiprob[iparset])) {
	    cerr << "t2kneutreweight Warning: NEfsiprob[iparset] is infinite, setting to 1" << endl;
	    NEfsiprob[iparset] = 1;
	  }
	}
	else NEfsiprob[iparset] = 1;
      }
      
      
      if (verbose) {
	cout << "Event " << ientry << ": EVDIFCRS = " << NEdifcrs[iparset] << ", EVPIPROB = " << NEfsiprob[iparset] << endl << endl;
      }

    }
   
    if (verbose) 
      cout << endl << endl;
    
    // Do not weight default set
    NEdifcrs[defset] = 1;
    NEfsiprob[defset] = 1;
    

    outTree->Fill();
    

  }

  outTree->Write();

}



void t2kneutreweight::fill_allparams(Int_t iparset) {
  fill_necard(iparset);
  fill_neutparams(iparset);
  fill_neutmodel(iparset);
    
}

void t2kneutreweight::fill_allevent(Int_t ivtx) {
  fill_nework(ivtx);
  fill_fsihist(ivtx);
}

void t2kneutreweight::fill_neutmodel(Int_t iparset) {
  // Documentation: See neutmodel.h and necard.h 

  nemdls_.xmaqe       =  pars[iparset][QEMA];
  nemdls_.xmvqe       =  pars[iparset][QEMV];
  nemdls_.kapp        =  pars[iparset][QEKAPP];

  nemdls_.xmaspi      =  pars[iparset][SPIMA];
  nemdls_.xmvspi      =  pars[iparset][SPIMV];

  nemdls_.xmacoh      =  pars[iparset][COHMA];

  neutdis_.nepdf      =  pars[iparset][DISPDF];
  neutdis_.nebodek    =  pars[iparset][DISBOD];

  neutcoh_.necohepi = 0; // To be implemented

  nefillmodel_(); // copy parameters to neutmodel common block
}

void t2kneutreweight::fill_necard(Int_t iparset) {
  // Documentation: See necard.h

  // Get target nucleus type 

  if (MCID==sk || MCID==piscat) {
    neuttarget_.numbndn =  8;
    neuttarget_.numbndp =  8;
    neuttarget_.numfrep =  2;
    neuttarget_.numatom = 16; 
    nesetfgparams_(); 
  }
  else if (MCID==nd280numc) {
    
   // PDG nuclear code format
    int nucID = abs(StdHepPdg[1]) - 1000000000;
    int nStrange = nucID/10000000;
    int Z = (nucID - nStrange*10000000)/10000;
    int A = (nucID - nStrange*10000000 - Z*10000)/10;
    int I = (nucID - nStrange*10000000 - Z*10000 - A*10);

    if (Z==1) {
      neuttarget_.numatom = 1; 
      neuttarget_.numbndp = 0;
      neuttarget_.numbndn = 0;
      neuttarget_.numfrep = 1; 
    }
    else if (Z>1){
      neuttarget_.numatom = A; 
      neuttarget_.numbndp = Z;
      neuttarget_.numbndn = A-Z;
      neuttarget_.numfrep = 0; 
      nesetfgparams_(); 
    }
    else {
      std::cerr << "t2kneutreweight Error: Unknown Z = " << Z << std::endl;
      exit (-1);
    }

    
  } else {
    
    std::cerr << "t2kneutreweight Error: MCID = " << MCID << " nucleus filling not yet implemented" << std::endl;
    exit (-1);

  }


  // FG parameter variations should go here in fractional form
  // but MC should have simulated the required phase space

  // Default paramaters
  neutcard_.nefrmflg  =  0;
  neutcard_.nepauflg  =  0;
  neutcard_.nenefo16  =  0;
  neutcard_.nemodflg  =  0;
  neutcard_.nenefmodl = 1;
  neutcard_.nenefmodh = 1;
  neutcard_.nenefkinh = 1;
  neutpiabs_.neabspiemit = 1;

  if (MCID==piscat)
    neutcard_.nusim = 0;
  else
    neutcard_.nusim = 1;
  
}   

void t2kneutreweight::fill_neutparams(Int_t iparset) {

  // Nuclear parameter //
  // Note: Nucleus dependant so set using nesetfgparams_() after nucleus is determined
  //nenupr_.pfsurf      = 0.225;
  //nenupr_.pfmax       = 0.225;
  //nenupr_.vnuini      = -1. * (sqrt(0.939 * 0.939 + 0.225 * 0.225) - 0.939);
  //nenupr_.vnufin      =  0;
  nenupr_.iformlen    =  1;  // On by default

  // Pion FSI Parameters
  neffpr_.fefqe   = pars[iparset][FSIQE];
  neffpr_.fefqeh  = pars[iparset][FSIQEH];
  neffpr_.fefinel = pars[iparset][FSIINEL];
  neffpr_.fefabs  = pars[iparset][FSIABS];
  neffpr_.fefcx   = pars[iparset][FSICX];
  neffpr_.fefcxh  = pars[iparset][FSICXH];
  
  // Unvaried values set to default
  neffpr_.fefcoh = 1;
  neffpr_.fefqehf = 1;
  neffpr_.fefcxhf = 0;
  neffpr_.fefcohf = 0;
}



void t2kneutreweight::fill_nework(Int_t ivtx)
{
  // Documentation: See nework.h


  // Common to all MC
  nework_.numne  = 4; // Only need first 4 particles 


  if (MCID==sk || MCID==piscat || nd280numc) {

    if (MCID==sk)
      nework_.modene = Neutmode;
    else if (MCID==nd280numc)
      nework_.modene = atoi(fString.Data());
    
    for ( int i = 0; i<nework_.numne; i++ ) { 
      nework_.ipne[i] = Ipvc[i];
      for ( int j = 0 ; j < 3 ; j++ ){
	nework_.pne[i][j] = (float)Pvc[i][j]/1000;  // VC(NE)WORK in M(G)eV
      }
    }
  } 


  else if (MCID==nd280anal) {
    nework_.modene = Vtx_NEneutmode[ivtx];
    
    for ( int i = 0; i<nework_.numne; i++ ) {  // Only need first 4 particles 
      nework_.ipne[i] = Vtx_NEipvc[ivtx][i];
      for ( int j = 0 ; j < 3 ; j++ ){
	nework_.pne[i][j] = (float)Vtx_NEpvc[ivtx][i][j]/1000;  // VC(NE)WORK in M(G)eV
      }
    }
  }

}

void t2kneutreweight::fill_fsihist(Int_t ivtx)
{

  if (MCID==sk || MCID==nd280numc || MCID==piscat) {
    fsihist_.nvert = Nvert;
    for (int ivert=0; ivert<Nvert; ivert++) {
      fsihist_.iflgvert[ivert] = Iflgvert[ivert];
      for (int j=0; j<3; j++)
	fsihist_.posvert[ivert][j] = (float)Posvert[ivert][j];
    }
  
    fsihist_.nvcvert = Nvcvert;
    for (int ip=0; ip<Nvcvert; ip++) {
      fsihist_.abspvert[ip] = (float)Abspvert[ip];
      fsihist_.abstpvert[ip] = (float)Abstpvert[ip];
      fsihist_.ipvert[ip] = Ipvert[ip];
      fsihist_.iverti[ip] = Iverti[ip];
      fsihist_.ivertf[ip] = Ivertf[ip];
      for (int j=0; j<3; j++)
	fsihist_.dirvert[ip][j] = (float)Dirvert[ip][j];
    }
  }
  
  else if (MCID==nd280anal) {
    fsihist_.nvert = Vtx_NEnvert[ivtx];
    for (int ivert=0; ivert<Vtx_NEnvert[ivtx]; ivert++) {
      fsihist_.iflgvert[ivert] = Vtx_NEiflgvert[ivtx][ivert];
      for (int j=0; j<3; j++)
	fsihist_.posvert[ivert][j] = (float)Vtx_NEposvert[ivtx][ivert][j];
    }
  
    fsihist_.nvcvert = Vtx_NEnvcvert[ivtx];
    for (int ip=0; ip<Vtx_NEnvcvert[ivtx]; ip++) {
      fsihist_.abspvert[ip] = (float)Vtx_NEabspvert[ivtx][ip];
      fsihist_.abstpvert[ip] = (float)Vtx_NEabstpvert[ivtx][ip];
      fsihist_.ipvert[ip] = Vtx_NEipvert[ivtx][ip];
      fsihist_.iverti[ip] = Vtx_NEiverti[ivtx][ip];
      fsihist_.ivertf[ip] = Vtx_NEivertf[ivtx][ip];
      for (int j=0; j<3; j++)
	fsihist_.dirvert[ip][j] = (float)Vtx_NEdirvert[ivtx][ip][j];
    }
  }
  
}


void t2kneutreweight::print_fsihist() {
  
  cout << endl << "-------------------------------" << endl;
  cout << "FSI History Common Block" << endl;
  cout << endl << "Nvert: " << fsihist_.nvert << endl;
  cout << "ivert  iflgvert posvert[3]" << endl;
  for (int ivert=0; ivert<fsihist_.nvert; ivert++) {
    cout << fsihist_.iflgvert[ivert] << " ";
    for (int j=0; j<3; j++)
      cout << fsihist_.posvert[ivert][j] << " ";
    cout << endl;
  }
    
  
  cout << endl << "Nvcvert: " << fsihist_.nvcvert << endl;
  cout << "ip ipvert abspvert iverti ivertf dirvert[3]" << endl;
  for (int ip=0; ip<fsihist_.nvcvert; ip++) {
    cout << fsihist_.ipvert[ip] << " " << fsihist_.abspvert[ip] << " " << fsihist_.iverti[ip] << " " << fsihist_.ivertf[ip] << " ";
    for (int j=0; j<3; j++)
      cout << fsihist_.dirvert[ip][j] << " ";
    cout << endl;
  }
 
  cout << endl;
}


void t2kneutreweight::print_nework() {
  
  cout << endl << "-------------------------------" << endl;
  cout << "NEWORK Common Block" << endl;
  cout << endl << "Mode = " << nework_.modene << ", Numne = " << nework_.numne << endl;
  cout << "ipne pne[3]" << endl;
  for (int ip=0; ip<nework_.numne; ip++) {
    cout << nework_.ipne[ip] << " ";
    for (int j=0; j<3; j++)
      cout << nework_.pne[ip][j] << " ";
    cout << endl;
  }
  
  cout << endl;
}

void t2kneutreweight::print_allevent() {
  print_nework();
  print_fsihist();
}


void t2kneutreweight::print_allparams() {

  cout << endl << endl;

  // Input parameters
  cout << endl << "-- User input variable parameters --" << endl;

  cout << "nemdls_.xmaqe = " <<         nemdls_.xmaqe   << endl;    
  cout << "nemdls_.xmvqe = " <<        nemdls_.xmvqe     << endl;  
  cout << "nemdls_.kapp = " <<          nemdls_.kapp      << endl;  
  cout << "nemdls_.xmaspi = " <<      nemdls_.xmaspi     << endl;
  cout << "nemdls_.xmvspi = " <<      nemdls_.xmvspi   << endl;
  cout << "nemdls_.xmacoh = " <<         nemdls_.xmacoh   << endl;    
  cout << "neutdis_.nepdf = " <<       neutdis_.nepdf   << endl;  
  cout << "neutdis_.nebodek = " <<     neutdis_.nebodek  << endl;  
  cout << "neffpr_.fefqe = " <<     neffpr_.fefqe   << endl;
  cout << "neffpr_.fefqeh = " <<    neffpr_.fefqeh  << endl;
  cout << "neffpr_.fefinel = " <<   neffpr_.fefinel << endl;
  cout << "neffpr_.fefabs = " <<    neffpr_.fefabs  << endl;
  cout << "neffpr_.fefcx = " <<     neffpr_.fefcx   << endl;
  cout << "neffpr_.fefcxh = " <<    neffpr_.fefcxh  << endl;

  // Fixed parameters
  cout << endl << "-- These depend on target nucleus --" << endl;
  
  cout << "neuttarget_.numbndn = " << neuttarget_.numbndn << endl;
  cout << "neuttarget_.numbndp = " <<   neuttarget_.numbndp << endl;
  cout << "neuttarget_.numfrep = " <<  neuttarget_.numfrep<< endl;
  cout << "neuttarget_.numatom = " <<  neuttarget_.numatom<< endl;

  cout << "nenupr_.pfsurf = " <<  nenupr_.pfsurf<< endl;
  cout << "nenupr_.pfmax = " <<  nenupr_.pfmax<< endl;
  cout << "nenupr_.vnuini = " << nenupr_.vnuini<< endl;
  cout << "nenupr_.vnufin = " << nenupr_.vnufin<< endl;
  

  cout << endl << "-- These are unvariable/default parameters [should = ()] --" << endl;
  cout << "nenupr_.iformlen = (1)" <<     nenupr_.iformlen   << endl;
  cout << "neutcard_.nefrmflg = (0)" <<  neutcard_.nefrmflg<< endl; 
  cout << "neutcard_.nepauflg = (0)" <<   neutcard_.nepauflg<< endl; 
  cout << "neutcard_.nenefo16 = (0)" <<  neutcard_.nenefo16 << endl;
  cout << "neutcard_.nemodflg = (0)" <<  neutcard_.nemodflg << endl;
  cout << "neutcard_.nenefmodl = (1)" <<  neutcard_.nenefmodl << endl;
  cout << "neutcard_.nenefmodh = (1)" <<  neutcard_.nenefmodh << endl;
  cout << "neutcard_.nenefkinh = (1)" <<  neutcard_.nenefkinh << endl;
  cout << "neutpiabs_.neabspiemit = (0)" <<  neutpiabs_.neabspiemit << endl;
  cout << "neutcard_.nusim = " <<  neutcard_.nusim << endl;

  cout << "neutcoh_.necohepi  = (0)" <<  neutcoh_.necohepi << endl;

  cout << "neffpr_.fefcoh = (1)" <<   neffpr_.fefcoh << endl;
  cout << "neffpr_.fefqehf = (1)" <<  neffpr_.fefqehf<< endl;
  cout << "neffpr_.fefcxhf = (0)" <<  neffpr_.fefcxhf<< endl;
  cout << "neffpr_.fefcohf  = (0)" << neffpr_.fefcohf<< endl;

  cout << endl;
}

void t2kneutreweight::openFile() {

  TFile *mcRwOutFile = new TFile(theFile,"UPDATE");
  
}

// Has some trouble with piscat ROOT files for some reason
void t2kneutreweight::copyFile()
{
  // Make new file name
  TString newFile = theFile;
  
  Int_t startIndex;
  startIndex = newFile.Index(".root",0);
  newFile.Insert(startIndex,".neutrw");
  
  cout << "Copying to new file: " << newFile << endl << endl;

  // Open input file to make copy of contents
  TFile *inFile = new TFile(theFile);
  TDirectory *dir = gDirectory;
  dir->ReadAll();
  
  // Make new file and copy contents
  mcRwOutFile = new TFile(newFile,"RECREATE");
  cout << dir->GetList()->Write() << endl;

  inFile->Close();
  
}

void t2kneutreweight::makeParTree()
{
  
  // Hack so that ROOT can read array into branch
  const int kParSets = nParSets;
  Double_t pars_temp[kParSets][nPars];
  for (int iparset=0; iparset<nParSets; iparset++) {
    cout << iparset << " ";
    for (int ipar=0; ipar<nPars; ipar++) {
      pars_temp[iparset][ipar] = pars[iparset][ipar];

      cout << pars_temp[iparset][ipar] << " " ;
    }
    cout << endl;
  }

  // Tree containing parameter sets
  TTree *parTree = new TTree("NEUTparams", "Parameter sets for NEUT reweighting");
  parTree->Branch("NEnparsets", &nParSets, "NEnparsets/I");
  parTree->Branch("NEpars", pars_temp,Form("NEpars[NEnparsets][%d]/D",nPars));
  parTree->Fill();
  parTree->Write(); 
  
}



// Initialize necessary branches for reweighting
void t2kneutreweight::Init()
{
  // Select tree name depending on MC type
  if (MCID==sk || MCID==piscat) treeKey = "h1";
  else if (MCID==nd280numc) treeKey = "nRooTracker";
  else if (MCID==nd280anal) treeKey = "TruthDir/NRooTrackerVtx";
  else {
    std::cerr << "Error in t2kneutreweight::Init(): Unknown MC ID = " << MCID << std::endl;
    exit(-1);
  }

  // Load tree (or chain of trees)
  fNeutChain = new TChain(treeKey.Data(),"");
  fNeutChain->Add(Form("%s/%s",theFile.Data(),treeKey.Data()));

  // Set branch addresses and branch pointers
  fCurrent = -1;
  fNeutChain->SetMakeClass(1);
  
  if (MCID==sk || MCID==piscat) {

    // New SK output with NEUT block (read in as floats)
    if (MCID!=piscat)
      fNeutChain->SetBranchAddress("Neutmode", &Neutmode, &b_Neutmode);
    fNeutChain->SetBranchAddress("Npvc", &Npvc, &b_Npvc);
    fNeutChain->SetBranchAddress("Ipvc", Ipvc, &b_Ipvc);
    fNeutChain->SetBranchAddress("Pvc", Pvc, &b_Pvc);
    fNeutChain->SetBranchAddress("Nvert", &Nvert, &b_Nvert);
    fNeutChain->SetBranchAddress("Posvert", Posvert, &b_Posvert);
    fNeutChain->SetBranchAddress("Iflgvert", Iflgvert, &b_Iflgvert);
    fNeutChain->SetBranchAddress("Nvcvert", &Nvcvert, &b_Nvcvert);
    fNeutChain->SetBranchAddress("Dirvert", Dirvert, &b_Dirvert);
    fNeutChain->SetBranchAddress("Abspvert", Abspvert, &b_Abspvert);
    fNeutChain->SetBranchAddress("Abstpvert", Abstpvert, &b_Abstpvert);
    fNeutChain->SetBranchAddress("Ipvert", Ipvert, &b_Ipvert);
    fNeutChain->SetBranchAddress("Iverti", Iverti, &b_Iverti);
    fNeutChain->SetBranchAddress("Ivertf", Ivertf, &b_Ivertf);
    /**/



    // Official 2010a SK output
    /*
    fNeutChain->SetBranchAddress("mode", &Neutmode, &b_Neutmode);
    fNeutChain->SetBranchAddress("numnu", &Npvc, &b_Npvc);
    fNeutChain->SetBranchAddress("ipnu", Ipvc, &b_Ipvc);
    fNeutChain->SetBranchAddress("pnu", fAbspvc, &b_Abspvc);
    fNeutChain->SetBranchAddress("dirnu", fPvc, &b_Pvc);
    /**/

  } else if (MCID==nd280numc) {
    
    // Read in as doubles

    fNeutChain->SetBranchAddress("fString", &fString, &b_Neutmode);
    
    fNeutChain->SetBranchAddress("StdHepN", &Npvc, &b_Npvc);
    fNeutChain->SetBranchAddress("StdHepPdg", StdHepPdg, &b_Ipvc);
    fNeutChain->SetBranchAddress("StdHepP4", StdHepP4, &b_Pvc);
       
    fNeutChain->SetBranchAddress("NEneutmode", &Neutmode, &b_Neutmode);
    fNeutChain->SetBranchAddress("NEnvc", &Npvc, &b_Npvc);
   
    fNeutChain->SetBranchAddress("NEipvc", Ipvc, &b_Ipvc);
    fNeutChain->SetBranchAddress("NEpvc", Pvc, &b_Pvc);
    fNeutChain->SetBranchAddress("NEnvert", &Nvert, &b_Nvert);
    fNeutChain->SetBranchAddress("NEposvert", Posvert, &b_Posvert);
    fNeutChain->SetBranchAddress("NEiflgvert", Iflgvert, &b_Iflgvert);
    fNeutChain->SetBranchAddress("NEnvcvert", &Nvcvert, &b_Nvcvert);
    fNeutChain->SetBranchAddress("NEdirvert", Dirvert, &b_Dirvert);
    fNeutChain->SetBranchAddress("NEabspvert", Abspvert, &b_Abspvert);
    fNeutChain->SetBranchAddress("NEabstpvert", Abstpvert, &b_Abstpvert);
    fNeutChain->SetBranchAddress("NEipvert", Ipvert, &b_Ipvert);
    fNeutChain->SetBranchAddress("NEiverti", Iverti, &b_Iverti);
    fNeutChain->SetBranchAddress("NEivertf", Ivertf, &b_Ivertf);
  } else if (MCID==nd280anal) {

    std::cerr << "t2kneutreweight Error: nd280anal input not yet implemented." << std::endl;
    exit (-1);

    fNeutChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
    fNeutChain->SetBranchAddress("Vtx.NEneutmode", Vtx_NEneutmode, &b_Neutmode);
    fNeutChain->SetBranchAddress("Vtx.NEnvc", Vtx_NEnvc, &b_Npvc);
    fNeutChain->SetBranchAddress("Vtx.NEipvc[100]", Vtx_NEipvc, &b_Ipvc);
    fNeutChain->SetBranchAddress("Vtx.NEpvc[100][3]", Vtx_NEpvc, &b_Pvc);
    fNeutChain->SetBranchAddress("Vtx.NEnvert", Vtx_NEnvert, &b_Nvert);
    fNeutChain->SetBranchAddress("Vtx.NEposvert[50][3]", Vtx_NEposvert, &b_Posvert);
    fNeutChain->SetBranchAddress("Vtx.NEiflgvert[50]", Vtx_NEiflgvert, &b_Iflgvert);
    fNeutChain->SetBranchAddress("Vtx.NEnvcvert", Vtx_NEnvcvert, &b_Nvcvert);
    fNeutChain->SetBranchAddress("Vtx.NEdirvert[300][3]", Vtx_NEdirvert, &b_Dirvert);
    fNeutChain->SetBranchAddress("Vtx.NEabspvert[300]", Vtx_NEabspvert, &b_Abspvert);
    fNeutChain->SetBranchAddress("Vtx.NEabstpvert[300]", Vtx_NEabstpvert, &b_Abstpvert);
    fNeutChain->SetBranchAddress("Vtx.NEipvert[300]", Vtx_NEipvert, &b_Ipvert);
    fNeutChain->SetBranchAddress("Vtx.NEiverti[300]", Vtx_NEiverti, &b_Iverti);
    fNeutChain->SetBranchAddress("Vtx.NEivertf[300]", Vtx_NEivertf, &b_Ivertf);
    
  }
}


Int_t t2kneutreweight::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fNeutChain) return 0;
  return fNeutChain->GetEntry(entry);
}


Long64_t t2kneutreweight::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fNeutChain) return -5;
  Long64_t centry = fNeutChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fNeutChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fNeutChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
  }
  return centry;
}
