#include <iostream>
#include <cstdlib>

#include "NTotCrs.h"
  
using namespace neut;
using namespace neut::rew;


NTotCrs * NTotCrs::fInstance = 0;
//____________________________________________________________________________
NTotCrs::NTotCrs()
{
}
//____________________________________________________________________________
NTotCrs::~NTotCrs()
{
  fInstance = 0;
}
//____________________________________________________________________________

NTotCrs * NTotCrs::Instance()
{
  if(fInstance == 0) {
    //std::cout << "NTotCrs late initialization" << '\n';
    static NTotCrs::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new NTotCrs;
  
    if (!getenv("NEUT_ROOT")) {
      std::cerr << "NTotCrs::LoadCCQE() Error: Set \"NEUT_ROOT\" environment variable (without the 'src')" << '\n';
      exit (1);
    }
    else {
      char temp[500];
      sprintf(temp,"%s/src/reweight/inputs",getenv("NEUT_ROOT"));
      fInstance->neut_folder = temp;
    }

    fInstance->LoadCCQE();
    fInstance->LoadRESSPI();

  }
  return fInstance;
}
//____________________________________________________________________________


void NTotCrs::LoadCCQE() {
  //TODO: Use TH1::SetDirectory(0) and TFile::GetObject() to replace the craziness and allow us to close the input file

  // Load CCQE cross section histogram
  ccqeCrsFile = new TFile(Form("%s/ccqe_tcrs.root",neut_folder.c_str()));

  if (ccqeCrsFile) {
    if (ccqeCrsFile->IsOpen()) {
      ccqe_crs[numu] = new TH3D(*(TH3D*)ccqeCrsFile->Get("h_ccqe_tcrs_numu"));
      ccqe_crs[numub] = new TH3D(*(TH3D*)ccqeCrsFile->Get("h_ccqe_tcrs_numub"));
      ccqe_crs[nue] = new TH3D(*(TH3D*)ccqeCrsFile->Get("h_ccqe_tcrs_nue"));
      ccqe_crs[nueb] = new TH3D(*(TH3D*)ccqeCrsFile->Get("h_ccqe_tcrs_nueb"));
    } else {
      std::cerr << "NTotCrs::LoadCCQE() Error: File " << Form("%s/ccqe_tcrs.root",neut_folder.c_str()) << " not open" << '\n';
      exit (-1);
    }
  } else {
    std::cerr << "NTotCrs::LoadCCQE() Error: ccqeCrsFile is NULL for file " << Form("%s/ccqe_tcrs.root",neut_folder.c_str()) << " not open" << '\n';
    exit (-1);
  }
}



void NTotCrs::LoadRESSPI() {
  //TODO: Use TH1::SetDirectory(0) and TFile::GetObject() to replace the craziness and allow us to close the input file

  resspiCrsFile = new TFile(Form("%s/resspi_tcrs.root",neut_folder.c_str()));

  if (resspiCrsFile) {
    if (resspiCrsFile->IsOpen()) {
      for (int imode=0; imode<7; imode++) {

	// Anti-neutrino NC modes are not in the same order as the neutrino NC modes
	// so make them the same order in the histogram array (MC generation used the 
	// incorrect cross section, so comment out the following swap and use 
	// the incorrect cross section consistently to match the MC)
	int imode_swap = imode;
	//if (imode == 3) imode_swap = 5;
	//else if (imode == 4) imode_swap = 6;
	//else if (imode == 5) imode_swap = 3;
	//else if (imode == 6) imode_swap = 4;	

	resspi_crs[numu][imode] = new TH2D(*(TH2D*)resspiCrsFile->Get(Form("hresspi_crs_numu_%d_0",imode+1)));
	resspi_crs[numub][imode_swap] = new TH2D(*(TH2D*)resspiCrsFile->Get(Form("hresspi_crs_numub_%d_0",imode+1)));
	resspi_crs[nue][imode] = new TH2D(*(TH2D*)resspiCrsFile->Get(Form("hresspi_crs_nue_%d_0",imode+1)));
	resspi_crs[nueb][imode_swap] = new TH2D(*(TH2D*)resspiCrsFile->Get(Form("hresspi_crs_nueb_%d_0",imode+1)));
      }
    } else {
      std::cerr << "NTotCrs::LoadRESSPI() Error: File " << Form("%s/resspi_tcrs.root",neut_folder.c_str()) << " not open" << '\n';
      exit (-1);
    }
  } else {
    std::cerr << "NTotCrs::LoadRESSPI() Error: resspiCrsFile is NULL for file " << Form("%s/resspi_tcrs.root",neut_folder.c_str()) << " not open" << '\n';
    exit (-1);
  }

}
