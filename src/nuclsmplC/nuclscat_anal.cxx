#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>

#include "TFile.h"
#include "TTree.h"

#include "nuclscat_output_blocks.hxx"

namespace {

  TFile* InputFile = NULL;
  TTree* InputTree = NULL;

  char const * InputFileName;
}

int main(int argc, char const * argv[]){
  if(argc != 2){
    std::cerr << "[ERROR]: Expect 1 CLI argument specifying the file to read."
      << std::endl;
    return 1;
  }
  InputFileName = argv[1];


  InputFile = TFile::Open(InputFileName);
  if(!InputFile || !InputFile->IsOpen()){
    std::cerr << "[ERROR]: Failed to open input file." << std::endl;
    return 2;
  }
  InputTree = dynamic_cast<TTree*>(InputFile->Get("h1"));

  if(!InputTree){
    std::cerr << "[ERROR]: Failed to open input tree." << std::endl;
    return 4;
  }

  if(!VectorCommon_SetBranchAddresses(InputTree) ||
     !FSIHistoryCommon_SetBranchAddresses(InputTree) ||
     !NucleusTargetCommon_SetBranchAddresses(InputTree) ||
     !NucleonFSIHistoryCommon_SetBranchAddresses(InputTree) ){
    std::cerr << "[ERROR]: Failed to set input branch addresses." << std::endl;
    return 8;
  }

  for(Long64_t entry = 0; entry < InputTree->GetEntries(); ++entry){
    std::cout << "[INFO]: Getting entry #" << entry << std::endl;
    InputTree->GetEntry(entry);

    std::cout << "[INFO]: Copying read variables to FORTRAN Common blocks."
      << std::endl;
    Copy_VectorCommon_C2F();
    Copy_FSIHistoryCommon_C2F();
    Copy_NucleusTargetCommon_C2F();
    Copy_NucleonFSIHistoryCommon_C2F();

    std::cout << "[INFO]: Found " << VectorCommon.NParVec
      << " | " << vcwork->nvc << " entries in this vector." << std::endl;
    for(Int_t i = 3; i < VectorCommon.NParVec; ++i){
      std::cout << "\t[" << i << "]IOrgVec: " << VectorCommon.IOrgVec[i]
        << " | " << vcwork->iorgvc[i] << std::endl;
      std::cout << "\t[" << i << "]IPVec: " << VectorCommon.IPVec[i]
        << " | " << vcwork->ipvc[i] << std::endl;
      std::cout << "\t[" << i << "]IChVec: " << VectorCommon.IChVec[i]
        << " | " << vcwork->icrnvc[i] << std::endl;
      std::cout << "\t[" << i << "]IFlgVec: " << VectorCommon.IFlgVec[i]
        << " | " << vcwork->iflgvc[i] << std::endl;

      std::cout << "\t[" << i << "]PVec[0]: " << VectorCommon.PVec[i][0]
        << " | " << vcwork->pvc[i][0] << std::endl;
      std::cout << "\t[" << i << "]PVec[1]: " << VectorCommon.PVec[i][1]
        << " | " << vcwork->pvc[i][1] << std::endl;
      std::cout << "\t[" << i << "]PVec[2]: " << VectorCommon.PVec[i][2]
        << " | " << vcwork->pvc[i][2] << std::endl;
    }

    std::cout << "[INFO]: Target description:" << std::endl;
    std::cout << "\tNumBndN: " << NucleusTargetCommon.NumBndN << " | "
      << neuttarget->numbndn << std::endl;;
    std::cout << "\tNumBndP: " << NucleusTargetCommon.NumBndP << " | "
      << neuttarget->numbndp << std::endl;;
    std::cout << "\tNumFreP: " << NucleusTargetCommon.NumFreP << " | "
      << neuttarget->numfrep << std::endl;;
    std::cout << "\tNumAtom: " << NucleusTargetCommon.NumAtom << " | "
      << neuttarget->numatom << std::endl;;

    std::cout << "[INFO]: Found " << FSIHistoryCommon.NVert << " | "
      << fsihist->nvert << " FSI vertices and " << FSIHistoryCommon.NVCVert
      << " | " << fsihist->nvcvert << " VC FSI vertices in this vector."
      << std::endl;

    for(Int_t IVert = 0; (IVert < fsihist->nvert) && (IVert < MaxVert);
      ++IVert){

      std::cout << "\t[" << IVert << "]IFlgVert: "
      << FSIHistoryCommon.IFlgVert[IVert] << " | "
        << fsihist->iflgvert[IVert] << std::endl;
      std::cout << "\t[" << IVert << "]PosVert[0]: "
        << FSIHistoryCommon.PosVert[IVert][0] << " | "
        << fsihist->posvert[IVert][0] << std::endl;
      std::cout << "\t[" << IVert << "]PosVert[1]: "
        << FSIHistoryCommon.PosVert[IVert][1] << " | "
        << fsihist->posvert[IVert][1] << std::endl;
      std::cout << "\t[" << IVert << "]PosVert[2]: "
        << FSIHistoryCommon.PosVert[IVert][2] << " | "
        << fsihist->posvert[IVert][2] << std::endl;
    }

    for(Int_t IVCVert = 0;
      (IVCVert < fsihist->nvcvert) && (IVCVert < MaxVCVert); ++IVCVert){

      std::cout << "\t[" << IVCVert << "]IPVert: "
      << FSIHistoryCommon.IPVert[IVCVert] << " | "
        << fsihist->ipvert[IVCVert] << std::endl;
      std::cout << "\t[" << IVCVert << "]IVertI: "
        << FSIHistoryCommon.IVertI[IVCVert] << " | "
        << fsihist->iverti[IVCVert] << std::endl;
      std::cout << "\t[" << IVCVert << "]IVertF: "
        << FSIHistoryCommon.IVertF[IVCVert] << " | "
        << fsihist->ivertf[IVCVert] << std::endl;
      std::cout << "\t[" << IVCVert << "]DirVert[0]: "
        << FSIHistoryCommon.DirVert[IVCVert][0] << " | "
        << fsihist->dirvert[IVCVert][0] << std::endl;
      std::cout << "\t[" << IVCVert << "]DirVert[1]: "
        << FSIHistoryCommon.DirVert[IVCVert][1] << " | "
        << fsihist->dirvert[IVCVert][1] << std::endl;
      std::cout << "\t[" << IVCVert << "]DirVert[2]: "
        << FSIHistoryCommon.DirVert[IVCVert][2] << " | "
        << fsihist->dirvert[IVCVert][2] << std::endl;
      std::cout << "\t[" << IVCVert << "]AbsPVert: "
        << FSIHistoryCommon.AbsPVert[IVCVert] << " | "
        << fsihist->abspvert[IVCVert] << std::endl;
      std::cout << "\t[" << IVCVert << "]AbsTPVert: "
        << FSIHistoryCommon.AbsTPVert[IVCVert] << " | "
        << fsihist->abstpvert[IVCVert] << std::endl;
    }
    std::cout << "\tFSIProb: "
      << FSIHistoryCommon.FSIProb << " | "
      << fsihist->fsiprob << std::endl;

    std::cout << "[INFO]: Found " << NucleonFSIHistoryCommon.NFNVert << " | "
      << nucleonfsihist->nfnvert << " nucleon FSI vertices and "
      << NucleonFSIHistoryCommon.NFNStep << " | " << nucleonfsihist->nfnstep
      << " nucleon FSI steps in this vector."
      << std::endl;

    for(Int_t IFNVert = 0;
      (IFNVert < nucleonfsihist->nfnvert) && (IFNVert < MaxNucleonVert);
      ++IFNVert){
      std::cout << "\t[" << IFNVert << "]NFIFlag: "
        << NucleonFSIHistoryCommon.NFIFlag[IFNVert] << " | "
        << nucleonfsihist->nfiflag[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFFirstStep: "
        << NucleonFSIHistoryCommon.NFFirstStep[IFNVert] << " | "
        << nucleonfsihist->nffirststep[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFx: "
        << NucleonFSIHistoryCommon.NFx[IFNVert] << " | "
        << nucleonfsihist->nfx[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFy: "
        << NucleonFSIHistoryCommon.NFy[IFNVert] << " | "
        << nucleonfsihist->nfy[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFz: "
        << NucleonFSIHistoryCommon.NFz[IFNVert] << " | "
        << nucleonfsihist->nfz[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFpx: "
        << NucleonFSIHistoryCommon.NFpx[IFNVert] << " | "
        << nucleonfsihist->nfpx[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFpy: "
        << NucleonFSIHistoryCommon.NFpy[IFNVert] << " | "
        << nucleonfsihist->nfpy[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFpz: "
        << NucleonFSIHistoryCommon.NFpz[IFNVert] << " | "
        << nucleonfsihist->nfpz[IFNVert] << std::endl;
      std::cout << "\t[" << IFNVert << "]NFe: "
        << NucleonFSIHistoryCommon.NFe[IFNVert] << " | "
        << nucleonfsihist->nfe[IFNVert] << std::endl;
    }

    for(Int_t IFNStep = 0;
      (IFNStep < nucleonfsihist->nfnstep) && (IFNStep < MaxNucleonStep);
      ++IFNStep){
      std::cout << "\t[" << IFNStep << "]NFptot: "
        << NucleonFSIHistoryCommon.NFptot[IFNStep] << " | "
        << nucleonfsihist->nfptot[IFNStep] << std::endl;
      std::cout << "\t[" << IFNStep << "]NFecms2: "
        << NucleonFSIHistoryCommon.NFecms2[IFNStep] << " | "
        << nucleonfsihist->nfecms2[IFNStep] << std::endl;
    }
  }

  InputFile->Close();
  return 0;
}
