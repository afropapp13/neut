#include "TTree.h"
#include "Rtypes.h"

//used for MaxVC macro definition
#include "vcworkC.h"
//used for MAXNE macro definition
#include "neworkC.h"
//used for MAXVERT and MAXVCVERT macro definitions
#include "fsihistC.h"
//used for MAXNUCLEONVERT and MAXNUCLEONSTEP macro definitions
#include "nucleonfsihistC.h"
//Provides neuttarget common block
#include "necardC.h"

namespace {
  const Int_t MaxNE = MAXNE;
  const Int_t MaxVC = MAXVC;
  const Int_t MaxVert = MAXVERT;
  const Int_t MaxVCVert = MAXVCVERT;
  const Int_t MaxNucleonVert = MAXNUCLEONVERT;
  const Int_t MaxNucleonStep = MAXNUCLEONSTEP;

  template<typename T>
  inline std::string tostr(T const &o){
    std::stringstream ss("");
    ss << o;
    return ss.str();
  }
}

template<typename T>
inline void Clear_CommonBlock(T &c){
  memset(&c, 0, sizeof(T));
}

struct NEUTEventCommon_{
  Int_t IEvent;
  Int_t Mode;
  Int_t NParNEUT;
  Int_t IPNEUT[MaxNE];

  Float_t AbsPNEUT[MaxNE];
  Float_t PNEUT[MaxNE][3];
} NEUTEventCommon;

inline Bool_t Add_NEUTEventCommon_Branches(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (tree->Branch("NEUTEventCommon.IEvent", &NEUTEventCommon.IEvent,
    "NEUTEventCommon.IEvent/I") != NULL);

  succ = succ &&
    (tree->Branch("NEUTEventCommon.Mode", &NEUTEventCommon.Mode,
    "NEUTEventCommon.Mode/I") != NULL);

  succ = succ &&
    (tree->Branch("NEUTEventCommon.NParNEUT", &NEUTEventCommon.NParNEUT,
    "NEUTEventCommon.NParNEUT/I") != NULL);

  succ = succ &&
    (tree->Branch("NEUTEventCommon.IPNEUT", NEUTEventCommon.IPNEUT,
    "NEUTEventCommon.IPNEUT[NEUTEventCommon.NParNEUT]/I") != NULL);

  succ = succ &&
    (tree->Branch("NEUTEventCommon.AbsPNEUT", NEUTEventCommon.AbsPNEUT,
    "NEUTEventCommon.AbsPNEUT[NEUTEventCommon.NParNEUT]/F") != NULL);

  succ = succ &&
    (tree->Branch("NEUTEventCommon.PNEUT", NEUTEventCommon.PNEUT,
   ( "NEUTEventCommon.PNEUT["+tostr(MaxNE)+"][3]/F").c_str() ) != NULL);

  return succ;
}

inline Bool_t NEUTEventCommon_SetBranchAddresses(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (!tree->SetBranchAddress("NEUTEventCommon.IEvent",
      &NEUTEventCommon.IEvent));

  succ = succ &&
    (!tree->SetBranchAddress("NEUTEventCommon.Mode",
      &NEUTEventCommon.Mode));

  succ = succ &&
    (!tree->SetBranchAddress("NEUTEventCommon.NParNEUT",
      &NEUTEventCommon.NParNEUT));

  succ = succ &&
    (!tree->SetBranchAddress("NEUTEventCommon.IPNEUT",
      NEUTEventCommon.IPNEUT));

  succ = succ &&
    (!tree->SetBranchAddress("NEUTEventCommon.AbsPNEUT",
      NEUTEventCommon.AbsPNEUT));

  succ = succ &&
    (!tree->SetBranchAddress("NEUTEventCommon.PNEUT",
      NEUTEventCommon.PNEUT));

  return succ;
}

inline void Copy_NEUTEventCommon_F2C(){
  NEUTEventCommon.Mode = nework->modene;
  NEUTEventCommon.NParNEUT = nework->numne;
  for(Int_t INe = 0; INe < (nework->numne) && (INe < MaxNE); ++INe){
    NEUTEventCommon.IPNEUT[INe] = nework->ipne[INe];
    NEUTEventCommon.AbsPNEUT[INe] = sqrt(nework->pne[INe][0]*nework->pne[INe][0]
                                  + nework->pne[INe][1]*nework->pne[INe][1]
                                  + nework->pne[INe][2]*nework->pne[INe][2]);
    NEUTEventCommon.PNEUT[INe][0] = nework->pne[INe][0];
    NEUTEventCommon.PNEUT[INe][1] = nework->pne[INe][1];
    NEUTEventCommon.PNEUT[INe][2] = nework->pne[INe][2];
  }
}

struct VectorCommon_{
  Int_t NParVec;
  Int_t IOrgVec[MaxVC];
  Int_t IPVec[MaxVC];
  Int_t IChVec[MaxVC];
  Int_t IFlgVec[MaxVC];

  Float_t AbsPVec[MaxVC];
  Float_t PVec[MaxVC][3];
  Float_t PosV[3];
} VectorCommon;

inline Bool_t Add_VectorCommon_Branches(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (tree->Branch("VectorCommon.NParVec", &VectorCommon.NParVec,
    "VectorCommon.NParVec/I") != NULL);

  succ = succ &&
    (tree->Branch("VectorCommon.IOrgVec", VectorCommon.IOrgVec,
    "VectorCommon.IOrgVec[VectorCommon.NParVec]/I") != NULL);
  succ = succ &&
    (tree->Branch("VectorCommon.IPVec", VectorCommon.IPVec,
    "VectorCommon.IPVec[VectorCommon.NParVec]/I") != NULL);
  succ = succ &&
    (tree->Branch("VectorCommon.IChVec", VectorCommon.IChVec,
    "VectorCommon.IChVec[VectorCommon.NParVec]/I") != NULL);
  succ = succ &&
    (tree->Branch("VectorCommon.IFlgVec", VectorCommon.IFlgVec,
    "VectorCommon.IFlgVec[VectorCommon.NParVec]/I") != NULL);

  succ = succ &&
    (tree->Branch("VectorCommon.AbsPVec", VectorCommon.AbsPVec,
    "VectorCommon.AbsPVec[VectorCommon.NParVec]/F") != NULL);
  succ = succ &&
    (tree->Branch("VectorCommon.PVec", VectorCommon.PVec,
   ( "VectorCommon.PVec["+tostr(MaxVC)+"][3]/F").c_str()) != NULL);
  succ = succ &&
    (tree->Branch("VectorCommon.PosV", VectorCommon.PosV,
    "VectorCommon.PosV[3]/F") != NULL);

  return succ;
}

inline Bool_t VectorCommon_SetBranchAddresses(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.NParVec", &VectorCommon.NParVec));

  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.IOrgVec", VectorCommon.IOrgVec));
  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.IPVec", VectorCommon.IPVec));
  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.IChVec", VectorCommon.IChVec));
  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.IFlgVec", VectorCommon.IFlgVec));

  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.AbsPVec", VectorCommon.AbsPVec));
  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.PVec", VectorCommon.PVec));
  succ = succ &&
    (!tree->SetBranchAddress("VectorCommon.PosV", VectorCommon.PosV));

  return succ;
}

inline void Copy_VectorCommon_F2C(){
  VectorCommon.NParVec = vcwork->nvc;
  for(Int_t n1 = 3; (n1 < vcwork->nvc) && (n1 < MaxVC); ++n1){
    VectorCommon.IOrgVec[n1] = vcwork->iorgvc[n1];
    VectorCommon.IPVec[n1] = vcwork->ipvc[n1];
    VectorCommon.IChVec[n1] = vcwork->icrnvc[n1];
    VectorCommon.IFlgVec[n1] = vcwork->iflgvc[n1];
    VectorCommon.AbsPVec[n1] = sqrt(vcwork->pvc[n1][0]*vcwork->pvc[n1][0]
                                  + vcwork->pvc[n1][1]*vcwork->pvc[n1][1]
                                  + vcwork->pvc[n1][2]*vcwork->pvc[n1][2]);
    VectorCommon.PVec[n1][0] = vcwork->pvc[n1][0];
    VectorCommon.PVec[n1][1] = vcwork->pvc[n1][1];
    VectorCommon.PVec[n1][2] = vcwork->pvc[n1][2];
  }
}

inline void Copy_VectorCommon_C2F(){
  vcwork->nvc = VectorCommon.NParVec;
  for(Int_t n1 = 3; (n1 < vcwork->nvc) && (n1 < MaxVC); ++n1){
    vcwork->iorgvc[n1] = VectorCommon.IOrgVec[n1];
    vcwork->ipvc[n1] = VectorCommon.IPVec[n1];
    vcwork->icrnvc[n1] = VectorCommon.IChVec[n1];
    vcwork->iflgvc[n1] = VectorCommon.IFlgVec[n1];
    vcwork->pvc[n1][0] = VectorCommon.PVec[n1][0];
    vcwork->pvc[n1][1] = VectorCommon.PVec[n1][1];
    vcwork->pvc[n1][2] = VectorCommon.PVec[n1][2];
  }
}

struct FSIHistoryCommon_{
  Int_t NVert;
  Int_t NVCVert;

  Int_t IFlgVert[MaxVert];
  Int_t IPVert[MaxVCVert];
  Int_t IVertI[MaxVCVert];
  Int_t IVertF[MaxVCVert];

  Float_t PosVert[MaxVert][3];
  Float_t DirVert[MaxVCVert][3];
  Float_t AbsPVert[MaxVCVert];
  Float_t AbsTPVert[MaxVCVert];
  Float_t FSIProb;
} FSIHistoryCommon;

inline Bool_t Add_FSIHistoryCommon_Branches(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.NVert", &FSIHistoryCommon.NVert,
    "FSIHistoryCommon.NVert/I") != NULL);
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.NVCVert", &FSIHistoryCommon.NVCVert,
    "FSIHistoryCommon.NVCVert/I") != NULL);

  succ = succ &&
    (tree->Branch("FSIHistoryCommon.IFlgVert", FSIHistoryCommon.IFlgVert,
    "FSIHistoryCommon.IFlgVert[FSIHistoryCommon.NVert]/I") != NULL);
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.IPVert", FSIHistoryCommon.IPVert,
    "FSIHistoryCommon.IPVert[FSIHistoryCommon.NVCVert]/I") != NULL);
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.IVertI", FSIHistoryCommon.IVertI,
    "FSIHistoryCommon.IVertI[FSIHistoryCommon.NVCVert]/I") != NULL);
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.IVertF", FSIHistoryCommon.IVertF,
    "FSIHistoryCommon.IVertF[FSIHistoryCommon.NVCVert]/I") != NULL);

  succ = succ &&
    (tree->Branch("FSIHistoryCommon.PosVert", FSIHistoryCommon.PosVert,
    ("FSIHistoryCommon.PosVert["+tostr(MaxVert)+"][3]/F").c_str() ) != NULL);
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.DirVert", FSIHistoryCommon.DirVert,
    ("FSIHistoryCommon.DirVert["+tostr(MaxVCVert)+"][3]/F").c_str() ) != NULL);
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.AbsPVert", FSIHistoryCommon.AbsPVert,
    "FSIHistoryCommon.AbsPVert[FSIHistoryCommon.NVCVert]/F") != NULL);
  succ = succ &&
    (tree->Branch("FSIHistoryCommon.AbsTPVert", FSIHistoryCommon.AbsTPVert,
    "FSIHistoryCommon.AbsTPVert[FSIHistoryCommon.NVCVert]/F") != NULL);

  succ = succ &&
    (tree->Branch("FSIHistoryCommon.FSIProb", &FSIHistoryCommon.FSIProb,
    "FSIHistoryCommon.FSIProb/F") != NULL);

  return succ;
}

inline Bool_t FSIHistoryCommon_SetBranchAddresses(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.NVert",
      &FSIHistoryCommon.NVert));
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.NVCVert",
      &FSIHistoryCommon.NVCVert));

  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.IFlgVert",
      FSIHistoryCommon.IFlgVert));
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.IPVert",
      FSIHistoryCommon.IPVert));
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.IVertI",
      FSIHistoryCommon.IVertI));
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.IVertF",
      FSIHistoryCommon.IVertF));

  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.PosVert",
      FSIHistoryCommon.PosVert));
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.DirVert",
      FSIHistoryCommon.DirVert));
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.AbsPVert",
      FSIHistoryCommon.AbsPVert));
  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.AbsTPVert",
      FSIHistoryCommon.AbsTPVert));

  succ = succ &&
    (!tree->SetBranchAddress("FSIHistoryCommon.FSIProb",
      &FSIHistoryCommon.FSIProb));

  return succ;
}

inline void Copy_FSIHistoryCommon_F2C(){
  FSIHistoryCommon.NVert = fsihist->nvert;
  for(Int_t IVert = 0; (IVert < fsihist->nvert) && (IVert < MaxVert);
    ++IVert){

    FSIHistoryCommon.IFlgVert[IVert] = fsihist->iflgvert[IVert];
    FSIHistoryCommon.PosVert[IVert][0] = fsihist->posvert[IVert][0];
    FSIHistoryCommon.PosVert[IVert][1] = fsihist->posvert[IVert][1];
    FSIHistoryCommon.PosVert[IVert][2] = fsihist->posvert[IVert][2];
  }

  FSIHistoryCommon.NVCVert = fsihist->nvcvert;
  for(Int_t IVCVert = 0;
    (IVCVert < fsihist->nvcvert) && (IVCVert < MaxVCVert); ++IVCVert){

    FSIHistoryCommon.IPVert[IVCVert] = fsihist->ipvert[IVCVert];
    FSIHistoryCommon.IVertI[IVCVert] = fsihist->iverti[IVCVert];
    FSIHistoryCommon.IVertF[IVCVert] = fsihist->ivertf[IVCVert];
    FSIHistoryCommon.DirVert[IVCVert][0] = fsihist->dirvert[IVCVert][0];
    FSIHistoryCommon.DirVert[IVCVert][1] = fsihist->dirvert[IVCVert][1];
    FSIHistoryCommon.DirVert[IVCVert][2] = fsihist->dirvert[IVCVert][2];
    FSIHistoryCommon.AbsPVert[IVCVert] = fsihist->abspvert[IVCVert];
    FSIHistoryCommon.AbsTPVert[IVCVert] = fsihist->abstpvert[IVCVert];
  }
  FSIHistoryCommon.FSIProb = fsihist->fsiprob;
}

inline void Copy_FSIHistoryCommon_C2F(){
  fsihist->nvert = FSIHistoryCommon.NVert;
  for(Int_t IVert = 0; (IVert < fsihist->nvert) && (IVert < MaxVert);
    ++IVert){

    fsihist->iflgvert[IVert] = FSIHistoryCommon.IFlgVert[IVert];
    fsihist->posvert[IVert][0] = FSIHistoryCommon.PosVert[IVert][0];
    fsihist->posvert[IVert][1] = FSIHistoryCommon.PosVert[IVert][1];
    fsihist->posvert[IVert][2] = FSIHistoryCommon.PosVert[IVert][2];
  }

  fsihist->nvcvert = FSIHistoryCommon.NVCVert;
  for(Int_t IVCVert = 0;
    (IVCVert < fsihist->nvcvert) && (IVCVert < MaxVCVert); ++IVCVert){

    fsihist->ipvert[IVCVert] = FSIHistoryCommon.IPVert[IVCVert];
    fsihist->iverti[IVCVert] = FSIHistoryCommon.IVertI[IVCVert];
    fsihist->ivertf[IVCVert] = FSIHistoryCommon.IVertF[IVCVert];
    fsihist->dirvert[IVCVert][0] = FSIHistoryCommon.DirVert[IVCVert][0];
    fsihist->dirvert[IVCVert][1] = FSIHistoryCommon.DirVert[IVCVert][1];
    fsihist->dirvert[IVCVert][2] = FSIHistoryCommon.DirVert[IVCVert][2];
    fsihist->abspvert[IVCVert] = FSIHistoryCommon.AbsPVert[IVCVert];
    fsihist->abstpvert[IVCVert] = FSIHistoryCommon.AbsTPVert[IVCVert];
  }
  fsihist->fsiprob = FSIHistoryCommon.FSIProb;
}

struct NucleusTargetCommon_{
  Int_t NumBndN;
  Int_t NumBndP;
  Int_t NumFreP;
  Int_t NumAtom;
} NucleusTargetCommon;

inline Bool_t Add_NucleusTargetCommon_Branches(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (tree->Branch("NucleusTargetCommon.NumBndN", &NucleusTargetCommon.NumBndN,
    "NucleusTargetCommon.NumBndN/I") != NULL);
  succ = succ &&
    (tree->Branch("NucleusTargetCommon.NumBndP", &NucleusTargetCommon.NumBndP,
    "NucleusTargetCommon.NumBndP/I") != NULL);
  succ = succ &&
    (tree->Branch("NucleusTargetCommon.NumFreP", &NucleusTargetCommon.NumFreP,
    "NucleusTargetCommon.NumFreP/I") != NULL);
  succ = succ &&
    (tree->Branch("NucleusTargetCommon.NumAtom", &NucleusTargetCommon.NumAtom,
    "NucleusTargetCommon.NumAtom/I") != NULL);

  return succ;
}

inline Bool_t NucleusTargetCommon_SetBranchAddresses(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (!tree->SetBranchAddress("NucleusTargetCommon.NumBndN",
      &NucleusTargetCommon.NumBndN));
  succ = succ &&
    (!tree->SetBranchAddress("NucleusTargetCommon.NumBndP",
      &NucleusTargetCommon.NumBndP));
  succ = succ &&
    (!tree->SetBranchAddress("NucleusTargetCommon.NumFreP",
      &NucleusTargetCommon.NumFreP));
  succ = succ &&
    (!tree->SetBranchAddress("NucleusTargetCommon.NumAtom",
      &NucleusTargetCommon.NumAtom));

  return succ;
}

inline void Copy_NucleusTargetCommon_F2C(){
  NucleusTargetCommon.NumBndN = neuttarget->numbndn;
  NucleusTargetCommon.NumBndP = neuttarget->numbndp;
  NucleusTargetCommon.NumFreP = neuttarget->numfrep;
  NucleusTargetCommon.NumAtom = neuttarget->numatom;
}

inline void Copy_NucleusTargetCommon_C2F(){
  neuttarget->numbndn = NucleusTargetCommon.NumBndN;
  neuttarget->numbndp = NucleusTargetCommon.NumBndP;
  neuttarget->numfrep = NucleusTargetCommon.NumFreP;
  neuttarget->numatom = NucleusTargetCommon.NumAtom;
}

struct NucleonFSIHistoryCommon_ {
  Int_t NFNVert;
  Int_t NFNStep;
  Int_t NFIFlag[MaxNucleonVert];
  Int_t NFFirstStep[MaxNucleonVert];

  Float_t NFx[MaxNucleonVert];
  Float_t NFy[MaxNucleonVert];
  Float_t NFz[MaxNucleonVert];
  Float_t NFpx[MaxNucleonVert];
  Float_t NFpy[MaxNucleonVert];
  Float_t NFpz[MaxNucleonVert];
  Float_t NFe[MaxNucleonVert];
  Float_t NFptot[MaxNucleonStep];
  Float_t NFecms2[MaxNucleonStep];
} NucleonFSIHistoryCommon;

inline Bool_t Add_NucleonFSIHistoryCommon_Branches(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFNVert",
      &NucleonFSIHistoryCommon.NFNVert,
    "NucleonFSIHistoryCommon.NFNVert/I") != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFNStep",
      &NucleonFSIHistoryCommon.NFNStep,
    "NucleonFSIHistoryCommon.NFNStep/I") != NULL);

  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFIFlag",
      NucleonFSIHistoryCommon.NFIFlag,
    "NucleonFSIHistoryCommon.NFIFlag[NucleonFSIHistoryCommon.NFNVert]/I")
    != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFFirstStep",
      NucleonFSIHistoryCommon.NFFirstStep,
    "NucleonFSIHistoryCommon.NFFirstStep[NucleonFSIHistoryCommon.NFNVert]/I")
    != NULL);

  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFx",
      NucleonFSIHistoryCommon.NFx,
    "NucleonFSIHistoryCommon.NFx[NucleonFSIHistoryCommon.NFNVert]/F") != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFy",
      NucleonFSIHistoryCommon.NFy,
    "NucleonFSIHistoryCommon.NFy[NucleonFSIHistoryCommon.NFNVert]/F") != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFz",
      NucleonFSIHistoryCommon.NFz,
    "NucleonFSIHistoryCommon.NFz[NucleonFSIHistoryCommon.NFNVert]/F") != NULL);


  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFpx",
      NucleonFSIHistoryCommon.NFpx,
    "NucleonFSIHistoryCommon.NFpx[NucleonFSIHistoryCommon.NFNVert]/F") != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFpy",
      NucleonFSIHistoryCommon.NFpy,
    "NucleonFSIHistoryCommon.NFpy[NucleonFSIHistoryCommon.NFNVert]/F") != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFpz",
      NucleonFSIHistoryCommon.NFpz,
    "NucleonFSIHistoryCommon.NFpz[NucleonFSIHistoryCommon.NFNVert]/F") != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFe",
      NucleonFSIHistoryCommon.NFe,
    "NucleonFSIHistoryCommon.NFe[NucleonFSIHistoryCommon.NFNVert]/F") != NULL);


  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFptot",
      NucleonFSIHistoryCommon.NFptot,
    "NucleonFSIHistoryCommon.NFptot[NucleonFSIHistoryCommon.NFNStep]/F") != NULL);
  succ = succ &&
    (tree->Branch("NucleonFSIHistoryCommon.NFecms2",
      NucleonFSIHistoryCommon.NFecms2,
    "NucleonFSIHistoryCommon.NFecms2[NucleonFSIHistoryCommon.NFNStep]/F") != NULL);

  return succ;
}

inline Bool_t NucleonFSIHistoryCommon_SetBranchAddresses(TTree * tree){
  Bool_t succ = true;
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFNVert",
      &NucleonFSIHistoryCommon.NFNVert) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFNStep",
      &NucleonFSIHistoryCommon.NFNStep) );

  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFIFlag",
      NucleonFSIHistoryCommon.NFIFlag) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFFirstStep",
      NucleonFSIHistoryCommon.NFFirstStep) );

  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFx",
      NucleonFSIHistoryCommon.NFx) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFy",
      NucleonFSIHistoryCommon.NFy) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFz",
      NucleonFSIHistoryCommon.NFz) );


  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFpx",
      NucleonFSIHistoryCommon.NFpx) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFpy",
      NucleonFSIHistoryCommon.NFpy) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFpz",
      NucleonFSIHistoryCommon.NFpz) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFe",
      NucleonFSIHistoryCommon.NFe) );


  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFptot",
      NucleonFSIHistoryCommon.NFptot) );
  succ = succ &&
    (!tree->SetBranchAddress("NucleonFSIHistoryCommon.NFecms2",
      NucleonFSIHistoryCommon.NFecms2) );

  return succ;
}

inline void Copy_NucleonFSIHistoryCommon_F2C(){
  NucleonFSIHistoryCommon.NFNVert = nucleonfsihist->nfnvert;
  for(Int_t IFNVert = 0;
    (IFNVert < nucleonfsihist->nfnvert) && (IFNVert < MaxNucleonVert);
    ++IFNVert){
    NucleonFSIHistoryCommon.NFIFlag[IFNVert] =
      nucleonfsihist->nfiflag[IFNVert];
    NucleonFSIHistoryCommon.NFFirstStep[IFNVert] =
      nucleonfsihist->nffirststep[IFNVert];
    NucleonFSIHistoryCommon.NFx[IFNVert] =
      nucleonfsihist->nfx[IFNVert];
    NucleonFSIHistoryCommon.NFy[IFNVert] =
      nucleonfsihist->nfy[IFNVert];
    NucleonFSIHistoryCommon.NFz[IFNVert] =
      nucleonfsihist->nfz[IFNVert];
    NucleonFSIHistoryCommon.NFpx[IFNVert] =
      nucleonfsihist->nfpx[IFNVert];
    NucleonFSIHistoryCommon.NFpy[IFNVert] =
      nucleonfsihist->nfpy[IFNVert];
    NucleonFSIHistoryCommon.NFpz[IFNVert] =
      nucleonfsihist->nfpz[IFNVert];
    NucleonFSIHistoryCommon.NFe[IFNVert] =
      nucleonfsihist->nfe[IFNVert];
  }

  NucleonFSIHistoryCommon.NFNStep = nucleonfsihist->nfnstep;
  for(Int_t IFNStep = 0;
    (IFNStep < nucleonfsihist->nfnstep) && (IFNStep < MaxNucleonStep);
    ++IFNStep){
    NucleonFSIHistoryCommon.NFptot[IFNStep] =
    nucleonfsihist->nfptot[IFNStep];
    NucleonFSIHistoryCommon.NFecms2[IFNStep] =
      nucleonfsihist->nfecms2[IFNStep];
  }
}

inline void Copy_NucleonFSIHistoryCommon_C2F(){
  nucleonfsihist->nfnvert = NucleonFSIHistoryCommon.NFNVert;
  for(Int_t IFNVert = 0;
    (IFNVert < nucleonfsihist->nfnvert) && (IFNVert < MaxNucleonVert);
    ++IFNVert){
    nucleonfsihist->nfiflag[IFNVert] =
      NucleonFSIHistoryCommon.NFIFlag[IFNVert];
    nucleonfsihist->nffirststep[IFNVert] =
      NucleonFSIHistoryCommon.NFFirstStep[IFNVert];
    nucleonfsihist->nfx[IFNVert] =
      NucleonFSIHistoryCommon.NFx[IFNVert];
    nucleonfsihist->nfy[IFNVert] =
      NucleonFSIHistoryCommon.NFy[IFNVert];
    nucleonfsihist->nfz[IFNVert] =
      NucleonFSIHistoryCommon.NFz[IFNVert];
    nucleonfsihist->nfpx[IFNVert] =
      NucleonFSIHistoryCommon.NFpx[IFNVert];
    nucleonfsihist->nfpy[IFNVert] =
      NucleonFSIHistoryCommon.NFpy[IFNVert];
    nucleonfsihist->nfpz[IFNVert] =
      NucleonFSIHistoryCommon.NFpz[IFNVert];
    nucleonfsihist->nfe[IFNVert] =
      NucleonFSIHistoryCommon.NFe[IFNVert];
  }

  nucleonfsihist->nfnstep = NucleonFSIHistoryCommon.NFNStep;
  for(Int_t IFNStep = 0;
    (IFNStep < nucleonfsihist->nfnstep) && (IFNStep < MaxNucleonStep);
    ++IFNStep){
    nucleonfsihist->nfptot[IFNStep] =
      NucleonFSIHistoryCommon.NFptot[IFNStep];
    nucleonfsihist->nfecms2[IFNStep] =
      NucleonFSIHistoryCommon.NFecms2[IFNStep];
  }
}
