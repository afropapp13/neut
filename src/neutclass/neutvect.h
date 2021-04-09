#ifndef _NEUTVECT_H_
#define _NEUTVECT_H_

#include <TArrayI.h>

#include <TObject.h>
#include <TLorentzVector.h>

#include "neutpart.h"
#include "neutfsipart.h"
#include "neutfsivert.h"
#include "neutnucfsivert.h"
#include "neutnucfsistep.h"

class NeutVect : public TObject {

 public:
  NeutVect(Int_t np = 0, Int_t nfsip = 0, Int_t nfsiv = 0);

  ~NeutVect();

  void ClearVars();

  void Dump();

  /* Generated vector information */

  Int_t          EventNo;  // Event number

  /* Target parameters */
  Int_t          TargetA;   // Atomic number of the target nucleus
  Int_t          TargetZ;   // Mass number of the target nucleus
  Int_t          TargetH;   // Number of free protons
  Int_t          Ibound;    // Interaction on nucleus (1) or free proton (0)

  /* Nuclear potentials */
  Double_t       VNuclIni;  // Nuclear potential for the initial state (MeV)
                            //
  Double_t       VNuclFin;  // Nuclear potential for the final state (MeV)


  /* Fermi momentum related */
  Double_t       PFSurf;    // Fermi surface momentum (MeV/c)
  Double_t       PFMax;     // Maximum value of the Fermi momentum (MeV/c)
  //

  /* Flux ID including detector code?*/
  Int_t          FluxID;    // Flux ID code
  //  0: Randomly generated

  // Interaction model
  // QE model
  Int_t          QEModel;
  // Axialform factor choice CCQE
  Int_t          QEAVForm;
  Int_t          QEVForm;
  // Radiative correction
  Int_t          RADcorr;
  // Single pion model
  Int_t          SPIModel;
  // Rein-Sehgal 1pi ejection
  Int_t          SPIRSEj;
  // Coherent model
  Int_t          COHModel;
  // DIS model
  Int_t          DISModel;

  // 1 pi & the other meson productions form factor type
  Int_t          SPIForm;
  Int_t          RESForm;

  // 1 pi & the other meson productions NR type ( form factor detail )
  Int_t          SPINRType;
  Int_t          RESNRType;

  // CCQE model parameters
  Float_t        QEMA;
  Float_t        QEMV;
  Float_t        KAPPA;

  // Single pion model parameters
  Float_t        SPIMA;
  Float_t        RESMA;
  Float_t        SPIMV;
  Float_t        RESMV;
  Float_t        SPICA5I;
  Float_t        SPIBGScale;
  Float_t        SPIMaBKGM;

  // Coherent model parameters
  Float_t        COHMA;
  Float_t        COHR0;
  Float_t        COHA1err;
  Float_t        COHb1err;
  // NEUT version
  Int_t          COREVer;
  Int_t          NUCEVer;
  Int_t          NUCCVer;

  // NEUT card configuration
  Int_t          FrmFlg;
  Int_t          PauFlg;
  Int_t          NefO16;
  Int_t          ModFlg;
  Int_t          SelMod;

  // Formation length
  Int_t          FormLen;
  // Pionless delta decay
  Int_t          IPilessDcy;
  Float_t        RPilessDcy;
  // Rescattering
  Int_t          NucScat;
  // Neutrino generating mode
  Float_t        NucFac;

  Float_t        NuceffKinVersion;

  Float_t        NuceffFactorPIQE;
  Float_t        NuceffFactorPIInel;
  Float_t        NuceffFactorPIAbs;
  Float_t        NuceffFactorPIQEH;
  Float_t        NuceffFactorPICX;
  Float_t        NuceffFactorPICXH;
  Float_t        NuceffFactorPICoh;
  Float_t        NuceffFactorPIAll;

  Float_t        NuceffFactorPIQEHKin;
  Float_t        NuceffFactorPIQELKin;
  Float_t        NuceffFactorPICXKin;

  /////////////////////////////////////////////////////////////////////////
  /* Nieves 1p1h model parameters */

  Int_t          NVQERFG;
  Int_t          NVQEBind;
  Int_t          NVQERPA;
  Float_t        XNVRPAFP0in;
  Float_t        XNVRPAPF0ex;
  Float_t        XNVRPAFstar;
  Float_t        XNVRPAF;
  Float_t        XNVRPAPILambda;
  Float_t        XNVRPACR0;
  Float_t        XNVRPARHOLambda;
  Float_t        XNVRPAGp;
  Float_t        XNVRPAXMPI;
  Float_t        XNVRPAXMRHO;
  Float_t        XNVRPAIrel;
  Int_t          FFTYPE;
  Int_t          NVBINDFermiCor;

  /////////////////////////////////////////////////////////////////////////
  /* Ineraction mode */
  Int_t          Mode;

  /* Cross sections */
  Float_t        Totcrs;

  Float_t        CrsEnergy;
  Float_t        DifCrsNE[8];
  Float_t        Crsx;
  Float_t        Crsy;
  Float_t        Crsz;
  Float_t        Crsphi;
  Float_t        Crsq2;

  /* Number of particles */
  Int_t          Npart(void)        const {return fNpart; }
  void           SetNpart(Int_t npart);

  /* Number of primary particles */
  Int_t          Nprimary(void)     const {return fNprimary; }
  void           SetNprimary(Int_t nprimary);

  /* Vertex ID */
  Int_t          VertexID(int idx)   const { return fVertexID[idx]; };
  void           SetVertexID(Int_t idx, Int_t VtxID);
  void           SetVertexID(Int_t np, Int_t *VtxID_array);

  /* Particle information (PID, 4 momentum etc.) */
  NeutPart       *PartInfo(Int_t idx)
  { return (NeutPart *)((fNpart>idx) ? (fPartInfo->At(idx)) : NULL );};
  void SetPartInfo(Int_t idx      , NeutPart PInfo);
  void SetPartInfo(Int_t npart   , NeutPart *PInfo_array);

  /* Parent particle's index */
  Int_t          ParentIdx(int idx)  const { return fParentIdx[idx]; };
  void           SetParentIdx(Int_t idx, Int_t ParentIdx);
  void           SetParentIdx(Int_t np, Int_t *ParentIdx_array);


  /* Number of FSI particles */
  Int_t          NfsiPart(void)        const {return fNfsiPart; }
  void           SetNfsiPart(Int_t npart);

  /* FSI Particle Information */
  NeutFsiPart       *FsiPartInfo(Int_t idx)
  { return (NeutFsiPart *)((fNfsiPart>idx) ? (fFsiPartInfo->At(idx)) : NULL );};
  void SetFsiPartInfo(Int_t idx      , NeutFsiPart PInfo);
  void SetFsiPartInfo(Int_t npart   , NeutFsiPart *PInfo_array);


  /* Number of FSI Vertices */
  Int_t          NfsiVert(void)        const {return fNfsiVert; }
  void           SetNfsiVert(Int_t nvert);

  /* FSI Vertex Information */
  NeutFsiVert       *FsiVertInfo(Int_t idx)
  { return (NeutFsiVert *)((fNfsiVert>idx) ? (fFsiVertInfo->At(idx)) : NULL );};
  void SetFsiVertInfo(Int_t idx     , NeutFsiVert VInfo);
  void SetFsiVertInfo(Int_t nvert   , NeutFsiVert *VInfo_array);

  Float_t Fsiprob;

  /* Number of Nucleon FSI Vertices */
  Int_t          NnucFsiVert(void)        const {return fNnucFsiVert; }
  void           SetNnucFsiVert(Int_t nnucvert);

  /* Nucleon FSI Vertex Information */
  NeutNucFsiVert       *NucFsiVertInfo(Int_t idx)
  { return (NeutNucFsiVert *)((fNnucFsiVert>idx) ? (fNucFsiVertInfo->At(idx)) : NULL );};
  void SetNucFsiVertInfo(Int_t idx      , NeutNucFsiVert VInfo);
  void SetNucFsiVertInfo(Int_t nnucvert , NeutNucFsiVert *VInfo_array);

  /* Number of Nucleon FSI Steps */
  Int_t          NnucFsiStep(void)        const {return fNnucFsiStep; }
  void           SetNnucFsiStep(Int_t nnucstep);

  /* Nucleon FSI Step Information */
  NeutNucFsiStep       *NucFsiStepInfo(Int_t idx)
  { return (NeutNucFsiStep *)((fNnucFsiStep>idx) ? (fNucFsiStepInfo->At(idx)) : NULL );};
  void SetNucFsiStepInfo(Int_t idx      , NeutNucFsiStep VInfo);
  void SetNucFsiStepInfo(Int_t nnucstep , NeutNucFsiStep *VInfo_array);

 private:

  Int_t          fNpart;     // Number of particles

  Int_t          fNprimary;  // Number of primary particles

  TObjArray      *fPartInfo; // ->

  TArrayI        fVertexID;  //         Vertex Index ID
  //         ( To refer NeutVertex )

  TArrayI        fParentIdx;  //        Idx of parent particle

  Int_t          fRandomSeed[16]; // Random seeds

  Int_t          fNfsiPart;  // Number of FSI particles
  Int_t          fNfsiVert;  // Number of FSI vertices

  TObjArray      *fFsiPartInfo; // ->
  TObjArray      *fFsiVertInfo; // ->

  Int_t          fNnucFsiVert;  // Number of Nucleon FSI vertices
  Int_t          fNnucFsiStep;  // Number of Nucleon FSI steps

  TObjArray      *fNucFsiVertInfo; // ->
  TObjArray      *fNucFsiStepInfo; // ->

  TLorentzVector fZeroVect;  // !

  ClassDef(NeutVect, 6)
};

#endif
