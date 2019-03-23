#include "neutread.h"

// STL headers
#include <iostream>
#include <iomanip>

// ROOT headers
#include <TTree.h>

// Headers containing TObjects from `${NEUT_ROOT}/src/neutclass`
#include "neutvect.h"
#include "neutvtx.h"

//#define NEUT_REWEIGHT_DEBUG

void readtree(TTree* tree, long long entry)
{
#ifdef NEUT_REWEIGHT_DEBUG
  std::cout << "[" << __FILE__ << ":L" << __LINE__ << "] "
            << "Don't use this function in combination with a NeutrootSingleton object as it will "
            << "invalidate the pointers in the other object. This function should only be used in "
            << "small scale debugging situations.\n"
            << "Instead, get the NeutVect pointer from the singleton and directly call readvect()\n";
#endif

  if (entry < 0 || tree->GetEntries() <= entry) {
#ifdef NEUT_REWEIGHT_DEBUG
    std::cerr << "[" << __FILE__ << ":L" << __LINE__ << "] "
              << " -!- Skipping illegal entry: " << entry
              << ", which is not in the allowed range [0, " << tree->GetEntries() << ")\n";
#endif
    return;
  }

  //
  // From TBranchElement::SetAddress(void* addobj)
  //
  //  ...
  //
  //    If addr is not zero and the pointer addr points at is
  //    also not zero, then the caller has allocated a branch
  //    object and is asking us to use it.  The caller owns it
  //    and must delete it when it is no longer needed.
  //
  //    Example:
  //
  //      Event* event = new Event();
  //      branch->SetAddress(&event);
  //      ... Do some work.
  //      delete event;
  //      event = 0;
  //
  //    These rules affect users of TTree::Branch(),
  //    TTree::SetBranchAddress(), and TChain::SetBranchAddress()
  //    as well because those routines call this one.
  //
  //  ...
  // 
  NeutVect* nvect = new NeutVect();
  NeutVtx* nvtx = new NeutVtx();

  tree->ResetBranchAddresses();
  tree->SetBranchStatus("vectorbranch", 1);
  tree->SetBranchStatus("vertexbranch", 1);
  tree->SetBranchAddress("vectorbranch", &nvect);
  tree->SetBranchAddress("vertexbranch", &nvtx);

  tree->GetEntry(entry);
  readvect(nvect, nvtx);
  readvtx(nvect, nvtx);

  delete nvect;
  delete nvtx;

  // Intentionally omit resetting the pointer to NULL,
  // because they go out of scope immediately and will
  // never be used again.
  //
  //nvect = nullptr;
  //nvtx = nullptr;
}


// Read the NeutVect object
void readvect(NeutVect* nvect, NeutVtx* nvtx)
{
  readnework(nvect, nvtx);
  readvcwork(nvect, nvtx);
  readneutparams(nvect, nvtx);
  readposinnuc(nvect, nvtx);
  readnecard(nvect, nvtx);
  readnrcard(nvect, nvtx);
  readneutmodel(nvect, nvtx);
  readver(nvect, nvtx);
  readneutcrs(nvect, nvtx);
  readfsihist(nvect, nvtx);
  readnucleonfsihist(nvect, nvtx);
}


// Read the NeutVtx object
void readvtx(NeutVect* nvect, NeutVtx* nvtx) {
  readvcvrtx(nvect, nvtx);
}


/**
 *     COMMON /NEWORK/
 *     =================================================================
 *     MODENE      : MODE OF INTERACTION
 *     NUMNE       : # OF PARTICLE
 *     IPNE(I)     : PARTICLE CODE OF I-TH PARTICLE
 *     PNE(3,I)    : MOMENTUM OF I-TH PARTICLE ( GEV/C )
 *     IORGNE(I)   : ID OF ORIGIN PARTICLE
 *     IFLGNE(I)   : FLAG OF FINAL STATE
 *                  -1 : initial particle
 *                   0 : DETERMINED LATER PROCEDURE
 *                   1 : DECAY TO OTHER PARTICLE
 *                   2 : ESCAPE FROM DETECTOR
 *                   3 : ABSORPTION
 *                   4 : CHARGE EXCHANGE
 *                   5 : STOP AND NOT CONSIDER IN M.C.
 *                   6 : E.M. SHOWER
 *     ICRNNE(I)  : FLAG OF TO CHASE OR NOT
 *                   0 : DO NOT CHASE
 *                   1 : CHASE
 *     =================================================================
 **/
void readnework(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  nework_.modene = nvect->Mode;
  nework_.numne = nvect->Nprimary();

  const float MeV_per_GeV = 1.e3;
  NeutPart* pinfo = nullptr;
  int npart = nvect->Npart();
  for (int i=0; i<npart; ++i) {
    pinfo = nvect->PartInfo(i);
    nework_.ipne[i] = pinfo->fPID;
    nework_.pne[i][0] = pinfo->fP.X() / MeV_per_GeV;
    nework_.pne[i][1] = pinfo->fP.Y() / MeV_per_GeV;
    nework_.pne[i][2] = pinfo->fP.Z() / MeV_per_GeV;
    nework_.iorgne[i] = nvect->ParentIdx(i);
    nework_.iflgne[i] = pinfo->fStatus;
    nework_.icrnne[i] = pinfo->fIsAlive;
  }
}


/**
 *     COMMON /VCWORK/
 *     =================================================================
 *     NVC         : # OF PARTICLE
 *     POSVC(3)    : VERTEX POSITION OF INITIAL INTERACTION
 *     IPVC(I)     : PARTICLE CODE OF I-TH PARTICLE
 *     AMASVC(I)   : MASS OF I-TH PARTICLE ( MEV/C**2 )
 *     PVC(3,I)    : MOMENTUM OF I-TH PARTICLE ( MEV/C )
 *     IORGVC(I)   : ID OF ORIGIN PARTICLE
 *     IFLGVC(I)   : FLAG OF FINAL STATE
 *                   0 : DETERMINED LATER PROCEDURE
 *                   1 : DECAY TO OTHER PARTICLE
 *                   2 : ESCAPE FROM DETECTOR
 *                   3 : ABSORPTION
 *                   4 : CHARGE EXCHANGE
 *                   5 : STOP AND NOT CONSIDER IN M.C.
 *                   6 : E.M. SHOWER
 *                   7 : HADRON PRODUCTION
 *                   8 : QUASI-ELASTIC SCATTER
 *                   9 : FORWARD (ELASTIC-LIKE) SCATTER
 *     ICRNVC(I)   : FLAG OF TO CHASE OR NOT
 *                   0 : DO NOT CHASE
 *                   1 : CHASE
 *     TIMVC(I)    : TIME OF VERTEX
 *     POSIVC(3,I) : INITIAL VERTEX POSITION
 *     IVTIVC(I)   : ID OF INITIAL VERTEX
 *     POSFVC(3,I) : FINAL VERTEX POSITION
 *     IVTFVC(I)   : ID OF FINAL VERTEX
 *     =================================================================
 **/
void readvcwork(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  NeutPart* pinfo = nullptr;
  vcwork_.nvc = nvect->Npart();
  for (int i=0; i<vcwork_.nvc; ++i) {
    vcwork_.ivtivc[i] = nvect->VertexID(i);
    vcwork_.iorgvc[i] = nvect->ParentIdx(i);

    pinfo = nvect->PartInfo(i);
    vcwork_.ipvc[i] = pinfo->fPID;
    vcwork_.amasvc[i] = pinfo->fMass;
    vcwork_.icrnvc[i] = pinfo->fIsAlive;
    vcwork_.iflgvc[i] = pinfo->fStatus;

    vcwork_.pvc[i][0] = pinfo->fP.X();
    vcwork_.pvc[i][1] = pinfo->fP.Y();
    vcwork_.pvc[i][2] = pinfo->fP.Z();

    vcwork_.posivc[i][0] = pinfo->fPosIni.X();
    vcwork_.posivc[i][1] = pinfo->fPosIni.Y();
    vcwork_.posivc[i][2] = pinfo->fPosIni.Z();

    vcwork_.posfvc[i][0] = pinfo->fPosFin.X();
    vcwork_.posfvc[i][1] = pinfo->fPosFin.Y();
    vcwork_.posfvc[i][2] = pinfo->fPosFin.Z();
  }
}


/**
 *     COMMON /NENUPR/
 *     =================================================================
 *     PFSURF : Fermi srface momentum in (GeV)
 *
 *     PFMAX  : Maximum value of the Fermi momentum (usually, this value is same as PFSURF.) in (GeV)
 *
 *     VNUINI : Nuclear potential (initial state) Negative VALUE (GeV)
 *
 *     VNUFIN : Nuclear potential (final state) Negative VALUE (GeV)
 *
 *     IFORMLEN: FORMATION LENGTH effect ON/OFF
 *               IFORMLEN=   1  : ALL ON (default)
 *               IFORMLEN=   0  : ALL OFF
 *               IFORMLEN= 110  : OFF for QE/Elastic
 *               IFORMLEN= 100  : ON  for mPi/DIS only
 *
 *     FZMU2   : FORMATION ZONE FREE PARAMETER (SCAT MODEL, DEFAULT 0.08)
 *     =================================================================
 *
 *     COMMON /NEFFPR/
 *     =================================================================
 *     FEFQE  : Correction Factor to the probability (p<500)
 *                               of quasielastic scattering ( single pi )
 *
 *     FEFQEH : Correction Factor to the probability (p>500)
 *                               of quasielastic scattering ( single pi )
 *  
 *     FEFINEL: Correction Factor to the probability (p>500)
 *                               of hadron prod. (inel. scatter)
 *  
 *     FEFABS : Correction Factor to the probability (p<500)
 *                               of absorption
 *  
 *     FEFCOH : Correction Factor to the probability (p>500)
 *                               of coherent(forward)-scattering
 *  
 *     FEFCX  : Correction Factor to the charge exchange amplitude (p<500)
 *  
 *     FEFCXH : Correction Factor to the charge exchange prob. (p>500)
 *  
 *     FEFQEHF: Portion of QE scattering that has isotropic resonance decay
 *                               kinematics (p>500)
 *  
 *     FEFCOHF: Amount of forward (elastic) scatter relative
 *                               to inelastic (p<500)
 *  
 *     FEFCXHF: Portion of inel. scattering that includes true CX (p>500)
 *  
 *     FEFCOUL: Coulomb correction to pion cascade trajectory (0: default, off)
 *  
 *     FEFALL:  Correction factor to the total mean free path (all momenta)
 *  
 **/
void readneutparams(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
    nenupr_.pfsurf = nvect->PFSurf;
    nenupr_.pfmax = nvect->PFMax;
    nenupr_.vnuini = nvect->VNuclIni;
    nenupr_.vnufin = nvect->VNuclFin;
    nenupr_.iformlen = nvect->FormLen;
    nenupr_.fzmu2 = 0.08; // (not filled)

    neffpr_.fefqe = nvect->NuceffFactorPIQE;
    neffpr_.fefqeh = nvect->NuceffFactorPIQEH;
    neffpr_.fefinel = nvect->NuceffFactorPIInel;
    neffpr_.fefabs = nvect->NuceffFactorPIAbs;
    neffpr_.fefcoh = nvect->NuceffFactorPICoh;
    neffpr_.fefqehf = nvect->NuceffFactorPIQEHKin;
    neffpr_.fefcohf = nvect->NuceffFactorPIQELKin;
    neffpr_.fefcx = nvect->NuceffFactorPICX;
    neffpr_.fefcxhf = nvect->NuceffFactorPICXKin;
    neffpr_.fefcxh = nvect->NuceffFactorPICXH;
    neffpr_.fefcoul = 0; // (not filled)

# if (NECOREVER == 5321) || (NECOREVER == 5322)
    // special case for 5.3.2.1 and 5.3.2.2
    // do nothing.
# elif (NECOREVER >= 535)
  if (nvect->Class_Version() >= 4) {
    neffpr_.fefall = nvect->NuceffFactorPIAll;
  }
# endif
}


/**
 *     COMMON /POSINNUC/
 *     =================================================================
 *     IBOUND : INTERACTION ON BOUND NUCLEON OR FREE NUCLEON
 *
 *     POSNUC(3,I) : INITIAL VERTEX POSITION
 **/
void readposinnuc(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  posinnuc_.ibound = nvect->Ibound;

  NeutPart* pinfo = nullptr;
  int npart = nvect->Npart();
  for (int i=0; i<npart; ++i) {
    pinfo = nvect->PartInfo(i);
    posinnuc_.posnuc[i][0] = pinfo->fPosIni.X();
    posinnuc_.posnuc[i][1] = pinfo->fPosIni.Y();
    posinnuc_.posnuc[i][2] = pinfo->fPosIni.Z();
  }
}


/**
 *     COMMON /NEUTMODEL/
 *     =================================================================
 *     MODELDIS = 70 ; GRV94 Original
 *                71 ; GRV94 Bodek
 *               120 ; GRV98 Original
 *               121 ; GRV98 Bodek
 *     
 *      MDLCOH =   0 ; Rein&Sehgal w/ lepton mass corr.
 *                 1 ; Kartavtsev et al.
 *                 2 ; Berger&Sehgal
 *     
 *     =================================================================
 *
 *     COMMON /NEMDLS/
 *     =================================================================
 *     MDLQE   = XXX1 ; CC : Smith-Moniz cross-section with classical kinematics
 *               XXX2 ; CC : Smith-Moniz cross-section with classical kinematics
 *                          ( BBBA05 vector form factor )
 *
 *             = XX0X ; NC : Original code
 *             = XX1X ; NC : Dipole
 *             = XX2X ; NC : BBBA05 calculation
 *
 *             = X0XX ; ALL; No correction on GM form factor
 *             = X1XX ; ALL; Correction on GM form factor ( Bodek et al )
 *
 *             = X4XX ; ALL; Use SF nuclear model for kinematics
 *
 *             = 1XXX ; ALL; Apply non-relativisitic RPA corrections
 *
 *     -- Axial form factor for QE
 *     MDLQEAF = 1 ; Simple dipole form factor
 *             = 0 ; BBBA07
 *
 *     XMAQE (XMVQE) : Axial (Vector) mass value of QE
 *
 *     KAPP : Kappa factor for QE (ala MiniBooNE)
 *
 *     PFSF : FermiMomentum used in SF model only (for pauli blocking)
 *              (if user doesn't define, uses default for nucleus)
 *
 *     SCCFV (SCCFA) : Second class form factors for CCQE (both 0.0 as default)
 *
 *     QEFP : Error term for pseudoscalar form factor for CCQE
 *
 *     -- Model for Single pion production
 *     MDLSPI = 1 ; Rein Sehgal
 *
 *     XMARSRES (XMVRSRES) : MA (MV) for Rein-Sehgal Single meson production
 *
 *     XMANFFRES (XMVNFFRES) : MA for New Form factor Single meson production
 *
 *     XMACOH  (Default = 1.0)
 *     RAD0NU  (Default = 1 fm)
 *     fA1COH  (Default = 0.0)
 *     fb1COH  (Default = 0.0)
 *     =================================================================
 **/
void readneutmodel(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  nemdls_.mdlqe = nvect->QEModel;
  nemdls_.mdlspi = nvect->SPIModel;
  nemdls_.mdldis = nvect->DISModel;
  nemdls_.mdlcoh = nvect->COHModel;
  nemdls_.mdlqeaf = nvect->QEVForm;
  nemdls_.xmaqe = nvect->QEMA;
  nemdls_.xmvqe = nvect->QEMV;
  nemdls_.xmaspi = nvect->SPIMA;
  nemdls_.xmvspi = nvect->SPIMV;
  nemdls_.kapp = nvect->KAPPA;
  nemdls_.xmacoh = nvect->COHMA;
  nemdls_.rad0nu = nvect->COHR0;
  nemdls_.iffspi = 1; // (not filled, value from neutcore/necard.F:NEIFF)
  nemdls_.nrtypespi = 0; // (not filled, value from neutcore/necard.F:NENRTYPE)
  nemdls_.rca5ispi = 1.01; // (not filled, value from neutcore/necard.F:RNECA5I)
  nemdls_.rbgsclspi = 1.30; // (not filled, value from neutcore/necard.F:RNEBGSCL)
  nemdls_.xmares = nvect->RESMA;
  nemdls_.xmvres = nvect->RESMV;
  nemdls_.sccfv = 0.; // (not filled, value from neutcore/necard.F)
  nemdls_.sccfa = 0.; // (not filled, value from neutcore/necard.F)
  nemdls_.pfsf = -1.; // (not filled, value from neutcore/necard.F)

  neutmodel_.modeldis = neutdis_.nepdf*10 + neutdis_.nebodek; // (not filled)
  neutmodel_.modelcoh = nemdls_.mdlcoh;

# if (NECOREVER == 5321) || (NECOREVER == 5322)
    // special case for 5.3.2.1 and 5.3.2.2
    // do nothing.
# elif (NECOREVER >= 533)
  if (nvect->Class_Version() >= 3) {
    nemdls_.fa1coh = nvect->COHA1err;
    nemdls_.fb1coh = nvect->COHb1err;
    nemdls_.fpqe = 1.; // (not filled)
  }
# endif
}


/*       NEFRMFLG : ( NEUT-FERM ) Flag on Fermi motion*/
/*                     0 : on ( default ) */
/*                     1 : off*/
/*       NEPAUFLG : ( NEUT-PAUL ) Flag on Pauli brokking*/
/*                     0 : on ( default ) */
/*                     1 : off*/
/*       NENEFO16 : ( NEUT-NEFF ) Flag on Nuclear effect in O16*/
/*                     0 : on ( default ) */
/*                     1 : off*/
/*       NENEFMODL: ( NEUT-MODL ) Select model for low energy pion*/
/*                     0 : default (Salcedo et. al)*/
/*                     1 : Tuned to pion scattering data*/
/*       NENEFMODH: ( NEUT-MODH ) Select model for high energy pion*/
/*                     0 : default*/
/*                     1 : proton/neutron separated*/
/*       NENEFKINH: ( NEUT-KINH ) Select kinematical model for */
/*                                   hi-nrg pion inelastic scattering*/
/*                     0 : Isotropic resonance decay*/
/*                     1 : SAID WI08 Phase Shifts*/
/*       neabspiemit : ( NEUT-neabspiemit ) Flag on Nucleon Ejection*/
/*                     0 : off */
/*                     1 : on ( default ) */
/*       NUSIM       : Neutrino Simulation Flag*/
/*                     1 : Neutrino Simulation   (default)*/
/*                     0 : Other (piscat, gampi)*/
/**/
/*       NEMODFLG : ( NEUT-MODE ) Flag on Interaction mode on neutrino int.*/
/*                     0 : normal ( default )*/
/*                    -1 : input cross section by CRSNEUT*/
/*                     n : sellect one mode ( n > 0 )   See nemodsel.F*/
/*                           n =  1 : charged current Q.E. */
/*                           n = 11,12,13*/
/*                                  : charged current Single pi production */
/*                           n = 21 : charged current Multi pi production*/
/*                           n = 31,32,33,34*/
/*                                  : neutral current Single pi production */
/*                           n = 41 : neutral current Multi pi production*/
/*                           n = 51,52 : neutral current elastic*/
/*                           n = 22,42,43 : single eta production */
/*                           n = 23,44,45 : single  K  production */
/* */
/*       CRSNEUT(28)   : ( NEUT-CRS ) Multiplied factor to cross section*/
/*                                    on each mode.   See nemodsel.F*/
/*       CRSNEUTB(28)  : ( NEUT-CRSB ) Multiplied factor to cross section*/
/*                                    on each mode.   See nemodsel.F*/
/**/
/*       NECOHEPI      : ( NEUT_COHEPI ) Select Coherent pi model*/
/*                            0 : Rein & Sehgal*/
/*                            1 : Kartavtsev */
/*                            2 : Berger & Sehgal*/
/**/
/*       NEPDF         : ( NEUT-Select Parton distribution function*/
/*                           n = 7 : GRV94*/
/*                           n =12 : GRV98*/
/**/
/*       NEBODEK       : ( NEUT- urn off/on Bodek-Yang correction*/
/*                       0 :   off*/
/*                       1 :   on */
/**/
/*       ITAUFLGCORE  : control Tau run mode(this is not set here)*/
/**/
/*       NUMBNDN  : total number of neutron      (e.g. H2O => 8 , Fe => 30)*/
/*       NUMBNDP  : total number of bound proton (e.g. H2O => 8 , Fe => 26)*/
/*       NUMFREP  : total number of free proton  (e.g. H2O => 2 , Fe =>  0)*/
/*       NUMATOM  : atomic number                (e.g.   O =>16 , Fe => 56)*/
/**/
/*       IPILESSDCY: Pion less delta deacy             ( 1 : on / 0: off )*/
/*       RPILESSDCY: Fraction of pion less delta deacy ( default 0.2 )*/
/**/
/*       NEIFF:    Tells which form factors are used (see rsdcrs)*/
/*       NENRTYPE: For R-S form factors, there is two separate*/
/**/
/*       IRADCORR:  Radiative correction on/off ( 1 : on / 0: off )*/
/*                 ( Currently off by default )*/
/**/
/*       QUIET : Screen output verbosity*/
/*         0 : Default (prints all initial state info)*/
/*         1 : Print only neutrino energy*/
/*         2 : Prints almost nothing (except PYTHIA output)*/
/**/
void readnecard(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  neutcard_.nefrmflg = nvect->FrmFlg;
  neutcard_.nepauflg = nvect->PauFlg;
  neutcard_.nenefo16 = nvect->NefO16;
  neutcard_.nenefmodl = 1; // (not filled, value from neutcore/necard.F)
  neutcard_.nenefmodh = 1; // (not filled, value from neutcore/necard.F)
  neutcard_.nenefkinh = 1; // (not filled, value from neutcore/necard.F)
  neutcard_.nemodflg = nvect->ModFlg;
  neutcard_.neselmod = nvect->SelMod;
  {
  const float mode_weight = (neutcard_.nemodflg == -1) ? 1. : 0.;
  for (size_t i=0; i!=28; ++i) {
    neutcard_.crsneut[i]  = mode_weight; // (not filled)
    neutcard_.crsneutb[i] = mode_weight; // (not filled) 
  }
  }
  neutcard_.itauflgcore = 0; // (not filled, value from neutcore/necard.F)
  neutcard_.nusim = 1; // (not filled, value from neutcore/necard.F)
  neutcard_.quiet = 0; // (not filled, value from neutcore/necard.F)

  neutdis_.nepdf = 12; // (not filled, value from neutcore/necard.F)
  neutdis_.nebodek = 1; // (not filled, value from neutcore/necard.F)

  neut1pi_.xmanffres = 0.95; // (not filled, value from neutcore/necard.F)
  neut1pi_.xmvnffres = 0.84; // (not filled, value from neutcore/necard.F)
  neut1pi_.xmarsres  = 1.21; // (not filled, value from neutcore/necard.F)
  neut1pi_.xmvrsres  = 0.84; // (not filled, value from neutcore/necard.F)
  neut1pi_.neiff = nvect->SPIForm;
  neut1pi_.nenrtype = nvect->SPINRType;
  neut1pi_.rneca5i = nvect->SPICA5I;
  neut1pi_.rnebgscl = nvect->SPIBGScale;

  neutcoh_.necohepi = 0; // (not filled, value from neutcore/necard.F)

  neuttarget_.numbndn = nvect->TargetA - nvect->TargetZ;
  neuttarget_.numbndp = nvect->TargetZ;
  neuttarget_.numfrep = nvect->TargetH;
  neuttarget_.numatom = nvect->TargetA;

  neutpiabs_.neabspiemit = 1; // (not filled, value from neutcore/necard.F)

  neutpiless_.ipilessdcy = nvect->IPilessDcy;
  neutpiless_.rpilessdcy = nvect->RPilessDcy;
}


void readnrcard(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  nucres_.nucrescat = nvect->NucScat;
  nucres_.xnucfact = nvect->NucFac;
  //nucres_.nucresflg; // (not filled)
}


void readver(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  neutversion_.corev = nvect->COREVer;
  neutversion_.nucev = nvect->NUCEVer;
  neutversion_.nuccv = nvect->NUCCVer;
}


void readneutcrs(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  neutcrscom_.crsenergy = nvect->CrsEnergy;
  neutcrscom_.totcrsne = nvect->Totcrs;
  for (int i=0; i<8; ++i) {
    neutcrscom_.difcrsne[i] = nvect->DifCrsNE[i];
  }
  neutcrscom_.crsx = nvect->Crsx;
  neutcrscom_.crsy = nvect->Crsy;
  neutcrscom_.crsz = nvect->Crsz;
  neutcrscom_.crsphi = nvect->Crsphi;
  neutcrscom_.crsq2 = nvect->Crsq2;
}


/**/
/*     COMMON/FSIHIST/ */
/**/
/*     NVERT       : # OF VERTICES*/
/*     POSVERT(3,I): POSITION OF I-TH VERTEX*/
/*     IFLGVERT(I) : INTERACTION TYPE OF I-TH VERTEX */
/*                     (*10 FOR HI-NRG)*/
/*                     (*100 for SKDETSIM non-NEFFECT interactions e.g. elastic SGPIEL;*/
/*                          +0 Free Hydrogen, +1 Oxygen)*/
/*                  -1 : ESCAPE*/
/*                   0 : INITIAL (or unmatched parent vertex if I>1)*/
/*                   1 :*/
/*                   2 : */
/*                   3 : ABSORPTION*/
/*                   4 : CHARGE EXCHANGE*/
/*                   5 : */
/*                   6 : */
/*                   7 : HADRON PRODUCTION (hi-nrg only, i.e. 70)*/
/*                   8 : QUASI-ELASTIC SCATTER*/
/*                   9 : FORWARD (ELASTIC-LIKE) SCATTER*/
/**/
/*     NVCVERT     : # OF INTERMEDIATE PARTICLES*/
/*     DIRVERT(3,J): DIRECTION OF J-TH PARTICLE*/
/*     ABSPVERT(J) : ABSOLUTE MOM. OF J-TH PART. in lab frame (MeV/c)*/
/*     ABSTPVERT(J): ABSOLUTE MOM. OF J-TH PART. in nucleon rest frame */
/*     IPVERT(J)   : PARTICLE CODE OF J-TH PARTICLE*/
/*     IVERTI(J)   : INDEX OF INITIAL VERTEX OF J-TH PARTICLE*/
/*     IVERTF(J)   : INDEX OF FINAL VERTEX OF J-TH PARTICLE*/
/**/
void readfsihist(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
  fsihist_.fsiprob = nvect->Fsiprob;

  NeutFsiPart* fsipinfo = nullptr;
  fsihist_.nvcvert = nvect->NfsiPart();
  if (fsihist_.nvcvert == 1) {
    // Check for dummy object
    fsipinfo = nvect->FsiPartInfo(0);
    if (   (fsipinfo->fPID == -1)
        && (fsipinfo->fMomLab = -1)
        && (fsipinfo->fMomNuc = -1)
        && (fsipinfo->fVertStart = -1)
        && (fsipinfo->fVertEnd = -1)
       )
    {
      fsihist_.nvcvert = 0;
    }
  }

  for (int i=0; i<fsihist_.nvcvert; ++i) {
    fsipinfo = nvect->FsiPartInfo(i);
    fsihist_.ipvert[i] = fsipinfo->fPID;
    fsihist_.dirvert[i][0] = fsipinfo->fDir.X();
    fsihist_.dirvert[i][1] = fsipinfo->fDir.Y();
    fsihist_.dirvert[i][2] = fsipinfo->fDir.Z();
    fsihist_.abspvert[i] = fsipinfo->fMomLab;
    fsihist_.abstpvert[i] = fsipinfo->fMomNuc;
    fsihist_.iverti[i] = fsipinfo->fVertStart;
    fsihist_.ivertf[i] = fsipinfo->fVertEnd;
  }

  NeutFsiVert* fsivinfo = nullptr;
  fsihist_.nvert = nvect->NfsiVert();
  if (fsihist_.nvert == 1) {
    fsivinfo = nvect->FsiVertInfo(0);
    if ((fsivinfo->fVertID == -99999) && (fsihist_.fsiprob == -1)) {
      fsihist_.nvert = 0;
    }
  }
  for (int i=0; i<fsihist_.nvert; ++i) {
    fsivinfo = nvect->FsiVertInfo(i);
    fsihist_.posvert[i][0] = fsivinfo->fPos.X();
    fsihist_.posvert[i][1] = fsivinfo->fPos.Y();
    fsihist_.posvert[i][2] = fsivinfo->fPos.Z();
    fsihist_.iflgvert[i] = fsivinfo->fVertID;
  }
}


/**/
/*       NFnvert        : number of vertices*/
/**/
/*       NFiflag(i)     : 4-digit flag for interaction type at i-th vertex, in the form "BNTP":*/
/*                         N: charge nucleon propagated through nucleus (0 = neutron, 1 = proton)*/
/*                         T: charge "target" nucleon the interaction is taking place on*/
/*                         P: scattering process:*/
/*                            P=0: start tracking of nucleon (i.e. gets "created")*/
/*                            P=1: elastic scattering*/
/*                            P=2: single pion production*/
/*                            P=3: double pion production*/
/*                            P=4: stop tracking of nucleon (i.e. leaves nucleus)*/
/*                         B: Pauli blocking flag (0 = not blocked, 1 = interaction was Pauli blocked*/
/*                            and actually did not take place)*/
/*                         Examples:*/
/*                          - 103 means double pion production when a proton scattered on a neutron*/
/*                          - 1011 means elastic scattering of a neutron on a proton did not take*/
/*                            place due to Pauli blocking*/
/*                         For P=0 and P=4, "T" is without meaning and always set to 0.*/
/*       NFx(i)         : x-component of i-th vertex position inside nucleus*/
/*       NFy(i)         : y-component of i-th vertex position inside nucleus*/
/*       NFz(i)         : z-component of i-th vertex position inside nucleus*/
/*       NFpx(i)        : x-component of momentum of nucleon leaving the i-th vertex*/
/*       NFpy(i)        : y-component of momentum of nucleon leaving the i-th vertex*/
/*       NFpz(i)        : z-component of momentum of nucleon leaving the i-th vertex*/
/*       NFe(i)         : energy of nucleon leaving the i-th vertex*/
/*       NFfirststep(i) : first step index of this track (to obtain the CMS energies for each step)*/
/**/
/*       NFnstep        : number of steps*/
/**/
/*       NFecms2(k)     : CMS energy squared of collision at k-th step (i.e. before interacting).*/
/*                        The sign of this value indicates the charge of the target nucleon:*/
/*                         NFecms2 > 0: proton,  NFecms2 < 0: neutron (same as "T" in NFiflag)*/
/*       NFptot(k)      : total probability at k-th step (for testing only, will be removed)*/
/**/
/*       Remarks:*/
/*        - a "vertex" is actually better described as a start, end or scattering point of a track*/
/*        - at each scattering point, the first nucleon will be followed in the same track, while the*/
/*          second one will create a new track*/
/*        - each track consists of a series of consecutive vertices. The first vertex has P=0, the*/
/*          last P=4. In between may be any number (including 0) vertices where an actual scattering*/
/*          took place (P=1,2,3).*/
/*        - it is not possible (and not needed) to connect the second track of a scattering vertex*/
/*          with the original one. Note that "first" and "second" is purely arbitrary. For nucleon*/
/*          FSI uncertainties, only the probabilities of the scattering processes have to be*/
/*          calculated, so it is not important to know which tracks belong to each other.*/
/**/
void readnucleonfsihist(NeutVect* nvect, NeutVtx* /*nvtx*/)
{
# if (NECOREVER == 5321)
  // special case for 5.3.2.1
  // do nothing.
# elif (NECOREVER == 5322) || (NECOREVER >= 535)
  //
  // special case for 5.3.2.2, which contains nucleon FSI
  //
  //if (nvect->Class_Version() < 4) {
  //  // Nucleon FSI is not stored until V4.
  //  return;
  //}

  NeutNucFsiVert* nucfsivinfo = nullptr;
  nucleonfsihist_.nfnvert = nvect->NnucFsiVert();
  if (nucleonfsihist_.nfnvert == 1) {
    nucfsivinfo = nvect->NucFsiVertInfo(0);
    if ((nucfsivinfo->fVertFlag == -1) && (nucfsivinfo->fVertFirstStep == -1)) {
      nucleonfsihist_.nfnvert = 0;
    }
  }
  for (int i=0; i<nucleonfsihist_.nfnvert; ++i) {
    nucfsivinfo = nvect->NucFsiVertInfo(i);
    nucleonfsihist_.nfiflag[i] = nucfsivinfo->fVertFlag;
    nucleonfsihist_.nffirststep[i] = nucfsivinfo->fVertFirstStep;
    nucleonfsihist_.nfx[i] = nucfsivinfo->fPos.X();
    nucleonfsihist_.nfy[i] = nucfsivinfo->fPos.Y();
    nucleonfsihist_.nfz[i] = nucfsivinfo->fPos.Z();
    nucleonfsihist_.nfpx[i] = nucfsivinfo->fMom.X();
    nucleonfsihist_.nfpy[i] = nucfsivinfo->fMom.Y();
    nucleonfsihist_.nfpz[i] = nucfsivinfo->fMom.Z();
    nucleonfsihist_.nfe[i] = nucfsivinfo->fMom.E();
  }

  NeutNucFsiStep* nucfsisinfo = nullptr;
  nucleonfsihist_.nfnstep = nvect->NnucFsiStep();
  if (nucleonfsihist_.nfnstep == 1) {
    nucfsisinfo = nvect->NucFsiStepInfo(0);
    if ((nucfsisinfo->fECMS2 == 0) && (nucfsisinfo->fProb == -1)) {
      nucleonfsihist_.nfnstep = 0;
    }
  }
  for (int i=0; i<nucleonfsihist_.nfnstep; ++i) {
    nucfsisinfo = nvect->NucFsiStepInfo(i);
    nucleonfsihist_.nfecms2[i] = nucfsisinfo->fECMS2;
    nucleonfsihist_.nfptot[i] = nucfsisinfo->fProb;
  }
# endif
}


/**/
/*     NVTXVC      : # OF VERTEX */
/*     PVTXVC(3,I) : POSITION OF VERTEX*/
/*     IFLVVC(I)   : KIND OF VERTEX*/
/*                   0 : DETERMINED LATER PROCEDURE*/
/*                   1 : PRIMARY POSITION*/
/*                   2 : DECAY*/
/*                   3 : STOP*/
/*     IPARVC(I)   : PARENT PARTICLE*/
/*     TIMVVC(I)   : TIME OF VERTEX*/
/**/
void readvcvrtx(NeutVect* /*nvect*/, NeutVtx* nvtx)
{
  TLorentzVector* pos = nullptr;
  vcvrtx_.nvtxvc = nvtx->Nvtx();
  for (int i=0; i<vcvrtx_.nvtxvc; ++i) {
    pos = nvtx->Pos(i);
    vcvrtx_.pvtxvc[i][0] = pos->X();
    vcvrtx_.pvtxvc[i][1] = pos->Y();
    vcvrtx_.pvtxvc[i][2] = pos->Z();
    vcvrtx_.timvvc[i]    = pos->T();
  }
}

void neut::PrintCommonBlocks() {
  using namespace std;
  cout << boolalpha << fixed << showpoint << noshowpos << setprecision(4);

  cout << "\n\n/COMMON NEWORK/\n";
  cout << setw(15)<< right << "nework_" << "." << setw(13) << left <<"modene"
       << "=" << setw(10) << right << nework_.modene << "\n";

  cout << setw(15)<< right << "nework_" << "." << setw(13) << left <<"numne"
       << "=" << setw(10) << right << nework_.numne << "\n";

  for (int i=0; i<vcwork_.nvc; ++i) {
    cout << "\n";
    cout << setw(15)<< right << "nework_" << "."
         << setw(6) << left <<"ipne" << "[" << setw(1) << right << i << setw(5) << left<< "]"
         << "=" << setw(10) << right << nework_.ipne[i] << "\n";

    cout << setw(15)<< right << "nework_" << "."
         << setw(6) << left <<"pne" <<"[" << setw(1) << right << i << setw(5) << left<< "][0]"
         << "=" << setw(10) << right << nework_.pne[i][0] << "\n";

    cout << setw(15)<< right << "nework_" << "."
         << setw(6) << left <<"pne" <<"[" << setw(1) << right << i << setw(5) << left<< "][1]"
         << "=" << setw(10) << right << nework_.pne[i][1] << "\n";

    cout << setw(15)<< right << "nework_" << "."
         << setw(6) << left <<"pne" <<"[" << setw(1) << right << i << setw(5) << left<< "][2]"
         << "=" << setw(10) << right << nework_.pne[i][2] << "\n";

    cout << setw(15)<< right << "nework_" << "."
         << setw(6) << left <<"iorgne" <<"[" << setw(1) << right << i << setw(5) << left<< "]"
         << "=" << setw(10) << right << nework_.iorgne[i] << "\n";

    cout << setw(15)<< right << "nework_" << "."
         << setw(6) << left <<"iflgne" <<"[" << setw(1) << right << i << setw(5) << left<< "]"
         << "=" << setw(10) << right << nework_.iflgne[i] << "\n";

    cout << setw(15)<< right << "nework_" << "."
         << setw(6) << left <<"icrnne" <<"[" << setw(1) << right << i << setw(5) << left<< "]"
         << "=" << setw(10) << right << nework_.icrnne[i] << "\n";
  }


  cout << "\n\n/COMMON NEMDLS/\n";
  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"mdlqe"
       << "=" << setw(10) << right << nemdls_.mdlqe << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"mdlqeaf"
       << "=" << setw(10) << right << nemdls_.mdlqeaf << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"sccfa"
       << "=" << setw(10) << right << nemdls_.sccfa << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"sccfv"
       << "=" << setw(10) << right << nemdls_.sccfv << "\n";

#if (NECOREVER != 5321) && (NECOREVER != 5322) && (NECOREVER >= 533)
  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"fpqe"
       << "=" << setw(10) << right << nemdls_.fpqe << "\n";
#endif

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"xmaqe"
       << "=" << setw(10) << right << nemdls_.xmaqe << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"xmvqe"
       << "=" << setw(10) << right << nemdls_.xmvqe << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"kapp"
       << "=" << setw(10) << right << nemdls_.kapp << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"xmaspi"
       << "=" << setw(10) << right << nemdls_.xmaspi << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"xmvspi"
       << "=" << setw(10) << right << nemdls_.xmvspi << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"xmares"
       << "=" << setw(10) << right << nemdls_.xmares << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"xmvres"
       << "=" << setw(10) << right << nemdls_.xmvres << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"xmacoh"
       << "=" << setw(10) << right << nemdls_.xmacoh << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"rad0nu"
       << "=" << setw(10) << right << nemdls_.rad0nu << "\n";

#if (NECOREVER != 5321) && (NECOREVER != 5322) && (NECOREVER >= 533)
  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"fa1coh"
       << "=" << setw(10) << right << nemdls_.fa1coh << "\n";

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"fb1coh"
       << "=" << setw(10) << right << nemdls_.fb1coh << "\n";
#endif

  cout << setw(15)<< right << "nemdls_" << "." << setw(13) << left <<"mdlcoh"
       << "=" << setw(10) << right << nemdls_.mdlcoh << "\n";


  cout << "\n\n/COMMON NENUPR/\n";
  cout << setw(15)<< right << "nenupr_" << "." << setw(13) << left <<"pfsurf"
       << "=" << setw(10) << right << nenupr_.pfsurf << "\n";

  cout << setw(15)<< right << "nenupr_" << "." << setw(13) << left <<"pfmax"
       << "=" << setw(10) << right << nenupr_.pfmax << "\n";

  cout << setw(15)<< right << "nenupr_" << "." << setw(13) << left <<"vnuini"
       << "=" << setw(10) << right << nenupr_.vnuini << "\n";

  cout << setw(15)<< right << "nenupr_" << "." << setw(13) << left <<"vnufin"
       << "=" << setw(10) << right << nenupr_.vnufin << "\n";

  cout << setw(15)<< right << "nenupr_" << "." << setw(13) << left <<"iformlen"
       << "=" << setw(10) << right << nenupr_.iformlen << "\n";

  cout << setw(15)<< right << "nenupr_" << "." << setw(13) << left <<"fzmu2"
       << "=" << setw(10) << right << nenupr_.fzmu2 << "\n";


  cout << "\n\n/COMMON NEFFPR/\n";
  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefqe"
       << "=" << setw(10) << right << neffpr_.fefqe << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefqeh"
       << "=" << setw(10) << right << neffpr_.fefqeh << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefinel"
       << "=" << setw(10) << right << neffpr_.fefinel << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefabs"
       << "=" << setw(10) << right << neffpr_.fefabs << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefcoh"
       << "=" << setw(10) << right << neffpr_.fefcoh << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefqehf"
       << "=" << setw(10) << right << neffpr_.fefqehf << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefcohf"
       << "=" << setw(10) << right << neffpr_.fefcohf << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefcx"
       << "=" << setw(10) << right << neffpr_.fefcx << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefcxhf"
       << "=" << setw(10) << right << neffpr_.fefcxhf << "\n";

  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefcxh"
       << "=" << setw(10) << right << neffpr_.fefcxh << "\n";

# if (NECOREVER != 5321) && (NECOREVER != 5322) && (NECOREVER >= 535)
  cout << setw(15)<< right << "neffpr_" << "." << setw(13) << left <<"fefall"
       << "=" << setw(10) << right << neffpr_.fefall << "\n";
# endif


  cout << "\n\n/COMMON NEUTTARGET/\n";
  cout << setw(15)<< right << "neuttarget_" << "." << setw(13) << left <<"numbndn"
       << "=" << setw(10) << right << neuttarget_.numbndn << "\n";

  cout << setw(15)<< right << "neuttarget_" << "." << setw(13) << left <<"numbndp"
       << "=" << setw(10) << right << neuttarget_.numbndp << "\n";

  cout << setw(15)<< right << "neuttarget_" << "." << setw(13) << left <<"numfrep"
       << "=" << setw(10) << right << neuttarget_.numfrep << "\n";

  cout << setw(15)<< right << "neuttarget_" << "." << setw(13) << left <<"numatom"
       << "=" << setw(10) << right << neuttarget_.numatom << "\n";


  // Input parameters
  cout << "\n\n-- Variable Model Parameters --\n";
  cout << "neut1pi_.xmanffres = " << neut1pi_.xmanffres << "\n";
  cout << "neut1pi_.xmvnffres = " << neut1pi_.xmvnffres << "\n";
  cout << "neut1pi_.xmarsres = "  << neut1pi_.xmarsres  << "\n";
  cout << "neut1pi_.xmvrsres = "  << neut1pi_.xmvrsres  << "\n";
  cout << "neut1pi_.neiff = "     << neut1pi_.neiff     << "\n";
  cout << "neut1pi_.nenrtype = "  << neut1pi_.nenrtype  << "\n";
  cout << "neut1pi_.rneca5i = "   << neut1pi_.rneca5i   << "\n";
  cout << "neut1pi_.rnebgscl = "  << neut1pi_.rnebgscl  << "\n";
  cout << "neutdis_.nepdf = "     << neutdis_.nepdf     << "\n";
  cout << "neutdis_.nebodek = "   << neutdis_.nebodek   << "\n";

  // Fixed parameters
  cout << "\n-- These are fixed parameters [should = ()] --\n";
  cout << "neutcard_.nefrmflg   =   (0)" << neutcard_.nefrmflg << "\n";
  cout << "neutcard_.nepauflg   =   (0)" << neutcard_.nepauflg << "\n";
  cout << "neutcard_.nenefo16   =   (0)" << neutcard_.nenefo16 << "\n";
  cout << "neutcard_.nemodflg   =   (0)" << neutcard_.nemodflg << "\n";
  cout << "neutcard_.nenefmodl  =   (1)" << neutcard_.nenefmodl << "\n";
  cout << "neutcard_.nenefmodh  =   (1)" << neutcard_.nenefmodh << "\n";
  cout << "neutcard_.nenefkinh  =   (1)" << neutcard_.nenefkinh << "\n";
  cout << "neutcard_.nusim      =   (1)" << neutcard_.nusim << "\n";

  cout << "neutpiabs_.neabspiemit = (1)"<< neutpiabs_.neabspiemit << "\n";

  cout << "neutcoh_.necohepi    =   (0)" <<  neutcoh_.necohepi << "\n";

  cout << "\n";
}

