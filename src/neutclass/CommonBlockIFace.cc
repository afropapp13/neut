
#include "CommonBlockIFace.h"

#include "cardnameC.h"

#include "fsihistC.h"
#include "neutcrsC.h"
#include "neworkC.h"
#include "nucleonfsihistC.h"
#include "vcworkC.h"

// Must come after vcworkC.h
#include "posinnucC.h"

#include <cstring>
#include <iostream>
#include <sstream>

extern "C" {
void necard_();
void necardev_();
}

namespace neut {

static CommonBlockIFace *gCommonBlockIFace = NULL;

void CommonBlockIFace::SetGenCard(std::string const &GenCardLocation) {
  if (!GenCardLocation.size()) {
    std::cout << "[ERROR]: CommonBlockIFace::SetGenCard passed empty string as "
                 "card location to read."
              << std::endl;
    throw;
  }
  // Set Card
  //Set the whole array to spaces first
  std::memset(necardname_.fcard, ' ', 1024);
  //Copy just the string characters without null termination as FORTRAN doesn't
  //do null termination...
  std::strncpy(necardname_.fcard, GenCardLocation.c_str(),
               std::min(size_t(1023), GenCardLocation.size()));
  necardname_.isset = true;
  // Read Card
  necard_();
  necardev_();

  // Call final common block setup functions
  NEUTSetParams();

  // Copy common blocks
  fneutcard_gen = neutcard_;
  fnuceffver_gen = nuceffver_;
  fneutdis_gen = neutdis_;
  fneut1pi_gen = neut1pi_;
  fneutdif_gen = neutdif_;
  fneutcoh_gen = neutcoh_;
  fneuttarget_gen = neuttarget_;
  fneutpiabs_gen = neutpiabs_;
  fneutpiless_gen = neutpiless_;
  fneutradcorr_gen = neutradcorr_;
  fnemdls_gen = nemdls_;
  fnenupr_gen = nenupr_;
  fneffpr_gen = neffpr_;
}
void CommonBlockIFace::Initialize(std::string const &GenCardLocation) {
  std::cout << "Initializing with " << GenCardLocation << std::endl;
  if (!gCommonBlockIFace) {
    gCommonBlockIFace = new CommonBlockIFace();
    gCommonBlockIFace->SetGenCard(GenCardLocation);
  }
}
CommonBlockIFace const &CommonBlockIFace::Get() {
  if (!gCommonBlockIFace) {
    std::cerr << "[ERROR]: Must manually call CommonBlockIFace::Initialize  "
                 "before first requesting an instance."
              << std::endl;
    abort();
  }
  return *gCommonBlockIFace;
}

void CommonBlockIFace::ResetGenValues() const {

  neutcard_ = fneutcard_gen;
  nuceffver_ = fnuceffver_gen;
  neutdis_ = fneutdis_gen;
  neut1pi_ = fneut1pi_gen;
  neutdif_ = fneutdif_gen;
  neutcoh_ = fneutcoh_gen;
  neuttarget_ = fneuttarget_gen;
  neutpiabs_ = fneutpiabs_gen;
  neutpiless_ = fneutpiless_gen;
  neutradcorr_ = fneutradcorr_gen;
  nemdls_ = fnemdls_gen;
  nenupr_ = fnenupr_gen;
  neffpr_ = fneffpr_gen;

  NEUTSetParams();
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
void CommonBlockIFace::ReadNEWORK(NeutVect *nvect) {
  nework_.modene = nvect->Mode;
  nework_.numne = nvect->Nprimary();

  const float MeV_per_GeV = 1.e3;
  NeutPart *pinfo = nullptr;
  int npart = nvect->Npart();
  for (int i = 0; i < npart; ++i) {
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

void CommonBlockIFace::ReadNEUTCRS(NeutVect *nvect) {
  neutcrscom_.crsenergy = nvect->CrsEnergy;
  neutcrscom_.totcrsne = nvect->Totcrs;
  for (int i = 0; i < 8; ++i) {
    neutcrscom_.difcrsne[i] = nvect->DifCrsNE[i];
  }
  neutcrscom_.crsx = nvect->Crsx;
  neutcrscom_.crsy = nvect->Crsy;
  neutcrscom_.crsz = nvect->Crsz;
  neutcrscom_.crsphi = nvect->Crsphi;
  neutcrscom_.crsq2 = nvect->Crsq2;
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
void CommonBlockIFace::ReadVCWORK(NeutVect *nvect) {
  NeutPart *pinfo = nullptr;
  vcwork_.nvc = nvect->Npart();
  for (int i = 0; i < vcwork_.nvc; ++i) {
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

/**/
/*     COMMON/FSIHIST/ */
/**/
/*     NVERT       : # OF VERTICES*/
/*     POSVERT(3,I): POSITION OF I-TH VERTEX*/
/*     IFLGVERT(I) : INTERACTION TYPE OF I-TH VERTEX */
/*                     (*10 FOR HI-NRG)*/
/*                     (*100 for SKDETSIM non-NEFFECT interactions e.g. elastic
 * SGPIEL;*/
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
void CommonBlockIFace::ReadFSIHIST(NeutVect *nvect) {
  fsihist_.fsiprob = nvect->Fsiprob;

  NeutFsiPart *fsipinfo = nullptr;
  fsihist_.nvcvert = nvect->NfsiPart();
  if (fsihist_.nvcvert == 1) {
    // Check for dummy object
    fsipinfo = nvect->FsiPartInfo(0);
    if ((fsipinfo->fPID == -1) && (fsipinfo->fMomLab = -1) &&
        (fsipinfo->fMomNuc = -1) && (fsipinfo->fVertStart = -1) &&
        (fsipinfo->fVertEnd = -1)) {
      fsihist_.nvcvert = 0;
    }
  }

  for (int i = 0; i < fsihist_.nvcvert; ++i) {
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

  NeutFsiVert *fsivinfo = nullptr;
  fsihist_.nvert = nvect->NfsiVert();
  if (fsihist_.nvert == 1) {
    fsivinfo = nvect->FsiVertInfo(0);
    if ((fsivinfo->fVertID == -99999) && (fsihist_.fsiprob == -1)) {
      fsihist_.nvert = 0;
    }
  }
  for (int i = 0; i < fsihist_.nvert; ++i) {
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
/*       NFiflag(i)     : 4-digit flag for interaction type at i-th vertex, in
 * the form "BNTP":*/
/*                         N: charge nucleon propagated through nucleus (0 =
 * neutron, 1 = proton)*/
/*                         T: charge "target" nucleon the interaction is taking
 * place on*/
/*                         P: scattering process:*/
/*                            P=0: start tracking of nucleon (i.e. gets
 * "created")*/
/*                            P=1: elastic scattering*/
/*                            P=2: single pion production*/
/*                            P=3: double pion production*/
/*                            P=4: stop tracking of nucleon (i.e. leaves
 * nucleus)*/
/*                         B: Pauli blocking flag (0 = not blocked, 1 =
 * interaction was Pauli blocked*/
/*                            and actually did not take place)*/
/*                         Examples:*/
/*                          - 103 means double pion production when a proton
 * scattered on a neutron*/
/*                          - 1011 means elastic scattering of a neutron on a
 * proton did not take*/
/*                            place due to Pauli blocking*/
/*                         For P=0 and P=4, "T" is without meaning and always
 * set to 0.*/
/*       NFx(i)         : x-component of i-th vertex position inside nucleus*/
/*       NFy(i)         : y-component of i-th vertex position inside nucleus*/
/*       NFz(i)         : z-component of i-th vertex position inside nucleus*/
/*       NFpx(i)        : x-component of momentum of nucleon leaving the i-th
 * vertex*/
/*       NFpy(i)        : y-component of momentum of nucleon leaving the i-th
 * vertex*/
/*       NFpz(i)        : z-component of momentum of nucleon leaving the i-th
 * vertex*/
/*       NFe(i)         : energy of nucleon leaving the i-th vertex*/
/*       NFfirststep(i) : first step index of this track (to obtain the CMS
 * energies for each step)*/
/**/
/*       NFnstep        : number of steps*/
/**/
/*       NFecms2(k)     : CMS energy squared of collision at k-th step (i.e.
 * before interacting).*/
/*                        The sign of this value indicates the charge of the
 * target nucleon:*/
/*                         NFecms2 > 0: proton,  NFecms2 < 0: neutron (same as
 * "T" in NFiflag)*/
/*       NFptot(k)      : total probability at k-th step (for testing only, will
 * be removed)*/
/**/
/*       Remarks:*/
/*        - a "vertex" is actually better described as a start, end or
 * scattering point of a track*/
/*        - at each scattering point, the first nucleon will be followed in the
 * same track, while the*/
/*          second one will create a new track*/
/*        - each track consists of a series of consecutive vertices. The first
 * vertex has P=0, the*/
/*          last P=4. In between may be any number (including 0) vertices where
 * an actual scattering*/
/*          took place (P=1,2,3).*/
/*        - it is not possible (and not needed) to connect the second track of a
 * scattering vertex*/
/*          with the original one. Note that "first" and "second" is purely
 * arbitrary. For nucleon*/
/*          FSI uncertainties, only the probabilities of the scattering
 * processes have to be*/
/*          calculated, so it is not important to know which tracks belong to
 * each other.*/
/**/
void CommonBlockIFace::ReadNUCLEONFSIHIST(NeutVect *nvect) {

  NeutNucFsiVert *nucfsivinfo = nullptr;
  nucleonfsihist_.nfnvert = nvect->NnucFsiVert();
  if (nucleonfsihist_.nfnvert == 1) {
    nucfsivinfo = nvect->NucFsiVertInfo(0);
    if ((nucfsivinfo->fVertFlag == -1) && (nucfsivinfo->fVertFirstStep == -1)) {
      nucleonfsihist_.nfnvert = 0;
    }
  }
  for (int i = 0; i < nucleonfsihist_.nfnvert; ++i) {
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

  NeutNucFsiStep *nucfsisinfo = nullptr;
  nucleonfsihist_.nfnstep = nvect->NnucFsiStep();
  if (nucleonfsihist_.nfnstep == 1) {
    nucfsisinfo = nvect->NucFsiStepInfo(0);
    if ((nucfsisinfo->fECMS2 == 0) && (nucfsisinfo->fProb == -1)) {
      nucleonfsihist_.nfnstep = 0;
    }
  }
  for (int i = 0; i < nucleonfsihist_.nfnstep; ++i) {
    nucfsisinfo = nvect->NucFsiStepInfo(i);
    nucleonfsihist_.nfecms2[i] = nucfsisinfo->fECMS2;
    nucleonfsihist_.nfptot[i] = nucfsisinfo->fProb;
  }
}

/**
 *     COMMON /POSINNUC/
 *     =================================================================
 *     IBOUND : INTERACTION ON BOUND NUCLEON OR FREE NUCLEON
 *
 *     POSNUC(3,I) : INITIAL VERTEX POSITION
 **/
void CommonBlockIFace::ReadPOSINNUC(NeutVect *nvect) {
  posinnuc_.ibound = nvect->Ibound;

  NeutPart *pinfo = nullptr;
  int npart = nvect->Npart();
  for (int i = 0; i < npart; ++i) {
    pinfo = nvect->PartInfo(i);
    posinnuc_.posnuc[i][0] = pinfo->fPosIni.X();
    posinnuc_.posnuc[i][1] = pinfo->fPosIni.Y();
    posinnuc_.posnuc[i][2] = pinfo->fPosIni.Z();
  }
}

void CommonBlockIFace::ReadVect(NeutVect *nvect) {
  ReadNEWORK(nvect);
  ReadVCWORK(nvect);
  ReadPOSINNUC(nvect);
  ReadNEUTCRS(nvect);
  ReadFSIHIST(nvect);
  ReadNUCLEONFSIHIST(nvect);
  ReadPARS(nvect);
}

  // Read parameters that may be important for reweight work and aren't filled in necard
void CommonBlockIFace::ReadPARS(NeutVect *nvect) {
  nemdls_.spidelta = nvect->SPIDelta;
}

std::string CommonBlockIFace::ParamsToString(bool isinstance) {
  std::stringstream ss("");

  neutcard_common const &neutcard =
      isinstance ? CommonBlockIFace::Get().fneutcard_gen : neutcard_;

  ss << "NEUTCARD: \n"
     << "\tnefrmflg:" << neutcard.nefrmflg << "\n"
     << "\tnepauflg:" << neutcard.nepauflg << "\n"
     << "\tnenefo16:" << neutcard.nenefo16 << "\n"
     << "\tnenefmodl:" << neutcard.nenefmodl << "\n"
     << "\tnenefmodh:" << neutcard.nenefmodh << "\n"
     << "\tnenefkinh:" << neutcard.nenefkinh << "\n"
     << "\tnemodflg:" << neutcard.nemodflg << "\n"
     << "\tneselmod:" << neutcard.neselmod << "\n"
     << "\tcrsneut: [";
  for (size_t i = 0; i < 30; ++i) {
    ss << neutcard.crsneut[i] << ((i != 29) ? ", " : " ");
  }
  ss << "]\n"
     << "\tcrsneutb: [";
  for (size_t i = 0; i < 30; ++i) {
    ss << neutcard.crsneutb[i] << ((i != 29) ? ", " : " ");
  }
  ss << "]\n"
     << "\titauflgcore:" << neutcard.itauflgcore << "\n"
     << "\tnusim:" << neutcard.nusim << "\n"
     << "\tquiet:" << neutcard.quiet << "\n";

  nuceffver_common const &nuceffver =
      isinstance ? CommonBlockIFace::Get().fnuceffver_gen : nuceffver_;

  ss << "NUCEFFVER:\n"
     << "\tnefkinver: " << nuceffver.nefkinver << "\n";

  neutdis_common const &neutdis =
      isinstance ? CommonBlockIFace::Get().fneutdis_gen : neutdis_;

  ss << "NEUTDIS:\n"
     << "\tnepdf: " << neutdis.nepdf << "\n"
     << "\tnebodek: " << neutdis.nebodek << "\n"
     << "\tnemult: " << neutdis.nemult << "\n";

  neut1pi_common const &neut1pi =
      isinstance ? CommonBlockIFace::Get().fneut1pi_gen : neut1pi_;

  ss << "NEUT1pi:\n"
     << "\txmanffres: " << neut1pi.xmanffres << "\n"
     << "\txmvnffres: " << neut1pi.xmvnffres << "\n"
     << "\txmarsres: " << neut1pi.xmarsres << "\n"
     << "\txmvrsres: " << neut1pi.xmvrsres << "\n"
     << "\tneiff: " << neut1pi.neiff << "\n"
     << "\tnenrtype: " << neut1pi.nenrtype << "\n"
     << "\trneca5i: " << neut1pi.rneca5i << "\n"
     << "\trnebgscl: " << neut1pi.rnebgscl << "\n";

  neutdif_common const &neutdif =
      isinstance ? CommonBlockIFace::Get().fneutdif_gen : neutdif_;

  ss << "NEUTDIF:\n"
     << "\tnedifpi: " << neutdif.nedifpi << "\n";

  neutcoh_common const &neutcoh =
      isinstance ? CommonBlockIFace::Get().fneutcoh_gen : neutcoh_;

  ss << "NEUTCOH:\n"
     << "\tnecohepi: " << neutcoh.necohepi << "\n";

  neuttarget_common const &neuttarget =
      isinstance ? CommonBlockIFace::Get().fneuttarget_gen : neuttarget_;

  ss << "NEUTTARGET:\n"
     << "\tnumbndn: " << neuttarget.numbndn << "\n"
     << "\tnumbndp: " << neuttarget.numbndp << "\n"
     << "\tnumfrep: " << neuttarget.numfrep << "\n"
     << "\tnumatom: " << neuttarget.numatom << "\n";

  neutpiabs_common const &neutpiabs =
      isinstance ? CommonBlockIFace::Get().fneutpiabs_gen : neutpiabs_;

  ss << "NEUTPIABS:\n"
     << "\tneabspiemit: " << neutpiabs.neabspiemit << "\n";

  neutpiless_common const &neutpiless =
      isinstance ? CommonBlockIFace::Get().fneutpiless_gen : neutpiless_;

  ss << "NEUTPILESS:\n"
     << "\tipilessdcy: " << neutpiless.ipilessdcy << "\n"
     << "\trpilessdcy: " << neutpiless.rpilessdcy << "\n";

  neutradcorr_common const &neutradcorr =
      isinstance ? CommonBlockIFace::Get().fneutradcorr_gen : neutradcorr_;

  ss << "neutradcorr:\n"
     << "\tiradcorr: " << neutradcorr.iradcorr << "\n";

  nemdls_common const &nemdls =
      isinstance ? CommonBlockIFace::Get().fnemdls_gen : nemdls_;

  ss << "nemdls:\n"
     << "\tmdlqe: " << nemdls.mdlqe << "\n"
     << "\tmdlspi: " << nemdls.mdlspi << "\n"
     << "\tmdlspiej: " << nemdls.mdlspiej << "\n"
     << "\tmdldis: " << nemdls.mdldis << "\n"
     << "\tmdlcoh: " << nemdls.mdlcoh << "\n"
     << "\tmdldif: " << nemdls.mdldif << "\n"
     << "\tmdlqeaf: " << nemdls.mdlqeaf << "\n"
     << "\txmaqe: " << nemdls.xmaqe << "\n"
     << "\txmaspi: " << nemdls.xmaspi << "\n"
     << "\txmvqe: " << nemdls.xmvqe << "\n"
     << "\txmvspi: " << nemdls.xmvspi << "\n"
     << "\tkapp: " << nemdls.kapp << "\n"
     << "\txmacoh: " << nemdls.xmacoh << "\n"
     << "\trad0nu: " << nemdls.rad0nu << "\n"
     << "\tfa1coh: " << nemdls.fa1coh << "\n"
     << "\tfb1coh: " << nemdls.fb1coh << "\n"
     << "\tiffspi: " << nemdls.iffspi << "\n"
     << "\tnrtypespi: " << nemdls.nrtypespi << "\n"
     << "\trca5ispi: " << nemdls.rca5ispi << "\n"
     << "\trbgsclspi: " << nemdls.rbgsclspi << "\n"
     << "\txmares: " << nemdls.xmares << "\n"
     << "\txmvres: " << nemdls.xmvres << "\n"
     << "\tsccfv: " << nemdls.sccfv << "\n"
     << "\tsccfa: " << nemdls.sccfa << "\n"
     << "\tfpqe: " << nemdls.fpqe << "\n"
     << "\tpfsf: " << nemdls.pfsf << "\n"
     << "\txmadif: " << nemdls.xmadif << "\n"
     << "\tnucvoldif: " << nemdls.nucvoldif << "\n"
     << "\taxffalpha: " << nemdls.axffalpha << "\n"
     << "\taxffgamma: " << nemdls.axffgamma << "\n"
     << "\taxfftheta: " << nemdls.axfftheta << "\n"
     << "\taxffbeta: " << nemdls.axffbeta << "\n"
     << "\taxzexpq4: " << nemdls.axzexpq4 << "\n"
     << "\taxzexpnt: " << nemdls.axzexpnt << "\n"
     << "\taxzexpt0: " << nemdls.axzexpt0 << "\n"
     << "\taxzexptc: " << nemdls.axzexptc << "\n"
     << "\taxzexpa0: " << nemdls.axzexpa0 << "\n"
     << "\taxzexpa1: " << nemdls.axzexpa1 << "\n"
     << "\taxzexpa2: " << nemdls.axzexpa2 << "\n"
     << "\taxzexpa3: " << nemdls.axzexpa3 << "\n"
     << "\taxzexpa4: " << nemdls.axzexpa4 << "\n"
     << "\taxzexpa5: " << nemdls.axzexpa5 << "\n"
     << "\taxzexpa6: " << nemdls.axzexpa6 << "\n"
     << "\taxzexpa7: " << nemdls.axzexpa7 << "\n"
     << "\taxzexpa8: " << nemdls.axzexpa8 << "\n"
     << "\taxzexpa9: " << nemdls.axzexpa9 << "\n"
     << "\tmdl2p2h: " << nemdls.mdl2p2h << "\n"
     << "\txmancel: " << nemdls.xmancel << "\n";

  nenupr_common const &nenupr =
      isinstance ? CommonBlockIFace::Get().fnenupr_gen : nenupr_;

  ss << "nenupr:\n"
     << "\tpfsurf: " << nenupr.pfsurf << "\n"
     << "\tpfmax: " << nenupr.pfmax << "\n"
     << "\tvnuini: " << nenupr.vnuini << "\n"
     << "\tvnufin: " << nenupr.vnufin << "\n"
     << "\tiformlen: " << nenupr.iformlen << "\n"
     << "\tfzmu2: " << nenupr.fzmu2 << "\n"
     << "\tsfebshift: " << nenupr.sfebshift << "\n";

  neffpr_common const &neffpr =
      isinstance ? CommonBlockIFace::Get().fneffpr_gen : neffpr_;

  ss << "neffpr:\n"
     << "\tfefqe: " << neffpr.fefqe << "\n"
     << "\tfefqeh: " << neffpr.fefqeh << "\n"
     << "\tfefinel: " << neffpr.fefinel << "\n"
     << "\tfefabs: " << neffpr.fefabs << "\n"
     << "\tfefcoh: " << neffpr.fefcoh << "\n"
     << "\tfefqehf: " << neffpr.fefqehf << "\n"
     << "\tfefcohf: " << neffpr.fefcohf << "\n"
     << "\tfefcx: " << neffpr.fefcx << "\n"
     << "\tfefcxhf: " << neffpr.fefcxhf << "\n"
     << "\tfefcxh: " << neffpr.fefcxh << "\n"
     << "\tfefcoul: " << neffpr.fefcoul << "\n"
     << "\tfefall: " << neffpr.fefall << "\n";

  return ss.str();
}

} // namespace neut
