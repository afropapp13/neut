#include "NeutRootHandlers.h"

#include "necardC.h"
#include "neutcrsC.h"
#include "neutmodelC.h"
#include "neworkC.h"
#include "vcvrtxC.h"
#include "vcworkC.h"

#include "fsihistC.h"
#include "nucleonfsihistC.h"

#include "posinnucC.h"

#include "neutmodesCPP.h"

#include "CommonBlockIFace.h"
#include "HepMCReader.h"

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"

#include "NuHepMC/Print.hxx"
#include "NuHepMC/WriterTools"

#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <type_traits>

std::map<int, std::unique_ptr<NuHepMC::WriterRootTree> > HepMCWriters;
std::map<int, size_t> EventCounts;
int id = 0;

extern "C" {
int hepmcopen_(char *filename, int *filename_len, double *FATC);
void hepmcclose_(int *fid);
void hepmcfill_(int *fid);
void mcmass_(int *pid, float *massout);
};

void nulltermstr(char *str, int len) {
  int i;
  for (i = 0; i < len; i++) {
    if (str[i] == ' ') {
      str[i] = '\0';
      break;
    }
  }
}

int hepmcopen_(char *filename, int *filename_len, double *FATC) {

  nulltermstr(filename, *filename_len);

  EventCounts[id] = 0;
  std::cout << "[INFO]: Opened file: " << filename << ", at ID: " << id
            << std::endl;

  std::shared_ptr<HepMC3::GenRunInfo> gri = NuHepMC::genruninfo::GRIFactory(
      "NEUT", VERSION, "The NEUT neutrino interaction generator");

  NuHepMC::genruninfo::SetVertexEnumStandard(gri, "1.0");
  NuHepMC::genruninfo::SetParticleEnumStandard(gri, "1.0");

  std::map<int, std::string> ModeNames;

  for (size_t i = 0; i < kNneutModes; ++i) {
    if (!neutModeID[i]) { // Ignore Total
      continue;
    }
    ModeNames[neutModeID[i]] = std::string("neutrino ") + neutModeTitle[i];
    ModeNames[-neutModeID[i]] =
        std::string("anti-neutrino ") + neutModeTitle[i];
  }

  NuHepMC::genruninfo::SetHardScatterModeDefinitions(gri, ModeNames);

  NuHepMC::genruninfo::SetExtraVertexEnumDefinitions(
      gri,
      {{NuHepMC::labels::e2i(neut::VertexState::kNEUTVertex), "NEUTVertex"},
       {NuHepMC::labels::e2i(neut::VertexState::kPionFSIVertex),
        "PionFSIVertex"},
       {NuHepMC::labels::e2i(neut::VertexState::kNucleonFSIVertex),
        "NucleonFSIVertex"}});

  NuHepMC::genruninfo::SetExtraParticleEnumDefinitions(
      gri,
      {
          // Particle stack and pion FSI history
          {NuHepMC::labels::e2i(neut::ParticleState::kChargeExchange),
           "ChargeExchange"},
          {NuHepMC::labels::e2i(neut::ParticleState::kBlocked), "Blocked"},
          {NuHepMC::labels::e2i(neut::ParticleState::kEMShower), "EMShower"},
          {NuHepMC::labels::e2i(neut::ParticleState::kHadronProduction),
           "HadronProduction"},
          {NuHepMC::labels::e2i(neut::ParticleState::kQEScatter), "QEScatter"},
          {NuHepMC::labels::e2i(neut::ParticleState::kForwardElasticScatter),
           "ForwardElasticScatter"},

          // Nucleon FSI history
          {NuHepMC::labels::e2i(neut::ParticleState::kNucleonElastic),
           "NucleonElastic"},
          {NuHepMC::labels::e2i(neut::ParticleState::kNucleonPiProd),
           "NucleonPiProd"},
          {NuHepMC::labels::e2i(neut::ParticleState::kNucleonDblPiProd),
           "NucleonDblPiProd"},

      });

  NuHepMC::genruninfo::SetFluxAveragedTotalCrossSection(gri, *FATC);

  // Add the common block 'card' values
  for (auto const &attr :
       neut::CommonBlockIFace::SerializeModelCommonBlocksToAttributeList()) {
    gri->add_attribute(attr.first, attr.second);
  }

  HepMCWriters[id] = std::unique_ptr<NuHepMC::WriterRootTree>(
      new NuHepMC::WriterRootTree(filename, gri));

  return id++;
}

void hepmcclose_(int *fid) {
  HepMCWriters[*fid]->close();
  HepMCWriters.erase(*fid);
}
/*
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
 */

int GetParticleStatusEnum(int status, int chase) {
  if (chase == 0) {
    switch (status) {
    case -1: {
      return NuHepMC::labels::e2i(
          NuHepMC::labels::ParticleState::kInitialState);
    }
    case 1: {
      return NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kDecayed);
    }
    case 3: {
      return NuHepMC::labels::e2i(neut::ParticleState::kAbsorption);
    }
    case 4: {
      return NuHepMC::labels::e2i(neut::ParticleState::kChargeExchange);
    }
    case 5: {
      return NuHepMC::labels::e2i(neut::ParticleState::kBlocked);
    }
    case 6: {
      return NuHepMC::labels::e2i(neut::ParticleState::kEMShower);
    }
    case 7: {
      return NuHepMC::labels::e2i(neut::ParticleState::kHadronProduction);
    }
    case 8: {
      return NuHepMC::labels::e2i(neut::ParticleState::kQEScatter);
    }
    case 9: {
      return NuHepMC::labels::e2i(neut::ParticleState::kForwardElasticScatter);
    }
    default: {
      return NuHepMC::labels::e2i(
          NuHepMC::labels::ParticleState::kIntermediate);
    }
    }
  } else if (chase == 1) {
    switch (status) {
    case 0:
    case 2: {
      return NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kFinalState);
    }
    default: {
      return NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kOther);
    }
    }

  } else {
    return NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kIntermediate);
  }
}

int GetNucleonFSIStatusEnum(int B, int P) {
  if (B) { // Pauli blocked
    return NuHepMC::labels::e2i(neut::ParticleState::kBlocked);
  }

  // P=0: start tracking of nucleon (i.e. gets "created")
  //  * !     P=1: elastic scattering
  //  * !     P=2: single pion production
  //  * !     P=3: double pion production
  //  * !     P=4: stop tracking of nucleon (i.e. leaves nucleus)

  switch (P) {
  case 0: {
    return NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kInitialState);
  }
  case 1: {
    return NuHepMC::labels::e2i(neut::ParticleState::kNucleonElastic);
  }
  case 2: {
    return NuHepMC::labels::e2i(neut::ParticleState::kNucleonPiProd);
  }
  case 3: {
    return NuHepMC::labels::e2i(neut::ParticleState::kNucleonDblPiProd);
  }
  case 4: {
    return NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kFinalState);
  }
  default: {
    return NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kIntermediate);
  }
  }
}

//let the compiler get the types correct
template <typename T, size_t N>
void FillVectorAttribute(HepMC3::GenEvent &evt, std::string const &name,
                         T (&t)[N], size_t n) {
  std::vector<T> vect(n);
  std::copy_n(t, n, vect.begin());
  evt.add_attribute(name, NuHepMC::AsAttribute(vect));
}

template <typename T, size_t N, size_t stride>
void FillVectorAttribute(HepMC3::GenEvent &evt, std::string const &name,
                         T (&t)[N][stride], size_t n) {
  std::vector<T> vect(n * stride);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < stride; ++j) {
      vect[(i * stride) + j] = t[i][j];
    }
  }

  evt.add_attribute(name, NuHepMC::AsAttribute(vect));
}

void hepmcfill_(int *fid) {

  HepMC3::GenEvent evt(HepMC3::Units::MEV, HepMC3::Units::CM);
  evt.set_event_number(EventCounts[*fid]++);

  NuHepMC::genevent::SetHardScatterMode(evt, nework_.modene);

  // Build the NuHepMC Standard format
  HepMC3::GenVertexPtr LabFrameVtx =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector::ZERO_VECTOR());
  LabFrameVtx->set_status(
      NuHepMC::labels::e2i(NuHepMC::labels::VertexState::kLabFrame));

  if (posinnuc_.ibound) {
    int A = neuttarget_.numbndp + neuttarget_.numbndn;
    LabFrameVtx->add_particle_in(
        NuHepMC::genevent::MakeNuclearParticle(neuttarget_.numbndp, A));
  } else {
    LabFrameVtx->add_particle_in(NuHepMC::genevent::MakeNuclearParticle(1, 1));
  }

  HepMC3::GenVertexPtr HardScatterVtx =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector::ZERO_VECTOR());
  HardScatterVtx->set_status(
      NuHepMC::labels::e2i(NuHepMC::labels::VertexState::kHardScatter));

  std::vector<HepMC3::GenVertexPtr> vertices;

  // Add the NEUT vector layout
  for (int vtx_it = 0; vtx_it < vcvrtx_.nvtxvc; vtx_it++) {
    vertices.push_back(std::make_shared<HepMC3::GenVertex>(
        HepMC3::FourVector(vcvrtx_.pvtxvc[vtx_it][0], vcvrtx_.pvtxvc[vtx_it][1],
                           vcvrtx_.pvtxvc[vtx_it][2], vcvrtx_.timvvc[vtx_it])));
    vertices.back()->set_status(
        NuHepMC::labels::e2i(neut::VertexState::kNEUTVertex));
  }

  for (int p_it = 0; p_it < vcwork_.nvc; p_it++) {

    // Fortran to C indexing
    int vtxid = vcwork_.ivtivc[p_it] - 1;
    if (vertices.size() <= vtxid) {
      std::cout << "[ERROR]: Trying to add particle to non-existent vertex ("
                << vtxid << "/" << vertices.size() << ")." << std::endl;
      abort();
    }

    double p_E =
        sqrt(pow(vcwork_.amasvc[p_it], 2) + pow(vcwork_.pvc[p_it][0], 2) +
             pow(vcwork_.pvc[p_it][1], 2) + pow(vcwork_.pvc[p_it][2], 2));

    int state =
        GetParticleStatusEnum(vcwork_.iflgvc[p_it], vcwork_.icrnvc[p_it]);

    HepMC3::GenParticlePtr p1 = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(vcwork_.pvc[p_it][0], vcwork_.pvc[p_it][1],
                           vcwork_.pvc[p_it][2], p_E),
        vcwork_.ipvc[p_it], state);
    p1->set_generated_mass(vcwork_.amasvc[p_it]);

    switch (state) {
    case NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kInitialState): {
      // Fill out hard scatter vertex.
      if (NuHepMC::pid::IsLepton(p1->pid())) {
        LabFrameVtx->add_particle_in(
            std::make_shared<HepMC3::GenParticle>(p1->data()));
      }
      // Only put 'primary' (NEWORK) particles into the hard scatter vtx
      if (p_it < nework_.numne) {
        HardScatterVtx->add_particle_in(
            std::make_shared<HepMC3::GenParticle>(p1->data()));
      }
      vertices[vtxid]->add_particle_in(p1);
      break;
    }
    case NuHepMC::labels::e2i(NuHepMC::labels::ParticleState::kFinalState): {
      LabFrameVtx->add_particle_out(
          std::make_shared<HepMC3::GenParticle>(p1->data()));
      if (p_it < nework_.numne) {
        HardScatterVtx->add_particle_out(
            std::make_shared<HepMC3::GenParticle>(p1->data()));
      }
      vertices[vtxid]->add_particle_out(p1);

      break;
    }
    default: {
      if (p_it < nework_.numne) {
        HardScatterVtx->add_particle_out(
            std::make_shared<HepMC3::GenParticle>(p1->data()));
      }
      vertices[vtxid]->add_particle_out(p1);

      break;
    }
    }
  }

  evt.add_vertex(LabFrameVtx);
  evt.add_vertex(HardScatterVtx);
  for (auto &vtx : vertices) {
    evt.add_vertex(vtx);
  }

  // Add generator specific passthrough for reweighting
  evt.add_attribute("Totcrs", NuHepMC::AsAttribute(neutcrscom_.totcrsne));
  evt.add_attribute("CrsEnergy", NuHepMC::AsAttribute(neutcrscom_.crsenergy));

  FillVectorAttribute(evt, "DifCrsNE", neutcrscom_.difcrsne, 8);

  evt.add_attribute("Crsx", NuHepMC::AsAttribute(neutcrscom_.crsx));
  evt.add_attribute("Crsy", NuHepMC::AsAttribute(neutcrscom_.crsy));
  evt.add_attribute("Crsz", NuHepMC::AsAttribute(neutcrscom_.crsz));
  evt.add_attribute("Crsphi", NuHepMC::AsAttribute(neutcrscom_.crsphi));
  evt.add_attribute("Crsq2", NuHepMC::AsAttribute(neutcrscom_.crsq2));

  evt.add_attribute("numbndn", NuHepMC::AsAttribute(neuttarget_.numbndn));
  evt.add_attribute("numbndp", NuHepMC::AsAttribute(neuttarget_.numbndp));
  evt.add_attribute("numfrep", NuHepMC::AsAttribute(neuttarget_.numfrep));
  evt.add_attribute("numatom", NuHepMC::AsAttribute(neuttarget_.numatom));

  evt.add_attribute("pfsurf", NuHepMC::AsAttribute(nenupr_.pfsurf));
  evt.add_attribute("pfmax", NuHepMC::AsAttribute(nenupr_.pfmax));
  evt.add_attribute("vnuini", NuHepMC::AsAttribute(nenupr_.vnuini));
  evt.add_attribute("vnufin", NuHepMC::AsAttribute(nenupr_.vnufin));

  evt.add_attribute("ibound", NuHepMC::AsAttribute(posinnuc_.ibound));

  // Want to write all of the event common block info straight into the event
  // attributes.

  /**/
  /*     COMMON /VCWORK/ */
  /*     ================================================================= */
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
  evt.add_attribute("nvtxvc", NuHepMC::AsAttribute(vcvrtx_.nvtxvc));
  FillVectorAttribute(evt, "pvtxvc", vcvrtx_.pvtxvc, vcvrtx_.nvtxvc);
  FillVectorAttribute(evt, "iflvvc", vcvrtx_.iflvvc, vcvrtx_.nvtxvc);
  FillVectorAttribute(evt, "iparvc", vcvrtx_.iparvc, vcvrtx_.nvtxvc);
  FillVectorAttribute(evt, "timvvc", vcvrtx_.timvvc, vcvrtx_.nvtxvc);

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
  evt.add_attribute("nvc", NuHepMC::AsAttribute(vcwork_.nvc));
  FillVectorAttribute(evt, "posvc", vcwork_.posvc, 3);
  FillVectorAttribute(evt, "ipvc", vcwork_.ipvc, vcwork_.nvc);
  FillVectorAttribute(evt, "amasvc", vcwork_.amasvc, vcwork_.nvc);
  FillVectorAttribute(evt, "pvc", vcwork_.pvc, vcwork_.nvc);
  FillVectorAttribute(evt, "iorgvc", vcwork_.iorgvc, vcwork_.nvc);
  FillVectorAttribute(evt, "iflgvc", vcwork_.iflgvc, vcwork_.nvc);
  FillVectorAttribute(evt, "icrnvc", vcwork_.icrnvc, vcwork_.nvc);
  FillVectorAttribute(evt, "timvc", vcwork_.timvc, vcwork_.nvc);
  FillVectorAttribute(evt, "posivc", vcwork_.posivc, vcwork_.nvc);
  FillVectorAttribute(evt, "ivtivc", vcwork_.ivtivc, vcwork_.nvc);
  FillVectorAttribute(evt, "posfvc", vcwork_.posfvc, vcwork_.nvc);
  FillVectorAttribute(evt, "ivtfvc", vcwork_.ivtfvc, vcwork_.nvc);

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
  evt.add_attribute("numne", NuHepMC::AsAttribute(nework_.numne));
  evt.add_attribute("modene", NuHepMC::AsAttribute(nework_.modene));

  FillVectorAttribute(evt, "ipne", nework_.ipne, nework_.numne);
  FillVectorAttribute(evt, "pne", nework_.pne, nework_.numne);
  FillVectorAttribute(evt, "iorgne", nework_.iorgne, nework_.numne);
  FillVectorAttribute(evt, "iflgne", nework_.iflgne, nework_.numne);
  FillVectorAttribute(evt, "icrnne", nework_.icrnne, nework_.numne);

  /**/
  /*     COMMON/FSIHIST/ */
  /**/
  /*     NVERT       : # OF VERTICES*/
  /*     POSVERT(3,I): POSITION OF I-TH VERTEX*/
  /*     IFLGVERT(I) : INTERACTION TYPE OF I-TH VERTEX */
  /*                     (*10 FOR HI-NRG)*/
  /*                     (*100 for SKDETSIM non-NEFFECT interactions e.g.
   * elastic SGPIEL;*/
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
  evt.add_attribute("nvert", NuHepMC::AsAttribute(fsihist_.nvert));

  FillVectorAttribute(evt, "posvert", fsihist_.posvert, fsihist_.nvert);
  FillVectorAttribute(evt, "iflgvert", fsihist_.iflgvert, fsihist_.nvert);

  evt.add_attribute("nvcvert", NuHepMC::AsAttribute(fsihist_.nvcvert));

  FillVectorAttribute(evt, "dirvert", fsihist_.dirvert, fsihist_.nvcvert);
  FillVectorAttribute(evt, "abspvert", fsihist_.abspvert, fsihist_.nvcvert);
  FillVectorAttribute(evt, "abstpvert", fsihist_.abstpvert, fsihist_.nvcvert);
  FillVectorAttribute(evt, "ipvert", fsihist_.ipvert, fsihist_.nvcvert);
  FillVectorAttribute(evt, "iverti", fsihist_.iverti, fsihist_.nvcvert);
  FillVectorAttribute(evt, "ivertf", fsihist_.ivertf, fsihist_.nvcvert);

  /**/
  /*       NFnvert        : number of vertices*/
  /**/
  /*       NFiflag(i)     : 4-digit flag for interaction type at i-th vertex, in
   * the form "BNTP":*/
  /*                         N: charge nucleon propagated through nucleus (0 =
   * neutron, 1 = proton)*/
  /*                         T: charge "target" nucleon the interaction is
   * taking place on*/
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
  /*       NFptot(k)      : total probability at k-th step (for testing only,
   * will be removed)*/
  /**/
  /*       Remarks:*/
  /*        - a "vertex" is actually better described as a start, end or
   * scattering point of a track*/
  /*        - at each scattering point, the first nucleon will be followed in
   * the same track, while the*/
  /*          second one will create a new track*/
  /*        - each track consists of a series of consecutive vertices. The first
   * vertex has P=0, the*/
  /*          last P=4. In between may be any number (including 0) vertices
   * where an actual scattering*/
  /*          took place (P=1,2,3).*/
  /*        - it is not possible (and not needed) to connect the second track of
   * a scattering vertex*/
  /*          with the original one. Note that "first" and "second" is purely
   * arbitrary. For nucleon*/
  /*          FSI uncertainties, only the probabilities of the scattering
   * processes have to be*/
  /*          calculated, so it is not important to know which tracks belong to
   * each other.*/
  /**/

  evt.add_attribute("nfnvert", NuHepMC::AsAttribute(nucleonfsihist_.nfnvert));
  FillVectorAttribute(evt, "nfiflag", nucleonfsihist_.nfiflag,
                      nucleonfsihist_.nfnvert);

  FillVectorAttribute(evt, "nfx", nucleonfsihist_.nfx, nucleonfsihist_.nfnvert);
  FillVectorAttribute(evt, "nfy", nucleonfsihist_.nfy, nucleonfsihist_.nfnvert);
  FillVectorAttribute(evt, "nfz", nucleonfsihist_.nfz, nucleonfsihist_.nfnvert);

  FillVectorAttribute(evt, "nfpz", nucleonfsihist_.nfpz,
                      nucleonfsihist_.nfnvert);
  FillVectorAttribute(evt, "nfpy", nucleonfsihist_.nfpy,
                      nucleonfsihist_.nfnvert);
  FillVectorAttribute(evt, "nfpz", nucleonfsihist_.nfpz,
                      nucleonfsihist_.nfnvert);
  FillVectorAttribute(evt, "nfe", nucleonfsihist_.nfe, nucleonfsihist_.nfnvert);
  FillVectorAttribute(evt, "nffirststep", nucleonfsihist_.nffirststep,
                      nucleonfsihist_.nfnvert);

  evt.add_attribute("nfnstep", NuHepMC::AsAttribute(nucleonfsihist_.nfnstep));

  FillVectorAttribute(evt, "nfecms2", nucleonfsihist_.nfecms2,
                      nucleonfsihist_.nfnstep);
  FillVectorAttribute(evt, "nfptot", nucleonfsihist_.nfptot,
                      nucleonfsihist_.nfnstep);

  // HepMC3::Print::content(evt);

  HepMCWriters[*fid]->write_event(evt);
}
