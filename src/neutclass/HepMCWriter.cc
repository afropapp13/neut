#include "NeutRootHandlers.h"

#include "neutmodelC.h"
#include "neworkC.h"
#include "vcvrtxC.h"
#include "vcworkC.h"

#include "posinnucC.h"

#include "neutmodesCPP.h"

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "HepMCNuEvtTools/WriterTools"

#include <cstdint>
#include <iostream>
#include <map>
#include <memory>

std::map<int, std::unique_ptr<HepMC3Nu::WriterRootTree> > HepMCWriters;
std::map<int, size_t> EventCounts;
int id = 0;

extern "C" {
int hepmcopen_(char *filename, int *filename_len, double *FATC);
void hepmcclose_(int *fid);
void hepmcfill_(int *fid);
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

  std::shared_ptr<HepMC3::GenRunInfo> gri = HepMC3Nu::genruninfo::GRIFactory(
      "NEUT", VERSION, "The NEUT neutrino interaction generator");

  HepMC3Nu::genruninfo::SetVertexEnumStandard(gri, "1.0");
  HepMC3Nu::genruninfo::SetParticleEnumStandard(gri, "1.0");

  std::map<int, std::string> ModeNames;

  for (size_t i = 0; i < kNneutModes; ++i) {
    if (!neutModeID[i]) { // Ignore Total
      continue;
    }
    ModeNames[neutModeID[i]] = std::string("neutrino ") + neutModeTitle[i];
    ModeNames[-neutModeID[i]] =
        std::string("anti-neutrino ") + neutModeTitle[i];
  }

  HepMC3Nu::genruninfo::SetHardScatterModeDefinitions(gri, ModeNames);

  HepMC3Nu::genruninfo::SetFluxAveragedTotalCrossSection(gri, *FATC);

  HepMCWriters[id] = std::unique_ptr<HepMC3Nu::WriterRootTree>(
      new HepMC3Nu::WriterRootTree(filename, gri));

  return id++;
}

void hepmcclose_(int *fid) {
  HepMCWriters[*fid]->close();
  HepMCWriters.erase(*fid);
}

HepMC3Nu::labels::ParticleState GetHepMCNuEvtParticleEnum(int status,
                                                          int chase) {
  if (chase == 0) {
    switch (status) {
    case -1: {
      return HepMC3Nu::labels::ParticleState::kInitialState;
    }
    case 1: {
      return HepMC3Nu::labels::ParticleState::kDecayed;
    }
    default: { return HepMC3Nu::labels::ParticleState::kIntermediate; }
    }
  } else if (chase == 1) {
    switch (status) {
    case 0:
    case 2: {
      return HepMC3Nu::labels::ParticleState::kFinalState;
    }
    default: { return HepMC3Nu::labels::ParticleState::kOther; }
    }

  } else {
    return HepMC3Nu::labels::ParticleState::kIntermediate;
  }
}

void hepmcfill_(int *fid) {

  HepMC3::GenEvent evt(HepMC3::Units::MEV, HepMC3::Units::CM);
  evt.set_event_number(EventCounts[*fid]++);

  evt.add_attribute("HardScatterMode",
                    std::make_shared<HepMC3::IntAttribute>(nework_.modene));

  HepMC3::GenVertexPtr LabFrameVtx =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector::ZERO_VECTOR());
  LabFrameVtx->set_status(
      HepMC3Nu::labels::e2i(HepMC3Nu::labels::VertexState::kLabFrame));

  HepMC3::GenVertexPtr HardScatterVtx =
      std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector::ZERO_VECTOR());
  HardScatterVtx->set_status(
      HepMC3Nu::labels::e2i(HepMC3Nu::labels::VertexState::kHardScatter));

  // std::vector<HepMC3::GenVertexPtr> vertices;

  // for (int vtx_it = 0; vtx_it < vcvrtx_.nvtxvc; vtx_it++) {
  //   vertices.push_back(std::make_shared<HepMC3::GenVertex>(
  //       HepMC3::FourVector(vcvrtx_.pvtxvc[vtx_it][0],
  //       vcvrtx_.pvtxvc[vtx_it][1],
  //                          vcvrtx_.pvtxvc[vtx_it][2],
  //                          vcvrtx_.timvvc[vtx_it])));
  //   vertices.back()->set_status(
  //       HepMC3Nu::labels::e2i(HepMC3Nu::labels::VertexState::kLabFrame));
  // }

  for (int p_it = 0; p_it < vcwork_.nvc; p_it++) {

    // Fortran to C indexing
    // int vtxid = vcwork_.ivtivc[p_it] - 1;
    // if (vertices.size() <= vtxid) {
    //   std::cout << "[ERROR]: Trying to add particle to non-existent vertex ("
    //             << vtxid << "/" << vertices.size() << ")." << std::endl;
    //   throw;
    // }
    double p_E =
        sqrt(pow(vcwork_.amasvc[p_it], 2) + pow(vcwork_.pvc[p_it][0], 2) +
             pow(vcwork_.pvc[p_it][1], 2) + pow(vcwork_.pvc[p_it][2], 2));

    HepMC3Nu::labels::ParticleState state =
        GetHepMCNuEvtParticleEnum(vcwork_.iflgvc[p_it], vcwork_.icrnvc[p_it]);

    HepMC3::GenParticlePtr p1 = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(vcwork_.pvc[p_it][0], vcwork_.pvc[p_it][1],
                           vcwork_.pvc[p_it][2], p_E),
        vcwork_.ipvc[p_it], HepMC3Nu::labels::e2i(state));
    p1->set_generated_mass(vcwork_.amasvc[p_it]);

    // std::cout << "[PART]: " << vcwork_.ipvc[p_it]
    //           << ", iflg: " << vcwork_.iflgvc[p_it]
    //           << ", icrn: " << vcwork_.icrnvc[p_it] << ", state = " << state
    //           << std::endl;

    // For now, only fill out lab frame vertex.
    switch (state) {
    case HepMC3Nu::labels::ParticleState::kInitialState: {
      // Fill out hard scatter vertex.
      if (HepMC3Nu::pid::IsLepton(p1->pid())) {
        LabFrameVtx->add_particle_in(
            std::make_shared<HepMC3::GenParticle>(p1->data()));
      }
      HardScatterVtx->add_particle_in(
          std::make_shared<HepMC3::GenParticle>(p1->data()));
      break;
    }
    case HepMC3Nu::labels::ParticleState::kFinalState: {
      // Fill out hard scatter vertex.
      if (HepMC3Nu::pid::IsLepton(p1->pid())) {
        LabFrameVtx->add_particle_out(
            std::make_shared<HepMC3::GenParticle>(p1->data()));
      }
      HardScatterVtx->add_particle_out(
          std::make_shared<HepMC3::GenParticle>(p1->data()));
      break;
    }
    default: {
      // Fill out FSI vertices.
      break;
    }
    }
  }

  evt.add_vertex(LabFrameVtx);
  evt.add_vertex(HardScatterVtx);

  // HepMC3::Print::listing(evt);

  HepMCWriters[*fid]->write_event(evt);
}
