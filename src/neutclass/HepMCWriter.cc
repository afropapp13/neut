#include "NeutRootHandlers.h"

#include "necardC.h"
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

#include "NuHepMC/WriterTools"

#include <cstdint>
#include <iostream>
#include <map>
#include <memory>

std::map<int, std::unique_ptr<NuHepMC::WriterRootTree> > HepMCWriters;
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

  NuHepMC::genruninfo::SetFluxAveragedTotalCrossSection(gri, *FATC);

  HepMCWriters[id] = std::unique_ptr<NuHepMC::WriterRootTree>(
      new NuHepMC::WriterRootTree(filename, gri));

  return id++;
}

void hepmcclose_(int *fid) {
  HepMCWriters[*fid]->close();
  HepMCWriters.erase(*fid);
}

NuHepMC::labels::ParticleState GetHepMCNuEvtParticleEnum(int status,
                                                         int chase) {
  if (chase == 0) {
    switch (status) {
    case -1: {
      return NuHepMC::labels::ParticleState::kInitialState;
    }
    case 1: {
      return NuHepMC::labels::ParticleState::kDecayed;
    }
    default: {
      return NuHepMC::labels::ParticleState::kIntermediate;
    }
    }
  } else if (chase == 1) {
    switch (status) {
    case 0:
    case 2: {
      return NuHepMC::labels::ParticleState::kFinalState;
    }
    default: {
      return NuHepMC::labels::ParticleState::kOther;
    }
    }

  } else {
    return NuHepMC::labels::ParticleState::kIntermediate;
  }
}

void hepmcfill_(int *fid) {

  HepMC3::GenEvent evt(HepMC3::Units::MEV, HepMC3::Units::CM);
  evt.set_event_number(EventCounts[*fid]++);

  NuHepMC::genevent::SetHardScatterMode(evt, nework_.modene);

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

  // std::vector<HepMC3::GenVertexPtr> vertices;

  // for (int vtx_it = 0; vtx_it < vcvrtx_.nvtxvc; vtx_it++) {
  //   vertices.push_back(std::make_shared<HepMC3::GenVertex>(
  //       HepMC3::FourVector(vcvrtx_.pvtxvc[vtx_it][0],
  //       vcvrtx_.pvtxvc[vtx_it][1],
  //                          vcvrtx_.pvtxvc[vtx_it][2],
  //                          vcvrtx_.timvvc[vtx_it])));
  //   vertices.back()->set_status(
  //       NuHepMC::labels::e2i(NuHepMC::labels::VertexState::kLabFrame));
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

    NuHepMC::labels::ParticleState state =
        GetHepMCNuEvtParticleEnum(vcwork_.iflgvc[p_it], vcwork_.icrnvc[p_it]);

    HepMC3::GenParticlePtr p1 = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(vcwork_.pvc[p_it][0], vcwork_.pvc[p_it][1],
                           vcwork_.pvc[p_it][2], p_E),
        vcwork_.ipvc[p_it], NuHepMC::labels::e2i(state));
    p1->set_generated_mass(vcwork_.amasvc[p_it]);

    // std::cout << "[PART]: " << vcwork_.ipvc[p_it]
    //           << ", iflg: " << vcwork_.iflgvc[p_it]
    //           << ", icrn: " << vcwork_.icrnvc[p_it] << ", state = " << state
    //           << std::endl;

    // For now, only fill out lab frame vertex.
    switch (state) {
    case NuHepMC::labels::ParticleState::kInitialState: {
      // Fill out hard scatter vertex.
      if (NuHepMC::pid::IsLepton(p1->pid())) {
        LabFrameVtx->add_particle_in(
            std::make_shared<HepMC3::GenParticle>(p1->data()));
      }
      HardScatterVtx->add_particle_in(
          std::make_shared<HepMC3::GenParticle>(p1->data()));
      break;
    }
    case NuHepMC::labels::ParticleState::kFinalState: {
      LabFrameVtx->add_particle_out(
          std::make_shared<HepMC3::GenParticle>(p1->data()));
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
