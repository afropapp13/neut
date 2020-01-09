#include "NeutRootHandlers.h"

#include "neutmodelC.h"
#include "neworkC.h"
#include "vcvrtxC.h"
#include "vcworkC.h"

#include "posinnucC.h"

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/WriterRootTree.h"

#include <cstdint>
#include <iostream>
#include <map>
#include <memory>

std::map<int, std::unique_ptr<HepMC3::WriterRootTree> > HepMCWriters;
std::map<int, size_t> EventCounts;
int id = 0;

extern "C" {
int hepmcopen_(char *filename, int filename_len);
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

int hepmcopen_(char *filename, int filename_len) {

  nulltermstr(filename, filename_len);

  EventCounts[id] = 0;
  std::cout << "[INFO]: Opened file: " << filename << ", at ID: " << id
            << std::endl;

  std::shared_ptr<HepMC3::GenRunInfo> gri =
      std::make_shared<HepMC3::GenRunInfo>();
  gri->tools().push_back({"NEUT", "NEUTCLASSVER", ""});
  gri->add_attribute("QEModel",
                     std::make_shared<HepMC3::IntAttribute>(nemdls_.mdlqe));
  gri->add_attribute("QEMA",
                     std::make_shared<HepMC3::DoubleAttribute>(nemdls_.xmaqe));

  HepMCWriters[id] = std::unique_ptr<HepMC3::WriterRootTree>(
      new HepMC3::WriterRootTree(filename, gri));

  return id++;
}

void hepmcclose_(int *fid) {
  HepMCWriters[*fid]->close();
  HepMCWriters.erase(*fid);
}

void hepmcfill_(int *fid) {

  HepMC3::GenEvent evt(HepMC3::Units::MEV, HepMC3::Units::MM);
  evt.set_event_number(EventCounts[*fid]++);

  evt.add_attribute("Mode",
                    std::make_shared<HepMC3::IntAttribute>(nework_.modene));

  std::vector<HepMC3::GenVertexPtr> vertices;

  for (int vtx_it = 0; vtx_it < vcvrtx_.nvtxvc; vtx_it++) {
    vertices.push_back(std::make_shared<HepMC3::GenVertex>(
        HepMC3::FourVector(vcvrtx_.pvtxvc[vtx_it][0], vcvrtx_.pvtxvc[vtx_it][1],
                           vcvrtx_.pvtxvc[vtx_it][2], vcvrtx_.timvvc[vtx_it])));
    vertices.back()->set_status(1);
  }

  for (int p_it = 0; p_it < vcwork_.nvc; p_it++) {

    int vtxid = vcwork_.ivtivc[p_it] - 1;
    if (vertices.size() <= vtxid) {
      std::cout << "[ERROR]: Trying to add particle to non-existent vertex ("
                << vtxid << "/" << vertices.size() << ")." << std::endl;
      throw;
    }
    double p_E =
        sqrt(pow(vcwork_.amasvc[p_it], 2) + pow(vcwork_.pvc[p_it][0], 2) +
             pow(vcwork_.pvc[p_it][1], 2) + pow(vcwork_.pvc[p_it][2], 2));

    HepMC3::GenParticlePtr p1 = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(vcwork_.pvc[p_it][0], vcwork_.pvc[p_it][1],
                           vcwork_.pvc[p_it][2], p_E),
        vcwork_.ipvc[p_it], 100 * vcwork_.icrnvc[p_it] + vcwork_.iflgvc[p_it]);
    p1->set_generated_mass(vcwork_.amasvc[p_it]);

    bool is_in = vcwork_.iflgvc[p_it] == -1;
    if (is_in) {
      vertices[vtxid]->add_particle_in(p1);
    } else {
      vertices[vtxid]->add_particle_out(p1);
    }
  }

  for (auto vtx : vertices) {
    evt.add_vertex(vtx);
  }

  HepMC3::Print::listing(evt);

  HepMCWriters[*fid]->write_event(evt);
}
