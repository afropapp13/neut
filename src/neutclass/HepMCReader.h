#pragma once

namespace HepMC3 {
class GenEvent;
class GenRunInfo;
} // namespace HepMC3
#include <memory>

namespace neut {
enum class VertexState {
  kNEUTVertex = 10,
  kPionFSIVertex = 13,
  kNucleonFSIVertex = 23
};

enum class ParticleState {
  kAbsorption = 13,
  kChargeExchange = 23,
  kBlocked = 19,
  kEMShower = 33,
  kHadronProduction = 43,
  kQEScatter = 53,
  kForwardElasticScatter = 63,
  kNucleonElastic = 73,
  kNucleonPiProd = 83,
  kNucleonDblPiProd = 93
};
} // namespace neut


void ReadGenRunInfo(std::shared_ptr<HepMC3::GenRunInfo const> gri);
void ReadNuHepMCEvent(std::shared_ptr<HepMC3::GenEvent> evt);
