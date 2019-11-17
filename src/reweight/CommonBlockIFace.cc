
#include "CommonBlockIFace.h"

#include <iostream>
#include <sstream>

namespace neut {
namespace rew {

static CommonBlockIFace *gCommonBlockIFace = NULL;

void CommonBlockIFace::SetGenCard(std::string const &GenCardLocation) {
  // Set Card
  // Read Card
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

std::string CommonBlockIFace::ToString() const {
  std::stringstream ss("");

  ss << "NEUTCARD: \n";

  return ss.str();
}

} // namespace rew
} // namespace neut
