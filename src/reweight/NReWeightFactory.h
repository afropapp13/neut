#pragma once

#include "NReWeight.h"

#include "CascadePion.h"
#include "NuXSecCCQE.h"
#include "NuXSecRES.h"

#include <memory>

namespace neut {
namespace rew {
inline std::unique_ptr<NReWeight> MakeNReWeightInstance() {
  std::unique_ptr<NReWeight> nrw(new NReWeight());

  nrw->AdoptWeightEngine("CascadePion", std::unique_ptr<NReWeightEngineI>(
                                            new CascadePionEngine()));
  nrw->AdoptWeightEngine(
      "NuXSecCCQE", std::unique_ptr<NReWeightEngineI>(new NuXSecCCQEEngine()));
  nrw->AdoptWeightEngine(
      "NuXSecRES", std::unique_ptr<NReWeightEngineI>(new NuXSecRESEngine()));

  return nrw;
}
} // namespace rew
} // namespace neut
