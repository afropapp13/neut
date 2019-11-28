#include "NReWeightFactory.h"

#include "NReWeightCasc.h"
#include "NReWeightNuXSecCCQE.h"
#include "NReWeightNuXSecRES.h"

namespace neut {
namespace rew {
NReWeight *DefaultConfiguredNReWeightFactory() {
  NReWeight *nrw = new NReWeight();

  nrw->AdoptWghtCalc("NReWeightNuXSecCCQE", new NReWeightNuXSecCCQE());
  nrw->AdoptWghtCalc("NReWeightNuXSecRES", new NReWeightNuXSecRES());
  nrw->AdoptWghtCalc("NReWeightCasc", new NReWeightCasc());

  return nrw;
}
} // namespace rew
} // namespace neut
