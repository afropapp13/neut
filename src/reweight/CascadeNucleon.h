#pragma once

#include <vector>

#include "NReWeightEngineI.h"
#include "NSyst.h"

namespace neut {
namespace rew {

class CascadeNucleonEngine : public NReWeightEngineI {
public:
  CascadeNucleonEngine();
  ~CascadeNucleonEngine() {}

  void Reconfigure(void);
  double CalcWeight();

  std::string EngineName() const { return "CascadeNucleonEngine"; }

private:
  NSyst_t NucleonFSI_Total;
  NSyst_t NucleonFSI_Elastic;
  NSyst_t NucleonFSI_Single;
  NSyst_t NucleonFSI_Double;

};

} // namespace rew
} // namespace neut
