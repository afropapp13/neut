#pragma once

#include <vector>

#include "NReWeightEngineI.h"
#include "NSyst.h"

namespace neut {
namespace rew {

class CascadePionEngine : public NReWeightEngineI {
public:
  CascadePionEngine();
  ~CascadePionEngine() {}

  void Reconfigure(void);
  double CalcWeight();

  std::string EngineName() const { return "CascadePionEngine"; }

private:
  NSyst_t FrAbs_pi;
  NSyst_t FrInelLow_pi;
  NSyst_t FrCExLow_pi;
  NSyst_t FrInelHigh_pi;
  NSyst_t FrCExHigh_pi;
  NSyst_t FrPiProd_pi;
};

} // namespace rew
} // namespace neut
