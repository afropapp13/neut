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
  NSyst_t PionFSI_AbsProb;
  NSyst_t PionFSI_QELowMomProb;
  NSyst_t PionFSI_QEHighMomProb;
  NSyst_t PionFSI_CExLowMomProb;
  NSyst_t PionFSI_CExHighMomProb;
  NSyst_t PionFSI_InelProb;

};

} // namespace rew
} // namespace neut
