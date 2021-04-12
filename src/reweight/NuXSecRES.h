#pragma once

#include "NReWeightEngineI.h"
#include "NSyst.h"

namespace neut {
namespace rew {

class NuXSecRESEngine : public NReWeightEngineI {
public:
  NuXSecRESEngine();
  ~NuXSecRESEngine() {}

  void Reconfigure(void);
  double CalcWeight();

  std::string EngineName() const { return "NuXSecRESEngine"; }

private:
  NSyst_t MaRES;
  NSyst_t MvRES;
  NSyst_t CA5RES;
  NSyst_t BgSclRES;
  NSyst_t MDLSPiEj;

  NSyst_t UseSeparateBgSclLMCPiBar;
  NSyst_t BgSclLMCPiBarRES;
};

} // namespace rew
} // namespace neut