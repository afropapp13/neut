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

  NSyst_t ROP3P3;
  NSyst_t ROP1P1;
  NSyst_t ROM3M3;
  NSyst_t ROM1M1;
  NSyst_t ROP3P1;
  NSyst_t ROM1M3;
  NSyst_t ROP3M1;
  NSyst_t ROP1M3;

};

} // namespace rew
} // namespace neut
