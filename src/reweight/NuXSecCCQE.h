#pragma once

#include "NReWeightEngineI.h"

namespace neut {
namespace rew {

class NuXSecCCQEEngine : public NReWeightEngineI {
public:
  NuXSecCCQEEngine();
  ~NuXSecCCQEEngine();

  void Reconfigure();
  double CalcWeight();

  std::string EngineName() { return "NuXSecCCQEEngine"; }

private:
  void Init();
  double CalcWeightMa();

  NSyst_t MaCCQE;
  NSyst_t AxlFFCCQE;
  NSyst_t SCCVecQE;
  NSyst_t SCCAxlQE;
  NSyst_t PsFF;
  NSyst_t FAxlCCQEAlpha;
  NSyst_t FAxlCCQEGamma;
  NSyst_t FAxlCCQEBeta;
  NSyst_t FAxlCCQETheta;
  // NSyst_t FAZExp_NTerms;
  // NSyst_t FAZExp_TCut;
  // NSyst_t FAZExp_T0;
  // NSyst_t FAZExp_Q4Cut;
  NSyst_t FAZExp_A1;
  NSyst_t FAZExp_A2;
  NSyst_t FAZExp_A3;
  NSyst_t FAZExp_A4;
};

} // namespace rew
} // namespace neut
