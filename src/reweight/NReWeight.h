#pragma once

#include "NSyst.h"

#include <map>
#include <memory>
#include <string>

namespace neut {
namespace rew {

class NReWeightEngineI;

class NReWeight {
public:
  NReWeight() {}
  ~NReWeight() {}

  void AdoptWeightEngine(std::string const &name,
                         std::unique_ptr<NReWeightEngineI> &&engine);

  NReWeightEngineI &WeightEngine(std::string const &name);

  void Reconfigure();

  double CalcWeight();

  bool DialIsHandled(NSyst_t) const;
  NSyst_t DialFromString(std::string const &) const;
  std::string DialAsString(NSyst_t) const;
  double GetDial_From_Value(NSyst_t) const;
  double GetDial_To_Value(NSyst_t) const;
  double GetDial_OneSigma(NSyst_t, double) const;
  void SetDial_NumberOfSigmas(NSyst_t, double);
  void SetDial_To_Value(NSyst_t, double);

  void Reset();

private:
  std::map<std::string, std::unique_ptr<NReWeightEngineI>> WeightEngines;
};

} // namespace rew
} // namespace neut
