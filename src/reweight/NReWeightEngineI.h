#pragma once

#include "NSyst.h"

#include <iostream>
#include <map>
#include <utility>

namespace neut {
namespace rew {

class NReWeightEngineI {

public:
#define ABORT_IF_UNHANDLED(syst)                                               \
  if (!DialIsHandled(syst)) {                                                  \
    std::cout << syst << " in engine: " << EngineName() << std::endl;          \
    abort();                                                                     \
  }

  virtual ~NReWeightEngineI(){};

  //! does the current engine handle the input nuisance param?
  bool DialIsHandled(NSyst_t syst) const { return fDials.count(syst); }

  double GetDial_From_Value(NSyst_t syst) const {
    ABORT_IF_UNHANDLED(syst);
    return fDials.at(syst).FromValue;
  }

  // Note that this only gets updated after calling Reconfigure
  double GetDial_To_Value(NSyst_t syst) const {
    ABORT_IF_UNHANDLED(syst);
    return fDials.at(syst).ToValue;
  }

  double GetDial_OneSigma(NSyst_t syst, double tweak) const {
    ABORT_IF_UNHANDLED(syst);
    return fDials.at(syst).GetOneSigma(tweak);
  }

  void SetDial_NumberOfSigmas(NSyst_t syst, double nsigmas) {
    ABORT_IF_UNHANDLED(syst);
    fDials[syst].Tweak = nsigmas;
  }

  void SetDial_To_Value(NSyst_t syst, double value) {
    ABORT_IF_UNHANDLED(syst);
    fDials[syst].Set_ToValue(value);
  }

  bool AnyTweaked() const {
    for (auto const &d : fDials) {
      if (d.second.IsTweaked()) {
        return true;
      }
    }
    return false;
  }

  //!  set all nuisance parameters to From values
  void Reset() {
    for (auto &d : fDials) {
      d.second.Tweak = 0;
      d.second.ApplyTweak();
    }
  }

  /// Any checks of the new values should be performed in here.
  virtual void Reconfigure() = 0;
  virtual double CalcWeight() = 0;

  virtual std::string EngineName() const = 0;

protected:
  void ApplyTweaks() {
    for (auto &d : fDials) {
      d.second.ApplyTweak();
    }
  }

  NSyst_t RegisterDial(std::string const &name, double FromValue,
                       double OneSigmaUp = 0, double OneSigmaDown = 0) {
    auto index_info_pair = NSyst::DeclareNewDialInstance(
        name, FromValue, OneSigmaUp, OneSigmaDown);
    int index = index_info_pair.first;
    fDials.emplace(std::move(index_info_pair));

    return index;
  }

  double CheckReturnWeight(double weight) {
    if (weight && !std::isnormal(weight)) {
      std::cout << "[WARN]: Abnormal weight being returned by " << EngineName()
                << std::endl;
    }
    return weight;
  }

  std::map<NSyst_t, NSyst::Info> fDials;

  NReWeightEngineI(){};

#undef ABORT_IF_UNHANDLED
};

} // namespace rew
} // namespace neut