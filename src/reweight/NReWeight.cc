#include "NReWeight.h"
#include "NReWeightEngineI.h"
#include "NSyst.h"

#include "neworkC.h"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace neut {
namespace rew {

void NReWeight::AdoptWeightEngine(std::string const &name,
                                  std::unique_ptr<NReWeightEngineI> &&wengine) {
  if (!wengine) {
    std::cout << "[ERROR]: nullptr passed as weightengine named: " << name
              << std::endl;
    throw std::invalid_argument(name);
  }
  WeightEngines[name] = std::move(wengine);
}

NReWeightEngineI &NReWeight::WeightEngine(std::string const &name) {
  if (WeightEngines.count(name)) {
    return *WeightEngines[name];
  }
  std::cout << "[ERROR]: Unknown weight engine requested: " << name
            << std::endl;
  throw std::invalid_argument(name);
}

void NReWeight::Reconfigure() {

  for (auto &wght_engine : WeightEngines) {
    wght_engine.second->Reconfigure();
  }
}

double NReWeight::CalcWeight() {
  if (nework_.modene ==
      0) { // This event is considered invalid and should be squished.
    return 0;
  }

  double weight = 1;
  for (auto &wght_engine : WeightEngines) {
    weight *= wght_engine.second->CalcWeight();
  }
  return weight;
}

bool NReWeight::DialIsHandled(NSyst_t syst) const {
  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      return true;
    }
  }
  return false;
}

NSyst_t NReWeight::DialFromString(std::string const &name) const {
  return NSyst::DialFromString(name);
}

std::string NReWeight::DialAsString(NSyst_t syst) const {
  return NSyst::DialAsString(syst);
}

void CheckEngineUnique(
    std::map<std::string, std::unique_ptr<NReWeightEngineI> > const
        &WeightEngines,
    NSyst_t syst) {
  int nweightengines = 0;
  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      nweightengines++;
    }
  }

  if (nweightengines > 1) {
    std::cout << "[ERROR]: Multiple weight engines handle dial "
              << NSyst::DialAsString(syst) << "\n";
    for (auto &wght_engine : WeightEngines) {
      if (wght_engine.second->DialIsHandled(syst)) {
        std::cout << "\t" << wght_engine.second->EngineName() << "\n";
      }
    }
    std::cout << "Please Use NReWeight::WeightEngine to interact "
                 "directly with the weight engine.\n"
              << std::endl;
    throw std::invalid_argument(NSyst::DialAsString(syst));
  }
}

double NReWeight::GetDial_From_Value(NSyst_t syst) const {
  CheckEngineUnique(WeightEngines, syst);

  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      return wght_engine.second->GetDial_From_Value(syst);
    }
  }
  std::cout << "[ERROR]: Dial: " << syst
            << " is unhandled by any known weight engine." << std::endl;
  std::cout << syst << " = " << NSyst::DialAsString(syst) << std::endl;
  throw std::invalid_argument(NSyst::DialAsString(syst));
}

double NReWeight::GetDial_To_Value(NSyst_t syst) const {
  CheckEngineUnique(WeightEngines, syst);

  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      return wght_engine.second->GetDial_To_Value(syst);
    }
  }
  std::cout << "[ERROR]: Dial: " << syst
            << " is unhandled by any known weight engine." << std::endl;
  std::cout << syst << " = " << NSyst::DialAsString(syst) << std::endl;
  throw std::invalid_argument(NSyst::DialAsString(syst));
}

double NReWeight::GetDial_Tweak(NSyst_t syst) const {
  CheckEngineUnique(WeightEngines, syst);

  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      return wght_engine.second->GetDial_Tweak(syst);
    }
  }
  std::cout << "[ERROR]: Dial: " << syst
            << " is unhandled by any known weight engine." << std::endl;
  std::cout << syst << " = " << NSyst::DialAsString(syst) << std::endl;
  throw std::invalid_argument(NSyst::DialAsString(syst));
}

double NReWeight::GetDial_OneSigma(NSyst_t syst, double tweak) const {
  CheckEngineUnique(WeightEngines, syst);

  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      return wght_engine.second->GetDial_OneSigma(syst, tweak);
    }
  }
  std::cout << "[ERROR]: Dial: " << syst
            << " is unhandled by any known weight engine." << std::endl;
  std::cout << syst << " = " << NSyst::DialAsString(syst) << std::endl;
  throw std::invalid_argument(NSyst::DialAsString(syst));
}

void NReWeight::SetDial_NumberOfSigmas(NSyst_t syst, double nsigmas) {
  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      return wght_engine.second->SetDial_NumberOfSigmas(syst, nsigmas);
    }
  }
}

void NReWeight::SetDial_To_Value(NSyst_t syst, double value) {
  for (auto &wght_engine : WeightEngines) {
    if (wght_engine.second->DialIsHandled(syst)) {
      return wght_engine.second->SetDial_To_Value(syst, value);
    }
  }
}

void NReWeight::Reset() {
  for (auto &wght_engine : WeightEngines) {
    wght_engine.second->Reset();
  }
}

} // namespace rew
} // namespace neut