#include "NSyst.h"

#include <iostream>
#include <stdexcept>

namespace neut {
namespace NSyst {

std::map<std::string, NSyst_t> Dials = {{"NullSystematic", kNullSystematic}};
std::map<NSyst_t, std::string> Dial_Documentation;

NSyst_t EnsureDialNameRegistered(std::string const &name) {
  if (Dials.count(name)) {
    return Dials[name];
  }
  NSyst_t index = Dials.size();
  Dials[name] = index;
  return index;
}

void DocumentDial(NSyst_t syst, std::string docstring) {
  Dial_Documentation[syst] = std::move(docstring);
}

std::string GetDialDocumentation(NSyst_t syst) {
  return Dial_Documentation[syst];
}

std::string DialAsString(NSyst_t syst) {
  for (auto const &d : Dials) {
    if (d.second == syst) {
      return d.first;
    }
  }

  std::cout << "Unknown dial index: " << syst << "\n";
  std::cout << "Known dials:"
            << "\n";
  for (auto const &d : Dials) {
    std::cout << "\t" << d.first << ": " << d.second << "\n";
  }
  std::cout << "!Did you instantiate an NReWeight instance before trying to "
               "access the global dial list?\n"
            << std::endl;
  throw std::invalid_argument(std::to_string(syst));
}

NSyst_t DialFromString(std::string const &name) {
  if (Dials.count(name)) {
    return Dials[name];
  }

  std::cout << "Unknown dial name: " << name << "\n";
  std::cout << "Known dials:"
            << "\n";
  for (auto const &d : Dials) {
    std::cout << "\t" << d.first << ": " << d.second << "\n";
  }
  std::cout << "!Did you instantiate an NReWeight instance before trying to "
               "access the global dial list?\n"
            << std::endl;
  throw std::invalid_argument(name);
}

std::pair<NSyst_t, Info> DeclareNewDialInstance(std::string const &name,
                                                double FromValue,
                                                double OneSigmaUp,
                                                double OneSigmaDown) {
  Info newinfo;
  newinfo.FromValue = FromValue;
  newinfo.ToValue = FromValue;
  newinfo.Tweak = 0;
  newinfo.OneSigmaUp = OneSigmaUp;
  newinfo.OneSigmaDown = OneSigmaDown;

  NSyst_t index = EnsureDialNameRegistered(name);

  return {index, newinfo};
}
} // namespace NSyst
} // namespace neut