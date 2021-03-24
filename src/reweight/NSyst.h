#pragma once

#include <cmath>
#include <map>
#include <string>

namespace neut {

using NSyst_t = int;

namespace NSyst {

const NSyst_t kNullSystematic = 0;

extern std::map<std::string, NSyst_t> Dials;

NSyst_t EnsureDialNameRegistered(std::string const &name);

std::string DialAsString(NSyst_t);
NSyst_t DialFromString(std::string const &);

struct Info {
  // This is the value of the dial in relevant units before a reweight, also
  // known as 'nominal', 'generated', or 'default'
  double FromValue;

  // This is the value of a parameter after reweighting
  double ToValue;

  // This is the dial variation, usually expressed in terms of number of
  // sigmas, except for dials with no defined uncertainty, where it is
  // expressed in an absolute variation from FromValue
  double Tweak;

  double OneSigmaUp;
  double OneSigmaDown;

  double GetOneSigma(double tweak) const {
    return (tweak >= 0) ? OneSigmaUp : OneSigmaDown;
  }

  void ApplyTweak() {
    double OneSigmaUncert = GetOneSigma(Tweak);
    if (std::isnormal(OneSigmaUncert)) {
      ToValue = FromValue + Tweak * OneSigmaUncert;
    } else {
      ToValue = FromValue + Tweak;
    }
  }

  void Set_ToValue(double value) {
    double uncert = GetOneSigma(FromValue - value);

    // For some dials the uncertainty may not be well defined, treak Tweak as
    // absolute
    if (!std::isnormal(uncert)) {
      Tweak = value - FromValue;
    } else { // Otherwise treat Tweak as nsigmas
      Tweak = (value - FromValue) / uncert;
    }
  }

  bool IsTweaked() const { return (std::fabs(Tweak) > 1E-8); }

  void Reset() {
    Tweak = 0;
    ApplyTweak();
  }
};

std::pair<NSyst_t, Info> DeclareNewDialInstance(std::string const &name,
                                                   double FromValue,
                                                   double OneSigmaUp = 0,
                                                   double OneSigmaDown = 0);

} // namespace neut

} // namespace neut
