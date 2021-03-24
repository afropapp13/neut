#pragma once

#include <cmath>

namespace neut {
namespace rew {

namespace NModeDefn {

inline bool isCC(int imode) {
  return ((1 <= std::abs(imode)) && (std::abs(imode) <= 26));
}

inline bool isNC(int imode) {
  return ((31 <= std::abs(imode)) && (std::abs(imode) <= 52));
}

inline bool isNCEL(int imode) {
  return ((std::abs(imode) == 51) || (std::abs(imode) == 52));
}

inline bool isCCQE(int imode) { return (std::abs(imode) == 1); }

inline bool isCCRES(int imode) {
  return (((11 <= std::abs(imode)) && (std::abs(imode) <= 13)) ||
          (std::abs(imode) == 17) || (std::abs(imode) == 22) ||
          (std::abs(imode) == 23));
}

inline bool isCC1PI(int imode) {
  return ((11 <= std::abs(imode)) && (std::abs(imode) <= 13));
}

inline bool isCC1PI0(int imode) { return (std::abs(imode) == 12); }

inline bool isCC1PIP(int imode) {
  return ((11 == std::abs(imode)) || (std::abs(imode) == 13));
}

inline bool isCCCOH(int imode) { return (std::abs(imode) == 16); }

inline bool isNCRES(int imode) {
  return (((31 <= std::abs(imode)) && (std::abs(imode) <= 34)) ||
          (std::abs(imode) == 38) || (std::abs(imode) == 39) ||
          ((42 <= std::abs(imode)) && (std::abs(imode) <= 45)));
}

inline bool isNC1PI(int imode) {
  return ((31 <= std::abs(imode)) && (std::abs(imode) <= 34));
}

inline bool isNC1PI0(int imode) {
  return ((31 <= std::abs(imode)) && (std::abs(imode) <= 32));
}

inline bool isNCCOH(int imode) { return (std::abs(imode) == 36); }

inline bool isCOH(int imode) {
  return ((std::abs(imode) == 16) || (std::abs(imode) == 36));
}

inline bool isDIS(int imode) {
  return ((std::abs(imode) == 26) || (std::abs(imode) == 46));
}

inline bool isMPI(int imode) {
  return ((std::abs(imode) == 21) || (std::abs(imode) == 41));
}

inline bool isRES(int imode) { return (isCCRES(imode) || isNCRES(imode)); }

inline bool is1PI(int imode) { return (isCC1PI(imode) || isNC1PI(imode)); }
} // namespace NModeDefn

} // namespace rew
} // namespace neut
