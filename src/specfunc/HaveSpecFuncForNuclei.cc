#include "necardC.h"

extern "C" {
bool havespecfuncfornuclei_() {

  if (((neuttarget_.numbndp == 6) && (neuttarget_.numatom == 12)) ||  // C12
      ((neuttarget_.numbndp == 8) && (neuttarget_.numatom == 16)) ||  // O16
      ((neuttarget_.numbndp == 26) && (neuttarget_.numatom == 56))) { // Fe56
    return true;
  }
  return false;
}
bool haveeffspecfuncfornuclei_() {

  if (((neuttarget_.numbndp == 1) && (neuttarget_.numatom == 2)) ||   // H2
      ((neuttarget_.numbndp == 2) && (neuttarget_.numatom == 3)) ||   // He3
      ((neuttarget_.numbndp == 2) && (neuttarget_.numatom == 4)) ||   // He4
      ((neuttarget_.numbndp == 6) && (neuttarget_.numatom == 12)) ||  // C12
      ((neuttarget_.numbndp == 8) && (neuttarget_.numatom == 16)) ||  // O16
      ((neuttarget_.numbndp == 10) && (neuttarget_.numatom == 20)) || // Ne20
      ((neuttarget_.numbndp == 13) && (neuttarget_.numatom == 27)) || // Al27
      ((neuttarget_.numbndp == 18) && (neuttarget_.numatom == 40)) || // Ar40
      ((neuttarget_.numbndp == 26) && (neuttarget_.numatom == 56)) || // Fe56
      ((neuttarget_.numbndp == 29) && (neuttarget_.numatom == 63)) || // Cu63
      ((neuttarget_.numbndp == 30) && (neuttarget_.numatom == 64)) || // Zn64
      ((neuttarget_.numbndp == 82) && (neuttarget_.numatom == 208))   // Pb208
  ) {
    return true;
  }
  return false;
}
}
