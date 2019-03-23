# include "CrossSection.hh"

void CrossSection::Initialize(){

  if (ITYPE == 12) {
    mass_lep = 0.000511;
  } else if (ITYPE == 14) {
    mass_lep = 0.1057;
  } else if (ITYPE == 16) {
    mass_lep = 1.7770;
  }

  if (CCNC == 1) {
    ccnc_ratio = 2.;
  } else if (CCNC == 0) {
    ccnc_ratio = 1.;
  } else {
    ccnc_ratio = 0.;
  }

  LeptonMass = 0;  // Lepton mass effect on(1), off(0)

}
