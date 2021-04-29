#include "neutnucfsistep.h"

NeutNucFsiStep::NeutNucFsiStep() {

  fECMS2 = 0.;
  fProb = -1.;

#ifdef NEUT_WRITE_NUCFSI
  fVertFlagStep = -999;
  fVertFsiRhon = -999;
  fStepPel = -999;
  fStepPsp = -999;
  fStepPdp = -999;
#endif
};

ClassImp(NeutNucFsiStep)
