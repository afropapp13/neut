#include "neutnucfsistep.h"

NeutNucFsiStep::NeutNucFsiStep() {

  fECMS2 = 0.;
  fProb = -1.;

#if defined(NEUT_WRITE_NUCFSI) or defined(NEUT_READ_NUCFSI)
  fVertFlagStep = -999;
  fVertFsiRhon = -999;
  fStepPel = -999;
  fStepPsp = -999;
  fStepPdp = -999;
#endif
};

ClassImp(NeutNucFsiStep)
