#include "neutnucfsistep.h"

NeutNucFsiStep::NeutNucFsiStep(){

  fECMS2     = 0.;     
  fProb      = -1.;  

  fPosStep.SetXYZT(0.,0.,0.,0.);
  fMomStep.SetPxPyPzE(0.,0.,0.,0.);
  //fVertFirstStep = -1;
  //fVertFlag      = -1;
  fVertFlagStep      = -999;
  fVertFsiRhon = -999;
  fStepPel = -999;
  fStepPsp = -999;
  fStepPdp = -999;
};

ClassImp(NeutNucFsiStep)
