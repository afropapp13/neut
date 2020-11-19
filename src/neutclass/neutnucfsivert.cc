#include "neutnucfsivert.h"

NeutNucFsiVert::NeutNucFsiVert(){

  fPos.SetXYZT(0.,0.,0.,0.);
  fMom.SetPxPyPzE(0.,0.,0.,0.);
  fVertFirstStep = -1;
  fVertFlag      = -1;

};

ClassImp(NeutNucFsiVert)
