#include "neutpart.h"

NeutPart::NeutPart(){

  fPID       = 0;
  fMass      = 0.0;
  fIsAlive   = false;
  fStatus    = -2;

  fP.SetXYZM(0.,0.,0.,0.);

  fPosIni.SetXYZT(0.,0.,0.,0.);
  fPosFin.SetXYZT(0.,0.,0.,0.);

};

ClassImp(NeutPart)
