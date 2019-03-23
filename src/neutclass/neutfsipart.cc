#include "neutfsipart.h"

NeutFsiPart::NeutFsiPart(){

  fPID       = 0;     
  fDir.SetXYZM(0.,0.,0.,0.);
  fMomLab    = 0;  
  fMomNuc    = 0;  
  fVertStart = -1;
  fVertEnd   = -1;  

};

ClassImp(NeutFsiPart)
