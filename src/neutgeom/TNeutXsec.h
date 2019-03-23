//////////////////////////////////////////////////////////////////////////////
// TNeutXsec.h
//
// Calls the neutxs function to determine the total cross section
// Control for obtaining cross sections for exclusive channels is through
// the card file (neut.card). If exclusive channels are specified there,
// getTotalXsec will return the cross section for the exclusive process
//
// The variable fExponentCm2 stores the log_10 of the units of the
// returned cross section in cm2, i.e. -38 means that the cros sections are
// returned in units of 10^{-38} cm^2.
//
// The function expects:
// 1. neutrino species +/- 12/14/16 for nue/numu/nutau
// 2. energy in CLHEP units (MeV)
// 3. Z/A of target nucleus
//
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include "SystemOfUnits.h"

using namespace std;

extern "C" {
  void neutxs_(int* nu_in,float* e_in,
               int* z_in, int* a_in, int* mode,
               float* totxs);
}

class TNeutXsec {

public:

  TNeutXsec() { 
    fExponentCm2 = -38;
  };

  ~TNeutXsec(){ };

  void  setExponentCm2(Int_t i) { fExponentCm2 = i;    }
  Int_t getExponentCm2()        { return fExponentCm2; };
  
  float getTotalXsec(int nutype, float e, int z, int a) 
  {  
    float xsec;
    float egev = e/CLHEP::GeV;
    int   ztemp = z;
    int   atemp = a;
    int   mode  = 0;

    neutxs_(&nutype, &egev, &ztemp, &atemp, &mode, &xsec);

    return xsec;
  };

  


private:

  int   fExponentCm2;

};
