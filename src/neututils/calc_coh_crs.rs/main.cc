#include <iostream>
#include <cmath>
#include <cstdio>

#include "CrossSection.hh"

using namespace std;

int main(int argc, char *argv[]) {

  int iEnergy, nloop;
  const int MaxiEnergy = 370;
  double Energy, xsec;

  FILE *output;
  output = fopen("cross.dat","w");

  //--- initialize
  int nu_type = atoi(argv[1]);
  int current = atoi(argv[2]); // CC(1), NC(0)

  CrossSection* COHCrossSection = new CrossSection(nu_type, current);
   
  //--- loop over neutrino energy
  nloop = 250;
  for (iEnergy=1; iEnergy<=MaxiEnergy; iEnergy++) {
      if (iEnergy < 100) {
	  Energy = double(iEnergy)*0.01;           // 0 - 1GeV
      } else if (iEnergy >= 100 && iEnergy < 190) {
	  Energy = double(iEnergy-100)*0.1 + 1.0;  // 1 - 10GeV
      } else if (iEnergy >= 190 && iEnergy < 280) {
	  Energy = double(iEnergy-190)*1. + 10.0;  // 10 - 100GeV
      } else if (iEnergy >= 280 && iEnergy <= 370) {
	  Energy = double(iEnergy-280)*10. + 100.; // 100 - 1000GeV
      }      

      if (Energy != 0.) {
	  xsec = COHCrossSection->TotalCrossSection(Energy, nloop);
      } else {
	  xsec = 0.;
      }

    fprintf(output, "%9.5e %9.5e\n", Energy, xsec);
  }

  fclose(output);
  delete COHCrossSection;

}
