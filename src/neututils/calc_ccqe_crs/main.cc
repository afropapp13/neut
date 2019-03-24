# include <stdio.h>
# include <iostream>
# include <cmath>

# include "CrossSection.hh"

using namespace std;

#include "neutmodelC.h"

struct neutmodel_common neutmodel_;
struct nemdls_common    nemdls_;

int main(int argc, char *argv[]) {

  int iEnergy;
  int MaxiEnergy, nloop;
  double Energy, xsec;

  FILE *output;
  output = fopen("cross.dat","w");

  //--- initialize
  int nu_type    = atoi(argv[1]);
  int nu_sign    = atoi(argv[2]);
  double Pfermi  = atof(argv[3]);
  double XMaxial = atof(argv[4]);

  CrossSection* QECrossSection = new CrossSection(Pfermi, nu_type, nu_sign, XMaxial);
   
  if (nu_type != 16) {
    MaxiEnergy = 210;
    nloop = 200; // precision for clenshaw-curtis's rule
  } else {
    MaxiEnergy = 125;
    nloop = 500;
  }
  //--- loop over neutrino energy
  for (iEnergy=0; iEnergy<MaxiEnergy; iEnergy++) {
    if (nu_type != 16) {
      Energy = energy_table1[iEnergy]*1000.; // in MeV
    } else {
      Energy = energy_table2[iEnergy]*1000.;
    }

    xsec = QECrossSection->TotalCrossSection(Energy, nloop);
    fprintf(output, "%5d %9.1f %8.5f\n", nu_type*nu_sign, Energy, xsec);
  }

  fclose(output);
  delete QECrossSection;

}
