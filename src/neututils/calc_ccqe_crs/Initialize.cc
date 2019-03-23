#include <iostream>
# include "CrossSection.hh"

using namespace std;

void CrossSection::Initialize(){

  double const XMN = 939.6;  
  
  //PFermi =  225.; // Fermi momentum (MeV/c)
  // Ebind  =   27.; // Binding energy (MeV/c)

  Ebind = sqrt(XMN*XMN+PFermi*PFermi) - XMN;
  
//  Maxial = 1310.; // Axial mass (MeV/c)
//  Maxial = 1210.; // Axial mass (MeV/c)
// Maxial = 1110.; // Axial mass (MeV/c)
  Maxial = 1030.; // Axial mass (MeV/c)
  Nnum = 8.;      // Number of neutron
  Znum = 8.;      // Number of proton
  //ITYPE = 12;     // Neutrino type
  //NUSIG = 1;      // Neutrino(1)/Anti-neutrino(-1)
  kappa = 1.0;    // kappa factor for pauli blocking effect(see PRL,100.032301)
  FF_model = 1;   // Axial form factor (0:dipole, 1:BBBA07)
  Trans_Corr = 0; // Transverse correction ( 0: default, 1: Corrected )

  cout << "PFermi = " << PFermi << " ( Ebind = " << Ebind << ")\n";
}
