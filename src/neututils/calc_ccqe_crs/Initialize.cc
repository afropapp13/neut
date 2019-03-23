#include <iostream>
# include "CrossSection.hh"

#define MAIN
#include "neutmodelC.h"

using namespace std;

void CrossSection::Initialize(double xamass){

  double const XMN = 939.6;  
  
  //PFermi =  225.; // Fermi momentum (MeV/c)
  // Ebind  =   27.; // Binding energy (MeV/c)

  Ebind = sqrt(XMN*XMN+PFermi*PFermi) - XMN;
  
// Maxial = 1310.; // Axial mass (MeV/c)
// Maxial = 1210.; // Axial mass (MeV/c)
// Maxial = 1110.; // Axial mass (MeV/c)
// Maxial = 1030.; // Axial mass (MeV/c)
  XMaxial = xamass; // Axial mass (MeV/c)
  Nnum = 8.;      // Number of neutron
  Znum = 8.;      // Number of proton
  //ITYPE = 12;     // Neutrino type
  //NUSIG = 1;      // Neutrino(1)/Anti-neutrino(-1)
  kappa = 1.0;    // kappa factor for pauli blocking effect(see PRL,100.032301)
//  FF_model = 1;   // Axial  form factor (0:dipole, 1:BBBA07)
  A_FF_model = 1;   // Axial  form factor (0:dipole, 1:BBBA07)
  V_FF_model = 2;   // Vector form factor (0:dipole, 1:BBBA05, 2: BBBA07)
//  V_FF_model = 0;   // Vector form factor (0:dipole, 1:BBBA05)
  Trans_Corr = 0; // Transverse correction ( 0: default, 1: Corrected )

  nemdls_.xmaqe = XMaxial/1.e3;

  cout << "PFermi = " << PFermi << " Ebind = " << Ebind 
	   << "Axial  = " << XMaxial<< " A_FF  = " << A_FF_model
	   << "V_FF   = " << V_FF_model << endl;
}
