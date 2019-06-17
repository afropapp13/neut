#ifndef __Nucleus2p2h__
#define __Nucleus2p2h__ 

#include <iostream>     // std::cout
#include <algorithm>    // std::max
#include <math.h>
#include <map>
#include <math.h>

#define hbarc 197.3269602
#define pi    3.14159265359
#define pi2   9.86960440109

typedef struct{
  double DZZ; 
  double DAA; 
  double DXP;
  double DXXP;
  double DA0P;
  double DXN;
  double DXXN;
  double DA0N;
  double klave;

  double NormDensityN;
  double NormDensityP; 
  double MaxDensityN;
  double MaxDensityP;
  
  double bindEnergyP;
  double bindEnergyN;
  double BindEnergyPP;
  double BindEnergyNN;
  double BindEnergyPN;
  
  double RMAXP; 
  double RMAXN;
  
  double FermiRFG;
} Nuclei; 


typedef std::map<int, Nuclei* > mapnuclei;

class Nucleus2p2h{

 private :
  mapnuclei N;
  bool isproton(int isospin) { return (isospin>0);}
  bool isneutron(int isospin) { return (isospin<0);}
  
 public :
  
  Nucleus2p2h(){;}

  void InitializeNucleus(int nuclei); 

  double GetMaxDensity(int nuclei,int isospin) { if( isospin > 0 ) return N[nuclei]->MaxDensityP; else return N[nuclei]->MaxDensityN; }

  double GetBind(int nuclei,int itype) {  if( itype > 0 ) return N[nuclei]->bindEnergyP; else return N[nuclei]->bindEnergyN;}

  double GetBind(int nuclei,int itype, int pn ) {
    if( itype > 0 ) {  // neutrino 
      if( pn ) 
	return N[nuclei]->BindEnergyPP;
      else
	return N[nuclei]->BindEnergyPN;
    }
    else {  // antineutrino 
      if( pn ) 
	return N[nuclei]->BindEnergyNN;
      else
	return N[nuclei]->BindEnergyPN;
    }
  }
  
  double GenerateR(int nuclei,int isospin); 

  double Density(int nuclei,int isospin, double R); 

  double GetFermiRFG(int nuclei) { return N[nuclei]->FermiRFG; }

  double GetFermiLFG(double R,int nuclei,int isospin){
    double kf = pow((3.*pi2*Density(nuclei,isospin,R))/R/R,(1./3.)); //in fm
	//    if( std::isnan(kf) ) std::cout << " >>>>>>>>>>> " <<  R << " " << Density(nuclei,isospin,R) << std::endl; 
	if( isnan(kf) ) std::cout << " >>>>>>>>>>> " <<  R << " " << Density(nuclei,isospin,R) << std::endl; 
    return kf*hbarc/1000.;  // fm to GeV  
  }

  double GetRmax(int nuclei) {
    if( N[nuclei]->RMAXP > N[nuclei]->RMAXN ) return N[nuclei]->RMAXP;
    return N[nuclei]->RMAXN;
  }
  
}; 

#endif 
