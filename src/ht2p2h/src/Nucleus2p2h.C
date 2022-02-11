#include "Nucleus2p2h.h"
#include <stdlib.h>
#include <math.h>

double RandomNuclei(void) {
  double a =  (double)random()/(double)RAND_MAX;
  return a; 
}

#define MAX(a,b) ( a > b ? a:b ) 

void Nucleus2p2h::InitializeNucleus(int nuclei){
  Nuclei *Nl = new Nuclei; 
  double DNCXP,DNCA0P,DNCXN,DNCA0N;

  if( N[nuclei] != (Nuclei*) 0 ) return; // Already initalized. 

  std::cout << " Initializing nucleus " << nuclei << std::endl; 

  N[nuclei] = Nl; 

  double uma = 0.9314941;
  
  if( nuclei == 12 ) {
    Nl->DZZ=6.; 
    Nl->DAA=12.; 
    DNCXP=DNCXN=1.692;
    DNCA0P=DNCA0N=1.082;
    Nl->klave = 1;
    Nl->FermiRFG = 0.225;
    // NP -> PP  12C -> 12N* -> 10B + 2 p
    Nl->BindEnergyPP = (10.012937-12)*uma+2.*0.938272; 
    // NN -> PN  12C -> 12N* -> 10C + p + n = 30.53
    Nl->BindEnergyPN = (10.01685-12)*uma+0.938272+0.9395654133;
    // NP -> NN  12C -> 12B* -> 10B + n + n
    Nl->BindEnergyNN = (10.012937-12)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if( nuclei == 16 ) {
    Nl->DZZ=8.;
    Nl->DAA=16.; 
    Nl->klave=1;
    DNCXP=DNCXN=1.833;
    DNCA0P=DNCA0N=1.544;
    Nl->FermiRFG = 0.225;

    // NP -> PP  16O -> 16F* -> 14N + 2 p
    Nl->BindEnergyPP = (14.00307-15.9949149)*uma+2.*0.938272; 
    // NN -> PN  16O -> 16F* -> 14O + p + n 
    Nl->BindEnergyPN = (14.00859-15.9949149)*uma+0.938272+0.9395654133;
    // NP -> NN  16O -> 16N* -> 14N + n + n
    Nl->BindEnergyNN = (14.00307-15.9949149)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 27 ) {
    Nl->DZZ=13.;
    Nl->DAA=27.;
    Nl->klave=0;       
    DNCXP=DNCXN=3.05;
    DNCA0P=DNCA0N=0.535;
    Nl->FermiRFG = 0.225;

   // NP -> PP  27Al -> 27Si* -> 25Mg + 2 p
    Nl->BindEnergyPP = (24.985837-26.98153863)*uma+2.*0.938272; 
    // NN -> PN  27Al -> 27Si* -> 25Al + p + n 
    Nl->BindEnergyPN = (24.9904281-26.98153863)*uma+0.938272+0.9395654133;
    // NP -> NN  27Al -> 27Mg* -> 25Mg + n + n
    Nl->BindEnergyNN = (24.985837-26.98153863)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 28 ) {
    Nl->DZZ=14.;
    Nl->DAA=28.;
    Nl->klave=0;
    DNCXP=DNCXN=2.93;
    DNCA0P=DNCA0N=0.569;
    Nl->FermiRFG = 0.225;

    // NP -> PP  28Si -> 28P* -> 26Al + 2 p
    Nl->BindEnergyPP = (25.98689-27.976926)*uma+2.*0.938272; 
    // NN -> PN  28Si -> 28P* -> 26Si + p + n 
    Nl->BindEnergyPN = (25.992330-27.976926)*uma+0.938272+0.9395654133;
    // NP -> NN  28Si -> 28P* -> 26Al + n + n
    Nl->BindEnergyNN = (25.98689-27.976926)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if( nuclei == 40 ) {
    Nl->DZZ=18.; 
    Nl->DAA=40.; 
    Nl->klave =0;    
    DNCXP=3.47;
    DNCA0P=0.569; 
    DNCXN=3.64;
    DNCA0N=0.569;
    Nl->FermiRFG = 0.225;

    // NP -> PP  40Ca -> 40Sc* -> 38K + 2 p
    Nl->BindEnergyPP = (37.96908-39.962591)*uma+2.*0.938272; 
    // NN -> PN  40Ca -> 40Sc* -> 38Ca + p + n 
    Nl->BindEnergyPN = (37.976318-39.962591)*uma+0.938272+0.9395654133;
    // NP -> NN  40Ca -> 40Sc* -> 38K + n + n
    Nl->BindEnergyNN = (37.96908-39.962591)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if( nuclei == 56 ) {
    Nl->DZZ=26.; 
    Nl->DAA=56.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  56Fe -> 56Co* -> 54Mn + 2 p
    Nl->BindEnergyPP = (53.9403589-55.938293)*uma+2.*0.938272; 
    // NN -> PN  56Fe -> 56Co* -> 54Fe + p + n 
    Nl->BindEnergyPN = (53.93961-55.938293)*uma+0.938272+0.9395654133;
    // NP -> NN  56Fe -> 56Co* -> 54Mn + n + n
    Nl->BindEnergyNN = (53.9403589-55.938293)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 112 ) {  // Copied from Lead.
    Nl->DZZ=56.; 
    Nl->DAA=112.;   
    Nl->klave =0;    
    DNCXP=6.624;
    DNCA0P=0.54; 
    DNCXN=6.89;
    DNCA0N=0.54;
    Nl->FermiRFG = 0.225;
    
    // NP -> PP  208Pb -> 208Bi* -> 206Tl + 2 p
    Nl->BindEnergyPP = (205.97611-207.9766521)*uma+2.*0.938272; 
    // NN -> PN  208Pb -> 208Bi* -> 206Pb + p + n 
    Nl->BindEnergyPN = (205.9744653-207.9766521)*uma+0.938272+0.9395654133;
    // NP -> NN  208Pb -> 208Bi* -> 206Tl + n + n
    Nl->BindEnergyNN = (205.97611-207.9766521)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 208 ) {
    Nl->DZZ=82.; 
    Nl->DAA=208.; 
    Nl->klave =0;    
    DNCXP=6.624;
    DNCA0P=0.54; 
    DNCXN=6.89;
    DNCA0N=0.54;
    Nl->FermiRFG = 0.225;

    // NP -> PP  208Pb -> 208Bi* -> 206Tl + 2 p
    Nl->BindEnergyPP = (205.97611-207.9766521)*uma+2.*0.938272; 
    // NN -> PN  208Pb -> 208Bi* -> 206Pb + p + n 
    Nl->BindEnergyPN = (205.9744653-207.9766521)*uma+0.938272+0.9395654133;
    // NP -> NN  208Pb -> 208Bi* -> 206Tl + n + n
    Nl->BindEnergyNN = (205.97611-207.9766521)*uma+2.*0.9395654133;
    
    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else {
    std::cout << " Nuclei " << nuclei << " is not implemented " << std::endl;
    exit(0); 
  }

 if( Nl->klave == 0 ) {
    Nl->DXP=DNCXP+3.45*DNCXP/(15.*DNCXP*DNCXP+7.*pi2*DNCA0P*DNCA0P);
    Nl->DA0P=sqrt((DNCXP*DNCXP*DNCXP+pi2*Nl->DXP*DNCA0P*DNCA0P-Nl->DXP*Nl->DXP*Nl->DXP)/(Nl->DXP*pi2));
    Nl->DXN=DNCXN+3.45*DNCXN/(15.*DNCXN*DNCXN+7.*pi2*DNCA0N*DNCA0N);
    Nl->DA0N=sqrt((DNCXN*DNCXN*DNCXN+pi2*Nl->DXN*DNCA0N*DNCA0N-Nl->DXN*Nl->DXN*Nl->DXN)/(Nl->DXN*pi2));
  }
  else  {
    double DXPP = 0.46;
    // DXPP ES IGUAL A 2/3 DEL RADIO CUADRATICO MEDIO DEL PROTON EN FM2
    Nl->DXP  = sqrt(DNCXP*DNCXP-DXPP);
    Nl->DXXP = DNCA0P*(DNCXP*DNCXP)/((1.+1.5*DNCA0P)*Nl->DXP*Nl->DXP);
    Nl->DA0P = 2.*Nl->DXXP/(2.-3.*Nl->DXXP);
    Nl->DXN  = sqrt(DNCXN*DNCXN-DXPP);
    Nl->DXXN = DNCA0N*(DNCXN*DNCXN)/((1.+1.5*DNCA0N)*Nl->DXN*Nl->DXN);
    Nl->DA0N = 2.*Nl->DXXN/(2.-3.*Nl->DXXN); 
  }

 // Normalization and maximum value.

 Nl->NormDensityN = 1.; // No normalization to start with.
 Nl->NormDensityP = 1.;


 for( int isospin = -1.; isospin <= 1; isospin += 2) {  //one -1 (neutron) and one 1 (proton)

   double RMAX =  100.;

   if (Nl->klave == 0) {
     RMAX= MAX(Nl->DXP,Nl->DXN)+9.25*MAX(Nl->DA0P,Nl->DA0N);
   }
   else {
     RMAX=sqrt(20.)*MAX(Nl->DXN,Nl->DXP);
   }

   if( isproton(isospin) ) 
     Nl->RMAXP = RMAX; 
   else
     Nl->RMAXN = RMAX; 

   std::cout << " R_max = " << RMAX << std::endl; 
   
   double step = RMAX/10000; 

   double density = 0.;
   double maxdensity = 0.; 
   
   for( double r = 0.; r < RMAX; r += step ){
     double ll = Density(nuclei,isospin,r);  
     density += ll*step;
     if( ll > maxdensity ) maxdensity = ll; 
   }

   if( isneutron(isospin) ) {
     Nl->NormDensityN = (4.*pi)*density/(Nl->DAA-Nl->DZZ);
     Nl->MaxDensityN = maxdensity/Nl->NormDensityN; 
   } else {
     Nl->NormDensityP = (4.*pi)*density/Nl->DZZ;
     Nl->MaxDensityP = maxdensity/Nl->NormDensityP; 
   }
 }

 std::cout << " Density " << Nl->MaxDensityP << "  " << Nl->MaxDensityN << std::endl; 
 
 return; 
}


double Nucleus2p2h::Density(int nuclei,int isospin, double R){

  // if( N[nuclei] != (Nuclei*) 0 ) return 0.; // Already initalized. 

  Nuclei *Nl = N[nuclei];
  
  double density = 0; 
  
  if( Nl->klave == 0 ) {
    if( isproton(isospin) ) {
      density = R*R/(exp((R-Nl->DXP)/Nl->DA0P)+1.);
      density /= Nl->NormDensityP;
    }
    else {
      density = R*R/(exp((R-Nl->DXN)/Nl->DA0N)+1.);
      density /= Nl->NormDensityN;
    }
  }
  else {
    if( isproton(isospin) ) {
      density = R*R*(1.+Nl->DA0P*(R/Nl->DXP)*(R/Nl->DXP))*exp(-(R/Nl->DXP)*(R/Nl->DXP));
      density /= Nl->NormDensityP;
    }
    else {
      density = R*R*(1.+Nl->DA0N*(R/Nl->DXN)*(R/Nl->DXN))*exp(-(R/Nl->DXN)*(R/Nl->DXN));
      density /= Nl->NormDensityN;
    }
  }

  if( density < 0. ) density = 0.;

  return density; 
}


double Nucleus2p2h::GenerateR(int nuclei, int isospin ) {
  double max,rmax;
  
  if( isneutron(isospin) ) {
    max = N[nuclei]->MaxDensityN;
    rmax = N[nuclei]->RMAXN;
  }
  else {
    max = N[nuclei]->MaxDensityP;
    rmax = N[nuclei]->RMAXP;
  }

  double r = 0.;
  double y = 0.;
  
  do{
    r = RandomNuclei()*rmax;
    y = max * RandomNuclei();
  } while( Density(nuclei,isospin,r) < y  ) ;
  
  return r;
}
