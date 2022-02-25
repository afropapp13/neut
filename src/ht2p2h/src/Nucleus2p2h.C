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
  if( nuclei == 10 ) { // 10B with proper parameters
    Nl->DZZ=5.;
    Nl->DAA=10.; 
    DNCXP=DNCXN=1.71;
    DNCA0P=DNCA0N=0.837;
    Nl->klave = 1;
    Nl->FermiRFG = 0.225;
    // NP -> PP  10B -> 10C* -> 8Be + 2 p
    Nl->BindEnergyPP = (8.00530510-10.0129370)*uma+2.*0.938272; 
    // NN -> PN  10B -> 10C* -> 8B + p + n = 30.53
    Nl->BindEnergyPN = (8.0246073-10.0129370)*uma+0.938272+0.9395654133;
    // NP -> NN  10B -> 10Be* -> 8Be + n + n
    Nl->BindEnergyNN = (8.00530510-10.0129370)*uma+2.*0.9395654133;
    
    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
     Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if( nuclei == 11 ) { // 11B with proper parameters
    Nl->DZZ=5.;
    Nl->DAA=11.; 
    DNCXP=DNCXN=1.692;
    DNCA0P=DNCA0N=1.082;
    Nl->klave = 1;
    Nl->FermiRFG = 0.225;
    // NP -> PP  11B -> 11C* -> 9Be + 2 p
    Nl->BindEnergyPP = (9.01218306-11.0093054)*uma+2.*0.938272; 
    // NN -> PN  11B -> 11C* -> 9B + p + n = 30.53
    Nl->BindEnergyPN = (9.0133296-11.0093054)*uma+0.938272+0.9395654133;
     // NP -> NN  11B -> 11Be* -> 9Be + n + n
    Nl->BindEnergyNN = (9.01218306-11.0093054)*uma+2.*0.9395654133;
    
    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if( nuclei == 12 ) {
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
  else if( nuclei == 13 ) {
    Nl->DZZ=6.; 
    Nl->DAA=13.; 
    DNCXP=DNCXN=1.692;
    DNCA0P=DNCA0N=1.082;
    Nl->klave = 1;
    Nl->FermiRFG = 0.225;
    // NP -> PP  13C -> 13N* -> 11B + 2 p
    Nl->BindEnergyPP = (11.009305167-13.0033548378)*uma+2.*0.938272; 
    // NN -> PN  13C -> 13N* -> 11C + p + n = 30.53
    Nl->BindEnergyPN = (11.01143260-13.0033548378)*uma+0.938272+0.9395654133;
    // NP -> NN  13C -> 13B* -> 11B + n + n
    Nl->BindEnergyNN = (11.009305167-13.0033548378)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if( nuclei == 14 ) {
    Nl->DZZ=7.; 
    Nl->DAA=14.; 
    DNCXP=DNCXN=1.729;
    DNCA0P=DNCA0N=1.;
    Nl->klave = 1;
    Nl->FermiRFG = 0.225;
    // NP -> PP  14N -> 14O* -> 12C + 2 p
    Nl->BindEnergyPP = (12.-14.0030740048)*uma+2.*0.938272; 
    // NN -> PN  14N -> 14O* -> 12O + p + n = 30.53
    Nl->BindEnergyPN = (12.034368-14.0030740048)*uma+0.938272+0.9395654133;
    // NP -> NN  14N -> 14C* -> 12C + n + n
    Nl->BindEnergyNN = (12.-14.0030740048)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if( nuclei == 15 ) {
    Nl->DZZ=7.; 
    Nl->DAA=15.; 
    DNCXP=DNCXN=1.729;
    DNCA0P=DNCA0N=1.;
    Nl->klave = 1;
    Nl->FermiRFG = 0.225;
    // NP -> PP  15N -> 15O* -> 13C + 2 p
    Nl->BindEnergyPP = (13.00335483534-14.0030740048)*uma+2.*0.938272; 
    // NN -> PN  15N -> 15O* -> 13O + p + n = 30.53
    Nl->BindEnergyPN = (13.024815-14.0030740048)*uma+0.938272+0.9395654133;
    // NP -> NN  15N -> 15C* -> 13C + n + n
    Nl->BindEnergyNN = (13.00335483534-14.0030740048)*uma+2.*0.9395654133;

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
 else if( nuclei == 17 ) {
    Nl->DZZ=8.;
    Nl->DAA=17.; 
    Nl->klave=1;
    DNCXP=1.798 ;
    DNCXN=1.798;
    DNCA0P=DNCA0N=1.498;
    Nl->FermiRFG = 0.225;

    // NP -> PP  17O -> 17F* -> 15N + 2 p
    Nl->BindEnergyPP = (15.0001088983-16.9991317560)*uma+2.*0.938272; 
    // NN -> PN  17O -> 17F* -> 15O + p + n 
    Nl->BindEnergyPN = (15.0030656-16.9991317560)*uma+0.938272+0.9395654133;
    // NP -> NN  17O -> 17N* -> 15N + n + n
    Nl->BindEnergyNN = (15.0001088983-16.9991317560)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 18 ) {
    Nl->DZZ=8.;
    Nl->DAA=18.; 
    Nl->klave=1;
    DNCXP=1.881 ;
    DNCXN=1.73;
    DNCA0P=DNCA0N=1.544;
    Nl->FermiRFG = 0.225;

    // NP -> PP  18O -> 18F* -> 16N + 2 p
    Nl->BindEnergyPP = (16.0061019-17.9991596121)*uma+2.*0.938272; 
    // NN -> PN  18O -> 18F* -> 16O + p + n 
    Nl->BindEnergyPN = (15.9949149-17.9991596121)*uma+0.938272+0.9395654133;
    // NP -> NN  18O -> 18N* -> 16N + n + n
    Nl->BindEnergyNN = (16.0061019-17.9991596121)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 19 ) {
    Nl->DZZ=9.;
    Nl->DAA=19.; 
    Nl->klave=0;
    DNCXP=2.59 ;
    DNCXN=2.25;
    DNCA0P=DNCA0N=0.564;
    Nl->FermiRFG = 0.225;

    // NP -> PP  19F -> 19Ne* -> 17O + 2 p
    Nl->BindEnergyPP = (16.9991317560-18.9984032)*uma+2.*0.938272; 
    // NN -> PN  19F -> 19Ne* -> 17F + p + n 
    Nl->BindEnergyPN = (17.00209524-18.9984032)*uma+0.938272+0.9395654133;
    // NP -> NN  19F -> 190* -> 17O + n + n
    Nl->BindEnergyNN = (16.9991317560-18.9984032)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 20 ) {
    Nl->DZZ=10.;
    Nl->DAA=20.; 
    Nl->klave=0;
    DNCXP=2.805;
    DNCXN=2.805;
    DNCA0P=DNCA0N=0.571;
    Nl->FermiRFG = 0.225;

    // NP -> PP  20Ne -> 20Na* -> 18F + 2 p
    Nl->BindEnergyPP = (18.0009373-19.9924401762)*uma+2.*0.938272; 
    // NN -> PN  20Ne -> 20Na* -> 18Ne + p + n 
    Nl->BindEnergyPN = (18.0057087-19.9924401762)*uma+0.938272+0.9395654133;
    // NP -> NN  20Ne -> 20F* -> 18F + n + n
    Nl->BindEnergyNN = (18.0009373-19.9924401762)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 23 ) {
    Nl->DZZ=11.;
    Nl->DAA=23.; 
    Nl->klave=0;
    DNCXP=2.773;
    DNCXN=2.776;
    DNCA0P= 0.540;
    DNCA0N=0.549;
    Nl->FermiRFG = 0.225;

    // NP -> PP  23Na -> 23Mg* -> 21Ne + 2 p
    Nl->BindEnergyPP = (20.99384669-22.9897692820)*uma+2.*0.938272; 
    // NN -> PN  23Na -> 23Mg* -> 21Na + p + n 
    Nl->BindEnergyPN = (20.99765470-22.9897692820)*uma+0.938272+0.9395654133;
    // NP -> NN  23Na -> 23Ne* -> 21Ne + n + n
    Nl->BindEnergyNN = (20.99384669-22.9897692820)*uma+2.*0.9395654133;

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
    DNCXP=DNCXN=3.14;
    DNCA0P=DNCA0N=0.537;
    Nl->FermiRFG = 0.225;

    // NP -> PP  28Si -> 28P* -> 26Al + 2 p
    Nl->BindEnergyPP = (25.98689-27.976926)*uma+2.*0.938272; 
    // NN -> PN  28Si -> 28P* -> 26Si + p + n 
    Nl->BindEnergyPN = (25.992330-27.976926)*uma+0.938272+0.9395654133;
    // NP -> NN  28Si -> 28Al* -> 26Al + n + n
    Nl->BindEnergyNN = (25.98689-27.976926)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 29 ) {
    Nl->DZZ=14.;
    Nl->DAA=30.;
    Nl->klave=0;
    DNCXP=DNCXN=3.17;
    DNCA0P=DNCA0N=0.52;
    Nl->FermiRFG = 0.225;

    // NP -> PP  29Si -> 29P* -> 27Al + 2 p
    Nl->BindEnergyPP = (26.98153863-28.976494700)*uma+2.*0.938272; 
    // NN -> PN  29Si -> 29P* -> 27Si + p + n 
    Nl->BindEnergyPN = (26.98670491-28.976494700)*uma+0.938272+0.9395654133;
    // NP -> NN  29Si -> 29Al* -> 27Al + n + n
    Nl->BindEnergyNN = (26.98153863-28.976494700)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 30 ) {
    Nl->DZZ=14.;
    Nl->DAA=30.;
    Nl->klave=0;
    DNCXP=DNCXN=3.17;
    DNCA0P=DNCA0N=0.52;
    Nl->FermiRFG = 0.225;

    // NP -> PP  30Si -> 30P* -> 28Al + 2 p
    Nl->BindEnergyPP = (27.98191009-29.973770137)*uma+2.*0.938272; 
    // NN -> PN  30Si -> 30P* -> 28Si + p + n 
    Nl->BindEnergyPN = (27.976926-29.973770137)*uma+0.938272+0.9395654133;
    // NP -> NN  30Si -> 30Al* -> 28Al + n + n
    Nl->BindEnergyNN = (27.98191009-29.973770137)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 35 ) {
    Nl->DZZ=17.;
    Nl->DAA=35.;
    Nl->klave=0;
    DNCXP=3.47;
    DNCXN=3.64;
    DNCA0P=DNCA0N=0.569;
    Nl->FermiRFG = 0.225;

    // NP -> PP  35Cl -> 35Ar* -> 33S + 2 p
    Nl->BindEnergyPP = (32.9714589099-34.968856268)*uma+2.*0.938272; 
    // NN -> PN  35Cl -> 35Ar* -> 33Cl + p + n 
    Nl->BindEnergyPN = (32.9774520-34.968856268)*uma+0.938272+0.9395654133;
    // NP -> NN  35Cl -> 35S* -> 33S + n + n
    Nl->BindEnergyNN = (32.9714589099-34.968856268)*uma+2.*0.9395654133;

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
    // NP -> NN  40Ca -> 40K* -> 38K + n + n
    Nl->BindEnergyNN = (37.96908-39.962591)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 46 ) {
    Nl->DZZ=22.; 
    Nl->DAA=46.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.04;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  46Ti -> 46V* -> 44Sc + 2 p
    Nl->BindEnergyPP = (43.9594028-45.9526316)*uma+2.*0.938272; 
    // NN -> PN  46Ti -> 46V* -> 44Ti + p + n 
    Nl->BindEnergyPN = (43.9596901-45.9526316)*uma+0.938272+0.9395654133;
    // NP -> NN  46Ti -> 46Sc* -> 44Sc + n + n
    Nl->BindEnergyNN = (43.9594028-45.9526316)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  } 
  else if( nuclei == 47 ) {
    Nl->DZZ=22.; 
    Nl->DAA=47.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  47Ti -> 47V* -> 45Sc + 2 p
    Nl->BindEnergyPP = (44.9559119-46.9517631)*uma+2.*0.938272; 
    // NN -> PN  47Ti -> 47V* -> 45Ti + p + n 
    Nl->BindEnergyPN = (44.9581256-46.9517631)*uma+0.938272+0.9395654133;
    // NP -> NN  47Ti -> 47Sc* -> 45Sc + n + n
    Nl->BindEnergyNN = (44.9559119-46.9517631)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 48 ) {
    Nl->DZZ=22.; 
    Nl->DAA=48.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  48Ti -> 48V* -> 46Sc + 2 p
    Nl->BindEnergyPP = (45.9551719-47.9479463)*uma+2.*0.938272; 
    // NN -> PN  48Ti -> 48V* -> 46Ti + p + n 
    Nl->BindEnergyPN = (45.9526316-47.9479463)*uma+0.938272+0.9395654133;
    // NP -> NN  48Ti -> 48Sc* -> 46Sc + n + n
    Nl->BindEnergyNN = (45.9551719-47.9479463)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 49 ) {
    Nl->DZZ=22.; 
    Nl->DAA=49.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  49Ti -> 49V* -> 47Sc + 2 p
    Nl->BindEnergyPP = (46.9524075-48.9478700)*uma+2.*0.938272; 
    // NN -> PN  49Ti -> 49V* -> 47Ti + p + n 
    Nl->BindEnergyPN = (46.9517631-48.9478700)*uma+0.938272+0.9395654133;
    // NP -> NN  49Ti -> 49Sc* -> 47Sc + n + n
    Nl->BindEnergyNN = (46.9524075-48.9478700)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 50 ) {
    Nl->DZZ=22.; 
    Nl->DAA=50.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  50Ti -> 50V* -> 48Sc + 2 p
    Nl->BindEnergyPP = (47.952231-49.9447912)*uma+2.*0.938272; 
    // NN -> PN  50Ti -> 50V* -> 48Ti + p + n 
    Nl->BindEnergyPN = (47.9479463-49.9447912)*uma+0.938272+0.9395654133;
    // NP -> NN  50Ti -> 50Sc* -> 48Sc + n + n
    Nl->BindEnergyNN = (47.952231-49.9447912)*uma+2.*0.9395654133;

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
  else if( nuclei == 59 ) {
    Nl->DZZ=27.; 
    Nl->DAA=59.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  59Co -> 59Ni* -> 57Fe + 2 p
    Nl->BindEnergyPP = (56.9353928-58.9331950)*uma+2.*0.938272; 
    // NN -> PN  59Co -> 59Ni* -> 57Co + p + n 
    Nl->BindEnergyPN = (56.9362914-58.9331950)*uma+0.938272+0.9395654133;
    // NP -> NN  59Co -> 59Fe* -> 57Fe + n + n
    Nl->BindEnergyNN = (56.9353928-58.9331950)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 63 ) {
    Nl->DZZ=29.; 
    Nl->DAA=63.; 
    Nl->klave =0;    
    DNCXP=4.214;
    DNCA0P=0.586; 
    DNCXN=4.214;
    DNCA0N=0.586;
    Nl->FermiRFG = 0.225;

    // NP -> PP  63Cu -> 63Zn* -> 61Ni + 2 p
    Nl->BindEnergyPP = (60.9310560-62.9295975)*uma+2.*0.938272; 
    // NN -> PN  63Cu -> 63Zn* -> 61Cu + p + n 
    Nl->BindEnergyPN = (60.9334578-62.9295975)*uma+0.938272+0.9395654133;
    // NP -> NN  63Cu -> 63Ni* -> 61Ni + n + n
    Nl->BindEnergyNN = (60.9310560-62.9295975)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 64 ) {
    Nl->DZZ=30.; 
    Nl->DAA=64.; 
    Nl->klave =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    Nl->FermiRFG = 0.225;

    // NP -> PP  64Zn -> 64Ga* -> 62Cu + 2 p
    Nl->BindEnergyPP = (61.932584-63.9291422)*uma+2.*0.938272; 
    // NN -> PN  64Zn -> 64Ga* -> 62Zn + p + n 
    Nl->BindEnergyPN = (61.934330-63.9291422)*uma+0.938272+0.9395654133;
    // NP -> NN  64Zn -> 64Cu* -> 62Cu + n + n
    Nl->BindEnergyNN = (61.932584-63.9291422)*uma+2.*0.9395654133;

    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
 else if( nuclei == 65 ) {
    Nl->DZZ=29.; 
    Nl->DAA=65.; 
    Nl->klave =0;    
    DNCXP=4.158;
    DNCA0P=0.632; 
    DNCXN=4.158;
    DNCA0N=0.632;
    Nl->FermiRFG = 0.225;

    // NP -> PP  65Cu -> 65Zn* -> 63Ni + 2 p
    Nl->BindEnergyPP = (62.9296694-64.9277895)*uma+2.*0.938272; 
    // NN -> PN  65Cu -> 65Zn* -> 63Cu + p + n 
    Nl->BindEnergyPN = (62.9295975-64.9277895)*uma+0.938272+0.9395654133;
    // NP -> NN  65Cu -> 65Ni* -> 63Ni + n + n
    Nl->BindEnergyNN = (62.9296694-64.9277895)*uma+2.*0.9395654133;

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
 else if ( nuclei == 207 ) {
    Nl->DZZ=82.; 
    Nl->DAA=207.; 
    Nl->klave =0;    
    DNCXP=6.62;
    DNCA0P=0.546; 
    DNCXN=6.62;
    DNCA0N=0.546;
    Nl->FermiRFG = 0.225;

    // NP -> PP  207Pb -> 207Bi* -> 205Tl + 2 p
    Nl->BindEnergyPP = (204.9744275-206.9758969)*uma+2.*0.938272; 
    // NN -> PN  207Pb -> 207Bi* -> 205Pb + p + n 
    Nl->BindEnergyPN = (204.9744818-206.9758969)*uma+0.938272+0.9395654133;
    // NP -> NN  207Pb -> 207Tl* -> 205Tl + n + n
    Nl->BindEnergyNN = (204.9744275-206.9758969)*uma+2.*0.9395654133;
    
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
    // NP -> NN  208Pb -> 208Tl* -> 206Tl + n + n
    Nl->BindEnergyNN = (205.97611-207.9766521)*uma+2.*0.9395654133;
    
    Nl->bindEnergyP = (Nl->BindEnergyPN+Nl->BindEnergyNN)/2.;
    Nl->bindEnergyN = (Nl->BindEnergyPP+Nl->BindEnergyPN)/2.;
  }
  else if ( nuclei == 209 ) {
    Nl->DZZ=83.; 
    Nl->DAA=209.; 
    Nl->klave =0;    
    DNCXP=6.64;
    DNCA0P=0.54; 
    DNCXN=6.87;
    DNCA0N=0.54;
    Nl->FermiRFG = 0.225;

    // NP -> PP  209Bi -> 209Po* -> 207Pb + 2 p
    Nl->BindEnergyPP = (206.9758969-208.9803987)*uma+2.*0.938272; 
    // NN -> PN  209Bi -> 209Po* -> 207Bi + p + n 
    Nl->BindEnergyPN = (206.9784707-208.9803987)*uma+0.938272+0.9395654133;
    // NP -> NN  209Bi -> 209Pb* -> 207Pb + n + n
    Nl->BindEnergyNN = (206.9758969-208.9803987)*uma+2.*0.9395654133;
    
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
