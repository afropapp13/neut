#ifndef __Nucleus__
#define __Nucleus__ 

#include <iostream>     // std::cout
#include <algorithm>    // std::max
#include <math.h> 

#include "f2cdeclarations.h"


class Nucleus{

 private :

  int DAA; 
  int DZZ;  
  int KLAVE;
  double DMA; 
  double DNCXP; 
  double DNCA0P; 
  double DNCXN; 
  double DNCA0N; 
  double DXP;
  double DA0P; 
  double DXN; 
  double DA0N;
  double DROP; 
  double DRON; 
  double DRO; 
  double DRO0; 
  double DRO0P;
  double DRO0N;
  double QP; 
  double QValue; 
  double qp_Neut;
  double qvalue_Neut;
  double qn_ANeut;
  double qvalue_ANeut;
  double vcd_Neut[2000];
  double vcd_ANeut[2000];
 
  int ILIN; 
  int NXRO;                     
  int NCXRO; 
  int IPOL; 
  int MN; 
  int MNR; 
  int NR;
  int NRP;
  double RMAX; 
  double DR[2000];
  double DRP[2000];
  double VCD[2000];
  double FermiRFG; 
  

 public :
  
  Nucleus(int type ); 
  
  void CopyNucleiStructure(int ieta); 

  double GetQvalueGeV(int id ){

    if( id < 0 )
      return qvalue_ANeut/1000.;
    else
      return qvalue_Neut/1000.;
    
  }

  double Density(double R);

  double FindRMaximumDensity(void);

  double GetMaximumFermi(int id,double R);

  void NoBindEnergy(void) { qvalue_.qvalue = 0; }

  double GetMaximumRadious(void) {return RMAX;}

  double VC(double R );
  
  double GetFermiRFG(void) { return FermiRFG; }

  void DumpInfo(void);
}; 

#endif 
