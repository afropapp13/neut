#include "Nucleus.h"

Nucleus::Nucleus(int type){

  double dpi = 3.141592653589793;
  
  double protonmass = 938.27208;
  double neutronmass = 939.56542;
  double uma =  931.4940954; 
  double electronmass = 0.5109989461;
  
  if(type == 11 ) {
    DZZ=5.; 
    DAA=11.; 
    DNCXP=1.692;
    DNCA0P=1.082;
    DNCXN=1.692;
    DNCA0N=1.082;
    KLAVE=1;
    
    //    qp_Neut =0.601;
    //    qvalue_Neut =16.827+qp_Neut;
    // B11 -> B10 + proton

    qvalue_Neut  = (10.0129370-11.0093054)*uma+protonmass; 

    //    qn_ANeut=3.370;
    //    qvalue_ANeut=13.880+qn_ANeut;
    // B11 -> Be10 + neutron

     qvalue_ANeut  = (10.0135338-11.0093054)*uma+neutronmass+electronmass; 

    FermiRFG = 0.225;
  }else
  if(type == 12 ) {
    DZZ=6.; 
    DAA=12.; 
    DNCXP=1.692;
    DNCA0P=1.082;
    DNCXN=1.692;
    DNCA0N=1.082;
    KLAVE=1;
    
    //    qp_Neut =0.601;
    //    qvalue_Neut =16.827+qp_Neut;
    // C12 -> C11 + proton

    qvalue_Neut  = (11.0114336-12.)*uma+protonmass; 

    //    qn_ANeut=3.370;
    //    qvalue_ANeut=13.880+qn_ANeut;
    // C12 -> B11 + neutron

     qvalue_ANeut  = (11.0093054-12.)*uma+neutronmass+electronmass; 

    FermiRFG = 0.225;
  }else
  if(type == 13 ) {
    DZZ=6.; 
    DAA=13.; 
    DNCXP=1.692;
    DNCA0P=1.082;
    DNCXN=1.692;
    DNCA0N=1.082;
    KLAVE=1;
    
    //    qp_Neut =0.601;
    //    qvalue_Neut =16.827+qp_Neut;
    // C13 -> C12 + proton

    qvalue_Neut  = (12-13.0033548378)*uma+protonmass; 

    //    qn_ANeut=3.370;
    //    qvalue_ANeut=13.880+qn_ANeut;
    // C13 -> B12 + neutron

     qvalue_ANeut  = (12.0143521-13.0033548378)*uma+neutronmass+electronmass; 

    FermiRFG = 0.225;
  }
  else if(type == 14 ) {
    DZZ=7.; 
    DAA=14.; 
    DNCXP=1.729;
    DNCA0P=1.291;
    DNCXN=1.729;
    DNCA0N=1.291;
    KLAVE=1;
    
    //    qp_Neut =0.601;
    //    qvalue_Neut =16.827+qp_Neut;
    // N14 -> N13 + proton

    qvalue_Neut  = (13.00573861-14.0030740048)*uma+protonmass; 

    //    qn_ANeut=3.370;
    //    qvalue_ANeut=13.880+qn_ANeut;
    // N14 -> C13 + neutron

     qvalue_ANeut  = (13.0033548378-14.0030740048)*uma+neutronmass+electronmass; 

    FermiRFG = 0.225;
  }
  else if(type == 15 ) { /* nitrogen isotope : tentative */
    DZZ=7.; 
    DAA=15.; 
    DNCXP=1.729;
    DNCA0P=1.291;
    DNCXN=1.729;
    DNCA0N=1.291;
    KLAVE=1;
    
    //    qp_Neut =0.601;
    //    qvalue_Neut =16.827+qp_Neut;
    // N14 -> N13 + proton

    qvalue_Neut  = (13.00573861-14.0030740048)*uma+protonmass; 

    //    qn_ANeut=3.370;
    //    qvalue_ANeut=13.880+qn_ANeut;
    // N14 -> C13 + neutron

     qvalue_ANeut  = (13.0033548378-14.0030740048)*uma+neutronmass+electronmass; 

    FermiRFG = 0.225;
  }
  else if(type == 16 ) {
    DZZ=8.;
    DAA=16.; 
    KLAVE=1;
    DNCXP=1.833;
    DNCA0P=1.544;
    DNCXN=1.833; 
    DNCA0N=1.544; 
    
    //    qp_Neut = -0.536; 
    //    qvalue_Neut=14.906+qp_Neut;
    // O16 -> O15 + proton

    qvalue_Neut= (15.0030656-15.99491461956)*uma+protonmass ;

    //    qn_ANeut=2.489; 
    //    qvalue_ANeut =10.931+qn_ANeut; 
    // O16 -> N15 + neutron

    qvalue_ANeut= (15.0001088982-15.99491461956)*uma+neutronmass+electronmass;

    FermiRFG = 0.225;    

  }
  else if(type == 17 ) {
    DZZ=9.;
    DAA=17.; 
    KLAVE=1;
    DNCXP=1.833;
    DNCA0P=1.544;
    DNCXN=1.833; 
    DNCA0N=1.544; 
    
    //    qp_Neut = -0.536; 
    //    qvalue_Neut=14.906+qp_Neut;
    // O16 -> O15 + proton

    qvalue_Neut= (15.0030656-15.99491461956)*uma+protonmass ;

    //    qn_ANeut=2.489; 
    //    qvalue_ANeut =10.931+qn_ANeut; 
    // O16 -> N15 + neutron

    qvalue_ANeut= (15.0001088982-15.99491461956)*uma+neutronmass+electronmass;

    FermiRFG = 0.225;    

  }
  else if(type == 18 ) {
    DZZ=10.;
    DAA=18.; 
    KLAVE=1;
    DNCXP=1.833;
    DNCA0P=1.544;
    DNCXN=1.833; 
    DNCA0N=1.544; 
    
    //    qp_Neut = -0.536; 
    //    qvalue_Neut=14.906+qp_Neut;
    // O16 -> O15 + proton

    qvalue_Neut= (15.0030656-15.99491461956)*uma+protonmass ;

    //    qn_ANeut=2.489; 
    //    qvalue_ANeut =10.931+qn_ANeut; 
    // O16 -> N15 + neutron

    qvalue_ANeut= (15.0001088982-15.99491461956)*uma+neutronmass+electronmass;

    FermiRFG = 0.225;    

  }
  else if(type == 19 ) { /* not correct parameters */
    DZZ=9.;
    DAA=19.; 
    KLAVE=1;
    DNCXP=1.833;
    DNCA0P=1.544;
    DNCXN=1.833; 
    DNCA0N=1.544; 
    
    //    qp_Neut = -0.536; 
    //    qvalue_Neut=14.906+qp_Neut;
    // O16 -> O15 + proton

    qvalue_Neut= (15.0030656-15.99491461956)*uma+protonmass ;

    //    qn_ANeut=2.489; 
    //    qvalue_ANeut =10.931+qn_ANeut; 
    // O16 -> N15 + neutron

    qvalue_ANeut= (15.0001088982-15.99491461956)*uma+neutronmass+electronmass;

    FermiRFG = 0.225;    

  }
  else if(type == 23 ) { /* not correct parameters */
    DZZ=11.;
    DAA=12.; 
    KLAVE=1;
    DNCXP=1.833;
    DNCA0P=1.544;
    DNCXN=1.833; 
    DNCA0N=1.544; 
    
    //    qp_Neut = -0.536; 
    //    qvalue_Neut=14.906+qp_Neut;
    // O16 -> O15 + proton

    qvalue_Neut= (15.0030656-15.99491461956)*uma+protonmass ;

    //    qn_ANeut=2.489; 
    //    qvalue_ANeut =10.931+qn_ANeut; 
    // O16 -> N15 + neutron

    qvalue_ANeut= (15.0001088982-15.99491461956)*uma+neutronmass+electronmass;

    FermiRFG = 0.225;    

  }
  else if( type == 27 ) {
    DZZ=13.;
    DAA=27.;
    KLAVE=0;       
    DNCXP=DNCXN=3.05;
    DNCA0P=DNCA0N=0.535;
    FermiRFG = 0.225;

    //    qp_Neut = -0.536; 
    //   qvalue_Neut=14.906+qp_Neut;
    // Al27 -> Al26 + proton
    qvalue_Neut = (25.98689169-26.98153863)*uma+protonmass; 

    // qn_ANeut=2.489; 
    // qvalue_ANeut =10.931+qn_ANeut; 
    // Al27 -> Mg26 + neutron
    qvalue_ANeut = (25.982592929-26.98153863)*uma+neutronmass+electronmass; 
  }
  else if ( type == 28 ) {
    DZZ=14.;
    DAA=28.;
    KLAVE=0;
    DNCXP=DNCXN=2.93;
    DNCA0P=DNCA0N=0.569;
    FermiRFG = 0.225;

    // Si28 -> Si27 + proton
    qvalue_Neut = (26.98670491-27.9769265325)*uma+protonmass; 
    
    // Si28 -> Al27 + neutron
    qvalue_ANeut = (26.98153863-27.9769265325)*uma+neutronmass+electronmass; 
    
  }
  else if(type == 35 ) {
    DZZ=17.; 
    DAA=18.; 
    KLAVE=0; 
    DNCXP=3.47;
    DNCA0P=0.569; 
    DNCXN=3.64;
    DNCA0N=0.569; 

    FermiRFG = 0.225;

    cteropnc_.dncxp = DNCXP;
    cteropnc_.dnca0p = DNCA0P;

    //    qp_Neut=7.582;
    //    qvalue_Neut=0.994+qp_Neut;
    // Cl35 -> Cl34 + proton
    qvalue_Neut = (33.97376282-34.968856268)*uma+protonmass; 
    
    //    qn_ANeut=5.830;
    //   qvalue_ANeut=7.991+qn_ANeut;
    // Cl35 -> S34 + neutron
    qvalue_ANeut = (33.96786690-34.968856268)*uma+neutronmass+electronmass; 

  }
  else if(type == 40 ) {
    DZZ=18.; 
    DAA=40.; 
    KLAVE=0; 
    DNCXP=3.47;
    DNCA0P=0.569; 
    DNCXN=3.64;
    DNCA0N=0.569; 

    FermiRFG = 0.225;

    cteropnc_.dncxp = DNCXP;
    cteropnc_.dnca0p = DNCA0P;

    //    qp_Neut=7.582;
    //    qvalue_Neut=0.994+qp_Neut;
    // Ar40 -> Ar39 + proton
    qvalue_Neut = (38.964313-39.9623831225)*uma+protonmass; 
    
    //    qn_ANeut=5.830;
    //   qvalue_ANeut=7.991+qn_ANeut;
    // Ar40 -> Cl39 + neutron
    qvalue_ANeut = (38.9680082-39.9623831225)*uma+neutronmass+electronmass; 

  }
  else if( type == 46 ) {
    DZZ=22.; 
    DAA=24.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;
	
    // Ti48 -> Ti47 + proton
    qvalue_Neut = (46.9517631-47.9479463)*uma+protonmass; 
	
    // Ti48 -> Sc47 + neutron
    qvalue_ANeut = (46.9524075-47.9479463)*uma+neutronmass+electronmass; 
  }
  else if( type == 47 ) {
    DZZ=22.; 
    DAA=25.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;
	
    // Ti48 -> Ti47 + proton
    qvalue_Neut = (46.9517631-47.9479463)*uma+protonmass; 
	
    // Ti48 -> Sc47 + neutron
    qvalue_ANeut = (46.9524075-47.9479463)*uma+neutronmass+electronmass; 
  }
  else if( type == 48 ) {
    DZZ=22.; 
    DAA=26.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;
	
    // Ti48 -> Ti47 + proton
    qvalue_Neut = (46.9517631-47.9479463)*uma+protonmass; 
	
    // Ti48 -> Sc47 + neutron
    qvalue_ANeut = (46.9524075-47.9479463)*uma+neutronmass+electronmass; 
  }
  else if( type == 49 ) {
    DZZ=22.; 
    DAA=27.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;
	
    // Ti48 -> Ti47 + proton
    qvalue_Neut = (46.9517631-47.9479463)*uma+protonmass; 
	
    // Ti48 -> Sc47 + neutron
    qvalue_ANeut = (46.9524075-47.9479463)*uma+neutronmass+electronmass; 
  }
  else if( type == 50 ) {
    DZZ=22.; 
    DAA=28.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;
	
    // Ti48 -> Ti47 + proton
    qvalue_Neut = (46.9517631-47.9479463)*uma+protonmass; 
	
    // Ti48 -> Sc47 + neutron
    qvalue_ANeut = (46.9524075-47.9479463)*uma+neutronmass+electronmass; 
  }
  else if( type == 56 ) {
    DZZ=26.; 
    DAA=56.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;
	
    // Fe56 -> Fe55 + proton
    qvalue_Neut = (54.9382934-55.9349363)*uma+protonmass; 
	
    // Fe56 -> Mn55 + neutron
    qvalue_ANeut = (54.9380451-55.9349363)*uma+neutronmass+electronmass; 
  }
  else if( type == 59 ) {
	DZZ=27.; 
    DAA=59.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;

    // Co59 -> Co58 + proton
    qvalue_Neut = (57.9357528-58.9331950)*uma+protonmass; 
    
    // Co59 -> Fe58 + neutron
    qvalue_ANeut = (57.9332744-58.9331950)*uma+neutronmass+electronmass; 
  }
  else if( type == 63 ) {
	DZZ=29.; 
    DAA=63.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;

    // Cu63 -> Cu62 + proton
    qvalue_Neut = (61.9334578-62.9295975)*uma+protonmass; 
    
    // Cu63 -> Ni62 + neutron
    qvalue_ANeut = (61.9283451-62.9295975)*uma+neutronmass+electronmass; 
  }
  else if( type == 64 ) {
	DZZ=30.; 
    DAA=64.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;

    // Zn64 -> Zn63 + proton
    qvalue_Neut = (62.9332116-63.9291422)*uma+protonmass; 
    
    // Zn64 -> Cu63 + neutron
    qvalue_ANeut = (62.9295975-63.9291422)*uma+neutronmass+electronmass; 
  }
  else if( type == 65 ) {
	DZZ=29.; 
    DAA=65.; 
    KLAVE =0;    
    DNCXP=3.971;
    DNCA0P=0.5935; 
    DNCXN=4.05;
    DNCA0N=0.5935;
    FermiRFG = 0.225;

    // Cu63 -> Cu62 + proton
    qvalue_Neut = (61.9334578-62.9295975)*uma+protonmass; 
    
    // Cu63 -> Ni62 + neutron
    qvalue_ANeut = (61.9283451-62.9295975)*uma+neutronmass+electronmass; 
  }
  else if ( type == 207 ) { /* copy of Pb208 */
    DZZ=82.; 
    DAA=207.; 
    KLAVE =0;    
    DNCXP=6.624;
    DNCA0P=0.54; 
    DNCXN=6.89;
    DNCA0N=0.54;
    FermiRFG = 0.225;

    // Pb208 -> Pb207 + proton
    qvalue_Neut = (206.9758969-207.9766521)*uma+protonmass; 
    
    // Pb208 -> Tl207 + neutron
    qvalue_ANeut = (206.977419-207.9766521)*uma+neutronmass+electronmass; 

  }
  else if ( type == 208 ) {
    DZZ=82.; 
    DAA=208.; 
    KLAVE =0;    
    DNCXP=6.624;
    DNCA0P=0.54; 
    DNCXN=6.89;
    DNCA0N=0.54;
    FermiRFG = 0.225;

    // Pb208 -> Pb207 + proton
    qvalue_Neut = (206.9758969-207.9766521)*uma+protonmass; 
    
    // Pb208 -> Tl207 + neutron
    qvalue_ANeut = (206.977419-207.9766521)*uma+neutronmass+electronmass; 

  }
  else {
    std::cout << " Nuclei type " << type << " is not implemented " << std::endl;
    exit(0); 
  }
  
  double hbarc = datos_.hbarc; 
  double xuma = datos2_.xuma;
  
  DMA = DAA*xuma; 
  
  MN=2;
  MNR=10;
  
  ILIN=2; 
  
  NXRO=0;                     
  NCXRO=0;                   
  
  int lc = 0; 
  
  convolucion_(&DNCXP,&DNCA0P,&DNCXN,&DNCA0N,&KLAVE,&lc,&DXP,&DA0P,&DXN,&DA0N);
      
  double double0 = 0.; 
  
  nxro_.nxro = NXRO; 
  ncxro_.ncxro = NCXRO; 
  nuc_.dzz = DZZ;  
  nuc_.daa = DAA;  
  roferos_.klave = KLAVE; 
  cteropnc_.dncxp = DNCXP;
  cteropnc_.dnca0p = DNCA0P;
  cteropn_.dxp = DXP; 
  cteropn_.dxn = DXN;
  cteropn_.da0p = DA0P; 
  cteropn_.da0n = DA0N; 
  cteropn_.da0n = DA0N; 

  xro_(&double0);

  DRO0P =  cteropn_.dro0n;
  DRO0N =  cteropn_.dro0n;

  
  if( KLAVE == 0 ) 
    RMAX = std::max(DXP,DXN) + 9.25* std::max(DA0P,DA0N); 
  else if ( KLAVE == 1 ) 
    RMAX=sqrt(20.)*std::max(DXN,DXP); 
  
  double df0[2000];

  double ipol = enregistre_.ipol; 

  if (ipol == 1) {
    //
    //    Finite Coulomb Potencial calculation 
    //       
    
    dsg20r_(&double0,&RMAX,&MNR,DR,&NR);
    
    double alpha = 1./137.036;
    
    for(int ir=0; ir < NR ; ir++ ) {
 
      int int5 = 5; 
      double r = DR[ir];
      
      dsg20r_(&double0,&r,&int5,DRP,&NRP);
      
      for(int irp = 0; irp< NRP; irp++ ) {
	double rp=DRP[irp];  
	df0[irp] = (rp*rp)*densq_(&rp)/r; 
      }                  
            
      double f1; 

      drg20r_(&double0,&r,&int5,df0,&f1);
      
      dsg20r_(&r,&RMAX,&int5,DRP,&NRP);
      
      for(int irp =  0; irp < NRP; irp++ ) {
	double rp = DRP[irp];
	df0[irp] = rp*densq_(&rp);
      }                  
     
      double f2;       

      drg20r_(&r,&RMAX,&int5,df0,&f2);
    
      VCD[ir]=-alpha*4.*dpi*(f1+f2);
      
    }
  }

  qvalue_.qvalue = qvalue_Neut; 

  DROP = vfpi_.drop; 
  DRON = vfpi_.dron; 
  DRO = densidad_.dro; 
  DRO0 = densidad_.dro0; 
  NXRO = nxro_.nxro; 
  NCXRO = ncxro_.ncxro;

  flags_.isinitialized = true; 

  if( type == 12 ) 
    DumpInfo();

  return; 
}

void Nucleus::CopyNucleiStructure(int ieta) {
  
  // Change potential according to neut or anti-neut 
  
  if( ieta == 1) 
    qvalue_.qvalue = qvalue_Neut; 
  else if ( ieta == -1 ) 
     qvalue_.qvalue = qvalue_ANeut; 
  else
    std::cout << " Ieta can be only 1. or -1. got: " << ieta << std::endl; 
  
  nuc_.dzz = DZZ;  
  nuc_.daa = DAA;  
  roferos_.klave = KLAVE; 
  cteropnc_.dncxp = DNCXP;
  cteropnc_.dnca0p = DNCA0P;

  cteropn_.dxp = DXP; 
  cteropn_.dxn = DXN;
  cteropn_.da0p = DA0P; 
  cteropn_.da0n = DA0N; 
  cteropn_.da0n = DA0N; 
  cteropn_.dro0p = DRO0P;
  cteropn_.dro0n = DRO0N;

  datos_.dma = DMA;   
  rinit_.rmax = RMAX; 

  nxro_.nxro = NXRO; 
  // ncxro_.ncxro = NCXRO; 

  densidad_.dro = DRO; 
  densidad_.dro0 = DRO0; 

  verif_.mnr = MNR; 
  verif_.mn = MN; 
  verif_.ilin = ILIN; 

  testieta_.ieta = ieta; 

  vfpi_.drop = DROP; 
  vfpi_.dron = DRON; 

  vcfto_.nr = NR; 
  for(int i = 0; i < NR; i++ )  {
    vcfto_.vcd[i] = VCD[i]*(double)ieta;
    vcfto_.dr[i] = DR[i];
  }

}


double Nucleus::Density(double R){
  return densq_(&R)*R*R;
}


double Nucleus::FindRMaximumDensity(void){

  double rmid = 0.; 
  double rmin = 0.;
  double rmax = RMAX; 
  double deltar = (rmax-rmin)/1e+6; 

  if( (Density(rmin+deltar)-Density(rmin))*(Density(rmax)-Density(rmax-deltar)) > 0 ) std::cout << " Error in finding maximal density " << std::endl; 

  std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl; 
  
  do {

    rmid = (rmax+rmin)/2.;
    
    if( (Density(rmin+deltar)-Density(rmin))*(Density(rmid)-Density(rmid-deltar)) < 0 ) {
      rmax = rmid; 
    } else {
      if( (Density(rmid+deltar)-Density(rmid))*(Density(rmax)-Density(rmax-deltar)) < 0 ) {
	rmin = rmid; 
      }
    }
    
    if( rmax-rmin < 1.e-4 ) break; 
    
    std::cout << rmin << "   " << rmax << std::endl; 
    
  } while(1); 
  
  std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << std::endl; 
  
  return (rmin+rmax)/2.;
}


double Nucleus::GetMaximumFermi(int id,double R) {

  double dpi = 3.141592653589793;
  
  xro_(&R); 
  
  if (id < 0 ) {
    double DKFP=pow((3.*dpi*dpi*vfpi_.drop),1./3.);
    return DKFP*datos_.hbarc/1000.;
  }
  else {
  double DKFN=pow((3.*dpi*dpi*vfpi_.dron),1./3.);
    return DKFN*datos_.hbarc/1000.;
  }
}


double Nucleus::VC(double R ){

  int IMIN = 0;
  
  int IMAX = NR-1;
  //
  // Find radious with Newton's method and interpolate the potential that is precomputed. 
  //
                              
  if( R < DR[IMIN] ){
    IMAX = IMIN+1;
  }
  else if( R > DR[IMAX] ) {
    IMIN = IMAX-1;
  }
  else{
    while ( (IMAX-IMIN) > 1 ){
                     
      int IMID = (IMIN+IMAX)/2;
	
      if( R < DR[IMID] ) {
	IMAX = IMID;
      }
      else {
	IMIN = IMID;
      }
    }
  }

  return (VCD[IMAX]-VCD[IMIN])/(DR[IMAX]-DR[IMIN])*(R-DR[IMIN])+VCD[IMIN];
}



void Nucleus::DumpInfo(void) {

  std::cout << " --------------------------------------- " << std::endl; 
  
  std::cout << " Z " << nuc_.dzz << std::endl; 
  std::cout << " A " << nuc_.daa << std::endl; 
  std::cout << " Q value " <<  qvalue_.qvalue << std::endl; 
  std::cout << " Klave " << roferos_.klave << std::endl; 
  std::cout << " DNCXP " << cteropnc_.dncxp << std::endl; 
  std::cout << " DNCA0P " << cteropnc_.dnca0p << std::endl;
  std::cout << " DMA " <<  datos_.dma << std::endl; 
  std::cout << " RMAX " <<  rinit_.rmax << std::endl; 
  std::cout << " NXRO " << nxro_.nxro << std::endl; 
  std::cout << " NCXRO " << ncxro_.ncxro << std::endl;
  std::cout << " Densidad DRO " << densidad_.dro << std::endl;
  std::cout << " Densidad DRO0 " << densidad_.dro0 << std::endl; 
  std::cout << " MNR " << verif_.mnr << std::endl;
  std::cout << " MN " <<  verif_.mn << std::endl; 
  std::cout << " ILIN " << verif_.ilin << std::endl;
  std::cout << " TEST IETA " << testieta_.ieta << std::endl;
  std::cout << " DROP " << vfpi_.drop << std::endl; 
  std::cout << " DRON " << vfpi_.dron << std::endl; 
  std::cout << " NR " << vcfto_.nr << std::endl; 
  
  for(int i = 0; i < vcfto_.nr; i++ )  {
    double r = vcfto_.dr[i];
    std::cout << "  [ "<< i << " ] = " << vcfto_.vcd[i] << "  " <<  vcfto_.dr[i] << " Pfermi =  " << GetMaximumFermi(-1,r) << "  " <<  GetMaximumFermi(1,r)  <<  std::endl; 
  }

}
