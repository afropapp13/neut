#ifndef _PARS_PISCAT_H
#define _PARS_PISCAT_H

// Pion Stuff
const int nPions=3;
char *pionName[nPions] = {
  "pi0", 
  "piP", 
  "piM"
};

char *pionTitle[nPions] = {
  "#pi^{0}",
  "#pi^{+}",
  "#pi^{-}"
};
float pionMass[nPions] = {134.9766,139.57018,139.57018};
enum pion_enum {pi0,pip,pim};

// Nucleon stuff
char *nucName[2] = {"n","p"};
char *nucTitle[2] = {"Neutron", "Proton"};
const float nucMass[2] = {939.565560,938.272013};

// Interaction channels
const int nInts = 9;
char *intName[nInts] = {
  "rxn",    
  "absb",   
  "scatt",  
  "scx",    
  "dcx",    
  "elas",   
  "hadpro", 
  "tot",    
  "normal" 
};
enum enum_iint {rxn,absb,scatt,scx,dcx,elas,hadpro,tot,normal};

char *intTitle[nInts] = {
  "Reactive",
  "Absorption",
  "Quasi-elastic",
  "Single CX",
  "Double CX",
  "Elastic",
  "#pi Production",
  "Total",
  "Normal"
};



// 2-pi production
const int nPiprods = 7;
char *piprodName[nPiprods] = {
  "pipiTot",  
  "pi0pi0",
  "piPpiM",
  "piPpi0",
  "piPpiP",
  "piMpi0",  
  "piMpiM"  
};

char *piprodTitle[nPiprods] = {
  "#pi#piX",
  "#pi^{0}#pi^{0}X",
  "#pi^{+}#pi^{-}X",
  "#pi^{+}#pi^{0}X",
  "#pi^{+}#pi^{+}X",
  "#pi^{-}#pi^{0}X",
  "#pi^{-}#pi^{-}X"
};

#endif
