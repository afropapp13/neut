#ifndef _PARS_PISCAT_H
#define _PARS_PISCAT_H

// Pion Stuff
#define nPions 3
#ifdef MAIN
const char *pionName[nPions] = {
  "pi0", 
  "piP", 
  "piM"
};

const char *pionTitle[nPions] = {
  "#pi^{0}",
  "#pi^{+}",
  "#pi^{-}"
};
static const float pionMass[nPions] = {134.9766,139.57018,139.57018};
#else
extern const char *pionName[nPions];
extern const char *pionTitle[nPions];
extern const float pionMass[nPions];
#endif
enum pion_enum {pi0,pip,pim};

// Nucleon stuff
#ifdef MAIN
const char *nucName[2] = {"n","p"};
const char *nucTitle[2] = {"Neutron", "Proton"};
const float nucMass[2] = {939.565560,938.272013};
#else
extern const char *nucName[2];
extern const char *nucTitle[2];
extern const float nucMass[2];
#endif

// Interaction channels
#define nInts 9
#ifdef MAIN
const char *intName[nInts] = {
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
const char *intTitle[nInts] = {
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
#else
extern const char *intName[nInts];
extern const char *intTitle[nInts];
#endif
enum enum_iint {rxn,absb,scatt,scx,dcx,elas,hadpro,tot,normal};

// 2-pi production
#define nPiprods 7
#ifdef MAIN
const char *piprodName[nPiprods] = {
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
#else
extern const char *piprodName[nPiprods];
extern char *piprodTitle[nPiprods];
#endif


#endif
