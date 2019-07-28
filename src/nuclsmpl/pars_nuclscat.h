#ifndef _PARS_PISCAT_H
#define _PARS_PISCAT_H

// Pion Stuff
#define nNucs 2
#ifdef MAIN
const char *nucName[nNucs] = {
  "proton\0",
  "neutron\0"
};

const char *nucTitle[nNucs] = {
  "proton\0",
  "neutron\0"
};

static const float nucMass[nNucs] = {938.272, 939.566};
#else
extern const char *nucName[nNucs];
extern const char *nucTitle[nNucs];
extern const float nucMass[nNucs];
#endif
enum nuc_enum {proton,neutron};

// Interaction channels
#define nInts 4
#ifdef MAIN
const char *intName[nInts] = {
  "rxn",    
  "elas",  
  "hadpro", 
  "normal" 
};
const char *intTitle[nInts] = {
  "Reactive",
  "Elastic",
  "Hadron Production",
  "Normal"
};
#else
extern const char *intName[nInts];
extern const char *intTitle[nInts];
#endif
enum enum_iint {rxn,elas,hadpro,tot,normal};

#endif
