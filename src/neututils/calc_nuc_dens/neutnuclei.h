#ifndef NEUTNUCLEI_H
#define NEUTNUCLEI_H

// Available NEUT nuclei (see neutcore/nesettarg.F)
const int nNuclei = 26;
int nBoundProtons[nNuclei] = {5,6,7,8,9,11,13,14,16,17,18,20,22,26,27,28,29,30,40,41,45,50,67,79,82,83};

float r_step = 0.1; // Defined by step size in eftrace.F
float max_r = 15;   // Must be greater than 2.5*eftarget_.c of largest atom

#endif
