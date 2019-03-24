#ifndef __HadronTensor__
#define __HadronTensor__

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h> 

class HadronTensor{
 private:
  int inucleus; 
  int nbinq3;
  int nbinq0;
  int icode; 
  double q3max;
  double q0max;
  double q3bin;
  double q0bin;
  double *Tensor; 
  bool isPN; 
  
  double htfull(int i, int iq0, int iq3); 
  
 public:

  HadronTensor(std::string filename); 
  HadronTensor();
  ~HadronTensor() { delete Tensor;} 

  void Initialize(std::string filename); 

  double Interpolate(int i, int j, double q3,double q0); 
  
  double Interpolate(int i, double q3,double q0); 

  int GetNuclei() { return inucleus; }

  bool IsPN() { return isPN;}

  double GetMaxQ0(void) { return q0max;}
  double GetMaxQ3(void) { return q3max;}

  int GetBinsQ3(void) { return nbinq3;}
  int GetBinsQ0(void) { return nbinq0;}

  double GetQ0Step(void) { return q0bin;}
  double GetQ3Step(void) { return q3bin;}
  
};

#endif
