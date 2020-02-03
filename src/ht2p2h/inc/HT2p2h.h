#ifndef _HT2p2h_
#define _HT2p2h_

#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <vector> 
#include <fstream>
#include <string>
#include "HadronTensor.h"
#include "Nucleus2p2h.h" 
#include "nieves2p2hC.h"

extern "C"
{
	float rlu_(int *);
}

typedef  std::pair<int,int>  Precompindx;  // It defines the precomputed table indeces.

class HT2p2h:Nucleus2p2h{

 private:
  double Enubin;
  double Emax;
  double facconv;
  double GFermi; 
  bool  RFG;

  int  HadKinMode; 

  double default_pn_nn_fraction; 
  
  double debug; 

  std::map< int , double > leptonmass;
  double protonmass;
  double neutronmass;

  double idneutron;
  double idproton; 

  // Hadron final state parameters
  double HFSparam1; 
  double HFSparam2;

  bool ApplyBindEnergy;
  
  std::map< int , int > daughter; 

  std::vector<int> neutrinoIdlist; 

  //  std::map< Precompindx , double > QvalueGeV;
  std::map< Precompindx , std::vector<double> > IntCrossSection;       // Contains the integrated cross-sections.
  std::map< Precompindx , std::vector<double> > MaxCrossSection;       // Contains the maximal cross-sections. 
  std::map< int , HadronTensor* > Tensor;                      // Contains the tensors.....
  std::map< int , bool > Tensor_init;                      // Initialized flag of the tensor
  std::map< int , std::string > Tensor_fname;                     // Filename of the tensor 
  void   ReadTensors(std::string dirname );
  void   ComputeIntegrals(int nuclei);
  int    ReadIntegrals(int nuclei);

  void   CheckNuclei(int id);
  void   CheckNuclei_2(int id);

 public:
  HT2p2h();
  HT2p2h(std::string dirname, bool ApplyBind, bool d=false );

  void Initialize(std::string dirname, bool bind = true, bool d = false );

  virtual ~HT2p2h() { }
  
  double IntegralCrossSection(int id,int nuclei,double Enu);
  
  void   CreateVectors(int id,int nuclei,double Enu);  

  int GenerateLeptonKinematics(int id,int nuclei,double Enu,double &TLepton,double &xcos,double &R);

  double DoubleDifferential(int id,int nuclei,double Enu,double TMu,double xcos, double R, bool pn = false);

  double GetHTElement(int nuclei,double q3,double q0);

  int GenerateHadronKinematics( int id, int nuclei, double frac, double q[4], double ni1[4], double ni2[4], double n1[4], double n2[4], int idNucleon[4], double &R );

  int GenerateVectors(int id, int nuclei, double pnu[4], double p[6][4], int idpart[6], int parent[6], double &R );

  double GetFraction(int id, int nuclei, double Enu,double TLepton,double coslepton,double R); 
  
  void SetDebug(bool a) { debug = a;}

  void SetRFG() { RFG = true; }
  void SetLFG() { RFG = false; }

  void SetApplyBind(){ ApplyBindEnergy = true;}
  void ResetApplyBind(){ ApplyBindEnergy = false;}
  
  void SetHFSparameters(double p1,double p2) { HFSparam1 = p1; HFSparam2 = p2; }

  void SetHadronKinMode(int val){  HadKinMode = val; }

  double GetQ0Fermivalue(int nuclei,int id,double R);

  double Random(void)
  {

	int dummy;
	double a = (double)(rlu_(&dummy));
	return a;

  }

}; 


#endif 


