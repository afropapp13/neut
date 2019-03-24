#ifndef _N1p1h_
#define _N1p1h_

#include <iostream>
#include <map>
#include <vector> 
#include <fstream>
#include <string>
#include "Nucleus.h" 

typedef  std::pair<int,int>  Precompindx;  // It defines the precomputed table indeces.

class N1p1h{

 private:
  double Enubin;
  double Emax;
  double facconv;
  double GFermi;  
  double unit;
 
  int NRC; //0->loop on R /  1->Random R  /  2->free nucleon
  int NPC; //0->dur     / 1->dur1(integral) /  2->durt(prand)

  char dirtable[256];
  
  int debug; 

  std::map< int , double > leptonmass;
  double protonmass;
  double neutronmass;

  double idneutron;
  double idproton; 

  // Hadron final state parameters
  double HFSparam1; 
  double HFSparam2;
  
  std::map< int , int > daughter; 

  std::vector<int> neutrinoIdlist; 

  //  std::map< Precompindx , double > QvalueGeV;
  std::map< Precompindx , std::vector<double> > IntCrossSection;       // Contains the integrated cross-sections.
  std::map< Precompindx , std::vector<double> > MaxCrossSection;       // Contains the maximal cross-sections. 
  std::map< int, Nucleus* > Nuclei;                                    // Contains the nuclei information. 
  
  void   ComputeIntegrals(int nuclei);

  void   ComputeHintegrals(int nuclei);
  
  void   CheckNuclei(int id);

  // Debugging variables 

  int numberofiterationslastevent; 
  double tmomminkept; 
  double TlepmaxKept; 
  double TlepminKept;

  // RPA values 

  double xmag; 
  int irpa;
  bool RFG; 
  bool ApplyBindEnergy;
  double RPAfp0in; 
  double RPAfp0ex;
  double RPAf; 
  double RPAfstar; 
  double RPApilambda; 
  double RPAcr0; 
  double RPArholambda; 
  double RPAgp; 
  double RPAxmpi; 
  double RPAxmrho; 
  int RPAirel;

  //Excited Nucleus formed
  double Ebind;

  //for Hydrogen only
  double cosH;
  
  // place to store params
  std::map<std::string,float> card_params;
  std::map<std::string,bool>  card_params_chk;
#define N1P1H_CARD_PARAMS 16
  char param_names[N1P1H_CARD_PARAMS][32];

 public:

  N1p1h();
  N1p1h(bool d);
  N1p1h(std::string directory, bool d = false );
  N1p1h(std::string directory, int nuclei,bool d = false );
  
  void Initialize( bool d = false );
  void Initialize( std::string directory, bool d = false );

  void InitializeAllNuclei(void); 

  void InitializeNucleus(int id, bool CInt);

  virtual ~N1p1h() { }
  
  double IntegralCrossSection(int id,int nuclei,double Enu);
  
  double CrossSectionmax(int id,int nuclei,double Enu);
  
  void  CreateVectors(int id,int nuclei,double Enu);  

  double DoubleDifferential(int id,int nuclei,double Enu,double TMu,double xcos);   
  double FourDifferential(int id,int nuclei,double Enu,double TMu,double xcos, double &R, double &dP, double &vc); 

  //  void GenerateVectors(int id, double pnu[4], double p[4][4], int idpart[4], int parent[4]);
  void GenerateVectors(int id,int nuclei,double pnu[4],  double p[4][4], int idpart[4], int parent[4], double &R, double &xsec);

    void GenerateHVectors(int id,int nuclei,double pnu[4],  double p[4][4], int idpart[4], int parent[4], double &R, double &xsec);

  double GetFraction(int id, int nuclei, double Enu,double TLepton,double coslepton); 
  
  void SetDebug(bool a) { debug = a;}
  
  void SetRFG() { RFG = true; }
  void SetLFG() { RFG = false; }

  void SetApplyBind(){ ApplyBindEnergy = true;}
  void ResetApplyBind(){ ApplyBindEnergy = false;}

  double GetEmax(void) { return Emax;}
  

  void SelectNucleusLepton(int leptonid,int nuclei); 

  //  double FindMaximumCrossSection(int nuclei,int id, double Enu);

  double EmuMax(double E,double xcos,double ml, int id,double pfermi, double vc,double qval);

  double EmuMin(double E,double xcos,double ml, int id,double pfermi, double vc,double qval);

  // Read and Write from the Tables to speed up the starting 

  bool ReadFromTable(int nuclei); 

  bool WriteTable(int nuclei);

  void RecomputeMaximum(int nuclei);

  int GetNumberofItLastEvent(void)  { return numberofiterationslastevent; }

  double GetMinTargetMomentum(void) { return tmomminkept; }

  double GetMaxTLepton(void) { return TlepmaxKept; }

  double GetMinTLepton(void) { return TlepminKept; }

  double GetMaxEnergy(void) { return Emax; }

  bool ReadConfig(const char *dir);

  bool IsNucleiInitialized(int id) { return (bool) Nuclei[id];}
   
}; 


#endif 
