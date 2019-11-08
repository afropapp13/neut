#ifndef _NEUTNUCFSISTEP_H_
#define _NEUTNUCFSISTEP_H_

#include <TObject.h>
#include <TObjArray.h>
#include <TLorentzVector.h>



class NeutNucFsiStep : public TObject {

 public:
  NeutNucFsiStep();

  ~NeutNucFsiStep(){};

  // Variables 
  Int_t          fVertFlagStep;       // Vertex flag       
  Double_t       fVertFsiRhon;
  Double_t       fECMS2; // CMS energy squared
  Double_t       fProb;  // Probability
  Double_t       fStepPel;
  Double_t       fStepPsp;
  Double_t       fStepPdp;

  //Int_t          fVertFlagStep;   // Vertex flag                                                                                                                 
  //Int_t          fVertFirstStepStep;  // ID of the First step                                                                                                      
  TLorentzVector fPosStep;            // Vertex position                                                                                                            
  TLorentzVector fMomStep;            // 4 momentum
  
  
  ClassDef(NeutNucFsiStep, 1)
};

#endif // _NEUTNUCFSISTEP_H_
