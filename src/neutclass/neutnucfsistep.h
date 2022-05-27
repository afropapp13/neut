#ifndef _NEUTNUCFSISTEP_H_
#define _NEUTNUCFSISTEP_H_

#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TObject.h>

class NeutNucFsiStep : public TObject {

public:
  NeutNucFsiStep();

  ~NeutNucFsiStep(){};

  // Variables

  Int_t fVertFlagStep; // Vertex flag
  Double_t fVertFsiRhon;
  Double_t fECMS2; // CMS energy squared
  Double_t fProb;  // Probability
  Double_t fStepPel;
  Double_t fStepPsp;
  Double_t fStepPdp;

  ClassDef(NeutNucFsiStep, 3)
};

#endif // _NEUTNUCFSISTEP_H_
