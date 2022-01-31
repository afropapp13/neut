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

#if defined(NEUT_WRITE_NUCFSI) or defined(NEUT_READ_NUCFSI)
  Int_t fVertFlagStep; // Vertex flag
  Double_t fVertFsiRhon;
#endif
  Double_t fECMS2; // CMS energy squared
  Double_t fProb;  // Probability
#if defined(NEUT_WRITE_NUCFSI) or defined(NEUT_READ_NUCFSI)
  Double_t fStepPel;
  Double_t fStepPsp;
  Double_t fStepPdp;

  ClassDef(NeutNucFsiStep, 3)
#else
  ClassDef(NeutNucFsiStep, 2)
#endif
};

#endif // _NEUTNUCFSISTEP_H_
