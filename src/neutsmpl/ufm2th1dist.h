
//  Class to convert a random number drawn from [0,1) uniform distribution
//  into a random number following a TH1D distribution.
//
//  2011.08: First version by S. Tobayama

#include <TH1D.h>

class Ufm2TH1dist {
private:
  Double_t *binEdgs;
  Double_t *binIntg;
  int nbin;
  
public:
  ~Ufm2TH1dist();
  void Init(TH1D* hist);
  double GetValue(double rnd);
};

