#include <iostream>
#include "skflux.h"
#include <TFile.h>

using namespace std;

extern "C"{
  void hlimit_(int *);
  void necard_();

};

struct PAWC{
  float H[5000000];
} pawc_;

int
main()
{
  int i;
  int    ierr;
  double enu;

  int hsiz = 5000000;

  TFile *tf;
  TH1D  *edist;

  necard_();

  tf    = new TFile("edist.root","NEW");
  if (tf->IsZombie()){
	cout << "Failed to open file.";
	exit(1);
  }

  edist = new TH1D("edist","Energy@SK",100,0.,10.);

  hlimit_(&hsiz);
  SKflux *flux_rand = new SKflux();

  for( i = 0 ; i < 1000000 ; i++){
	enu = flux_rand->rand_enuflux(14,&ierr);
	edist->Fill(enu);
	cout <<  i << " " << enu << "\n";
  }

  edist->Write();
  tf->Close();
  
}
