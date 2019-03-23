#include <iostream>
#include <TH1D.h>
#include <TUnuran.h>
#include <TUnuranEmpDist.h>
#include <TRandom3.h>

int
main(int argc, char **argv)
{

  TH1D *histo;
  int i;
  Double_t r,w;
  
  TRandom3 *tr3;
  TUnuranEmpDist *unr_emp;
  TUnuran        *unr;

  tr3 = new TRandom3();

  histo = new TH1D("sample","sampletitle",200,0,10);

  for ( i = 0 ; i < 100 ; i++){
	r = tr3->Rndm();
	w = tr3->Rndm();
	histo->Fill(r,w);
  }

  unr_emp = new TUnuranEmpDist(histo);
  unr     = new TUnuran();
  unr->Init(*unr_emp,"empk");
  

}


