#include <iostream>
#include <TFile.h>
#include <stdlib.h>
#include "t2kflux_sk.h"

using namespace std;

extern "C"{
  void necard_();
  void necardbm_();

#ifdef gFortran
  void _gfortran_set_args(int, char **);
#else
  void f_setarg(int, char **);
#endif
};

int
main(int argc, char **argv)
{
  int i,j;
  int    ierr;

  char fname_in_histo[1024];
  char fname_out_histo[1024];

  int hsiz = 5000000;

  int idflux  = 14;
  int idcrs   = 14;
  double enusk;

  TFile *f_histo_in,*f_histo_out 

#ifdef gFortran
  _gfortran_set_args(argc, argv);
#else
  f_setarg(argc, argv);
#endif

  if ( argc < 1 ){
	cout << "Usage : " 
		 << argv[0] 
		 << "[card file] " 
		 << "[input root flux file] " 
		 << "[output root flux file]" 
		 << endl;
	exit(1);
  }
  necard_();
  necardbm_();

  strncpy(fname,argv[1],sizeof(fname_in_histo));
  strncpy(fname,argv[2],sizeof(fname_out_histo));

  f_histo_in    = new TFile(fname_in_histo);
  if (f_histo_in->IsZombie()){
	cout << "Failed to open file." << fname_in_histo << endl;
	exit(1);
  }
  
  f_histo_out    = new TFile(fname_out_histo,"NEW");
  if (f_histo_out->IsZombie()){
	cout << "Failed to open file." << fname_out_histo << endl;
	exit(1);
  }
  
  for ( i = 0 ; i < 4 ; i++ ){
	(skflux->fluxhisto[i])->Write();
	for ( j = 0 ; j < 4 ; j++ ){
	  (skflux->ratehisto[i][j])->Write();
	}
  }
  tf->Close();
  
}
