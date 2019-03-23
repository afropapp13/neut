#include <iostream>
#include <TFile.h>
#include <stdlib.h>
#include "t2kflux_sk.h"

TRandom3 *global_tr3;

using namespace std;

extern "C"{
  void hlimit_(int *);
  void necard_();
  void necardbm_();

  /*
  void rmarin_(int *,int *,int *);
  void rmarin_(int *i,int *j, int *k)
  {
	printf("rmarin was called.\n");
	return;
  }
  */
#ifdef gFortran
  void _gfortran_set_args(int, char **);
#else
  void f_setarg(int, char **);
#endif
};

struct PAWC{
  float H[5000000];
} pawc_;

int
main(int argc, char **argv)
{
  int i,j;
  int    ierr;

  char fname[1024];

  int hsiz = 5000000;

  global_tr3 = new TRandom3();

  int idflux  = 14;
  int idcrs   = 14;
  double enusk;

  TFile *tf;

#ifdef gFortran
  _gfortran_set_args(argc, argv);
#else
  f_setarg(argc, argv);
#endif

  if ( argc < 1 ){
	cout << "Usage : " << argv[0] << " [output root flux file]" << endl;
	exit(1);
  }
  necard_();
  necardbm_();

  strncpy(fname,argv[2],sizeof(fname));

  tf    = new TFile(fname,"NEW");
  if (tf->IsZombie()){
	cout << "Failed to open file." << fname << endl;
	exit(1);
  }

  hlimit_(&hsiz);
  T2Kflux_SK *skflux = new T2Kflux_SK(global_tr3);

  enusk = skflux->rand_enu(idflux,idcrs,&ierr);
	   
  for ( i = 0 ; i < 4 ; i++ ){
	(skflux->fluxhisto[i])->Write();
	for ( j = 0 ; j < 4 ; j++ ){
	  (skflux->ratehisto[i][j])->Write();
	}
  }
  tf->Close();
  
}
