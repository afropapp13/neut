#include <iostream>
#include <necardbmC.h>
#include <necardC.h>
#include <math.h>
#include <cstdlib>

#define MAIN
#include "t2kneut_sk.h"
#include "t2kflux_sk.h"

#include "skheadC.h"


#define ISIZE 5000000
#define LUNO  20
#define PI    3.1415926535

using namespace std;

struct PAWC{
  float h[ISIZE];
} pawc_;

extern "C"
{
  void necard_();
  void necardbm_();

  void kzinit_();
  void kzeclr_();
  void kzwrit_(int *);
  void kzend_(int *);

  void hlimit_(int *);

  void skopenf_(int *, int *, char *, int *, int);

  void nevect_(int *, float *, float *, int *);

  void rnpos_(float *, float *);

  void nemknebk_(float *);
  void nemkcrsbk_();
  void nemkmodelbk_();
  void nemkfsibk_();
  void nemknpos_();
  void nemknetarg_();
  
  void vcmkvc_();
  void vcmkvx_();

  void mcmkhd_();
  void mcmkmh_();
  void mcmkwt_();
#ifdef gFortran
  void _gfortran_set_args(int, char **);
#else
  void f_setarg(int, char **);
#endif

}

int
main(int argc, char **argv)
{

  int isize   = -1 * ISIZE;
  int luno    = LUNO;
  int file_no = 1;
  char chopt[]="Z\0";

  int idflux  = 14;
  int idcrs   = 14;

  float dwall = -50;

  float pos[3],pmom[3];

  int i,ierr;
  
  T2Kflux_SK *skflux;

  float enusk;

  /*******************************************************/
#ifdef gFortran
  _gfortran_set_args(argc, argv);
#else
  f_setarg(argc, argv);
#endif

  global_tr3 = new TRandom3();

  necard_();
  necardbm_();

  idflux = nubmgen_.idbmflx;
  idcrs  = nubmgen_.idbmpid;
	

  kzinit_();
  hlimit_(&isize);

  skflux = new T2Kflux_SK(global_tr3);

  skopenf_(&luno, &file_no, chopt, &ierr, strlen(chopt));
  if (ierr != 0){
	cout << "Failed to open file (LUN=" << luno << ".\n";
	exit(1);
  }
  
  for ( i = 0 ; i < nubmgen_.nevbm ; i++){
	kzeclr_();
	do{
           ierr = 0;
	   enusk = skflux->rand_enu(idflux,idcrs,&ierr);
	   if (ierr != 0){
	     cout << "Failed to fix energy for event #" << i << ".\n";
	     exit(2);
	   }else{
	     if (neutcard_.quiet==0)
	       cout << "Event #" << i << " : Energy = " << enusk << "GeV.\n";
	   }
	
	   rnpos_(pos,&dwall);

	   enusk = enusk * 1000.;
	   pmom[2] = enusk * sin(1.388*PI/180.);
	   pmom[0] = sqrt(enusk*enusk-pmom[2]*pmom[2]);
	   pmom[1] = 0;

	   pmom[1] = sin(312.064*PI/180.)*pmom[0];
	   pmom[0] = sqrt(pmom[0]*pmom[0]-pmom[1]*pmom[1]);
	
	   nevect_(&idcrs, pos, pmom, &ierr);
	}while(ierr != 0);

	nemknebk_(pos);
	nemkmodelbk_();
	nemkfsibk_();
	nemkcrsbk_();
	vcmkvc_();
	vcmkvx_();
	nemknpos_();

	skheadg_.sk_geometry = 4;
	nemknetarg_();

	mcmkhd_();
	mcmkmh_();
	mcmkwt_();

	kzwrit_(&luno);
  }
  kzend_(&luno);

  exit(0);
}
