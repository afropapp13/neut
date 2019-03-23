#include <iostream>
#include <TRandom3.h>
#include <t2kneut_sk.h>

using namespace std;

extern "C"
{
  void ranlux_(float *, int *);

  void rluxgo_(int *, int *, int *, int *);

  void rluxat_(int *, int *, int *, int *);

};

void
ranlux_(float *dr, int *len)
{

  Int_t i;

  for ( i = 0 ; i < *len ; i++ ){
	dr[i] = global_tr3->Rndm();
  }

}

void
rluxgo_(int *lux,int *iseed, int *k1, int *k2)
{
  
  Int_t itseed;

  itseed = *iseed;

  global_tr3->SetSeed(itseed);
  /* reset k1 and k2 */

  cout << "Setting random seed for TRandom3 : seed = " 
	   << *iseed
	   << "\n";

}

void
rluxat_(int *lux,int *iseed, int *k1, int *k2)
{
  
  Int_t itseed;

  itseed = global_tr3->GetSeed();

  *iseed = itseed;
  *k1    = 0;
  *k2    = 0;

  //  cout << "Current random seed for TRandom3 : seed = " << *iseed;

}

