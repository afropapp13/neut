#ifndef SEEN_NEUTCINCLUDES_HXX
#define SEEN_NEUTCINCLUDES_HXX
//Provides necard common block
#include "necardC.h"
//Provides eftarget common block
#include "efpionC.h"
//Provides vcwork common block
#include "vcworkC.h"
//Provides posinnuc common block
#include "posinnucC.h"


#ifdef __cplusplus
extern "C" {
#endif

extern struct nevccard_common {
  Int_t nectecvt;
  Int_t idptevct;
  Int_t mposevct;

  Float_t posevct[3];
  Float_t radevct;

  Int_t mdirevct;
  Float_t direvct[3];

  Int_t mpvevct;
  Float_t pvevct[2];

  char filenmevct[80];
  char histnmevct[80];

  Int_t inmececvt;

} nevccard_;
#ifdef STATIC_COMMON_POINTERS
static struct nevccard_common *nevccard = &nevccard_;
#endif

  //Defined in neutcore/necard.F
  void necard_();
  //Defined in skmcsvc/necardev.F
  void necardev_();
  float rlu_(int*);

  // skmcsvc/vcclcm.F
  // Clears VCWORK and VCVRTX common block arrays
  void vcclcm_();

  // nuceff/efclfsi.F
  // Clears FSIHIST common block arrays
  void efclfsi_();

  //nuccorspl/nrintr.F
  //chases and simulates nucleons from the VC work block
  void nrintr_();

  //skmcsvc/mcmass.F
  void mcmass_(int*, float*);

  //nuceff/evpiprob.F
  float evpiprob_();

  //nuceff/efpqeabh.F
  void efpqeabh_(float*,float*,float*,int*,int*,float*,float*,float*,float*,float*);
  
  void nesettarg_();

  //nuccorspl/nrfermi.F
  void nrfermi_(float*);

  //nuccorspl/nrnuc.F
  void nrnuc_(int*);
  
  //nuccorspl/nrprton.F
  void nrprton_(float*, float*, float*, int*, float*, int*, float*, int*, int*, int*);




#ifdef __cplusplus
} // closes extern "C"
#endif

template<size_t N>
std::ostream &operator<<(std::ostream& os, Float_t (&v)[N] ){
  os << "[" << std::flush;

  for(size_t i = 0; i < (N-1); ++i){
    os << v[i] << ", " << std::flush;
  }
  os << v[N-1] << "]" << std::flush;
  return os;
}


void DebugVals(){
  std::cout << "eftarget: {" << std::endl;
  std::cout << "\tan: " << eftarget->an << std::endl;
  std::cout << "\tzz: " << eftarget->zz << std::endl;
  std::cout << "\tc: " << eftarget->c << std::endl;
  std::cout << "\tcnn: " << eftarget->cnn << std::endl;
  std::cout << "\taf: " << eftarget->af << std::endl;
  std::cout << "\tcc2: " << eftarget->cc2 << std::endl;
  std::cout << "\tdfact: " << eftarget->dfact << std::endl;
  std::cout << "\trmsrad: " << eftarget->rmsrad << std::endl;
  std::cout << "\twparm: " << eftarget->wparm << std::endl;
  std::cout << "}" << std::endl;

  std::cout << "nevccard: {" << std::endl;
  std::cout << "\tnectecvt: " << nevccard->nectecvt << std::endl;
  std::cout << "\tidptevct: " << nevccard->idptevct << std::endl;
  std::cout << "\tmposevct: " << nevccard->mposevct << std::endl;
  std::cout << "\tmdirevct: " << nevccard->mdirevct << std::endl;
  std::cout << "\tmpvevct: " << nevccard->mpvevct << std::endl;
  std::cout << "\tposevct: " << nevccard->posevct << std::endl;
  std::cout << "\tradevct: " << nevccard->radevct << std::endl;
  std::cout << "\tdirevct: " << nevccard->direvct << std::endl;
  std::cout << "\tpvevct: " << nevccard->pvevct << std::endl;
  std::cout << "\tinmececvt: " << nevccard->inmececvt << std::endl;
  std::cout << "}" << std::endl;
}

#endif
