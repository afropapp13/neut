#ifndef _N_TOTCRS_H_
#define _N_TOTCRS_H_

#include <string>

#include "NSyst.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"


namespace neut {
namespace rew   {

 class NTotCrs
 { 
 public:
  static NTotCrs* Instance (void);
  
  NTotCrs();
  //NTotCrs(const NTotCrs & err);
  ~NTotCrs();

  static NTotCrs * fInstance;
  
  std::string neut_folder;
  void LoadCCQE();
  TFile *ccqeCrsFile;
  TH3D *ccqe_crs[4];    ///< CCQE cross section histogram as a function of Enu, PFSurf, MA (for MA-shape variations)

  void LoadRESSPI();
  TFile *resspiCrsFile;
  TH2D *resspi_crs[4][7]; ///< ResSpi cross section histograms as a function of Enu, MA (for MA-shape variations)
  enum neutResSpiModes {neutmode11,neutmode12,neutmode13,neutmode32,neutmode34,neutmode31,neutmode33};

  enum neutTotCrsNuTypes {numu,numub,nue,nueb};

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (NTotCrs::fInstance !=0) {
            delete NTotCrs::fInstance;
            NTotCrs::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;

 };

} // rew   namespace
} // neut namespace

#endif

