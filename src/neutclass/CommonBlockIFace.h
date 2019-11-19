#ifndef _N_COMMONBLOCKIFACE_H_
#define _N_COMMONBLOCKIFACE_H_

#include "necardC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"

#include "neutvect.h"

#include <string>

extern "C" {
void nesetfgparams_();
void nefillmodel_();
void zexpconfig_();
double evdifcrs_();
double evpiprob_();
}

namespace neut {

inline double NEUTGetXSec() { return evdifcrs_(); }
inline double NEUTGetPiCascProb() { return evpiprob_(); }
inline void NEUTSetParams() {
  nefillmodel_();
  nesetfgparams_();
  zexpconfig_();
}

class CommonBlockIFace {

public:
  // necardC.h
  neutcard_common fneutcard_gen;
  nuceffver_common fnuceffver_gen;
  neutdis_common fneutdis_gen;
  neut1pi_common fneut1pi_gen;
  neutdif_common fneutdif_gen;
  neutcoh_common fneutcoh_gen;
  neuttarget_common fneuttarget_gen;
  neutpiabs_common fneutpiabs_gen;
  neutpiless_common fneutpiless_gen;
  neutradcorr_common fneutradcorr_gen;

  // neutmodelC.h
  nemdls_common fnemdls_gen;
  // neutparamsC.h
  nenupr_common fnenupr_gen;
  neffpr_common fneffpr_gen;

private:
  void SetGenCard(std::string const &GenCardLocation);

public:
  static void Initialize(std::string const &GenCardLocation);
  static CommonBlockIFace const &Get();

  void ResetGenValues() const;

  // Reads event properties straight into the common blocks
  static void ReadNEWORK(NeutVect *);
  static void ReadVCWORK(NeutVect *);
  static void ReadPOSINNUC(NeutVect *);
  static void ReadNEUTCRS(NeutVect *);
  static void ReadFSIHIST(NeutVect *);
  static void ReadNUCLEONFSIHIST(NeutVect *);
  static void ReadVect(NeutVect *);

  static std::string ParamsToString(bool isinstance = false);
};

} // namespace neut

#endif
