#ifndef _N_COMMONBLOCKIFACE_H_
#define _N_COMMONBLOCKIFACE_H_

#include "necardC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"

#include "cardnameC.h"
#include "necardevC.h"
#include "neutfilepathC.h"
#include "nieves1p1hC.h"
#include "nrcardC.h"

#include "neutvect.h"

#ifdef USE_HEPMC
namespace HepMC3 {
class GenEvent;
class GenRunInfo;
class Attribute;
} // namespace HepMC3
#include <memory>
#endif

#include <map>
#include <string>

extern "C" {
void nesetfgparams_();
void nefillmodel_();
void zexpconfig_();
double evdifcrs_();
double evpiprob_();
double evpiw_();
void silencelibgfortran_();
void unsilencelibgfortran_();
}

namespace neut {

inline double NEUTGetXSec() { return evdifcrs_(); }
inline double NEUTGetPiCascProb() { return evpiprob_(); }
inline double NEUTGetPiEj() { return evpiw_(); }
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
  neutpiabs_common fneutpiabs_gen;
  neutpiless_common fneutpiless_gen;
  neutradcorr_common fneutradcorr_gen;

  // neutmodelC.h
  nemdls_common fnemdls_gen;
  neutmodel_common fneutmodel_gen;

  // neutparamsC.h
  nenupr_common fnenupr_gen;
  neffpr_common fneffpr_gen;

  // cardnameC.h
  necardname_common fnecardname_gen;

  // necardevC.h
  nevccard_common fnevccard_gen;

  // nrcardC.h
  nucres_common fnucres_gen;

  // neutfilepathC.h
  neutfilepath_common fneutfilepath_gen;

  // nieves1p1hC.h
  nievesqepar_common fnievesqepar_gen;

private:
  void SetGenCard(std::string const &GenCardLocation);

public:
  static void Initialize(std::string const &GenCardLocation,
                         bool quiet = false);
  static CommonBlockIFace const &Get();

#ifdef USE_HEPMC
  static void Initialize(std::shared_ptr<HepMC3::GenRunInfo const> gri);
  static void ReadEvent(std::shared_ptr<HepMC3::GenEvent> evt);

  static std::map<std::string, std::shared_ptr<HepMC3::Attribute> >
  SerializeModelCommonBlocksToAttributeList();
#endif

  void ResetGenValues() const;

  // Reads event properties straight into the common blocks
  static void ReadNEWORK(NeutVect *);
  static void ReadVCWORK(NeutVect *);
  static void ReadPOSINNUC(NeutVect *);
  static void ReadNEUTCRS(NeutVect *);
  static void ReadFSIHIST(NeutVect *);
  static void ReadNUCLEONFSIHIST(NeutVect *);
  static void ReadNEUTTARGET(NeutVect *);
  static void ReadNENUPR(NeutVect *);

  static void ReadVect(NeutVect *);
  static void ReadEVENTPARS(NeutVect *);

  static std::string ParamsToString(bool isinstance = false);
};

} // namespace neut

#endif
