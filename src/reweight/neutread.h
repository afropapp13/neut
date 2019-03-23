#ifndef NREWEIGHT_NEUTREAD_H_INCLUDED
#define NREWEIGHT_NEUTREAD_H_INCLUDED

// Headers containing common blocks from `${NEUT_ROOT}/include`
#include "neworkC.h"
#include "vcworkC.h"
#include "vcvrtxC.h"
#include "neutparamsC.h"
#include "posinnucC.h"
#include "neutmodelC.h"
#include "necardC.h"
#include "nrcardC.h"
#include "nefillverC.h"
#include "neutcrsC.h"
#include "fsihistC.h"
# if (NECOREVER == 5321)
  // special case for 5.3.2.1
  // do nothing.
# elif (NECOREVER == 5322) || (NECOREVER >= 535)
  //
  // special case for 5.3.2.2, which contains nucleon FSI
  //
#include "nucleonfsihistC.h"
#endif

#include "neutfilepathC.h"

class TTree;
class NeutVect;
class NeutVtx;

#ifdef __cplusplus
extern "C" {
#endif
  void readtree(TTree* tree, long long event);
  void readvect(NeutVect*, NeutVtx*);
  void readvtx(NeutVect*, NeutVtx*);

  void readnework(NeutVect*, NeutVtx*);
  void readvcwork(NeutVect*, NeutVtx*);
  void readposinnuc(NeutVect*, NeutVtx*);
  void readneutparams(NeutVect*, NeutVtx*);
  void readneutmodel(NeutVect*, NeutVtx*);
  void readnecard(NeutVect*, NeutVtx*);
  void readnrcard(NeutVect*, NeutVtx*);
  void readver(NeutVect*, NeutVtx*);
  void readneutcrs(NeutVect*, NeutVtx*);
  void readfsihist(NeutVect*, NeutVtx*);
  void readnucleonfsihist(NeutVect*, NeutVtx*);

  void readvcvrtx(NeutVect*, NeutVtx*);

#ifdef __cplusplus
};
#endif

namespace neut {
void PrintCommonBlocks();
}

#endif /* NREWEIGHT_NEUTREAD_H_INCLUDED */
