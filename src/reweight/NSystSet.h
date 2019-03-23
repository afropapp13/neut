//____________________________________________________________________________
/*!

\class    neut::rew::NSystSet

\brief    Set of systematics to be considered by the reweighting package.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_SET_OF_SYSTEMATICS_H_
#define _N_SET_OF_SYSTEMATICS_H_

#include <string>
#include <map>
#include <vector>

#include "NSyst.h"

using std::string;
using std::map;
using std::vector;

namespace neut {
namespace rew   {

class NSystInfo;

class NSystSet {

public:  
  NSystSet();
  NSystSet(const NSystSet & syst_set);
 ~NSystSet();

  void Init    (NSyst_t syst, double init=0., double min=-1., double max=+1., double step=0.05);
  void Remove  (NSyst_t syst);
  void Set     (NSyst_t syst, double current_value);

  int  Size    (void) const;
  bool Added   (NSyst_t syst) const;
  void Print   (void);
  void Copy    (const NSystSet & syst_set);

  const NSystInfo * Info(NSyst_t syst) const;

  vector<neut::rew::NSyst_t> AllIncluded (void);

private:
  
  map<NSyst_t, NSystInfo *>  fSystematics;  
};

class NSystInfo {

public:
  NSystInfo() : 
     CurValue(0), InitValue(0), MinValue(0), MaxValue(0), Step(0) 
  { 

  }
  NSystInfo(double init, double min, double max, double step) : 
     CurValue(init), InitValue(init), MinValue(min), MaxValue(max), Step(step) 
  {
 
  }
 ~NSystInfo() 
  { 

  }

  double CurValue;
  double InitValue;
  double MinValue;
  double MaxValue;
  double Step;
};

} // rew   namespace
} // neut namespace

#endif 

