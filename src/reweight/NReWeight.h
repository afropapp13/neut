//____________________________________________________________________________
/*!

\class    neut::rew::NReWeight

\brief    Interface to the NEUT event reweighting engines

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

\created  Apr 5, 2011

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_H_
#define _N_REWEIGHT_H_

#include <string>
#include <map>

#include "NSystSet.h"
#include "NReWeightI.h"

using std::string;
using std::map;

namespace neut {


namespace rew   {

 class NReWeight
 {
 public:
   NReWeight();
  ~NReWeight();

   void        AdoptWghtCalc (string name, NReWeightI* wcalc);   ///< add concrete weight calculator, transfers ownership
   NReWeightI* WghtCalc      (string name);                      ///< access a weight calculator by name
   NSystSet &  Systematics   (void);                             ///< set of enabled systematic params & values
   void        Reconfigure   (void);                             ///< reconfigure weight calculators with new params
   double      CalcWeight    ();                                 ///< calculate weight for input event
   double      CalcChisq     (void);                             ///< calculate penalty chisq for current values of tweaking dials
   void        Print         (void);                             ///< print

   double      GetWeightXsec () {return fWeightXsec;};           ///< weight from cross section variations only
  private:

   void CleanUp (void);

   NSystSet                  fSystSet;   ///< set of enabled nuisance parameters
   map<string, NReWeightI *> fWghtCalc;  ///< concrete weight calculators

   double fWeightXsec;                   ///< store weight from cross section variations only
 };

} // rew   namespace
} // neut namespace

#endif

