//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightNuclPiless

\brief    Reweight NEUT CC resonance neutrino-production

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_NUCL_PILESS_H_
#define _N_REWEIGHT_NUCL_PILESS_H_

#include "NReWeightI.h"
#include "NFortFns.h"
#include "NModeDefn.h"

namespace neut {
namespace rew   {

 class NReWeightNuclPiless : public NReWeightI 
 {
 public:

   NReWeightNuclPiless();
  ~NReWeightNuclPiless();

   // implement the NReWeightI interface
   bool   IsHandled      (NSyst_t syst);
   void   SetSystematic  (NSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     ();
   double CalcChisq      (void);

   // various config options
   void SetModeEnu  (double aModeEnu) { fModeEnu = aModeEnu; }
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }

 private:

   void   Init                (void);
   double CalcWeightPilessDcy ();

   double fModeEnu;      ///< 0: Apply to all energies and weight all events (Default)
                         ///< fModeEnu>0: Apply to Enu<fModeEnu[GeV] and weight all events
                         ///< fModeEnu<0: Apply to Enu<-fModeEnu[GeV] and weight only PDD events
   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?

   double fPilessDcyTwkDial;
   double fPilessDcyDef;      ///<
   double fPilessDcyCurr;     ///<

   NFortFns * fortFns;
   NModeDefn modeDefn;
 };

} // rew   namespace
} // neut namespace

#endif

