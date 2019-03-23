//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightNuXSecCOH

\brief    Reweighting NEUT coherent neutrino-nucleus cross sections

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

#ifndef _N_REWEIGHT_NU_XSEC_COH_H_
#define _N_REWEIGHT_NU_XSEC_COH_H_

#include "NReWeightI.h"
#include "NFortFns.h"
#include "NModeDefn.h"

namespace neut {
namespace rew   {

 class NReWeightNuXSecCOH : public NReWeightI 
 {
 public:
   NReWeightNuXSecCOH();
  ~NReWeightNuXSecCOH();

   // implement the NReWeightI interface
   bool   IsHandled      (NSyst_t syst);
   void   SetSystematic  (NSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     ();
   double CalcChisq      (void);

   // various config options
   void RewNue      (bool tf ) { fRewNue     = tf; }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf; }
   void RewNumu     (bool tf ) { fRewNumu    = tf; }
   void RewNumubar  (bool tf ) { fRewNumubar = tf; }
   void RewCC       (bool tf ) { fRewCC      = tf; }
   void RewNC       (bool tf ) { fRewNC      = tf; }

 private:

   void Init (void);

   bool   fRewNue;       ///< reweight nu_e?
   bool   fRewNuebar;    ///< reweight nu_e_bar?
   bool   fRewNumu;      ///< reweight nu_mu?
   bool   fRewNumubar;   ///< reweight nu_mu_bar?
   bool   fRewCC;        ///< reweight CC?
   bool   fRewNC;        ///< reweight NC?

   double fNECOHEPITwkDial;    ///<
   double fNECOHEPIDef;        ///<
   double fNECOHEPICurr;       ///<

   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<

   double fR0TwkDial;    ///<
   double fR0Def;        ///<
   double fR0Curr;       ///<

   double fA1TwkDial;    ///<
   double fA1Def;        ///<
   double fA1Curr;       ///<  

   double fb1TwkDial;    ///<
   double fb1Def;        ///<
   double fb1Curr;       ///<

   NFortFns * fortFns;
   NModeDefn modeDefn;
 };

} // rew   namespace
} // neut namespace

#endif

