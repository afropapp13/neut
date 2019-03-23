//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightNuXSecNCEL

\brief    Reweighting NCEL NEUT neutrino cross sections

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

\created  Nov 25, 2010

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_NU_XSEC_NCEL_H_
#define _N_REWEIGHT_NU_XSEC_NCEL_H_


#include "NReWeightI.h"
#include "NFortFns.h"
#include "NModeDefn.h"

namespace neut {
namespace rew   {

 class NReWeightNuXSecNCEL : public NReWeightI 
 {
 public:
   NReWeightNuXSecNCEL();
  ~NReWeightNuXSecNCEL();

   // implement the NReWeightI interface
   bool   IsHandled      (NSyst_t syst);
   void   SetSystematic  (NSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     ();
   double CalcChisq      (void);

   // various config options
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }

 private:

   void Init(void);

   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<
   //double fEtaTwkDial;    ///<
   //double fEtaDef;        ///<
   //double fEtaCurr;       ///<

   NFortFns * fortFns;
   NModeDefn modeDefn;
   
 };

} // rew   namespace
} // neut namespace

#endif

