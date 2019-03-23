//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightNuXSecDIS

\brief    Reweighting GENIE DIS neutrino-nucleus cross sections

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

#ifndef _N_REWEIGHT_NU_XSEC_DIS_H_
#define _N_REWEIGHT_NU_XSEC_DIS_H_

#include "NReWeightI.h"
#include "NFortFns.h"
#include "NModeDefn.h"

namespace neut {
namespace rew   {

 class NReWeightNuXSecDIS : public NReWeightI 
 {
 public:
   //static const int kModeABCV12u      = 0;
   //static const int kModeABCV12uShape = 1;

   NReWeightNuXSecDIS();
  ~NReWeightNuXSecDIS();

   // implement the NReWeightI interface
   bool   IsHandled      (NSyst_t syst);
   void   SetSystematic  (NSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     ();
   double CalcChisq      (void);

   // various config options
   void SetMode      (int    m )  { fMode       = m;  }
   void RewNue       (bool   tf)  { fRewNue     = tf; }
   void RewNuebar    (bool   tf)  { fRewNuebar  = tf; }
   void RewNumu      (bool   tf)  { fRewNumu    = tf; }
   void RewNumubar   (bool   tf)  { fRewNumubar = tf; }
   void RewCC        (bool   tf)  { fRewCC      = tf; }
   void RewNC        (bool   tf)  { fRewNC      = tf; }

 private:

   void   Init                   (void);
   double CalcWeightNorm     ();
   double CalcWeightBYOnOff     ();

   int    fMode;                 ///< 0: DIS+MPI, 1: DIS-only, 2: MPI-only
   bool   fRewNue;               ///< reweight nu_e?
   bool   fRewNuebar;            ///< reweight nu_e_bar?
   bool   fRewNumu;              ///< reweight nu_mu?
   bool   fRewNumubar;           ///< reweight nu_mu_bar?
   bool   fRewCC;                ///< reweight CC?
   bool   fRewNC;                ///< reweight NC?

   double fNormTwkDial;  ///<
   double fNormDef;      ///<
   double fNormCurr;     ///<

   double fBYOnOffTwkDial;  ///<
   double fBYOnOffDef;      ///<
   double fBYOnOffCurr;     ///<

   NFortFns * fortFns;
   NModeDefn modeDefn;
   
 };

} // rew   namespace
} // neut namespace

#endif

