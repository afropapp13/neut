//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightNuXSecNCEL

\brief    Reweighting NCEL NEUT neutrino cross sections

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

#ifndef _N_REWEIGHT_NU_XSEC_NCEL_H_
#define _N_REWEIGHT_NU_XSEC_NCEL_H_

#include <iostream>

#include "NReWeightI.h"
#include "NFortFns.h"
#include "NModeDefn.h"
#include "NTotCrs.h"

namespace neut {
namespace rew   {

 class NReWeightNuXSecNCEL : public NReWeightI 
 {
 public:
   static const int kModeMa             = 0;
   static const int kModeNormAndMaShape = 1;
   static const int kMode1overMa2 = 2;

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
   void SetMode     (int mode) { fMode       = mode; }
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }

 private:

   void   Init              (void);
   double CalcWeightNorm    ();
   double CalcWeightMaShape ();
   double CalcWeightMa      ();


   int    fMode;         ///< 0: Ma, 1: Norm and MaShape, 2: 1/Ma^2
   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?

   double fNormTwkDial;  ///<
   double fNormDef;      ///<
   double fNormCurr;     ///<

   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<

   float fAxlFFTwkDial;    ///<
   float fAxlFFDef;        ///<
   float fAxlFFCurr;       ///<

   float fVecFFTwkDial;    ///<
   float fVecFFDef;        ///<
   float fVecFFCurr;       ///<

   NFortFns * fortFns;
   NModeDefn modeDefn;

   NTotCrs *neutTotCrs;
 };

} // rew   namespace
} // neut namespace

#endif

