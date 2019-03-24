//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightCasc

\brief    Reweighting NEUT FSI Cascade Model

          See T2K-TN-033

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_CASC_H_
#define _N_REWEIGHT_CASC_H_

#include <vector>

#include "NReWeightI.h"

#include "NFortFns.h"

using namespace neut::rew;
using namespace neut;

namespace neut {
namespace rew   {

 class NReWeightCasc : public NReWeightI 
 {
 public:
   NReWeightCasc();
  ~NReWeightCasc();

   // implement the NReWeightI interface
   bool   IsHandled      (NSyst_t syst);
   void   SetSystematic  (NSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     ();
   double CalcChisq      (void);

   vector<double> GetCurrParVals(void);

 private:
   
   void   Init              (void);

   double fCExLowTwkDial;    ///<
   double fInelLowTwkDial;   ///<
   double fAbsTwkDial;       ///<
   double fPiProdTwkDial;    ///<
   double fCExHighTwkDial;   ///<
   double fInelHighTwkDial;  ///<
   double fAllTwkDial;       ///<
   
   double fCExLowCurr;       ///<
   double fInelLowCurr;      ///<
   double fAbsCurr;          ///<
   double fPiProdCurr;       ///<
   double fCExHighCurr;      ///<
   double fInelHighCurr;     ///<
   double fAllCurr;          ///<
     
   double fCExLowDef;        ///<
   double fInelLowDef;       ///<
   double fAbsDef;           ///<
   double fPiProdDef;        ///<
   double fCExHighDef;       ///<
   double fInelHighDef;      ///< 
   double fAllDef;           ///< 

   NFortFns * fortFns;

 };

} // rew
} // neut

#endif

