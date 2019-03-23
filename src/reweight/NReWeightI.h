//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightI

\brief    NEUT event reweighting engine ABC 

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_ABC_H_
#define _N_REWEIGHT_ABC_H_

#include "NSyst.h"
#include "NSystSet.h"

namespace neut {

namespace rew   {

 class NReWeightI 
 {
 public:
  virtual ~NReWeightI() 
  { 

  } 

  //
  // define the NReWeightI interface
  //

  //! does the current weight calculator handle the input nuisance param?
  virtual bool IsHandled (NSyst_t syst) = 0; 

  //! update the value for the specified nuisance param
  virtual void SetSystematic (NSyst_t syst, double val) = 0; 

  //!  set all nuisance parameters to default values
  virtual void Reset (void) = 0;            

  //! propagate updated nuisance parameter values to actual MC, etc
  virtual void Reconfigure (void) = 0;            
  
  //! calculate a weight for the input event using the current nuisance param values
  virtual double CalcWeight () = 0;  

  //! calculate penalty factors
  virtual double CalcChisq (void) = 0;        

 protected:

   NReWeightI() 
   { 
   
   }
 };

} // rew   namespace
} // neut namespace

#endif

