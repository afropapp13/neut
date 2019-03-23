//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightBindingEnergy

\brief    Reweighting CCQE NEUT neutrino cross sections
          applying corrections exclusively to the binding energy

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Patrick de Perio <pdeperio \at physics.utoronto.ca>
          University of Toronto

          Matt Dunkman <mdunkman \at msu.edu>
          Michigan State University

\created  Aug 1, 2009
\modified  Jun 30, 2016

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_BINDING_ENERGY_H_
#define _N_REWEIGHT_BINDING_ENERGY_H_

#include <iostream>
#include <string>

#include "NReWeightI.h"
#include "NModeDefn.h"
#include "NSyst.h"

class TTree;
class TH3F;
class TH2F;

namespace neut {
namespace rew   {

namespace binding_energy {
 class TargetMaterial {
  public:
   const int atomic_number;
   const double default_value;
   const double fractional_uncertainty;

  protected:
   explicit TargetMaterial(int atomic_number_, double default_value_, double uncertainty_);
   virtual ~TargetMaterial() { }
 };

 class CarbonTarget : public TargetMaterial {
  public:
   explicit CarbonTarget(double default_value_=25., double uncertainty_=9.);
   ~CarbonTarget() { }
 };

 class OxygenTarget : public TargetMaterial {
  public:
   explicit OxygenTarget(double default_value_=27., double uncertainty_=9.);
   ~OxygenTarget() { }
 };

 class AluminumTarget : public TargetMaterial {
  public:
   explicit AluminumTarget(double default_value_=28., double uncertainty_=9.);
   ~AluminumTarget() { }
 };

 class IronTarget : public TargetMaterial {
  public:
   explicit IronTarget(double default_value_=33., double uncertainty_=9.);
   ~IronTarget() { }
 };

 class CopperTarget : public TargetMaterial {
  public:
   explicit CopperTarget(double default_value_=35., double uncertainty_=9.);
   ~CopperTarget() { }
 };

 class ZincTarget : public TargetMaterial {
  public:
   explicit ZincTarget(double default_value_=35., double uncertainty_=9.);
   ~ZincTarget() { }
 };

 class LeadTarget : public TargetMaterial {
  public:
   explicit LeadTarget(double default_value_=44., double uncertainty_=9.);
   ~LeadTarget() { }
 };

}

 class NReWeightBindingEnergy: public NReWeightI
 {
 public:
   explicit NReWeightBindingEnergy(
       const binding_energy::TargetMaterial& target=binding_energy::CarbonTarget(),
       int desired_mdlqe=2, double weight_cap=100., bool debug=false, bool trace=false
       );
   ~NReWeightBindingEnergy();

   // implement the NReWeightI interface
   bool IsHandled(NSyst_t syst);
   void SetSystematic(NSyst_t syst, double val);
   void Reset();
   void Reconfigure();
   double CalcWeight(TTree* tree, long long event);
   double CalcWeight();
   double CalcChisq();

   //TH3F GenerateLookupTable(TTree* tree);
   TH2F GenerateLookupTable(TTree* tree);
   TH3F GenerateLookupTable();

   // various config options
   void RewNue(bool tf) { fRewNue = tf; }
   void RewNuebar(bool tf) { fRewNuebar = tf; }
   void RewNumu(bool tf) { fRewNumu = tf; }
   void RewNumubar(bool tf) { fRewNumubar = tf; }
   void Absolute(bool value) { fAbsolute = value; }

 private:

   void Init();
   NSyst_t GetSystematic() const;
   bool IsQEModelValid() const;
   bool GoodNeutrino() const;
   bool EventIsOkayForDesiredMDLQE() const;
   size_t GetIndexOfOutgoingLepton() const;
   double CalcWeightI();
   double CalcXsecSF();
   double CalcXsecSF(const double k[], const double p[], const double in_nucleon_mass,
                     const double kPrime[], const double pPrime[], const double out_nucleon_mass);
   double CalcXsecRFG();


   bool fRewNue;       ///< reweight nu_e CC?
   bool fRewNuebar;    ///< reweight nu_e_bar CC?
   bool fRewNumu;      ///< reweight nu_mu CC?
   bool fRewNumubar;   ///< reweight nu_mu_bar CC?
   bool fAbsolute;     ///< is the value passed into SetSystematic() absolute, or in sigma?
   bool fDebug;        ///< print messages when unexpected things happen?
   bool fTrace;        ///< print a whole bunch of junk?


   int fDesiredMDLQE;   ///<
   const double fWeightCap;   ///< maximum value to return by CalcWeight()
   float fEbTwkDial;    ///<
   float fEbCurr;       ///<
   float fEbDef;        ///<
   float fEbUncertainty;///<
   int fNumAtom;        ///< Atomic number of target nucleon
   NSyst_t fSystematic; ///< enum member

   //NFortFns * fortFns;
   NModeDefn modeDefn;

   //NTotCrs *neutTotCrs;
 };


} // rew   namespace
} // neut namespace

#endif

