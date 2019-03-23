//____________________________________________________________________________
/*!

\class    neut::rew::NReWeightINuke

\brief    Reweighting NEUT Nucleon FSI Cascade Model

\author   Martin Hierholzer <martin@hierholzer.info>
          University of Bern

\created  Nov 04, 2014
*/
//____________________________________________________________________________

#ifndef _N_REWEIGHT_INUKE_H_
#define _N_REWEIGHT_INUKE_H_

#include <vector>

#include "NReWeightI.h"

#include "NFortFns.h"

using namespace neut::rew;
using namespace neut;

namespace neut {
namespace rew   {

 class NReWeightINuke : public NReWeightI 
 {
 public:
   NReWeightINuke();
  ~NReWeightINuke();

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

   double fMFPDial;     ///< mean free path

   double fMFPCurr;     ///< mean free path
     
   double fMFPDef;      ///< mean free path

   // some function converted from Fortran used to calculate the probability for rescattering
   void nrsettarg();
   double nrROXY(double r);
   void nrfermi(double r, double *pin, int N);
   void nrhis(double x, double rbin, int num, int &ilo, int &ihi, double &erem);
   void nrnuc(int ichint);


   // some nuclear properties and variables changing per step or interaction
   int A;                   // number of nucleons
   int Z;                   // number of protons
   double NRRMSRAD;
   double NRC;
   double NRCNN;
   double NRAF;
   double NRWPARM;
   double NRDFACT;
   double NRPFSURF;         // fermi momentum (in GeV)
       
   int nbin;                // number of bins for nuclear density
   double rhotab[200];      // nuclear density in NBIN bins
       
   double pnorm;
       
   double ecms2;            // CMS energy for the current (possible) interaction - set for each step!

   double unucl;            // ??? from nrfermi()
   double prat;             // faction of protons in nucleos, from nrfermi()
   double rhon;             // local nuclear density, from nrfermi()

   double ptot;             // total interaction probability at current step
   double ppel;
   double ppsp;
   double ppdp;
   double pnel;
   double pnsp;
   double pndp;

   double ptot_s;           // total interaction probability at current step
   double ppel_s;
   double ppsp_s;
   double ppdp_s;
   double pnel_s;
   double pnsp_s;
   double pndp_s;
       
   int ichtrgt;             // charge of target nucleon in collision

 };

} // rew
} // neut

#endif

