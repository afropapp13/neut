//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory 

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ May 18, 2010 - CA
   AdoptWghtCalc(string,NReWeightI*) allows user to decide which weight 
   calculator to include. Weight calculators are owned by NReWeight and are 
   identified by a name. Weight calculators can be retrieved via the 
   WghtCalc(string) method and their reweighting options can be fine-tuned.
 @ Apr 6, 2011 - PD
   Implemented in NEUT
*/
//____________________________________________________________________________

#include <iostream>
#include <vector>

#include <TMath.h>
#include <TString.h>

#include "NReWeight.h"

using std::vector;
using std::cout;
using std::endl;

//#define _NEUT_REWEIGHT_DEBUG_

using namespace neut;
using namespace neut::rew;

//____________________________________________________________________________
NReWeight::NReWeight()
{

}
//____________________________________________________________________________
NReWeight::~NReWeight()
{
  this->CleanUp();
}
//____________________________________________________________________________
void NReWeight::AdoptWghtCalc(string name, NReWeightI* wcalc)
{
  if(!wcalc) return;

  fWghtCalc.insert(map<string, NReWeightI*>::value_type(name,wcalc));
}
//____________________________________________________________________________
NReWeightI* NReWeight::WghtCalc(string name)
{ 
  map<string, NReWeightI*>::iterator iter = fWghtCalc.find(name);
  if(iter != fWghtCalc.end()) return iter->second;
  
  return 0;
}
//____________________________________________________________________________
NSystSet & NReWeight::Systematics(void)
{ 
  return fSystSet; 
}
//____________________________________________________________________________
void NReWeight::Reconfigure(void)
{
#ifdef _NEUT_REWEIGHT_DEBUG_
  cout << "NReWeight: Reconfiguring ...";
#endif 

  vector<neut::rew::NSyst_t> svec = fSystSet.AllIncluded();

  map<string, NReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {

      NReWeightI * wcalc = it->second;

      vector<neut::rew::NSyst_t>::const_iterator parm_iter = svec.begin();
      for( ; parm_iter != svec.end(); ++parm_iter) {
          NSyst_t syst = *parm_iter;
          double val = fSystSet.Info(syst)->CurValue;
          wcalc->SetSystematic(syst, val);
      }//params

      wcalc->Reconfigure();

  }//weight calculators

#ifdef _NEUT_REWEIGHT_DEBUG_
  cout << "Done reconfiguring" << endl;;
#endif
}
//____________________________________________________________________________
double NReWeight::CalcWeight() 
{
// calculate weight for all tweaked physics parameters
//
  double weight = 1.0;
  fWeightXsec = 1;

  map<string, NReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {

    NReWeightI * wcalc = it->second;
    double w = wcalc->CalcWeight();

 #ifdef _NEUT_REWEIGHT_DEBUG_
    cout 
      << "Calculator: " << it->first << " => wght = " << w << endl;	
#endif

    // Store weight from neutrino cross section variations only, 
    // not FSI (casc) or nuclear effects (e.g. piless decay).
    // This assumes the mapping is named with "xsec" as in 
    // ${T2KREWEIGHT}/src/T2KNeutReWeight.cxx.
    // This is used by ${NIWG}/neut/ana/NInputXsec for reweighting
    // the input (bare) cross sections.
    if ( (it->first).find("xsec") != string::npos )
      fWeightXsec *=w;

    weight *= w;

  }
  return weight;
}
//____________________________________________________________________________
double NReWeight::CalcChisq(void) 
{
// calculate the sum of penalty terms for all tweaked physics parameters
//
  double tot_chisq = 0.0;

  map<string, NReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    NReWeightI * wcalc = it->second;
    double chisq = wcalc->CalcChisq(); 
#ifdef _NEUT_REWEIGHT_DEBUG_
    //cout << "Calculator: " << it->first << " => chisq = " << chisq;	
#endif
    tot_chisq *= chisq;
  }
  return tot_chisq;
}
//____________________________________________________________________________
void NReWeight::CleanUp(void)
{
  map<string, NReWeightI *>::iterator it = fWghtCalc.begin();
  for( ; it != fWghtCalc.end(); ++it) {
    NReWeightI * rw = it->second;
    if(rw) {
      delete rw;
      rw=0;
    }
  }
  fWghtCalc.clear();
}
//____________________________________________________________________________
void NReWeight::Print()
{
  vector<neut::rew::NSyst_t> syst_vec = this->Systematics().AllIncluded();
  int vec_size = syst_vec.size();

  cout << "NReWeight: Current set of systematic params:" << endl;;	
  for(int i = 0 ; i < vec_size ; i ++){
     cout 
        << " --o "  << NSyst::AsString(syst_vec[i])
        << " is set at " << this->Systematics().Info(syst_vec[i])->CurValue << endl;;
  }		       	        

  double chi2val = this->CalcChisq();

  cout << "Chisq_{penalty} = " << chi2val << endl;
}
//____________________________________________________________________________


