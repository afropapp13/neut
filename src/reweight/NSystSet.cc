//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.neut-mc.org
 or see $NEUT/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include <iostream>

#include "NSystSet.h"

using std::cout;
using std::endl;

using namespace neut;
using namespace neut::rew;

//_______________________________________________________________________________________
NSystSet::NSystSet()
{

}
//_______________________________________________________________________________________
NSystSet::NSystSet(const NSystSet & syst)
{
  this->Copy(syst);
}
//_______________________________________________________________________________________
NSystSet::~NSystSet()
{
  fSystematics.clear();
}
//_______________________________________________________________________________________
void NSystSet::Init(NSyst_t syst, double init, double min, double max, double step)
{
  if(syst == kNullSystematic) return;

  if(this->Added(syst)) {
    this->Remove(syst);    
  }

  NSystInfo * syst_info = new NSystInfo(init,min,max,step);
  fSystematics.insert( map<NSyst_t, NSystInfo*>::value_type(syst, syst_info) );     
}
//_______________________________________________________________________________________
void NSystSet::Remove(NSyst_t syst)
{
  fSystematics.erase(syst);
}
//_______________________________________________________________________________________
int NSystSet::Size(void) const
{
  return fSystematics.size();
}
//_______________________________________________________________________________________
bool NSystSet::Added(NSyst_t syst) const
{
  return (fSystematics.find(syst) != fSystematics.end());
}
//_______________________________________________________________________________________
vector<neut::rew::NSyst_t> NSystSet::AllIncluded(void)
{
  vector<NSyst_t> svec;

  map<NSyst_t, NSystInfo*>::const_iterator it = fSystematics.begin();
  for( ; it != fSystematics.end(); ++it) {
    NSyst_t syst = it->first;
    svec.push_back(syst);
  }
  return svec;
}
//_______________________________________________________________________________________
const NSystInfo * NSystSet::Info(NSyst_t syst) const
{
  if ( this->Added(syst) ) {
    return fSystematics.find(syst)->second;
  }
  return 0;
}
//_______________________________________________________________________________________
void NSystSet::Set(NSyst_t syst, double val)
{
  if ( this->Added(syst) ) {
    fSystematics[syst]->CurValue = val;
  }
  else {
    this->Init(syst);
    this->Set(syst,val);
  }
}
//_______________________________________________________________________________________
void NSystSet::Print(void)
{
  cout 
     << "Considering " << this->Size() << " systematics";
				    
  vector<neut::rew::NSyst_t> svec = this->AllIncluded();

  unsigned int i=0;
  vector<neut::rew::NSyst_t>::const_iterator it = svec.begin();
  for( ; it != svec.end(); ++it) {
    NSyst_t syst = *it;
    cout << "(" << i++ << ") : " << NSyst::AsString(syst);
  }
}
//_______________________________________________________________________________________
void NSystSet::Copy(const NSystSet & syst_set)
{
  return fSystematics.clear();

  map<NSyst_t, NSystInfo*>::const_iterator it = syst_set.fSystematics.begin();
  for( ; it != syst_set.fSystematics.end(); ++it) {
    NSyst_t     syst       = it->first;
    NSystInfo * syst_info  = it->second;

    double cur  = syst_info->CurValue;
    double init = syst_info->InitValue;
    double min  = syst_info->MinValue;
    double max  = syst_info->MaxValue;
    double step = syst_info->Step;

    this->Init(syst,init,min,max,step);
    this->Set(syst,cur);
  }
}
//_______________________________________________________________________________________
