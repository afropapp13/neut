//____________________________________________________________________________ 
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 11, 2009 - CA
   Was first added in v2.5.1. Adapted from the T2K-specific version of the
   GENIE reweighting tool.
 @ Dec 17, 2010 - JD
   Added method to calculate weight for a modified formation zone. 
 @ Apr 6, 2011 - PD
   Implemented for NEUT
*/
//____________________________________________________________________________

#include "NReWeightUtils.h"

using namespace neut;
using namespace neut::rew;

//____________________________________________________________________________
int neut::utils::rew::Sign(double twkdial)
{
  if(twkdial < 0.) return -1;
  if(twkdial > 0.) return +1;
  return 0;
}
//____________________________________________________________________________

