/********************************************************************************************
*                                                                                           *
*     ( purpose )                                                                           *
*       Set constants for MK-model                                                          *
*                                                                                           *
*     ( creation date and author )                                                          *
*     Minoo Kabirnezhad, March 2019                                                         *
*                                                                                           *
*     ( comment )                                                                           *
*       Reference: M. Kabirnezhad, Phys.Rev.D 97 (2018) no.1, 013002                        *
*                  M. Kabirnezhad thesis:                                                   *
*                  https://www.ncbj.gov.pl/sites/default/files/m.kabirnezhad_thesis_0.pdf   *
*                  Rein and Sehgal,Ann. of Phys. 133(1981),79-153                           *
*                  D. Rein ,Z.Phys.C 35(1987),43-64                                         *
*                  S. L. Adler, Ann. Phys. (N.Y.) 50, 189 (1968).                           *
********************************************************************************************/
#ifndef __MKCONS_IS_DEFINED__
#define __MKCONS_IS_DEFINED__
// C++ includes
#include <iostream>
#include <cmath>

// NEUT include for parameters
#include "neutmodelC.h"

// An anonymous namespace for Minoo's shared constants for the mk_imodeXY.cc
// Implementation intentionally made to mimic Rein-Sehgal code in NEUT
// The parameters are now read from the NEUT card file instead of this namespace (e.g. MARES, CA5(0))

namespace {

  // Make these the same as in rscons.h to avoid running into numerical issues
  // Mass of nucleon (GeV)
  const double M = 0.938919;
  //const double M= 0.9389186795; //GeV
  // Mass of pion (GeV)
  const double m_pi = 0.138;
  //const double m_pi= 0.13804;//GeV

  // ******************
  // Variables of model which are tunable
  // Mass of resonance background, set in card file now
  // MAR
  // CA5(0) axial form factor: now set with NEUT card file for each imode
  // CA5[0]
  // Axial mass: now set with NEUT card file for each imode
  // M_A
  // Vector mass: now set with NEUT card file for each imode
  // M_V
  // ******************

  // Weinberg angle
  const double sin2Wein=0.231;

  double MAR = 1.06;
  double M_A = 1.03;
  double M_V=0.84;

  // Mass of rho
  const double m_rho= 0.7758;
  const double F0_A= 1.;
  const double F0_1=1.0;
  const double F0_pi= 1.0;
  const double F0_2=1.85;
  const double F0_rho= 1.;
  const double f_pi=0.093;
  const double g_A= 1.26;
  const double g_NN=13.50;//g_NN\pi
  const double G_F= 1.1664e-5 ;//GeV-2
  const double Omega= 1.05 ;//GeV+2
  //const double Z= 3./4. ;
  const double Z= 0.76;
  //const double pi = 3.141592653589793;
  const double pi = 3.1415926;
  const double Gvcm= 5.06;

  // Number of resonances
  const int nRes = 17;

  // Resonance masses
  const double MR[nRes] = {1.232, 1.430, 1.515, 1.535, 1.63, 1.655, 1.675, 1.685, 1.7, 1.7, 1.71, 1.72, 1.88, 1.89, 1.92, 1.93, 1.6};

  // Resonance widths
  const double RWA[nRes] = {0.117, 0.35, 0.115, 0.15, 0.14, 0.14, 0.15, 0.13, 0.15, 0.3, 0.1, 0.25, 0.33, 0.28, 0.26, 0.285, 0.32};

  // Resonance braching ratios
  const double BRA[nRes] = {1., 0.65, 0.6, 0.45, 0.25, 0.7, 0.4, 0.67, 0.12, 0.15, 0.12, 0.11, 0.12, 0.22, 0.12, 0.4, 0.18};

  // Resonance normalizations
  const double NA[nRes] = {1.01762, 0.841642, 0.990136, 0.936202, 0.947601, 0.949313, 0.98815, 1.04593, 0.98157, 1.15457, 0.858881, 0.927106 ,1.24422, 0.915585, 0.903074, 1.17743, 1.};  

  // Rein-Sehgal signs
  const int Dsgn[nRes] = {1, 1,-1,1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1};
  const int Cjsgn_plus[nRes]  ={1,-1,-1, 1, 1, 1, 1,-1, -1,-1,-1, 1, -1, -1, 1, 1,1};
  const int Cjsgn_minus[nRes] ={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1};

  const int Jsgn[nRes] ={1,-1,-1, 1, 1, 1, 1,-1, -1,-1,-1, 1, -1, -1, 1, 1,1};

  // 2l+1
  const int LP[nRes] = {3, 3, 5, 1, 1, 1, 5, 7, 5, 5, 3, 3, 7, 3, 3, 7, 3};
  // 2j+1/2
  const int JP[nRes] = {2, 1, 2, 1, 1, 1, 3, 3, 2, 2, 1, 2, 3, 1, 2, 4, 2};

  const int Del[nRes] ={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ,1 ,1 ,1, 1, 1};

  // Different CA5_0 form factor for different resonances
  double C50[nRes] = {1.02, 1., 0.5, 1.53, 1., 1., 0.6, 1.3 ,1.0, 1., 1.2,   1.1 , 0.6, 0.6, 0.6, 1.4,1.}; 
  const double CV[nRes] = {1.12, 1., 1., 1., 1., 1., 1. , 1., 1., 1., 1. ,1. ,1. ,1., 1.,1.,1.};
  const double qA[nRes] ={2.97, 0.64,  0.,  0.93, 0., 0., 0.,  0.,  0.0, 0., 0.,  0. , 0. ,0. ,0., 0., 0.};
  const double qV[nRes] ={-2.86, 0.,  2.498,  0., 0., pi, 0.,  0.,  0.0, 0., 0.,  0. , 0. ,0. ,0., 0., 0.};

  inline double sqrt_chk(double number) {
    if (number >= 0) return sqrt(number);
    if (number + 1E-5 < 0) {
      std::cerr << "Had to add more than 1E-5 to make positive" << std::endl;
      std::cerr << "Number is " << number << std::endl;
      std::cerr << "Returning sqrt(" << number << ")=0" << std::endl;
    }
    return 0;
  }
}


#endif
