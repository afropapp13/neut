#ifndef __FA3COMP_H__
#define __FA3COMP_H__
//////////////////////////////////////////////////////////////
// This macro contains functions for the 2 and 3 component
// axial form factors
// Author  P.Stowell (05/08/2016)
//////////////////////////////////////////////////////////////

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
#include "neutmodelC.h"

using namespace std;

extern "C"
{
  double fa2COMP(double x);
  double fa2comp_(double *x);
  double fa3COMP(double x);
  double fa3comp_(double *x);
}

// Axial Coupling
static const double GA = 1.267;

double fa2COMP(double x){

  // Read Params
  double q2 = fabs(x);
  double axial_ff_gamma = nemdls_.axffgamma;
  double axial_ff_alpha = nemdls_.axffalpha;  
  if (axial_ff_gamma < 0.0) axial_ff_gamma = 0.0;

  //  std::cout << "Alpha, Gamma, Beta, Theta = " << nemdls_.axffalpha << " " << nemdls_.axffgamma << std::endl;

  // Set Ma_Axl as axial meson
  double ma_axl = 1.230;

  // Gamma term
  double gterm = 1.0 / pow(1.0 + axial_ff_gamma * q2, 2);

  // Alpha Term
  double aterm = (1.0 - axial_ff_alpha + \
		  (axial_ff_alpha * (ma_axl * ma_axl) / (ma_axl*ma_axl + q2)));

  //  cout << " Returning 2 Comp " << -1.0 * GA * gterm * aterm << " " << gterm << " " << aterm << endl;
  double FA = -1.0 * GA * gterm * aterm;
  if (FA < 0.0) return FA;
  else return 0.0;
}

double fa3COMP(double x){

  // std::cout << "Alpha, Gamma, Beta, Theta = " << nemdls_.axffalpha << " " << nemdls_.axffgamma << " " << nemdls_.axffbeta << " " << nemdls_.axfftheta << std::endl;

  // Read Params
  double q2 = fabs(x);
  double axial_ff_theta = nemdls_.axfftheta;
  double axial_ff_beta  = nemdls_.axffbeta;
  if (axial_ff_beta < 0.0) axial_ff_beta = 0.0;

  // Get First 2Comp Term
  double comp2_term = fa2COMP(q2);

  // Get Exponential
  double thetaprime = sqrt( fabs(axial_ff_theta) * axial_ff_beta );
  if (axial_ff_theta < 0.0) thetaprime *= -1.0;

  double exp_term = -1.0 * GA * q2 * thetaprime * exp(- 1.0 * axial_ff_beta * q2);

  //  cout << " Returning 3 Comp " << q2 << " " << comp2_term + exp_term << endl;
  double FA = comp2_term + exp_term;
  if (FA < 0.0) return FA;
  else return 0.0;
}

// Awful extra function that neut seems to need...
double fa2comp_(double *x){
  return fa2COMP(*x);
}

// Awful extra function that neut seems to need...
double fa3comp_(double *x){
  return fa3COMP(*x);
}

#endif
