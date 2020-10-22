#ifndef __FAZEXP_H__
#define __FAZEXP_H__
//////////////////////////////////////////////////////////////
// Original comments start from here
//////////////////////////////////////////////////////////////
// This macro contains the calculations for z expansion axial
// form factors.
// Author. P. Stowell (05/08/2016)
//////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "neutmodelC.h"

using namespace std;

extern "C" {
double faZEXP(double x);
double fazexp_(double *x);
void zexpconfig_();
}

// Axial Coupling
static const double GA = 1.267;
static const int kmax = 15;
static double zexp_aterms[15];
static int zexp_nterms;
static double zexp_tc;
static double zexp_t0;

double GetZ(double q2) {

  /*
  std::cout <<  "Getting Z "
            << " TC = " << zexp_tc
            << " T0 = " << zexp_t0
            << " q2 = " << q2
            << std::endl;
  */

  // T Cut
  double num = sqrt(zexp_tc - q2) - sqrt(zexp_tc - zexp_t0);
  double den = sqrt(zexp_tc - q2) + sqrt(zexp_tc - zexp_t0);

  // cout << "Z = " << num/den << endl;
  return num / den;
}

double faZEXP(double x) {

  // Read Params
  double q2 = -fabs(x);

  // Calculate z
  double z = GetZ(q2);
  double FA = 0.0;

  int ncount = zexp_nterms;
  if (nemdls_.axzexpq4 > 0)
    ncount += 4;
  // cout << "Nterms = " << ncount << endl;

  // std::cout << "faZEXP: n: " << zexp_nterms << ", tc: " << zexp_tc
  //           << ", t0: " << zexp_t0 << ", t[0] " << zexp_aterms[0] << ", t[1] "
  //           << zexp_aterms[1] << ", t[2] " << zexp_aterms[2] << ", t[3] "
            // << zexp_aterms[3] << std::endl;

  for (int i = 0; i <= ncount; i++) {

    /*
    cout << i << " Adding " << pow(z,i) << " * "
         << zexp_aterms[i] << " = "
         << pow(z,i) * zexp_aterms[i] << endl;
    */

    FA += pow(z, i) * zexp_aterms[i];
  }

  return FA;
}

// Awful extra function that neut seems to need...
double fazexp_(double *x) { return faZEXP(*x); }

void PrintZExpTerms(bool showFA) {

  cout << " ZEXP State! " << endl;
  cout << " ------------------" << endl;
  cout << " T0 = " << zexp_t0 << endl;
  cout << " TC = " << zexp_tc << endl;

  int ncount = zexp_nterms;
  if (nemdls_.axzexpq4 > 0)
    ncount += 4;
  for (int i = 0; i <= ncount; i++) {
    cout << "ZEXP A" << i << " = " << zexp_aterms[i] << endl;
  }

  if (showFA) {
    cout << " FA Values " << endl;
    cout << " FAZ(0.00) = " << faZEXP(0.00) << endl;
    cout << " FAZ(0.25) = " << faZEXP(0.25) << endl;
    cout << " FAZ(0.50) = " << faZEXP(0.50) << endl;
    cout << " FAZ(0.75) = " << faZEXP(0.75) << endl;
    cout << " FAZ(1.00) = " << faZEXP(1.00) << endl;
    cout << " FAZ(1.50) = " << faZEXP(1.50) << endl;
    cout << " FAZ(2.00) = " << faZEXP(2.00) << endl;
    cout << " FAZ(3.00) = " << faZEXP(3.00) << endl;
  }
}

void ReadZExpBlock() {

  // Setup Params
  zexp_nterms = nemdls_.axzexpnt;

  zexp_tc = nemdls_.axzexptc;
  zexp_t0 = nemdls_.axzexpt0;

  zexp_aterms[0] = nemdls_.axzexpa0;
  zexp_aterms[1] = nemdls_.axzexpa1;
  zexp_aterms[2] = nemdls_.axzexpa2;
  zexp_aterms[3] = nemdls_.axzexpa3;
  zexp_aterms[4] = nemdls_.axzexpa4;
  zexp_aterms[5] = nemdls_.axzexpa5;
  zexp_aterms[6] = nemdls_.axzexpa6;
  zexp_aterms[7] = nemdls_.axzexpa7;
  zexp_aterms[8] = nemdls_.axzexpa8;
  zexp_aterms[9] = nemdls_.axzexpa9;
}

void UpdateZExpBlock() {

  nemdls_.axzexpnt = zexp_nterms;

  nemdls_.axzexptc = zexp_tc;
  nemdls_.axzexpt0 = zexp_t0;

  // Truncate Terms above Max
  int ncount = zexp_nterms;
  if (nemdls_.axzexpq4 > 0)
    ncount += 4;

  for (int i = ncount + 1; i < kmax; i++) {
    zexp_aterms[i] = 0.0;
  }

  // Fill Block of Coeff
  nemdls_.axzexpa0 = zexp_aterms[0];
  nemdls_.axzexpa1 = zexp_aterms[1];
  nemdls_.axzexpa2 = zexp_aterms[2];
  nemdls_.axzexpa3 = zexp_aterms[3];
  nemdls_.axzexpa4 = zexp_aterms[4];
  nemdls_.axzexpa5 = zexp_aterms[5];
  nemdls_.axzexpa6 = zexp_aterms[6];
  nemdls_.axzexpa7 = zexp_aterms[7];
  nemdls_.axzexpa8 = zexp_aterms[8];
  nemdls_.axzexpa9 = zexp_aterms[9];
}

void zexp_applysumrules() {
  //  cout << " ZEXP: Applying Sum Rules" << endl;
  //  PrintZExpTerms(false);

  // The Code below is from private correspondence
  // with Aaron Meyer on the calculation of sum Rules.
  // - P. Stowell

  // Gives the Q^-4 format at high Q^2
  double k0 = (double)zexp_nterms;
  double z0 = GetZ(0.0);

  double k1 = (double)zexp_nterms + 1;
  double z1 = pow(z0, (int)k1);

  double k2 = (double)zexp_nterms + 2;
  double z2 = pow(z0, (int)k2);

  double k3 = (double)zexp_nterms + 3;
  double z3 = pow(z0, (int)k3);

  double k4 = (double)zexp_nterms + 4;
  double z4 = pow(z0, (int)k4);

  // Get Delta (z shifts through terms)
  double del = 6.0 - 1.0 * k4 * k3 * k2 * z1 + 3.0 * k4 * k3 * z2 * k1 -
               3.0 * k4 * z3 * k2 * k1 + 1.0 * z4 * k3 * k2 * k1;

  // Setup Starting Parameters
  double b0 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;
  double b0z = 1.267;

  for (int ki = 1; ki <= zexp_nterms; ki++) {
    b0 += zexp_aterms[ki];
    b1 += ki * zexp_aterms[ki];
    b2 += ki * (ki - 1) * zexp_aterms[ki];
    b3 += ki * (ki - 1) * (ki - 2) * zexp_aterms[ki];

    b0z += zexp_aterms[ki] * pow(z0, ki);
  }

  // Setup A Terms ----
  // Copied Verbatim from ZExp PDF Aaron Sent

  /*
  // A0
  zexp_aterms[0] =
    (- 6*b0z - b0*(del-6.) + b3*(-z1 + 3.*z2 - 3*z3 + z4)
     + b2 * (3*(N+2)*z1 - 3*(3*N+5)*z2 + 3*(3*N+4)*z3 - 3*(N+1)*z4)
     + b1 * (-3*(N+3)*(N+2)*z1 + 3*(N+3)*(3N+4)*z2
             -3*(N+1)*(3N+8)*z3 + 3*(N+2)*(N+1)*z4) ) / (del);

  // AN+1
  zexp_aterms[zexp_nterms+1] = \
    (- (b0-b0z)*(N+4)*(N+3)*(N+2) \
     + b3*(1 - 0.5*(N+4)*(N+3)*z2 + (N+4)*(N+2)*z3 - 0.5*(N+3)*(N+2)*z4) \
     + b2*(-3*(N+2) + (N+4)*(N+3)*(N+2)*z2 \
           - (N+4)*(N+2)*(2N+3)*z3 + (N+3)*(N+2)*(N+1)*z4) \
     + b1*(3*(N+3)*(N+2) - 0.5*(N+4)*(N+3)*(N+3)*(N+2)*z2 \
           + (N+4)*(N+3)*(N+2)*(N+1)*z3 - 0.5*(N+3)*(N+2)*(N+2)*(N+1)*z4) \
     ) / (del) ;

  // A2
  zexp_aterms[zexp_nterms+2] =
    ( + 3.*(b0-b0z)*(N+4)*(N+3)*(N+1) \
      + b3*(-3. + 0.5*(N+4)*(N+3)*z1 - (3./2.)*(N+4)*(N+1)*z3 + (N+3)*(N+1)*z4)
      + b2*(3.*(3.*N+5) - (N+4)*(N+3)*(N+2)*z1 + 3*(N+4)*(N+1)*(N+1)*z3
            - (N+3)*(N+1)*(2N+1)*z4)
      + b1*(-3*(N+3)*(3N+4) + 0.5*(N+4)*(N+3)*(N+3)*(N+2)*z1
            -(3/2)*(N+4)*(N+3)*(N+1)*N*z3 + (N+3)*(N+2)*(N+1)*N*z4)
      ) / (del);

  // A3
  zexp_aterms[zexp_nterms+3] =
    (- 3*(b0-bz)*(N+4)*(N+2)*(N+1)
     + b3*(3 - (N+4)*(N+2)*z1 + (3/2)*(N+4)*(N+1)*z2 - 0.5*(N+2)*(N+1)*z4)
     + b2*(-3*(3N+4) + (N+4)*(N+2)*(2N+3)*z1
           - 3*(N+4)*(N+1)*(N+1)*z2 + (N+2)*(N+1)*N*z4)
     + b1*(3*(N+1)*(3N+8) - (N+4)*(N+3)*(N+2)*(N+1)*z1
           +(3/2)*(N+4)*(N+3)*(N+1)*N*z2 - (1/2)*(N+2)*(N+1)*(N+1)*N*z4)
     ) / (del);

  // A4
  zexp_aterms[zexp_nterms+4] =
    ( + (b0-b0z)*(N+3)*(N+2)*(N+1)
      + b3*(-1 + (1./2.)*(N+3)*(N+2)*z1 - (N+3)*(N+1)*z2 + 0.5*(N+2)*(N+1)*z3)
      + b2*(3*(N+1) - (N+3)*(N+2)*(N+1)*z1 + (N+3)*(N+1)*(2N+1)*z2
            - (N+2)*(N+1)*N*z3)
      + b1*(-3*(N+2)*(N+1) + 0.5*(N+3)*(N+2)*(N+2)*(N+1)*z1
            -(N+3)*(N+2)*(N+1)*N*z2 + 0.5*(N+2)*(N+1)*(N+1)*N*z3)
      ) / (del);
  */

  // Below is the same as above but with some replacements
  // to make life easier...

  // A0
  zexp_aterms[0] =
      (-6. * b0z - b0 * (del - 6.) + b3 * (-z1 + 3. * z2 - 3. * z3 + z4) +
       b2 * (3. * k2 * z1 - 3. * (3. * k0 + 5.) * z2 +
             3. * (3. * k0 + 4.) * z3 - 3. * k1 * z4) +
       b1 * (-3. * k3 * k2 * z1 + 3. * k3 * (3. * k0 + 4.) * z2 -
             3. * k1 * (3. * k0 + 8.) * z3 + 3. * k2 * k1 * z4)) /
      (del);

  // A1
  zexp_aterms[(int)k1] =
      (-(b0 - b0z) * k4 * k3 * k2 +
       b3 * (1. - 0.5 * k4 * k3 * z2 + k4 * k2 * z3 - 0.5 * k3 * k2 * z4) +
       b2 * (-3. * k2 + k4 * k3 * k2 * z2 - k4 * k2 * (2. * k0 + 3.) * z3 +
             k3 * k2 * k1 * z4) +
       b1 * (3. * k3 * k2 - 0.5 * k4 * k3 * k3 * k2 * z2 +
             k4 * k3 * k2 * k1 * z3 - 0.5 * k3 * k2 * k2 * k1 * z4)) /
      (del);

  // A2
  zexp_aterms[(int)k2] =
      (+3. * (b0 - b0z) * k4 * k3 * k1 +
       b3 * (-3. + 0.5 * k4 * k3 * z1 - (3. / 2.) * k4 * k1 * z3 +
             k3 * k1 * z4) +
       b2 * (3. * (3. * k0 + 5) - k4 * k3 * k2 * z1 + 3 * k4 * k1 * k1 * z3 -
             k3 * k1 * (2. * k0 + 1.) * z4) +
       b1 * (-3. * k3 * (3. * k0 + 4.) + 0.5 * k4 * k3 * k3 * k2 * z1 -
             (3. / 2.) * k4 * k3 * k1 * k0 * z3 + k3 * k2 * k1 * k0 * z4)) /
      (del);

  // A3
  zexp_aterms[(int)k3] =
      (-3. * (b0 - b0z) * k4 * k2 * k1 +
       b3 *
           (3. - k4 * k2 * z1 + (3. / 2.) * k4 * k1 * z2 - 0.5 * k2 * k1 * z4) +
       b2 * (-3. * (3. * k0 + 4.) + k4 * k2 * (2. * k0 + 3.) * z1 -
             3. * k4 * k1 * k1 * z2 + k2 * k1 * k0 * z4) +
       b1 * (3. * k1 * (3. * k0 + 8.) - k4 * k3 * k2 * k1 * z1 +
             (3. / 2.) * k4 * k3 * k1 * k0 * z2 -
             (1. / 2.) * k2 * k1 * k1 * k0 * z4)) /
      (del);

  // A4
  zexp_aterms[(int)k4] =
      (+(b0 - b0z) * k3 * k2 * k1 +
       b3 * (-1. + (1. / 2.) * k3 * k2 * z1 - k3 * k1 * z2 +
             0.5 * k2 * k1 * z3) +
       b2 * (3. * k1 - k3 * k2 * k1 * z1 + k3 * k1 * (2. * k0 + 1.) * z2 -
             k2 * k1 * k0 * z3) +
       b1 * (-3. * k2 * k1 + 0.5 * k3 * k2 * k2 * k1 * z1 -
             k3 * k2 * k1 * k0 * z2 + 0.5 * k2 * k1 * k1 * k0 * z3)) /
      (del);

  // Set NTerms to be higher now and update blocks
  UpdateZExpBlock();

  // PRINT
  //  cout << " Updated Terms! " << endl;
  //  PrintZExpTerms(true);
  //  sleep(10);

  return;
}

// Fix A0
void zexp_applyq0limit() {
  //  cout << " ZEXP: Fixing A0 Limits" << endl;
  //  PrintZExpTerms(false);

  double z = GetZ(0.0);
  double FA = 0.0;

  for (int i = 1; i <= zexp_nterms; i++) {
    FA = FA + pow(z, i) * zexp_aterms[i];
  }
  zexp_aterms[0] = -1.267 - FA;

  // Update
  UpdateZExpBlock();

  // PRINT
  //  PrintZExpTerms(true);
  //  sleep(10);

  return;
}

// Fix Z expansion Co-efficients
void zexpconfig_() {

  // Read Info From NEMDLS
  ReadZExpBlock();

  // Truncate at 10 terms
  if (zexp_nterms > 10)
    zexp_nterms = 10;

  // Sort Parameter Fixing
  if (nemdls_.axzexpq4) {
    zexp_applysumrules();
  } else {
    zexp_applyq0limit();
  }

  return;
}

#endif
