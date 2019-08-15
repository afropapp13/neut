#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "necardC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"
#include "rsconsC.h"

#ifdef __cplusplus
extern "C" {
double elaxff_(double *Q2_GeV, bool *nc, int *sgn);
}
#else
double elaxff_(double *Q2_GeV, _Bool *nc, int *sgn);
#endif

double qedifcrs_(int *ip, int *iscc, float *energy_GeV, float *elep_GeV,
                 float *coslep) {
  /* All in GeV */

  //--- Parameters

  double const HC = 0.197326; // GeV fm
  double const GF = 1.1664e-11;
  double const COSB2 = .9494;
  double const XMNc = XMN; // p/n average
  double const DMN = 0.;
  double const EBF = 0.;

  double const TRCORA = 6.0;
  double const TRCORB = 0.34;

  int mec = 0;

  float Maxial = nemdls_.xmaqe;
  float Mvec = nemdls_.xmvqe;
  float kappa = nemdls_.kapp;

  float fperr = nemdls_.fpqe;

  float PFermi = nenupr_.pfsurf;
  float Ebind = -nenupr_.vnuini;

  float Nnum = neuttarget_.numbndn;
  float Znum = neuttarget_.numbndp;

  float sccv = nemdls_.sccfv;
  float scca = nemdls_.sccfa;

  //--- Local variables
  double W_GeV, Q_GeV, Q2_GeV2, PLEP, WEF, QVEC2, AP, BP, Q2EF;
  double FGE, FGM, F1, F2, FA, FP, FV3, FA3;
  double GPE, GPM, GNE, GNM;
  double T1, T2, TA, TB, T8;
  double W1, W2, WA, WB, W8;
  double A1, A2, A3, A4, A5, A6, A7;
  double B0, B1, B2;
  double C, D, X0, DBLS, S1, S2, EU, EL;
  double DBLQ2;
  double COX, DSIG;
  double PFn;

  int nupdgsign;

  double Mlep_GeV = 0;

  DSIG = 0.;
  EL = 0.;

#ifdef __cplusplus
  bool isnc;
#else
  _Bool isnc;
#endif

  if (*iscc) {
    isnc = 0;
    if (abs(*ip) == 12) {
      Mlep_GeV = XME;
    } else if (abs(*ip) == 14) {
      Mlep_GeV = XMMU;
    } else if (abs(*ip) == 16) {
      Mlep_GeV = XMTAU;
    }
  } else {
    isnc = 1;
    Mlep_GeV = 0.;
  }

  if (*ip > 0) {
    nupdgsign = 1;
  } else {
    nupdgsign = -1;
  }

  if ((*elep_GeV) < Mlep_GeV) {
    return 0;
  }

  //--- Kinematical variables
  W_GeV = *energy_GeV - (*elep_GeV);
  if (W_GeV <= 0) {
    return 0;
  }
  PLEP = sqrt(((*elep_GeV)) * ((*elep_GeV)) - Mlep_GeV * Mlep_GeV);
  EU = sqrt((PFermi) * (PFermi) + XMNc * XMNc);
  WEF = W_GeV + EBF - (Ebind);
  AP = (Ebind) * (1. + EBF / W_GeV);
  BP = EBF * (1. - (Ebind) / W_GeV);
  DBLQ2 = ((*energy_GeV) * (*energy_GeV) + PLEP * PLEP -
           2. * (*energy_GeV) * PLEP * (*coslep));
  if (DBLQ2 < 0. && DBLQ2 >= -1.e-7)
    DBLQ2 = 0.;
  Q_GeV = sqrt(DBLQ2);
  QVEC2 = Q_GeV * Q_GeV;
  Q2_GeV2 = Q_GeV * Q_GeV - W_GeV * W_GeV;
  Q2EF = Q_GeV * Q_GeV - WEF * WEF - DMN;

  //--- Form factor
  if ((nemdls_.mdlqe % 10) == 2) {
    // printf("qedifcrs.c Q2_GeV2 = %g\n",Q2_GeV2);
    bbba05_(&Q2_GeV2, &GPE, &GPM, &GNE, &GNM);
    // printf("qedifcrs.c GPE GPM GNE GNM = %g %g %g %g\n",GPE   , GPM  ,  GNE,
    // GNM);
    FGE = GPE - GNE;
    FGM = GPM - GNM;
  } else if ((nemdls_.mdlqe % 10) == 3) {
    // printf("qedifcrs.c Q2_GeV2 = %g\n",Q2_GeV2);
    bbba07_(&Q2_GeV2, &GPE, &GPM, &GNE, &GNM);
    // printf("qedifcrs.c GPE GPM GNE GNM = %g %g %g %g\n",GPE   , GPM  ,  GNE,
    // GNM);
    FGE = GPE - GNE;
    FGM = GPM - GNM;
  } else if ((nemdls_.mdlqe % 10) == 1) {
    FGE = pow(1. + Q2_GeV2 / ((Mvec) * (Mvec)), -2.);
    FGM = (1. + 3.71) * FGE;
  } else {
    printf("qedifcrs: Unknown MDLQE = %d\n", nemdls_.mdlqe);
  }

  // Meson Exchange Current, also in dnelsq2.F
  mec = ((nemdls_.mdlqe) / 100);
  // printf("mec = %d\n",mec);
  if (mec == 1)
    FGM = FGM * sqrt(1 + TRCORA * Q2_GeV2 * exp(-1 * Q2_GeV2 / TRCORB));

  F1 = (FGE + Q2_GeV2 * FGM / (XMNc * XMNc) / 4.) /
       (1. + Q2_GeV2 / (XMNc * XMNc) / 4.);
  F2 = 0.5 * (FGE - FGM) / XMNc / (1. + Q2_GeV2 / (XMNc * XMNc) / 4.);

  FA = elaxff_(&Q2_GeV2, &isnc, &nupdgsign);

  // printf("qedifcrs: FA = %f\n",FA);

  // FP = 2.*XMNc*FA/(Q2_GeV2 + 139.57*139.57);
  FP = fperr * 2. * XMNc * FA / (Q2_GeV2 + 0.13957 * 0.13957);

  // Second-class currents
  FV3 = sccv * (-1.) * pow(1. + Q2_GeV2 / ((Mvec) * (Mvec)), -2.);
  FA3 = scca * FA;

  //--- T parameters,still need to add SCC here
  T1 = 0.5 * Q2_GeV2 * pow(F1 - 2. * XMNc * F2, 2.) +
       (2. * XMNc * XMNc + .5 * Q2_GeV2) * FA * FA;
  T2 = 2. * XMNc * XMNc * (F1 * F1 + Q2_GeV2 * (F2 * F2 + FA3 * FA3) + FA * FA);
  //  TA = XMNc*XMNc*(2.*XMNc*F1*F2 + (.5*Q2_GeV2 - 2.*XMNc*XMNc)*F2*F2
  //				-2.*XMNc*FA*FP + .5*Q2_GeV2*FP*FP);
  // add second class currents
  // TA = XMNc*XMNc*(2.*XMNc*F1*F2 + (.5*Q2_GeV2 - 2.*XMNc*XMNc)*(F2*F2 -
  // FV3*FV3) 				-2.*XMNc*FA*FP + .5*Q2_GeV2*FP*FP
  //              + FA3*(2*XMNc*FA - Q2_GeV2*FP) + FV3*(2*XMNc*F1 +
  //              Q2_GeV2*F2));
  TA = XMNc * XMNc *
       (2 * XMNc * F1 * F2 + (.5 * Q2_GeV2 - 2 * XMNc * XMNc) * F2 * F2 +
        FV3 * FV3 * (0.5 * Q2_GeV2 + 2 * XMNc * XMNc) - 2. * XMNc * FA * FP +
        0.5 * Q2_GeV2 * FP * FP + FA3 * (2 * XMNc * FA - Q2_GeV2 * FP) +
        FV3 * (2 * XMNc * F1 + Q2_GeV2 * F2) + 0.5 * Q2_GeV2 * FA3 * FA3);
  TB = -.5 * T2 - XMNc * XMNc * FV3 * (2 * XMNc * F1 + F2 * Q2_GeV2) -
       XMNc * XMNc * FA3 * (2 * XMNc * FA - Q2_GeV2 * FP);
  T8 = 2. * XMNc * XMNc * FA * (F1 - 2. * XMNc * F2);

  //--- B parameters
  //  if (Q_GeV == 0.){
  //    printf("returning 0 because Q_GeV == 0.");
  //    return 0.;
  //  }

  C = -WEF / Q_GeV;
  D = Q2EF / (2. * Q_GeV * XMNc);
  PFn = pow(2. * (Nnum) / ((Znum) + (Nnum)), 1. / 3.) * (PFermi);
  X0 = 3. * (Nnum) / (4. * Q_GeV * pow(PFn, 3.));
  if (WEF > 0.) {
    DBLS = 1. - C * C + D * D;
    if (DBLS <= 0.)
      DBLS = 0.;
    S1 = XMNc * (C * D + sqrt(DBLS)) / (1. - C * C);
    S2 = kappa * (EU - WEF);
    EL = fmax(S1, S2);
  }

  //  if (WEF < 0. || EU < EL){
  //
  //    if (WEF < 0.)
  //      printf("returning 0 because WEF < 0\n");

  //    if ( EU < EL)
  //      printf("return 0 because EU < EL\n");
  //
  //    return 0.;

  //  }

  B0 = X0 * (EU - EL + AP * log((EU - (Ebind)) / (EL - (Ebind))) +
             BP * log((EU - (Ebind) + W_GeV) / (EL - (Ebind) + W_GeV)));
  B1 = X0 / XMNc *
       (.5 * (EU * EU - EL * EL) +
        AP * (EU - EL + (Ebind)*log((EU - (Ebind)) / (EL - (Ebind)))) +
        BP * ((EU - EL) + ((Ebind)-W_GeV) * log((EU - (Ebind) + W_GeV) /
                                                (EL - (Ebind) + W_GeV))));
  B2 = X0 / (XMNc * XMNc) *
       (1. / 3. * (pow(EU, 3.) - pow(EL, 3.)) +
        AP * (.5 * (EU * EU - EL * EL) + (Ebind) * (EU - EL) +
              (Ebind) * (Ebind)*log((EU - (Ebind)) / (EL - (Ebind)))) +
        BP * (.5 * (EU * EU - EL * EL) + ((Ebind)-W_GeV) * (EU - EL) +
              pow((Ebind)-W_GeV, 2.) *
                  log((EU - (Ebind) + W_GeV) / (EL - (Ebind) + W_GeV))));

  //--- A parameters
  A1 = B0;
  A2 = B2 - B0;
  A3 = C * C * B2 + 2. * C * D * B1 + D * D * B0;
  A4 = B2 - 2. * (Ebind) / XMNc * B1 + (Ebind) * (Ebind) / (XMNc * XMNc) * B0;
  A5 = C * B2 + (D - (Ebind)*C / XMNc) * B1 - (Ebind)*D / XMNc * B0;
  A6 = C * B1 + D * B0;
  A7 = B1 - (Ebind) / XMNc * B0;

  //--- W_GeV parameters
  W1 = A1 * T1 + .5 * (A2 - A3) * T2;
  W2 = (A4 + 2 * W_GeV / Q_GeV * A5 + W_GeV * W_GeV / (Q_GeV * Q_GeV) * A3 +
        .5 * Q2_GeV2 / (Q_GeV * Q_GeV) * (A2 - A3)) *
       T2;
  WA = 1. / (Q_GeV * Q_GeV) * (1.5 * A3 - .5 * A2) * T2 +
       1. / (XMNc * XMNc) * A1 * TA + 2. / (XMNc * Q_GeV) * A6 * TB;
  WB = 1. / XMNc * (A7 + W_GeV / Q_GeV * A6) * TB +
       W_GeV / (Q_GeV * Q_GeV) * (1.5 * A3 - .5 * A2 + Q_GeV / W_GeV * A5) * T2;
  W8 = nupdgsign / XMNc * (A7 + W_GeV / Q_GeV * A6) * T8;

  //--- Differential cross section
  COX = PLEP / (*elep_GeV) * (*coslep);
  DSIG = COSB2 * pow(GF * HC, 2.) * (*elep_GeV) * PLEP / PI * .5 *
         (W2 * (1. + COX) + (2. * W1 + Mlep_GeV * Mlep_GeV * WA) * (1. - COX) -
          2. * W8 * ((*energy_GeV) + (*elep_GeV)) * (1. - COX) +
          (WB + W8) * 2. * Mlep_GeV * Mlep_GeV / (*elep_GeV));

  return DSIG;
}
