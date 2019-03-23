#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "neutmodelC.h"
#include "neutparamsC.h"
#include "necardC.h"
#include "rsconsC.h"



/*  for 64 bits OS & g77 compiler -- M. Fechner*/
//#if defined(f2cFortran)&&!defined(gFortran)
//#define FUNCTION_RETURN double
//#else
//#define FUNCTION_RETURN float
//#endif

//extern "C" {
  //FUNCTION_RETURN bbba05_(double* Q2,double* GPE,double* GPM,double* GNE,double* GNM);
  double bbba05_(double* Q2,double* GPE,double* GPM,double* GNE,double* GNM);
  double bbba07_(double* Q2,double* GPE,double* GPM,double* GNE,double* GNM);
//}

/*
extern struct nemdls_common {
  int    mdlqe;
  int    mdlspi;
  int    mdldis;
  int    mdlcoh;
  int    mdlqeaf;
  float  xmaqe;
  float  xmaspi;
  float  xmvqe;
  float  xmvspi;
  float  kapp;
} nemdls_;

extern struct nenupr_common {
  float  pfsurf;
  float  pfmax;
  float  vnuini;
  float  vnufin;
  int    iformlen;
} nenupr_;

extern struct neuttarget_common {
  int    numbndn;
  int    numbndp;
  int    numfrep;
  int    numatom;
} neuttarget_;
*/

extern double faBBBA07(const double);

double
qedifcrs_(int *ip, int *ccnc, float *energy, float *elep, float *coslep)
{
  /* All in MeV */

  //--- Parameters

  //  double const PI = 3.1415926;
  //  double const HC = 197.326;  // MeV fm
  double const HC = 0.197326;   // GeV fm
  double const GF = 1.1664e-11;
  double const COSB2 = .9494;
  //  double const XMNc = XMN*1000;  // p/n average
  double const XMNc = XMN;  // p/n average
  //double XMNc;
  double const DMN = 0.;
  double const EBF = 0.;

  double const TRCORA = 6.0;
  double const TRCORB = 0.34;

  int mec = 0;

  //  float Maxial = nemdls_.xmaqe*1000;
  //  float Mvec   = nemdls_.xmvqe*1000;
  float Maxial = nemdls_.xmaqe;
  float Mvec   = nemdls_.xmvqe;
  float kappa  = nemdls_.kapp;

  int FF_model = nemdls_.mdlqeaf;

  float fperr = nemdls_.fpqe;

  //  float PFermi = nenupr_.pfsurf*1000;
  //  float Ebind  = nenupr_.vnuini*(-1000);
  float PFermi = nenupr_.pfsurf;
  float Ebind  = -nenupr_.vnuini;
  
  float Nnum   = neuttarget_.numbndn;
  float Znum   = neuttarget_.numbndp;
  
  float sccv = nemdls_.sccfv;
  float scca = nemdls_.sccfa;
  
  //--- Local variables
  double W,Q,Q2,PLEP,WEF,QVEC2,AP,BP,Q2EF;
  double FGE,FGM,F1,F2,FA,FP,FV3,FA3;
  double  GPE,GPM,GNE,GNM;
  double T1,T2,TA,TB,T8;
  double W1,W2,WA,WB,W8;
  double A1,A2,A3,A4,A5,A6,A7;
  double B0,B1,B2;
  double C,D,X0,DBLS,S1,S2,EU,EL;
  double DBLQ2;
  double COX, DSIG;
  double PFn;

  int    NUSIG;

  double XM=0;

  DSIG = 0.;
  EL = 0.;

  if (*ccnc != 0){
	if (abs(*ip)==12){
	  //	  XM = XME*1000;
	  XM = XME;
	}else if (abs(*ip)==14){
	  //	  XM = XMMU*1000;
	  XM = XMMU;
	}else if (abs(*ip)==16){
	  //	  XM = XMTAU*1000;
	  XM = XMTAU;
	}
  }else{
	XM = 0.;
  }

  if (*ip>0){
	NUSIG=1;
	//XMNc = XMNE;
  }else{
	NUSIG=-1;
	//XMNc = XMP;
  }

  if ((*elep) < XM) return 0.;

  //--- Kinematical variables
  W = *energy - (*elep);
  if (W <= 0.) return 0.;
  PLEP = sqrt(((*elep))*((*elep)) - XM*XM);
  EU = sqrt((PFermi)*(PFermi) + XMNc*XMNc);
  WEF = W + EBF - (Ebind);
  AP = (Ebind)*(1. + EBF/W);
  BP = EBF*(1. - (Ebind)/W);
  DBLQ2 = ((*energy)*(*energy) + PLEP*PLEP - 2.*(*energy)*PLEP*(*coslep));
  if ( DBLQ2 < 0. && DBLQ2 >= -1.e-7) DBLQ2 = 0.;
  Q  = sqrt(DBLQ2);
  QVEC2 = Q*Q;
  Q2 = Q*Q - W*W;
  Q2EF = Q*Q - WEF*WEF - DMN;

  //--- Form factor
  if ((nemdls_.mdlqe%10)==2) {
    //printf("qedifcrs.c Q2 = %g\n",Q2);
    bbba05_(&Q2,&GPE,&GPM,&GNE,&GNM);
    //printf("qedifcrs.c GPE GPM GNE GNM = %g %g %g %g\n",GPE   , GPM  ,  GNE, GNM);
    FGE=GPE-GNE;
    FGM=GPM-GNM;
  }
  else if ((nemdls_.mdlqe%10)==3) {
    //printf("qedifcrs.c Q2 = %g\n",Q2);
    bbba07_(&Q2,&GPE,&GPM,&GNE,&GNM);
    //printf("qedifcrs.c GPE GPM GNE GNM = %g %g %g %g\n",GPE   , GPM  ,  GNE, GNM);
    FGE=GPE-GNE;
    FGM=GPM-GNM;
  }
  else if ((nemdls_.mdlqe%10)==1) {
    FGE = pow(1. + Q2/((Mvec)*(Mvec)), -2.);
    FGM = (1. + 3.71)*FGE;
  } else {
    printf("qedifcrs: Unknown MDLQE = %d\n",nemdls_.mdlqe);
  }

  //Meson Exchange Current, also in dnelsq2.F
  mec = ((nemdls_.mdlqe)/100);
  // printf("mec = %d\n",mec);
  if (mec==1) FGM = FGM * sqrt(1+TRCORA * Q2 * exp( -1*Q2 / TRCORB ));

  F1 = (FGE + Q2*FGM/(XMNc*XMNc)/4.) / (1. + Q2/(XMNc*XMNc)/4.);
  F2 = 0.5*(FGE-FGM)/XMNc/(1.+Q2/(XMNc*XMNc)/4.);

  if (FF_model == 1) {
    printf("qedifcrs: axial dipole FF_model = %d\n",FF_model);
    // FA = -1.232*pow(1.+Q2/(Maxial*Maxial), -2.);
    FA = -1.267*pow(1.+Q2/((Maxial)*(Maxial)), -2.);
  } else if (FF_model == 2) {
    printf("qedifcrs: bbba07 BBBA07 FF_model = %d\n",FF_model);
    // FA = faBBBA07(Q2/1.e+6);
    FA = faBBBA07(Q2);
  }else{
    printf("qedifcrs: Unknown FF_model = %d\n",FF_model);
    exit(1);
  }
  //FP = 2.*XMNc*FA/(Q2 + 139.57*139.57);
  FP = fperr * 2.*XMNc*FA/(Q2 + 0.13957*0.13957); 

  //Second-class currents
  FV3 = sccv*(-1.)*pow(1. + Q2/((Mvec)*(Mvec)), -2.);
  FA3 = scca * FA ;


  //--- T parameters,still need to add SCC here
  T1 = 0.5*Q2*pow(F1- 2.*XMNc*F2, 2.) + (2.*XMNc*XMNc + .5*Q2)*FA*FA;
  T2 = 2.*XMNc*XMNc*(F1*F1 + Q2*(F2*F2 + FA3*FA3) + FA*FA);
//  TA = XMNc*XMNc*(2.*XMNc*F1*F2 + (.5*Q2 - 2.*XMNc*XMNc)*F2*F2 
//				-2.*XMNc*FA*FP + .5*Q2*FP*FP);
// add second class currents
  //TA = XMNc*XMNc*(2.*XMNc*F1*F2 + (.5*Q2 - 2.*XMNc*XMNc)*(F2*F2 - FV3*FV3)
  //				-2.*XMNc*FA*FP + .5*Q2*FP*FP
  //              + FA3*(2*XMNc*FA - Q2*FP) + FV3*(2*XMNc*F1 + Q2*F2));
  TA = XMNc*XMNc*(2*XMNc*F1*F2 + (.5*Q2-2*XMNc*XMNc)*F2*F2 
		  + FV3*FV3*(0.5*Q2+2*XMNc*XMNc) -2.*XMNc*FA*FP + 0.5*Q2*FP*FP
              + FA3*(2*XMNc*FA - Q2*FP) + FV3*(2*XMNc*F1 + Q2*F2)
		  + 0.5*Q2*FA3*FA3);
  TB = -.5*T2 - XMNc*XMNc*FV3*(2*XMNc*F1 + F2*Q2) - XMNc*XMNc*FA3*(2*XMNc*FA - Q2*FP);
  T8 = 2.*XMNc*XMNc*FA*(F1- 2.*XMNc*F2);

  //--- B parameters
  if (Q == 0.) return 0.;

  C = -WEF/Q;
  D = Q2EF/(2.*Q*XMNc);
  PFn = pow(2.*(Nnum)/((Znum)+(Nnum)), 1./3.)*(PFermi);
  X0 = 3.*(Nnum)/(4.*Q*pow(PFn,3.));
  if (WEF > 0.) {
	DBLS = 1.-C*C+D*D;
	if (DBLS <= 0.) DBLS = 0.;
	S1 = XMNc*(C*D+sqrt(DBLS))/(1.-C*C);
	S2 = kappa*(EU - WEF);
	EL = fmax(S1, S2);
  }

  if (WEF < 0. || EU < EL) return 0.;

  B0 = X0*(EU - EL+AP*log((EU-(Ebind))/(EL-(Ebind)))
		   +BP*log((EU-(Ebind)+W)/(EL-(Ebind)+W)));
  B1 = X0/XMNc*(.5*(EU*EU-EL*EL)+
			   AP*(EU-EL+(Ebind)*log((EU- (Ebind))/(EL-(Ebind))))
			   +BP*((EU-EL)
					+((Ebind)-W)*log((EU-(Ebind)+W)/(EL-(Ebind)+W))));
  B2 = X0/(XMNc*XMNc)*(1./3.*(pow(EU,3.)-pow(EL,3.))
					 +AP*(.5*(EU*EU-EL*EL)+(Ebind)*(EU-EL)
						  +(Ebind)*(Ebind)*log((EU-(Ebind))/(EL-(Ebind))))
					 +BP*(.5*(EU*EU-EL*EL)+((Ebind)-W)*
						  (EU-EL)+pow((Ebind)-W,2.)
						  *log((EU-(Ebind)+W)/(EL-(Ebind)+W))));
  
  //--- A parameters
  A1 = B0;
  A2 = B2 - B0;
  A3 = C*C*B2 + 2.*C*D*B1+D*D*B0;
  A4 = B2 - 2.*(Ebind)/XMNc*B1+(Ebind)*(Ebind)/(XMNc*XMNc)*B0;
  A5 = C*B2+(D-(Ebind)*C/XMNc)*B1-(Ebind)*D/XMNc*B0;
  A6 = C*B1+D*B0;
  A7 = B1-(Ebind)/XMNc*B0;

  //--- W parameters
  W1 = A1*T1+.5*(A2-A3)*T2;
  W2 = (A4+2*W/Q*A5+W*W/(Q*Q)*A3+.5*Q2/(Q*Q)
		*(A2-A3))*T2;
  WA = 1./(Q*Q)*(1.5*A3-.5*A2)*T2
	+1./(XMNc*XMNc)*A1*TA + 2./(XMNc*Q)*A6*TB;
  WB = 1./XMNc * (A7+W/Q*A6)*TB
	+ W/(Q*Q)*(1.5*A3-.5*A2+Q/W*A5)*T2;
  W8 = NUSIG/XMNc*(A7+W/Q*A6)*T8;

  //--- Differential cross section
  COX = PLEP/(*elep)*(*coslep);
  DSIG = COSB2*pow(GF*HC,2.)*(*elep)*PLEP/PI *.5 *
	(W2*(1.+COX)+(2.*W1+XM*XM*WA)*(1.-COX)
	 -2.*W8*((*energy)+(*elep))*(1.-COX)
	 +(WB+W8)*2.*XM*XM/(*elep));

  return DSIG;

}
