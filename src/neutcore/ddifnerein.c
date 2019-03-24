#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "neutmodelC.h"
#include "neutparamsC.h"
#include "necardC.h"
#include "rsconsC.h"


double ddifnerein_(float *energy, int *ip, int *imode, float *x, float *y, float *t)
{


  const double M = XMP;//M is the proton mass
  const double pie = PI;//fundamental constant

  double Q2,Epi,El,ml,mpi;
  double plab,qlab,nu;
  double A1,A2,A3;
  double A1min,A2min,A3min;
  double A1max,A2max,A3max;
  double f1,f2,xsecdif;
  double tmin,tmax;
  double ymin,ymax;
  double ymine,ymaxe;
  double sTot,prop;

  double factor;
  float E,tt;

  float Ma = nemdls_.xmadif;
  float b = nemdls_.nucvoldif;//units of GeV^-2

  if (abs(*imode)==15) {//CC
    mpi = 0.13957;
    factor = 1.0;
  } else if (abs(*imode)==35) {//NC
    mpi = 0.1349766;
    factor = 0.5;//NC interactions have extra factor of 0.5
  } else {
    printf("ddifnerein: unknown interaction mode!\n");
    return 0.;
  }

  if (abs(*ip)==12){
    ml = XME;
  }else if (abs(*ip)==14){
    ml = XMMU;
  }else if (abs(*ip)==16){
    ml = XMTAU;
  }else{
    ml = 0.;
  }

  E = (*energy);//neutrino energy
  tt = (*t);

  const double G2 = (1.16637)*(1.16637)*pow(10.,-10.)*M/(pie*pie);//Actually, GF^2*M/pi^2; note GF is in GeV^-2, total units GeV^-3
  //pion forward scattering parameter squared below:
  const double fp2 = 0.93*0.93*0.13957*0.13957;//units GeV^2

  ymine = mpi/E;//this should be the absolute minimum (some E needed to create pion)
  ymaxe = (E-ml-mpi)/E;//this should be maximum, some energy needed to create pion & muon


  if (ymine>ymaxe) return 0.;
  if ((*y)>ymaxe) return 0.;
  if ((*y)<ymine) return 0.;
  if ((*y)<0.) return 0.;


  El = E*(1-(*y));//lepton energy
  Q2 = 2*M*E*(*x)*(*y);//= -q2
  nu = E*(*y);//E-El
  if (Q2<0.) return 0.;//sanity check
  if (Q2>2.*M*E) return 0.;//sanity check
  if (El<ml) return 0.;//sanity check
  if (nu<mpi) return 0.;//sanity chack
  if (fabs((El-(Q2+ml*ml)/(2.*E))/sqrt(El*El-ml*ml))>1) return 0.;

  prop = (Ma*Ma/(Ma*Ma+Q2))*(Ma*Ma/(Ma*Ma+Q2));//propagator

  //below are values for getting tmin & tmax, well, tmin & t|xsi=0
  //tmin is set by Rein paper, tmax is assumed from xsi being [-1,1]
  A1min = 1+(2*nu/M)+nu*nu/(M*M)-(nu*nu+Q2)/(M*M);
  A2min = (1+(nu/M))*(mpi*mpi-Q2-2*nu*nu)+2*nu*(nu*nu+Q2)/M;
  A3min = (mpi*mpi-Q2-2*nu*nu)*(mpi*mpi-Q2-2*nu*nu)
    -4*(nu*nu+Q2)*(nu*nu-mpi*mpi);

  A1max = 1+(2*nu/M)+nu*nu/(M*M);
  A2max = (1+(nu/M))*(mpi*mpi-Q2-2*nu*nu);
  A3max = (mpi*mpi-Q2-2*nu*nu)*(mpi*mpi-Q2-2*nu*nu);

  if (A2min*A2min-A1min*A3min<0) return 0.;//don't want imaginary terms
  if (A2max*A2max-A1max*A3max<0) return 0.;
  tmin = (A2min+sqrt(A2min*A2min-A1min*A3min))/A1min;
  tmax = (A2max+sqrt(A2max*A2max-A1max*A3max))/A1max;
  if (tmax>tmin) return 0.;//sanity check; remember, though, t is negative

  //is t within its absolute limits?
  if (fabs(tt)<fabs(tmin)) return 0;
  if (fabs(tt)>fabs(tmax)) return 0;


  Epi =  nu + tt/(2.*M);//from equation 8 in Rein's paper, for EpiLab, which I think is reasonable to use here
  //otherwise will need to change to Epi = E*y;
  if (Epi<mpi) return 0.;
  if (E-El-Epi<0) return 0.;//double check to make sure have at least a proton's rest mass of energy left
  //check to see from current set of variables if cos theta_qpi 
  //(the angle of the pion WRT the lepton scattering plane) is between -1 & 1
  //use the Epi from above
  if (fabs((tt-mpi*mpi+Q2+2*Epi*(E-El))/(2.*sqrt(Epi*Epi-mpi*mpi)*sqrt(E*E+sqrt(El*El-ml*ml)*sqrt(El*El-ml*ml)-2.*E*sqrt(El*El-ml*ml)*((El-(Q2+ml*ml)/(2.*E))/sqrt(El*El-ml*ml)))))>1.) return 0.;


  //pieN xsec term below, using Regge parameterization
  sTot = 24. + 12./sqrt(E*(*y));//in mb, need to divide by .389 to get GeV^-2

  //d3sigma/dx/dy/d|t| units = GeV^-3 * GeV^2 * GeV * mb * GeV-2 => mb/GeV^2
  xsecdif = factor*G2*fp2*E*(1-(*y))*prop*sTot*(sTot/0.389)*exp(-1.*b*fabs(tt))/(16.*pie);
  xsecdif *= pow(10.,11.);//factors of 10^-27 && 10^38 to get in 10^-38 cm^2

  return xsecdif;
}
