# include "CrossSection.hh"
# include "clenshaw_curtis.cc"
# include "bbba07.cc"

//#define DEBUG

CrossSection::CrossSection(double Pfermi, int nu_type, int nu_sign)
  : PFermi(Pfermi), ITYPE(nu_type), NUSIG(nu_sign)
{
#ifdef DEBUG
  cout << "Neutrino type= " << ITYPE << " : " << NUSIG 
       << ", Fermi momentum= " << PFermi << endl;
#endif
  Initialize();
}
CrossSection::~CrossSection(){}

double CrossSection::DifferentialCrossSection(const double energy) {

  //--- Parameters
  double const PI = 3.1415926;
  double const HC = 197.326;
  double const GF = 1.1664e-11;
  double const COSB2 = .9494;
  double const XMN = 939.6;
  double const DMN = 0.;
  double const EBF = 0.;
  
  //--- Local variables
  double W,Q,Q2,PLEP,WEF,QVEC2,AP,BP,Q2EF;
  double FGE,FGM,F1,F2,FA,FP;
  double T1,T2,TA,TB,T8;
  double W1,W2,WA,WB,W8;
  double A1,A2,A3,A4,A5,A6,A7;
  double B0,B1,B2;
  double C,D,X0,DBLS,S1,S2,EU,EL;
  double DBLQ2;
  double COX, DSIG;
  double PFn;

  DSIG = 0.;
  EL = 0.;

  if (elep < XM) return 0.;

  //--- Kinematical variables
  W = energy - elep;
  if (W <= 0.) return 0.;
  PLEP = sqrt(elep*elep - XM*XM);
  EU = sqrt(PFermi*PFermi + XMN*XMN);
  WEF = W + EBF - Ebind;
  AP = Ebind*(1. + EBF/W);
  BP = EBF*(1. - Ebind/W);
  DBLQ2 = (energy*energy + PLEP*PLEP - 2.*energy*PLEP*coslep);
  if ( DBLQ2 < 0. && DBLQ2 >= -1.e-7) DBLQ2 = 0.;
  Q  = sqrt(DBLQ2);
  QVEC2 = Q*Q;
  Q2 = Q*Q - W*W;
  Q2EF = Q*Q - WEF*WEF - DMN;

  //--- Form factor
  FGE = pow(1. + Q2/7.1e+5, -2.);
  FGM = (1. + 3.71)*FGE;

  if (Trans_Corr == 1){
	FGM = FGM * sqrt( 1 + 6.0 * Q2 /1.e+6 * exp( -1.*Q2/1.e+6/.34 ) );
  }

  F1 = (FGE + Q2*FGM/(XMN*XMN)/4.) / (1. + Q2/(XMN*XMN)/4.);
  F2 = 0.5*(FGE-FGM)/XMN/(1.+Q2/(XMN*XMN)/4.);

  if (FF_model == 0) {
	//    FA = -1.232*pow(1.+Q2/(Maxial*Maxial), -2.);
	FA = -1.267*pow(1.+Q2/(Maxial*Maxial), -2.);
  } else if (FF_model == 1) {
    FA = faBBBA07(Q2/1.e+6);
  }else{
	exit(1);
  }
  FP = 2.*XMN*FA/(Q2 + 139.57*139.57);
      
  //--- T parameters
  T1 = 0.5*Q2*pow(F1- 2.*XMN*F2, 2.) + (2.*XMN*XMN + .5*Q2)*FA*FA;
  T2 = 2.*XMN*XMN*(F1*F1 + Q2*F2*F2 + FA*FA);
  TA = XMN*XMN*(2.*XMN*F1*F2 + (.5*Q2 - 2.*XMN*XMN)*F2*F2
                -2.*XMN*FA*FP + .5*Q2*FP*FP);
  TB = -.5*T2;
  T8 = 2.*XMN*XMN*FA*(F1- 2.*XMN*F2);
    
  //--- B parameters
  if (Q == 0.) return 0.;

  C = -WEF/Q;
  D = Q2EF/(2.*Q*XMN);
  PFn = pow(2.*Nnum/(Znum+Nnum), 1./3.)*PFermi;
  X0 = 3.*Nnum/(4.*Q*pow(PFn,3.));
  if (WEF > 0.) {
    DBLS = 1.-C*C+D*D;
    if (DBLS <= 0.) DBLS = 0.;
    S1 = XMN*(C*D+sqrt(DBLS))/(1.-C*C);
    S2 = kappa*(EU - WEF);
    EL = max(S1, S2);
  }

  if (WEF < 0. || EU < EL) return 0.;

  B0 = X0*(EU - EL+AP*log((EU-Ebind)/(EL-Ebind))
         +BP*log((EU-Ebind+W)/(EL-Ebind+W)));
  B1 = X0/XMN*(.5*(EU*EU-EL*EL)+
               AP*(EU-EL+Ebind*log((EU- Ebind)/(EL-Ebind)))
               +BP*((EU-EL)+(Ebind-W)*log((EU-Ebind+W)/(EL-Ebind+W))));
  B2 = X0/(XMN*XMN)*(1./3.*(pow(EU,3.)-pow(EL,3.))
                     +AP*(.5*(EU*EU-EL*EL)+Ebind*(EU-EL)
                          +Ebind*Ebind*log((EU-Ebind)/(EL-Ebind)))
                     +BP*(.5*(EU*EU-EL*EL)+(Ebind-W)*
                          (EU-EL)+pow(Ebind-W,2.)
                          *log((EU-Ebind+W)/(EL-Ebind+W))));

  //--- A parameters
  A1 = B0;
  A2 = B2 - B0;
  A3 = C*C*B2 + 2.*C*D*B1+D*D*B0;
  A4 = B2 - 2.*Ebind/XMN*B1+Ebind*Ebind/(XMN*XMN)*B0;
  A5 = C*B2+(D-Ebind*C/XMN)*B1-Ebind*D/XMN*B0;
  A6 = C*B1+D*B0;
  A7 = B1-Ebind/XMN*B0;

  //--- W parameters
  W1 = A1*T1+.5*(A2-A3)*T2;
  W2 = (A4+2*W/Q*A5+W*W/(Q*Q)*A3+.5*Q2/(Q*Q)
        *(A2-A3))*T2;
  WA = 1./(Q*Q)*(1.5*A3-.5*A2)*T2
    +1./(XMN*XMN)*A1*TA + 2./(XMN*Q)*A6*TB;
  WB = 1./XMN * (A7+W/Q*A6)*TB
    + W/(Q*Q)*(1.5*A3-.5*A2+Q/W*A5)*T2;
  W8 = NUSIG/XMN*(A7+W/Q*A6)*T8;

  //--- Differential cross section
  COX = PLEP/elep*coslep;
  DSIG = COSB2*pow(GF*HC,2.)*elep*PLEP/PI *.5 *
    (W2*(1.+COX)+(2.*W1+XM*XM*WA)*(1.-COX)
     -2.*W8*(energy+elep)*(1.-COX) 
     +(WB+W8)*2.*XM*XM/elep);

  return DSIG;

}

///////////////////////////////////////////////////
double CrossSection::TotalCrossSection(const double energy, const int nloop) {

  int dim;
  int dim_num;
  int order;
  int *order_1d;
  int order_nd;
  double *point;
  double *weight;
  double weight_sum;
  double X[2];

  dim_num = 2;

  order_1d = new int[dim_num];
  for ( dim = 0; dim < dim_num; dim++ ) {
    order_1d[dim] = nloop;
  }
  
  order_nd = i4vec_product ( dim_num, order_1d ); // order_1d[0]*order_1d[1]
  
  point  = new double[dim_num*order_nd];
  weight = new double[order_nd];
  
  clenshaw_curtis_compute_nd ( dim_num, order_1d, point, weight );

  for ( order = 0; order < order_nd; order++ ) {
    for ( dim = 0; dim < dim_num; dim++ ) {
      X[dim] = point[dim+order*dim_num];
    }
    if (SetParameters(energy, X[0], X[1]) != 0) {
      weight[order] = 0.;
    } else {
      weight[order] *= DifferentialCrossSection(energy)*(energy-mass_lep)*2.; // 1.-(-1.) = 2.
    }
  }

  weight_sum = r8vec_sum ( order_nd, weight ); // Sigma_{i}^{order_nd}weight[i]

  delete [] order_1d;
  delete [] point;
  delete [] weight;
  
  return weight_sum*1.e+12/8./4.; // /#of nucleon = 8/[-1,1]^2 = 4

}
