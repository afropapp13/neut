# include "CrossSection.hh"
# include "clenshaw_curtis.hh"

//#define DEBUG

///////////////////////////////////////////////////
CrossSection::CrossSection(int nu_type, int current)
  : ITYPE(nu_type), CCNC(current)
{
#ifdef DEBUG
  cout << "Neutrino type = " << ITYPE << ", Interaction = " << CCNC << endl;
#endif
  Initialize();
}
CrossSection::~CrossSection(){}


///////////////////////////////////////////////////
double CrossSection::Coefficient(const double energy) {

  double bte, coeff;

  bte = pow(RAD0NU,2.)*pow(ATOMN,2./3.)*pow(energy,2.)/(3.*pow(HBARC,2.));

  coeff = PRMASS*pow(ATOMN*PIMASS*GFERMI,2.)*pow(energy,3.)*pow(FPIFAC,2.)
    /pow(PI,2.)/pow(HBARC,8.)*exp(bte*pow(PIMASS/energy,2.))*pow(FMTCTM,2.);

  return coeff;
}


///////////////////////////////////////////////////
double CrossSection::DifferentialCrossSection(const double energy) {

  double dsigdt;
  double qsmx,qsmn,xmax,xmin,y2mpe2,bte,soln,soln_tmp,epitot,epikin;
  double arg1,arg2,arg3,femto1,femto2;
  double ymin,ymax,q2min,q2max,q2,correction;
  double r, bess;

  if(y < PIMASS/energy) return 0.;
  if(y > 1./(1.+PRMASS*x/(2.*energy))) return 0.;

  qsmx = 4.*energy*energy*(1.-y);
  qsmn = 0.;
  if(qsmx < 0.) return 0.;

  xmax = qsmx/(2.*PRMASS*y*energy);
  xmin = qsmn/(2.*PRMASS*y*energy);
  if (x > xmax) return 0.;
  if (x < xmin) return 0.;

  //--- Lepton mass effects
  ymin = PIMASS/energy;
  ymax = 1.-mass_lep/energy;

  q2min = pow(mass_lep,2.)*y/(1-y);
  q2max = 2.*PRMASS*energy*(1.-mass_lep/energy);
  q2 = 2.*PRMASS*energy*x*y;

  correction = 1.-q2min/(2.*(q2+pow(PIMASS,2.)))
    + y/4.*q2min*(q2-q2min)/pow(q2+pow(PIMASS,2.),2.);

  if (y < ymin || y > ymax || q2 < q2min || q2 > q2max) correction = 0.;
  //-----

  y2mpe2 = pow(y,2.) - pow(PIMASS/energy,2.);
  bte = pow(RAD0NU,2.)*pow(ATOMN,2./3.)*pow(energy,2.)/(3.*pow(HBARC,2.));

  soln_tmp = sqrt(abs(y*(y+2.*PRMASS*x/energy)*abs(y2mpe2)))
    *(1.-y)/pow(1.+2.*PRMASS*energy*x*y/pow(AM,2.),2.);
  epitot = energy*y;
  epikin = energy*y-PIMASS;
  if(epikin < 0.) return 0.;

  arg1 = -2.*bte*y*(y+PRMASS*x/energy-z*(1+PRMASS*x/energy)*sqrt(abs(y2mpe2)));
  arg2 = -2.*bte*sqrt(abs(((1.-pow(z,2.))*y2mpe2*(2.*PRMASS/energy)*x*y
			   *(1.-y*(1.+PRMASS*x/(2.*energy))))));

  femto1 = crosin_(epikin);
  arg3 = -9.*pow(ATOMN,1./3.)*femto1/(16.*PI*pow(RAD0NU,2.));

  soln = soln_tmp*exp(arg1+arg3);

  femto2 = crosto_(epikin);
  r = getr_(epikin,epitot,femto2);
  dsigdt = P2OS*pow(femto2,2.)*(1.+pow(r,2.))/(4.*PI);

  if (arg2 < -50.) arg2 = -50.;
  bess = dbesi0_(arg2);

  if (CCNC != 1 || LeptonMass != 1) correction = 1.;
  return ccnc_ratio*soln*bess*dsigdt*correction;

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
  double X[3];

  dim_num = 3;

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
    if (SetParameters(energy, X[0], X[1], X[2]) != 0) {
      weight[order] = 0.;
    } else {
      weight[order] *= DifferentialCrossSection(energy)
	*(1./(1.+PRMASS*x/(2.*energy))-PIMASS/energy)*2.; // 1.-(-1.) = 2.
    }
  }

  weight_sum = r8vec_sum ( order_nd, weight ); // Sigma_{i}^{order_nd}weight[i]

  delete [] order_1d;
  delete [] point;
  delete [] weight;

  return weight_sum*Coefficient(energy)*1.E+40/8.; // /8 = [-1,1]^3

}
