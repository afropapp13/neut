# include "CrossSection.hh"
//#define VERBOSE

int CrossSection::SetParameters(const double energy,
				const double x1, const double x2, const double x3) {
  
  double xhig, xlow;
  double yhig, ylow;
  double zhig, zlow;

  // z; minimum, maximum
  zhig =  1.;
  zlow = -1.;
  z = (zhig - zlow)*(x1 + 1.)/2. + zlow;

  // x; minimum, maximum
  xhig = 1.;
  xlow = 0.;
  x = (xhig - xlow)*(x2 + 1.)/2. + xlow;
  if (x == 0.) x = 0.0001;
  if (x == 1.) x = 0.9999;

  // y; minimum, maximum
  yhig = 1./(1.+PRMASS*x/(2.*energy));
  ylow = PIMASS/energy;
  y = (yhig - ylow)*(x3 + 1.)/2. + ylow;

#ifdef VERBOSE  
  printf("x = %f, y = %f, z = %f\n", x, y, z);
#endif
  return 0;

}
