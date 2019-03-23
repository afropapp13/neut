# include "CrossSection.hh"
//#define VERBOSE

int CrossSection::SetParameters(const double energy, const double x, const double y) {
  
  double elepmin;
  double coslmax, coslmin;

  // lepton mass
  if (ITYPE == 12) {
    XM = 0.511;
  } else if (ITYPE == 14) {
    XM = 105.7;
  } else if (ITYPE == 16) {
    XM = 1784.1;
  }
  mass_lep = XM;

  // lepton momentum; minimum, maximum
  elepmin = mass_lep;
  if (energy < elepmin) return 1;
  elep = (energy - elepmin)*(x + 1.)/2. + elepmin;

  // lepton scattering angle; minimum, maximum
  coslmax =  1.;
  coslmin = -1.;
  coslep = (coslmax - coslmin)*(y + 1.)/2. + coslmin;

#ifdef VERBOSE
  printf("elep = %f, cosl = %f\n", elep, coslep);
#endif
  return 0;

}
