# include <stdlib.h>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

class CrossSection {

public:

  CrossSection(double Pfermi, int nu_type, int nu_sign, double XAmass);
  ~CrossSection();

  int    SetParameters(const double energy, const double x, const double y);
  double TotalCrossSection(const double energy, const int nloop);
  double DifferentialCrossSection(const double energy);
  void   Initialize(double xamass);

private:

  double mass_lep, PFermi, Ebind;
  double elep, coslep, XMaxial;
  double Nnum, Znum;
  double XM;
  int ITYPE, NUSIG;
  double kappa;
  int A_FF_model,V_FF_model;
  int Trans_Corr;
};

// Energy table for nue or numu (in GeV)
const double energy_table1[210] = {
  0.025, 0.075, 0.125, 0.175, 0.225, 
  0.275, 0.325, 0.375, 0.425, 0.475,  
  0.525, 0.575, 0.625, 0.675, 0.725,  
  0.775, 0.825, 0.875, 0.925, 0.975,  
  1.025, 1.075, 1.125, 1.175, 1.225,
  1.275, 1.325, 1.375, 1.425, 1.475,  
  1.525, 1.575, 1.625, 1.675, 1.725,  
  1.775, 1.825, 1.875, 1.925, 1.975,  
  2.025, 2.075, 2.125, 2.175, 2.225,  
  2.275, 2.325, 2.375, 2.425, 2.475, 
  2.525, 2.575, 2.625, 2.675, 2.725,  
  2.775, 2.825, 2.875, 2.925, 2.975,  
  3.025, 3.075, 3.125, 3.175, 3.225,  
  3.275, 3.325, 3.375, 3.425, 3.475,  
  3.525, 3.575, 3.625, 3.675, 3.725,  
  3.775, 3.825, 3.875, 3.925, 3.975,  
  4.025, 4.075, 4.125, 4.175, 4.225,  
  4.275, 4.325, 4.375, 4.425, 4.475,  
  4.525, 4.575, 4.625, 4.675, 4.725,  
  4.775, 4.825, 4.875, 4.925, 4.975, 
  5.025, 5.075, 5.125, 5.175, 5.225,  
  5.275, 5.325, 5.375, 5.425, 5.475,  
  5.525, 5.575, 5.625, 5.675, 5.725,  
  5.775, 5.825, 5.875, 5.925, 5.975,  
  6.025, 6.075, 6.125, 6.175, 6.225,  
  6.275, 6.325, 6.375, 6.425, 6.475,  
  6.525, 6.575, 6.625, 6.675, 6.725,  
  6.775, 6.825, 6.875, 6.925, 6.975,  
  7.025, 7.075, 7.125, 7.175, 7.225,  
  7.275, 7.325, 7.375, 7.425, 7.475, 
  7.525, 7.575, 7.625, 7.675, 7.725,  
  7.775, 7.825, 7.875, 7.925, 7.975,  
  8.025, 8.075, 8.125, 8.175, 8.225,  
  8.275, 8.325, 8.375, 8.425, 8.475,  
  8.525, 8.575, 8.625, 8.675, 8.725,  
  8.775, 8.825, 8.875, 8.925, 8.975,  
  9.025, 9.075, 9.125, 9.175, 9.225,  
  9.275, 9.325, 9.375, 9.425, 9.475,  
  9.525, 9.575, 9.625, 9.675, 9.725,  
  9.775, 9.825, 9.875, 9.925, 9.975, 
  11.00, 12.00, 13.00, 14.00, 15.00, 
  16.00, 17.00, 18.00, 19.00, 20.00
};

// Energy table for nutau (in GeV)
const double energy_table2[125] = {
  .270E+01, .280E+01, .290E+01, .300E+01, .310E+01, 
  .320E+01, .330E+01, .340E+01, .350E+01, .360E+01, 
  .370E+01, .380E+01, .390E+01, .400E+01, .410E+01, 
  .420E+01, .430E+01, .440E+01, .450E+01, .460E+01, 
  .470E+01, .480E+01, .490E+01, .500E+01, .510E+01, 
  .520E+01, .530E+01, .540E+01, .550E+01, .560E+01, 
  .570E+01, .580E+01, .590E+01, .600E+01, .610E+01, 
  .620E+01, .630E+01, .640E+01, .650E+01, .660E+01, 
  .670E+01, .680E+01, .690E+01, .700E+01, .710E+01, 
  .720E+01, .730E+01, .740E+01, .750E+01, .760E+01, 
  .770E+01, .780E+01, .790E+01, .800E+01, .810E+01, 
  .820E+01, .830E+01, .840E+01, .850E+01, .860E+01, 
  .870E+01, .880E+01, .890E+01, .900E+01, .910E+01, 
  .920E+01, .930E+01, .940E+01, .950E+01, .960E+01, 
  .970E+01, .980E+01, .990E+01, .100E+02, .100E+02, 
  .110E+02, .120E+02, .130E+02, .140E+02, .150E+02, 
  .160E+02, .170E+02, .180E+02, .190E+02, .200E+02, 
  .210E+02, .220E+02, .230E+02, .240E+02, .250E+02, 
  .260E+02, .270E+02, .280E+02, .290E+02, .300E+02, 
  .310E+02, .320E+02, .330E+02, .340E+02, .350E+02, 
  .360E+02, .370E+02, .380E+02, .390E+02, .400E+02, 
  .410E+02, .420E+02, .430E+02, .440E+02, .450E+02, 
  .460E+02, .470E+02, .480E+02, .490E+02, .500E+02, 
  .510E+02, .520E+02, .530E+02, .540E+02, .550E+02, 
  .560E+02, .570E+02, .580E+02, .590E+02, .600E+02
};
