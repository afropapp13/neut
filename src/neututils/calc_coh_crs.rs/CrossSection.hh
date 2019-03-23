# include <stdlib.h>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

const double PI     = 3.14159265;
const double HBARC  = 0.19732705;
const double GFERMI = 0.89719167E-7;
const double FPIFAC = 0.93;
const double PIMASS = 0.135;
const double PRMASS = 0.939;
const double AM     = 1.;
const double ATOMN  = 16.;
const double RAD0NU = 1.;
const double BRNTC2 = 1.E-24;
const double BRNTF2 = 1.E+2;
const double FMTCTM = 1.E-13;
const double P2OS   = 1./4.;

class CrossSection {

public:

  CrossSection(int nu_type, int current);
  ~CrossSection();

  int    SetParameters(const double energy, const double x1, const double x2, const double x3);
  double TotalCrossSection(const double energy, const int nloop);
  double DifferentialCrossSection(const double energy);
  double Coefficient(const double energy);
  void   Initialize();

private:

  int nu_type;
  int current;
  double mass_lep;
  double x, y, z;

  int ITYPE, CCNC;
  double ccnc_ratio;
  int LeptonMass;
};

extern "C" { double crosin_(double&); }
extern "C" { double crosto_(double&); }
extern "C" { double getr_(double&, double&, double&); }
extern "C" { double dbesi0_(double&); }
