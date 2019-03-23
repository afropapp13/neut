# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

//
// This routine is based on the ROOT macro made by A.Bodek
//

using namespace std;

// proton mass, GeV
double const MP = 0.9383;

// proton magnetic moment
double const MuP = 2.7928;

// neutron mass, GeV
double const MN = 0.9396;

// neutron magnetic moment
double const MuN = -1.913;

// nucleon mass for axial ff (MP+MN)/2
double const MNucl = (MP+MN)/2.;

// standard value of LambdaD for dipole
double const LambdaD = 0.71;

// squared axial mass
//double const MA = 1.015;
double const MA = 1.21;
double const MA2 = MA*MA;

// 
//double const GA = -1.267;
double const GA = -1.232;

// Nachtman variable xi
double xi(double qsq, double m)
{
    if (qsq==0) {
        return 0.;
    } else {
        return 2./(1.+pow(1.+4.*m*m/qsq, 0.5));
    }
}

// dipole parameterization
double dipole(const double& x, const double* par)
{
    return par[1]*pow(1. + x/par[0],-2);
}

double const ParDipoleA[2] = { MA2, GA };


// Lagrange parameterization
double lagrange(const double& x, const double* par)
{
  int const N = 7;
    double const nodes[N] = { 0., 1./6., 2./6., 3./6., 4./6., 5./6., 1. };

    double qsq = x;
    double mass = par[7];
    double xi1;
    xi1 = xi(qsq, mass);

    double sum = 0.;
    for (int i = 0; i<N; ++i) {
        double part = par[i];
        for (int j = 0; j<N; ++j) {
            if ( i != j ) {
                part *= (xi1 - nodes[j])/(nodes[i]-nodes[j]);
            }
        }
        sum += part;
    }
    return sum;
}

// fa
double const ParFaLagrange[8] = {1., 0.9133, 0.9955, 1.1043, 1.1753, 1.3912, 0.7443, MNucl};

// BBBA07
// Fa
double faBBBA07(const double& x)
{  
  double lagr;
  double dip;
  lagr = lagrange(x, ParFaLagrange);
  dip = dipole(x, ParDipoleA);
  return lagr*dip;
}
