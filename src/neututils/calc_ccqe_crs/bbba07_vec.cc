//////////////////////////////////////////////////////////////
//
// 
// 
// 
// Original comments start from here
//////////////////////////////////////////////////////////////
// This ROOT macro contains BBBA07 parameterizations of the
// elastic nucleon form factors
// version 0.2
// 07/30/07

// Usage example 
// $ root
// root [0] .L BBBA07.C
// root [1] drawGepBBBA07()
//////////////////////////////////////////////////////////////

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

extern "C"
{
  double gepBBBA07(double *x);
  double gmpBBBA07(double *x);
  double genBBBA07(double *x);
  double gmnBBBA07(double *x);

  int    bbba07_(double *x, 
				 double *gep, double *gmp, double *gen, double *gmn);
}

						  
// proton mass, GeV
static const double MP = 0.9383;

// proton magnetic moment
static const double MuP = 2.7928;

// neutron mass, GeV
static const double MN = 0.9396;

// neutron magnetic moment
static const double MuN = -1.913;

// nucleon mass for axial ff (MP+MN)/2
static const double MNucl = (MP+MN)/2.;

// standard value of LambdaD for dipole
static const double LambdaD = 0.71;

// squared axial mass
//static const double MA = 1.015;
// Use value from the card
/*
double const MA = nemdls_.xmaqe;
static const double MA2 = MA*MA;
*/


// 
static const double GA = 1.267;

// define interval of Q2
static const double Q2min = 0.0; // GeV2
static const double Q2max = 20.; // GeV2

// tau = Q2/4M**2
double tau(double qsq, double m)
{
    return qsq/(4.*m*m);
}

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
double dipole(double* x, const double* par)
{
    return par[1]*pow(1. + x[0]/par[0],-2);
}

// Galster parameterization
double galster(double* x, double* par)
{
    double tau1 = tau(x[0], par[4]);
    return par[0]*tau1/(1. + par[1]*tau1)*par[3]*pow(1. + x[0]/par[2],-2);
}

// "Galster factor" a*tau/(1+b*tau) 
double galsterFactor(double* x, const double* par)
{
    double tau1 = tau(x[0], par[2]);
    return par[0]*tau1/(1. + par[1]*tau1);
}


static const double ParDipole[2] = { LambdaD, 1. };

static const double ParGalster[5] = { 1.7, 3.3, LambdaD, 1., MN };

static const double ParGalsterFactor[3] = { 1.7, 3.3, MN };

// Kelly parameterization
double kelly(double* x, const double* par)
{
    double tau1 = tau(x[0], par[4]);
    double numerator = 1.;
    double denominator = 1.;
    numerator += par[0]*tau1;
    denominator += par[1]*tau1;
    denominator += par[2]*tau1*tau1;
    denominator += par[3]*tau1*tau1*tau1;
    return numerator/denominator;
}

static const double ParGepKelly[5] = { -0.24, 10.98, 12.82, 21.97, MP };
static const double ParGmpKelly[5] = { 0.1717, 11.26, 19.32, 8.33, MP };

// Lagrange parameterization
double lagrange(double* x, const double* par)
{
    static const int N = 7;
    static const double nodes[N] = { 0., 1./6., 2./6., 3./6., 4./6., 5./6., 1. };

    double qsq = x[0];
    double mass = par[7];
    double xi1 = xi(qsq, mass);

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

// Kelly Parameters
// proton electric
static const double ParGepLagrange[8] = { 1., .9927, .9898, .9975, .9812, .9340, 1., MP };
// proton magnetic
static const double ParGmpLagrange[8] = { 1., 1.0011, .9992, .9974, 1.0010, 1.0003,1., MP };

// Lagrange Parameters
// neutron magnetic
static const double ParGmnLagrange_25[8] = { 1., .9958, .9877, 1.0193, 1.0350, .9164, .7300, MN };
static const double ParGmnLagrange_43[8] = { 1., .9959, .9851, 1.0187, 1.0307, .9080, .9557, MN };

// neutron electric
static const double ParGenLagrange_25[8] = {1., 1.1011, 1.1392, 1.0203, 1.1093, 1.5429, 0.9706, MN};
static const double ParGenLagrange_43[8] = {1., 1.1019, 1.1387, 1.0234, 1.1046, 1.5395, 1.2708, MN };

// BBBA07

// Gep
double gep_BBBA07(double* x, const double* par)
{
    double lagr = lagrange(x, par);
    double kelly1 = kelly(x, ParGepKelly);
    return lagr*kelly1;
}

double gepBBBA07(double* x)
{
  double result;
  result = gep_BBBA07(x, ParGepLagrange);
  return result;
}

// Gmp
double gmp_BBBA07(double* x, const double* par)
{
    double lagr = lagrange(x, par);
    double kelly1 = kelly(x, ParGmpKelly);
    return lagr*kelly1;
}

double gmpBBBA07(double* x)
{
  double result;
  result = gmp_BBBA07(x,ParGmpLagrange);
  return result;
}
  

// Gmn
double gmn_BBBA07(double* x, const double* par)
{
    double lagr = lagrange(x, par);
    double gmp = gmp_BBBA07(x, ParGmpLagrange);
    return lagr*gmp;
}

double gmnBBBA07_25(double* x)
{
  double result;
  result = gmn_BBBA07(x, ParGmnLagrange_25);
  return result;
}

double gmnBBBA07_43(double* x)
{
  double result;
  result = gmn_BBBA07(x, ParGmnLagrange_43);
  return result;
}

double gmnBBBA07(double* x)
{
  double result;
  result = gmnBBBA07_25(x);
  return result;
}


// Gen
double gen_BBBA07(double* x, const double* par)
{
    double lagr = lagrange(x, par);
    double galst1 = galsterFactor(x, ParGalsterFactor);
    double gep = gep_BBBA07(x, ParGepLagrange);
    return lagr*galst1*gep;
}

double genBBBA07_25(double* x)
{
  double result;
  result = gen_BBBA07(x, ParGenLagrange_25);
  return result;
}

double genBBBA07_43(double* x)
{
  double result;
  result = gen_BBBA07(x, ParGenLagrange_43);
  return result;
}

double genBBBA07(double* x)
{
  double result;
  result = genBBBA07_25(x);
  return result;
}


// combined
int bbba07_(double *x, 
			double *gep, double *gmp, double *gen, double *gmn)
{
    *gep = gepBBBA07(x);
	*gmp = gmpBBBA07(x) * 2.7928;
	*gen = genBBBA07(x);
    *gmn = gmnBBBA07(x) * -1.9130;
	return 0;
}

#ifdef HAS_ROOT
// print values 
void printParam(TF1* param, double* q2, int n)
{
    for (int i=0; i<n; ++i) {
        cout << "// " << q2[i] << "     " << param->Eval(q2[i]) << endl;
    }
}


void printValues()
{
    TF1* gep_BBBA07 = new TF1("gep_BBBA07", gep_BBBA07, Q2min, Q2max, 8);
    gep_BBBA07->SetParameters(ParGepLagrange);

    TF1* gmp_BBBA07 = new TF1("gmp_BBBA07", gmp_BBBA07, Q2min, Q2max, 8);
    gmp_BBBA07->SetParameters(ParGmpLagrange);

    TF1* gmn_BBBA07_25 = new TF1("gmn_BBBA07_25", gmn_BBBA07, Q2min, Q2max, 8);
    gmn_BBBA07_25->SetParameters(ParGmnLagrange_25);

    TF1* gmn_BBBA07_43 = new TF1("gmn_BBBA07_43", gmn_BBBA07, Q2min, Q2max, 8);
    gmn_BBBA07_43->SetParameters(ParGmnLagrange_43);

    TF1* gen_BBBA07_25 = new TF1("gen_BBBA07_25", gen_BBBA07, Q2min, Q2max, 8);
    gen_BBBA07_25->SetParameters(ParGenLagrange_25);

    TF1* gen_BBBA07_43 = new TF1("gen_BBBA07_43", gen_BBBA07, Q2min, Q2max, 8);
    gen_BBBA07_43->SetParameters(ParGenLagrange_43);

    const double q2[4] = {0.1, 0.5, 1., 2.};

    cout << "// Q**2, GeV**2        value" << endl;
    cout << "// Gep:" << endl;
    printParam(gep_BBBA07, q2, 4);
    cout << "// Gmp:" << endl;
    printParam(gmp_BBBA07, q2, 4);
    cout << "// Gmn25:" << endl;
    printParam(gmn_BBBA07_25, q2, 4);
    cout << "// Gmn43:" << endl;
    printParam(gmn_BBBA07_43, q2, 4);
    cout << "// Gen25:" << endl; 
    printParam(gen_BBBA07_25, q2, 4);
    cout << "// Gen43:" << endl;
    printParam(gen_BBBA07_43, q2, 4);
}


// drawing examples
// BBBA07
void drawGep_BBBA07()
{
    TF1* gep_BBBA07 = new TF1("gep_BBBA07", gep_BBBA07, Q2min, Q2max, 8);
    gep_BBBA07->SetParameters(ParGepLagrange);
    gep_BBBA07->Draw();
}

#endif

// table of values of BBBA07 as printed by printValues()
// Q**2, GeV**2        value
// Gep:
// 0.1     0.742447
// 0.5     0.334549
// 1     0.162618
// 2     0.0537109
// Gmp:
// 0.1     0.752648
// 0.5     0.339299
// 1     0.176445
// 2     0.0725707
// Gmn25:
// 0.1     0.742646
// 0.5     0.347636
// 1     0.183311
// 2     0.0723006
// Gmn43:
// 0.1     0.741354
// 0.5     0.347576
// 1     0.182908
// 2     0.0714289
// Gen25:
// 0.1     0.0376581
// 0.5     0.0554865
// 1     0.0430783
// 2     0.0237363
// Gen43:
// 0.1     0.0376437
// 0.5     0.0556677
// 1     0.0430036
// 2     0.0234982
// Fa:
// 0.1     1.00519
// 0.5     0.63817
// 1     0.377064
// 2     0.186606

/* just for check
int
main(int argc, char **argv)
{
  double q2tbl[4]={0.1,0.5,1.0,2.0};
  int i;

  for ( i = 0 ; i < 4 ; i++ ){
	cout << gepBBBA07(&(q2tbl[i])) << endl;
  }
  for ( i = 0 ; i < 4 ; i++ ){
	cout << gmpBBBA07(&(q2tbl[i])) << endl;
  }
  for ( i = 0 ; i < 4 ; i++ ){
	cout << genBBBA07(&(q2tbl[i])) << endl;
  }
  for ( i = 0 ; i < 4 ; i++ ){
	cout << gmnBBBA07(&(q2tbl[i])) << endl;
  }

}
*/
